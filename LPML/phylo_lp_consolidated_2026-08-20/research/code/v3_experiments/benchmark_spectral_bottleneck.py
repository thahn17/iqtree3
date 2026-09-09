#!/usr/bin/env python3
from __future__ import annotations
import sys, math, random, time, csv
from pathlib import Path
import numpy as np
from scipy.linalg import expm
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phylo_substitution import gtr_f, set_active_model
from fixed_flow_unit_generator_variable_lengths import LibraryBuilder
import multisite_pruned_general as pg
import multisite_variable_length_compressed as vl
import multisite_lowerarm_taxonflow as low
import test_rootonly_pruning as tr
import fractional_tree_hierarchical_units as hd
from phylo_branch_opt import TreeLikelihood
RATES=(1.0,2.3,0.6,0.85,1.5,0.95); PI=np.array([0.30,0.20,0.25,0.25],float); BASE='ACGT'
def taxa(t): return (t,) if isinstance(t,int) else tuple(sorted(taxa(t[0])+taxa(t[1])))
def canon(t):
    if isinstance(t,int): return t
    a,b=canon(t[0]),canon(t[1]); return (a,b) if taxa(a)<=taxa(b) else (b,a)
def key(t): return repr(canon(t))
def random_tree(n,rng):
    q=list(range(n))
    while len(q)>1:
        i,j=sorted(rng.sample(range(len(q)),2),reverse=True); a=q.pop(i);b=q.pop(j);q.append(canon((a,b)))
    return canon(q[0])
def edge_paths(t,path=()):
    if isinstance(t,int): return []
    out=[]
    for i,ch in enumerate(t): p=path+(i,);out.append(p);out.extend(edge_paths(ch,p))
    return out
def simulate(t,S,model,rng,regime):
    paths=edge_paths(t)
    if regime=='moderate': lengths={p:rng.uniform(0.05,0.25) for p in paths}
    else:
        grid=list(np.geomspace(1e-4,10.0,len(paths)));rng.shuffle(grid);lengths={p:float(grid[i]) for i,p in enumerate(paths)}
    P={p:expm(model.q*lengths[p]) for p in paths};n=len(taxa(t));seq=[[] for _ in range(n)]
    for _ in range(S):
        rs=int(rng.choices(range(4),weights=model.pi,k=1)[0])
        def rec(node,state,path=()):
            if isinstance(node,int): seq[node].append(BASE[state]);return
            for i,ch in enumerate(node):
                p=path+(i,);ns=int(rng.choices(range(4),weights=P[p][state],k=1)[0]);rec(ch,ns,p)
        rec(t,rs)
    seqs=[''.join(x) for x in seq];cols=[tuple(seqs[i][s] for i in range(n)) for s in range(S)];return seqs,cols,lengths
def build_lib(n,root):
    d=Path(root)/f'lib{n}';f=d/'library.json'
    if not f.exists():
        d.mkdir(parents=True,exist_ok=True);b=LibraryBuilder([f't{i}' for i in range(n)],variable_length_paths=True,long_path_edges=2,variable_length_scope='branches');b.build();b.write_library(d)
    return str(f)
def configure(model):
    set_active_model(model);vl.configure_substitution_model(model);low.configure_substitution_model(model);tr.base.configure_substitution_model(model);tr.PI=model.pi;pg.base.configure_substitution_model(model);pg.PI=model.pi
def empirical_model(cols):
    cnt=np.zeros(4)
    for c in cols:
        for ch in c: cnt[BASE.index(ch)]+=1
    return gtr_f((cnt+0.5)/(cnt.sum()+2.0),RATES)
def rescore(t,labels,seqs,model,decoded=None,frac=None):
    return TreeLikelihood(t,labels,seqs,model,decoded=decoded,frac=frac,short_time=pg.base.branch_times,long_ratio=2.0).optimize(multistart=4,maxiter=180,ftol=1e-9).log_likelihood
def one(n,S,regime,rep,mode,root,topk=40,support=5):
    seed=810000+n*1000+S*10+rep+(0 if regime=='moderate' else 100000);rng=random.Random(seed);sim=gtr_f(PI,RATES);true=random_tree(n,rng);seqs,cols,lengths=simulate(true,S,sim,rng,regime);model=empirical_model(cols);configure(model);lib=build_lib(n,root);labels=[f't{i}' for i in range(n)]
    t=time.time();M=pg.PrunedGeneralModel(lib,cols,preprocess_mode='fast',substitution_model=model,branch_relaxation=mode,spectral_support_points=support);build=time.time()-t;t=time.time();score,res,_=M.solve();solve=time.time()-t
    t=time.time();cands,frac=hd.rank_integral_tree_candidates(M,res,topk=topk,root_beam=256 if n<=12 else 384,leaf_beam=64 if n<=12 else 96,child_branching=16,leaf_branching=10,merge_branching=128,max_expansions=500000 if n<=12 else 800000);ident=time.time()-t
    vals=[];trank=9999;t=time.time()
    for i,z in enumerate(cands,1):
        if key(z.topology)==key(true):trank=i
        vals.append((rescore(z.topology,labels,seqs,model,z,frac),i,z.proportion))
    bt=time.time()-t;vals.sort(reverse=True);best=vals[0] if vals else (-math.inf,0,0);true_ll=rescore(true,labels,seqs,model)
    return dict(n=n,S=S,regime=regime,rep=rep,branch_mode=mode,topk=topk,support_points=support,split_units=len(M.splits),lp_vars=len(M.lp.names),build_s=build,solve_s=solve,tree_id_s=ident,branch_rescore_s=bt,lp_score=score,candidates=len(cands),true_structural_rank=trank,best_structural_rank=best[1],best_bottleneck=best[2],best_logL=best[0],true_topology_logL=true_ll,best_minus_true=best[0]-true_ll,true_min_branch=min(lengths.values()),true_max_branch=max(lengths.values()))
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--S',type=int,required=True);ap.add_argument('--regime',choices=['moderate','wide'],required=True);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--mode',choices=['two-point','spectral'],required=True);ap.add_argument('--out',required=True);ap.add_argument('--root',default='/mnt/data/spectral_bottleneck_libs');a=ap.parse_args();r=one(a.n,a.S,a.regime,a.rep,a.mode,a.root);print(r,flush=True)
 exists=Path(a.out).exists();
 with open(a.out,'a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(r));
  if not exists:w.writeheader()
  w.writerow(r)
