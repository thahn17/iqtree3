#!/usr/bin/env python3
from __future__ import annotations
import math,time,random,copy,csv,sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm
sys.path.insert(0,str(Path(__file__).resolve().parent))
import benchmark_spectral_bottleneck as bench
import multisite_pruned_general as pg
import multisite_lowerarm_taxonflow as base
import fractional_tree_hierarchical_units as hd
from phylo_branch_opt import TreeLikelihood
BASE='ACGT'

def comp_expr(M,u,port,pat,b):
    m,s,_=M._req_slot(u,port);d={}
    for t,ch in enumerate(pat):
        if ch==b:
            x=M.tmem[m,0,s,t];d[x]=d.get(x,0.0)+1.0
    return d

def ref_matrix(kind='single',tref=0.1,spread=2.0):
    if kind=='identity':return np.eye(4)
    if kind=='single':return expm(base.Q*float(tref))
    if kind=='logmean':
        lo=max(float(tref)/float(spread),1e-8);hi=float(tref)*float(spread)
        ts=np.geomspace(lo,hi,7);return np.mean(np.stack([expm(base.Q*t) for t in ts]),axis=0)
    if kind=='wide':
        ts=np.r_[0.0,np.geomspace(0.01,1.0,7),5.0];return np.mean(np.stack([expm(base.Q*t) for t in ts]),axis=0)
    raise ValueError(kind)

def unit_weight(M,u,scheme):
    m=int(M.units[u]['descendant_count']);a=int(M.units[u]['left_flow']);b=int(M.units[u]['right_flow'])
    if scheme=='equal':return 1.0
    if scheme=='inverse_flow':return 1/max(m,1)
    if scheme=='sqrt_inverse':return 1/math.sqrt(max(m,1))
    if scheme=='copy_normalized':return 1/max(len(M.by[m]),1)
    if scheme=='flow_capacity':return 1/max(M.L//max(m,1),1)
    if scheme=='balanced':return 4*a*b/(m*m)
    if scheme=='balanced_sqrt':return math.sqrt(max(4*a*b/(m*m),0.0))
    raise ValueError(scheme)

def add_overlap_scores(M,anchors=8,weight='equal',kind='single',tref=0.1,spread=2.0,baseline_center=True,score_mode='overlap'):
    P=ref_matrix(kind,tref,spread);npat=len(M.patterns);rank=sorted(range(npat),key=lambda p:(-M.weights[p],-len(set(M.patterns[p])),p));chosen=rank[:min(anchors,npat)]
    totalw=max(sum(float(M.weights[p]) for p in chosen),1.0);obj={};vars0=len(M.lp.names);rows0=len(M.lp.le)
    # baseline stationary-like score per active unit, subtracted to avoid rewarding uninformative fragments.
    m0=P@np.asarray(base.PI,float); baseline=float(np.asarray(base.PI)@m0) if baseline_center else 0.0
    for p in chosen:
        pat=M.patterns[p];pw=float(M.weights[p]/totalw)
        for u in M.splits:
            a=float(M.units[u]['left_flow']);b=float(M.units[u]['right_flow']);A=M.act[u];uw=unit_weight(M,u,weight)
            # transformed normalized child empirical state profiles; h_i <= both sides.
            for i in range(4):
                h=M.lp.var(f'ov:{p}:{u}:{i}',0,1)
                row={h:1.0}
                for bb in range(4):
                    coef=float(P[i,bb]/a)
                    if coef:
                        for v,z in comp_expr(M,u,'left',pat,BASE[bb]).items():row[v]=row.get(v,0.0)-coef*z
                M.lp.add_le(row,0)
                row={h:1.0}
                for bb in range(4):
                    coef=float(P[i,bb]/b)
                    if coef:
                        for v,z in comp_expr(M,u,'right',pat,BASE[bb]).items():row[v]=row.get(v,0.0)-coef*z
                M.lp.add_le(row,0)
                obj[h]=obj.get(h,0.0)+(1.0 if score_mode=='overlap' else -1.0)*pw*uw*float(base.PI[i])
            if baseline:
                obj[A]=obj.get(A,0.0)+((-baseline) if score_mode=='overlap' else baseline)*pw*uw
    return obj,dict(score_vars=len(M.lp.names)-vars0,score_rows=len(M.lp.le)-rows0,anchors=len(chosen),baseline=baseline)

def mix_result(primary,secondary,alpha):
    rr=copy.copy(secondary);x=secondary.x.copy();n=len(primary.x);x[:n]=(1-alpha)*primary.x+alpha*secondary.x[:n];rr.x=x;return rr

def rescore(M,res,seqs,model,true,K=25,multistart=1):
    labels=[f't{i}' for i in range(M.L)];C,F=hd.rank_integral_tree_candidates(M,res,topk=K,root_beam=384,leaf_beam=96,child_branching=18,leaf_branching=12,merge_branching=160,max_expansions=700000)
    vals=[];tr=9999
    for i,z in enumerate(C,1):
        if bench.key(z.topology)==bench.key(true):tr=i
        ll=TreeLikelihood(z.topology,labels,seqs,model,decoded=z,frac=F,short_time=pg.base.branch_times,long_ratio=2).optimize(multistart=multistart,maxiter=130,ftol=1e-8).log_likelihood
        vals.append((ll,i,z.proportion))
    vals.sort(reverse=True);return (vals[0] if vals else (-math.inf,0,0)),len(C),tr

def run(n=8,S=40,rep=0,regime='moderate',delta=2,weight='equal',anchors=8,kind='single',tref=0.1,spread=2.0,score_mode='overlap',alphas=(0,0.5,1)):
    seed=810000+n*1000+S*10+rep+(100000 if regime=='wide' else 0);rng=random.Random(seed);sim=bench.gtr_f(bench.PI,bench.RATES);true=bench.random_tree(n,rng);seqs,cols,_=bench.simulate(true,S,sim,rng,regime);model=bench.empirical_model(cols);bench.configure(model);lib=bench.build_lib(n,'/mnt/data/spectral_bottleneck_libs')
    tb=time.time();M=pg.PrunedGeneralModel(lib,cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='spectral',spectral_support_points=7);build=time.time()-tb;z,r,ps=M.solve();basebest,nc,tr=rescore(M,r,seqs,model,true)
    ts=time.time();obj,stats=add_overlap_scores(M,anchors,weight,kind,tref,spread,True,score_mode);scorebuild=time.time()-ts;raw=sum(v*r.x[j] for j,v in M.obj.items());M.lp.add_ge(M.obj,raw-delta);sv,r2,ss=M.lp.solve(obj,True);true_ll=bench.rescore(true,[f't{i}' for i in range(n)],seqs,model)
    rows=[]
    for a in alphas:
        rr=mix_result(r,r2,float(a));best,cands,trank=rescore(M,rr,seqs,model,true)
        rows.append(dict(n=n,S=S,rep=rep,regime=regime,delta=delta,weight=weight,kind=kind,tref=tref,spread=spread,score_mode=score_mode,alpha=float(a),lp_score=z,secondary_score=sv,build_s=build,primary_solve_s=ps,score_build_s=scorebuild,secondary_solve_s=ss,baseline_best=basebest[0],best_logL=best[0],true_logL=true_ll,best_minus_true=best[0]-true_ll,candidates=cands,true_rank=trank,**stats))
    return rows
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=8);ap.add_argument('--S',type=int,default=40);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--regime',default='moderate');ap.add_argument('--delta',type=float,default=2);ap.add_argument('--weight',default='equal');ap.add_argument('--anchors',type=int,default=8);ap.add_argument('--kind',default='single');ap.add_argument('--tref',type=float,default=.1);ap.add_argument('--spread',type=float,default=2);ap.add_argument('--score-mode',default='overlap');ap.add_argument('--out',default='/mnt/data/overlap_local_score.csv');a=ap.parse_args();rows=run(a.n,a.S,a.rep,a.regime,a.delta,a.weight,a.anchors,a.kind,a.tref,a.spread,a.score_mode);print(rows,flush=True);p=Path(a.out);ex=p.exists()
 with p.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));
  if not ex:w.writeheader()
  w.writerows(rows)
