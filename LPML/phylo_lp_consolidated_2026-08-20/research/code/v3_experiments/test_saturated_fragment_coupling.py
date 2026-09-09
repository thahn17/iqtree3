#!/usr/bin/env python3
from __future__ import annotations
import sys,time,random,math,csv
from pathlib import Path
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import benchmark_spectral_bottleneck as bench
import multisite_pruned_general as pg
import multisite_lowerarm_taxonflow as base
import fractional_tree_hierarchical_units as hd
from phylo_branch_opt import TreeLikelihood

def count_expr(M,u,port,pat,b):
    m,s,_=M._req_slot(u,port);d={}
    for t,ch in enumerate(pat):
        if ch==b:
            v=M.tmem[m,0,s,t];d[v]=d.get(v,0)+1
    return d,m

def build_satvar(M,u,port,pat,b,T,name):
    d,flow=count_expr(M,u,port,pat,b);A=M.act[u]
    if flow<=T:
        # exact count; return expression rather than new variable
        return d
    v=M.lp.var(name,0,T)
    # v <= count; v <= T A; v >= (T/flow) count
    row={v:1};
    for j,c in d.items():row[j]=row.get(j,0)-c
    M.lp.add_le(row,0);M.lp.add_le({v:1,A:-T},0)
    row={v:-1}
    for j,c in d.items():row[j]=row.get(j,0)+(T/flow)*c
    M.lp.add_le(row,0)
    return {v:1.0}

def cloud(M,I,g,a,b,T):
    raw=I['raw'];VA=M._branch_vertices(a);VB=M._branch_vertices(b);X=[];Y=[];null_total=M.L-sum(g)
    for (ma,cL),AL in raw.items():
      if ma!=a:continue
      for (mb,cR),AR in raw.items():
       if mb!=b:continue
       ct=tuple(cL[k]+cR[k] for k in range(4))
       if any(ct[k]>g[k] for k in range(4)) or (a+b)-sum(ct)>null_total:continue
       for zA,PA in VA:
        LA=PA@AL[:,1]
        for zB,PB in VB:
         q=float(base.PI@(LA*(PB@AR[:,1])))
         feat=[min(cL[k],T) for k in range(4)]+[min(cR[k],T) for k in range(4)]+list(zA)+list(zB)
         X.append(feat);Y.append(q)
    Y=np.asarray(Y,float);mx=max(float(Y.max()),1e-300);return np.asarray(X,float),Y/mx

def add_scores(M,T=2,anchors=6,facets=3):
    P=len(M.patterns);rank=sorted(range(P),key=lambda p:(-M.weights[p],-len(set(M.patterns[p])),p));chosen=rank[:anchors]
    cache={};obj={};vars0=len(M.lp.names);rows0=len(M.lp.le);totalw=sum(M.weights[p] for p in chosen)
    for p in chosen:
      pat=M.patterns[p];I=M.info[p];g=tuple(I['counts'])
      for u in M.splits:
       a=int(M.units[u]['left_flow']);b=int(M.units[u]['right_flow']);key=(g,a,b,T,facets)
       if key not in cache:
        X,y=cloud(M,I,g,a,b,T);cache[key]=pg.shifted_majorant_planes(X,y,K=facets,fit_cap=1600)
       fs=cache[key]
       q=M.lp.var(f'satfrag:{T}:{p}:{u}',0,1);sat=[]
       for port in ('left','right'):
        for bb in 'ACGT':sat.append(build_satvar(M,u,port,pat,bb,T,f'sat:{T}:{p}:{u}:{port}:{bb}'))
       lv=M._branch_vars(u,'left');rv=M._branch_vars(u,'right')
       for f in fs:
        a0,*co=f;row={q:1,M.act[u]:-a0};off=0
        for d in sat:
         c=co[off];off+=1
         if c:
          for v,z in d.items():row[v]=row.get(v,0)-c*z
        for v in lv+rv:
         c=co[off];off+=1
         if c:row[v]=row.get(v,0)-c
        M.lp.add_le(row,0)
       obj[q]=float(M.weights[p]/totalw)
    return obj,dict(vars=len(M.lp.names)-vars0,rows=len(M.lp.le)-rows0,chosen=chosen,classes=len(cache))

def rescore(M,res,seqs,model,true,K=15):
 labels=[f't{i}' for i in range(M.L)];C,F=hd.rank_integral_tree_candidates(M,res,topk=K,root_beam=256,leaf_beam=64,child_branching=16,leaf_branching=10,merge_branching=128,max_expansions=500000);vals=[]
 for i,z in enumerate(C,1):
  ll=TreeLikelihood(z.topology,labels,seqs,model,decoded=z,frac=F,short_time=pg.base.branch_times,long_ratio=2).optimize(multistart=1,maxiter=100,ftol=1e-8).log_likelihood;vals.append((ll,i))
 vals.sort(reverse=True);return vals[0] if vals else (-math.inf,0),len(C)

def one(n=8,S=40,regime='moderate',rep=0,T=2,delta=2,anchors=6,facets=3):
 seed=810000+n*1000+S*10+rep+(0 if regime=='moderate' else 100000);rng=random.Random(seed);sim=bench.gtr_f(bench.PI,bench.RATES);true=bench.random_tree(n,rng);seqs,cols,_=bench.simulate(true,S,sim,rng,regime);model=bench.empirical_model(cols);bench.configure(model);lib=bench.build_lib(n,'/mnt/data/spectral_bottleneck_libs')
 M=pg.PrunedGeneralModel(lib,cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='spectral',spectral_support_points=7);z,r,_=M.solve();basebest,nc=rescore(M,r,seqs,model,true)
 t=time.time();obj,st=add_scores(M,T,anchors,facets);bt=time.time()-t;raw=sum(v*r.x[j] for j,v in M.obj.items());M.lp.add_ge(M.obj,raw-delta);sv,r2,sd=M.lp.solve(obj,True);best,nc2=rescore(M,r2,seqs,model,true);tl=bench.rescore(true,[f't{i}' for i in range(n)],seqs,model)
 return dict(n=n,S=S,regime=regime,rep=rep,T=T,delta=delta,baseline=basebest[0],best=best[0],improvement=best[0]-basebest[0],true_ll=tl,gap=best[0]-tl,candidates=nc2,score=sv,build_score_s=bt,solve_score_s=sd,**st)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=8);p.add_argument('--S',type=int,default=40);p.add_argument('--regime',default='moderate');p.add_argument('--rep',type=int,default=0);p.add_argument('--T',type=int,default=2);p.add_argument('--delta',type=float,default=2);p.add_argument('--out',default='/mnt/data/saturated_fragment.csv');a=p.parse_args();r=one(a.n,a.S,a.regime,a.rep,a.T,a.delta);print(r);fn=Path(a.out);ex=fn.exists();
 with fn.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(r));
  if not ex:w.writeheader();
  w.writerow(r)
