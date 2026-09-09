#!/usr/bin/env python3
from __future__ import annotations
import math,time,random,copy,csv,sys
from pathlib import Path
import numpy as np
sys.path[:0]=['/mnt/data/work_v3','/mnt/data/_v3','/mnt/data']
import benchmark_spectral_bottleneck as bench
import multisite_pruned_general as pg
import multisite_lowerarm_taxonflow as base
import fractional_tree_hierarchical_units as hd
from phylo_branch_opt import TreeLikelihood
from test_overlap_local_score import ref_matrix, comp_expr, mix_result
BASE='ACGT'

def pair_dist(seqs):
    X=np.frombuffer(''.join(seqs).encode(),dtype='S1').reshape(len(seqs),len(seqs[0]))
    n=len(seqs); D=np.zeros((n,n),float)
    for i in range(n):
        # vectorized against all j>i
        if i+1<n:
            vals=np.mean(X[i+1:]!=X[i],axis=1)
            D[i,i+1:]=vals; D[i+1:,i]=vals
    return D

def unit_pair_modifiers(M,res,D,eta=0.5):
    vals={}; bytype={}
    for u in M.splits:
        act=float(res.x[M.act[u]])
        m=int(M.units[u]['descendant_count']); a=int(M.units[u]['left_flow']); b=int(M.units[u]['right_flow'])
        if act<1e-7:
            vals[u]=None; continue
        xl=np.array([res.x[M.tmem[M._req_slot(u,'left')[0],0,M._req_slot(u,'left')[1],t]] for t in range(M.L)],float)/act
        xr=np.array([res.x[M.tmem[M._req_slot(u,'right')[0],0,M._req_slot(u,'right')[1],t]] for t in range(M.L)],float)/act
        cross=float(xl@D@xr/max(a*b,1))
        if a>1:
            wa=float((xl@D@xl)/(a*(a-1))) # D diag=0, ordered pairs denominator
        else: wa=0.0
        if b>1:
            wb=float((xr@D@xr)/(b*(b-1)))
        else: wb=0.0
        s=cross-max(wa,wb)
        vals[u]=s; bytype.setdefault((m,a,b),[]).append(s)
    allv=np.array([v for v in vals.values() if v is not None],float)
    med=float(np.median(allv)) if len(allv) else 0.0
    mad=float(np.median(np.abs(allv-med))) if len(allv) else 1.0
    scale=max(1.4826*mad,float(np.std(allv)) if len(allv) else 0.0,1e-6)
    out={}
    for u in M.splits:
        v=vals[u]
        if v is None: z=0.0
        else: z=np.clip((v-med)/scale,-2.0,2.0)
        out[u]=float(np.exp(eta*z))
    return out,dict(pair_median=med,pair_scale=scale,pair_active=sum(v is not None for v in vals.values()))

def structural_weight(M,u,beta=0.0,gamma=0.0):
    m=float(M.units[u]['descendant_count']); a=float(M.units[u]['left_flow']); b=float(M.units[u]['right_flow'])
    B=max(4*a*b/(m*m),1e-6)
    return (m**beta)*(B**gamma)

def add_weighted_overlap(M,pairmods=None,beta=0.0,gamma=0.0,anchors=8,tref=.1):
    P=ref_matrix('single',tref,2.0); npat=len(M.patterns)
    rank=sorted(range(npat),key=lambda p:(-M.weights[p],-len(set(M.patterns[p])),p)); chosen=rank[:min(anchors,npat)]
    totalw=max(sum(float(M.weights[p]) for p in chosen),1.0); obj={}; vars0=len(M.lp.names); rows0=len(M.lp.le)
    m0=P@np.asarray(base.PI,float); baseline=float(np.asarray(base.PI)@m0)
    # normalize structural weights geometrically over split units so scale is comparable across beta/gamma.
    sw=np.array([structural_weight(M,u,beta,gamma) for u in M.splits],float)
    gmean=float(np.exp(np.mean(np.log(np.clip(sw,1e-12,None))))) if len(sw) else 1.0
    smap={u:structural_weight(M,u,beta,gamma)/gmean for u in M.splits}
    for p in chosen:
        pat=M.patterns[p]; pw=float(M.weights[p]/totalw)
        for u in M.splits:
            a=float(M.units[u]['left_flow']); b=float(M.units[u]['right_flow']); A=M.act[u]
            uw=smap[u]*(pairmods[u] if pairmods is not None else 1.0)
            for i in range(4):
                h=M.lp.var(f'wov:{p}:{u}:{i}',0,1)
                row={h:1.0}
                for bb in range(4):
                    coef=float(P[i,bb]/a)
                    if coef:
                        for v,z in comp_expr(M,u,'left',pat,BASE[bb]).items(): row[v]=row.get(v,0.0)-coef*z
                M.lp.add_le(row,0)
                row={h:1.0}
                for bb in range(4):
                    coef=float(P[i,bb]/b)
                    if coef:
                        for v,z in comp_expr(M,u,'right',pat,BASE[bb]).items(): row[v]=row.get(v,0.0)-coef*z
                M.lp.add_le(row,0)
                obj[h]=obj.get(h,0.0)+pw*uw*float(base.PI[i])
            obj[A]=obj.get(A,0.0)-baseline*pw*uw
    return obj,dict(score_vars=len(M.lp.names)-vars0,score_rows=len(M.lp.le)-rows0,anchors=len(chosen))

def rescore_best(M,res,seqs,model,true,K=40):
    if M.L>=20: K=min(K,15)
    labels=[f't{i}' for i in range(M.L)]
    C,F=hd.rank_integral_tree_candidates(M,res,topk=K,root_beam=384,leaf_beam=96,child_branching=18,leaf_branching=12,merge_branching=160,max_expansions=800000)
    vals=[]; tr=9999
    for i,z in enumerate(C,1):
        if bench.key(z.topology)==bench.key(true): tr=i
        ll=TreeLikelihood(z.topology,labels,seqs,model,decoded=z,frac=F,short_time=pg.base.branch_times,long_ratio=2).optimize(multistart=1,maxiter=90,ftol=1e-7).log_likelihood
        vals.append((ll,i,z.proportion))
    vals.sort(reverse=True)
    return (vals[0] if vals else (-math.inf,0,0)),len(C),tr

def run(n=8,S=40,rep=0,regime='moderate',delta=2.0,beta=-.2,gamma=-.2,eta=.5,anchors=8,alphas=(0,.5,1)):
    seed=810000+n*1000+S*10+rep+(100000 if regime=='wide' else 0); rng=random.Random(seed)
    sim=bench.gtr_f(bench.PI,bench.RATES); true=bench.random_tree(n,rng); seqs,cols,_=bench.simulate(true,S,sim,rng,regime); model=bench.empirical_model(cols); bench.configure(model); lib=bench.build_lib(n,'/mnt/data/spectral_bottleneck_libs')
    # old two-point baseline
    M0=pg.PrunedGeneralModel(lib,cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='two-point',spectral_support_points=7); z0,r0,_=M0.solve(); b0,c0,t0=rescore_best(M0,r0,seqs,model,true)
    # spectral primary
    M=pg.PrunedGeneralModel(lib,cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='spectral',spectral_support_points=7); z,r,_=M.solve(); bp,cp,tp=rescore_best(M,r,seqs,model,true)
    D=pair_dist(seqs); mods,pstats=unit_pair_modifiers(M,r,D,eta)
    obj,stats=add_weighted_overlap(M,mods,beta,gamma,anchors,.1)
    raw=sum(v*r.x[j] for j,v in M.obj.items()); M.lp.add_ge(M.obj,raw-delta); sv,r2,ss=M.lp.solve(obj,True)
    rows=[]
    best_union=(-math.inf,0,0); union_alpha=None
    for a in alphas:
        rr=mix_result(r,r2,float(a)); best,cands,tr=rescore_best(M,rr,seqs,model,true)
        if best[0]>best_union[0]: best_union=best; union_alpha=float(a)
        rows.append(dict(n=n,S=S,rep=rep,regime=regime,delta=delta,beta=beta,gamma=gamma,eta=eta,alpha=float(a),two_point_best=b0[0],spectral_primary_best=bp[0],weighted_best=best[0],best_union_so_far=best_union[0],best_union_alpha=union_alpha,two_point_candidates=c0,spectral_candidates=cp,weighted_candidates=cands,two_point_true_rank=t0,spectral_true_rank=tp,weighted_true_rank=tr,lp_two_point=z0,lp_spectral=z,secondary_score=sv,**stats,**pstats))
    return rows

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser(); ap.add_argument('--n',type=int,default=8);ap.add_argument('--S',type=int,default=40);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--regime',default='moderate');ap.add_argument('--delta',type=float,default=2);ap.add_argument('--beta',type=float,default=-.2);ap.add_argument('--gamma',type=float,default=-.2);ap.add_argument('--eta',type=float,default=.5);ap.add_argument('--anchors',type=int,default=8);ap.add_argument('--out',default='/mnt/data/weighted_secondary_coupling.csv');a=ap.parse_args();rows=run(a.n,a.S,a.rep,a.regime,a.delta,a.beta,a.gamma,a.eta,a.anchors);print(rows,flush=True);p=Path(a.out);ex=p.exists()
 with p.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));
  if not ex:w.writeheader()
  w.writerows(rows)
