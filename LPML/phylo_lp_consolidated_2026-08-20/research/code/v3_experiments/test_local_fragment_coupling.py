#!/usr/bin/env python3
from __future__ import annotations
import sys, math, time, random, csv
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.linalg import expm
sys.path.insert(0,str(Path(__file__).resolve().parent))
import benchmark_spectral_bottleneck as bench
import multisite_pruned_general as pg
import multisite_lowerarm_taxonflow as base
import fractional_tree_hierarchical_units as hd
from phylo_branch_opt import TreeLikelihood

BASE='ACGT'

def _comp_expr(M,u,port,pat):
    m,s,_=M._req_slot(u,port)
    out=[]
    for b in 'ACG':
        d={}
        for t,ch in enumerate(pat):
            if ch==b:
                v=M.tmem[m,0,s,t];d[v]=d.get(v,0.0)+1.0
        out.append(d)
    return out


def _physical_score_vertices(M,flow,grid_n=7):
    if base.branch_relaxation_mode()!='spectral': return M._branch_vertices(flow)
    S=base.spectral_relaxation();mu=np.asarray(S.mu,float)
    # Expected-rate-normalized Q: cover near-zero, ordinary, long, and stationary endpoint.
    core=np.geomspace(1e-3,5.0,max(int(grid_n)-2,1)) if grid_n>2 else np.array([])
    times=np.r_[0.0,core,np.inf]
    out=[]
    for t in times:
        if np.isinf(t): z=np.zeros(len(mu));P=np.ones((4,1))*base.PI[None,:]
        else: z=np.exp(-mu*t);P=expm(base.Q*t)
        out.append((np.asarray(z,float),np.asarray(P,float)))
    return out

def _cloud_for_type(M,I,g,a,b,physical_grid=None):
    raw=I['raw']; VA=_physical_score_vertices(M,a,physical_grid) if physical_grid else M._branch_vertices(a); VB=_physical_score_vertices(M,b,physical_grid) if physical_grid else M._branch_vertices(b)
    FA=np.asarray([z for z,P in VA],float); FB=np.asarray([z for z,P in VB],float)
    PA=np.asarray([P for z,P in VA],float); PB=np.asarray([P for z,P in VB],float)
    na,nb=len(VA),len(VB); X=[];Y=[]
    null_total=M.L-sum(g)
    for (ma,cL),AL in raw.items():
        if ma!=a: continue
        for (mb,cR),AR in raw.items():
            if mb!=b: continue
            ctot=tuple(cL[k]+cR[k] for k in range(4))
            if any(ctot[k]>g[k] for k in range(4)): continue
            if (a+b)-sum(ctot)>null_total: continue
            av=AL[:,1]; bv=AR[:,1]
            LA=np.einsum('vij,j->vi',PA,av); LB=np.einsum('vij,j->vi',PB,bv)
            q=(LA*base.PI[None,:])@LB.T
            for ia in range(na):
                for ib in range(nb):
                    X.append([cL[0],cL[1],cL[2],cR[0],cR[1],cR[2],*FA[ia].tolist(),*FB[ib].tolist()])
                    Y.append(float(q[ia,ib]))
    if not Y: return np.empty((0,6+2*base.branch_feature_dim())),np.empty(0)
    Y=np.asarray(Y,float); mx=max(float(Y.max()),1e-300); return np.asarray(X,float),Y/mx


def _sample_cloud_for_type(M,I,g,a,b,samples=256,seed=0):
    raw=I['raw']; VA=M._branch_vertices(a); VB=M._branch_vertices(b)
    left=[(c,A) for (m,c),A in raw.items() if m==a]; right=[(c,A) for (m,c),A in raw.items() if m==b]
    if not left or not right:return np.empty((0,6+2*base.branch_feature_dim())),np.empty(0)
    rng=np.random.default_rng(seed); X=[];Y=[];null_total=M.L-sum(g);tries=0;target=int(samples)
    while len(Y)<target and tries<target*30:
        tries+=1;cL,AL=left[int(rng.integers(len(left)))];cR,AR=right[int(rng.integers(len(right)))]
        ctot=tuple(cL[k]+cR[k] for k in range(4))
        if any(ctot[k]>g[k] for k in range(4)) or (a+b)-sum(ctot)>null_total:continue
        ia=int(rng.integers(len(VA)));ib=int(rng.integers(len(VB)));zA,PA=VA[ia];zB,PB=VB[ib]
        q=float(base.PI@((PA@AL[:,1])*(PB@AR[:,1])))
        X.append([cL[0],cL[1],cL[2],cR[0],cR[1],cR[2],*np.asarray(zA).tolist(),*np.asarray(zB).tolist()]);Y.append(q)
    if not Y:return np.empty((0,6+2*base.branch_feature_dim())),np.empty(0)
    Y=np.asarray(Y,float);mx=max(float(Y.max()),1e-300);return np.asarray(X,float),Y/mx

def add_local_fragment_scores(M, anchors=8, facets=4, size_weight='equal', max_flow=None, sampled_points=None, physical_grid=None):
    # deterministic informative pattern selection: multiplicity then number of distinct states, with spread.
    P=len(M.patterns)
    rank=sorted(range(P),key=lambda p:(-M.weights[p],-len(set(M.patterns[p])),p))
    if anchors is None or anchors>=P: chosen=rank
    else:
        # include high-weight first half, then evenly spread through remaining ranks
        h=max(1,anchors//2); chosen=rank[:h]
        rem=rank[h:]
        if len(chosen)<anchors and rem:
            idx=np.linspace(0,len(rem)-1,anchors-len(chosen),dtype=int)
            chosen += [rem[i] for i in idx]
        chosen=list(dict.fromkeys(chosen))[:anchors]
    cache={}; qvars={}; obj={}; rows0=len(M.lp.le); vars0=len(M.lp.names)
    totalw=max(float(sum(M.weights[p] for p in chosen)),1.0)
    for p in chosen:
        pat=M.patterns[p]; I=M.info[p]; g=tuple(I['counts'])
        for u in M.splits:
            m=int(M.units[u]['descendant_count']); a=int(M.units[u]['left_flow']); b=int(M.units[u]['right_flow'])
            if max_flow is not None and m>max_flow: continue
            key=(g,a,b,base.branch_relaxation_mode(),base.SPECTRAL_SUPPORT_POINTS,facets)
            if key not in cache:
                if sampled_points:
                    seed=(hash((g,a,b,p)) & 0xffffffff);X,y=_sample_cloud_for_type(M,I,g,a,b,sampled_points,seed)
                else:X,y=_cloud_for_type(M,I,g,a,b,physical_grid)
                cache[key]=pg.shifted_majorant_planes(X,y,K=facets,fit_cap=1600)
            fs=cache[key]
            if not fs: continue
            q=M.lp.var(f'localfrag:{p}:{u}',-1 if sampled_points else 0,1);qvars[p,u]=q
            le=_comp_expr(M,u,'left',pat); re=_comp_expr(M,u,'right',pat)
            lv=M._branch_vars(u,'left'); rv=M._branch_vars(u,'right')
            for f in fs:
                a0,*co=f; row={q:1,M.act[u]:-float(a0)}; off=0
                for d in le+re:
                    c=float(co[off]);off+=1
                    if c:
                        for v,z in d.items(): row[v]=row.get(v,0.0)-c*z
                for v in lv+rv:
                    c=float(co[off]);off+=1
                    if c: row[v]=row.get(v,0.0)-c
                M.lp.add_le(row,0)
            sw=1.0
            if size_weight=='inverse_flow': sw=1.0/max(m,1)
            elif size_weight=='inverse_internal': sw=1.0/max(m-1,1)
            obj[q]=float(M.weights[p]/totalw)*sw
    return obj,dict(anchor_ids=chosen,variables=len(M.lp.names)-vars0,rows=len(M.lp.le)-rows0,qvars=len(qvars),cache_classes=len(cache))

def solve_secondary(M, primary_res, local_obj, delta=0.0):
    raw=sum(v*primary_res.x[j] for j,v in M.obj.items())
    M.lp.add_ge(M.obj,raw-float(delta))
    return M.lp.solve(local_obj,True)

def rescore_candidates(M,res,seqs,model,true,topk=20,multistart=2):
    labels=[f't{i}' for i in range(M.L)]
    cands,frac=hd.rank_integral_tree_candidates(M,res,topk=topk,root_beam=256 if M.L<=12 else 384,leaf_beam=64 if M.L<=12 else 96,child_branching=16,leaf_branching=10,merge_branching=128,max_expansions=500000)
    vals=[];trank=9999
    for i,z in enumerate(cands,1):
        if bench.key(z.topology)==bench.key(true): trank=i
        ll=TreeLikelihood(z.topology,labels,seqs,model,decoded=z,frac=frac,short_time=pg.base.branch_times,long_ratio=2.0).optimize(multistart=multistart,maxiter=120,ftol=1e-8).log_likelihood
        vals.append((ll,i,z.proportion))
    vals.sort(reverse=True)
    return (vals[0] if vals else (-math.inf,0,0)),len(cands),trank

def run_case(n=8,S=40,regime='moderate',rep=0,mode='spectral',delta=0.0,anchors=8,facets=4,size_weight='equal',topk=20,multistart=2,max_flow=None,sampled_points=None,physical_grid=None):
    seed=810000+n*1000+S*10+rep+(0 if regime=='moderate' else 100000);rng=random.Random(seed)
    sim=bench.gtr_f(bench.PI,bench.RATES);true=bench.random_tree(n,rng);seqs,cols,lengths=bench.simulate(true,S,sim,rng,regime);model=bench.empirical_model(cols);bench.configure(model);lib=bench.build_lib(n,'/mnt/data/spectral_bottleneck_libs')
    t=time.time();M=pg.PrunedGeneralModel(lib,cols,preprocess_mode='fast',substitution_model=model,branch_relaxation=mode,spectral_support_points=7); build=time.time()-t
    t=time.time();lp_score,res,solve= M.solve();
    basebest,nc0,tr0=rescore_candidates(M,res,seqs,model,true,topk,multistart)
    t=time.time();locobj,stats=add_local_fragment_scores(M,anchors,facets,size_weight,max_flow,sampled_points,physical_grid);couple_build=time.time()-t
    t=time.time();locval,res2,secsolve=solve_secondary(M,res,locobj,delta); secsolve_wall=time.time()-t
    best,nc,tr=rescore_candidates(M,res2,seqs,model,true,topk,multistart)
    true_ll=bench.rescore(true,[f't{i}' for i in range(n)],seqs,model)
    return dict(n=n,S=S,regime=regime,rep=rep,mode=mode,delta=delta,anchors=anchors,facets=facets,size_weight=size_weight,max_flow=max_flow,sampled_points=sampled_points,physical_grid=physical_grid,
                lp_score=lp_score,primary_build_s=build,primary_solve_s=solve,baseline_best=basebest[0],baseline_candidates=nc0,baseline_true_rank=tr0,
                local_score=locval,couple_build_s=couple_build,couple_solve_s=secsolve,couple_wall_s=secsolve_wall,best_logL=best[0],best_structural_rank=best[1],best_bottleneck=best[2],candidates=nc,true_rank=tr,true_logL=true_ll,best_minus_true=best[0]-true_ll,**stats)

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=8);ap.add_argument('--S',type=int,default=40);ap.add_argument('--regime',default='moderate');ap.add_argument('--rep',type=int,default=0);ap.add_argument('--mode',default='spectral');ap.add_argument('--delta',type=float,default=0);ap.add_argument('--anchors',type=int,default=8);ap.add_argument('--facets',type=int,default=4);ap.add_argument('--size-weight',default='equal');ap.add_argument('--physical-grid',type=int,default=None);ap.add_argument('--sampled-points',type=int,default=None);ap.add_argument('--max-flow',type=int,default=None);ap.add_argument('--topk',type=int,default=20);ap.add_argument('--multistart',type=int,default=2);ap.add_argument('--out',default='/mnt/data/local_fragment_coupling.csv');a=ap.parse_args()
 r=run_case(a.n,a.S,a.regime,a.rep,a.mode,a.delta,a.anchors,a.facets,a.size_weight,a.topk,a.multistart,a.max_flow,a.sampled_points,a.physical_grid);print(r,flush=True)
 p=Path(a.out);ex=p.exists()
 with p.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(r));
  if not ex:w.writeheader()
  w.writerow(r)
