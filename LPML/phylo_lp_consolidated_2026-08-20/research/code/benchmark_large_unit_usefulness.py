#!/usr/bin/env python3
from __future__ import annotations
import random,math,time,csv,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.linalg import expm
sys.path.insert(0,str(Path(__file__).resolve().parent))
import benchmark_spectral_bottleneck as bench

BASE='ACGT'; BI={b:i for i,b in enumerate(BASE)}

def set_sub(t,path,val):
    if not path:return val
    i=path[0];q=list(t);q[i]=set_sub(q[i],path[1:],val);return bench.canon(tuple(q))

def nni_neighbors(t):
    out={}
    def rec(node,path=()):
        if isinstance(node,int):return
        for ci in (0,1):
            C=node[ci];S=node[1-ci]
            if not isinstance(C,int):
                for cj in (0,1):
                    A=C[cj];B=C[1-cj];C2=bench.canon((S,B));P2=bench.canon((A,C2));z=bench.canon(set_sub(t,path,P2));out[bench.key(z)]=z
            rec(C,path+(ci,))
    rec(t);out.pop(bench.key(t),None);return list(out.values())

def walk(t,steps,rng):
    cur=t
    for _ in range(steps):
        ns=nni_neighbors(cur)
        if not ns:break
        cur=rng.choice(ns)
    return cur

def clades(t):
    out=[]
    def rec(z):
        if isinstance(z,int):return frozenset([z])
        A=rec(z[0]);B=rec(z[1]);C=A|B
        if len(C)>1:out.append(C)
        return C
    root=rec(t);return set(c for c in out if c!=root)

def split_units(t):
    out=[]
    def rec(z):
        if isinstance(z,int):return frozenset([z])
        A=rec(z[0]);B=rec(z[1]);C=A|B
        aa,bb=(A,B) if (len(A),tuple(sorted(A))) <= (len(B),tuple(sorted(B))) else (B,A)
        out.append((C,aa,bb));return C
    rec(t);return out

def rf(t,true):
    A=clades(t);B=clades(true);return len(A^B)

def generate_candidates(true,ncand,rng):
    n=len(bench.taxa(true)); d={bench.key(true):true}
    # dense local neighborhoods plus farther walks and a few random trees
    frontier=[true]
    for _ in range(2):
        nf=[]
        for z in frontier:
            for q in nni_neighbors(z):
                if len(d)>=ncand//2:break
                k=bench.key(q)
                if k not in d:d[k]=q;nf.append(q)
            if len(d)>=ncand//2:break
        frontier=nf or frontier
    while len(d)<int(ncand*.9):
        st=rng.randint(2,max(3,n//3));q=walk(true,st,rng);d.setdefault(bench.key(q),q)
    while len(d)<ncand:
        q=bench.random_tree(n,rng);d.setdefault(bench.key(q),q)
    return list(d.values())[:ncand]

def onehot_sequences(seqs):
    n=len(seqs);S=len(seqs[0]);X=np.zeros((n,S,4),float)
    for t,s in enumerate(seqs):
        for j,ch in enumerate(s):X[t,j,BI[ch]]=1
    return X

def unit_scores(candidates,X,model,tref=0.1,kind='single'):
    n,S,_=X.shape
    if kind=='single':P=expm(model.q*tref)
    elif kind=='logmean':
        ts=np.geomspace(max(tref/2,1e-5),tref*2,7);P=np.mean(np.stack([expm(model.q*t) for t in ts]),axis=0)
    elif kind=='wide':
        ts=np.r_[0,np.geomspace(.01,1,7),5];P=np.mean(np.stack([expm(model.q*t) for t in ts]),axis=0)
    else:raise ValueError(kind)
    total=X.sum(0)/n; MG=total@P.T; baseline=(MG*model.pi[None,:]).sum(1)
    cache={}; meta={}
    all_units={}
    for T in candidates:
        for C,A,B in split_units(T):all_units[(A,B)]=(len(C),len(A),len(B))
    for (A,B),(m,a,b) in all_units.items():
        fA=X[list(A)].mean(0);fB=X[list(B)].mean(0);MA=fA@P.T;MB=fB@P.T
        ov=(np.minimum(MA,MB)*model.pi[None,:]).sum(1)
        raw=float(np.mean(ov));center=float(np.mean(ov-baseline));cache[(A,B)]=(raw,center);meta[(A,B)]=(m,a,b)
    # flow-size normalization based on all candidate units (diagnostic; identifies size drift).
    bym=defaultdict(list)
    for k,(raw,cen) in cache.items():bym[meta[k][0]].append(cen)
    norm={}
    for m,v in bym.items():
        arr=np.asarray(v);med=float(np.median(arr));mad=float(np.median(np.abs(arr-med)));scale=max(1.4826*mad,float(np.std(arr)),1e-8)
        for k in [q for q in cache if meta[q][0]==m]:norm[k]=(cache[k][1]-med)/scale
    return cache,norm,meta

def tree_scores(T,cache,norm,meta,n):
    units=split_units(T); vals=np.array([cache[(A,B)][1] for C,A,B in units],float); z=np.array([norm[(A,B)] for C,A,B in units],float)
    ms=np.array([len(C) for C,A,B in units],float); ab=np.array([4*len(A)*len(B)/(len(C)**2) for C,A,B in units],float)
    def wm(w):return float(np.sum(w*vals)/max(np.sum(np.abs(w)),1e-12))
    # Shift z for bottleneck only; higher is better.
    return {
      'sum_equal':float(vals.mean()),
      'sum_invflow':wm(1/ms),
      'sum_sqrtinv':wm(1/np.sqrt(ms)),
      'sum_balance':wm(ab),
      'sum_balance_sqrt':wm(np.sqrt(ab)),
      'sum_bal_sqrtinv':wm(np.sqrt(ab)/np.sqrt(ms)),
      'sum_bal_invflow':wm(ab/np.sqrt(ms)),
      'flow_z_mean':float(z.mean()),
      'bottleneck_raw':float(vals.min()),
      'bottleneck_z':float(z.min()),
      'q10_z':float(np.quantile(z,.10)),
      'q20_z':float(np.quantile(z,.20)),
    }

def run(n=20,S=1000,rep=0,ncand=300,tref=.1,kind='single'):
    seed=6010000+n*10000+S*10+rep;rng=random.Random(seed);sim=bench.gtr_f(bench.PI,bench.RATES);true=bench.random_tree(n,rng);seqs,cols,_=bench.simulate(true,S,sim,rng,'moderate');model=bench.empirical_model(cols)
    t=time.time();C=generate_candidates(true,ncand,rng);gen=time.time()-t;X=onehot_sequences(seqs);t=time.time();cache,norm,meta=unit_scores(C,X,model,tref,kind);scoretime=time.time()-t
    metrics=list(tree_scores(true,cache,norm,meta,n));rows=[];truekey=bench.key(true)
    allscores={k:[] for k in metrics}
    for T in C:
        sc=tree_scores(T,cache,norm,meta,n);d=rf(T,true)
        for k,v in sc.items():allscores[k].append((v,d,bench.key(T)==truekey))
    for k,arr in allscores.items():
        order=sorted(arr,key=lambda q:q[0],reverse=True);rank=next(i+1 for i,q in enumerate(order) if q[2]);top=order[0]
        vals=np.array([q[0] for q in arr]);dists=np.array([q[1] for q in arr],float)
        corr=float(np.corrcoef(vals,-dists)[0,1]) if np.std(vals)>0 and np.std(dists)>0 else 0
        rows.append(dict(n=n,S=S,rep=rep,ncand=len(C),tref=tref,kind=kind,metric=k,true_rank=rank,top_rf=top[1],corr_neg_rf=corr,candidate_gen_s=gen,score_s=scoretime,unique_units=len(cache)))
    return rows

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=20);ap.add_argument('--S',type=int,default=1000);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--ncand',type=int,default=300);ap.add_argument('--tref',type=float,default=.1);ap.add_argument('--kind',default='single');ap.add_argument('--out',default='/mnt/data/large_unit_usefulness.csv');a=ap.parse_args();rows=run(a.n,a.S,a.rep,a.ncand,a.tref,a.kind);print(rows,flush=True);p=Path(a.out);ex=p.exists()
 with p.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));
  if not ex:w.writeheader()
  w.writerows(rows)
