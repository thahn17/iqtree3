#!/usr/bin/env python3
from __future__ import annotations
import random,math,time,csv,sys
from pathlib import Path
import numpy as np
sys.path[:0]=['/mnt/data/work_v3','/mnt/data/_v3','/mnt/data']
import benchmark_spectral_bottleneck as bench
import benchmark_large_unit_usefulness as b

def collect_moves(t,path=()):
    out=[]
    if isinstance(t,int): return out
    for ci in (0,1):
        C=t[ci]
        if not isinstance(C,int):
            out.append((path,ci,0)); out.append((path,ci,1))
            out.extend(collect_moves(C,path+(ci,)))
    return out

def apply_move(t,move):
    path,ci,cj=move
    def rec(node,p=()):
        if p==path:
            C=node[ci]; S=node[1-ci]; A=C[cj]; B=C[1-cj]
            C2=bench.canon((S,B)); P2=bench.canon((A,C2)); return P2
        if isinstance(node,int): return node
        return bench.canon((rec(node[0],p+(0,)),rec(node[1],p+(1,))))
    return bench.canon(rec(t))

def random_nni(t,rng,steps):
    cur=t
    for _ in range(steps):
        mv=collect_moves(cur)
        if not mv: break
        cur=apply_move(cur,rng.choice(mv))
    return cur

def candidates(true,ncand,rng):
    d={bench.key(true):true}
    depths=[1,2,3,5,8,12,20,30]
    i=0
    while len(d)<int(.9*ncand):
        q=random_nni(true,rng,depths[i%len(depths)]); d.setdefault(bench.key(q),q); i+=1
        if i>ncand*30: break
    while len(d)<ncand:
        q=bench.random_tree(len(bench.taxa(true)),rng); d.setdefault(bench.key(q),q)
    return list(d.values())[:ncand]

def pair_dist(seqs):
    n=len(seqs); S=len(seqs[0]); X=np.frombuffer(''.join(seqs).encode(),dtype='S1').reshape(n,S); D=np.zeros((n,n),np.float32)
    for i in range(n-1):
        v=np.mean(X[i+1:]!=X[i],axis=1,dtype=np.float64).astype(np.float32)
        D[i,i+1:]=v; D[i+1:,i]=v
    return D

def mp(D,A,B=None):
    A=np.fromiter(A,dtype=np.int32)
    if B is None:
        if len(A)<2:return 0.0
        Q=D[np.ix_(A,A)]; return float(Q.sum()/(len(A)*(len(A)-1)))
    B=np.fromiter(B,dtype=np.int32)
    return float(D[np.ix_(A,B)].mean())

def build_unit_scores(C,D):
    dat={}
    for T in C:
        for Cc,A,B in b.split_units(T):
            k=(A,B)
            if k in dat:continue
            cross=mp(D,A,B); wa=mp(D,A); wb=mp(D,B); s=cross-max(wa,wb)
            m=len(Cc); a=len(A); bb=len(B); bal=4*a*bb/(m*m)
            dat[k]=(s,m,bal)
    vals=np.array([v[0] for v in dat.values()]); ref=float(np.quantile(vals,.9))
    return dat,ref

def rf(T,true):return b.rf(T,true)

def score(T,dat,ref,beta,gamma):
    U=b.split_units(T); s=np.array([dat[(A,B)][0] for C,A,B in U]); m=np.array([dat[(A,B)][1] for C,A,B in U],float); B=np.array([dat[(A,B)][2] for C,A,B in U],float)
    if beta==0 and gamma==0:return float(s.min())
    d=ref-s; w=(m**beta)*(np.clip(B,1e-6,None)**gamma); w/=math.exp(float(np.mean(np.log(np.clip(w,1e-12,None)))))
    return -float(np.max(w*d))

def run(n=100,S=10000,rep=0,ncand=180,regime='moderate'):
    seed=13370000+n*10000+S+rep; rng=random.Random(seed); sim=bench.gtr_f(bench.PI,bench.RATES); true=bench.random_tree(n,rng)
    t=time.time(); seqs,cols,_=bench.simulate(true,S,sim,rng,regime); tsim=time.time()-t
    t=time.time(); C=candidates(true,ncand,rng); tc=time.time()-t
    t=time.time(); D=pair_dist(seqs); td=time.time()-t
    t=time.time(); dat,ref=build_unit_scores(C,D); tu=time.time()-t
    tk=bench.key(true); rows=[]
    grid=[(0,0),(-.05,-.05),(-.1,-.1),(-.15,-.15),(-.2,-.2),(-.3,-.3),(-.2,0),(0,-.2),(.1,0),(0,.1)]
    for beta,gamma in grid:
        arr=[(score(T,dat,ref,beta,gamma),rf(T,true),bench.key(T)==tk) for T in C]; arr.sort(reverse=True,key=lambda q:q[0]); rank=next(i+1 for i,q in enumerate(arr) if q[2]); rows.append(dict(n=n,S=S,rep=rep,regime=regime,ncand=len(C),beta=beta,gamma=gamma,true_rank=rank,top_rf=arr[0][1],sim_s=tsim,candidate_s=tc,pairdist_s=td,unit_score_s=tu,unique_units=len(dat),ref=ref))
    return rows
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=100);ap.add_argument('--S',type=int,default=10000);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--ncand',type=int,default=180);ap.add_argument('--regime',default='moderate');ap.add_argument('--out',default='/mnt/data/unit_weight_100taxa.csv');a=ap.parse_args();rows=run(a.n,a.S,a.rep,a.ncand,a.regime);print(rows,flush=True);p=Path(a.out);ex=p.exists()
 with p.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));
  if not ex: w.writeheader()
  w.writerows(rows)
