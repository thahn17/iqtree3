#!/usr/bin/env python3
from __future__ import annotations
import random,math,time,csv,sys
from pathlib import Path
import numpy as np
sys.path[:0]=['/mnt/data/work_v3','/mnt/data/_v3','/mnt/data']
import benchmark_spectral_bottleneck as bench
import benchmark_large_unit_usefulness as b
from benchmark_100taxa_weighting import candidates
BASE='ACGT'; BI={c:i for i,c in enumerate(BASE)}

def cov_matrix(seqs,pi):
    n=len(seqs);S=len(seqs[0]); enc=np.fromiter((BI[c] for s in seqs for c in s),dtype=np.int8,count=n*S).reshape(n,S)
    K=np.zeros((n,n),np.float32); inv=1/np.asarray(pi,float)
    for i in range(n-1):
        a=enc[i]; rest=enc[i+1:]
        same=(rest==a)
        # baseline -1 each site; same base b adds 1/pi_b.
        vals=-np.ones((len(rest),S),dtype=np.float32)
        # avoid large 3d tensors: fill same positions using broadcast inv[a]
        vals[same]=np.broadcast_to(inv[a],vals.shape)[same]-1.0
        v=vals.mean(axis=1,dtype=np.float64).astype(np.float32);K[i,i+1:]=v;K[i+1:,i]=v
    np.fill_diagonal(K,3.0)
    return K

def mean_cross(K,A,B):
    A=np.fromiter(A,dtype=np.int32);B=np.fromiter(B,dtype=np.int32);return float(K[np.ix_(A,B)].mean())
def mean_within(K,A):
    A=np.fromiter(A,dtype=np.int32)
    if len(A)<2:return 3.0
    Q=K[np.ix_(A,A)];return float((Q.sum()-np.trace(Q))/(len(A)*(len(A)-1)))
def build(C,K):
    dat={}
    for T in C:
        for Cc,A,B in b.split_units(T):
            k=(A,B)
            if k in dat:continue
            wa=mean_within(K,A);wb=mean_within(K,B);cross=mean_cross(K,A,B)
            # both children should be internally correlated relative to between-child correlation.
            s=min(wa,wb)-cross
            m=len(Cc);aa=len(A);bb=len(B);Bf=4*aa*bb/(m*m)
            dat[k]=(s,m,Bf)
    vals=np.array([x[0] for x in dat.values()]); ref=float(np.quantile(vals,.9));return dat,ref

def score(T,dat,ref,beta,gamma):
    U=b.split_units(T);s=np.array([dat[(A,B)][0] for C,A,B in U]);m=np.array([dat[(A,B)][1] for C,A,B in U],float);B=np.array([dat[(A,B)][2] for C,A,B in U],float)
    if beta==0 and gamma==0:return float(s.min())
    d=ref-s;w=(m**beta)*(np.clip(B,1e-6,None)**gamma);w/=math.exp(float(np.mean(np.log(np.clip(w,1e-12,None)))));return -float(np.max(w*d))
def run(n=100,S=10000,rep=0,ncand=180,regime='moderate'):
    seed=13370000+n*10000+S+rep;rng=random.Random(seed);sim=bench.gtr_f(bench.PI,bench.RATES);true=bench.random_tree(n,rng);t=time.time();seqs,cols,_=bench.simulate(true,S,sim,rng,regime);ts=time.time()-t;model=bench.empirical_model(cols);t=time.time();C=candidates(true,ncand,rng);tc=time.time()-t;t=time.time();K=cov_matrix(seqs,model.pi);tk=time.time()-t;t=time.time();dat,ref=build(C,K);tu=time.time()-t;truek=bench.key(true);rows=[]
    for beta,gamma in [(0,0),(-.05,-.05),(-.1,-.1),(-.15,-.15),(-.2,-.2),(-.3,-.3),(-.1,0),(0,-.1)]:
        arr=[(score(T,dat,ref,beta,gamma),b.rf(T,true),bench.key(T)==truek) for T in C];arr.sort(reverse=True,key=lambda q:q[0]);rank=next(i+1 for i,q in enumerate(arr) if q[2]);rows.append(dict(n=n,S=S,rep=rep,regime=regime,ncand=len(C),beta=beta,gamma=gamma,true_rank=rank,top_rf=arr[0][1],sim_s=ts,cand_s=tc,cov_s=tk,units_s=tu,unique_units=len(dat),ref=ref))
    return rows
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=100);ap.add_argument('--S',type=int,default=10000);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--ncand',type=int,default=180);ap.add_argument('--regime',default='moderate');ap.add_argument('--out',default='/mnt/data/paircov_100taxa.csv');a=ap.parse_args();rows=run(a.n,a.S,a.rep,a.ncand,a.regime);print(rows,flush=True);p=Path(a.out);ex=p.exists()
 with p.open('a',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));
  if not ex:w.writeheader()
  w.writerows(rows)
