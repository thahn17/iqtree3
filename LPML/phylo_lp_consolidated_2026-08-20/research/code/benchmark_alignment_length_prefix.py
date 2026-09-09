#!/usr/bin/env python3
from __future__ import annotations
import random,time,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import benchmark_spectral_bottleneck as bench
import benchmark_large_unit_usefulness as lu

def run(n=30,Smax=30000,rep=0,ncand=180,prefixes=(1000,5000,10000,30000)):
    seed=7710000+n*10000+rep;rng=random.Random(seed);sim=bench.gtr_f(bench.PI,bench.RATES);true=bench.random_tree(n,rng);seqs,cols,_=bench.simulate(true,Smax,sim,rng,'moderate');C=lu.generate_candidates(true,ncand,rng);rows=[]
    for S in prefixes:
        seqp=[s[:S] for s in seqs];colp=cols[:S];model=bench.empirical_model(colp);X=lu.onehot_sequences(seqp);cache,norm,meta=lu.unit_scores(C,X,model,.1,'single');tk=bench.key(true);alls={}
        for T in C:
            sc=lu.tree_scores(T,cache,norm,meta,n);d=lu.rf(T,true);k=bench.key(T)
            for name,v in sc.items():alls.setdefault(name,[]).append((v,d,k==tk))
        for name,arr in alls.items():
            order=sorted(arr,key=lambda q:q[0],reverse=True);rank=next(i+1 for i,q in enumerate(order) if q[2]);rows.append(dict(n=n,S=S,rep=rep,ncand=len(C),metric=name,true_rank=rank,top_rf=order[0][1]))
    return rows
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=30);p.add_argument('--Smax',type=int,default=30000);p.add_argument('--rep',type=int,default=0);p.add_argument('--ncand',type=int,default=180);p.add_argument('--out',default='/mnt/data/alignment_length_prefix.csv');a=p.parse_args();rows=run(a.n,a.Smax,a.rep,a.ncand);print(rows,flush=True);f=Path(a.out);ex=f.exists();
 with f.open('a',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]));
  if not ex:w.writeheader()
  w.writerows(rows)
