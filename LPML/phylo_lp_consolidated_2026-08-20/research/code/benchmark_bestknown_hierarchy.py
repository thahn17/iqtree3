#!/usr/bin/env python3
import sys,random,math,time
sys.path.insert(0,'/mnt/data')
import numpy as np,pandas as pd
from collections import defaultdict
import benchmark_similarity_hierarchy as b
import benchmark_fractional_composition_grammar as g

def nni_neighbors(t):
 out={}
 def rec(node,path=()):
  if isinstance(node,int):return
  for ci in (0,1):
   C=node[ci];S=node[1-ci]
   if not isinstance(C,int):
    for cj in (0,1):
     A=C[cj];B=C[1-cj];C2=b.canon((S,B));P2=b.canon((A,C2));z=b.canon(b.set_sub(t,path,P2));out[b.key(z)]=z
   rec(C,path+(ci,))
 rec(t);out.pop(b.key(t),None);return list(out.values())

def search_alignment(true,seq,rng,restarts=4,maxit=10):
 cache={}
 def score(t):
  k=b.key(t)
  if k not in cache:cache[k]=(t,b.loglik(t,seq))
  return cache[k][1]
 starts=[true]
 for r in range(restarts-1):starts.append(b.mutate(true,1+r*2,rng) if r<3 else b.random_tree(len(b.taxa(true)),rng))
 local=[]
 for st in starts:
  cur=st;cl=score(cur)
  for _ in range(maxit):
   ns=nni_neighbors(cur);vals=[score(z) for z in ns]
   if not vals:break
   j=int(np.argmax(vals))
   if vals[j]<=cl+1e-9:break
   cur,cl=ns[j],vals[j]
  local.append((cl,cur))
 bestll,best=max(local,key=lambda q:q[0])
 vals=sorted(cache.values(),key=lambda q:q[1],reverse=True)
 return best,bestll,vals

def make_ensemble(visited,eps,rng,topn=32,noise_n=32,temp=2.0):
 top=visited[:min(topn,len(visited))];rest=visited[min(topn,len(visited)):]
 x=np.array([ll for _,ll in top]);w=np.exp(np.clip((x-x.max())/temp,-50,0));w/=w.sum();w=(1-eps)*w
 noise=rng.sample(rest,min(noise_n,len(rest))) if rest else []
 trees=[t for t,_ in top]+[t for t,_ in noise];weights=list(w)+([eps/len(noise)]*len(noise) if noise else [])
 if not noise:weights=list(np.array(weights)/sum(weights))
 return trees,weights

def supports(trees,weights):
 cs=defaultdict(float);ss=defaultdict(float);pair=defaultdict(float)
 for T,w in zip(trees,weights):
  cl,sp=b.clades_and_splits(T)
  for C in cl:cs[C]+=w
  for C,A,B in sp:
   aa,bb=(A,B) if tuple(sorted(A))<=tuple(sorted(B)) else (B,A);ss[(C,aa,bb)]+=w
 return cs,ss,pair

def run(n,S,rep,eps,seed):
 rng=random.Random(seed);true=b.random_tree(n,rng);seq=b.simulate_alignment(true,S,rng)
 t0=time.time();best,bestll,visited=search_alignment(true,seq,rng,restarts=4 if n<100 else 3,maxit=8 if n<100 else 6);search_time=time.time()-t0
 trees,w=make_ensemble(visited,eps,rng);root,g0,g1,g2=g.grammar_counts(trees,w);cs,ss,pair=supports(trees,w)
 # G2 candidate generation wide enough for difficult cases, then robust bottleneck rerank.
 tk=g.topk_g2(n,root,g2,K=500,beam=2048);cand=[t for t,_ in tk]
 ordered=sorted(cand,key=lambda t:(b.tree_metrics(t,cs,ss,pair,n)['bottleneck'],b.tree_metrics(t,cs,ss,pair,n)['split_log']),reverse=True)
 keys=[b.key(t) for t in ordered];rank=keys.index(b.key(best))+1 if b.key(best) in keys else 9999
 # composition leakage from fixed local fragments
 unique={b.key(t):t for t in trees};orig=sum(g.prob_tree(t,root,g0,g1,g2,2) for t in unique.values())
 return dict(n=n,S=S,rep=rep,eps=eps,rank=rank,top1=rank<=1,top3=rank<=3,top5=rank<=5,top10=rank<=10,top20=rank<=20,top50=rank<=50,top200=rank<=200,best_is_true=int(b.key(best)==b.key(true)),delta_vs_true=bestll-b.loglik(true,seq),visited=len(visited),ensemble_eff=1/sum(x*x for x in w),original_mass=min(orig,1),hybrid_leakage=max(0,1-orig),search_time=search_time)

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int);ap.add_argument('--S',type=int);ap.add_argument('--rep',type=int,default=0);ap.add_argument('--eps',type=float,default=.15);ap.add_argument('--out',default='/mnt/data/bestknown_hierarchy.csv');a=ap.parse_args();seed=2000000+int(a.eps*1000)*100000+a.n*10000+a.S*10+a.rep
 r=run(a.n,a.S,a.rep,a.eps,seed);print(r);import os
 if os.path.exists(a.out):df=pd.read_csv(a.out);df=pd.concat([df,pd.DataFrame([r])],ignore_index=True)
 else:df=pd.DataFrame([r])
 df.to_csv(a.out,index=False)
