#!/usr/bin/env python3
exec(open('/mnt/data/_v3/benchmark_nearopt_face2.py').read().split("if __name__=='__main__':")[0])
import copy
from types import SimpleNamespace

def center_solve(M,raw,delta,root_tmem=False):
 L=copy.deepcopy(M.lp);L.add_ge(M.obj,raw-delta);obj={}
 for k,us in M.by.items():
  wt=1.0/max(1,len(us))
  for u in us:
   s=L.var(f'cact:{u}',0,.5);L.add_le({s:1,M.act[u]:-1},0);L.add_le({s:1,M.act[u]:1},1);obj[s]=wt
 if root_tmem:
  # s <= x and s <= act-x for root-side memberships only; normalize by roots*taxa.
  wt=.5/max(1,len(M.roots)*M.L)
  for r in M.roots:
   A=M.act[r]
   for port in ('left','right'):
    mm,slot,_=M._req_slot(r,port)
    for tax in range(M.L):
     x=M.tmem[mm,0,slot,tax];s=L.var(f'ct:{r}:{port}:{tax}',0,.5);L.add_le({s:1,x:-1},0);L.add_le({s:1,x:1,A:-1},0);obj[s]=wt
 v,res,dt=L.solve(obj,True);return res,dt

def run(n,S,seed):
 rng=random.Random(seed);simmod=gtr_f(PI,RATES);true=random_tree(n,rng);seq,cols=sim(true,S,simmod,rng);cnt=np.array([sum(s.count(b) for s in seq) for b in BASE]);model=gtr_f((cnt+.5)/(cnt.sum()+2),RATES);cfg(model);M=pg.PrunedGeneralModel(lib(n),cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='spectral',spectral_support_points=5);raw,res,sv=M.lp.solve(M.obj,True);truell=ll(true,seq,model)
 print('CASE',n,S,'true',truell,'splits',len(M.splits))
 for d in [1,2,4,8]:
  for rt in [False,True]:
   rr,dt=center_solve(M,raw,d,rt);c,_=rank_integral_tree_candidates(M,rr,topk=50,root_beam=256,leaf_beam=64,max_expansions=180000);best,_=cand_ll(c,seq,model,8);print('delta',d,'rt',rt,'ll',best,'solve',dt,'n',len(c))
run(8,60,20008)
run(12,20,20012)
