#!/usr/bin/env python3
exec(open('/mnt/data/_v3/benchmark_nearopt_face2.py').read().split("if __name__=='__main__':")[0])
import copy
n,S,seed=8,60,18008;rng=random.Random(seed);simmod=gtr_f(PI,RATES);true=random_tree(n,rng);seq,cols=sim(true,S,simmod,rng);cnt=np.array([sum(s.count(b) for s in seq) for b in BASE]);model=gtr_f((cnt+.5)/(cnt.sum()+2),RATES);cfg(model);M=pg.PrunedGeneralModel(lib(n),cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='spectral',spectral_support_points=5);raw,res,_=M.lp.solve(M.obj,True)
rows=[]
for k in sorted(M.by):
 if k==n:continue
 for u in M.by[k]:
  L=copy.deepcopy(M.lp);L.add_eq({M.act[u]:1},1)
  try:v,rr,dt=L.solve(M.obj,True);drop=raw-v
  except Exception:drop=float('inf')
  rows.append((k,u,float(res.x[M.act[u]]),drop))
for k in sorted(set(r[0] for r in rows)):
 z=[r for r in rows if r[0]==k];ds=np.array([r[3] for r in z]);print('flow',k,'units',len(z),'active_now',sum(r[2]>1e-8 for r in z),'drop min/median/max',np.min(ds),np.median(ds),np.max(ds),'zero<=1e-7',sum(ds<=1e-7))
print('all',len(rows),'zero',sum(r[3]<=1e-7 for r in rows))
