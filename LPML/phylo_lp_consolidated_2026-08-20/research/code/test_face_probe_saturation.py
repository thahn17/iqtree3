#!/usr/bin/env python3
exec(open('/mnt/data/_v3/benchmark_nearopt_face2.py').read().split("if __name__=='__main__':")[0])
from types import SimpleNamespace
n,S,seed=8,60,15008;rng=random.Random(seed);simmod=gtr_f(PI,RATES);true=random_tree(n,rng);seq,cols=sim(true,S,simmod,rng);cnt=np.array([sum(s.count(b) for s in seq) for b in BASE]);model=gtr_f((cnt+.5)/(cnt.sum()+2),RATES);cfg(model);M=pg.PrunedGeneralModel(lib(n),cols,preprocess_mode='fast',substitution_model=model,branch_relaxation='spectral',spectral_support_points=5);raw,res,_=M.lp.solve(M.obj,True);X,pt=sketch_points(M,raw,res,4.0,32,seed+77);truell=ll(true,seq,model);print('true',truell,'split units',len(M.splits),'root profiles',len(M.roots),'probe_time',pt)
base_roots=1+len(M.roots)
for P in [0,2,4,8,12,16,24,32]:
 m=1 if P==0 else min(len(X),base_roots+P)
 Z=X[:m];x=np.quantile(Z,.9,axis=0) if P else Z[0]
 c,_=rank_integral_tree_candidates(M,SimpleNamespace(x=x),topk=80,root_beam=320,leaf_beam=80,max_expansions=280000);best,_=cand_ll(c,seq,model,8)
 print('P',P,'points',m,'bestll',best,'n',len(c))
