import importlib.util,math,numpy as np
from collections import Counter
spec=importlib.util.spec_from_file_location('c','/mnt/data/multisite_compressed_onepass.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
h=c.h;base=c.base;BASES=c.BASES;BETA=c.BETA
COLS=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
cnt=Counter(COLS);P=list(cnt);W=np.array([cnt[p] for p in P],float)

def rootdata(pat,M):
 counts=tuple(pat.count(b) for b in BASES);out,_=h.build_composition_boxes_restricted(counts);vals={}
 for u in M.roots:vals[u]=h.root_split_upper(counts,out,int(M.units[u]['split_a']))[0]
 U=max(vals.values());L=max(min(base.root_split_lower(counts,out,int(M.units[u]['split_a'])) for u in M.roots),1e-300)
 return vals,U,L

def run(anchor):
 M=c.Model('/mnt/data/multisite_units8/library.json',[P[anchor]],1.1)
 # replace default anchor weight 1 with actual multiplicity
 M.obj={k:v*W[anchor] for k,v in M.obj.items()}; const=W[anchor]*(math.log(M.info[0]['sroot'])+8*math.log(BETA));err=W[anchor]*c.chord_gap(1.1)
 for p,pat in enumerate(P):
  if p==anchor:continue
  vals,U,L=rootdata(pat,M);q=M.lp.var(f'roq:{p}',L/U,1);M.lp.add_eq({q:1,**{M.act[u]:-vals[u]/U for u in M.roots}},0)
  pts=c.geometric_breaks(L/U,1.1);lam=[]
  for k,x in enumerate(pts):
   v=M.lp.var(f'rolam:{p}:{k}',0,1);lam.append(v);M.obj[v]=W[p]*math.log(x)
  M.lp.add_eq({v:1 for v in lam},1);M.lp.add_eq({q:1,**{lam[k]:-pts[k] for k in range(len(pts))}},0)
  const+=W[p]*(math.log(U)+8*math.log(BETA));err+=W[p]*c.chord_gap(1.1)
 val,res,dt=M.lp.solve(M.obj,True);return val+const,len(M.lp.names),len(M.lp.eq),len(M.lp.le),dt,err
if __name__=='__main__':
 exact,_,_=base.exact_shared_best(P,W);print('exact',exact)
 for i,p in enumerate(P):print(i,''.join(p),run(i))
