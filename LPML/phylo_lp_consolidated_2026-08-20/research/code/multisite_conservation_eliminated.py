#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,math,numpy as np
from collections import Counter,defaultdict
from pathlib import Path
import json
spec=importlib.util.spec_from_file_location('cp','/mnt/data/multisite_conservation_projection.py');cp=importlib.util.module_from_spec(spec);spec.loader.exec_module(cp)
ms=cp.ms;h=cp.h;base=cp.base;PI=h.PI;BETA=h.BETA;BASES=h.BASES

class EliminatedModel:
 def __init__(self,libpath,patterns,weights):
  self.lib=json.loads(Path(libpath).read_text());self.L=int(self.lib['summary']['leaf_count']);self.patterns=patterns;self.weights=np.array(weights,float);self.P=len(patterns);self.lp=base.SparseLP();self.units=self.lib['units'];self.roots=set(self.lib['root_candidates']);self.splits=[u for u,x in self.units.items() if x['kind']=='split'];self.tax=[u for u,x in self.units.items() if x['kind']=='taxon'];self.by=defaultdict(list)
  for u in self.splits:self.by[int(self.units[u]['descendant_count'])].append(u)
  self.reqdef=defaultdict(list)
  for u in self.splits:
   for port,f in [('left',int(self.units[u]['left_flow'])),('right',int(self.units[u]['right_flow']))]:
    if f<self.L:self.reqdef[f].append((u,port))
  self.act={u:self.lp.var('a:'+u,0,1) for u in self.splits};self.lp.add_eq({self.act[u]:1 for u in self.roots},1)
  # Keep the shared full topology transportation only once; site likelihood no longer has cell variables.
  self.outs={1:list(self.tax)}
  for m in range(2,self.L):self.outs[m]=list(self.by[m])
  self.y={}
  for m in range(1,self.L):
   R=self.reqdef[m];O=self.outs[m]
   for r,(pu,port) in enumerate(R):
    row={}
    for o,cu in enumerate(O):v=self.lp.var(f'y:{m}:{r}:{o}',0,1);self.y[m,r,o]=v;row[v]=1
    self.lp.add_eq({**row,self.act[pu]:-1},0)
   for o,cu in enumerate(O):
    col={self.y[m,r,o]:1 for r in range(len(R))}
    if m==1:self.lp.add_eq(col,1)
    else:self.lp.add_eq({**col,self.act[cu]:-1},0)
  self.info=[]
  for pat in patterns:
   counts=tuple(pat.count(b) for b in BASES);out,_=h.build_composition_boxes_restricted(counts);cb={};scale={}
   for m in range(1,self.L):
    Bs=[B for (mm,c),B in out.items() if mm==m];lo=np.min(np.stack([B[:,0] for B in Bs]),axis=0);hi=np.max(np.stack([B[:,1] for B in Bs]),axis=0);cb[m]=(lo,hi);scale[m]=float(hi.max())
   cbn={m:(cb[m][0]/scale[m],cb[m][1]/scale[m]) for m in cb};roots={u:(base.root_split_lower(counts,out,int(self.units[u]['split_a'])),h.root_split_upper(counts,out,int(self.units[u]['split_a']))[0]) for u in self.roots};rsc={u:scale[int(self.units[u]['left_flow'])]*scale[int(self.units[u]['right_flow'])] for u in self.roots};sroot=max(rsc.values());rootsn={u:(roots[u][0]/sroot,roots[u][1]/sroot) for u in self.roots}
   self.info.append({'cb':cbn,'s':scale,'roots':rootsn,'rootcoef':{u:rsc[u]/sroot for u in self.roots},'sroot':sroot})
  self.req={};self.z={};self.q={};self.logv={};self.build()
 def post_expr(self,p,u,i):
  # normalized post-edge X_i = coef * sum_j P_ij Z_j
  I=self.info[p];aa=int(self.units[u]['left_flow']);bb=int(self.units[u]['right_flow']);m=int(self.units[u]['descendant_count']);coef=I['s'][aa]*I['s'][bb]/I['s'][m];P=h.transition_for_count(m)
  return {self.z[p,u][j]:coef*float(P[i,j]) for j in range(4)}
 def build(self):
  # first create all request and product variables, because post-edge expressions reference z
  for p,pat in enumerate(self.patterns):
   I=self.info[p]
   for m in range(1,self.L):
    lo,hi=I['cb'][m]
    for r,(pu,port) in enumerate(self.reqdef[m]):
     a=self.act[pu];rv=[]
     for i in range(4):
      v=self.lp.var(f'req:{p}:{m}:{r}:{i}',0,float(hi[i]));rv.append(v);self.lp.add_ge({v:1,a:-float(lo[i])},0);self.lp.add_le({v:1,a:-float(hi[i])},0)
     self.req[p,pu,port]=rv
   for u in self.splits:
    aa=int(self.units[u]['left_flow']);bb=int(self.units[u]['right_flow']);a=self.act[u];A=self.req[p,u,'left'];B=self.req[p,u,'right'];lA,uA=I['cb'][aa];lB,uB=I['cb'][bb];zz=[]
    for i in range(4):
     z=self.lp.var(f'z:{p}:{u}:{i}',0,min(float(uA[i]),float(uB[i])));zz.append(z);base.perspective_mcc(self.lp,z,A[i],B[i],a,float(lA[i]),float(uA[i]),float(lB[i]),float(uB[i]));self.lp.add_le({z:1,A[i]:-1},0);self.lp.add_le({z:1,B[i]:-1},0)
    self.z[p,u]=zz
   # post-edge bounds and count conservation, with X substituted out
   for m in range(1,self.L):
    R=self.reqdef[m];lo,hi=I['cb'][m]
    if m>1:
     for u in self.outs[m]:
      a=self.act[u]
      for i in range(4):
       e=self.post_expr(p,u,i)
       self.lp.add_ge({**e,a:-float(lo[i])},0);self.lp.add_le({**e,a:-float(hi[i])},0)
    for i in range(4):
     d={self.req[p,pu,port][i]:1 for pu,port in R}
     if m==1:
      total=sum((h.scaled_leaf_vector(pat[int(self.units[uid]['label'].split('_')[-1])-1])/I['s'][1])[i] for uid in self.tax);self.lp.add_eq(d,total)
     else:
      for u in self.outs[m]:
       for v,c in self.post_expr(p,u,i).items():d[v]=d.get(v,0)-c
      self.lp.add_eq(d,0)
   q=self.lp.var(f'q:{p}',0,1);self.q[p]=q;eq={q:1}
   for u in self.roots:
    coef=I['rootcoef'][u];z=self.z[p,u];a=self.act[u]
    for i in range(4):eq[z[i]]=eq.get(z[i],0)-coef*float(PI[i])
    lo,hi=I['roots'][u];expr={z[i]:coef*float(PI[i]) for i in range(4)};self.lp.add_le({**expr,a:-hi},0);self.lp.add_ge({**expr,a:-lo},0)
   self.lp.add_eq(eq,0);t=self.lp.var(f'log:{p}',-1000,10);self.logv[p]=t
   for x0 in (1,.1,.01,.001,.0001):self.addtan(p,x0)
 def addtan(self,p,x0):self.lp.add_le({self.logv[p]:1,self.q[p]:-1/x0},math.log(x0)-1)
 def solve(self,tol=1e-4,maxit=12):
  obj={self.logv[p]:self.weights[p] for p in range(self.P)};tot=0;hist=[]
  for it in range(maxit):
   val,res,dt=self.lp.solve(obj,True);tot+=dt
   if not res.success:return None,res,tot,hist
   gaps=[]
   for p in range(self.P):q=max(res.x[self.q[p]],1e-12);gaps.append(res.x[self.logv[p]]-math.log(q))
   hist.append((it,max(gaps),dt))
   if max(gaps)<=tol:break
   for p,g in enumerate(gaps):
    if g>tol:self.addtan(p,max(res.x[self.q[p]],1e-12))
  const=sum(self.weights[p]*(math.log(self.info[p]['sroot'])+self.L*math.log(BETA)) for p in range(self.P));return val+const,res,tot,hist

def main():
 cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')];cnt=Counter(cols);pats=list(cnt);w=[cnt[p] for p in pats];exact,_,_=base.exact_shared_best(pats,np.array(w,float));M=EliminatedModel('/mnt/data/multisite_units8/library.json',pats,w);print('size',len(M.lp.names),len(M.lp.eq),len(M.lp.le));u,r,t,hist=M.solve();print('exact',exact,'ub',u,'gap',u-exact,'time',t,'hist',hist)
if __name__=='__main__':main()
