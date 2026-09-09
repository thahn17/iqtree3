#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,json,math,time
from collections import defaultdict
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('base','/mnt/data/multisite_shared_tree_lp_v3.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
h=base.h;BASES=h.BASES;PI=h.PI;BETA=h.BETA
scaled_leaf_vector=h.scaled_leaf_vector;transition_for_count=h.transition_for_count
build_composition_boxes_restricted=h.build_composition_boxes_restricted
root_split_upper=h.root_split_upper;root_split_lower=base.root_split_lower
SparseLP=base.SparseLP;perspective_mcc=base.perspective_mcc

def kv(vs):
 out=np.array([1.0])
 for v in vs:out=np.kron(out,v)
 return out

def km(ms):
 out=np.array([[1.0]])
 for M in ms:out=np.kron(out,M)
 return out

class Model:
 def __init__(self,libpath,blocks):
  self.lib=json.loads(Path(libpath).read_text());self.L=int(self.lib['summary']['leaf_count']);self.blocks=blocks;self.lp=SparseLP();self.units=self.lib['units'];self.roots=set(self.lib['root_candidates'])
  self.splits=[u for u,x in self.units.items() if x['kind']=='split'];self.tax=[u for u,x in self.units.items() if x['kind']=='taxon'];self.by=defaultdict(list)
  for u in self.splits:self.by[int(self.units[u]['descendant_count'])].append(u)
  self.reqdef=defaultdict(list)
  for u in self.splits:
   for port,f in [('left',int(self.units[u]['left_flow'])),('right',int(self.units[u]['right_flow']))]:
    if f<self.L:self.reqdef[f].append((u,port))
  self.outs={1:list(self.tax)}
  for m in range(2,self.L):self.outs[m]=list(self.by[m])
  self.act={u:self.lp.var('a:'+u,0,1) for u in self.splits};self.lp.add_eq({self.act[u]:1 for u in self.roots},1);self.y={}
  for m in range(1,self.L):
   R=self.reqdef[m];O=self.outs[m]
   for r,(pu,port) in enumerate(R):
    row={}
    for o,cu in enumerate(O):
     v=self.lp.var(f'y:{m}:{r}:{o}',0,1);self.y[m,r,o]=v;row[v]=1
    self.lp.add_eq({**row,self.act[pu]:-1},0)
   for o,cu in enumerate(O):
    col={self.y[m,r,o]:1 for r in range(len(R))}
    if m==1:self.lp.add_eq(col,1)
    else:self.lp.add_eq({**col,self.act[cu]:-1},0)
  self.info=[];self.x={};self.w={};self.req={};self.z={};self.q={};self.logv={};self.prep();self.build()
 def prep(self):
  for blk in self.blocks:
   sinfo=[]
   for pat in blk:
    counts=tuple(pat.count(b) for b in BASES);out,_=build_composition_boxes_restricted(counts);cb={};scale={}
    for m in range(1,self.L):
     Bs=[B for (mm,c),B in out.items() if mm==m];lo=np.min(np.stack([B[:,0] for B in Bs]),axis=0);hi=np.max(np.stack([B[:,1] for B in Bs]),axis=0);cb[m]=(lo,hi);scale[m]=float(hi.max())
    roots={u:(root_split_lower(counts,out,int(self.units[u]['split_a'])),root_split_upper(counts,out,int(self.units[u]['split_a']))[0]) for u in self.roots}
    sinfo.append({'cb':cb,'s':scale,'roots':roots})
   # block scale is product of marginal scales; normalized block component bounds are product bounds divided by this scale
   sblk={m:float(np.prod([s['s'][m] for s in sinfo])) for m in range(1,self.L)};cbblk={}
   for m in range(1,self.L):
    lo=kv([s['cb'][m][0] for s in sinfo])/sblk[m];hi=kv([s['cb'][m][1] for s in sinfo])/sblk[m];cbblk[m]=(lo,hi)
   root_scales={u: sblk[int(self.units[u]['left_flow'])]*sblk[int(self.units[u]['right_flow'])] for u in self.roots};sroot=max(root_scales.values())
   roots={u:(float(np.prod([s['roots'][u][0] for s in sinfo]))/sroot,float(np.prod([s['roots'][u][1] for s in sinfo]))/sroot) for u in self.roots}
   self.info.append({'sinfo':sinfo,'s':sblk,'cb':cbblk,'roots':roots,'rootcoef':{u:root_scales[u]/sroot for u in self.roots},'sroot':sroot,'pi':kv([PI]*len(blk)),'d':4**len(blk)})
 def build(self):
  for b,blk in enumerate(self.blocks):
   I=self.info[b];d=I['d']
   for u in self.splits:
    if u in self.roots:continue
    m=int(self.units[u]['descendant_count']);lo,hi=I['cb'][m];a=self.act[u];vv=[]
    for i in range(d):
     v=self.lp.var(f'x:{b}:{u}:{i}',0,float(hi[i]));vv.append(v);self.lp.add_ge({v:1,a:-float(lo[i])},0);self.lp.add_le({v:1,a:-float(hi[i])},0)
    self.x[b,u]=vv
   for m in range(1,self.L):
    R=self.reqdef[m];O=self.outs[m];lo,hi=I['cb'][m]
    for r,(pu,port) in enumerate(R):
     for o,cu in enumerate(O):
      y=self.y[m,r,o];wv=[]
      if m==1:
       taxidx=int(self.units[cu]['label'].split('_')[-1])-1;lv=kv([scaled_leaf_vector(pat[taxidx]) for pat in blk])/I['s'][1]
       for i in range(d):
        w=self.lp.var(f'w:{b}:{m}:{r}:{o}:{i}',0,float(lv[i]));wv.append(w);self.lp.add_eq({w:1,y:-float(lv[i])},0)
      else:
       for i in range(d):
        w=self.lp.var(f'w:{b}:{m}:{r}:{o}:{i}',0,float(hi[i]));wv.append(w);self.lp.add_ge({w:1,y:-float(lo[i])},0);self.lp.add_le({w:1,y:-float(hi[i])},0)
      self.w[b,m,r,o]=wv
    if m>1:
     for o,cu in enumerate(O):
      xv=self.x[b,cu]
      for i in range(d):self.lp.add_eq({**{self.w[b,m,r,o][i]:1 for r in range(len(R))},xv[i]:-1},0)
    for r,(pu,port) in enumerate(R):
     a=self.act[pu];rv=[]
     for i in range(d):
      v=self.lp.var(f'req:{b}:{m}:{r}:{i}',0,float(hi[i]));rv.append(v);self.lp.add_eq({v:1,**{self.w[b,m,r,o][i]:-1 for o in range(len(O))}},0);self.lp.add_ge({v:1,a:-float(lo[i])},0);self.lp.add_le({v:1,a:-float(hi[i])},0)
     self.req[b,pu,port]=rv
   for u in self.splits:
    aa=int(self.units[u]['left_flow']);bb=int(self.units[u]['right_flow']);m=int(self.units[u]['descendant_count']);act=self.act[u];A=self.req[b,u,'left'];B=self.req[b,u,'right'];lA,uA=I['cb'][aa];lB,uB=I['cb'][bb];zz=[]
    for i in range(d):
     z=self.lp.var(f'z:{b}:{u}:{i}',0,min(float(uA[i]),float(uB[i])));zz.append(z);perspective_mcc(self.lp,z,A[i],B[i],act,float(lA[i]),float(uA[i]),float(lB[i]),float(uB[i]));self.lp.add_le({z:1,A[i]:-1},0);self.lp.add_le({z:1,B[i]:-1},0)
    self.z[b,u]=zz
    if u not in self.roots:
     coef=I['s'][aa]*I['s'][bb]/I['s'][m];P=km([transition_for_count(m)]*len(blk));xv=self.x[b,u]
     for i in range(d):self.lp.add_eq({xv[i]:1,**{zz[j]:-coef*float(P[i,j]) for j in range(d)}},0)
   # normalized root likelihood q = physical beta-scaled block likelihood / sroot
   q=self.lp.var(f'q:{b}',0,1);self.q[b]=q;eq={q:1}
   for u in self.roots:
    z=self.z[b,u];a=self.act[u];coef=I['rootcoef'][u]
    for i in range(d):eq[z[i]]=eq.get(z[i],0)-coef*float(I['pi'][i])
    lo,hi=I['roots'][u];expr={z[i]:coef*float(I['pi'][i]) for i in range(d)};self.lp.add_le({**expr,a:-hi},0);self.lp.add_ge({**expr,a:-lo},0)
   self.lp.add_eq(eq,0)
   t=self.lp.var(f'log:{b}',-1000,10);self.logv[b]=t
   for x0 in (1,.1,.01,.001,.0001):self.addtan(b,x0)
 def addtan(self,b,x0):self.lp.add_le({self.logv[b]:1,self.q[b]:-1/x0},math.log(x0)-1)
 def solve(self,tol=2e-4,maxit=12):
  if len(self.blocks)==1:
   val,res,dt=self.lp.solve({self.q[0]:1},True);phys=math.log(max(val,1e-300))+math.log(self.info[0]['sroot'])+sum(len(bl) for bl in self.blocks)*self.L*math.log(BETA);return phys,res,dt,[]
  obj={self.logv[b]:1 for b in range(len(self.blocks))};tot=0;hist=[]
  for it in range(maxit):
   val,res,dt=self.lp.solve(obj,True);tot+=dt
   if not res.success:return None,res,tot,hist
   gaps=[]
   for b in range(len(self.blocks)):
    q=max(res.x[self.q[b]],1e-12);gaps.append(res.x[self.logv[b]]-math.log(q))
   hist.append((it,max(gaps),dt))
   if max(gaps)<=tol:break
   for b,g in enumerate(gaps):
    if g>tol:self.addtan(b,max(res.x[self.q[b]],1e-12))
  const=sum(math.log(I['sroot']) for I in self.info)+sum(len(bl) for bl in self.blocks)*self.L*math.log(BETA)
  return val+const,res,tot,hist

def main():
 pats=[tuple('AAAACCCC'),tuple('AACCGGTT'),tuple('ACGTACGT')];exact,_,_=base.exact_shared_best(pats,np.ones(3));print('exact',exact)
 for name,blocks in [('1+1+1',[[pats[0]],[pats[1]],[pats[2]]]),('2+1',[[pats[1],pats[2]],[pats[0]]]),('2+1alt',[[pats[0],pats[1]],[pats[2]]]),('3',[[pats[0],pats[1],pats[2]]])]:
  t=time.time();M=Model('/mnt/data/multisite_units8/library.json',blocks);print(name,'size',len(M.lp.names),len(M.lp.eq),len(M.lp.le),'build',time.time()-t,'scales',[(I['sroot'],max(I['rootcoef'].values())) for I in M.info],flush=True)
  ub,res,st,h=M.solve();print(' ->',res.status,res.message,'ub',ub,'gap',None if ub is None else ub-exact,'solve',st,'hist',h,flush=True)
if __name__=='__main__':main()
