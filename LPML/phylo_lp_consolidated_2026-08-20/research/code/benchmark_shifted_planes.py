#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,sys,time,math,csv
import numpy as np
spec=importlib.util.spec_from_file_location('t','/mnt/data/test_rootonly_pruning.py');t=importlib.util.module_from_spec(spec);sys.modules['t']=t;spec.loader.exec_module(t)
base=t.base;PI=base.PI

def shifted_planes(X,y,K):
    X=np.asarray(X,float);y=np.asarray(y,float);n,d=X.shape
    # normalized coordinates for distances
    mu=X.mean(0);sd=X.std(0);sd[sd<1e-9]=1;Z=(X-mu)/sd
    # representative centers: farthest-point sampling, deterministic
    centers=[int(np.argmax(np.sum((Z-Z.mean(0))**2,axis=1)))] if n else []
    mind=np.sum((Z-Z[centers[0]])**2,axis=1) if centers else np.array([])
    while len(centers)<min(K,n):
        j=int(np.argmax(mind));centers.append(j);mind=np.minimum(mind,np.sum((Z-Z[j])**2,axis=1))
    A=np.column_stack([np.ones(n),X]);planes=[]
    # bandwidth chosen so each local fit sees a meaningful fraction of points
    for j in centers:
        dist=np.sum((Z-Z[j])**2,axis=1);q=np.quantile(dist,0.35) if n>4 else max(float(np.max(dist)),1.0);q=max(float(q),1e-6)
        w=np.exp(-dist/q);Aw=A*np.sqrt(w[:,None]);yw=y*np.sqrt(w)
        coef,*_=np.linalg.lstsq(Aw,yw,rcond=None)
        pred=A@coef;coef[0]+=float(np.max(y-pred))+1e-12
        key=tuple(np.round(coef,12))
        if key not in [tuple(np.round(p,12)) for p in planes]:planes.append(coef.copy())
    # always include global LS shifted plane
    coef,*_=np.linalg.lstsq(A,y,rcond=None);pred=A@coef;coef[0]+=float(np.max(y-pred))+1e-12;planes.append(coef)
    # dedup
    out=[];seen=set()
    for p in planes:
        k=tuple(np.round(p,11))
        if k not in seen:seen.add(k);out.append(tuple(map(float,p)))
    return out

class RootOnlyShift(t.RootOnly):
    def __init__(self,*args,K=8,**kwargs):self._K=K;super().__init__(*args,**kwargs)
    def _add_root_length_partition_facets(self):
        self.root_length_facets=0;self.shift_fit_s=0;self.root_points=0
        for p,pat in enumerate(self.patterns):
            I=self.info[p];counts=I['counts'];raw=I['raw']
            for r in self.roots:
                a=int(self.units[r]['left_flow']);b=self.L-a;PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
                X=[];vals=[]
                for (mm,c1),A in raw.items():
                    if mm!=a:continue
                    c2=tuple(counts[k]-c1[k] for k in range(4))
                    if min(c2)<0 or (b,c2) not in raw:continue
                    B=raw[(b,c2)]
                    for la,PA in enumerate((PsA,PlA)):
                        for lb,PB in enumerate((PsB,PlB)):
                            X.append((c1[0],c1[1],c1[2],la,lb));vals.append(float(PI@((PA@A[:,1])*(PB@B[:,1])))/I['scale'][self.L])
                self.root_points+=len(X);st=time.perf_counter();fs=shifted_planes(X,vals,self._K);self.shift_fit_s+=time.perf_counter()-st
                m,s,_=self._req_slot(r,'left');Aact=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right']
                for f in fs:
                    a0,dA,dC,dG,dL,dR=f;row={self.qr[p,r]:1,Aact:-a0,eL:-dL,eR:-dR}
                    for tax,ch in enumerate(pat):
                        c={'A':dA,'C':dC,'G':dG,'T':0}[ch]
                        if c:
                            x=self.tmem[m,0,s,tax];row[x]=row.get(x,0)-c
                    self.lp.add_le(row,0);self.root_length_facets+=1

if __name__=='__main__':
 for L,lib,cols in [(8,'/mnt/data/test_match_units8/library.json',[tuple('ACGTACGT')]),(12,'/mnt/data/test_match_units12/library.json',[t.random_columns(12,1,702)[0]]),(16,'/mnt/data/test_match_units16/library.json',[t.random_columns(16,1,703)[0]])]:
  for K in [1,4,8,16]:
   st=time.perf_counter();M=RootOnlyShift(lib,cols,1.1,2.0,max_first_facets=0,max_quad_facets=1,K=K);build=time.perf_counter()-st
   v,_,sol=M.solve();print(L,K,v,build,sol,len(M.lp.names),len(M.lp.le),M.root_length_facets,M.shift_fit_s,M.root_points,flush=True)
