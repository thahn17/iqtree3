#!/usr/bin/env python3
"""Compatibility/root-only layer used by multisite_pruned_general.py.

This file reconstructs the compact root-only research base from the retained lower-arm
model.  It is intentionally small and keeps all topology/taxon-flow variables shared;
per-site variables are root likelihood contributions plus log interpolation variables.
"""
from __future__ import annotations
import math, random
from collections import Counter, defaultdict
from pathlib import Path
import json
import numpy as np
from scipy.optimize import linprog
import multisite_lowerarm_taxonflow as base

PI=base.PI


def random_columns(L,n,seed=1):
    rng=random.Random(seed); b='ACGT'
    return [tuple(rng.choice(b) for _ in range(L)) for _ in range(n)]


def facetfun(features, vals, max_facets=8):
    X=np.column_stack([np.ones(len(features)),np.asarray(features,float)]); y=np.asarray(vals,float)
    if len(X)==0:return []
    if len(X)<=max_facets:idx=list(range(len(X)))
    else:
        order=np.argsort(y);idx=set([int(order[0]),int(order[-1])])
        for z in np.linspace(0,len(order)-1,max_facets,dtype=int):idx.add(int(order[z]))
        idx=sorted(idx)
    fs=[]
    for ii in idx:
        res=linprog(X[ii],A_ub=-X,b_ub=-y,bounds=[(None,None)]*X.shape[1],method='highs')
        if res.success and np.all(X@res.x+1e-8>=y):
            f=tuple(np.round(res.x,11))
            if f not in fs:fs.append(f)
    return fs


class RootOnly(base.BoundaryTaxonModel):
    def __init__(self,libpath,columns,log_ratio=1.1,long_ratio=2.0,
                 pair_moments=True,max_first_facets=20,max_quad_facets=3,pair_bounds='full'):
        cnt=Counter(tuple(c) for c in columns);self.patterns=list(cnt);self.weights=np.array([cnt[p] for p in self.patterns],float)
        self.lib=json.loads(Path(libpath).read_text());self.L=int(self.lib['summary']['leaf_count']);self.units=self.lib['units'];self.roots=set(self.lib['root_candidates'])
        self.splits=[u for u,x in self.units.items() if x['kind']=='split'];self.tax=[u for u,x in self.units.items() if x['kind']=='taxon']
        self.by=defaultdict(list)
        for u in self.splits:self.by[int(self.units[u]['descendant_count'])].append(u)
        self.reqdef=defaultdict(list)
        for u in self.splits:
            for port in ('left','right'):
                self.reqdef[int(self.units[u][port+'_flow'])].append((u,port))
        self.log_ratio=log_ratio;self.long_ratio=long_ratio;self.variable_lengths=True;self.taxon_flow=True;self.root_partition=False
        self.max_first_facets=max_first_facets;self.max_quad_facets=max_quad_facets;self.pair_bounds=pair_bounds
        self.lp=base.SparseLP();self.act={u:self.lp.var('act:'+u,0,1) for u in self.splits};self.lp.add_eq({self.act[u]:1 for u in self.roots},1)
        self.wireact={};self.sw={};self.routermeta={};self._build_topology()
        self.branch_relaxation=base.branch_relaxation_mode();self.ell={};self.branch_effect={}
        if self.branch_relaxation=='spectral':
            S=base.spectral_relaxation();self.spectral_mu=np.asarray(S.mu,float)
            for u in self.splits:
                for port in ('left','right'):
                    q=[self.lp.var(f'smode:{u}:{port}:{r}',0,1) for r in range(len(S.mu))];self.branch_effect[u,port]=q
                    Aact=self.act[u]
                    for ii in range(len(S.b)):
                        row={q[r]:float(S.A[ii,r]) for r in range(len(q)) if abs(S.A[ii,r])>1e-15}
                        row[Aact]=row.get(Aact,0)-float(S.b[ii]);self.lp.add_le(row,0)
        else:
            for u in self.splits:
                for port in ('left','right'):
                    e=self.lp.var(f'long:{u}:{port}',0,1);self.ell[u,port]=e;self.lp.add_le({e:1,self.act[u]:-1},0);self.branch_effect[u,port]=[e]
        self.tmem={};self._build_taxon_flow()
        self.info=[];self._prep()
        self.req={};self.hlong={};self.z={};self.q={};self.lam={};self.obj={};self.log_error_bound=0;self.qr={};self.rpair={}
        self._build_sites()
        self._add_root_length_partition_facets()
        if pair_moments:
            self._add_pairs();self._add_quad_facets()

    def _build_sites(self):
        self.qr={};self.lam={};self.obj={};self.log_error_bound=0
        pts=base.vl.geometric_breaks(1e-12,self.log_ratio)
        for p,pat in enumerate(self.patterns):
            I=self.info[p]
            for r in self.roots:
                hi=max(I['rootvals'][r].values())/I['scale'][self.L]
                self.qr[p,r]=self.lp.var(f'qr:{p}:{r}',0,float(hi))
            lam=[]
            for kk,x in enumerate(pts):
                v=self.lp.var(f'lam:{p}:{kk}',0,1);lam.append(v);self.obj[v]=self.weights[p]*math.log(x)
            self.lam[p]=lam;self.lp.add_eq({v:1 for v in lam},1)
            row={self.qr[p,r]:1 for r in self.roots}
            for k,v in enumerate(lam):row[v]=row.get(v,0)-pts[k]
            self.lp.add_eq(row,0);self.log_error_bound+=self.weights[p]*base.vl.chord_gap(self.log_ratio)

    def _first_facets(self,p,r):
        I=self.info[p];counts=I['counts'];a=int(self.units[r]['left_flow']);b=self.L-a;raw=I['raw']
        PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio);X=[];y=[]
        for (mm,c1),A in raw.items():
            if mm!=a:continue
            c2=tuple(counts[k]-c1[k] for k in range(4))
            if min(c2)<0 or (b,c2) not in raw:continue
            for la,PA in enumerate((PsA,PlA)):
                for lb,PB in enumerate((PsB,PlB)):
                    X.append([c1[0],c1[1],c1[2],la,lb]);y.append(float(base.PI@((PA@A[:,1])*(PB@raw[(b,c2)][:,1])))/I['scale'][self.L])
        return facetfun(X,y,self.max_first_facets)

    def _add_root_length_partition_facets(self):
        self.root_length_facets=0
        for p,pat in enumerate(self.patterns):
            for r in self.roots:
                m,s,_=self._req_slot(r,'left');A=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right']
                for f in self._first_facets(p,r):
                    a0,dA,dC,dG,dL,dR=f;row={self.qr[p,r]:1,A:-a0,eL:-dL,eR:-dR}
                    for tax,ch in enumerate(pat):
                        c={'A':dA,'C':dC,'G':dG,'T':0.0}.get(ch,0.0)
                        if c:
                            x=self.tmem[m,0,s,tax];row[x]=row.get(x,0)-c
                    self.lp.add_le(row,0);self.root_length_facets+=1

    def _add_pairs(self):
        for r in self.roots:
            a=int(self.units[r]['left_flow'])
            if a<=1:continue
            A=self.act[r];m,s,_=self._req_slot(r,'left')
            for u in range(self.L):
                for v in range(u+1,self.L):
                    y=self.lp.var(f'rpair:{r}:{u}:{v}',0,1);self.rpair[r,u,v]=y
                    if self.pair_bounds=='full' and a>=3:
                        xu=self.tmem[m,0,s,u];xv=self.tmem[m,0,s,v]
                        self.lp.add_le({y:1,xu:-1},0);self.lp.add_le({y:1,xv:-1},0);self.lp.add_ge({y:1,xu:-1,xv:-1,A:1},0)
            for u in range(self.L):
                row={}
                for v in range(self.L):
                    if v==u:continue
                    row[self.rpair[r,min(u,v),max(u,v)]]=1
                row[self.tmem[m,0,s,u]]=-(a-1);self.lp.add_eq(row,0)

    def _moment_exprs(self,r,pat):
        m,s,_=self._req_slot(r,'left')
        C=[]
        for b in 'ACG':
            d={}
            for t,ch in enumerate(pat):
                if ch==b:d[self.tmem[m,0,s,t]]=1
            C.append(d)
        out=list(C); bases='ACGT'
        for bi,b in enumerate(bases):
            for ci,c in enumerate(bases[bi:],start=bi):
                d={}
                if b==c:
                    for t,ch in enumerate(pat):
                        if ch==b:d[self.tmem[m,0,s,t]]=d.get(self.tmem[m,0,s,t],0)+1
                    ids=[t for t,ch in enumerate(pat) if ch==b]
                    for ii,u in enumerate(ids):
                        for v in ids[ii+1:]:
                            y=self.rpair.get((r,min(u,v),max(u,v)))
                            if y is not None:d[y]=d.get(y,0)+2
                else:
                    ib=[t for t,ch in enumerate(pat) if ch==b];ic=[t for t,ch in enumerate(pat) if ch==c]
                    for u in ib:
                        for v in ic:
                            y=self.rpair.get((r,min(u,v),max(u,v)))
                            if y is not None:d[y]=d.get(y,0)+1
                out.append(d)
        return out

    def _quad_facets(self,p,r):
        I=self.info[p];counts=I['counts'];a=int(self.units[r]['left_flow']);b=self.L-a;raw=I['raw'];X=[];y=[]
        PsA,PlA=base.trans_pair(a,self.long_ratio);PsB,PlB=base.trans_pair(b,self.long_ratio)
        for (mm,c1),A in raw.items():
            if mm!=a:continue
            c2=tuple(counts[k]-c1[k] for k in range(4))
            if min(c2)<0 or (b,c2) not in raw:continue
            feat=[c1[0],c1[1],c1[2]]+[c1[i]*c1[j] for i in range(4) for j in range(i,4)]
            for la,PA in enumerate((PsA,PlA)):
                for lb,PB in enumerate((PsB,PlB)):
                    X.append(feat+[la,lb]);y.append(float(base.PI@((PA@A[:,1])*(PB@raw[(b,c2)][:,1])))/I['scale'][self.L])
        return facetfun(X,y,self.max_quad_facets)

    def _add_quad_facets(self):
        self.root_moment_facets=0
        for p,pat in enumerate(self.patterns):
            for r in self.roots:
                a=int(self.units[r]['left_flow'])
                if a<=1:continue
                ex=self._moment_exprs(r,pat);A=self.act[r];eL=self.ell[r,'left'];eR=self.ell[r,'right']
                for f in self._quad_facets(p,r):
                    a0,*co=f;row={self.qr[p,r]:1,A:-a0,eL:-co[len(ex)],eR:-co[len(ex)+1]}
                    for h,d in enumerate(ex):
                        if co[h]:
                            for v,c in d.items():row[v]=row.get(v,0)-co[h]*c
                    self.lp.add_le(row,0);self.root_moment_facets+=1

    def solve(self):
        val,res,dt=self.lp.solve(self.obj,True)
        const=sum(self.weights[p]*math.log(self.info[p]['scale'][self.L]) for p in range(len(self.patterns)))
        return val+const,res,dt
