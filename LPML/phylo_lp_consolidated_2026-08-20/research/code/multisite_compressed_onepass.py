#!/usr/bin/env python3
"""Compressed one-pass shared-tree multisite LP research harness.

Key reductions:
- shared physical Beneš topology: split activities + router wire activities + switch fractions
- no site-specific router likelihood variables
- no site-specific post-edge X variables (X=PZ substituted)
- site-specific request vectors + split product vectors only
- upper-only count/composition boxes (lower bounds dropped except root likelihood floor)
- duplicate site-pattern compression
- one-pass log objective via geometric chord interpolation (pure LP, no iterative OA)

This is an LP relaxation / research harness, not a finished production solver.
"""
from __future__ import annotations
import json,math,importlib.util
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('em','/mnt/data/multisite_conservation_eliminated.py')
em=importlib.util.module_from_spec(spec);spec.loader.exec_module(em)
base=em.base;h=em.h;PI=h.PI;BETA=h.BETA;BASES=h.BASES


def add_switch_hull(lp,x0,x1,y0,y1,s):
    # Activity-only 2x2 switch hull, all activities in [0,1].
    lp.add_eq({y0:1,y1:1,x0:-1,x1:-1},0)
    # p=x0-y0=s*(x0-x1), with d in [-1,1]
    lp.add_ge({x0:1,y0:-1,s:1},0)            # p >= -s
    lp.add_ge({y0:-1,x1:1,s:-1},-1)          # p >= d+s-1
    lp.add_le({x0:1,y0:-1,s:-1},0)           # p <= s
    lp.add_le({y0:-1,x1:1,s:1},1)            # p <= d-s+1


def geometric_breaks(lo,ratio):
    lo=max(float(lo),1e-300)
    pts=[1.0];x=1.0
    while x/ratio > lo*(1+1e-12):
        x/=ratio;pts.append(x)
    if pts[-1] > lo*(1+1e-12):pts.append(lo)
    return sorted(set(pts))


def chord_gap(r):
    c=math.log(r)/(r-1);y=(r-1)/math.log(r)
    return math.log(y)-(y-1)*c


class Model:
    def __init__(self,libpath,columns,log_ratio=1.1):
        cnt=Counter(tuple(c) for c in columns)
        self.patterns=list(cnt);self.weights=np.array([cnt[p] for p in self.patterns],float)
        self.lib=json.loads(Path(libpath).read_text());self.L=int(self.lib['summary']['leaf_count']);self.units=self.lib['units'];self.roots=set(self.lib['root_candidates'])
        self.splits=[u for u,x in self.units.items() if x['kind']=='split'];self.tax=[u for u,x in self.units.items() if x['kind']=='taxon'];self.by=defaultdict(list)
        for u in self.splits:self.by[int(self.units[u]['descendant_count'])].append(u)
        self.reqdef=defaultdict(list)
        for u in self.splits:
            for port,f in [('left',int(self.units[u]['left_flow'])),('right',int(self.units[u]['right_flow']))]:
                if f<self.L:self.reqdef[f].append((u,port))
        self.lp=base.SparseLP();self.act={u:self.lp.var('act:'+u,0,1) for u in self.splits};self.lp.add_eq({self.act[u]:1 for u in self.roots},1)
        self.wireact={};self.sw={};self._build_physical_topology()
        self.info=[];self._prep()
        self.req={};self.z={};self.q={};self.lam={};self.log_ratio=log_ratio;self._build_sites()
    def _build_physical_topology(self):
        for mstr,r in self.lib['routers'].items():
            m=int(mstr);W=int(r['padded_width']);dims=list(r['dimensions']);D=len(dims)
            if D==0:
                rr=r['input_requests'][0];oo=r['output_targets'][0];a=self.act[rr['parent_unit']]
                if m==1:self.lp.add_eq({a:1},1)
                else:self.lp.add_eq({a:1,self.act[oo['unit_id']]:-1},0)
                continue
            # internal layers only are variables; boundary values are represented by helper fixed/bound references
            # For compactness in this research implementation we still create all layer activities; site variables do not.
            for l in range(D+1):
                for w in range(W):self.wireact[m,l,w]=self.lp.var(f'wa:{m}:{l}:{w}',0,1)
            rin={int(x['slot']):x for x in r['input_requests']};rout={int(x['slot']):x for x in r['output_targets']}
            for w in range(W):
                v=self.wireact[m,0,w]
                if w in rin:self.lp.add_eq({v:1,self.act[rin[w]['parent_unit']]:-1},0)
                else:self.lp.add_eq({v:1},0)
                v=self.wireact[m,D,w]
                if w in rout:
                    if m==1:self.lp.add_eq({v:1},1)
                    else:self.lp.add_eq({v:1,self.act[rout[w]['unit_id']]:-1},0)
                else:self.lp.add_eq({v:1},0)
            for st,dim in enumerate(dims):
                seen=set()
                for w in range(W):
                    b=min(w,w^(1<<dim))
                    if b in seen:continue
                    seen.add(b);o=b^(1<<dim);s=self.lp.var(f'sw:{m}:{st}:{b}',0,1);self.sw[m,st,b]=s
                    add_switch_hull(self.lp,self.wireact[m,st,b],self.wireact[m,st,o],self.wireact[m,st+1,b],self.wireact[m,st+1,o],s)
    def _prep(self):
        for pat in self.patterns:
            counts=tuple(pat.count(b) for b in BASES);out,_=h.build_composition_boxes_restricted(counts);cb={};scale={}
            for m in range(1,self.L):
                Bs=[B for (mm,c),B in out.items() if mm==m];hi=np.max(np.stack([B[:,1] for B in Bs]),axis=0);cb[m]=hi;scale[m]=float(hi.max())
            cbn={m:cb[m]/scale[m] for m in cb};roots={u:(base.root_split_lower(counts,out,int(self.units[u]['split_a'])),h.root_split_upper(counts,out,int(self.units[u]['split_a']))[0]) for u in self.roots};rsc={u:scale[int(self.units[u]['left_flow'])]*scale[int(self.units[u]['right_flow'])] for u in self.roots};sroot=max(rsc.values())
            self.info.append({'u':cbn,'s':scale,'roots':{u:(roots[u][0]/sroot,roots[u][1]/sroot) for u in self.roots},'rootcoef':{u:rsc[u]/sroot for u in self.roots},'sroot':sroot})
    def post_expr(self,p,u,i):
        I=self.info[p];a=int(self.units[u]['left_flow']);b=int(self.units[u]['right_flow']);m=int(self.units[u]['descendant_count']);coef=I['s'][a]*I['s'][b]/I['s'][m];P=h.transition_for_count(m)
        return {self.z[p,u][j]:coef*float(P[i,j]) for j in range(4)}
    def _build_sites(self):
        self.obj={};self.log_error_bound=0.0
        for p,pat in enumerate(self.patterns):
            I=self.info[p]
            # request vectors, upper-only
            for m in range(1,self.L):
                hi=I['u'][m]
                for r,(pu,port) in enumerate(self.reqdef[m]):
                    a=self.act[pu];rv=[]
                    for i in range(4):
                        v=self.lp.var(f'req:{p}:{m}:{r}:{i}',0,float(hi[i]));rv.append(v);self.lp.add_le({v:1,a:-float(hi[i])},0)
                    self.req[p,pu,port]=rv
            # split products with [0,u] perspective hull: 3 inequalities/state
            for u in self.splits:
                aa=int(self.units[u]['left_flow']);bb=int(self.units[u]['right_flow']);a=self.act[u];A=self.req[p,u,'left'];B=self.req[p,u,'right'];uA=I['u'][aa];uB=I['u'][bb];zz=[]
                for i in range(4):
                    z=self.lp.var(f'z:{p}:{u}:{i}',0,min(float(uA[i]),float(uB[i])));zz.append(z)
                    self.lp.add_ge({z:1,A[i]:-float(uB[i]),B[i]:-float(uA[i]),a:float(uA[i]*uB[i])},0)
                    self.lp.add_le({z:1,A[i]:-float(uB[i])},0);self.lp.add_le({z:1,B[i]:-float(uA[i])},0)
                self.z[p,u]=zz
            # router likelihood projection: only total conservation per fixed count, post-edge X substituted
            for m in range(1,self.L):
                for i in range(4):
                    d={self.req[p,pu,port][i]:1 for pu,port in self.reqdef[m]}
                    if m==1:
                        total=sum((h.scaled_leaf_vector(pat[int(self.units[uid]['label'].split('_')[-1])-1])/I['s'][1])[i] for uid in self.tax);self.lp.add_eq(d,total)
                    else:
                        for u in self.by[m]:
                            for v,c in self.post_expr(p,u,i).items():d[v]=d.get(v,0)-c
                        self.lp.add_eq(d,0)
            # root likelihood and split-aware upper cuts
            q=self.lp.var(f'q:{p}',0,1);self.q[p]=q;eq={q:1}
            for u in self.roots:
                coef=I['rootcoef'][u];z=self.z[p,u];a=self.act[u]
                for i in range(4):eq[z[i]]=eq.get(z[i],0)-coef*float(PI[i])
                hi=I['roots'][u][1];expr={z[i]:coef*float(PI[i]) for i in range(4)};self.lp.add_le({**expr,a:-hi},0)
            self.lp.add_eq(eq,0)
            # one-pass chord log, using a valid physical root lower bound only for interpolation domain
            qlo=max(min(v[0] for v in I['roots'].values()),1e-300);pts=geometric_breaks(qlo,self.log_ratio);lam=[]
            for k,x in enumerate(pts):
                v=self.lp.var(f'lam:{p}:{k}',0,1);lam.append(v);self.obj[v]=self.weights[p]*math.log(x)
            self.lam[p]=lam;self.lp.add_eq({v:1 for v in lam},1);self.lp.add_eq({q:1,**{lam[k]:-pts[k] for k in range(len(pts))}},0)
            self.log_error_bound += self.weights[p]*chord_gap(self.log_ratio)
    def solve(self):
        val,res,dt=self.lp.solve(self.obj,True)
        const=sum(self.weights[p]*(math.log(self.info[p]['sroot'])+self.L*math.log(BETA)) for p in range(len(self.patterns)))
        return val+const,res,dt

if __name__=='__main__':
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    M=Model('/mnt/data/multisite_units8/library.json',cols,1.1);v,r,t=M.solve();print('obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'solve',t,'PWLerr<=',M.log_error_bound)
