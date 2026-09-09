#!/usr/bin/env python3
"""Compressed shared-tree multisite LP with independent short/long branch choices.

Research harness. Pure continuous LP.

Length representation:
- one shared long-activity variable ell_e for each likelihood-bearing edge;
- leaves: message is affine in ell_e because the observed base is fixed;
- non-root split unit u: one four-state long-part vector H_{site,u}, with
    0 <= H <= Z,
    H <= U_Z * ell_u,
    Z-H <= U_Z * (act_u-ell_u),
  followed by
    post = coef * [P_short Z + (P_long-P_short) H].

Thus branch-length flexibility adds only 4 variables per non-root split per
unique site pattern, plus one shared scalar per possible biological edge.
"""
from __future__ import annotations

import importlib.util
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

# Shared configurable substitution model (JC69 by default; GTR+F supported).
import sys
_HERE=Path(__file__).resolve().parent
if str(_HERE) not in sys.path: sys.path.insert(0,str(_HERE))
from phylo_substitution import get_active_model, BASES, BASE_INDEX
_m=get_active_model(); PI=_m.pi.copy(); Q=_m.q.copy()

def configure_substitution_model(model):
    global PI,Q
    PI=np.asarray(model.pi,float).copy(); Q=np.asarray(model.q,float).copy()


class SparseLP:
    def __init__(self):
        self.names=[]; self.bounds=[]
        self.le=[]; self.lerhs=[]; self.eq=[]; self.eqrhs=[]
    def var(self,name,lb=0.0,ub=None):
        j=len(self.names); self.names.append(name); self.bounds.append((lb,ub)); return j
    def add_le(self,d,rhs): self.le.append(dict(d)); self.lerhs.append(float(rhs))
    def add_ge(self,d,rhs): self.le.append({j:-v for j,v in d.items()}); self.lerhs.append(float(-rhs))
    def add_eq(self,d,rhs): self.eq.append(dict(d)); self.eqrhs.append(float(rhs))
    def _mat(self,rows):
        if not rows: return None
        data=[]; rr=[]; cc=[]
        for r,d in enumerate(rows):
            for j,v in d.items():
                if abs(v)>0:
                    rr.append(r); cc.append(j); data.append(v)
        return csr_matrix((data,(rr,cc)),shape=(len(rows),len(self.names)))
    def solve(self,obj,maximize=True):
        c=np.zeros(len(self.names))
        for j,v in obj.items(): c[j]=-v if maximize else v
        t=time.time()
        res=linprog(
            c,A_ub=self._mat(self.le),b_ub=np.asarray(self.lerhs) if self.le else None,
            A_eq=self._mat(self.eq),b_eq=np.asarray(self.eqrhs) if self.eq else None,
            bounds=self.bounds,method="highs",
            options={"dual_feasibility_tolerance":1e-8,"primal_feasibility_tolerance":1e-8},
        )
        dt=time.time()-t
        if not res.success:
            raise RuntimeError(res.message)
        val=sum(v*res.x[j] for j,v in obj.items())
        return val,res,dt


def add_switch_hull(lp,x0,x1,y0,y1,s):
    lp.add_eq({y0:1,y1:1,x0:-1,x1:-1},0)
    lp.add_ge({x0:1,y0:-1,s:1},0)
    lp.add_ge({y0:-1,x1:1,s:-1},-1)
    lp.add_le({x0:1,y0:-1,s:-1},0)
    lp.add_le({y0:-1,x1:1,s:1},1)


def geometric_breaks(lo,ratio):
    lo=max(float(lo),1e-300)
    pts=[1.0]; x=1.0
    while x/ratio > lo*(1+1e-12):
        x/=ratio; pts.append(x)
    if pts[-1] > lo*(1+1e-12): pts.append(lo)
    return sorted(set(pts))


def chord_gap(r):
    c=math.log(r)/(r-1); y=(r-1)/math.log(r)
    return math.log(y)-(y-1)*c


def short_time(m:int)->float:
    return 0.10 if m == 1 else 0.055 + 0.018 * math.log2(m)


def transition_pair(m:int,long_ratio:float):
    ts=short_time(m)
    return expm(Q*ts), expm(Q*(ts*long_ratio))


def build_boxes(global_counts:Tuple[int,int,int,int],long_ratio:float,variable_lengths:bool):
    """Composition-restricted component boxes allowing either length at each edge."""
    g=tuple(map(int,global_counts)); L=sum(g)
    P1s,P1l=transition_pair(1,long_ratio)
    beta=max(float(P1s.max()), float(P1l.max())) if variable_lengths else float(P1s.max())
    out={}; partial={}; by=defaultdict(list)
    for b in range(4):
        if g[b]<=0: continue
        c=tuple(1 if i==b else 0 for i in range(4))
        vs=P1s[:,b]/beta
        if variable_lengths:
            vl=P1l[:,b]/beta; lo=np.minimum(vs,vl); hi=np.maximum(vs,vl)
        else:
            lo=hi=vs
        out[(1,c)]=np.stack([lo,hi],axis=1); by[1].append(c)
    for m in range(2,L+1):
        acc_o={}; acc_z={}; Ps,Pl=transition_pair(m,long_ratio)
        for a in range(1,m//2+1):
            b=m-a
            if not by[a] or not by[b]: continue
            for c1 in by[a]:
                B1=out[(a,c1)]; lo1,hi1=B1[:,0],B1[:,1]
                for c2 in by[b]:
                    c=tuple(c1[i]+c2[i] for i in range(4))
                    if any(c[i]>g[i] for i in range(4)): continue
                    B2=out[(b,c2)]
                    zlo=lo1*B2[:,0]; zhi=hi1*B2[:,1]
                    los=[Ps@zlo]; his=[Ps@zhi]
                    if variable_lengths:
                        los.append(Pl@zlo); his.append(Pl@zhi)
                    olo=np.min(np.stack(los),axis=0); ohi=np.max(np.stack(his),axis=0)
                    if c not in acc_o:
                        acc_o[c]=[olo.copy(),ohi.copy()]; acc_z[c]=[zlo.copy(),zhi.copy()]
                    else:
                        acc_o[c][0]=np.minimum(acc_o[c][0],olo); acc_o[c][1]=np.maximum(acc_o[c][1],ohi)
                        acc_z[c][0]=np.minimum(acc_z[c][0],zlo); acc_z[c][1]=np.maximum(acc_z[c][1],zhi)
        by[m]=list(acc_o)
        for c,(lo,hi) in acc_o.items(): out[(m,c)]=np.stack([lo,hi],axis=1)
        for c,(lo,hi) in acc_z.items(): partial[(m,c)]=np.stack([lo,hi],axis=1)
    return out,partial,beta


def root_split_bounds(global_counts,out_boxes,split_a:int):
    g=tuple(map(int,global_counts)); L=sum(g); a=min(int(split_a),L-int(split_a)); b=L-a
    by=defaultdict(list)
    for m,c in out_boxes: by[m].append(c)
    lo_best=math.inf; hi_best=-math.inf
    for c1 in by[a]:
        c2=tuple(g[i]-c1[i] for i in range(4))
        if min(c2)<0 or (b,c2) not in out_boxes: continue
        B1=out_boxes[(a,c1)]; B2=out_boxes[(b,c2)]
        lo=float(PI@(B1[:,0]*B2[:,0])); hi=float(PI@(B1[:,1]*B2[:,1]))
        lo_best=min(lo_best,lo); hi_best=max(hi_best,hi)
    if not math.isfinite(lo_best): lo_best=1e-300
    return lo_best,hi_best


class Model:
    def __init__(self,libpath,columns,log_ratio=1.1,variable_lengths=True):
        cnt=Counter(tuple(c) for c in columns)
        self.patterns=list(cnt); self.weights=np.array([cnt[p] for p in self.patterns],float)
        self.lib=json.loads(Path(libpath).read_text()); self.L=int(self.lib['summary']['leaf_count'])
        self.units=self.lib['units']; self.roots=set(self.lib['root_candidates'])
        self.splits=[u for u,x in self.units.items() if x['kind']=='split']
        self.tax=[u for u,x in self.units.items() if x['kind']=='taxon']
        self.nonroot=[u for u in self.splits if u not in self.roots]
        self.by=defaultdict(list)
        for u in self.splits:self.by[int(self.units[u]['descendant_count'])].append(u)
        self.reqdef=defaultdict(list)
        for u in self.splits:
            for port,f in [('left',int(self.units[u]['left_flow'])),('right',int(self.units[u]['right_flow']))]:
                if f<self.L:self.reqdef[f].append((u,port))
        groups=self.lib.get('length_choice_groups',{})
        ratios={float(g['length_ratio']) for g in groups.values() if g.get('likelihood_bearing')}
        self.long_ratio=next(iter(ratios)) if ratios else 2.0
        if len(ratios)>1: raise ValueError('Research model currently expects one short/long path ratio')
        self.variable_lengths=bool(variable_lengths)
        self.lp=SparseLP()
        self.act={u:self.lp.var('act:'+u,0,1) for u in self.splits}
        self.lp.add_eq({self.act[u]:1 for u in self.roots},1)
        self.wireact={}; self.sw={}; self._build_physical_topology()
        # Shared branch-length long activities. One per possible biological edge only.
        self.ell={}
        if self.variable_lengths:
            for uid in self.tax:
                self.ell[uid]=self.lp.var('long:'+uid,0,1)
            for uid in self.nonroot:
                e=self.lp.var('long:'+uid,0,1); self.ell[uid]=e
        self.info=[]; self._prep()
        self.req={}; self.z={}; self.hlong={}; self.q={}; self.lam={}; self.log_ratio=log_ratio
        self._build_sites()

    def _build_physical_topology(self):
        for mstr,r in self.lib['routers'].items():
            m=int(mstr); W=int(r['padded_width']); dims=list(r['dimensions']); D=len(dims)
            if D==0:
                rr=r['input_requests'][0]; oo=r['output_targets'][0]; a=self.act[rr['parent_unit']]
                if m==1:self.lp.add_eq({a:1},1)
                else:self.lp.add_eq({a:1,self.act[oo['unit_id']]:-1},0)
                continue
            for l in range(D+1):
                for w in range(W):self.wireact[m,l,w]=self.lp.var(f'wa:{m}:{l}:{w}',0,1)
            rin={int(x['slot']):x for x in r['input_requests']}; rout={int(x['slot']):x for x in r['output_targets']}
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
                    if b in seen: continue
                    seen.add(b); o=b^(1<<dim)
                    s=self.lp.var(f'sw:{m}:{st}:{b}',0,1); self.sw[m,st,b]=s
                    add_switch_hull(self.lp,self.wireact[m,st,b],self.wireact[m,st,o],self.wireact[m,st+1,b],self.wireact[m,st+1,o],s)

    def _prep(self):
        betas=[]
        for pat in self.patterns:
            counts=tuple(pat.count(b) for b in BASES)
            out,partial,beta=build_boxes(counts,self.long_ratio,self.variable_lengths); betas.append(beta)
            cb={}; scale={}
            for m in range(1,self.L):
                Bs=[B for (mm,c),B in out.items() if mm==m]
                hi=np.max(np.stack([B[:,1] for B in Bs]),axis=0)
                cb[m]=hi; scale[m]=float(hi.max())
            cbn={m:cb[m]/scale[m] for m in cb}
            roots={u:root_split_bounds(counts,out,int(self.units[u]['split_a'])) for u in self.roots}
            rsc={u:scale[int(self.units[u]['left_flow'])]*scale[int(self.units[u]['right_flow'])] for u in self.roots}
            sroot=max(rsc.values())
            self.info.append({'u':cbn,'s':scale,'roots':{u:(roots[u][0]/sroot,roots[u][1]/sroot) for u in self.roots},'rootcoef':{u:rsc[u]/sroot for u in self.roots},'sroot':sroot,'beta':beta})
        if max(betas)-min(betas)>1e-12: raise AssertionError('beta should not depend on site')
        self.beta=betas[0]

    def post_expr(self,p,u,i):
        I=self.info[p]; aa=int(self.units[u]['left_flow']); bb=int(self.units[u]['right_flow']); m=int(self.units[u]['descendant_count'])
        coef=I['s'][aa]*I['s'][bb]/I['s'][m]; Ps,Pl=transition_pair(m,self.long_ratio)
        d={self.z[p,u][j]:coef*float(Ps[i,j]) for j in range(4)}
        if self.variable_lengths:
            H=self.hlong[p,u]
            for j in range(4): d[H[j]]=d.get(H[j],0)+coef*float(Pl[i,j]-Ps[i,j])
        return d

    def _build_sites(self):
        self.obj={}; self.log_error_bound=0.0
        P1s,P1l=transition_pair(1,self.long_ratio)
        for p,pat in enumerate(self.patterns):
            I=self.info[p]
            # requests
            for m in range(1,self.L):
                hi=I['u'][m]
                for r,(pu,port) in enumerate(self.reqdef[m]):
                    a=self.act[pu]; rv=[]
                    for i in range(4):
                        v=self.lp.var(f'req:{p}:{m}:{r}:{i}',0,float(hi[i])); rv.append(v)
                        self.lp.add_le({v:1,a:-float(hi[i])},0)
                    self.req[p,pu,port]=rv
            # split product vectors
            for u in self.splits:
                aa=int(self.units[u]['left_flow']); bb=int(self.units[u]['right_flow']); a=self.act[u]
                A=self.req[p,u,'left']; B=self.req[p,u,'right']; uA=I['u'][aa]; uB=I['u'][bb]; zz=[]
                for i in range(4):
                    ub=min(float(uA[i]),float(uB[i])); z=self.lp.var(f'z:{p}:{u}:{i}',0,ub); zz.append(z)
                    self.lp.add_ge({z:1,A[i]:-float(uB[i]),B[i]:-float(uA[i]),a:float(uA[i]*uB[i])},0)
                    self.lp.add_le({z:1,A[i]:-float(uB[i])},0); self.lp.add_le({z:1,B[i]:-float(uA[i])},0)
                self.z[p,u]=zz
                # Long-path share only for non-root upper branches.
                if self.variable_lengths and u in self.nonroot:
                    ell=self.ell[u]; hh=[]
                    for i,z in enumerate(zz):
                        ub=min(float(uA[i]),float(uB[i])); hv=self.lp.var(f'h:{p}:{u}:{i}',0,ub); hh.append(hv)
                        self.lp.add_le({hv:1,z:-1},0)                         # H <= Z
                        self.lp.add_le({hv:1,ell:-ub},0)                     # H <= ub*ell
                        self.lp.add_le({z:1,hv:-1,a:-ub,ell:ub},0)          # Z-H <= ub*(a-ell)
                    self.hlong[p,u]=hh
            # count-flow likelihood conservation
            for m in range(1,self.L):
                for i in range(4):
                    d={self.req[p,pu,port][i]:1 for pu,port in self.reqdef[m]}
                    if m==1:
                        rhs=0.0
                        for uid in self.tax:
                            idx=int(self.units[uid]['label'].split('_')[-1])-1; b=BASE_INDEX[pat[idx]]
                            vs=P1s[:,b]/self.beta/I['s'][1]
                            rhs += float(vs[i])
                            if self.variable_lengths:
                                vl=P1l[:,b]/self.beta/I['s'][1]
                                delta=float(vl[i]-vs[i]); e=self.ell[uid]
                                d[e]=d.get(e,0)-delta
                        self.lp.add_eq(d,rhs)
                    else:
                        for u in self.by[m]:
                            for v,c in self.post_expr(p,u,i).items(): d[v]=d.get(v,0)-c
                        self.lp.add_eq(d,0)
            # root q and split-aware bounds
            q=self.lp.var(f'q:{p}',0,1); self.q[p]=q; eq={q:1}
            for u in self.roots:
                coef=I['rootcoef'][u]; z=self.z[p,u]; a=self.act[u]
                for i in range(4): eq[z[i]]=eq.get(z[i],0)-coef*float(PI[i])
                hi=I['roots'][u][1]; expr={z[i]:coef*float(PI[i]) for i in range(4)}
                self.lp.add_le({**expr,a:-hi},0)
            self.lp.add_eq(eq,0)
            qlo=max(min(v[0] for v in I['roots'].values()),1e-300); pts=geometric_breaks(qlo,self.log_ratio); lam=[]
            for k,x in enumerate(pts):
                v=self.lp.var(f'lam:{p}:{k}',0,1); lam.append(v); self.obj[v]=self.weights[p]*math.log(x)
            self.lam[p]=lam; self.lp.add_eq({v:1 for v in lam},1); self.lp.add_eq({q:1,**{lam[k]:-pts[k] for k in range(len(pts))}},0)
            self.log_error_bound += self.weights[p]*chord_gap(self.log_ratio)

    def solve(self):
        val,res,dt=self.lp.solve(self.obj,True)
        const=sum(self.weights[p]*(math.log(self.info[p]['sroot'])+self.L*math.log(self.beta)) for p in range(len(self.patterns)))
        return val+const,res,dt

    def length_solution(self,res,threshold=1e-8):
        rows=[]
        if not self.variable_lengths: return rows
        for uid,e in self.ell.items():
            if uid in self.tax:
                a=1.0
            else:
                a=float(res.x[self.act[uid]])
            ell=float(res.x[e]); frac=ell/a if a>threshold else float('nan')
            m=int(self.units[uid]['descendant_count']); ts=short_time(m); tl=ts*self.long_ratio
            teff=ts+(0 if not math.isfinite(frac) else frac)*(tl-ts)
            rows.append((uid,a,ell,frac,ts,tl,teff))
        return rows


def demo(libpath='/mnt/data/varlen_units8/library.json'):
    cols=[tuple('AAAACCCC'),tuple('ACGTACGT'),tuple('AACCGGTT'),tuple('AGCTAGCT'),tuple('AAAACCGT'),tuple('AACCGGTT')]
    for var in (False,True):
        M=Model(libpath,cols,1.1,var); v,r,t=M.solve()
        print(('variable' if var else 'fixed'), 'obj',v,'vars',len(M.lp.names),'eq',len(M.lp.eq),'le',len(M.lp.le),'solve',t)
        if var:
            sol=[x for x in M.length_solution(r) if (x[1]>.05)]
            print('active length sample',sol[:12])

if __name__=='__main__': demo()
