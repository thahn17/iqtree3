#!/usr/bin/env python3
from __future__ import annotations
import importlib.util,sys,time,math,itertools,csv
from collections import defaultdict
import numpy as np

spec=importlib.util.spec_from_file_location('base','/mnt/data/multisite_lowerarm_taxonflow.py')
base=importlib.util.module_from_spec(spec);sys.modules['base']=base;spec.loader.exec_module(base)
PI=base.PI


def sig_leaf(b,T):
    s=[0]*4; s[b]=1 if T>1 else T
    return tuple(s)

def add_status(x,y,T):
    # values 0..T-1 exact, T means >=T
    if x==T or y==T: return T
    z=x+y
    return z if z<T else T

def combine_sig(a,b,T): return tuple(add_status(a[i],b[i],T) for i in range(4))

def status_lohi(s,T,g):
    lo=[];hi=[]
    for i,x in enumerate(s):
        if x<T:lo.append(x);hi.append(x)
        else:lo.append(T);hi.append(g[i])
    return lo,hi

def sig_possible(s,m,T,g):
    lo,hi=status_lohi(s,T,g)
    if any(lo[i]>g[i] for i in range(4)): return False
    return sum(lo)<=m<=sum(hi)

def pair_compatible(s1,s2,a,b,T,g):
    lo1,hi1=status_lohi(s1,T,g);lo2,hi2=status_lohi(s2,T,g)
    if not(sum(lo1)<=a<=sum(hi1) and sum(lo2)<=b<=sum(hi2)):return False
    for i in range(4):
        if g[i] < lo1[i]+lo2[i] or g[i] > hi1[i]+hi2[i]: return False
    return True

def saturated_boxes(global_counts,T,long_ratio=2.0):
    g=tuple(map(int,global_counts));L=sum(g)
    raw={};by=defaultdict(list)
    for b in range(4):
        if g[b]<=0: continue
        s=sig_leaf(b,T);v=np.zeros(4);v[b]=1.0
        raw[(1,s)]=v.copy();
        if s not in by[1]:by[1].append(s)
    pair_ops=0
    for m in range(2,L+1):
        acc={}
        for a in range(1,m//2+1):
            b=m-a
            if not by[a] or not by[b]:continue
            PaS,PaL=base.trans_pair(a,long_ratio);PbS,PbL=base.trans_pair(b,long_ratio)
            # cache transformed upper vectors by signature for this child count
            TA={s:np.maximum(PaS@raw[(a,s)],PaL@raw[(a,s)]) for s in by[a]}
            TB={s:np.maximum(PbS@raw[(b,s)],PbL@raw[(b,s)]) for s in by[b]}
            for s1 in by[a]:
                x=TA[s1]
                for s2 in by[b]:
                    pair_ops+=1
                    s=combine_sig(s1,s2,T)
                    if not sig_possible(s,m,T,g):continue
                    z=x*TB[s2]
                    if s not in acc:acc[s]=z.copy()
                    else:acc[s]=np.maximum(acc[s],z)
        by[m]=list(acc)
        for s,v in acc.items():raw[(m,s)]=v
    return raw,by,pair_ops

def saturated_root_upper(g,raw,by,T,long_ratio=2.0):
    L=sum(g);best=-1;best_a=None;root_pairs=0
    for a in range(1,L//2+1):
        b=L-a;PaS,PaL=base.trans_pair(a,long_ratio);PbS,PbL=base.trans_pair(b,long_ratio)
        TA={s:[PaS@raw[(a,s)],PaL@raw[(a,s)]] for s in by[a]}
        TB={s:[PbS@raw[(b,s)],PbL@raw[(b,s)]] for s in by[b]}
        for s1 in by[a]:
            for s2 in by[b]:
                if not pair_compatible(s1,s2,a,b,T,g):continue
                root_pairs+=1
                for x in TA[s1]:
                    for y in TB[s2]:
                        v=float(PI@(x*y))
                        if v>best:best=v;best_a=a
    return best,best_a,root_pairs

def exact_root_upper(g,raw,long_ratio=2.0):
    L=sum(g);by=defaultdict(list)
    for (m,c) in raw:by[m].append(c)
    best=-1;ba=None;root_pairs=0
    for a in range(1,L//2+1):
        b=L-a;PaS,PaL=base.trans_pair(a,long_ratio);PbS,PbL=base.trans_pair(b,long_ratio)
        for c1 in by[a]:
            c2=tuple(g[i]-c1[i] for i in range(4))
            if min(c2)<0 or (b,c2) not in raw:continue
            root_pairs+=1
            for PA in (PaS,PaL):
                x=PA@raw[(a,c1)][:,1]
                for PB in (PbS,PbL):
                    y=PB@raw[(b,c2)][:,1]
                    v=float(PI@(x*y))
                    if v>best:best=v;ba=a
    return best,ba,root_pairs

def run_case(g,Ts,do_exact=True):
    rows=[];exact=None
    if do_exact:
        st=time.perf_counter();raw=base.raw_boxes(g,2.0,True);tdp=time.perf_counter()-st
        st=time.perf_counter();ub,ba,rp=exact_root_upper(g,raw);tr=time.perf_counter()-st
        exact=ub
        rows.append(dict(L=sum(g),counts=str(g),method='exact_composition',T='',prep_dp_s=tdp,root_eval_s=tr,total_pre_s=tdp+tr,upper=ub,ratio_to_exact=1.0,states=len(raw),pair_ops='',root_pairs=rp,best_split=ba))
    for T in Ts:
        st=time.perf_counter();raw,by,ops=saturated_boxes(g,T);tdp=time.perf_counter()-st
        st=time.perf_counter();ub,ba,rp=saturated_root_upper(g,raw,by,T);tr=time.perf_counter()-st
        rows.append(dict(L=sum(g),counts=str(g),method=f'saturated_T{T}',T=T,prep_dp_s=tdp,root_eval_s=tr,total_pre_s=tdp+tr,upper=ub,ratio_to_exact=(ub/exact if exact else ''),states=len(raw),pair_ops=ops,root_pairs=rp,best_split=ba))
    return rows

if __name__=='__main__':
    rows=[]
    cases=[(8,(2,2,2,2),True),(16,(4,4,4,4),True),(24,(6,6,6,6),True),(32,(8,8,8,8),True),
           (32,(16,16,0,0),True),(64,(16,16,16,16),False),(100,(25,25,25,25),False)]
    for L,g,ex in cases:
        print('CASE',L,g,flush=True)
        rs=run_case(g,[1,2,3,4],ex);rows.extend(rs)
        for r in rs:print(r,flush=True)
    fn='/mnt/data/preprocess_tradeoff_saturated.csv'
    with open(fn,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('WROTE',fn)
