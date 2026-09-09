#!/usr/bin/env python3
"""Top-down / leaves-up decoding directly from fractional subtree-unit frequencies.

This decoder implements the hierarchical interpretation:

* The virtual root itself is a 100% match.  A root candidate contributes its solved
  activity frequency.
* For an internal lower request of flow k>1, any unused count-k split unit is a
  structurally valid child.  The child's residual activity is its frequency support.
  We do NOT require a literal request->unit fractional edge or an already resolved
  leaf set.  Taxon-membership overlap is only a secondary matching score.
* At flow 1, the child is a named taxon and the boundary taxon-membership variable is
  the actual frequency support.  Every taxon is used exactly once in a decoded tree.

Thus a fractional mixture of internal subtrees receives full structural credit until
it is recursively resolved.  This is exactly the case where strict named-taxon
matching at every internal port is unnecessarily strong.

A complete tree component has extractable frequency

  theta = min(root/split-unit residual activities, selected flow-1 taxon masses).

When theta is subtracted from every selected split unit and selected leaf attachment,
fixed-flow activity conservation remains balanced: a full integral tree removes the
same number of count-k request and count-k child frequencies for every k.  Hence the
component proportions are meaningful residual frequencies and sum to at most one.

Both root-down and leaves-up beam searches use theta as a monotone bound.  No NNI/SPR
walk and no mixture LP/QP are used.
"""
from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union
import math, random
import numpy as np

Tree=Union[int,Tuple['Tree','Tree']]
Request=Tuple[str,str]

def _taxa(t:Tree):
    if isinstance(t,int):return (t,)
    return tuple(sorted(_taxa(t[0])+_taxa(t[1])))
def _canon(t:Tree):
    if isinstance(t,int):return t
    a,b=_canon(t[0]),_canon(t[1]);return (a,b) if _taxa(a)<=_taxa(b) else (b,a)

@dataclass
class HierarchicalFractionalTree:
    labels:Tuple[str,...]
    activity:Dict[str,float]
    roots:Tuple[str,...]
    split_type:Dict[str,Tuple[int,int,int]]
    children_requests:Dict[str,Tuple[Request,Request]]
    by_flow:Dict[int,Tuple[str,...]]
    leaf_mass:Dict[Tuple[Request,int],float]
    port_taxon_distribution:Dict[Request,np.ndarray]
    unit_taxon_distribution:Dict[str,np.ndarray]
    conditional_long_intensity:Dict[Request,float]
    initial_branch_length:Dict[Request,float]

@dataclass
class Decoded:
    topology:Tree
    newick:str
    root_unit:str
    child_unit:Dict[Request,str]
    leaf_child:Dict[Request,int]
    selected_units:Tuple[str,...]
    proportion:float
    secondary_score:float
    direction:str

@dataclass
class Mixture:
    components:List[Decoded]
    extracted_mass:float
    residual_root_mass:float
    root_down_components:int
    leaves_up_components:int
    fractional_tree:HierarchicalFractionalTree
    notes:str=''


def build_hierarchical_fractional_tree(model,result):
    x=result.x;L=model.L;labels=tuple(model.lib.get('taxa',[f'taxon_{i+1}' for i in range(L)]))
    act={u:max(float(x[model.act[u]]),0.0) for u in model.splits};roots=tuple(sorted(model.roots))
    st={};cr={};by=defaultdict(list);leaf={};pd={};ud={};h={};initlen={}
    for u in model.splits:
        rec=model.units[u];m=int(rec['subtree_size']);a=int(rec['left_flow']);b=int(rec['right_flow']);st[u]=(m,a,b);by[m].append(u);cr[u]=((u,'left'),(u,'right'))
        whole=np.zeros(L)
        for p,flow in [('left',a),('right',b)]:
            mm,s,_=model._req_slot(u,p);v=np.array([max(float(x[model.tmem[mm,0,s,t]]),0.0) for t in range(L)])
            whole+=v
            den=flow*act[u]
            q=v/den if den>1e-14 else np.zeros(L);z=q.sum();pd[u,p]=q/z if z>0 else q
            if flow==1:
                for t in range(L):leaf[(u,p),t]=float(v[t])
            e=getattr(model,'ell',{});raw=max(float(x[e[u,p]]),0.0) if (u,p) in e else 0.0
            h[u,p]=0 if act[u]<=1e-14 else min(max(raw/act[u],0.0),1.0)
            # Spectral branch modes are retained by the LP, but in the pruned/root-only
            # model most internal mode variables are not likelihood-coupled.  Do not use
            # their arbitrary LP values to initialize final fixed-tree branch optimization.
        den=m*act[u];q=whole/den if den>1e-14 else np.zeros(L);z=q.sum();ud[u]=q/z if z>0 else q
    return HierarchicalFractionalTree(labels,act,roots,st,cr,{k:tuple(sorted(v)) for k,v in by.items()},leaf,pd,ud,h,initlen)

class Residual:
    def __init__(self,F):self.F=F;self.act=dict(F.activity);self.leaf=dict(F.leaf_mass)
    def root_total(self):return sum(self.act[r] for r in self.F.roots)
    def subtract(self,z:Decoded,tol=1e-10):
        th=z.proportion
        for u in z.selected_units:
            self.act[u]-=th
            if self.act[u]<-tol:raise RuntimeError('negative activity')
            self.act[u]=max(self.act[u],0)
        for r,t in z.leaf_child.items():
            self.leaf[r,t]-=th
            if self.leaf[r,t]<-tol:raise RuntimeError('negative leaf membership')
            self.leaf[r,t]=max(self.leaf[r,t],0)

def _overlap(a,b):return float(np.minimum(a,b).sum()) if len(a) else 0.0

def _nw(t,lab):
    if isinstance(t,int):return str(lab[t])
    return f'({_nw(t[0],lab)},{_nw(t[1],lab)})'

def _topology(F,root,cu,lc):
    def rec(u):
        vals=[]
        for r in F.children_requests[u]:
            flow=F.split_type[u][1] if r[1]=='left' else F.split_type[u][2]
            if flow==1:vals.append(lc[r])
            else:vals.append(rec(cu[r]))
        return _canon((vals[0],vals[1]))
    return rec(root)

@dataclass
class _TD:
    root:str;pending:Tuple[Request,...];cu:Dict[Request,str];lc:Dict[Request,int];used_units:frozenset;used_taxa:frozenset;theta:float;sec:float

def root_down(F,R,beam_width=24,child_branching=6,leaf_branching=4,max_expansions=100000):
    beam=[]
    for r in F.roots:
        a=R.act[r]
        if a>1e-12:beam.append(_TD(r,F.children_requests[r],{}, {},frozenset({r}),frozenset(),a,math.log(max(a,1e-300))))
    beam.sort(key=lambda s:(-s.theta,-s.sec));beam=beam[:beam_width];best=None;exp=0
    while beam and exp<max_expansions:
        nxt=[]
        for s in beam:
            if not s.pending:
                if len(s.used_taxa)==len(F.labels) and (best is None or (s.theta,s.sec)>(best.theta,best.sec)):best=s
                continue
            r=max(s.pending,key=lambda q:F.split_type[q[0]][1] if q[1]=='left' else F.split_type[q[0]][2]);rest=list(s.pending);rest.remove(r)
            u,p=r;m,a,b=F.split_type[u];flow=a if p=='left' else b
            if flow==1:
                opts=[]
                for t in range(len(F.labels)):
                    if t in s.used_taxa:continue
                    mass=R.leaf.get((r,t),0.0)
                    if mass>1e-12:opts.append((mass,t))
                opts.sort(reverse=True)
                for mass,t in opts[:leaf_branching]:
                    lc=dict(s.lc);lc[r]=t;ut=set(s.used_taxa);ut.add(t)
                    nxt.append(_TD(s.root,tuple(rest),s.cu,lc,s.used_units,frozenset(ut),min(s.theta,mass),s.sec+math.log(max(mass/max(R.act[u],1e-300),1e-300))))
                    exp+=1
            else:
                opts=[]
                for c in F.by_flow.get(flow,()):
                    if c in s.used_units:continue
                    mass=R.act[c]
                    if mass<=1e-12:continue
                    sim=_overlap(F.port_taxon_distribution[r],F.unit_taxon_distribution[c])
                    opts.append((mass,sim,c))
                opts.sort(key=lambda z:(-z[0],-z[1],z[2]))
                for mass,sim,c in opts[:child_branching]:
                    cu=dict(s.cu);cu[r]=c;uu=set(s.used_units);uu.add(c);pend=list(rest);pend.extend(F.children_requests[c])
                    nxt.append(_TD(s.root,tuple(pend),cu,s.lc,frozenset(uu),s.used_taxa,min(s.theta,mass),s.sec+math.log(max(mass,1e-300))+sim))
                    exp+=1
            if exp>=max_expansions:break
        if not nxt:break
        nxt.sort(key=lambda s:(-s.theta,-s.sec,len(s.pending)));beam=nxt[:beam_width]
    if best is None:
        for s in beam:
            if not s.pending and len(s.used_taxa)==len(F.labels) and (best is None or (s.theta,s.sec)>(best.theta,best.sec)):best=s
    if best is None:return None
    t=_topology(F,best.root,best.cu,best.lc)
    return Decoded(t,_nw(t,F.labels)+';',best.root,best.cu,best.lc,tuple(sorted(best.used_units)),best.theta,best.sec,'root-down')

@dataclass(frozen=True)
class _C:
    root:Union[str,int];tree:Tree;flow:int
@dataclass
class _BU:
    clusters:Tuple[_C,...];cu:Dict[Request,str];lc:Dict[Request,int];used:frozenset;theta:float;sec:float;root:Optional[str]

def leaves_up(F,R,beam_width=20,merge_branching=48,max_expansions=100000):
    L=len(F.labels);beam=[_BU(tuple(_C(t,t,1) for t in range(L)),{}, {},frozenset(),1.0,0.0,None)];best=None;exp=0
    while beam and exp<max_expansions:
        nxt=[]
        for s in beam:
            if len(s.clusters)==1:
                c=s.clusters[0]
                if isinstance(c.root,str) and c.root in F.roots:
                    th=min(s.theta,R.act[c.root]);z=_BU(s.clusters,s.cu,s.lc,s.used,th,s.sec+math.log(max(R.act[c.root],1e-300)),c.root)
                    if best is None or (z.theta,z.sec)>(best.theta,best.sec):best=z
                continue
            C=s.clusters;by=defaultdict(list)
            for i,c in enumerate(C):by[c.flow].append((i,c))
            opts=[]
            for u,(m,a,b) in F.split_type.items():
                if u in s.used or R.act[u]<=1e-12:continue
                if m==L and len(C)!=2:continue
                if not by[a] or not by[b]:continue
                for i,A in by[a]:
                    for j,B in by[b]:
                        if i==j:continue
                        ors=[(A,B)] if a!=b else [(A,B),(B,A)]
                        for LC,RC in ors:
                            lreq,rreq=F.children_requests[u];th=min(s.theta,R.act[u]);sec=s.sec+math.log(max(R.act[u],1e-300));ok=True;cuadd={};lcadd={}
                            for req,cl in [(lreq,LC),(rreq,RC)]:
                                flow=cl.flow
                                if flow==1:
                                    t=int(cl.root);mass=R.leaf.get((req,t),0.0)
                                    if mass<=1e-12:ok=False;break
                                    th=min(th,mass);sec+=math.log(max(mass/max(R.act[u],1e-300),1e-300));lcadd[req]=t
                                else:
                                    child=str(cl.root);cuadd[req]=child;sec+=_overlap(F.port_taxon_distribution[req],F.unit_taxon_distribution[child])
                            if ok and th>1e-12:opts.append((th,sec,u,i,j,LC,RC,cuadd,lcadd))
            opts.sort(key=lambda z:(-z[0],-z[1]))
            for th,sec,u,i,j,LC,RC,ca,la in opts[:merge_branching]:
                cu=dict(s.cu);cu.update(ca);lc=dict(s.lc);lc.update(la);tr=_canon((LC.tree,RC.tree));nc=_C(u,tr,LC.flow+RC.flow)
                cls=[c for q,c in enumerate(C) if q not in (i,j)]+[nc];cls.sort(key=lambda c:(c.flow,str(c.root)))
                nxt.append(_BU(tuple(cls),cu,lc,frozenset(set(s.used)|{u}),th,sec,u if nc.flow==L else s.root));exp+=1
                if exp>=max_expansions:break
            if exp>=max_expansions:break
        if not nxt:break
        nxt.sort(key=lambda s:(-s.theta,-s.sec,len(s.clusters)));beam=nxt[:beam_width]
    if best is None:
        for s in beam:
            if len(s.clusters)==1 and isinstance(s.clusters[0].root,str) and s.clusters[0].root in F.roots:
                th=min(s.theta,R.act[s.clusters[0].root]);z=_BU(s.clusters,s.cu,s.lc,s.used,th,s.sec,s.clusters[0].root)
                if best is None or (z.theta,z.sec)>(best.theta,best.sec):best=z
    if best is None:return None
    t=best.clusters[0].tree
    return Decoded(t,_nw(t,F.labels)+';',str(best.root),best.cu,best.lc,tuple(sorted(best.used)),best.theta,best.sec,'leaves-up')


def decompose(model,result,max_components=20,root_beam=24,leaf_beam=20,child_branching=6,leaf_branching=4,merge_branching=48,min_mass=1e-4,use_root_down=True,use_leaves_up=True):
    F=build_hierarchical_fractional_tree(model,result);R=Residual(F);out=[];rd=bu=0
    for _ in range(max_components):
        cand=[]
        if use_root_down:
            z=root_down(F,R,root_beam,child_branching,leaf_branching)
            if z:cand.append(z)
        if use_leaves_up:
            z=leaves_up(F,R,leaf_beam,merge_branching)
            if z:cand.append(z)
        if not cand:break
        z=max(cand,key=lambda q:(q.proportion,q.secondary_score))
        if z.proportion<min_mass:break
        R.subtract(z);out.append(z);rd+=z.direction=='root-down';bu+=z.direction=='leaves-up'
        if R.root_total()<min_mass:break
    return Mixture(out,sum(z.proportion for z in out),R.root_total(),rd,bu,F,
        notes=('Internal requests of flow k>1 match any unused count-k subtree unit; unit activity is frequency support and taxon overlap is only a tie-breaker. '
               'Named taxon membership is enforced only at flow 1. Component proportions are monotone extractable frequencies; root-down/leaves-up search avoids whole-tree-space traversal.'))

# ---------------------------------------------------------------------------
# Non-destructive top-K candidate ranking for final integral-tree rescoring.
# ---------------------------------------------------------------------------

def root_down_topk(F,R,topk=10,beam_width=128,child_branching=12,leaf_branching=8,max_expansions=500000):
    beam=[];completed={};exp=0
    for r in F.roots:
        a=R.act[r]
        if a>1e-12:beam.append(_TD(r,F.children_requests[r],{}, {},frozenset({r}),frozenset(),a,math.log(max(a,1e-300))))
    beam.sort(key=lambda s:(-s.theta,-s.sec));beam=beam[:beam_width]
    while beam and exp<max_expansions:
        nxt=[]
        for s in beam:
            if not s.pending:
                if len(s.used_taxa)==len(F.labels):
                    t=_topology(F,s.root,s.cu,s.lc);k=repr(t)
                    z=Decoded(t,_nw(t,F.labels)+';',s.root,s.cu,s.lc,tuple(sorted(s.used_units)),s.theta,s.sec,'root-down')
                    old=completed.get(k)
                    if old is None or (z.proportion,z.secondary_score)>(old.proportion,old.secondary_score):completed[k]=z
                continue
            r=max(s.pending,key=lambda q:F.split_type[q[0]][1] if q[1]=='left' else F.split_type[q[0]][2]);rest=list(s.pending);rest.remove(r)
            u,p=r;m,a,b=F.split_type[u];flow=a if p=='left' else b
            if flow==1:
                opts=[]
                for t in range(len(F.labels)):
                    if t in s.used_taxa:continue
                    mass=R.leaf.get((r,t),0.0)
                    if mass>1e-12:opts.append((mass,t))
                opts.sort(reverse=True)
                for mass,t in opts[:leaf_branching]:
                    lc=dict(s.lc);lc[r]=t;ut=set(s.used_taxa);ut.add(t)
                    nxt.append(_TD(s.root,tuple(rest),s.cu,lc,s.used_units,frozenset(ut),min(s.theta,mass),s.sec+math.log(max(mass/max(R.act[u],1e-300),1e-300))));exp+=1
            else:
                opts=[]
                for c in F.by_flow.get(flow,()):
                    if c in s.used_units:continue
                    mass=R.act[c]
                    if mass<=1e-12:continue
                    sim=_overlap(F.port_taxon_distribution[r],F.unit_taxon_distribution[c]);opts.append((mass,sim,c))
                opts.sort(key=lambda z:(-z[0],-z[1],z[2]))
                for mass,sim,c in opts[:child_branching]:
                    cu=dict(s.cu);cu[r]=c;uu=set(s.used_units);uu.add(c);pend=list(rest);pend.extend(F.children_requests[c])
                    nxt.append(_TD(s.root,tuple(pend),cu,s.lc,frozenset(uu),s.used_taxa,min(s.theta,mass),s.sec+math.log(max(mass,1e-300))+sim));exp+=1
            if exp>=max_expansions:break
        if not nxt:break
        nxt.sort(key=lambda s:(-s.theta,-s.sec,len(s.pending)));beam=nxt[:beam_width]
        if len(completed)>=topk and beam and beam[0].theta < sorted((z.proportion for z in completed.values()),reverse=True)[topk-1]-1e-12:break
    return sorted(completed.values(),key=lambda z:(-z.proportion,-z.secondary_score,z.newick))[:topk]


def leaves_up_topk(F,R,topk=10,beam_width=32,merge_branching=96,max_expansions=500000):
    L=len(F.labels);beam=[_BU(tuple(_C(t,t,1) for t in range(L)),{}, {},frozenset(),1.0,0.0,None)];completed={};exp=0
    while beam and exp<max_expansions:
        nxt=[]
        for s in beam:
            if len(s.clusters)==1:
                c=s.clusters[0]
                if isinstance(c.root,str) and c.root in F.roots:
                    th=min(s.theta,R.act[c.root]);t=c.tree;k=repr(t)
                    z=Decoded(t,_nw(t,F.labels)+';',str(c.root),s.cu,s.lc,tuple(sorted(s.used)),th,s.sec+math.log(max(R.act[c.root],1e-300)),'leaves-up')
                    old=completed.get(k)
                    if old is None or (z.proportion,z.secondary_score)>(old.proportion,old.secondary_score):completed[k]=z
                continue
            C=s.clusters;by=defaultdict(list)
            for i,c in enumerate(C):by[c.flow].append((i,c))
            opts=[]
            for u,(m,a,b) in F.split_type.items():
                if u in s.used or R.act[u]<=1e-12:continue
                if m==L and len(C)!=2:continue
                if not by[a] or not by[b]:continue
                for i,A in by[a]:
                    for j,B in by[b]:
                        if i==j:continue
                        ors=[(A,B)] if a!=b else [(A,B),(B,A)]
                        for LC,RC in ors:
                            lreq,rreq=F.children_requests[u];th=min(s.theta,R.act[u]);sec=s.sec+math.log(max(R.act[u],1e-300));ok=True;cuadd={};lcadd={}
                            for req,cl in ((lreq,LC),(rreq,RC)):
                                if cl.flow==1:
                                    t=int(cl.root);mass=R.leaf.get((req,t),0.0)
                                    if mass<=1e-12:ok=False;break
                                    th=min(th,mass);sec+=math.log(max(mass/max(R.act[u],1e-300),1e-300));lcadd[req]=t
                                else:
                                    child=str(cl.root);cuadd[req]=child;sec+=_overlap(F.port_taxon_distribution[req],F.unit_taxon_distribution[child])
                            if ok and th>1e-12:opts.append((th,sec,u,i,j,LC,RC,cuadd,lcadd))
            opts.sort(key=lambda z:(-z[0],-z[1]))
            for th,sec,u,i,j,LC,RC,ca,la in opts[:merge_branching]:
                cu=dict(s.cu);cu.update(ca);lc=dict(s.lc);lc.update(la);tr=_canon((LC.tree,RC.tree));nc=_C(u,tr,LC.flow+RC.flow)
                cls=[c for q,c in enumerate(C) if q not in (i,j)]+[nc];cls.sort(key=lambda c:(c.flow,str(c.root)))
                nxt.append(_BU(tuple(cls),cu,lc,frozenset(set(s.used)|{u}),th,sec,u if nc.flow==L else s.root));exp+=1
                if exp>=max_expansions:break
            if exp>=max_expansions:break
        if not nxt:break
        nxt.sort(key=lambda s:(-s.theta,-s.sec,len(s.clusters)));beam=nxt[:beam_width]
    return sorted(completed.values(),key=lambda z:(-z.proportion,-z.secondary_score,z.newick))[:topk]


def rank_integral_tree_candidates(model,result,topk=10,root_beam=128,leaf_beam=32,
                                  child_branching=12,leaf_branching=8,merge_branching=96,
                                  max_expansions=500000):
    """Return a non-destructive top-K list of integral structural matches.

    Primary score is the bottleneck support (`Decoded.proportion`); secondary score is
    accumulated local support/taxon overlap.  No mass is removed and no NNI/SPR walk is used.
    """
    F=build_hierarchical_fractional_tree(model,result);R=Residual(F);allc={}
    for z in root_down_topk(F,R,max(topk,topk//2+1),root_beam,child_branching,leaf_branching,max_expansions):
        allc[repr(z.topology)]=z
    for z in leaves_up_topk(F,R,max(topk,topk//2+1),leaf_beam,merge_branching,max_expansions):
        k=repr(z.topology);old=allc.get(k)
        if old is None or (z.proportion,z.secondary_score)>(old.proportion,old.secondary_score):allc[k]=z
    ranked=sorted(allc.values(),key=lambda z:(-z.proportion,-z.secondary_score,z.newick))[:topk]
    return ranked,F
