#!/usr/bin/env python3
"""Fractional-tree reconstruction and hierarchical integral-tree decoding.

The projected likelihood LP stores shared split activities and named-taxon marginals
at every fixed-flow lower request, but it intentionally does not retain a dense
request->child likelihood-routing vector per alignment column.  This module converts
those marginals into an explicit *fractional biological tree* before decoding.

For each descendant flow k<L:
  * rows are lower requests of flow k, with mass equal to their parent split activity;
  * columns are count-k child split units (or named taxa for k=1), with mass equal to
    child activity (1 for each taxon);
  * a taxon-compatible fractional transport matrix couples rows to columns.

The coupling is obtained by entropy-regularized matrix scaling (Sinkhorn), not an LP.
Its preference kernel is the overlap of the row and column taxon-membership
probability vectors.  Row/column sums therefore reproduce the solved fractional
unit frequencies exactly (up to numerical tolerance), while the coupling picks the
most taxon-compatible realization of the projected marginals.

An integral tree is then a one-to-one choice of one child output for every selected
lower request.  Root-down and leaves-up beam decoders operate directly on this DAG.
For a partial tree, the largest extractable component mass is the minimum residual
root/edge mass used so far.  Adding choices can only decrease this value, so it is a
monotone branch-and-bound score unlike NNI/SPR similarity.

After a complete tree is chosen with mass theta, subtract theta from its root and
all selected fractional DAG edges.  This preserves row/column conservation and gives
an actual residual fractional tree.  Repeating yields multiple component trees whose
reported proportions are literal extractable masses from this reconstructed DAG.

Branch length is required in the likelihood LP but excluded from topology matching.
For each lower arm e=(u,p): frequency=act[u], conditional long intensity=ell[e]/act[u].
The latter can be attached after topology recovery or ignored in favor of fast
integral-tree branch-length optimization.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union
import json, math, random
import numpy as np

Tree = Union[int, Tuple["Tree", "Tree"]]
Request = Tuple[str, str]
Output = str  # split-unit id or "taxon:<index>"


def _canon(t: Tree) -> Tree:
    if isinstance(t, int): return int(t)
    a,b=_canon(t[0]),_canon(t[1])
    return (a,b) if _taxa(a)<=_taxa(b) else (b,a)

def _taxa(t: Tree) -> Tuple[int,...]:
    if isinstance(t,int): return (t,)
    return tuple(sorted(_taxa(t[0])+_taxa(t[1])))


@dataclass
class FractionalDAGEdge:
    request: Request
    child: Output
    flow: int
    mass: float
    conditional: float
    taxon_similarity: float


@dataclass
class FractionalTreeDAG:
    taxa_labels: Tuple[str,...]
    root_mass: Dict[str,float]
    request_mass: Dict[Request,float]
    output_mass: Dict[Output,float]
    edge_mass: Dict[Tuple[Request,Output],float]
    edge_similarity: Dict[Tuple[Request,Output],float]
    request_flow: Dict[Request,int]
    output_flow: Dict[Output,int]
    split_children: Dict[str,Tuple[Request,Request]]
    split_type: Dict[str,Tuple[int,int,int]]
    arm_long_intensity: Dict[Request,float]
    sinkhorn_error: Dict[int,float] = field(default_factory=dict)

    def adjacency(self, tol=1e-12):
        out=defaultdict(list)
        for (r,c),m in self.edge_mass.items():
            if m>tol:
                rm=max(self.request_mass.get(r,0.0),1e-300)
                out[r].append(FractionalDAGEdge(r,c,self.request_flow[r],m,m/rm,self.edge_similarity.get((r,c),0.0)))
        for r in out:
            out[r].sort(key=lambda e:(-e.mass,-e.taxon_similarity,e.child))
        return out

    def reverse_adjacency(self, tol=1e-12):
        out=defaultdict(list)
        for egs in self.adjacency(tol).values():
            for e in egs: out[e.child].append(e)
        return out

    def to_dict(self):
        return {
            'taxa':list(self.taxa_labels),
            'root_mass':self.root_mass,
            'request_mass':{f'{u}:{p}':v for (u,p),v in self.request_mass.items()},
            'output_mass':self.output_mass,
            'edges':[{'parent_unit':r[0],'parent_port':r[1],'child':c,'flow':self.request_flow[r],
                      'mass':m,'conditional':m/max(self.request_mass[r],1e-300),
                      'taxon_similarity':self.edge_similarity.get((r,c),0.0)}
                     for (r,c),m in self.edge_mass.items() if m>1e-12],
            'arm_long_intensity':{f'{u}:{p}':v for (u,p),v in self.arm_long_intensity.items()},
            'sinkhorn_error':self.sinkhorn_error,
        }

    def write_json(self,path):
        with open(path,'w',encoding='utf-8') as f: json.dump(self.to_dict(),f,indent=2)


@dataclass
class DecodedTree:
    topology: Tree
    newick: str
    root_unit: str
    child_of_request: Dict[Request,Output]
    proportion: float
    log_match_score: float
    direction: str


@dataclass
class DAGMixtureResult:
    components: List[DecodedTree]
    extracted_mass: float
    residual_root_mass: float
    fractional_dag: FractionalTreeDAG
    root_down_components: int
    leaves_up_components: int
    notes: str=''


def _taxon_distribution_from_request(model,result,u,p,flow,activity):
    L=model.L
    if activity<=1e-14: return np.zeros(L)
    m,s,_=model._req_slot(u,p)
    v=np.array([max(float(result.x[model.tmem[m,0,s,t]]),0.0) for t in range(L)],float)
    denom=flow*activity
    if denom>0: v/=denom
    z=v.sum()
    return v/z if z>0 else np.full(L,1.0/L)


def _taxon_distribution_from_unit(model,result,u,flow,activity):
    L=model.L
    if activity<=1e-14: return np.zeros(L)
    v=np.zeros(L)
    for p in ('left','right'):
        m,s,_=model._req_slot(u,p)
        v += np.array([max(float(result.x[model.tmem[m,0,s,t]]),0.0) for t in range(L)],float)
    denom=flow*activity
    if denom>0: v/=denom
    z=v.sum()
    return v/z if z>0 else np.full(L,1.0/L)


def _sinkhorn_transport(row_mass,col_mass,similarity,tau=0.08,max_iter=1000,tol=1e-10):
    r=np.asarray(row_mass,float);c=np.asarray(col_mass,float);S=np.asarray(similarity,float)
    sr=r.sum();sc=c.sum()
    if sr<=0 or sc<=0: return np.zeros_like(S),0.0
    if abs(sr-sc)>1e-8*max(sr,sc,1.0):
        # Solver tolerances / projected marginals should make these equal.  A tiny
        # rescaling of columns preserves their relative solved frequencies.
        c=c*(sr/sc)
    # Stabilized positive kernel.  Similarity is in [0,1].
    K=np.exp((S-S.max())/max(tau,1e-6)) + 1e-15
    u=np.ones(len(r));v=np.ones(len(c))
    for it in range(max_iter):
        Kv=K@v; u=np.divide(r,Kv,out=np.zeros_like(r),where=Kv>0)
        KTu=K.T@u; v=np.divide(c,KTu,out=np.zeros_like(c),where=KTu>0)
        if it%10==0 or it==max_iter-1:
            Y=(u[:,None]*K)*v[None,:]
            err=max(np.max(np.abs(Y.sum(axis=1)-r)) if len(r) else 0,
                    np.max(np.abs(Y.sum(axis=0)-c)) if len(c) else 0)
            if err<tol: break
    Y=(u[:,None]*K)*v[None,:]
    err=max(np.max(np.abs(Y.sum(axis=1)-r)) if len(r) else 0,
            np.max(np.abs(Y.sum(axis=0)-c)) if len(c) else 0)
    return Y,float(err)


def build_fractional_tree_dag(model,result,*,tau=0.08,activity_tol=1e-10) -> FractionalTreeDAG:
    x=result.x;L=int(model.L);labels=tuple(model.lib.get('taxa',[f'taxon_{i+1}' for i in range(L)]))
    act={u:max(float(x[model.act[u]]),0.0) for u in model.splits}
    root_mass={r:act[r] for r in model.roots if act[r]>activity_tol}
    request_mass={};request_flow={};split_children={};split_type={};armh={}
    by_flow_req=defaultdict(list)
    for u in model.splits:
        rec=model.units[u];m=int(rec['subtree_size']);a=int(rec['left_flow']);b=int(rec['right_flow'])
        split_type[u]=(m,a,b);split_children[u]=((u,'left'),(u,'right'))
        for p,flow in (('left',a),('right',b)):
            r=(u,p);request_mass[r]=act[u];request_flow[r]=flow;by_flow_req[flow].append(r)
            e=getattr(model,'ell',{});raw=max(float(x[e[u,p]]),0.0) if (u,p) in e else 0.0
            armh[r]=0.0 if act[u]<=1e-14 else min(max(raw/act[u],0.0),1.0)

    output_mass={};output_flow={};by_flow_out=defaultdict(list)
    for t in range(L):
        o=f'taxon:{t}';output_mass[o]=1.0;output_flow[o]=1;by_flow_out[1].append(o)
    for u in model.splits:
        m=int(model.units[u]['subtree_size'])
        if m<L and act[u]>activity_tol:
            output_mass[u]=act[u];output_flow[u]=m;by_flow_out[m].append(u)

    edge_mass={};edge_similarity={};errors={}
    for k in range(1,L):
        rows=[r for r in by_flow_req.get(k,[]) if request_mass[r]>activity_tol]
        cols=[o for o in by_flow_out.get(k,[]) if output_mass[o]>activity_tol]
        if not rows and not cols: continue
        if not rows or not cols:
            raise RuntimeError(f'flow {k}: unbalanced nonempty request/output sets')
        rm=np.array([request_mass[r] for r in rows],float);cm=np.array([output_mass[o] for o in cols],float)
        RP=np.vstack([_taxon_distribution_from_request(model,result,r[0],r[1],k,request_mass[r]) for r in rows])
        CP=[]
        for o in cols:
            if o.startswith('taxon:'):
                q=np.zeros(L);q[int(o.split(':')[1])]=1.0
            else:q=_taxon_distribution_from_unit(model,result,o,k,output_mass[o])
            CP.append(q)
        CP=np.vstack(CP)
        # overlap coefficient sum min(p,q), vectorized via L small (<=100 target)
        S=np.zeros((len(rows),len(cols)))
        for i in range(len(rows)):
            S[i,:]=np.minimum(RP[i][None,:],CP).sum(axis=1)
        Y,err=_sinkhorn_transport(rm,cm,S,tau=tau);errors[k]=err
        for i,r in enumerate(rows):
            for j,o in enumerate(cols):
                y=float(Y[i,j])
                if y>1e-14:
                    edge_mass[r,o]=y;edge_similarity[r,o]=float(S[i,j])
    return FractionalTreeDAG(labels,root_mass,request_mass,output_mass,edge_mass,edge_similarity,
                             request_flow,output_flow,split_children,split_type,armh,errors)


class ResidualDAG:
    def __init__(self,dag:FractionalTreeDAG):
        self.dag=dag;self.root=dict(dag.root_mass);self.edge=dict(dag.edge_mass)
    def row_mass(self,r): return sum(m for (rr,_),m in self.edge.items() if rr==r and m>0)
    def adjacency(self,tol=1e-12):
        out=defaultdict(list)
        row=defaultdict(float)
        for (r,c),m in self.edge.items():
            if m>tol: row[r]+=m
        for (r,c),m in self.edge.items():
            if m>tol:
                out[r].append((c,m,m/max(row[r],1e-300),self.dag.edge_similarity.get((r,c),0.0)))
        for r in out: out[r].sort(key=lambda z:(-z[1],-z[3],z[0]))
        return out
    def subtract(self,tr:DecodedTree,tol=1e-10):
        th=tr.proportion;self.root[tr.root_unit]=self.root.get(tr.root_unit,0.0)-th
        if self.root[tr.root_unit]<-tol: raise RuntimeError('negative root residual')
        self.root[tr.root_unit]=max(self.root[tr.root_unit],0.0)
        for r,c in tr.child_of_request.items():
            k=(r,c);self.edge[k]=self.edge.get(k,0.0)-th
            if self.edge[k]<-tol: raise RuntimeError(f'negative edge residual {k} {self.edge[k]}')
            self.edge[k]=max(self.edge[k],0.0)
    def root_total(self):return sum(self.root.values())


def _topology_from_mapping(dag:FractionalTreeDAG,root:str,childmap:Dict[Request,Output]) -> Tree:
    def rec(u):
        vals=[]
        for r in dag.split_children[u]:
            c=childmap[r]
            if c.startswith('taxon:'):vals.append(int(c.split(':')[1]))
            else:vals.append(rec(c))
        return _canon((vals[0],vals[1]))
    return rec(root)

def _newick(t:Tree,labels):
    if isinstance(t,int):return str(labels[t])
    return f'({_newick(t[0],labels)},{_newick(t[1],labels)})'


@dataclass
class _TD:
    root:str
    pending:Tuple[Request,...]
    childmap:Dict[Request,Output]
    used: frozenset
    theta:float
    logscore:float


def decode_root_down_dag(dag:FractionalTreeDAG,res:ResidualDAG,*,beam_width=32,edge_branching=6,max_expansions=100000):
    adj=res.adjacency();beam=[]
    for r,m in res.root.items():
        if m<=1e-12:continue
        pending=tuple(dag.split_children[r]);beam.append(_TD(r,pending,{},frozenset(),m,math.log(max(m,1e-300))))
    beam.sort(key=lambda s:(-s.theta,-s.logscore));beam=beam[:beam_width]
    best=None;exp=0
    while beam and exp<max_expansions:
        nxt=[]
        for st in beam:
            if not st.pending:
                if best is None or (st.theta,st.logscore)>(best.theta,best.logscore):best=st
                continue
            # largest-flow pending request first = literal root-down traversal
            r=max(st.pending,key=lambda q:dag.request_flow[q]);rest=list(st.pending);rest.remove(r)
            for c,m,cond,sim in adj.get(r,[])[:edge_branching]:
                if c in st.used:continue
                th=min(st.theta,m)
                if th<=1e-12:continue
                if best is not None and th+1e-12<best.theta:continue
                cm=dict(st.childmap);cm[r]=c;used=set(st.used);used.add(c);pend=list(rest)
                if not c.startswith('taxon:'):pend.extend(dag.split_children[c])
                nxt.append(_TD(st.root,tuple(pend),cm,frozenset(used),th,st.logscore+math.log(max(cond,1e-300))))
                exp+=1
                if exp>=max_expansions:break
            if exp>=max_expansions:break
        if not nxt:break
        nxt.sort(key=lambda s:(-s.theta,-s.logscore,len(s.pending)));beam=nxt[:beam_width]
    if best is None:
        for st in beam:
            if not st.pending and (best is None or (st.theta,st.logscore)>(best.theta,best.logscore)):best=st
    if best is None:return None
    topo=_topology_from_mapping(dag,best.root,best.childmap)
    return DecodedTree(topo,_newick(topo,dag.taxa_labels)+';',best.root,best.childmap,float(best.theta),float(best.logscore),'root-down')


@dataclass(frozen=True)
class _Cl:
    output:Output
    tree:Tree
    flow:int
@dataclass
class _BU:
    clusters:Tuple[_Cl,...]
    used:frozenset
    childmap:Dict[Request,Output]
    theta:float
    logscore:float
    root:Optional[str]


def decode_leaves_up_dag(dag:FractionalTreeDAG,res:ResidualDAG,*,beam_width=24,merge_branching=48,max_expansions=100000):
    adj=res.adjacency();em={(r,c):(m,cond,sim) for r,ls in adj.items() for c,m,cond,sim in ls}
    L=len(dag.taxa_labels);init=tuple(_Cl(f'taxon:{t}',t,1) for t in range(L));beam=[_BU(init,frozenset(),{},1.0,0.0,None)]
    best=None;exp=0
    units=list(dag.split_type)
    while beam and exp<max_expansions:
        nxt=[]
        for st in beam:
            if len(st.clusters)==1:
                cl=st.clusters[0]
                if cl.flow==L and cl.output in dag.root_mass:
                    root=cl.output;rm=res.root.get(root,0.0);th=min(st.theta,rm)
                    cand=_BU(st.clusters,st.used,st.childmap,th,st.logscore+math.log(max(rm,1e-300)),root)
                    if best is None or (cand.theta,cand.logscore)>(best.theta,best.logscore):best=cand
                continue
            candidates=[];C=st.clusters
            byflow=defaultdict(list)
            for i,c in enumerate(C):byflow[c.flow].append((i,c))
            for u in units:
                if u in st.used:continue
                m,a,b=dag.split_type[u]
                if m>L or not byflow.get(a) or not byflow.get(b):continue
                # size-L units are roots; they should only be used for the final merge
                if m==L and len(C)!=2:continue
                if m<L and u in dag.root_mass:continue
                lreq,rreq=dag.split_children[u]
                for i,A in byflow[a]:
                    for j,B in byflow[b]:
                        if i==j:continue
                        orientations=[(A,B)]
                        if a==b:orientations=[(A,B),(B,A)]
                        for LC,RC in orientations:
                            z1=em.get((lreq,LC.output));z2=em.get((rreq,RC.output))
                            if z1 is None or z2 is None:continue
                            th=min(st.theta,z1[0],z2[0]);log=st.logscore+math.log(max(z1[1],1e-300))+math.log(max(z2[1],1e-300))
                            if m==L:
                                rm=res.root.get(u,0.0)
                                if rm<=1e-12:continue
                                th=min(th,rm);log+=math.log(max(rm,1e-300))
                            if th<=1e-12 or (best is not None and th+1e-12<best.theta):continue
                            candidates.append((th,log,u,i,j,LC,RC))
            candidates.sort(key=lambda z:(-z[0],-z[1]))
            for th,log,u,i,j,LC,RC in candidates[:merge_branching]:
                cm=dict(st.childmap);lreq,rreq=dag.split_children[u];cm[lreq]=LC.output;cm[rreq]=RC.output
                tr=_canon((LC.tree,RC.tree));new=_Cl(u,tr,LC.flow+RC.flow)
                cls=[c for q,c in enumerate(C) if q not in (i,j)]+[new];cls.sort(key=lambda c:(c.flow,c.output))
                nxt.append(_BU(tuple(cls),frozenset(set(st.used)|{u}),cm,th,log,u if new.flow==L else st.root))
                exp+=1
                if exp>=max_expansions:break
            if exp>=max_expansions:break
        if not nxt:break
        nxt.sort(key=lambda s:(-s.theta,-s.logscore,len(s.clusters)));beam=nxt[:beam_width]
    if best is None:
        for st in beam:
            if len(st.clusters)==1 and st.clusters[0].flow==L and st.clusters[0].output in dag.root_mass:
                root=st.clusters[0].output;rm=res.root.get(root,0.0);th=min(st.theta,rm)
                z=_BU(st.clusters,st.used,st.childmap,th,st.logscore+math.log(max(rm,1e-300)),root)
                if best is None or (z.theta,z.logscore)>(best.theta,best.logscore):best=z
    if best is None:return None
    topo=best.clusters[0].tree
    return DecodedTree(topo,_newick(topo,dag.taxa_labels)+';',str(best.root),best.childmap,float(best.theta),float(best.logscore),'leaves-up')


def decompose_fractional_dag(model,result,*,max_components=20,tau=0.08,root_beam=32,leaf_beam=24,
                             edge_branching=6,merge_branching=48,min_mass=1e-4,use_root_down=True,use_leaves_up=True):
    dag=build_fractional_tree_dag(model,result,tau=tau);res=ResidualDAG(dag);comps=[];rd=bu=0
    for _ in range(max_components):
        cand=[]
        if use_root_down:
            z=decode_root_down_dag(dag,res,beam_width=root_beam,edge_branching=edge_branching)
            if z:cand.append(z)
        if use_leaves_up:
            z=decode_leaves_up_dag(dag,res,beam_width=leaf_beam,merge_branching=merge_branching)
            if z:cand.append(z)
        if not cand:break
        best=max(cand,key=lambda z:(z.proportion,z.log_match_score))
        if best.proportion<min_mass:break
        res.subtract(best);comps.append(best)
        rd += best.direction=='root-down';bu += best.direction=='leaves-up'
        if res.root_total()<min_mass:break
    return DAGMixtureResult(comps,float(sum(c.proportion for c in comps)),float(res.root_total()),dag,rd,bu,
        notes=('The fractional DAG is a taxon-compatible Sinkhorn coupling of solved fixed-flow request and child-unit marginals. '
               'Tree proportions are literal masses subtracted from root and selected DAG edges. Root-down and leaves-up '
               'beam search use a monotone extractable-mass bound and do not traverse NNI/SPR space. Branch length remains '
               'in the likelihood LP but is excluded from topology matching via frequency/intensity separation.'))


# ---------------------------------------------------------------------------
# Preferred fractional-tree reconstruction: disaggregate the solved physical
# Beneš switch activities directly, without a transport LP.
# ---------------------------------------------------------------------------

def _switch_transport_from_solution(x0,x1,y0,y1,s):
    """Return a feasible 2x2 input->output flow matrix for one solved switch.

    Rows are inputs (0,1), columns outputs (0,1).  Row/column sums are exactly
    x/y.  The activity hull leaves one balanced straight/cross degree of freedom;
    choose it so total crossed mass is as close as possible to s*(x0+x1).
    This agrees with the integral endpoints s=0 (straight) and s=1 (cross).
    """
    x0=max(float(x0),0.0);x1=max(float(x1),0.0);y0=max(float(y0),0.0);y1=max(float(y1),0.0);s=min(max(float(s),0.0),1.0)
    # Correct tiny conservation error symmetrically through y1.
    total=x0+x1
    if abs((y0+y1)-total)>1e-8*max(1.0,total):
        scale=total/max(y0+y1,1e-300);y0*=scale;y1*=scale
    p=x0-y0  # c0-c1
    bp=max(p,0.0);bn=max(-p,0.0)
    qmax=max(0.0,min(x0-bp,y1-bp,x1-bn,y0-bn))
    desired=s*total
    q=min(max((desired-abs(p))/2.0,0.0),qmax)
    c0=bp+q;c1=bn+q
    F=np.array([[x0-c0,c0],[c1,x1-c1]],float)
    F[F<0]=0
    return F


def build_fractional_tree_dag_physical(model,result,*,activity_tol=1e-12) -> FractionalTreeDAG:
    """Collapse solved physical Beneš switches into biological request->child edges.

    This is the default 'actual fractional tree' representation because it uses
    only topology variables that were already in the LP.  No optimization is
    performed during collapse.  Named-taxon boundary marginals are used only to
    annotate each collapsed edge with a taxon-similarity score.
    """
    x=result.x;L=int(model.L);labels=tuple(model.lib.get('taxa',[f'taxon_{i+1}' for i in range(L)]))
    act={u:max(float(x[model.act[u]]),0.0) for u in model.splits}
    root_mass={r:act[r] for r in model.roots if act[r]>activity_tol}
    request_mass={};request_flow={};split_children={};split_type={};armh={}
    request_taxdist={}
    for u in model.splits:
        rec=model.units[u];m=int(rec['subtree_size']);a=int(rec['left_flow']);b=int(rec['right_flow'])
        split_type[u]=(m,a,b);split_children[u]=((u,'left'),(u,'right'))
        for p,flow in (('left',a),('right',b)):
            r=(u,p);request_mass[r]=act[u];request_flow[r]=flow
            request_taxdist[r]=_taxon_distribution_from_request(model,result,u,p,flow,act[u])
            e=getattr(model,'ell',{});raw=max(float(x[e[u,p]]),0.0) if (u,p) in e else 0.0
            armh[r]=0.0 if act[u]<=1e-14 else min(max(raw/act[u],0.0),1.0)
    output_mass={};output_flow={};output_taxdist={}
    # taxon output lookup from physical taxon units
    tax_index={}
    for uid,rec in model.units.items():
        if rec.get('kind')=='taxon':
            label=rec.get('label','')
            try: tax_index[uid]=labels.index(label)
            except ValueError:
                # generator labels are normally leaf_1 ...
                try: tax_index[uid]=int(str(label).split('_')[-1])-1
                except Exception: pass
    for t in range(L):
        o=f'taxon:{t}';output_mass[o]=1.0;output_flow[o]=1;q=np.zeros(L);q[t]=1;output_taxdist[o]=q
    for u in model.splits:
        m=int(model.units[u]['subtree_size'])
        if m<L:
            output_mass[u]=act[u];output_flow[u]=m;output_taxdist[u]=_taxon_distribution_from_unit(model,result,u,m,act[u])

    edge_mass={};edge_similarity={};errors={}
    routers=model.lib['routers']
    for mstr,meta in routers.items():
        m=int(mstr);W=int(meta['padded_width']);dims=list(meta['dimensions']);D=len(dims)
        inputs={int(z['slot']):(z['parent_unit'],z['parent_port']) for z in meta['input_requests']}
        outs={int(z['slot']):z['unit_id'] for z in meta['output_targets']}
        # map real output unit ids to biological output ids used by decoder
        bout={}
        for slot,uid in outs.items():
            if m==1:
                t=tax_index.get(uid)
                if t is None: continue
                bout[slot]=f'taxon:{t}'
            else: bout[slot]=uid
        if D==0:
            if inputs and bout:
                slot=next(iter(inputs));r=inputs[slot];c=bout.get(slot,next(iter(bout.values())))
                mass=request_mass[r];edge_mass[r,c]=mass
                edge_similarity[r,c]=float(np.minimum(request_taxdist[r],output_taxdist[c]).sum())
            errors[m]=0.0;continue
        # source distributions on wires, keyed by biological request id.
        src=[defaultdict(float) for _ in range(W)]
        for slot,r in inputs.items():
            mass=request_mass[r]
            if mass>activity_tol: src[slot][r]=mass
        for st,dim in enumerate(dims):
            nxt=[defaultdict(float) for _ in range(W)];seen=set()
            for w in range(W):
                b=min(w,w^(1<<dim))
                if b in seen:continue
                seen.add(b);o=b^(1<<dim)
                x0=float(x[model.wireact[m,st,b]]);x1=float(x[model.wireact[m,st,o]])
                y0=float(x[model.wireact[m,st+1,b]]);y1=float(x[model.wireact[m,st+1,o]])
                sval=float(x[model.sw[m,st,b]])
                F=_switch_transport_from_solution(x0,x1,y0,y1,sval)
                for inp,wire in enumerate((b,o)):
                    denom=(x0,x1)[inp]
                    if denom<=activity_tol:continue
                    for req,mass in src[wire].items():
                        for outidx,outwire in enumerate((b,o)):
                            q=F[inp,outidx]/denom
                            if q>0:nxt[outwire][req]+=mass*q
            src=nxt
        # final wire source decomposition -> biological outputs
        colerr=0.0
        for slot,c in bout.items():
            tot=0.0
            for r,mass in src[slot].items():
                if mass>1e-14:
                    edge_mass[r,c]=edge_mass.get((r,c),0.0)+float(mass);tot+=mass
                    edge_similarity[r,c]=float(np.minimum(request_taxdist[r],output_taxdist[c]).sum())
            colerr=max(colerr,abs(tot-output_mass.get(c,0.0)))
        # row error too
        for r in inputs.values():
            got=sum(v for (rr,_),v in edge_mass.items() if rr==r)
            colerr=max(colerr,abs(got-request_mass[r]))
        errors[m]=float(colerr)
    return FractionalTreeDAG(labels,root_mass,request_mass,output_mass,edge_mass,edge_similarity,
                             request_flow,output_flow,split_children,split_type,armh,errors)


def decompose_physical_fractional_tree(model,result,*,max_components=20,root_beam=32,leaf_beam=24,
                                       edge_branching=6,merge_branching=48,min_mass=1e-4,
                                       use_root_down=True,use_leaves_up=True):
    dag=build_fractional_tree_dag_physical(model,result);res=ResidualDAG(dag);comps=[];rd=bu=0
    for _ in range(max_components):
        cand=[]
        if use_root_down:
            z=decode_root_down_dag(dag,res,beam_width=root_beam,edge_branching=edge_branching)
            if z:cand.append(z)
        if use_leaves_up:
            z=decode_leaves_up_dag(dag,res,beam_width=leaf_beam,merge_branching=merge_branching)
            if z:cand.append(z)
        if not cand:break
        best=max(cand,key=lambda z:(z.proportion,z.log_match_score))
        if best.proportion<min_mass:break
        res.subtract(best);comps.append(best);rd+=best.direction=='root-down';bu+=best.direction=='leaves-up'
        if res.root_total()<min_mass:break
    return DAGMixtureResult(comps,float(sum(c.proportion for c in comps)),float(res.root_total()),dag,rd,bu,
        notes=('Preferred decoder: fractional biological DAG is obtained by locally disaggregating the solved physical Beneš '
               'switch activities, using the solved switch selector to resolve the one balanced straight/cross degree of freedom. '
               'No optimization and no tree-space traversal is used to create the DAG. Integral trees are then extracted root-down '
               'and/or leaves-up with a monotone residual-mass score.'))
