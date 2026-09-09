#!/usr/bin/env python3
"""Continuous branch-length ML optimization for recovered integral trees.

Uses Felsenstein pruning with an analytic gradient and SciPy L-BFGS-B on log branch
lengths.  DNA models are reversible JC69 or GTR+F from phylo_substitution.DNAmodel.
"""
from __future__ import annotations
from dataclasses import dataclass
from collections import Counter
import math
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

IUPAC={
 'A':(1,0,0,0),'C':(0,1,0,0),'G':(0,0,1,0),'T':(0,0,0,1),
 'U':(0,0,0,1),'R':(1,0,1,0),'Y':(0,1,0,1),'S':(0,1,1,0),'W':(1,0,0,1),
 'K':(0,0,1,1),'M':(1,1,0,0),'B':(0,1,1,1),'D':(1,0,1,1),'H':(1,1,0,1),
 'V':(1,1,1,0),'N':(1,1,1,1),'X':(1,1,1,1),'?':(1,1,1,1),'-':(1,1,1,1),'.':(1,1,1,1)
}

@dataclass
class OptimizedTree:
    log_likelihood: float
    newick: str
    branch_lengths: np.ndarray
    success: bool
    iterations: int
    evaluations: int
    message: str
    min_length_bound: float
    max_length_bound: float


def _quote_label(x):
    s=str(x)
    if all(ch.isalnum() or ch in '_.-' for ch in s):return s
    return "'"+s.replace("'","''")+"'"


def compress_alignment(sequences):
    if not sequences:return [],np.array([],float)
    L=len(sequences);n=len(sequences[0])
    if any(len(s)!=n for s in sequences):raise ValueError('Alignment sequences have unequal lengths')
    cnt=Counter(tuple(sequences[t][j].upper() for t in range(L)) for j in range(n))
    pats=list(cnt);weights=np.asarray([cnt[p] for p in pats],float)
    leaf=[]
    for t in range(L):
        leaf.append(np.asarray([IUPAC.get(p[t],IUPAC['N']) for p in pats],float))
    return leaf,weights


class _TreeLayout:
    def __init__(self,topology,labels,initial_by_request=None,decoded=None,frac=None,short_time=None,long_ratio=2.0):
        self.labels=labels;self.children={};self.edge=[];self.edge_init=[];self.leaf_node={};self._next=len(labels)
        # If decoded mapping is available, build from physical units so each edge gets the LP intensity initialization.
        if decoded is not None and frac is not None and short_time is not None:
            def rec_unit(u):
                nid=self._next;self._next+=1;chs=[]
                m,a,b=frac.split_type[u]
                for port,flow in (('left',a),('right',b)):
                    req=(u,port)
                    if flow==1:
                        tax=decoded.leaf_child[req];cid=tax;self.leaf_node[tax]=tax
                    else:
                        cid=rec_unit(decoded.child_unit[req])
                    e=len(self.edge);self.edge.append((nid,cid))
                    if req in getattr(frac,'initial_branch_length',{}):
                        init=float(frac.initial_branch_length[req])
                    else:
                        h=frac.conditional_long_intensity.get(req,0.0);ts=float(short_time(flow));init=ts*(1.0+h*(long_ratio-1.0))
                    self.edge_init.append(max(init,1e-8));chs.append((cid,e))
                self.children[nid]=chs;return nid
            self.root=rec_unit(decoded.root_unit)
        else:
            def rec(t):
                if isinstance(t,int):self.leaf_node[t]=t;return t
                nid=self._next;self._next+=1;chs=[]
                for c in t:
                    cid=rec(c);e=len(self.edge);self.edge.append((nid,cid));self.edge_init.append(0.1);chs.append((cid,e))
                self.children[nid]=chs;return nid
            self.root=rec(topology)
        self.n_edges=len(self.edge)
        self.post=[]
        def post(u):
            if u in self.children:
                for c,e in self.children[u]:post(c)
                self.post.append(u)
        post(self.root)
        self.pre=[]
        def pre(u):
            if u in self.children:
                self.pre.append(u)
                for c,e in self.children[u]:pre(c)
        pre(self.root)

    def newick(self,lengths):
        edge_len={(u,c):lengths[e] for e,(u,c) in enumerate(self.edge)}
        def rec(u):
            if u < len(self.labels):return _quote_label(self.labels[u])
            parts=[]
            for c,e in self.children[u]:parts.append(f'{rec(c)}:{edge_len[u,c]:.10g}')
            return '('+','.join(parts)+')'
        return rec(self.root)+';'


class TreeLikelihood:
    def __init__(self,topology,labels,sequences,model,decoded=None,frac=None,short_time=None,long_ratio=2.0):
        self.model=model;self.pi=np.asarray(model.pi,float);self.Q=np.asarray(model.q,float)
        self.leaf,self.weights=compress_alignment(sequences);self.Pn=len(self.weights)
        self.layout=_TreeLayout(topology,labels,decoded=decoded,frac=frac,short_time=short_time,long_ratio=long_ratio)

    def value_grad(self,log_lengths):
        length=np.exp(np.asarray(log_lengths,float));E=self.layout.n_edges
        P=[expm(self.Q*t) for t in length];dP=[self.Q@x for x in P]
        partial={i:self.leaf[i] for i in range(len(self.leaf))};msg={};dmsg={}
        for u in self.layout.post:
            prod=np.ones((self.Pn,4),float)
            for c,e in self.layout.children[u]:
                m=partial[c]@P[e].T;msg[e]=m;dmsg[e]=partial[c]@dP[e].T;prod*=m
            partial[u]=prod
        rootp=partial[self.layout.root];site=rootp@self.pi
        if np.any(site<=0) or not np.isfinite(site).all():return 1e300,np.zeros(E)
        ll=float(np.dot(self.weights,np.log(site)))
        adj={self.layout.root:self.weights[:,None]*self.pi[None,:]/site[:,None]};gt=np.zeros(E,float)
        for u in self.layout.pre:
            a=adj[u];chs=self.layout.children[u]
            for idx,(c,e) in enumerate(chs):
                sibling=np.ones((self.Pn,4),float)
                for j,(c2,e2) in enumerate(chs):
                    if j!=idx:sibling*=msg[e2]
                am=a*sibling
                gt[e]=float(np.sum(am*dmsg[e]))
                if c in self.layout.children:adj[c]=am@P[e]
        glog=gt*length
        return -ll,-glog

    def saturation_length(self,tol=1e-12):
        """Length beyond which every nonstationary CTMC eigenmode is <= tol.

        Q is normalized to one expected substitution/site.  For ordinary models this is
        O(10--30), while extreme GTR parameterizations can automatically receive a much
        larger range.  Lengths beyond this point are likelihood-indistinguishable from
        an infinite branch to approximately ``tol`` in the slowest spectral mode.
        """
        ev=np.linalg.eigvals(self.Q)
        rates=[-float(z.real) for z in ev if z.real < -1e-12]
        if not rates:return 1e6
        gap=min(rates)
        return max(10.0,-math.log(float(tol))/gap)

    def optimize(self,min_length=1e-12,max_length=None,maxiter=200,ftol=1e-10,saturation_tol=1e-12,multistart=4):
        """Optimize all branch lengths with bounded log-length L-BFGS-B.

        A small deterministic multistart set is used by default.  This is especially important
        after a fractional LP solve: an LP branch-effect variable can sit near a saturation
        endpoint even when that value is only weakly identified, and starting a continuous
        optimizer there can produce near-zero gradients.
        """
        if min_length <= 0:min_length=1e-12
        effective_max=self.saturation_length(saturation_tol) if max_length is None else float(max_length)
        if not math.isfinite(effective_max) or effective_max <= min_length:raise ValueError('Maximum branch length must exceed the minimum')
        init=np.clip(np.asarray(self.layout.edge_init,float),min_length,effective_max)
        starts=[init]
        n=max(1,int(multistart))
        if n>1:
            # Finite, nonsaturated starts spanning short through moderately long branches.
            lo=max(1e-3,min_length*10);hi=min(1.5,effective_max*0.25)
            if hi<=lo:hi=min(effective_max*0.5,max(lo*10,0.1))
            for v in np.geomspace(lo,hi,n-1):starts.append(np.full_like(init,float(v)))
        bounds=[(math.log(min_length),math.log(effective_max))]*len(init)
        best=None;tot_nit=0;tot_nfev=0;messages=[]
        for st in starts:
            x0=np.log(np.clip(st,min_length,effective_max))
            res=minimize(lambda x:self.value_grad(x),x0,jac=True,method='L-BFGS-B',bounds=bounds,
                         options={'maxiter':int(maxiter),'ftol':float(ftol),'gtol':1e-7,'maxls':40})
            tot_nit+=int(getattr(res,'nit',0));tot_nfev+=int(getattr(res,'nfev',0));messages.append(str(res.message))
            if best is None or float(res.fun)<float(best.fun):best=res
        lengths=np.exp(best.x);ll=-float(best.fun)
        msg=str(best.message) if n==1 else f'best of {len(starts)} starts: {best.message}'
        return OptimizedTree(ll,self.layout.newick(lengths),lengths,bool(best.success),tot_nit,tot_nfev,msg,float(min_length),float(effective_max))
