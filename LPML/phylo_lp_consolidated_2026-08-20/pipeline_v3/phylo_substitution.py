#!/usr/bin/env python3
"""Shared DNA substitution-model configuration for the LP and final tree scorer.

For non-JC reversible models the LP can use a 3-mode spectral branch relaxation.  A real
branch has z_r=exp(-mu_r t), r=1..3.  The relaxation outer-approximates the convex hull of
that one-dimensional curve with linear inequalities and keeps transition matrices
nonnegative.  This covers every t in [0,infinity) without requiring a literal LP length.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
import numpy as np
from scipy.spatial import ConvexHull,HalfspaceIntersection
from scipy.optimize import brentq,linprog

BASES = "ACGT"
BASE_INDEX = {b:i for i,b in enumerate(BASES)}

@dataclass(frozen=True)
class DNAmodel:
    name: str
    pi: np.ndarray
    q: np.ndarray
    rates: tuple[float, float, float, float, float, float] | None = None

@dataclass(frozen=True)
class SpectralRelaxation:
    pi: np.ndarray
    mu: np.ndarray                 # three positive decay rates, slow -> fast
    projectors: tuple[np.ndarray, ...]  # P(z)=1*pi^T + sum z_r B_r
    A: np.ndarray                  # A z <= b outer polytope
    b: np.ndarray
    vertices: np.ndarray           # vertices of A z <= b
    support_points: int

    def transition(self, z) -> np.ndarray:
        z=np.asarray(z,float)
        P=np.tile(self.pi,(4,1)).astype(float)
        for r,B in enumerate(self.projectors):P += float(z[r])*B
        return P


def _normalize_q(q: np.ndarray, pi: np.ndarray) -> np.ndarray:
    q=np.asarray(q,float).copy(); pi=np.asarray(pi,float)
    rate=-float(np.dot(pi,np.diag(q)))
    if not np.isfinite(rate) or rate<=0:
        raise ValueError("Substitution generator has non-positive expected rate")
    return q/rate


def jc69() -> DNAmodel:
    pi=np.full(4,0.25,float)
    q=np.full((4,4),1.0,float); np.fill_diagonal(q,0.0)
    for i in range(4): q[i,i]=-q[i].sum()
    q=_normalize_q(q,pi)
    return DNAmodel("JC69",pi,q,(1.,1.,1.,1.,1.,1.))


def gtr_f(empirical_frequencies, rates) -> DNAmodel:
    """GTR+F with rates ordered AC, AG, AT, CG, CT, GT."""
    pi=np.asarray(empirical_frequencies,float)
    if pi.shape!=(4,) or np.any(pi<0) or not np.isfinite(pi).all() or pi.sum()<=0:
        raise ValueError("Need four finite nonnegative empirical frequencies")
    pi=pi/pi.sum()
    # Reversible spectral transform requires positive stationary masses.  This is also what
    # production +F callers already enforce with a tiny pseudocount for absent bases.
    if np.any(pi<=0):
        raise ValueError("Spectral GTR relaxation requires positive +F frequencies; use a small frequency pseudocount")
    rr=tuple(float(x) for x in rates)
    if len(rr)!=6 or any((not np.isfinite(x) or x<=0) for x in rr):
        raise ValueError("GTR exchangeabilities AC,AG,AT,CG,CT,GT must be six positive numbers")
    R=np.zeros((4,4),float)
    pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for (i,j),r in zip(pairs,rr):R[i,j]=R[j,i]=r
    q=np.zeros((4,4),float)
    for i in range(4):
        for j in range(4):
            if i!=j:q[i,j]=R[i,j]*pi[j]
        q[i,i]=-q[i].sum()
    q=_normalize_q(q,pi)
    return DNAmodel("GTR+F",pi,q,rr)

_ACTIVE=jc69()
_SPECTRAL_CACHE={}

def set_active_model(model: DNAmodel) -> DNAmodel:
    global _ACTIVE
    _ACTIVE=model
    return model

def get_active_model() -> DNAmodel:
    return _ACTIVE


def is_jc_like(model:DNAmodel|None=None, tol=1e-9) -> bool:
    model=_ACTIVE if model is None else model
    pi=np.asarray(model.pi,float);q=np.asarray(model.q,float)
    if np.max(np.abs(pi-0.25))>tol:return False
    off=q.copy();np.fill_diagonal(off,np.nan)
    vals=off[np.isfinite(off)]
    return float(np.max(vals)-np.min(vals))<=tol*max(1.0,float(np.max(np.abs(vals))))


def spectral_decomposition(model:DNAmodel|None=None):
    """Return (mu, B) with P(t)=1*pi^T+sum_r exp(-mu_r t) B_r for reversible Q."""
    model=_ACTIVE if model is None else model
    pi=np.asarray(model.pi,float);Q=np.asarray(model.q,float)
    if np.any(pi<=0):raise ValueError('Positive stationary frequencies required')
    sq=np.sqrt(pi); inv=1.0/sq
    S=(sq[:,None]*Q)*inv[None,:]
    S=(S+S.T)/2.0
    lam,U=np.linalg.eigh(S)
    # stationary eigenvalue is the one closest to zero; remaining sorted slow -> fast
    iz=int(np.argmin(np.abs(lam)))
    ids=[i for i in range(4) if i!=iz]
    ids=sorted(ids,key=lambda i:-lam[i])
    mu=np.array([-lam[i] for i in ids],float)
    Bs=[]
    for i in ids:
        u=U[:,i]
        B=(inv[:,None]*np.outer(u,u))*sq[None,:]
        Bs.append(B)
    # numerical reconstruction checks at 0 and infinity
    P0=np.tile(pi,(4,1)).astype(float)
    for B in Bs:P0+=B
    if np.max(np.abs(P0-np.eye(4)))>5e-8:
        raise ValueError('Failed reversible spectral reconstruction')
    return mu,tuple(Bs)


def _polytope_vertices(A,b):
    """Vertices of A z<=b via a Chebyshev-center halfspace intersection."""
    A=np.asarray(A,float);b=np.asarray(b,float);d=A.shape[1];norm=np.linalg.norm(A,axis=1)
    res=linprog(np.r_[np.zeros(d),-1.0],A_ub=np.column_stack([A,norm]),b_ub=b,
                bounds=[(None,None)]*(d+1),method='highs')
    if not res.success or res.x[d]<=1e-12:raise ValueError('Could not find an interior point for spectral branch polytope')
    if d==1:
        lo=max((-b[i]/A[i,0] for i in range(len(b)) if A[i,0]<0),default=-np.inf)
        hi=min((b[i]/A[i,0] for i in range(len(b)) if A[i,0]>0),default=np.inf)
        return np.asarray([[lo],[hi]],float)
    H=HalfspaceIntersection(np.column_stack([A,-b]),res.x[:d]);V=np.asarray(H.intersections,float)
    return np.unique(np.round(V,12),axis=0)


def _merge_equal_modes(mu,Bs,tol=1e-9):
    out_mu=[];out_B=[]
    for m,B in zip(mu,Bs):
        if out_mu and abs(float(m)-out_mu[-1])<=tol*max(1.0,abs(float(m)),abs(out_mu[-1])):
            out_B[-1]=out_B[-1]+B
        else:out_mu.append(float(m));out_B.append(np.asarray(B,float).copy())
    return np.asarray(out_mu,float),tuple(out_B)


def _curve_support_halfspaces(mu, support_points):
    """Valid outer support halfspaces for z_r=x^(mu_r/mu_1), x in [0,1]."""
    d=len(mu);n=max(5,int(support_points))
    if d==1:return np.asarray([[1.0],[-1.0]]),np.asarray([1.0,0.0])
    ratios=np.asarray(mu,float)/float(mu[0])
    if d==2:
        r=float(ratios[1]);A=[[1,0],[-1,0],[-1,1]];b=[1,0,0] # 0<=x<=1, y<=x
        k=np.arange(n);xs=0.5*(1-np.cos(np.pi*k/(n-1)))
        for x0 in xs[1:]:
            slope=r*x0**(r-1);inter=(1-r)*x0**r
            A.append([slope,-1]);b.append(-inter+2e-12)
        return np.asarray(A,float),np.asarray(b,float)
    if d!=3:raise ValueError('DNA spectral relaxation supports at most three unique nonstationary modes')
    r2=float(ratios[1]);r3=float(ratios[2])
    k=np.arange(n);x=0.5*(1-np.cos(np.pi*k/(n-1)));X=np.column_stack([x,x**r2,x**r3])
    H=ConvexHull(X,qhull_options='QJ');A=[];bb=[]
    grid=np.r_[np.geomspace(1e-14,1e-3,80),np.linspace(1e-3,1,320)]
    for eq in H.equations:
        a=np.asarray(eq[:3],float);norm=float(np.linalg.norm(a))
        if norm<1e-14:continue
        a=a/norm
        def f(z):return float(a[0]*z+a[1]*z**r2+a[2]*z**r3)
        def df(z):return float(a[0]+a[1]*r2*z**(r2-1)+a[2]*r3*z**(r3-1))
        cand=[0.0,1.0];dv=np.array([df(z) for z in grid])
        for i in range(len(grid)-1):
            if dv[i]==0:cand.append(float(grid[i]))
            elif dv[i]*dv[i+1]<0:
                try:cand.append(float(brentq(df,float(grid[i]),float(grid[i+1]),xtol=2e-14)))
                except ValueError:pass
        rhs=max(f(z) for z in cand)+2e-12
        if not any(np.max(np.abs(a-aa))<1e-8 and abs(rhs-b0)<1e-8 for aa,b0 in zip(A,bb)):
            A.append(a);bb.append(rhs)
    return np.asarray(A,float),np.asarray(bb,float)


def build_spectral_relaxation(model:DNAmodel|None=None, support_points=7) -> SpectralRelaxation:
    """Tight outer polyhedral relaxation of all possible non-JC reversible branch effects.

    A real branch has z_r=exp(-mu_r t), t>=0.  We build supporting planes around that
    complete one-dimensional curve and add transition-matrix nonnegativity.  Every physical
    length from zero through infinity is included.  More support_points tighten the outer
    convex relaxation but increase root-envelope preprocessing.
    """
    model=_ACTIVE if model is None else model
    n=max(5,int(support_points));key=(tuple(np.round(model.pi,14)),tuple(np.round(model.q.ravel(),14)),n)
    if key in _SPECTRAL_CACHE:return _SPECTRAL_CACHE[key]
    mu0,Bs0=spectral_decomposition(model);mu,Bs=_merge_equal_modes(mu0,Bs0)
    A,b=_curve_support_halfspaces(mu,n)
    # Stochastic row sums are automatic in the spectral basis. Enforce P_ij(z)>=0 so every
    # relaxed branch effect remains a valid stochastic transition matrix.
    extraA=[];extrab=[]
    for i in range(4):
        for j in range(4):
            extraA.append([-Bs[r][i,j] for r in range(len(Bs))]);extrab.append(float(model.pi[j]))
    A=np.vstack([A,np.asarray(extraA,float)]);b=np.r_[b,np.asarray(extrab,float)]
    V=_polytope_vertices(A,b)
    V[np.abs(V)<2e-10]=0;V[np.abs(V-1)<2e-10]=1
    # Continuous-curve inclusion check across an extreme range of t.
    for t in np.r_[0.0,np.logspace(-8,6,200)]:
        z=np.exp(-mu*t)
        if np.max(A@z-b)>5e-8:raise ValueError('Physical branch effect excluded by spectral relaxation')
    stationary=np.tile(model.pi,(4,1))
    for z in V:
        P=stationary.copy()
        for r,B in enumerate(Bs):P+=z[r]*B
        if P.min()<-2e-7 or np.max(np.abs(P.sum(1)-1))>2e-7:
            raise ValueError('Spectral relaxation vertex produced invalid transition matrix')
    out=SpectralRelaxation(np.asarray(model.pi,float).copy(),mu,Bs,A,b,V,n)
    _SPECTRAL_CACHE[key]=out
    return out

