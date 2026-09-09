#!/usr/bin/env python3
"""Shared DNA substitution-model configuration for the LP and final tree scorer.

The active model is process-global on purpose: research LP modules are imported dynamically
and historically used module-level PI/Q constants.  Keeping one shared source prevents the
LP preprocessing and post-solve likelihood scorer from silently using different models.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

BASES = "ACGT"
BASE_INDEX = {b:i for i,b in enumerate(BASES)}

@dataclass(frozen=True)
class DNAmodel:
    name: str
    pi: np.ndarray
    q: np.ndarray
    rates: tuple[float, float, float, float, float, float] | None = None


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
    """GTR+F with rates ordered AC, AG, AT, CG, CT, GT.

    Frequencies are supplied by the selected alignment columns (+F).  Tiny values are
    allowed but the caller should normally regularize an exactly absent base first.
    """
    pi=np.asarray(empirical_frequencies,float)
    if pi.shape!=(4,) or np.any(pi<0) or not np.isfinite(pi).all() or pi.sum()<=0:
        raise ValueError("Need four finite nonnegative empirical frequencies")
    pi=pi/pi.sum()
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

def set_active_model(model: DNAmodel) -> DNAmodel:
    global _ACTIVE
    _ACTIVE=model
    return model

def get_active_model() -> DNAmodel:
    return _ACTIVE
