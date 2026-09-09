#!/usr/bin/env python3
"""Search the full rooted-binary tree space to decompose a fractional tree LP.

This is post-solve only.  It does NOT add integer variables to the likelihood LP.

Core distinction
----------------
For a split unit u and lower arm e=(u,port):

    frequency(e) = activity(u)
    conditional_long_intensity(e) = ell[e] / activity(u)   (when activity>0)

The raw LP long activity ell[e] = frequency * conditional length intensity, so it
must not be used as a topology-frequency feature.  Branch length remains essential
inside the LP likelihood calculation, but topology recovery uses only branch
frequency/taxon-state structure.  Length intensity can be attached after a tree is
identified or re-estimated by a fast integral-tree routine.

Tree-space search
-----------------
The previous prototype generated a finite candidate pool and then fit simplex
weights.  This module instead uses a Frank-Wolfe / matching-pursuit decomposition
whose linear oracle searches the rooted-binary tree space itself.

For L <= exact_taxa the oracle enumerates every rooted binary labeled topology.
For larger L it uses multi-start rooted-NNI search plus randomized NNI walks.  The
rooted-NNI graph is connected, so the search domain is the complete tree space;
there is no pre-generated candidate set.  The large-L oracle is heuristic rather
than a proof of the global maximizer, which is unavoidable for practical 100-taxon
postprocessing unless one uses a much more expensive exact combinatorial solver.

For a fixed topology, assignment of its internal nodes to the physical split-unit
copies is optimized exactly by a Hungarian assignment within each fixed split type.
Equal-size child orientation is also chosen to maximize the current oracle score.

The outer Frank-Wolfe updates are analytic (no LP/QP mixture fit): each newly found
integral tree is added with an exact line-search step.  The resulting component
weights are the reported approximate proportions.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union
import math
import random

import numpy as np
from scipy.optimize import linear_sum_assignment

Tree = Union[int, Tuple["Tree", "Tree"]]
FeatureKey = Tuple
TaxonSet = Tuple[int, ...]


@dataclass
class SearchedIntegralTree:
    topology: Tree
    newick: str
    root_unit: str
    selected_units: Tuple[str, ...]
    port_taxa: Dict[Tuple[str, str], TaxonSet]
    feature_dict: Dict[FeatureKey, float]
    oracle_score: float


@dataclass
class SearchComponent:
    proportion: float
    newick: str
    root_unit: str
    selected_units: Tuple[str, ...]


@dataclass
class TreeSpaceMixtureResult:
    components: List[SearchComponent]
    topology_explained_fraction: float
    topology_rmse: float
    frank_wolfe_iterations: int
    oracle_calls: int
    unique_trees_seen: int
    exact_oracle: bool
    arm_frequency: Dict[Tuple[str, str], float] = field(default_factory=dict)
    conditional_long_intensity: Dict[Tuple[str, str], float] = field(default_factory=dict)
    reconstructed_arm_frequency: Dict[Tuple[str, str], float] = field(default_factory=dict)
    notes: str = ""


def _taxa(t: Tree) -> Tuple[int, ...]:
    if isinstance(t, int):
        return (t,)
    return tuple(sorted(_taxa(t[0]) + _taxa(t[1])))


def _canon(t: Tree) -> Tree:
    if isinstance(t, int):
        return int(t)
    a, b = _canon(t[0]), _canon(t[1])
    ka, kb = _taxa(a), _taxa(b)
    return (a, b) if ka <= kb else (b, a)


def _key(t: Tree):
    if isinstance(t, int):
        return ("L", t)
    return ("N", _key(t[0]), _key(t[1]))


def _random_tree(leaves: Sequence[int], rng: random.Random) -> Tree:
    q: List[Tree] = [int(x) for x in leaves]
    while len(q) > 1:
        i = rng.randrange(len(q)); a = q.pop(i)
        j = rng.randrange(len(q)); b = q.pop(j)
        q.append(_canon((a, b)))
    return q[0]


def _nni_neighbors(t: Tree) -> List[Tree]:
    """All rooted-NNI neighbors, deduplicated canonically."""
    if isinstance(t, int):
        return []
    A, B = t
    out: Dict[Tuple, Tree] = {}
    # NNI moves on the edge root--A/root--B.
    if not isinstance(A, int):
        a1, a2 = A
        for z in (_canon((_canon((a1, B)), a2)), _canon((_canon((a2, B)), a1))):
            out[_key(z)] = z
    if not isinstance(B, int):
        b1, b2 = B
        for z in (_canon((_canon((b1, A)), b2)), _canon((_canon((b2, A)), b1))):
            out[_key(z)] = z
    # NNI moves internal to each child subtree.
    for na in _nni_neighbors(A):
        z = _canon((na, B)); out[_key(z)] = z
    for nb in _nni_neighbors(B):
        z = _canon((A, nb)); out[_key(z)] = z
    out.pop(_key(_canon(t)), None)
    return list(out.values())


@lru_cache(maxsize=None)
def _enumerate_subset(mask: int, L: int) -> Tuple[Tree, ...]:
    bits = [i for i in range(L) if mask >> i & 1]
    if len(bits) == 1:
        return (bits[0],)
    first = bits[0]
    restmask = mask & ~(1 << first)
    ans: Dict[Tuple, Tree] = {}
    sub = restmask
    while True:
        leftmask = sub | (1 << first)
        rightmask = mask ^ leftmask
        if rightmask:
            for a in _enumerate_subset(leftmask, L):
                for b in _enumerate_subset(rightmask, L):
                    z = _canon((a, b)); ans[_key(z)] = z
        if sub == 0:
            break
        sub = (sub - 1) & restmask
    return tuple(ans.values())


def enumerate_rooted_binary_trees(L: int) -> Iterable[Tree]:
    return _enumerate_subset((1 << L) - 1, L)


class TreeSpaceDecomposer:
    def __init__(self, model, result, include_pairs: bool = True):
        self.model = model
        self.result = result
        self.x = np.asarray(result.x, float)
        self.L = int(model.L)
        self.units = model.units
        self.splits = list(model.splits)
        self.roots = set(model.roots)
        self.taxa_labels = list(model.lib.get("taxa", [f"taxon_{i+1}" for i in range(self.L)]))
        self.include_pairs = bool(include_pairs and getattr(model, "rpair", None))

        self.act = {u: float(self.x[model.act[u]]) for u in self.splits}
        self.ell = {k: float(self.x[idx]) for k, idx in getattr(model, "ell", {}).items()}
        self.rpair = {k: float(self.x[idx]) for k, idx in getattr(model, "rpair", {}).items()}
        self.tmem_target = {k: float(self.x[idx]) for k, idx in model.tmem.items()}

        self.by_type: Dict[Tuple[int, int], List[str]] = {}
        for u in self.splits:
            m = int(self.units[u]["subtree_size"])
            a = int(self.units[u]["split_a"])
            self.by_type.setdefault((m, a), []).append(u)
        for v in self.by_type.values():
            v.sort()

        self.counts = {
            "act": max(len(self.act), 1),
            "tm": max(len(self.tmem_target), 1),
        }
        if self.include_pairs:
            self.counts["pair"] = max(len(self.rpair), 1)
        gw = {"act": 1.0, "tm": 3.0, "pair": 1.25}
        self.scale = {g: gw[g] / self.counts[g] for g in self.counts}
        self.target_norm2 = (
            self.scale["act"] * sum(v*v for v in self.act.values())
            + self.scale["tm"] * sum(v*v for v in self.tmem_target.values())
            + (self.scale.get("pair", 0.0) * sum(v*v for v in self.rpair.values()) if self.include_pairs else 0.0)
        )
        self.target_cache: Dict[FeatureKey, float] = {}
        self.eval_cache: Dict[Tuple, SearchedIntegralTree] = {}
        self.unique_seen: set = set()

    # ---- frequency / length normalization ----
    def arm_frequency_and_length_intensity(self):
        freq: Dict[Tuple[str, str], float] = {}
        h: Dict[Tuple[str, str], float] = {}
        for u in self.splits:
            a = max(self.act[u], 0.0)
            for p in ("left", "right"):
                k = (u, p)
                freq[k] = a
                e = self.ell.get(k, 0.0)
                h[k] = 0.0 if a <= 1e-14 else min(max(e/a, 0.0), 1.0)
        return freq, h

    def _target(self, k: FeatureKey) -> float:
        if k in self.target_cache:
            return self.target_cache[k]
        typ = k[0]
        if typ == "act": v = self.act.get(k[1], 0.0)
        elif typ == "tm": v = self.tmem_target.get((k[1], k[2], k[3], k[4]), 0.0)
        elif typ == "pair": v = self.rpair.get((k[1], k[2], k[3]), 0.0)
        else: raise KeyError(k)
        self.target_cache[k] = float(v)
        return float(v)

    def _w(self, k: FeatureKey) -> float:
        return self.scale[k[0]]

    # ---- map an abstract topology to physical split-unit copies ----
    def _collect_nodes(self, t: Tree, is_root=True):
        out = []
        def rec(z: Tree, rootflag: bool):
            if isinstance(z, int): return
            A, B = z; ta, tb = _taxa(A), _taxa(B)
            if len(ta) <= len(tb): small, large = A, B
            else: small, large = B, A; ta, tb = tb, ta
            out.append({"tree": z, "m": len(ta)+len(tb), "a": len(ta), "small": small, "large": large,
                        "equal": len(ta)==len(tb), "root": rootflag})
            rec(A, False); rec(B, False)
        rec(t, is_root)
        return out

    def _orientation_options(self, node):
        s, l = node["small"], node["large"]
        if node["equal"]:
            return [(_taxa(s), _taxa(l)), (_taxa(l), _taxa(s))]
        return [(_taxa(s), _taxa(l))]

    def _local_features(self, u: str, Ltax: TaxonSet, Rtax: TaxonSet, is_root: bool) -> Dict[FeatureKey, float]:
        f: Dict[FeatureKey, float] = {("act", u): 1.0}
        for port, taxa in (("left", Ltax), ("right", Rtax)):
            m, slot, _ = self.model._req_slot(u, port)
            for t in taxa:
                f[("tm", m, 0, slot, t)] = 1.0
        if is_root and self.include_pairs:
            ss = set(Ltax)
            for key in self.rpair:
                r, a, b = key
                if r == u and a in ss and b in ss:
                    f[("pair", r, a, b)] = 1.0
        return f

    @staticmethod
    def _dot_coeff(f: Dict[FeatureKey, float], coeff: Callable[[FeatureKey], float]) -> float:
        return sum(v * coeff(k) for k, v in f.items())

    def map_and_score(self, topology: Tree, coeff: Callable[[FeatureKey], float], use_cache=False) -> Optional[SearchedIntegralTree]:
        topology = _canon(topology); kk = _key(topology); self.unique_seen.add(kk)
        nodes = self._collect_nodes(topology)
        groups: Dict[Tuple[int, int], List[dict]] = {}
        for n in nodes: groups.setdefault((n["m"], n["a"]), []).append(n)
        selected_units=[]; port_taxa={}; allfeat={}; total=0.0; root_unit=None
        for typ, ns in groups.items():
            us = self.by_type.get(typ, [])
            if len(us) < len(ns): return None
            M=np.full((len(ns),len(us)),-1e100); choices={}
            for i,n in enumerate(ns):
                for j,u in enumerate(us):
                    best=None
                    for lt,rt in self._orientation_options(n):
                        f=self._local_features(u,lt,rt,n["root"])
                        sc=self._dot_coeff(f,coeff)
                        if best is None or sc>best[0]: best=(sc,lt,rt,f)
                    M[i,j]=best[0];choices[i,j]=best
            rows,cols=linear_sum_assignment(-M)
            for i,j in zip(rows,cols):
                n=ns[i];u=us[j];sc,lt,rt,f=choices[i,j]
                total+=float(sc);selected_units.append(u);port_taxa[u,"left"]=lt;port_taxa[u,"right"]=rt
                for k,v in f.items(): allfeat[k]=allfeat.get(k,0.0)+v
                if n["root"]: root_unit=u
        if root_unit is None: return None
        def nw(z):
            if isinstance(z,int): return str(self.taxa_labels[z])
            return f"({nw(z[0])},{nw(z[1])})"
        return SearchedIntegralTree(topology,nw(topology)+";",root_unit,tuple(sorted(selected_units)),port_taxa,allfeat,total)

    # ---- full-tree-space oracle ----
    def _hill_climb(self, start: Tree, coeff, max_steps: int, rng: random.Random,
                    random_walk_steps: int = 0, walk_temperature: float = 0.05):
        cur=_canon(start); best=self.map_and_score(cur,coeff); cache={_key(cur):best}
        if best is None: return None
        for _ in range(max_steps):
            nb=_nni_neighbors(cur); bestn=best; bestt=cur
            for z in nb:
                kz=_key(z); ev=cache.get(kz)
                if ev is None: ev=self.map_and_score(z,coeff);cache[kz]=ev
                if ev is not None and ev.oracle_score>bestn.oracle_score+1e-12:
                    bestn=ev;bestt=z
            if bestt is cur or _key(bestt)==_key(cur): break
            cur,best=bestt,bestn
        # Random NNI walk followed by another climb gives an inexpensive basin escape.
        if random_walk_steps>0:
            z=cur;ze=best
            for _ in range(random_walk_steps):
                nb=_nni_neighbors(z)
                if not nb: break
                cand=rng.choice(nb);ev=self.map_and_score(cand,coeff)
                if ev is None: continue
                delta=ev.oracle_score-ze.oracle_score
                if delta>=0 or rng.random()<math.exp(max(delta,-50.0)/max(walk_temperature,1e-9)):
                    z,ze=cand,ev
            fin=self._hill_climb(z,coeff,max_steps,rng,0,walk_temperature)
            if fin is not None and fin.oracle_score>best.oracle_score: best=fin
        return best

    def oracle(self, coeff: Callable[[FeatureKey], float], *, exact_taxa: int = 6,
               restarts: int = 10, hill_steps: int = 60, random_walk_steps: int = 18,
               seed: int = 1, warm_start: Optional[Tree] = None):
        if self.L <= exact_taxa:
            best=None
            for t in enumerate_rooted_binary_trees(self.L):
                ev=self.map_and_score(t,coeff)
                if ev is not None and (best is None or ev.oracle_score>best.oracle_score): best=ev
            return best, True
        rng=random.Random(seed);starts=[]
        if warm_start is not None: starts.append(warm_start)
        # A deterministic caterpillar start plus random full-space starts.
        cat=0
        for t in range(1,self.L): cat=_canon((cat,t))
        starts.append(cat)
        while len(starts)<max(restarts,2): starts.append(_random_tree(range(self.L),rng))
        best=None
        for st in starts:
            ev=self._hill_climb(st,coeff,hill_steps,rng,random_walk_steps,0.05)
            if ev is not None and (best is None or ev.oracle_score>best.oracle_score): best=ev
        return best, False

    # ---- sparse weighted algebra ----
    def _dot_target_feature(self, f): return sum(self._w(k)*self._target(k)*v for k,v in f.items())
    def _norm_feature(self, f): return sum(self._w(k)*v*v for k,v in f.items())
    def _dot_sparse(self,a,b):
        if len(a)>len(b): a,b=b,a
        return sum(self._w(k)*v*b.get(k,0.0) for k,v in a.items())
    def _dot_target_sparse(self,a): return sum(self._w(k)*self._target(k)*v for k,v in a.items())
    def _norm_sparse(self,a): return sum(self._w(k)*v*v for k,v in a.items())

    def decompose(self, max_components: int = 20, *, exact_taxa: int = 6,
                  restarts: int = 10, hill_steps: int = 60, random_walk_steps: int = 18,
                  seed: int = 1, min_step: float = 1e-4, min_proportion: float = 1e-3) -> TreeSpaceMixtureResult:
        # Best single tree for weighted squared distance: maximize sum w(2x-1)f for binary features.
        coeff0=lambda k: self._w(k)*(2.0*self._target(k)-1.0)
        first, exact=self.oracle(coeff0,exact_taxa=exact_taxa,restarts=restarts,hill_steps=hill_steps,
                                random_walk_steps=random_walk_steps,seed=seed)
        if first is None:
            return TreeSpaceMixtureResult([],0.0,float("inf"),0,1,0,exact,notes="No valid integral tree found.")
        trees={_key(first.topology):first}; weights={_key(first.topology):1.0}; mix=dict(first.feature_dict)
        oracle_calls=1; warm=first.topology

        for it in range(1,max_components):
            coeff=lambda k,mix=mix: self._w(k)*(self._target(k)-mix.get(k,0.0))
            tr, ex2=self.oracle(coeff,exact_taxa=exact_taxa,restarts=restarts,hill_steps=hill_steps,
                                random_walk_steps=random_walk_steps,seed=seed+104729*it,warm_start=warm)
            oracle_calls+=1; exact=exact and ex2
            if tr is None: break
            fk=tr.feature_dict
            mixnorm=self._norm_sparse(mix); fnorm=self._norm_feature(fk); mf=self._dot_sparse(mix,fk)
            txf=self._dot_target_feature(fk); txm=self._dot_target_sparse(mix)
            numer=txf-txm-mf+mixnorm
            denom=fnorm+mixnorm-2.0*mf
            gamma=0.0 if denom<=1e-15 else min(max(numer/denom,0.0),1.0)
            if gamma<min_step: break
            for k in list(weights): weights[k]*=(1.0-gamma)
            kt=_key(tr.topology); weights[kt]=weights.get(kt,0.0)+gamma;trees[kt]=tr
            keys=set(mix)|set(fk)
            mix={k:(1.0-gamma)*mix.get(k,0.0)+gamma*fk.get(k,0.0) for k in keys}
            mix={k:v for k,v in mix.items() if abs(v)>1e-14};warm=tr.topology

        # prune tiny weights for display only; reported reconstruction uses the full mixture.
        mixnorm=self._norm_sparse(mix); txm=self._dot_target_sparse(mix)
        residual2=max(self.target_norm2+mixnorm-2.0*txm,0.0)
        explained=1.0-residual2/max(self.target_norm2,1e-15)
        comps=[]
        for k,w in sorted(weights.items(),key=lambda kv:-kv[1]):
            if w>=min_proportion:
                tr=trees[k];comps.append(SearchComponent(float(w),tr.newick,tr.root_unit,tr.selected_units))

        freq,h=self.arm_frequency_and_length_intensity();recon={k:0.0 for k in freq}
        for kk,w in weights.items():
            sel=set(trees[kk].selected_units)
            for e in recon:
                if e[0] in sel: recon[e]+=w
        return TreeSpaceMixtureResult(
            comps,float(explained),math.sqrt(residual2),len(weights),oracle_calls,len(self.unique_seen),exact,
            arm_frequency=freq,conditional_long_intensity=h,reconstructed_arm_frequency=recon,
            notes=("Integral components obtained by Frank-Wolfe over a full-tree-space rooted-NNI oracle. "
                   "Branch lengths remain in the LP but are excluded from topology frequency matching: "
                   "frequency=a_u and conditional length intensity=ell/(a_u).")
        )


def search_integral_tree_mixture(model,result,**kwargs):
    return TreeSpaceDecomposer(model,result).decompose(**kwargs)
