#!/usr/bin/env python3
"""Post-solve decomposition of a fractional universal-tree LP into integral trees.

This module is deliberately separate from the likelihood relaxation.  It does
not add integrality to the optimization model.  After a continuous LP solve it:

1. reads shared topology quantities (split-unit activity and named-taxon
   membership on each lower port),
2. generates a diverse pool of valid integral rooted binary trees by randomized
   rounding/beam-like local choices,
3. fits simplex weights theta_j >= 0, sum(theta)=1 so the weighted integral-tree
   feature vectors approximate the fractional topology, and
4. optionally accounts for the shared lower-arm branch-length fractions.

The default matching features are topology-only:
  * split-unit activity,
  * named-taxon membership at every selected split port,
  * root taxon-pair moments when the model contains them.

Branch lengths are intentionally handled *after* topology matching by default.
For a physical lower arm e with fractional long activity ell_e and reconstructed
integral-tree activity a_hat_e, the best continuous post-fit long fraction is

    h_e = clip(ell_e / a_hat_e, 0, 1).

This cannot distort topology proportions and is sufficient when branch lengths
will later be re-estimated by a fast integral-tree routine.  A binary-length
"before" mode is included only for comparison/testing.

The generated trees are embedded in the physical split-unit library: a split
unit is used at most once, its fixed left/right descendant counts are obeyed,
and the child subtrees have the requested counts.  Since the fixed-flow Beneš
router is universal for distinct count-m outputs, these selected units can be
routed integrally by the existing generator/routing machinery.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import math
import random

import numpy as np
from scipy.optimize import minimize


TaxonSet = Tuple[int, ...]
FeatureKey = Tuple


@dataclass
class IntegralNode:
    unit_id: Optional[str]
    taxa: TaxonSet
    left: Optional["IntegralNode"] = None
    right: Optional["IntegralNode"] = None
    leaf_taxon: Optional[int] = None

    @property
    def is_leaf(self) -> bool:
        return self.leaf_taxon is not None


@dataclass
class IntegralTreeCandidate:
    root_unit: str
    root: IntegralNode
    heuristic_score: float
    # (unit,port) -> tuple of taxa on that lower side
    port_taxa: Dict[Tuple[str, str], TaxonSet] = field(default_factory=dict)
    selected_units: Tuple[str, ...] = ()
    newick: str = ""
    binary_long: Dict[Tuple[str, str], float] = field(default_factory=dict)


@dataclass
class MixtureComponent:
    proportion: float
    newick: str
    root_unit: str
    selected_units: Tuple[str, ...]
    heuristic_score: float


@dataclass
class TreeMixtureResult:
    components: List[MixtureComponent]
    topology_explained_fraction: float
    topology_rmse: float
    candidate_count: int
    length_postfit_rmse: Optional[float] = None
    length_explained_fraction: Optional[float] = None
    inferred_long_fraction: Dict[Tuple[str, str], float] = field(default_factory=dict)
    inferred_length_multiplier: Dict[Tuple[str, str], float] = field(default_factory=dict)
    notes: str = ""


def _gumbel(rng: random.Random) -> float:
    u = min(max(rng.random(), 1e-12), 1 - 1e-12)
    return -math.log(-math.log(u))


def _stable_soft_choice(items, scores, rng: random.Random, temperature: float):
    if len(items) == 1 or temperature <= 1e-12:
        return items[int(np.argmax(scores))]
    s = np.asarray(scores, float)
    z = (s - float(np.max(s))) / max(float(temperature), 1e-9)
    z = np.clip(z, -60, 0)
    p = np.exp(z); p /= p.sum()
    u = rng.random(); c = 0.0
    for item, pp in zip(items, p):
        c += float(pp)
        if u <= c:
            return item
    return items[-1]


class FractionalTreeDecomposer:
    """Extract high-mass integral trees from a solved fractional tree model.

    The model is expected to expose the objects used by the current research LP:
      model.act[u], model.tmem[m,0,slot,t], model.ell[u,port], model.rpair,
      model.splits, model.roots, model.units, model.by, model._req_slot(), and
      model.lp.names.  ``result`` is the scipy/HiGHS result returned by solve().
    """

    def __init__(self, model, result):
        self.model = model
        self.result = result
        self.x = np.asarray(result.x, float)
        self.L = int(model.L)
        self.units = model.units
        self.splits = list(model.splits)
        self.roots = list(model.roots)
        self.by_count = {int(k): list(v) for k, v in model.by.items()}
        self.taxa_labels = list(model.lib.get("taxa", [f"taxon_{i+1}" for i in range(self.L)]))

        self.act = {u: float(self.x[model.act[u]]) for u in self.splits}
        self.port_mem: Dict[Tuple[str, str], np.ndarray] = {}
        self.unit_mem: Dict[str, np.ndarray] = {}
        for u in self.splits:
            tot = np.zeros(self.L)
            for port in ("left", "right"):
                m, slot, _ = model._req_slot(u, port)
                v = np.array([self.x[model.tmem[m, 0, slot, t]] for t in range(self.L)], float)
                self.port_mem[u, port] = v
                tot += v
            self.unit_mem[u] = tot
        self.ell = {}
        for key, idx in getattr(model, "ell", {}).items():
            self.ell[key] = float(self.x[idx])
        self.rpair_target = {}
        for key, idx in getattr(model, "rpair", {}).items():
            self.rpair_target[key] = float(self.x[idx])

    # ---------- candidate generation ----------
    def _unit_fit(self, u: str, taxa: Sequence[int]) -> float:
        a = max(self.act[u], 1e-10)
        p = np.clip(self.unit_mem[u] / a, 0.0, 1.0)
        overlap = float(p[list(taxa)].sum()) / max(len(taxa), 1)
        # Activity is important but not allowed to overwhelm membership agreement.
        return 2.5 * overlap + 0.35 * math.log(a + 1e-8)

    def _choose_unit(self, m: int, taxa: Sequence[int], used: set, rng: random.Random,
                     temperature: float, pool: int) -> Optional[str]:
        cand = [u for u in self.by_count.get(m, []) if u not in used]
        if not cand:
            return None
        scores = np.array([self._unit_fit(u, taxa) for u in cand], float)
        order = np.argsort(-scores)[:max(1, min(pool, len(cand)))]
        cc = [cand[i] for i in order]; ss = [scores[i] for i in order]
        return _stable_soft_choice(cc, ss, rng, temperature)

    def _partition(self, u: str, taxa: Sequence[int], rng: random.Random,
                   temperature: float) -> Tuple[TaxonSet, TaxonSet, float]:
        left_n = int(self.units[u]["left_flow"])
        act = max(self.act[u], 1e-10)
        pL = np.clip(self.port_mem[u, "left"] / act, 1e-9, 1.0)
        pR = np.clip(self.port_mem[u, "right"] / act, 1e-9, 1.0)
        taxa = list(taxa)
        # Conditional log-odds is a more stable rounding statistic than pL-pR.
        z = np.array([math.log(pL[t]) - math.log(pR[t]) for t in taxa], float)
        if temperature > 1e-12:
            z = z + np.array([temperature * _gumbel(rng) for _ in taxa])
        idx = np.argsort(-z)
        Lset = tuple(sorted(taxa[i] for i in idx[:left_n]))
        Rset = tuple(sorted(t for t in taxa if t not in set(Lset)))
        # Port agreement in [roughly 0,2]; used only to rank/generate candidates.
        agree = (float(self.port_mem[u, "left"][list(Lset)].sum()) +
                 float(self.port_mem[u, "right"][list(Rset)].sum())) / act
        agree /= max(len(taxa), 1)
        return Lset, Rset, agree

    def _build_node(self, taxa: TaxonSet, used: set, rng: random.Random,
                    temperature: float, pool: int, forced_unit: Optional[str] = None):
        m = len(taxa)
        if m == 1:
            return IntegralNode(None, taxa, leaf_taxon=taxa[0]), 0.0, {}, []
        u = forced_unit or self._choose_unit(m, taxa, used, rng, temperature, pool)
        if u is None or u in used:
            return None
        used.add(u)
        Lset, Rset, part_score = self._partition(u, taxa, rng, temperature)
        left = self._build_node(Lset, used, rng, temperature, pool)
        if left is None:
            used.remove(u); return None
        right = self._build_node(Rset, used, rng, temperature, pool)
        if right is None:
            # unwind all units consumed in left as well by using a local used copy at caller level.
            return None
        lnode, ls, lp, lu = left; rnode, rs, rp, ru = right
        ports = {(u, "left"): Lset, (u, "right"): Rset}; ports.update(lp); ports.update(rp)
        units = [u] + lu + ru
        score = self._unit_fit(u, taxa) + part_score + ls + rs
        return IntegralNode(u, taxa, lnode, rnode), score, ports, units

    def _node_newick(self, n: IntegralNode) -> str:
        if n.is_leaf:
            label = self.taxa_labels[n.leaf_taxon]
            return str(label)
        a = self._node_newick(n.left); b = self._node_newick(n.right)
        # Canonical text only; preserve model left/right in port_taxa separately.
        return f"({a},{b})"

    def generate_candidates(self, n_samples: int = 128, seed: int = 1,
                            temperature: float = 0.18, candidate_pool: int = 4,
                            include_greedy: bool = True,
                            binary_length_sampling: bool = False) -> List[IntegralTreeCandidate]:
        rng = random.Random(seed)
        roots = [r for r in self.roots if self.act[r] > 1e-10]
        if not roots:
            roots = list(self.roots)
        rscore = [math.log(self.act[r] + 1e-8) for r in roots]
        unique: Dict[str, IntegralTreeCandidate] = {}

        attempts = n_samples + (len(roots) if include_greedy else 0)
        for k in range(attempts):
            local_rng = random.Random(rng.randrange(1 << 62))
            if include_greedy and k < len(roots):
                root = roots[k]
                temp = 0.0
            else:
                root = _stable_soft_choice(roots, rscore, local_rng, max(temperature, 0.05))
                temp = temperature
            used = set()
            # isolate failed recursion so it cannot leave stale used units
            out = self._build_node(tuple(range(self.L)), used, local_rng, temp, candidate_pool, root)
            if out is None:
                continue
            node, score, ports, units = out
            nw = self._node_newick(node) + ";"
            cand = IntegralTreeCandidate(root, node, score, ports, tuple(sorted(units)), nw)
            if binary_length_sampling:
                for u in units:
                    a = max(self.act[u], 1e-12)
                    for port in ("left", "right"):
                        h = min(max(self.ell.get((u, port), 0.0) / a, 0.0), 1.0)
                        cand.binary_long[u, port] = 1.0 if local_rng.random() < h else 0.0
            old = unique.get(nw)
            if old is None or cand.heuristic_score > old.heuristic_score:
                unique[nw] = cand
        return sorted(unique.values(), key=lambda c: c.heuristic_score, reverse=True)

    # ---------- feature matching ----------
    def _candidate_features(self, c: IntegralTreeCandidate, include_pairs: bool,
                            include_binary_lengths: bool) -> Dict[FeatureKey, float]:
        f: Dict[FeatureKey, float] = {}
        for u in c.selected_units:
            f[("act", u)] = 1.0
        for (u, port), taxa in c.port_taxa.items():
            m, slot, _ = self.model._req_slot(u, port)
            for t in taxa:
                f[("tm", m, slot, t)] = 1.0
        if include_pairs and self.rpair_target:
            Ltax = set(c.port_taxa[c.root_unit, "left"])
            for (r, a, b) in self.rpair_target:
                if r == c.root_unit and a in Ltax and b in Ltax:
                    f[("pair", r, a, b)] = 1.0
        if include_binary_lengths:
            for key, v in c.binary_long.items():
                if key[0] in c.selected_units:
                    f[("ell", key[0], key[1])] = float(v)
        return f

    def _target_value(self, key: FeatureKey) -> float:
        typ = key[0]
        if typ == "act": return self.act[key[1]]
        if typ == "tm":
            _, m, slot, t = key
            return float(self.x[self.model.tmem[m, 0, slot, t]])
        if typ == "pair": return self.rpair_target[(key[1], key[2], key[3])]
        if typ == "ell": return self.ell.get((key[1], key[2]), 0.0)
        raise KeyError(key)

    def _group_counts(self, include_pairs: bool, include_binary_lengths: bool):
        out = {"act": max(len(self.act), 1), "tm": max(len(self.model.tmem), 1)}
        if include_pairs and self.rpair_target: out["pair"] = len(self.rpair_target)
        if include_binary_lengths and self.ell: out["ell"] = len(self.ell)
        return out

    def _group_weights(self, include_pairs: bool, include_binary_lengths: bool):
        # Membership carries most topology identity; activities and pair moments supplement it.
        w = {"act": 1.0, "tm": 3.0}
        if include_pairs and self.rpair_target: w["pair"] = 1.25
        if include_binary_lengths and self.ell: w["ell"] = 0.75
        return w

    def _target_norm2(self, counts, weights, include_pairs, include_binary_lengths):
        total = 0.0
        sc = {g: weights[g] / counts[g] for g in counts}
        total += sc["act"] * sum(v*v for v in self.act.values())
        total += sc["tm"] * sum(float(self.x[idx])**2 for idx in self.model.tmem.values())
        if include_pairs and self.rpair_target:
            total += sc["pair"] * sum(v*v for v in self.rpair_target.values())
        if include_binary_lengths and self.ell:
            total += sc["ell"] * sum(v*v for v in self.ell.values())
        return total

    def fit_mixture(self, candidates: Sequence[IntegralTreeCandidate],
                    include_pairs: bool = True,
                    include_binary_lengths: bool = False,
                    min_proportion: float = 1e-4) -> TreeMixtureResult:
        if not candidates:
            return TreeMixtureResult([], 0.0, float("inf"), 0, notes="No integral candidates generated.")
        feats = [self._candidate_features(c, include_pairs, include_binary_lengths) for c in candidates]
        counts = self._group_counts(include_pairs, include_binary_lengths)
        weights = self._group_weights(include_pairs, include_binary_lengths)
        scale = {g: weights[g] / counts[g] for g in counts}
        n = len(candidates)
        G = np.zeros((n, n)); b = np.zeros(n)
        # Cache scaled sparse dictionaries and target dots.
        sf = []
        for d in feats:
            dd = {k: v * math.sqrt(scale[k[0]]) for k, v in d.items() if k[0] in scale}
            sf.append(dd)
        for i, d in enumerate(sf):
            b[i] = sum(v * math.sqrt(scale[k[0]]) * self._target_value(k) for k, v in d.items())
            G[i, i] = sum(v*v for v in d.values())
            for j in range(i):
                a, c = (d, sf[j]) if len(d) < len(sf[j]) else (sf[j], d)
                z = sum(v * c.get(k, 0.0) for k, v in a.items())
                G[i, j] = G[j, i] = z
        def fun(th): return 0.5 * float(th @ G @ th) - float(b @ th)
        def jac(th): return G @ th - b
        x0 = np.full(n, 1.0/n)
        opt = minimize(fun, x0, jac=jac, method="SLSQP",
                       bounds=[(0.0, 1.0)]*n,
                       constraints=[{"type":"eq", "fun":lambda th: float(th.sum()-1.0),
                                     "jac":lambda th: np.ones_like(th)}],
                       options={"ftol":1e-11, "maxiter":1000, "disp":False})
        theta = np.maximum(opt.x if opt.success else x0, 0.0)
        theta /= max(theta.sum(), 1e-15)
        target_norm2 = self._target_norm2(counts, weights, include_pairs, include_binary_lengths)
        residual2 = max(float(theta @ G @ theta - 2*b @ theta + target_norm2), 0.0)
        explained = 1.0 - residual2/max(target_norm2, 1e-15)
        comps=[]
        for j in np.argsort(-theta):
            if theta[j] < min_proportion: continue
            c=candidates[j]
            comps.append(MixtureComponent(float(theta[j]), c.newick, c.root_unit,
                                          c.selected_units, c.heuristic_score))
        out = TreeMixtureResult(comps, float(explained), math.sqrt(residual2), len(candidates))
        # Topology-first post-fit branch lengths.
        self._postfit_lengths(out, candidates, theta)
        return out

    def _postfit_lengths(self, out: TreeMixtureResult, candidates, theta):
        if not self.ell:
            return
        recon_act={k:0.0 for k in self.ell}
        for j,c in enumerate(candidates):
            if theta[j] <= 0: continue
            sel=set(c.selected_units)
            for key in self.ell:
                if key[0] in sel: recon_act[key] += float(theta[j])
        inferred={}; err2=0.0; norm2=0.0
        for key,target in self.ell.items():
            den=recon_act[key]
            h=0.0 if den<=1e-15 else min(max(target/den,0.0),1.0)
            inferred[key]=h
            rec=den*h
            err2+=(target-rec)**2; norm2+=target*target
        out.inferred_long_fraction=inferred
        ratio=float(getattr(self.model, "long_ratio", 2.0))
        out.inferred_length_multiplier={k:1.0+h*(ratio-1.0) for k,h in inferred.items()}
        out.length_postfit_rmse=math.sqrt(err2/max(len(self.ell),1))
        out.length_explained_fraction=1.0-err2/max(norm2,1e-15)

    def decompose(self, n_samples: int = 128, seed: int = 1,
                  temperature: float = 0.18, candidate_pool: int = 4,
                  include_pairs: bool = True,
                  length_strategy: str = "after",
                  min_proportion: float = 1e-4,
                  max_fit_candidates: int = 160) -> TreeMixtureResult:
        """Generate candidates and fit mixture proportions.

        length_strategy:
          "after"  (default): topology/pair features determine theta; continuous
                    branch-length fractions are fitted afterwards.
          "ignore": topology only, with no post-fit interpretation requested.
          "binary-before": sample short/long endpoint indicators for candidates
                    and include them in the mixture fit.  This is mainly an
                    evaluation mode; it is not recommended when branch lengths
                    will later be re-estimated on recovered integral trees.
        """
        binary = length_strategy == "binary-before"
        cands = self.generate_candidates(n_samples, seed, temperature, candidate_pool,
                                         include_greedy=True, binary_length_sampling=binary)
        if max_fit_candidates and len(cands) > max_fit_candidates:
            cands = cands[:max_fit_candidates]
        out = self.fit_mixture(cands, include_pairs=include_pairs,
                               include_binary_lengths=binary,
                               min_proportion=min_proportion)
        if length_strategy == "ignore":
            out.length_postfit_rmse=None;out.length_explained_fraction=None;out.inferred_long_fraction={};out.inferred_length_multiplier={}
        out.notes = ("Topology-first matching; branch lengths fitted after proportions." if length_strategy=="after"
                     else "Binary short/long endpoint features included before topology mixture fitting." if binary
                     else "Branch lengths ignored.")
        return out


def extract_integral_tree_mixture(model, result, **kwargs) -> TreeMixtureResult:
    """Convenience function requested by the main model/API."""
    return FractionalTreeDecomposer(model, result).decompose(**kwargs)
