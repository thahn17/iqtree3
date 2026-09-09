#!/usr/bin/env python3
"""Recommended taxon/likelihood matching lift for the lower-arm multisite LP.

Adds two structural/tightening layers to multisite_lowerarm_taxonflow.py:

1. Explicit fixed-flow taxon packing certificate
      sum_{u: descendant_count(u)=k} x_{u,t} <= 1
   where x_{u,t} is total taxon-t membership in split unit u.

2. Shared root taxon-pair co-membership moments y_{r,t,v}.  For a root
   candidate r whose left side has fixed size a,
      sum_{v != t} y_{r,t,v} = (a-1) x_{r,t}.
   These moments make every site's first and quadratic base-composition
   moments linear images of the SAME fractional root partition.

Site/root likelihood contributions are then constrained by affine majorants
in first + second composition moments and the two lower-arm length selectors.
The majorants are precomputed from count/composition DP values only; no subtree
or taxon-subset enumeration is required for these likelihood facets.

Pure continuous LP.
"""
from __future__ import annotations
import importlib.util, sys
from pathlib import Path
import numpy as np

BASE_PATH = '/mnt/data/multisite_lowerarm_taxonflow.py'
spec = importlib.util.spec_from_file_location('base_matching', BASE_PATH)
base = importlib.util.module_from_spec(spec)
sys.modules['base_matching'] = base
spec.loader.exec_module(base)

PI = base.PI
BASE_INDEX = base.BASE_INDEX


def affine_majorants_features(X, y, max_facets=3):
    """Small set of valid affine majorants over finitely many feature points."""
    from scipy.optimize import linprog
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    A = np.column_stack([np.ones(len(X)), X])
    if not len(A):
        return []
    if len(A) <= max_facets:
        idx = list(range(len(A)))
    else:
        order = np.argsort(y)
        idx = set([int(order[0]), int(order[-1])])
        for z in np.linspace(0, len(order)-1, max_facets, dtype=int):
            idx.add(int(order[z]))
        idx = sorted(idx)
    out = []
    for ii in idx:
        res = linprog(A[ii], A_ub=-A, b_ub=-y,
                      bounds=[(None, None)] * A.shape[1], method='highs')
        if res.success and np.all(A @ res.x + 1e-8 >= y):
            key = tuple(np.round(res.x, 11))
            if key not in out:
                out.append(key)
    return out


class TaxonStateMatchingModel(base.BoundaryTaxonModel):
    def __init__(self, libpath, columns, log_ratio=1.1, long_ratio=2.0,
                 max_root_moment_facets=3, pair_bounds='full'):
        self.max_root_moment_facets = max_root_moment_facets
        self.pair_bounds = pair_bounds
        super().__init__(libpath, columns, log_ratio, long_ratio)
        self.flow_packing_cuts = 0
        self._add_flow_taxon_packing()
        self.root_pair = {}
        self._add_root_pair_moments()
        self.root_moment_facets = 0
        self._add_root_moment_likelihood_facets()

    def _unit_taxon_expr(self, u, t):
        row = {}
        for port in ('left', 'right'):
            m, slot, _ = self._req_slot(u, port)
            v = self.tmem[m, 0, slot, t]
            row[v] = row.get(v, 0.0) + 1.0
        return row

    def _add_flow_taxon_packing(self):
        # Valid antichain inequality: a taxon can be in at most one count-k
        # subtree in an integral rooted tree.
        for k in range(2, self.L):
            for t in range(self.L):
                row = {}
                for u in self.by[k]:
                    for v, c in self._unit_taxon_expr(u, t).items():
                        row[v] = row.get(v, 0.0) + c
                if row:
                    self.lp.add_le(row, 1.0)
                    self.flow_packing_cuts += 1

    def _add_root_pair_moments(self):
        # Only the canonical left side is needed; right composition is the
        # global composition minus left composition.
        for r in self.roots:
            a = int(self.units[r]['left_flow'])
            if a <= 1:
                continue
            rho = self.act[r]
            m, slot, _ = self._req_slot(r, 'left')
            for t in range(self.L):
                for v in range(t + 1, self.L):
                    y = self.lp.var(f'root_pair:{r}:{t}:{v}', 0, 1)
                    self.root_pair[r, t, v] = y
                    if self.pair_bounds in ('upper', 'full'):
                        xt = self.tmem[m, 0, slot, t]
                        xv = self.tmem[m, 0, slot, v]
                        self.lp.add_le({y: 1, xt: -1}, 0)
                        self.lp.add_le({y: 1, xv: -1}, 0)
                    if self.pair_bounds == 'full':
                        xt = self.tmem[m, 0, slot, t]
                        xv = self.tmem[m, 0, slot, v]
                        self.lp.add_ge({y: 1, xt: -1, xv: -1, rho: 1}, 0)
            # Exact degree identity for an a-subset, preserved under convex mixtures.
            for t in range(self.L):
                row = {}
                for v in range(self.L):
                    if v == t:
                        continue
                    key = (r, min(t, v), max(t, v))
                    row[self.root_pair[key]] = row.get(self.root_pair[key], 0) + 1
                xt = self.tmem[m, 0, slot, t]
                row[xt] = row.get(xt, 0) - (a - 1)
                self.lp.add_eq(row, 0)

    def _site_moment_exprs(self, r, pattern):
        """3 first moments + 10 raw quadratic composition moments."""
        m, slot, _ = self._req_slot(r, 'left')
        expr = []
        for b in range(3):  # A,C,G; T determined by side size
            d = {}
            for t, ch in enumerate(pattern):
                if BASE_INDEX[ch] == b:
                    d[self.tmem[m, 0, slot, t]] = 1
            expr.append(d)
        for b in range(4):
            for c in range(b, 4):
                d = {}
                if b == c:
                    # n_b^2 = n_b + 2 * number of same-base taxon pairs.
                    for t, ch in enumerate(pattern):
                        if BASE_INDEX[ch] == b:
                            x = self.tmem[m, 0, slot, t]
                            d[x] = d.get(x, 0) + 1
                    for t in range(self.L):
                        if BASE_INDEX[pattern[t]] != b:
                            continue
                        for v in range(t + 1, self.L):
                            if BASE_INDEX[pattern[v]] == b and (r, t, v) in self.root_pair:
                                y = self.root_pair[r, t, v]
                                d[y] = d.get(y, 0) + 2
                else:
                    for t in range(self.L):
                        bt = BASE_INDEX[pattern[t]]
                        if bt not in (b, c):
                            continue
                        for v in range(t + 1, self.L):
                            bv = BASE_INDEX[pattern[v]]
                            if {bt, bv} == {b, c} and (r, t, v) in self.root_pair:
                                y = self.root_pair[r, t, v]
                                d[y] = d.get(y, 0) + 1
                expr.append(d)
        return expr

    def _add_root_moment_likelihood_facets(self):
        # Coefficients depend only on global site base counts, root split size,
        # branch endpoint matrices, and composition; not on taxon identities.
        for p, pattern in enumerate(self.patterns):
            I = self.info[p]
            raw = I['raw']
            counts = I['counts']
            for r in self.roots:
                a = int(self.units[r]['left_flow'])
                b = self.L - a
                rho = self.act[r]
                eL = self.ell[r, 'left']
                eR = self.ell[r, 'right']
                PsA, PlA = base.trans_pair(a, self.long_ratio)
                PsB, PlB = base.trans_pair(b, self.long_ratio)
                X, vals = [], []
                for (mm, c1), BA in raw.items():
                    if mm != a:
                        continue
                    c2 = tuple(counts[k] - c1[k] for k in range(4))
                    if min(c2) < 0 or (b, c2) not in raw:
                        continue
                    feat = [c1[0], c1[1], c1[2]]
                    for bb in range(4):
                        for cc in range(bb, 4):
                            feat.append(c1[bb] * c1[cc])
                    for la, PA in enumerate((PsA, PlA)):
                        for lb, PB in enumerate((PsB, PlB)):
                            val = float(PI @ ((PA @ BA[:, 1]) *
                                              (PB @ raw[(b, c2)][:, 1]))) / I['scale'][self.L]
                            X.append(feat + [la, lb])
                            vals.append(val)
                moment_expr = self._site_moment_exprs(r, pattern)
                qrow = {self.z[p, r][i]: float(PI[i]) for i in range(4)}
                for f in affine_majorants_features(X, vals, self.max_root_moment_facets):
                    a0, *coef = f
                    row = dict(qrow)
                    row[rho] = row.get(rho, 0) - a0
                    row[eL] = row.get(eL, 0) - coef[len(moment_expr)]
                    row[eR] = row.get(eR, 0) - coef[len(moment_expr) + 1]
                    for h, ex in enumerate(moment_expr):
                        c = coef[h]
                        if c:
                            for v, cc in ex.items():
                                row[v] = row.get(v, 0) - c * cc
                    self.lp.add_le(row, 0)
                    self.root_moment_facets += 1


if __name__ == '__main__':
    cols = [tuple('AAAACCCC'), tuple('ACGTACGT'), tuple('AACCGGTT'),
            tuple('AGCTAGCT'), tuple('AAAACCGT'), tuple('AACCGGTT')]
    M = TaxonStateMatchingModel('/mnt/data/test_match_units8/library.json', cols)
    v, res, sec = M.solve()
    print('objective', v)
    print('vars', len(M.lp.names), 'eq', len(M.lp.eq), 'le', len(M.lp.le))
    print('flow packing cuts', M.flow_packing_cuts,
          'root moment facets', M.root_moment_facets, 'solve', sec)
