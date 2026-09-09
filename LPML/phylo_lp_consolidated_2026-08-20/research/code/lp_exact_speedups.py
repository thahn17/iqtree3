#!/usr/bin/env python3
"""Exact-equivalence solve utilities for the sparse phylogeny LPs.

The transformations here do not change the feasible polytope or optimizer:

* compile the Python row dictionaries once and reuse the CSR matrices;
* scale each constraint by a positive row constant;
* scale the complete objective by a positive constant (normally site count);
* optionally discard byte-identical duplicate rows.

The returned objective is always evaluated with the original, unscaled
coefficients.  These helpers deliberately do not project the routing network or
drop fingerprint/site layers; those are approximations and must be benchmarked
separately.
"""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_array, csr_matrix, vstack


def _deduplicate(rows, rhs):
    """Remove only exactly identical sparse rows with exactly identical RHS."""
    kept_rows = []
    kept_rhs = []
    seen = set()
    for row, value in zip(rows, rhs):
        key = (tuple(sorted((int(j), float(v)) for j, v in row.items() if v != 0.0)), float(value))
        if key in seen:
            continue
        seen.add(key)
        kept_rows.append(row)
        kept_rhs.append(value)
    return kept_rows, kept_rhs


def _matrix(rows, nvars):
    if not rows:
        return None
    data = []
    row_index = []
    column_index = []
    for row_id, row in enumerate(rows):
        for variable, coefficient in row.items():
            if coefficient != 0.0:
                row_index.append(row_id)
                column_index.append(variable)
                data.append(coefficient)
    return csr_matrix(
        (np.asarray(data, dtype=np.float64), (row_index, column_index)),
        shape=(len(rows), nvars),
    )


def _scale_rows(matrix, rhs):
    """Max-norm row equilibration using positive constants only."""
    if matrix is None:
        return None, None
    matrix = matrix.copy().tocsr()
    rhs = np.asarray(rhs, dtype=np.float64).copy()
    maxima = np.zeros(matrix.shape[0], dtype=np.float64)
    for row in range(matrix.shape[0]):
        start, stop = matrix.indptr[row : row + 2]
        if stop > start:
            maxima[row] = np.max(np.abs(matrix.data[start:stop]))
    scale = np.ones_like(maxima)
    nonzero = maxima > 0.0
    scale[nonzero] = 1.0 / maxima[nonzero]
    matrix.data *= np.repeat(scale, np.diff(matrix.indptr))
    rhs *= scale
    return matrix, rhs


@dataclass
class SolveTiming:
    compile_s: float
    solve_s: float


class CompiledSparseLP:
    """A reusable SciPy/HiGHS view of the project's dictionary-based LP."""

    def __init__(self, sparse_lp, *, row_scale=False, deduplicate=False):
        started = time.perf_counter()
        le_rows, le_rhs = sparse_lp.le, sparse_lp.lerhs
        eq_rows, eq_rhs = sparse_lp.eq, sparse_lp.eqrhs
        if deduplicate:
            le_rows, le_rhs = _deduplicate(le_rows, le_rhs)
            eq_rows, eq_rhs = _deduplicate(eq_rows, eq_rhs)
        self.a_ub = _matrix(le_rows, len(sparse_lp.names))
        self.b_ub = np.asarray(le_rhs, dtype=np.float64) if le_rows else None
        self.a_eq = _matrix(eq_rows, len(sparse_lp.names))
        self.b_eq = np.asarray(eq_rhs, dtype=np.float64) if eq_rows else None
        if row_scale:
            self.a_ub, self.b_ub = _scale_rows(self.a_ub, self.b_ub)
            self.a_eq, self.b_eq = _scale_rows(self.a_eq, self.b_eq)
        self.bounds = tuple(sparse_lp.bounds)
        self.nvars = len(sparse_lp.names)
        self.compile_s = time.perf_counter() - started
        self.original_inequalities = len(sparse_lp.le)
        self.original_equalities = len(sparse_lp.eq)
        self.inequalities = 0 if self.a_ub is None else self.a_ub.shape[0]
        self.equalities = 0 if self.a_eq is None else self.a_eq.shape[0]

    def solve(
        self,
        objective,
        *,
        maximize=True,
        objective_scale=1.0,
        method="highs",
        feasibility_tolerance=1e-8,
        simplex_dual_edge_weight_strategy=None,
        extra_options=None,
    ):
        if not np.isfinite(objective_scale) or objective_scale <= 0:
            raise ValueError("objective_scale must be finite and positive")
        vector = np.zeros(self.nvars, dtype=np.float64)
        direction = -1.0 if maximize else 1.0
        for variable, coefficient in objective.items():
            vector[variable] = direction * float(coefficient) / objective_scale
        options = {
            "dual_feasibility_tolerance": feasibility_tolerance,
            "primal_feasibility_tolerance": feasibility_tolerance,
            "presolve": True,
        }
        if simplex_dual_edge_weight_strategy is not None:
            options["simplex_dual_edge_weight_strategy"] = simplex_dual_edge_weight_strategy
        if extra_options:
            options.update(extra_options)
        started = time.perf_counter()
        result = linprog(
            vector,
            A_ub=self.a_ub,
            b_ub=self.b_ub,
            A_eq=self.a_eq,
            b_eq=self.b_eq,
            bounds=self.bounds,
            method=method,
            options=options,
        )
        solve_s = time.perf_counter() - started
        if not result.success:
            raise RuntimeError(result.message)
        value = sum(float(coefficient) * result.x[variable] for variable, coefficient in objective.items())
        return value, result, SolveTiming(self.compile_s, solve_s)


class DirectHighsLP:
    """Cache SciPy's final combined CSC matrix before calling HiGHS.

    ``linprog`` rebuilds and converts ``[A_ub; A_eq]`` on every call.  This
    adapter follows SciPy's HiGHS translation exactly but caches that mechanical
    work.  It uses SciPy's private wrapper and is therefore an optional pinned-
    SciPy optimization, not a portable public API.
    """

    def __init__(self, sparse_lp):
        from scipy.optimize import _linprog_highs as scipy_highs

        started = time.perf_counter()
        compiled = CompiledSparseLP(sparse_lp)
        empty_ub = csr_matrix((0, compiled.nvars)) if compiled.a_ub is None else compiled.a_ub
        empty_eq = csr_matrix((0, compiled.nvars)) if compiled.a_eq is None else compiled.a_eq
        self.matrix = csc_array(vstack((empty_ub, empty_eq)))
        n_ub = empty_ub.shape[0]
        ub_rhs = np.empty(0) if compiled.b_ub is None else compiled.b_ub
        eq_rhs = np.empty(0) if compiled.b_eq is None else compiled.b_eq
        self.lhs = np.concatenate((np.full(n_ub, -np.inf), eq_rhs))
        self.rhs = np.concatenate((ub_rhs, eq_rhs))
        bounds = np.asarray(
            [
                (-np.inf if lower is None else lower, np.inf if upper is None else upper)
                for lower, upper in compiled.bounds
            ],
            dtype=float,
        )
        self.lower = bounds[:, 0].copy()
        self.upper = bounds[:, 1].copy()
        self.lhs = scipy_highs._replace_inf(self.lhs)
        self.rhs = scipy_highs._replace_inf(self.rhs)
        self.lower = scipy_highs._replace_inf(self.lower)
        self.upper = scipy_highs._replace_inf(self.upper)
        self.nvars = compiled.nvars
        self.compile_s = time.perf_counter() - started
        self._highs = scipy_highs

    def solve(self, objective, *, maximize=True, feasibility_tolerance=1e-8):
        h = self._highs
        vector = np.zeros(self.nvars, dtype=float)
        direction = -1.0 if maximize else 1.0
        for variable, coefficient in objective.items():
            vector[variable] = direction * float(coefficient)
        options = {
            "presolve": True,
            "sense": h.ObjSense.kMinimize,
            "solver": None,
            "time_limit": None,
            "highs_debug_level": h.HighsDebugLevel.kHighsDebugLevelNone,
            "dual_feasibility_tolerance": feasibility_tolerance,
            "ipm_optimality_tolerance": None,
            "log_to_console": False,
            "mip_max_nodes": None,
            "output_flag": False,
            "primal_feasibility_tolerance": feasibility_tolerance,
            "simplex_dual_edge_weight_strategy": None,
            "simplex_strategy": h.s_c.SimplexStrategy.kSimplexStrategyDual,
            "ipm_iteration_limit": None,
            "simplex_iteration_limit": None,
            "mip_rel_gap": None,
        }
        started = time.perf_counter()
        raw = h._highs_wrapper(
            vector,
            self.matrix.indptr,
            self.matrix.indices,
            self.matrix.data,
            self.lhs,
            self.rhs,
            self.lower,
            self.upper,
            np.empty(0, dtype=np.uint8),
            options,
        )
        solve_s = time.perf_counter() - started
        if raw.get("status") != h.HighsModelStatus.kOptimal:
            raise RuntimeError(raw.get("message", str(raw.get("status"))))
        result = type("DirectHighsResult", (), {})()
        result.x = np.asarray(raw["x"], dtype=float)
        result.success = True
        result.message = raw.get("message", "optimal")
        result.nit = raw.get("simplex_nit", 0)
        value = sum(float(coefficient) * result.x[variable] for variable, coefficient in objective.items())
        return value, result, SolveTiming(self.compile_s, solve_s)
