#!/usr/bin/env python3
"""Scalable saturated-composition root envelope for the actual 100-taxon LP.

The topology and named-taxon layers are unchanged: all physical split units,
fixed-flow routers, and boundary taxon memberships are instantiated.  Only the
composition preprocessing used to build root likelihood majorants is replaced by
a fixed-state saturated DP.  Values 0..T-1 are exact and T means ``>= T``.

The affine root planes are valid for the saturated DP bound.  Slopes are fitted on
a representative cloud, then each intercept is shifted against every feasible exact
root-side composition and every vertex of the active branch-effect polytope.
"""
from __future__ import annotations

from collections import defaultdict
import math
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PIPELINE = ROOT / "pipeline_v3"
sys.path.insert(0, str(PIPELINE))

import multisite_lowerarm_taxonflow as branch  # noqa: E402
import multisite_pruned_general as pruned  # noqa: E402
import test_rootonly_pruning as rootonly  # noqa: E402
import fractional_tree_hierarchical_units as hierarchy  # noqa: E402


def _leaf_signature(base_index, threshold):
    sig = [0] * 4
    sig[base_index] = 1 if threshold > 1 else threshold
    return tuple(sig)


def _add_status(left, right, threshold):
    if left == threshold or right == threshold:
        return threshold
    total = left + right
    return total if total < threshold else threshold


def _combine(left, right, threshold):
    return tuple(_add_status(left[i], right[i], threshold) for i in range(4))


def _status_bounds(signature, threshold, global_counts):
    lower = []
    upper = []
    for i, value in enumerate(signature):
        if value < threshold:
            lower.append(value)
            upper.append(value)
        else:
            lower.append(threshold)
            upper.append(global_counts[i])
    return tuple(lower), tuple(upper)


def _signature_possible(signature, size, threshold, global_counts):
    lower, upper = _status_bounds(signature, threshold, global_counts)
    return all(lower[i] <= global_counts[i] for i in range(4)) and sum(lower) <= size <= sum(upper)


def saturated_boxes(global_counts, threshold):
    """Return componentwise upper partials indexed by (subtree size, signature)."""
    counts = tuple(int(x) for x in global_counts)
    leaf_count = sum(counts)
    raw = {}
    by_size = defaultdict(list)
    for base_index in range(4):
        if counts[base_index] <= 0:
            continue
        signature = _leaf_signature(base_index, threshold)
        value = np.zeros(4)
        value[base_index] = 1.0
        raw[1, signature] = value
        if signature not in by_size[1]:
            by_size[1].append(signature)

    for size in range(2, leaf_count + 1):
        accumulated = {}
        for left_size in range(1, size // 2 + 1):
            right_size = size - left_size
            if not by_size[left_size] or not by_size[right_size]:
                continue
            left_vertices = branch.branch_transition_vertices(left_size)
            right_vertices = branch.branch_transition_vertices(right_size)
            left_messages = {
                signature: np.max(
                    np.stack([matrix @ raw[left_size, signature] for _, matrix in left_vertices]), axis=0
                )
                for signature in by_size[left_size]
            }
            right_messages = {
                signature: np.max(
                    np.stack([matrix @ raw[right_size, signature] for _, matrix in right_vertices]), axis=0
                )
                for signature in by_size[right_size]
            }
            for left_sig in by_size[left_size]:
                left_value = left_messages[left_sig]
                for right_sig in by_size[right_size]:
                    signature = _combine(left_sig, right_sig, threshold)
                    if not _signature_possible(signature, size, threshold, counts):
                        continue
                    value = left_value * right_messages[right_sig]
                    if signature not in accumulated:
                        accumulated[signature] = value.copy()
                    else:
                        accumulated[signature] = np.maximum(accumulated[signature], value)
        by_size[size] = sorted(accumulated)
        for signature, value in accumulated.items():
            raw[size, signature] = value
    return raw, by_size


def feasible_compositions(global_counts, left_size):
    """Enumerate exact four-base counts on one root side; <=456,976 total at L=100."""
    g = tuple(int(x) for x in global_counts)
    rows = []
    lo_a = max(0, left_size - g[1] - g[2] - g[3])
    hi_a = min(g[0], left_size)
    for a in range(lo_a, hi_a + 1):
        lo_c = max(0, left_size - a - g[2] - g[3])
        hi_c = min(g[1], left_size - a)
        for c in range(lo_c, hi_c + 1):
            lo_g = max(0, left_size - a - c - g[3])
            hi_g = min(g[2], left_size - a - c)
            for gg in range(lo_g, hi_g + 1):
                t = left_size - a - c - gg
                if 0 <= t <= g[3]:
                    rows.append((a, c, gg, t))
    return np.asarray(rows, dtype=float)


def signature_of_counts(counts, threshold):
    return tuple(min(int(value), threshold) for value in counts)


class SaturatedPrunedModel(rootonly.RootOnly):
    def __init__(
        self,
        library,
        columns,
        substitution_model=None,
        branch_relaxation="two-point",
        spectral_support_points=5,
        saturation_threshold=2,
        root_planes=4,
        log_ratio=1.1,
        log_floor=1e-12,
        long_ratio=2.0,
    ):
        self.saturation_threshold = int(saturation_threshold)
        self.root_plane_count = int(root_planes)
        self.log_floor = float(log_floor)
        if not (0 < self.log_floor < 1):
            raise ValueError("log_floor must lie strictly between zero and one")
        self._saturated_cache = {}
        self._composition_cache = {}
        self._root_state_cache = {}
        self._plane_cache = {}
        if substitution_model is not None:
            branch.configure_substitution_model(substitution_model)
        branch.configure_branch_relaxation(branch_relaxation, spectral_support_points)
        super().__init__(
            library,
            columns,
            log_ratio=log_ratio,
            long_ratio=long_ratio,
            pair_moments=False,
            max_first_facets=root_planes,
            max_quad_facets=0,
        )
        self.branch_relaxation = branch.branch_relaxation_mode()

    def _prep(self):
        cache = {}
        for pattern in self.patterns:
            counts = tuple(pattern.count(base) for base in branch.BASES)
            if counts not in cache:
                raw, by_size = saturated_boxes(counts, self.saturation_threshold)
                scale = {}
                for size in range(1, self.L + 1):
                    values = [raw[size, signature] for signature in by_size[size]]
                    scale[size] = max(float(np.max(np.stack(values))), 1e-300)
                cache[counts] = {
                    "counts": counts,
                    "sat_raw": raw,
                    "sat_by": by_size,
                    "scale": scale,
                }
            self.info.append(cache[counts])
        self.preprocess_count_classes = len(cache)

    def _compositions(self, counts, left_size):
        key = (tuple(counts), int(left_size))
        if key not in self._composition_cache:
            self._composition_cache[key] = feasible_compositions(counts, left_size)
        return self._composition_cache[key]

    def _root_states(self, info, left_size):
        counts = tuple(info["counts"])
        key = (counts, int(left_size), branch.branch_relaxation_mode(), branch.SPECTRAL_SUPPORT_POINTS)
        if key in self._root_state_cache:
            return self._root_state_cache[key]
        right_size = self.L - left_size
        compositions = self._compositions(counts, left_size)
        groups = defaultdict(list)
        for row in compositions.astype(int):
            complement = tuple(counts[i] - int(row[i]) for i in range(4))
            left_sig = signature_of_counts(row, self.saturation_threshold)
            right_sig = signature_of_counts(complement, self.saturation_threshold)
            groups[left_sig, right_sig].append(tuple(int(x) for x in row))

        left_vertices = branch.branch_transition_vertices(left_size, self.long_ratio)
        right_vertices = branch.branch_transition_vertices(right_size, self.long_ratio)
        states = []
        upper = 0.0
        norm = info["scale"][self.L]
        for (left_sig, right_sig), exact_rows in groups.items():
            left_raw = info["sat_raw"][left_size, left_sig]
            right_raw = info["sat_raw"][right_size, right_sig]
            left_z = np.stack([np.asarray(z, float) for z, _ in left_vertices])
            right_z = np.stack([np.asarray(z, float) for z, _ in right_vertices])
            left_msg = np.stack([matrix @ left_raw for _, matrix in left_vertices])
            right_msg = np.stack([matrix @ right_raw for _, matrix in right_vertices])
            score = ((left_msg * np.asarray(branch.PI)[None, :]) @ right_msg.T) / norm
            upper = max(upper, float(np.max(score)))
            states.append(
                {
                    "compositions": np.asarray(exact_rows, float),
                    "left_z": left_z,
                    "right_z": right_z,
                    "score": score,
                }
            )
        result = (states, upper)
        self._root_state_cache[key] = result
        return result

    def _build_sites(self):
        self.qr = {}
        self.lam = {}
        self.obj = {}
        self.log_error_bound = 0.0
        points = branch.vl.geometric_breaks(self.log_floor, self.log_ratio)
        for pattern_id, _ in enumerate(self.patterns):
            info = self.info[pattern_id]
            for root in self.roots:
                left_size = int(self.units[root]["left_flow"])
                _, upper = self._root_states(info, left_size)
                self.qr[pattern_id, root] = self.lp.var(f"qr:{pattern_id}:{root}", 0, upper)
            lambdas = []
            for point_id, point in enumerate(points):
                variable = self.lp.var(f"lam:{pattern_id}:{point_id}", 0, 1)
                lambdas.append(variable)
                self.obj[variable] = self.weights[pattern_id] * math.log(point)
            self.lam[pattern_id] = lambdas
            self.lp.add_eq({variable: 1 for variable in lambdas}, 1)
            row = {self.qr[pattern_id, root]: 1 for root in self.roots}
            for point_id, variable in enumerate(lambdas):
                row[variable] = row.get(variable, 0) - points[point_id]
            self.lp.add_eq(row, 0)
            self.log_error_bound += self.weights[pattern_id] * branch.vl.chord_gap(self.log_ratio)

    @staticmethod
    def _sample_cloud(states, max_compositions=160, max_branch_pairs=96):
        features = []
        values = []
        for state in states:
            comps = state["compositions"]
            if len(comps) > max_compositions:
                comp_ids = np.unique(np.linspace(0, len(comps) - 1, max_compositions, dtype=int))
                comps = comps[comp_ids]
            nl, nr = state["score"].shape
            pair_ids = np.arange(nl * nr)
            if len(pair_ids) > max_branch_pairs:
                pair_ids = np.unique(np.linspace(0, len(pair_ids) - 1, max_branch_pairs, dtype=int))
            for pair_id in pair_ids:
                li, ri = divmod(int(pair_id), nr)
                branch_features = np.concatenate([state["left_z"][li], state["right_z"][ri]])
                block = np.column_stack(
                    [comps[:, :3], np.tile(branch_features, (len(comps), 1))]
                )
                features.append(block)
                values.append(np.full(len(comps), float(state["score"][li, ri])))
        return np.vstack(features), np.concatenate(values)

    @staticmethod
    def _shift_intercept(slopes, states):
        comp_slope = np.asarray(slopes[:3], float)
        dim = (len(slopes) - 3) // 2
        left_slope = np.asarray(slopes[3 : 3 + dim], float)
        right_slope = np.asarray(slopes[3 + dim :], float)
        intercept = -math.inf
        for state in states:
            comp_residual = -state["compositions"][:, :3] @ comp_slope
            branch_residual = (
                state["score"]
                - (state["left_z"] @ left_slope)[:, None]
                - (state["right_z"] @ right_slope)[None, :]
            )
            intercept = max(intercept, float(np.max(comp_residual)) + float(np.max(branch_residual)))
        return intercept + 1e-11

    def _planes(self, pattern_id, root):
        info = self.info[pattern_id]
        left_size = int(self.units[root]["left_flow"])
        key = (
            tuple(info["counts"]),
            left_size,
            branch.branch_relaxation_mode(),
            branch.SPECTRAL_SUPPORT_POINTS,
            self.saturation_threshold,
            self.root_plane_count,
        )
        if key in self._plane_cache:
            return self._plane_cache[key]
        states, _ = self._root_states(info, left_size)
        features, values = self._sample_cloud(states)
        fitted = pruned.shifted_majorant_planes(
            features,
            values,
            K=max(1, self.root_plane_count - 1),
            fit_cap=2500,
        )
        slope_vectors = [np.zeros(features.shape[1])]
        slope_vectors.extend(np.asarray(plane[1:], float) for plane in fitted)
        planes = []
        seen = set()
        for slopes in slope_vectors:
            intercept = self._shift_intercept(slopes, states)
            plane = tuple([intercept] + slopes.tolist())
            marker = tuple(np.round(plane, 11))
            if marker not in seen:
                seen.add(marker)
                planes.append(plane)
        self._plane_cache[key] = planes
        return planes

    def _add_root_length_partition_facets(self):
        self.root_length_facets = 0
        dim = branch.branch_feature_dim()
        for pattern_id, pattern in enumerate(self.patterns):
            for root in self.roots:
                count, slot, _ = self._req_slot(root, "left")
                activity = self.act[root]
                for plane in self._planes(pattern_id, root):
                    intercept, d_a, d_c, d_g, *branch_coef = plane
                    row = {self.qr[pattern_id, root]: 1, activity: -intercept}
                    for j, variable in enumerate(self.branch_effect[root, "left"]):
                        row[variable] = row.get(variable, 0) - branch_coef[j]
                    for j, variable in enumerate(self.branch_effect[root, "right"]):
                        row[variable] = row.get(variable, 0) - branch_coef[dim + j]
                    for taxon, base_value in enumerate(pattern):
                        coefficient = {"A": d_a, "C": d_c, "G": d_g, "T": 0.0}.get(base_value, 0.0)
                        if coefficient:
                            variable = self.tmem[count, 0, slot, taxon]
                            row[variable] = row.get(variable, 0) - coefficient
                    self.lp.add_le(row, 0)
                    self.root_length_facets += 1


def select_fingerprint_patterns(columns, count=4, offset=0):
    """Greedily choose columns that refine taxon signatures.

    The criterion minimizes the number of still-colliding taxon pairs.  ``offset``
    selects a different high-entropy starting column for independent multi-view runs.
    """
    unique = list(dict.fromkeys(tuple(column) for column in columns))
    if not unique:
        raise ValueError("At least one alignment column is required")
    taxa = len(unique[0])
    signatures = [tuple() for _ in range(taxa)]
    chosen = []
    available = list(range(len(unique)))
    for step in range(min(int(count), len(unique))):
        scored = []
        for index in available:
            groups = defaultdict(int)
            column = unique[index]
            for taxon in range(taxa):
                groups[signatures[taxon] + (column[taxon],)] += 1
            collisions = sum(size * (size - 1) // 2 for size in groups.values())
            scored.append((collisions, index))
        scored.sort()
        pick_position = min(int(offset), len(scored) - 1) if step == 0 else 0
        _, selected = scored[pick_position]
        chosen.append(unique[selected])
        available.remove(selected)
        signatures = [signatures[taxon] + (unique[selected][taxon],) for taxon in range(taxa)]
    return chosen


def _category(values):
    code = 0
    for value in values:
        code = 4 * code + branch.BASE_INDEX[value]
    return code


class FingerprintSaturatedModel(SaturatedPrunedModel):
    """Full split-unit LP with joint-column taxon fingerprints at router boundaries."""

    def __init__(
        self,
        library,
        columns,
        fingerprint_columns=4,
        fingerprint_offset=0,
        anchor_weighting="uniform",
        **kwargs,
    ):
        anchors = select_fingerprint_patterns(columns, fingerprint_columns, fingerprint_offset)
        self.anchor_columns = tuple(anchors)
        self.fingerprint_blocks = tuple(
            tuple(range(start, min(start + 2, len(anchors))))
            for start in range(0, len(anchors), 2)
        )
        self.pattern_to_anchor = {tuple(pattern): index for index, pattern in enumerate(anchors)}
        taxa = len(anchors[0])
        self.block_taxon_category = {}
        self.block_categories = {}
        self.block_category_counts = {}
        for block_id, block in enumerate(self.fingerprint_blocks):
            taxon_codes = tuple(_category(tuple(anchors[pos][taxon] for pos in block)) for taxon in range(taxa))
            categories = tuple(sorted(set(taxon_codes)))
            self.block_taxon_category[block_id] = taxon_codes
            self.block_categories[block_id] = categories
            self.block_category_counts[block_id] = {
                category: taxon_codes.count(category) for category in categories
            }
        super().__init__(library, anchors, **kwargs)
        self.anchor_weighting = anchor_weighting
        if anchor_weighting not in {"uniform", "nearest"}:
            raise ValueError("anchor_weighting must be 'uniform' or 'nearest'")
        if anchor_weighting == "nearest":
            # Assign every alignment column to its closest selected anchor.  The
            # LP still has only len(anchors) site layers.  Dividing by alignment
            # size makes the complete primary objective a per-site average, so
            # coefficient magnitudes remain O(1) as the alignment grows.
            weights = np.zeros(len(anchors), dtype=float)
            anchor_array = np.asarray(anchors)
            for column in columns:
                distances = np.count_nonzero(anchor_array != np.asarray(column), axis=1)
                weights[int(np.argmin(distances))] += 1.0
            total = max(float(weights.sum()), 1.0)
            self.weights = weights / total
            points = branch.vl.geometric_breaks(self.log_floor, self.log_ratio)
            self.obj = {
                variable: float(self.weights[pattern_id]) * math.log(points[point_id])
                for pattern_id, variables in self.lam.items()
                for point_id, variable in enumerate(variables)
            }
            self.log_error_bound = branch.vl.chord_gap(self.log_ratio)

    def _build_taxon_flow(self):
        self.tmem = {}
        self.fmem = {}
        self.leaf_tmem = {}
        lp = self.lp
        # Joint-category memberships on every real request boundary.
        for block_id in range(len(self.fingerprint_blocks)):
            categories = self.block_categories[block_id]
            counts = self.block_category_counts[block_id]
            for unit in self.splits:
                activity = self.act[unit]
                for port in ("left", "right"):
                    flow, slot, _ = self._req_slot(unit, port)
                    row = {activity: -float(flow)}
                    for category in categories:
                        capacity = float(counts[category])
                        variable = lp.var(f"fm:{block_id}:{flow}:{slot}:{category}", 0, capacity)
                        self.fmem[block_id, flow, slot, category] = variable
                        lp.add_le({variable: 1, activity: -capacity}, 0)
                        row[variable] = 1
                    lp.add_eq(row, 0)

            # A subtree cannot contain more members of a category than exist globally.
            # Root units contain every category exactly at their root activity.
            for unit in self.splits:
                activity = self.act[unit]
                for category in categories:
                    row = {activity: -float(counts[category])}
                    for port in ("left", "right"):
                        flow, slot, _ = self._req_slot(unit, port)
                        variable = self.fmem[block_id, flow, slot, category]
                        row[variable] = row.get(variable, 0) + 1
                    if unit in self.roots:
                        lp.add_eq(row, 0)
                    else:
                        lp.add_le(row, 0)

            # Projected fixed-flow conservation for every joint category.
            for flow in range(1, self.L):
                for category in categories:
                    row = {}
                    for unit, port in self.reqdef[flow]:
                        _, slot, _ = self._req_slot(unit, port)
                        variable = self.fmem[block_id, flow, slot, category]
                        row[variable] = row.get(variable, 0) + 1
                    if flow == 1:
                        lp.add_eq(row, float(counts[category]))
                    else:
                        for child in self.by[flow]:
                            for port in ("left", "right"):
                                child_flow, slot, _ = self._req_slot(child, port)
                                variable = self.fmem[block_id, child_flow, slot, category]
                                row[variable] = row.get(variable, 0) - 1
                        lp.add_eq(row, 0)

        # Flow-1 ports are few (432 at L=100), so retain exact named-taxon
        # membership there.  This couples the separate fingerprint blocks and
        # makes the final leaf assignment exact without paying for named taxa on
        # every internal port.
        leaf_requests = list(self.reqdef[1])
        for unit, port in leaf_requests:
            _, slot, _ = self._req_slot(unit, port)
            activity = self.act[unit]
            row = {activity: -1.0}
            for taxon in range(self.L):
                variable = lp.var(f"leaf_tm:{slot}:{taxon}", 0, 1)
                self.leaf_tmem[slot, taxon] = variable
                lp.add_le({variable: 1, activity: -1}, 0)
                row[variable] = 1
            lp.add_eq(row, 0)
            for block_id in range(len(self.fingerprint_blocks)):
                codes = self.block_taxon_category[block_id]
                for category in self.block_categories[block_id]:
                    link = {self.fmem[block_id, 1, slot, category]: 1}
                    for taxon in range(self.L):
                        if codes[taxon] == category:
                            link[self.leaf_tmem[slot, taxon]] = -1
                    lp.add_eq(link, 0)
        for taxon in range(self.L):
            lp.add_eq(
                {
                    self.leaf_tmem[self._req_slot(unit, port)[1], taxon]: 1
                    for unit, port in leaf_requests
                },
                1,
            )

    def _anchor_location(self, pattern):
        anchor_index = self.pattern_to_anchor[tuple(pattern)]
        for block_id, block in enumerate(self.fingerprint_blocks):
            if anchor_index in block:
                return block_id, block.index(anchor_index)
        raise KeyError(pattern)

    def _add_root_length_partition_facets(self):
        self.root_length_facets = 0
        dim = branch.branch_feature_dim()
        for pattern_id, pattern in enumerate(self.patterns):
            block_id, position = self._anchor_location(pattern)
            block = self.fingerprint_blocks[block_id]
            categories = self.block_categories[block_id]
            for root in self.roots:
                flow, slot, _ = self._req_slot(root, "left")
                activity = self.act[root]
                for plane in self._planes(pattern_id, root):
                    intercept, d_a, d_c, d_g, *branch_coef = plane
                    row = {self.qr[pattern_id, root]: 1, activity: -intercept}
                    for j, variable in enumerate(self.branch_effect[root, "left"]):
                        row[variable] = row.get(variable, 0) - branch_coef[j]
                    for j, variable in enumerate(self.branch_effect[root, "right"]):
                        row[variable] = row.get(variable, 0) - branch_coef[dim + j]
                    width = len(block)
                    divisor = 4 ** (width - position - 1)
                    for category in categories:
                        base_index = (category // divisor) % 4
                        coefficient = (d_a, d_c, d_g, 0.0)[base_index]
                        if coefficient:
                            variable = self.fmem[block_id, flow, slot, category]
                            row[variable] = row.get(variable, 0) - coefficient
                    self.lp.add_le(row, 0)
                    self.root_length_facets += 1


class ProjectedFingerprintSaturatedModel(FingerprintSaturatedModel):
    """Project out physical Beneš switches but retain every split-unit activity.

    Units of a fixed descendant count are structurally interchangeable to the
    production decoder.  The projected layer keeps the exact count-wise supply and
    demand equations used by the routing network, while avoiding hundreds of
    thousands of wire/switch variables that do not identify a parent-child path.
    """

    def _build_topology(self):
        self.routermeta = {}
        for size_text, router in self.lib["routers"].items():
            size = int(size_text)
            width = int(router["padded_width"])
            depth = len(router["dimensions"])
            self.routermeta[size] = (router, width, depth)
            row = {}
            for request in router["input_requests"]:
                unit = request["parent_unit"]
                row[self.act[unit]] = row.get(self.act[unit], 0) + 1
            if size == 1:
                self.lp.add_eq(row, float(self.L))
            else:
                for target in router["output_targets"]:
                    unit = target["unit_id"]
                    row[self.act[unit]] = row.get(self.act[unit], 0) - 1
                self.lp.add_eq(row, 0)


def build_fingerprint_fractional_tree(model, result, numerical_zero=1e-12):
    """Map joint-category marginals to taxon plausibilities for beam decoding."""
    x = np.asarray(result.x, dtype=float).copy()
    if numerical_zero is not None and numerical_zero > 0:
        x[np.abs(x) < float(numerical_zero)] = 0.0
    labels = tuple(model.lib.get("taxa", [f"leaf_{i+1}" for i in range(model.L)]))
    activity = {unit: max(float(x[model.act[unit]]), 0.0) for unit in model.splits}
    roots = tuple(sorted(model.roots))
    split_type = {}
    child_requests = {}
    by_flow = defaultdict(list)
    leaf_mass = {}
    port_distribution = {}
    unit_distribution = {}
    conditional = {}
    initial = {}

    for unit in model.splits:
        record = model.units[unit]
        total = int(record["subtree_size"])
        left = int(record["left_flow"])
        right = int(record["right_flow"])
        split_type[unit] = (total, left, right)
        by_flow[total].append(unit)
        child_requests[unit] = ((unit, "left"), (unit, "right"))
        whole = np.zeros(model.L)
        for port, flow in (("left", left), ("right", right)):
            request = (unit, port)
            _, slot, _ = model._req_slot(unit, port)
            per_view = []
            for block_id in range(len(model.fingerprint_blocks)):
                codes = model.block_taxon_category[block_id]
                values = np.asarray(
                    [
                        max(float(x[model.fmem[block_id, flow, slot, codes[taxon]]]), 0.0)
                        for taxon in range(model.L)
                    ],
                    float,
                )
                # Correct for category multiplicity so common signatures do not dominate.
                counts = model.block_category_counts[block_id]
                values = np.asarray([values[t] / counts[codes[t]] for t in range(model.L)])
                per_view.append(values)
            plausibility = np.min(np.stack(per_view), axis=0)
            if plausibility.sum() > 0:
                distribution = plausibility / plausibility.sum()
            else:
                distribution = plausibility
            port_distribution[request] = distribution
            whole += flow * activity[unit] * distribution
            if flow == 1:
                for taxon in range(model.L):
                    leaf_mass[request, taxon] = max(float(x[model.leaf_tmem[slot, taxon]]), 0.0)
                exact = np.asarray([leaf_mass[request, taxon] for taxon in range(model.L)])
                if exact.sum() > 0:
                    port_distribution[request] = exact / exact.sum()
            if (unit, port) in getattr(model, "ell", {}):
                raw = max(float(x[model.ell[unit, port]]), 0.0)
                conditional[request] = 0.0 if activity[unit] <= 1e-14 else min(raw / activity[unit], 1.0)
            else:
                conditional[request] = 0.0
        unit_distribution[unit] = whole / whole.sum() if whole.sum() > 0 else whole

    return hierarchy.HierarchicalFractionalTree(
        labels,
        activity,
        roots,
        split_type,
        child_requests,
        {flow: tuple(sorted(units)) for flow, units in by_flow.items()},
        leaf_mass,
        port_distribution,
        unit_distribution,
        conditional,
        initial,
    )


def rank_fingerprint_candidates(model, result, topk=20, **kwargs):
    fractional = build_fingerprint_fractional_tree(
        model, result, numerical_zero=kwargs.get("numerical_zero", 1e-12)
    )
    residual = hierarchy.Residual(fractional)
    all_candidates = {}
    root_beam = kwargs.get("root_beam", 128)
    leaf_beam = kwargs.get("leaf_beam", 32)
    child_branching = kwargs.get("child_branching", 12)
    leaf_branching = kwargs.get("leaf_branching", 8)
    merge_branching = kwargs.get("merge_branching", 96)
    max_expansions = kwargs.get("max_expansions", 500000)
    for candidate in hierarchy.root_down_topk(
        fractional, residual, topk, root_beam, child_branching, leaf_branching, max_expansions
    ):
        all_candidates[repr(candidate.topology)] = candidate
    for candidate in hierarchy.leaves_up_topk(
        fractional, residual, topk, leaf_beam, merge_branching, max_expansions
    ):
        key = repr(candidate.topology)
        old = all_candidates.get(key)
        if old is None or (candidate.proportion, candidate.secondary_score) > (old.proportion, old.secondary_score):
            all_candidates[key] = candidate
    ranked = sorted(
        all_candidates.values(), key=lambda item: (-item.proportion, -item.secondary_score, item.newick)
    )[:topk]
    return ranked, fractional
