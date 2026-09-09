#!/usr/bin/env python3
"""
fixed_flow_unit_generator.py

Generate a compact count-stratified universal library of connected subtree units
for all rooted binary trees on L distinct taxa, with these properties:

  * Exactly one physical leaf/taxon unit per supplied taxon label.
  * Every physical unit has a known, fixed upper flow (= descendant leaf count).
  * Every logical edge inside every unit has a known, fixed leaf flow.
  * Every selectable inter-unit lower connection has at most two specific choices
    (k <= 2).
  * Every rooted binary labeled tree on the supplied taxa can be routed through
    the library.
  * The construction is O(L^2 log L) in the number of physical units.

The implementation minimizes units within this count-stratified/Beneš family by:

  1. Using the provably necessary multiplicity floor(L/m) for each fixed-flow
     split type m = a + b.
  2. Omitting a universal ROOT_ENTRY and the count-L router: a realization
     directly selects the appropriate size-L split unit as its root.
  3. Compacting each padded Beneš router by eliminating both boundary wire
     layers as standalone units. The first switch is represented directly by
     the requesting lower port's two connections; the final switch connects
     the last internal router unit directly to destination units.
  4. Omitting routers of width 1 entirely.

This is not claimed to be a globally minimum universal library over every
possible construction, but it removes the obvious redundant units while
retaining a simple constructive proof and polynomial-time routing.

Optional variable-length internal paths
---------------------------------------
With --variable-length-paths, every logical edge A--B inside every unit is
represented by two alternative paths with the SAME fixed flow:

    short: A -- B                         (1 edge)
    long : A -- h1 -- ... -- h_(r-1) -- B  (r edges, r >= 2)

They share the same endpoints and hybrid_group. A displayed tree selects one
of the two alternatives, varying the topological A-B path length without
changing connectivity endpoints or leaf flow. These helper vertices/edges are
inside a unit and do NOT count as additional subtree units.

Inputs
------
Library only:
    python fixed_flow_unit_generator.py A B C D --out-dir units4

Library + route a target rooted binary Newick tree:
    python fixed_flow_unit_generator.py \
        --tree '((A,B),(C,(D,E)));' --out-dir routed5

Enable alternate longer paths inside units:
    python fixed_flow_unit_generator.py \
        --tree '((A,B),(C,(D,E)));' \
        --variable-length-paths --long-path-edges 3 \
        --active-length-choice long \
        --out-dir routed5_long

Outputs
-------
Always:
    units.csv
    connections.csv
    unit_logical_edges.csv
    unit_internal_paths.csv
    unit_internal_connections.csv
    library.json
    summary.txt

With --tree / --tree-file:
    target_nodes.csv
    routes.csv
    active_units.csv
    active_connections.csv
    router_switches.csv
    active_internal_paths.csv
    active_internal_connections.csv
    routing.json
    routing_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------


@dataclass
class Unit:
    unit_id: str
    kind: str
    descendant_count: int
    upper_flow: int

    label: str = ""
    subtree_size: Optional[int] = None
    split_a: Optional[int] = None
    split_b: Optional[int] = None
    copy_index: Optional[int] = None

    left_flow: Optional[int] = None
    right_flow: Optional[int] = None
    down_flow: Optional[int] = None

    router_count: Optional[int] = None
    router_layer: Optional[int] = None
    router_wire: Optional[int] = None
    router_width: Optional[int] = None

    notes: str = ""


@dataclass
class Connection:
    parent_unit: str
    parent_port: str
    child_unit: str
    option_index: int
    edge_flow: int
    parent_upper_flow: int
    child_upper_flow: int
    connection_kind: str
    router_count: Optional[int] = None
    router_stage: Optional[int] = None
    switch_id: str = ""
    switch_state: str = ""


@dataclass
class LogicalEdge:
    unit_id: str
    logical_edge_id: str
    endpoint_a: str
    endpoint_b: str
    edge_flow: int


@dataclass
class InternalPath:
    unit_id: str
    logical_edge_id: str
    hybrid_group: str
    path_option: str
    path_length_edges: int
    edge_flow: int
    endpoint_a: str
    endpoint_b: str
    path_nodes: str


@dataclass
class InternalConnection:
    unit_id: str
    logical_edge_id: str
    hybrid_group: str
    path_option: str
    segment_index: int
    segment_from: str
    segment_to: str
    edge_flow: int


# ---------------------------------------------------------------------------
# Target-tree parser
# ---------------------------------------------------------------------------


@dataclass
class TargetNode:
    label: Optional[str] = None
    left: Optional["TargetNode"] = None
    right: Optional["TargetNode"] = None
    node_id: str = ""
    descendant_count: int = 0
    descendant_labels: Tuple[str, ...] = ()
    assigned_unit: str = ""

    @property
    def is_leaf(self) -> bool:
        return self.left is None and self.right is None


class NewickParser:
    """Small rooted-binary Newick parser; branch lengths/internal labels ignored."""

    def __init__(self, text: str):
        self.s = text.strip()
        self.i = 0

    def _skip_ws(self):
        while self.i < len(self.s) and self.s[self.i].isspace():
            self.i += 1

    def _peek(self) -> str:
        self._skip_ws()
        return self.s[self.i] if self.i < len(self.s) else ""

    def _consume(self, ch: str):
        self._skip_ws()
        if self.i >= len(self.s) or self.s[self.i] != ch:
            got = self.s[self.i:self.i + 24]
            raise ValueError(f"Expected {ch!r} near {got!r}")
        self.i += 1

    def _parse_label(self) -> str:
        self._skip_ws()
        if self.i >= len(self.s):
            return ""

        if self.s[self.i] == "'":
            self.i += 1
            out = []
            while self.i < len(self.s):
                if self.s[self.i] == "'":
                    if self.i + 1 < len(self.s) and self.s[self.i + 1] == "'":
                        out.append("'")
                        self.i += 2
                        continue
                    self.i += 1
                    return "".join(out)
                out.append(self.s[self.i])
                self.i += 1
            raise ValueError("Unterminated quoted Newick label")

        start = self.i
        while self.i < len(self.s) and self.s[self.i] not in "(),:;":
            self.i += 1
        return self.s[start:self.i].strip()

    def _skip_branch_length(self):
        self._skip_ws()
        if self.i < len(self.s) and self.s[self.i] == ":":
            self.i += 1
            while self.i < len(self.s) and self.s[self.i] not in ",);":
                self.i += 1

    def _parse_subtree(self) -> TargetNode:
        self._skip_ws()
        if self._peek() == "(":
            self._consume("(")
            left = self._parse_subtree()
            self._consume(",")
            right = self._parse_subtree()
            if self._peek() == ",":
                raise ValueError("Target Newick is not rooted-binary")
            self._consume(")")
            if self._peek() not in ("", ";", ",", ")", ":"):
                self._parse_label()  # internal label ignored
            self._skip_branch_length()
            return TargetNode(left=left, right=right)

        label = self._parse_label()
        if not label:
            raise ValueError("Leaf without a label in Newick")
        self._skip_branch_length()
        return TargetNode(label=label)

    def parse(self) -> TargetNode:
        root = self._parse_subtree()
        self._skip_ws()
        if self._peek() == ";":
            self.i += 1
        self._skip_ws()
        if self.i != len(self.s):
            raise ValueError(f"Unexpected trailing Newick text: {self.s[self.i:]!r}")
        return root


def annotate_target_tree(root: TargetNode) -> List[TargetNode]:
    """Canonicalize child order and assign counts, leaf sets, and IDs."""
    post: List[TargetNode] = []

    def visit(node: TargetNode):
        if node.is_leaf:
            assert node.label is not None
            node.descendant_count = 1
            node.descendant_labels = (node.label,)
        else:
            assert node.left is not None and node.right is not None
            visit(node.left)
            visit(node.right)
            lk = (node.left.descendant_count, node.left.descendant_labels)
            rk = (node.right.descendant_count, node.right.descendant_labels)
            if lk > rk:
                node.left, node.right = node.right, node.left
            node.descendant_count = node.left.descendant_count + node.right.descendant_count
            node.descendant_labels = tuple(
                sorted(node.left.descendant_labels + node.right.descendant_labels)
            )
        post.append(node)

    visit(root)
    for i, node in enumerate(post):
        node.node_id = f"T{i:05d}"
    return post


def leaf_labels_in_order(root: TargetNode) -> List[str]:
    labels: List[str] = []

    def walk(n: TargetNode):
        if n.is_leaf:
            assert n.label is not None
            labels.append(n.label)
        else:
            assert n.left is not None and n.right is not None
            walk(n.left)
            walk(n.right)

    walk(root)
    return labels


# ---------------------------------------------------------------------------
# Beneš helpers
# ---------------------------------------------------------------------------


def next_power_of_two(x: int) -> int:
    if x <= 1:
        return 1
    return 1 << (x - 1).bit_length()


def benes_dimensions(width: int) -> List[int]:
    if width == 1:
        return []
    if width & (width - 1):
        raise ValueError("Beneš width must be a power of two")
    d = int(math.log2(width))
    return list(range(d)) + list(range(d - 2, -1, -1))


def benes_route_switches(permutation: List[int]) -> List[Dict[int, int]]:
    """
    Route a permutation through the dimension-sequence power-of-two Beneš
    network. Returns one dict per switch stage: pair_base -> 0/1
    (straight/cross).
    """
    n = len(permutation)
    if n == 0 or (n & (n - 1)):
        raise ValueError("Permutation length must be a positive power of two")
    if sorted(permutation) != list(range(n)):
        raise ValueError("Input is not a permutation")
    if n == 1:
        return []

    d_total = int(math.log2(n))
    stages: List[Dict[int, int]] = [dict() for _ in range(2 * d_total - 1)]

    def recurse(p: List[int], depth: int, prefix: int, stage0: int):
        nloc = len(p)
        if nloc == 1:
            return
        dloc = int(math.log2(nloc))

        input_pair_edges: List[List[int]] = [[] for _ in range(nloc // 2)]
        output_pair_edges: List[List[int]] = [[] for _ in range(nloc // 2)]
        for inp, out in enumerate(p):
            input_pair_edges[inp // 2].append(inp)
            output_pair_edges[out // 2].append(inp)

        color: List[Optional[int]] = [None] * nloc
        for start in range(nloc):
            if color[start] is not None:
                continue
            color[start] = 0
            stack = [start]
            while stack:
                edge = stack.pop()
                assert color[edge] is not None
                in_pair = input_pair_edges[edge // 2]
                out_pair = output_pair_edges[p[edge] // 2]
                for pair in (in_pair, out_pair):
                    if len(pair) != 2:
                        raise RuntimeError("Malformed Beneš pair graph")
                    other = pair[0] if pair[1] == edge else pair[1]
                    want = 1 - int(color[edge])
                    if color[other] is None:
                        color[other] = want
                        stack.append(other)
                    elif color[other] != want:
                        raise RuntimeError("Beneš cycle coloring failed")

        last_stage = stage0 + 2 * dloc - 2

        for pair_index in range(nloc // 2):
            even_edge = 2 * pair_index
            actual_even_wire = (2 * pair_index << depth) | prefix
            pair_base = min(actual_even_wire, actual_even_wire ^ (1 << depth))
            stages[stage0][pair_base] = 0 if color[even_edge] == 0 else 1

        inverse = [0] * nloc
        for inp, out in enumerate(p):
            inverse[out] = inp
        for pair_index in range(nloc // 2):
            edge_to_even_output = inverse[2 * pair_index]
            actual_even_wire = (2 * pair_index << depth) | prefix
            pair_base = min(actual_even_wire, actual_even_wire ^ (1 << depth))
            stages[last_stage][pair_base] = 0 if color[edge_to_even_output] == 0 else 1

        for c in (0, 1):
            subperm: List[Optional[int]] = [None] * (nloc // 2)
            for inp, out in enumerate(p):
                if color[inp] == c:
                    subperm[inp // 2] = out // 2
            if any(x is None for x in subperm):
                raise RuntimeError("Incomplete Beneš sub-permutation")
            recurse(
                [int(x) for x in subperm],
                depth + 1,
                prefix | (c << depth),
                stage0 + 1,
            )

    recurse(permutation, depth=0, prefix=0, stage0=0)
    return stages


def simulate_benes_path(
    input_wire: int,
    dimensions: List[int],
    switch_stages: List[Dict[int, int]],
) -> Tuple[List[int], List[str]]:
    wire = input_wire
    wires = [wire]
    states: List[str] = []
    for stage, dim in enumerate(dimensions):
        pair_base = min(wire, wire ^ (1 << dim))
        state = switch_stages[stage][pair_base]
        states.append("cross" if state else "straight")
        if state:
            wire ^= 1 << dim
        wires.append(wire)
    return wires, states


# ---------------------------------------------------------------------------
# Library construction
# ---------------------------------------------------------------------------


class LibraryBuilder:
    def __init__(
        self,
        leaves: List[str],
        *,
        variable_length_paths: bool = False,
        long_path_edges: int = 2,
    ):
        if not leaves:
            raise ValueError("At least one taxon is required")
        if len(set(leaves)) != len(leaves):
            raise ValueError("Taxon labels must be distinct")
        if long_path_edges < 2:
            raise ValueError("long_path_edges must be >= 2")

        self.leaves = leaves
        self.L = len(leaves)
        self.variable_length_paths = variable_length_paths
        self.long_path_edges = long_path_edges

        self.leaf_unit_by_label = {
            label: f"LEAF_{i:04d}" for i, label in enumerate(leaves)
        }

        self.units: Dict[str, Unit] = {}
        self.connections: List[Connection] = []
        self.logical_edges: List[LogicalEdge] = []
        self.internal_paths: List[InternalPath] = []
        self.internal_connections: List[InternalConnection] = []

        self.split_units_by_m: Dict[int, List[str]] = defaultdict(list)
        self.split_units_by_type: Dict[Tuple[int, int], List[str]] = defaultdict(list)
        self.requests_by_count: Dict[int, List[Tuple[str, str]]] = defaultdict(list)
        self.routers: Dict[int, dict] = {}
        self.root_candidates: List[str] = []

    # ----- unit / geometry helpers ----------------------------------------

    def add_unit(self, unit: Unit):
        if unit.unit_id in self.units:
            raise RuntimeError(f"Duplicate unit ID {unit.unit_id}")
        self.units[unit.unit_id] = unit

    def add_logical_edge(self, unit_id: str, logical_id: str, a: str, b: str, flow: int):
        le = LogicalEdge(unit_id, logical_id, a, b, flow)
        self.logical_edges.append(le)
        self._materialize_internal_paths(le)

    def _materialize_internal_paths(self, le: LogicalEdge):
        group = f"{le.unit_id}:{le.logical_edge_id}"

        def add_path(option: str, nodes: List[str]):
            self.internal_paths.append(
                InternalPath(
                    unit_id=le.unit_id,
                    logical_edge_id=le.logical_edge_id,
                    hybrid_group=group,
                    path_option=option,
                    path_length_edges=len(nodes) - 1,
                    edge_flow=le.edge_flow,
                    endpoint_a=le.endpoint_a,
                    endpoint_b=le.endpoint_b,
                    path_nodes="|".join(nodes),
                )
            )
            for i in range(len(nodes) - 1):
                self.internal_connections.append(
                    InternalConnection(
                        unit_id=le.unit_id,
                        logical_edge_id=le.logical_edge_id,
                        hybrid_group=group,
                        path_option=option,
                        segment_index=i,
                        segment_from=nodes[i],
                        segment_to=nodes[i + 1],
                        edge_flow=le.edge_flow,
                    )
                )

        add_path("short" if self.variable_length_paths else "base", [le.endpoint_a, le.endpoint_b])

        if self.variable_length_paths:
            helper_count = self.long_path_edges - 1
            helpers = [
                f"{le.logical_edge_id}__H{i+1}" for i in range(helper_count)
            ]
            add_path("long", [le.endpoint_a] + helpers + [le.endpoint_b])

    def add_connection(
        self,
        parent: str,
        port: str,
        child: str,
        option_index: int,
        flow: int,
        kind: str,
        *,
        router_count: Optional[int] = None,
        router_stage: Optional[int] = None,
        switch_id: str = "",
        switch_state: str = "",
    ):
        self.connections.append(
            Connection(
                parent_unit=parent,
                parent_port=port,
                child_unit=child,
                option_index=option_index,
                edge_flow=flow,
                parent_upper_flow=self.units[parent].upper_flow,
                child_upper_flow=self.units[child].upper_flow,
                connection_kind=kind,
                router_count=router_count,
                router_stage=router_stage,
                switch_id=switch_id,
                switch_state=switch_state,
            )
        )

    def build_leaf_units(self):
        for i, label in enumerate(self.leaves):
            uid = f"LEAF_{i:04d}"
            self.add_unit(
                Unit(
                    unit_id=uid,
                    kind="taxon",
                    descendant_count=1,
                    upper_flow=1,
                    label=label,
                    subtree_size=1,
                    notes="Unique physical taxon unit; exactly one unit for this label.",
                )
            )
            self.add_logical_edge(uid, "taxon_edge", "upper", f"taxon:{label}", 1)

    def build_split_units_and_requests(self):
        for m in range(2, self.L + 1):
            q_m = self.L // m
            for a in range(1, m // 2 + 1):
                b = m - a
                for j in range(q_m):
                    uid = f"S_m{m:04d}_a{a:04d}_j{j:04d}"
                    self.add_unit(
                        Unit(
                            unit_id=uid,
                            kind="split",
                            descendant_count=m,
                            upper_flow=m,
                            subtree_size=m,
                            split_a=a,
                            split_b=b,
                            copy_index=j,
                            left_flow=a,
                            right_flow=b,
                            notes=(
                                f"Fixed-flow Y unit: upper={m}, left={a}, right={b}. "
                                f"Multiplicity floor(L/m) is the maximum simultaneous need."
                            ),
                        )
                    )
                    self.split_units_by_m[m].append(uid)
                    self.split_units_by_type[(m, a)].append(uid)
                    self.requests_by_count[a].append((uid, "left"))
                    self.requests_by_count[b].append((uid, "right"))

                    self.add_logical_edge(uid, "upper_arm", "upper", "branch", m)
                    self.add_logical_edge(uid, "left_arm", "branch", "left", a)
                    self.add_logical_edge(uid, "right_arm", "branch", "right", b)

        # q_L = 1, so each possible root split has one physical candidate.
        self.root_candidates = list(self.split_units_by_m[self.L]) if self.L >= 2 else []

    def real_output_targets(self, m: int) -> List[str]:
        if m == 1:
            return [self.leaf_unit_by_label[label] for label in self.leaves]
        return list(self.split_units_by_m[m])

    @staticmethod
    def _switch_id(m: int, stage: int, wire: int, dim: int) -> str:
        pair_base = min(wire, wire ^ (1 << dim))
        return f"R{m}_S{stage}_P{pair_base}"

    def build_compact_router(self, m: int):
        """
        Build a compact padded Beneš router for count m.

        Boundary wire layers are not physical units:
          request port --(stage 0 choice)--> first internal layer
          last internal layer --(final stage choice)--> target unit

        For width 2, there are no router units at all: each request port has at
        most two direct target choices, coupled by the single 2x2 switch state.
        """
        requests = list(self.requests_by_count[m])
        outputs = self.real_output_targets(m)
        if not requests:
            return

        width = next_power_of_two(max(len(requests), len(outputs), 1))
        dims = benes_dimensions(width)
        D = len(dims)

        # input and output slot metadata used by target routing
        input_requests = [
            {"slot": slot, "parent_unit": p, "parent_port": port}
            for slot, (p, port) in enumerate(requests)
        ]
        output_targets = [
            {"slot": slot, "unit_id": uid}
            for slot, uid in enumerate(outputs)
        ]

        internal_unit_count = 0

        if width == 1:
            # Only one possible request and one possible destination.
            if len(requests) != 1 or len(outputs) != 1:
                raise RuntimeError("Width-1 router has nontrivial I/O")
            p, port = requests[0]
            self.add_connection(
                p, port, outputs[0], 0, m, "direct_router",
                router_count=m, router_stage=0,
            )
        elif D == 1:
            # One 2x2 switch; request ports connect directly to real targets.
            dim = dims[0]
            for inp, (parent, port) in enumerate(requests):
                sid = self._switch_id(m, 0, inp, dim)
                for state in (0, 1):
                    out = inp ^ ((1 << dim) if state else 0)
                    if out < len(outputs):
                        self.add_connection(
                            parent, port, outputs[out], state, m,
                            "compact_benes_direct",
                            router_count=m,
                            router_stage=0,
                            switch_id=sid,
                            switch_state="cross" if state else "straight",
                        )
        else:
            # Physical internal wire units only for layers after stage 0 and
            # before the final stage: compact layers 1..D-1 inclusive.
            wire_id: Dict[Tuple[int, int], str] = {}
            for layer in range(1, D):
                for wire in range(width):
                    uid = f"R_m{m:04d}_L{layer:03d}_W{wire:05d}"
                    wire_id[(layer, wire)] = uid
                    self.add_unit(
                        Unit(
                            unit_id=uid,
                            kind="router",
                            descendant_count=m,
                            upper_flow=m,
                            subtree_size=m,
                            down_flow=m,
                            router_count=m,
                            router_layer=layer,
                            router_wire=wire,
                            router_width=width,
                            notes=(
                                f"Compact count-{m} Beneš wire unit; upper/down flow={m}."
                            ),
                        )
                    )
                    self.add_logical_edge(uid, "through", "upper", "down", m)
                    internal_unit_count += 1

            # First stage directly from requesting unit ports to layer 1.
            dim0 = dims[0]
            for inp, (parent, port) in enumerate(requests):
                sid = self._switch_id(m, 0, inp, dim0)
                for state in (0, 1):
                    w1 = inp ^ ((1 << dim0) if state else 0)
                    child = wire_id[(1, w1)]
                    self.add_connection(
                        parent, port, child, state, m,
                        "compact_benes_enter",
                        router_count=m,
                        router_stage=0,
                        switch_id=sid,
                        switch_state="cross" if state else "straight",
                    )

            # Middle stages: layer s -> layer s+1 for s=1..D-2.
            for stage in range(1, D - 1):
                dim = dims[stage]
                for wire in range(width):
                    parent = wire_id[(stage, wire)]
                    sid = self._switch_id(m, stage, wire, dim)
                    same = wire_id[(stage + 1, wire)]
                    cross = wire_id[(stage + 1, wire ^ (1 << dim))]
                    self.add_connection(
                        parent, "down", same, 0, m,
                        "compact_benes_straight",
                        router_count=m,
                        router_stage=stage,
                        switch_id=sid,
                        switch_state="straight",
                    )
                    self.add_connection(
                        parent, "down", cross, 1, m,
                        "compact_benes_cross",
                        router_count=m,
                        router_stage=stage,
                        switch_id=sid,
                        switch_state="cross",
                    )

            # Final stage directly into real destination units.
            final_stage = D - 1
            dimf = dims[final_stage]
            for wire in range(width):
                parent = wire_id[(D - 1, wire)]
                sid = self._switch_id(m, final_stage, wire, dimf)
                for state in (0, 1):
                    out = wire ^ ((1 << dimf) if state else 0)
                    if out < len(outputs):
                        self.add_connection(
                            parent, "down", outputs[out], state, m,
                            "compact_benes_exit",
                            router_count=m,
                            router_stage=final_stage,
                            switch_id=sid,
                            switch_state="cross" if state else "straight",
                        )

        self.routers[m] = {
            "count": m,
            "real_inputs": len(requests),
            "real_outputs": len(outputs),
            "padded_width": width,
            "switch_stages": D,
            "dimensions": dims,
            "internal_wire_layers": max(D - 1, 0),
            "router_units": internal_unit_count,
            "input_requests": input_requests,
            "output_targets": output_targets,
        }

    def build_routers(self):
        # No count-L router: root is selected directly from root_candidates.
        for m in range(1, self.L):
            self.build_compact_router(m)

    def build(self):
        self.build_leaf_units()
        if self.L >= 2:
            self.build_split_units_and_requests()
            self.build_routers()
        self.verify()

    # ----- verification ----------------------------------------------------

    def port_flow(self, unit: Unit, port: str) -> int:
        if unit.kind == "split":
            if port == "left":
                assert unit.left_flow is not None
                return unit.left_flow
            if port == "right":
                assert unit.right_flow is not None
                return unit.right_flow
        if unit.kind == "router" and port == "down":
            assert unit.down_flow is not None
            return unit.down_flow
        raise ValueError(f"Unknown lower port {unit.unit_id}:{port}")

    def verify(self):
        # Exactly one unit per supplied taxon.
        taxon_units = [u for u in self.units.values() if u.kind == "taxon"]
        if len(taxon_units) != self.L:
            raise AssertionError("Taxon unit count is not exactly L")
        labels = [u.label for u in taxon_units]
        if sorted(labels) != sorted(self.leaves) or len(set(labels)) != self.L:
            raise AssertionError("Taxon units are not one-to-one with supplied labels")

        # Unit fixed-flow invariants.
        for u in self.units.values():
            if u.upper_flow != u.descendant_count:
                raise AssertionError(f"Upper flow mismatch on {u.unit_id}")
            if u.kind == "split":
                assert u.left_flow is not None and u.right_flow is not None
                if u.left_flow + u.right_flow != u.upper_flow:
                    raise AssertionError(f"Flow conservation failed on {u.unit_id}")
            elif u.kind == "router":
                if u.down_flow != u.upper_flow:
                    raise AssertionError(f"Router flow mismatch on {u.unit_id}")
            elif u.kind == "taxon":
                if u.upper_flow != 1:
                    raise AssertionError(f"Taxon flow mismatch on {u.unit_id}")

        # Every logical edge has positive known flow and materialized path(s).
        if any(le.edge_flow <= 0 for le in self.logical_edges):
            raise AssertionError("Nonpositive logical edge flow")
        paths_by_logical: Dict[Tuple[str, str], List[InternalPath]] = defaultdict(list)
        for p in self.internal_paths:
            paths_by_logical[(p.unit_id, p.logical_edge_id)].append(p)
        for le in self.logical_edges:
            ps = paths_by_logical[(le.unit_id, le.logical_edge_id)]
            expected = 2 if self.variable_length_paths else 1
            if len(ps) != expected:
                raise AssertionError(f"Wrong path alternative count for {le}")
            if any(p.edge_flow != le.edge_flow for p in ps):
                raise AssertionError(f"Internal path flow mismatch for {le}")

        # k <= 2 and inter-unit edge flow exactly matches both lower port and
        # child upper flow.
        by_port: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        for c in self.connections:
            if c.parent_unit not in self.units or c.child_unit not in self.units:
                raise AssertionError("Connection references unknown unit")
            by_port[(c.parent_unit, c.parent_port)].add(c.child_unit)
            if c.child_upper_flow != c.edge_flow:
                raise AssertionError(f"Child upper flow mismatch on connection {c}")
            expected = self.port_flow(self.units[c.parent_unit], c.parent_port)
            if expected != c.edge_flow:
                raise AssertionError(f"Parent port flow mismatch on connection {c}")
        for key, dests in by_port.items():
            if len(dests) > 2:
                raise AssertionError(f"k>2 at {key}: {sorted(dests)}")

    # ----- summaries -------------------------------------------------------

    def max_k(self) -> int:
        d: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        for c in self.connections:
            d[(c.parent_unit, c.parent_port)].add(c.child_unit)
        return max((len(v) for v in d.values()), default=0)

    def summary(self) -> dict:
        kinds: Dict[str, int] = defaultdict(int)
        for u in self.units.values():
            kinds[u.kind] += 1
        denom = self.L * self.L * max(math.log2(max(self.L, 2)), 1.0)
        return {
            "leaf_count": self.L,
            "total_units": len(self.units),
            "units_by_kind": dict(kinds),
            "split_units": kinds.get("split", 0),
            "router_units": kinds.get("router", 0),
            "taxon_units": kinds.get("taxon", 0),
            "total_connections": len(self.connections),
            "max_k": self.max_k(),
            "root_candidate_count": len(self.root_candidates),
            "variable_length_paths": self.variable_length_paths,
            "long_path_edges": self.long_path_edges if self.variable_length_paths else None,
            "logical_internal_edges": len(self.logical_edges),
            "internal_path_alternatives": len(self.internal_paths),
            "internal_connection_segments": len(self.internal_connections),
            "normalized_units_over_L2log2L": len(self.units) / denom if denom else 0.0,
        }

    # ----- target assignment / routing ------------------------------------

    def assign_target_units(self, root: TargetNode, nodes: List[TargetNode]):
        target_labels = sorted(n.label for n in nodes if n.is_leaf and n.label is not None)
        if target_labels != sorted(self.leaves):
            raise ValueError(
                "Target tree leaf labels must match the supplied library taxa exactly"
            )
        if root.descendant_count != self.L:
            raise ValueError("Target tree leaf count does not match library")

        next_copy: Dict[Tuple[int, int], int] = defaultdict(int)
        for node in nodes:
            if node.is_leaf:
                assert node.label is not None
                node.assigned_unit = self.leaf_unit_by_label[node.label]
                continue

            assert node.left is not None and node.right is not None
            m = node.descendant_count
            a = node.left.descendant_count
            b = node.right.descendant_count
            if a > b:
                raise AssertionError("Target child orientation is not canonical")
            pool = self.split_units_by_type[(m, a)]
            idx = next_copy[(m, a)]
            if idx >= len(pool):
                raise RuntimeError(
                    f"Insufficient physical copies for split {m}={a}+{b}; "
                    f"needed > {len(pool)}"
                )
            node.assigned_unit = pool[idx]
            next_copy[(m, a)] += 1

        if root.assigned_unit not in self.root_candidates:
            raise AssertionError("Assigned root is not a root candidate")

    def _build_demands(self, root: TargetNode, nodes: List[TargetNode]) -> Dict[int, List[dict]]:
        demands: Dict[int, List[dict]] = defaultdict(list)
        route_counter = 0
        for node in nodes:
            if node.is_leaf:
                continue
            assert node.left is not None and node.right is not None
            for port, child in (("left", node.left), ("right", node.right)):
                route_id = f"Q{route_counter:05d}"
                route_counter += 1
                m = child.descendant_count
                demands[m].append(
                    {
                        "route_id": route_id,
                        "parent_target_node": node.node_id,
                        "child_target_node": child.node_id,
                        "parent_unit": node.assigned_unit,
                        "parent_port": port,
                        "child_unit": child.assigned_unit,
                        "requested_count": m,
                        "edge_flow": m,
                    }
                )
        return demands

    def _connection_exists(
        self, parent: str, port: str, child: str, state: Optional[str] = None
    ) -> bool:
        for c in self.connections:
            if c.parent_unit == parent and c.parent_port == port and c.child_unit == child:
                if state is None or c.switch_state == state or not c.switch_state:
                    return True
        return False

    def route_target(self, root: TargetNode, nodes: List[TargetNode]) -> dict:
        self.assign_target_units(root, nodes)
        demands_by_count = self._build_demands(root, nodes)

        active_connections: List[dict] = []
        route_records: List[dict] = []
        router_switches: Dict[int, List[dict]] = {}
        active_unit_target_nodes: Dict[str, Set[str]] = defaultdict(set)
        active_unit_routes: Dict[str, Set[str]] = defaultdict(set)

        for node in nodes:
            active_unit_target_nodes[node.assigned_unit].add(node.node_id)

        for m in range(1, self.L):
            router = self.routers.get(m)
            if router is None:
                if demands_by_count.get(m):
                    raise RuntimeError(f"No router exists for active count {m}")
                continue

            width = router["padded_width"]
            dims = list(router["dimensions"])
            D = len(dims)

            input_slot_by_request = {
                (r["parent_unit"], r["parent_port"]): r["slot"]
                for r in router["input_requests"]
            }
            output_slot_by_unit = {
                r["unit_id"]: r["slot"] for r in router["output_targets"]
            }

            active_map: Dict[int, int] = {}
            demand_by_input: Dict[int, dict] = {}
            for demand in demands_by_count.get(m, []):
                key = (demand["parent_unit"], demand["parent_port"])
                if key not in input_slot_by_request:
                    raise RuntimeError(f"R_{m} lacks input request {key}")
                if demand["child_unit"] not in output_slot_by_unit:
                    raise RuntimeError(f"R_{m} lacks output {demand['child_unit']}")
                inp = input_slot_by_request[key]
                out = output_slot_by_unit[demand["child_unit"]]
                if inp in active_map:
                    raise RuntimeError(f"Duplicate active input R_{m}:{inp}")
                if out in active_map.values():
                    raise RuntimeError(f"Duplicate active output R_{m}:{out}")
                active_map[inp] = out
                demand_by_input[inp] = demand

            permutation = [-1] * width
            for inp, out in active_map.items():
                permutation[inp] = out
            unused_inputs = [i for i, x in enumerate(permutation) if x < 0]
            used_outputs = set(active_map.values())
            unused_outputs = [o for o in range(width) if o not in used_outputs]
            for inp, out in zip(unused_inputs, unused_outputs):
                permutation[inp] = out

            switches = benes_route_switches(permutation)
            switch_rows: List[dict] = []
            for stage, states in enumerate(switches):
                dim = dims[stage]
                for pair_base in sorted(states):
                    switch_rows.append(
                        {
                            "stage": stage,
                            "dimension": dim,
                            "pair_base_wire": pair_base,
                            "switch_id": f"R{m}_S{stage}_P{pair_base}",
                            "state": "cross" if states[pair_base] else "straight",
                            "edge_flow": m,
                        }
                    )
            router_switches[m] = switch_rows

            for inp in sorted(demand_by_input):
                demand = demand_by_input[inp]
                route_id = demand["route_id"]
                wires, states = simulate_benes_path(inp, dims, switches)
                expected_out = active_map[inp]
                if wires[-1] != expected_out:
                    raise RuntimeError(
                        f"Beneš route failure R_{m}: input {inp}->{wires[-1]}, expected {expected_out}"
                    )

                internal_units: List[str] = []
                if D >= 2:
                    internal_units = [
                        f"R_m{m:04d}_L{layer:03d}_W{wires[layer]:05d}"
                        for layer in range(1, D)
                    ]
                    for uid in internal_units:
                        active_unit_routes[uid].add(route_id)

                path_units = [demand["parent_unit"]] + internal_units + [demand["child_unit"]]

                # Materialize active compact connections stage-by-stage.
                if D == 0:
                    # width 1 direct router
                    parent = demand["parent_unit"]
                    child = demand["child_unit"]
                    if not self._connection_exists(parent, demand["parent_port"], child):
                        raise AssertionError("Missing direct compact connection")
                    active_connections.append(
                        self._active_connection_record(
                            route_id, parent, demand["parent_port"], child, 0,
                            "direct_router", m, 0, "fixed"
                        )
                    )
                elif D == 1:
                    state = states[0]
                    option = 1 if state == "cross" else 0
                    parent = demand["parent_unit"]
                    child = demand["child_unit"]
                    if not self._connection_exists(parent, demand["parent_port"], child, state):
                        raise AssertionError("Missing compact direct switch connection")
                    active_connections.append(
                        self._active_connection_record(
                            route_id, parent, demand["parent_port"], child, option,
                            "compact_benes_direct", m, 0, state
                        )
                    )
                else:
                    # stage 0: parent -> layer 1
                    state = states[0]
                    option = 1 if state == "cross" else 0
                    active_connections.append(
                        self._active_connection_record(
                            route_id,
                            demand["parent_unit"],
                            demand["parent_port"],
                            internal_units[0],
                            option,
                            "compact_benes_enter",
                            m,
                            0,
                            state,
                        )
                    )

                    # middle stages 1..D-2
                    for stage in range(1, D - 1):
                        state = states[stage]
                        option = 1 if state == "cross" else 0
                        active_connections.append(
                            self._active_connection_record(
                                route_id,
                                internal_units[stage - 1],
                                "down",
                                internal_units[stage],
                                option,
                                "compact_benes_cross" if option else "compact_benes_straight",
                                m,
                                stage,
                                state,
                            )
                        )

                    # final stage D-1: last internal -> target
                    state = states[-1]
                    option = 1 if state == "cross" else 0
                    active_connections.append(
                        self._active_connection_record(
                            route_id,
                            internal_units[-1],
                            "down",
                            demand["child_unit"],
                            option,
                            "compact_benes_exit",
                            m,
                            D - 1,
                            state,
                        )
                    )

                route_records.append(
                    {
                        **demand,
                        "total_descendant_leaves": m,
                        "router_input_slot": inp,
                        "router_output_slot": expected_out,
                        "router_width": width,
                        "router_switch_stages": D,
                        "switch_states": "|".join(states),
                        "internal_router_units": "|".join(internal_units),
                        "full_unit_path": "|".join(path_units),
                    }
                )

        active_ids = set(active_unit_target_nodes) | set(active_unit_routes)
        active_units = [
            self._active_unit_record(
                uid,
                active_unit_target_nodes.get(uid, set()),
                active_unit_routes.get(uid, set()),
            )
            for uid in sorted(active_ids)
        ]

        return {
            "selected_root_unit": root.assigned_unit,
            "target_nodes": [self._target_node_record(n) for n in nodes],
            "routes": route_records,
            "active_connections": active_connections,
            "active_units": active_units,
            "router_switches": router_switches,
        }

    def _target_node_record(self, node: TargetNode) -> dict:
        if node.is_leaf:
            a = b = ""
        else:
            assert node.left is not None and node.right is not None
            a = node.left.descendant_count
            b = node.right.descendant_count
        return {
            "target_node_id": node.node_id,
            "kind": "taxon" if node.is_leaf else "internal",
            "label": node.label or "",
            "descendant_leaf_count": node.descendant_count,
            "descendant_leaf_labels": "|".join(node.descendant_labels),
            "split_a": a,
            "split_b": b,
            "assigned_unit": node.assigned_unit,
            "assigned_unit_upper_flow": self.units[node.assigned_unit].upper_flow,
        }

    def _active_unit_record(self, uid: str, target_nodes: Set[str], route_ids: Set[str]) -> dict:
        u = self.units[uid]
        return {
            "unit_id": uid,
            "kind": u.kind,
            "label": u.label,
            "total_descendant_leaves": u.descendant_count,
            "upper_flow": u.upper_flow,
            "left_flow": u.left_flow if u.left_flow is not None else "",
            "right_flow": u.right_flow if u.right_flow is not None else "",
            "down_flow": u.down_flow if u.down_flow is not None else "",
            "target_node_ids": "|".join(sorted(target_nodes)),
            "route_ids": "|".join(sorted(route_ids)),
        }

    def _active_connection_record(
        self,
        route_id: str,
        parent: str,
        port: str,
        child: str,
        option: int,
        kind: str,
        flow: int,
        stage: int,
        state: str,
    ) -> dict:
        return {
            "route_id": route_id,
            "parent_unit": parent,
            "parent_port": port,
            "child_unit": child,
            "option_index": option,
            "connection_kind": kind,
            "edge_flow": flow,
            "parent_upper_flow": self.units[parent].upper_flow,
            "child_upper_flow": self.units[child].upper_flow,
            "router_count": flow,
            "router_stage": stage,
            "switch_state": state,
        }

    # ----- output ----------------------------------------------------------

    @staticmethod
    def _write_rows(path: Path, rows: List[dict], default_fields: List[str]):
        with path.open("w", newline="", encoding="utf-8") as f:
            fields = list(rows[0].keys()) if rows else default_fields
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    def write_library(self, out_dir: Path) -> List[Path]:
        out_dir.mkdir(parents=True, exist_ok=True)

        units_path = out_dir / "units.csv"
        connections_path = out_dir / "connections.csv"
        logical_path = out_dir / "unit_logical_edges.csv"
        internal_paths_path = out_dir / "unit_internal_paths.csv"
        internal_conn_path = out_dir / "unit_internal_connections.csv"
        json_path = out_dir / "library.json"
        summary_path = out_dir / "summary.txt"

        unit_rows = [asdict(self.units[uid]) for uid in sorted(self.units)]
        conn_rows = [asdict(c) for c in self.connections]
        logical_rows = [asdict(x) for x in self.logical_edges]
        ipath_rows = [asdict(x) for x in self.internal_paths]
        iconn_rows = [asdict(x) for x in self.internal_connections]

        self._write_rows(units_path, unit_rows, list(asdict(Unit("", "", 1, 1)).keys()))
        self._write_rows(
            connections_path, conn_rows,
            list(asdict(Connection("", "", "", 0, 1, 1, 1, "")).keys()),
        )
        self._write_rows(logical_path, logical_rows, list(asdict(LogicalEdge("", "", "", "", 1)).keys()))
        self._write_rows(
            internal_paths_path, ipath_rows,
            list(asdict(InternalPath("", "", "", "", 1, 1, "", "", "")).keys()),
        )
        self._write_rows(
            internal_conn_path, iconn_rows,
            list(asdict(InternalConnection("", "", "", "", 0, "", "", 1)).keys()),
        )

        payload = {
            "description": (
                "Compact O(L^2 log L), k<=2 universal unit library with unique taxa, "
                "fixed unit upper flow, and fixed flow on every logical/internal/inter-unit edge."
            ),
            "taxa": self.leaves,
            "root_candidates": self.root_candidates,
            "summary": self.summary(),
            "routers": self.routers,
            "units": {uid: asdict(u) for uid, u in self.units.items()},
            "connections": conn_rows,
            "logical_edges": logical_rows,
            "internal_paths": ipath_rows,
            "internal_connections": iconn_rows,
        }
        json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        s = self.summary()
        lines = [
            f"L = {self.L}",
            f"total units = {s['total_units']}",
            f"taxon units = {s['taxon_units']} (exactly one per supplied taxon)",
            f"split units = {s['split_units']}",
            f"router units = {s['router_units']}",
            f"total inter-unit connections = {s['total_connections']}",
            f"max k = {s['max_k']}",
            f"root candidates = {s['root_candidate_count']}",
            f"variable-length internal paths = {s['variable_length_paths']}",
            f"units/(L^2 log2 L) = {s['normalized_units_over_L2log2L']:.6f}",
            "",
            "Every unit has upper_flow == descendant_count.",
            "Every connection row has a fixed edge_flow.",
            "Every internal logical edge and every alternative path segment has a fixed edge_flow.",
        ]
        summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        return [
            units_path,
            connections_path,
            logical_path,
            internal_paths_path,
            internal_conn_path,
            json_path,
            summary_path,
        ]

    def write_routing(
        self,
        out_dir: Path,
        routing: dict,
        tree_text: str,
        *,
        active_length_choice: str = "short",
    ) -> List[Path]:
        target_nodes_path = out_dir / "target_nodes.csv"
        routes_path = out_dir / "routes.csv"
        active_units_path = out_dir / "active_units.csv"
        active_conn_path = out_dir / "active_connections.csv"
        switches_path = out_dir / "router_switches.csv"
        active_ipaths_path = out_dir / "active_internal_paths.csv"
        active_iconn_path = out_dir / "active_internal_connections.csv"
        routing_json_path = out_dir / "routing.json"
        routing_summary_path = out_dir / "routing_summary.txt"

        self._write_rows(target_nodes_path, routing["target_nodes"], ["target_node_id"])
        self._write_rows(routes_path, routing["routes"], ["route_id"])
        self._write_rows(active_units_path, routing["active_units"], ["unit_id"])
        self._write_rows(active_conn_path, routing["active_connections"], ["route_id"])

        switch_rows = []
        for m, rows in sorted(routing["router_switches"].items()):
            for row in rows:
                switch_rows.append({"router_count": m, **row})
        self._write_rows(switches_path, switch_rows, ["router_count", "stage", "state"])

        active_ids = {r["unit_id"] for r in routing["active_units"]}
        if self.variable_length_paths:
            chosen = active_length_choice
        else:
            chosen = "base"
        active_ipaths = [
            asdict(p) for p in self.internal_paths
            if p.unit_id in active_ids and p.path_option == chosen
        ]
        active_iconns = [
            asdict(c) for c in self.internal_connections
            if c.unit_id in active_ids and c.path_option == chosen
        ]
        self._write_rows(active_ipaths_path, active_ipaths, ["unit_id", "path_option"])
        self._write_rows(active_iconn_path, active_iconns, ["unit_id", "path_option"])

        payload = {
            "target_newick": tree_text,
            "selected_root_unit": routing["selected_root_unit"],
            "active_length_choice": chosen,
            **routing,
            "active_internal_paths": active_ipaths,
            "active_internal_connections": active_iconns,
        }
        routing_json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        lines = [
            f"target = {tree_text}",
            f"L = {self.L}",
            f"selected root unit = {routing['selected_root_unit']}",
            f"active units = {len(routing['active_units'])}",
            f"active inter-unit connections = {len(routing['active_connections'])}",
            f"target parent->child routes = {len(routing['routes'])}",
            f"active internal path choice = {chosen}",
            "",
            "Every active connection has fixed edge_flow equal to the child subtree's descendant count.",
        ]
        routing_summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        return [
            target_nodes_path,
            routes_path,
            active_units_path,
            active_conn_path,
            switches_path,
            active_ipaths_path,
            active_iconn_path,
            routing_json_path,
            routing_summary_path,
        ]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Build a compact O(L^2 log L), k<=2 universal subtree-unit library "
            "with one unit per taxon and fixed flow on every unit edge."
        )
    )
    p.add_argument("leaves", nargs="*", help="Distinct taxon labels, e.g. A B C D")
    p.add_argument(
        "--leaf-count", type=int, default=None,
        help="Generate taxa leaf_1 ... leaf_L instead of listing labels.",
    )
    tg = p.add_mutually_exclusive_group()
    tg.add_argument("--tree", default=None, help="Rooted binary Newick tree to route")
    tg.add_argument("--tree-file", default=None, help="File containing one rooted binary Newick tree")
    p.add_argument("--out-dir", default="fixed_flow_units", help="Output directory")
    p.add_argument(
        "--variable-length-paths", action="store_true",
        help="Add a second longer A-B path for every logical edge inside every unit.",
    )
    p.add_argument(
        "--long-path-edges", type=int, default=2,
        help="Number of edges in each optional long internal A-B path (>=2; default 2).",
    )
    p.add_argument(
        "--active-length-choice", choices=["short", "long"], default="short",
        help="When routing a target and variable paths are enabled, which internal path to mark active.",
    )
    return p.parse_args()


def main():
    args = parse_args()

    tree_text: Optional[str] = None
    parsed_root: Optional[TargetNode] = None
    if args.tree_file:
        tree_text = Path(args.tree_file).read_text(encoding="utf-8").strip()
    elif args.tree:
        tree_text = args.tree.strip()

    if tree_text:
        parsed_root = NewickParser(tree_text).parse()
        tree_leaves = leaf_labels_in_order(parsed_root)
        if len(set(tree_leaves)) != len(tree_leaves):
            raise SystemExit("Target tree contains duplicate taxon labels")
    else:
        tree_leaves = []

    if args.leaf_count is not None:
        if args.leaves:
            raise SystemExit("Use either positional labels or --leaf-count, not both")
        if args.leaf_count < 1:
            raise SystemExit("--leaf-count must be >=1")
        leaves = [f"leaf_{i+1}" for i in range(args.leaf_count)]
    elif args.leaves:
        leaves = args.leaves
    elif tree_text:
        leaves = tree_leaves
    else:
        raise SystemExit("Provide taxon labels, --leaf-count, or --tree/--tree-file")

    if tree_text and sorted(leaves) != sorted(tree_leaves):
        raise SystemExit("Supplied taxon labels do not match target tree labels")

    builder = LibraryBuilder(
        leaves,
        variable_length_paths=args.variable_length_paths,
        long_path_edges=args.long_path_edges,
    )
    builder.build()

    out_dir = Path(args.out_dir)
    lib_paths = builder.write_library(out_dir)

    routing_paths: List[Path] = []
    if parsed_root is not None and tree_text is not None:
        nodes = annotate_target_tree(parsed_root)
        routing = builder.route_target(parsed_root, nodes)
        routing_paths = builder.write_routing(
            out_dir,
            routing,
            tree_text,
            active_length_choice=args.active_length_choice,
        )

    s = builder.summary()
    print(f"Built compact fixed-flow universal library for L={builder.L}")
    print(f"  units:        {s['total_units']}")
    print(f"  taxon units:  {s['taxon_units']} (exactly one per taxon)")
    print(f"  split units:  {s['split_units']}")
    print(f"  router units: {s['router_units']}")
    print(f"  max k:        {s['max_k']}")
    print(f"  root choices: {s['root_candidate_count']}")
    print(f"  variable internal lengths: {s['variable_length_paths']}")
    if parsed_root is not None:
        print("  target tree routed successfully")
    print("Outputs:")
    for path in lib_paths + routing_paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
