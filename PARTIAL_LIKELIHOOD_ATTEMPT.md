# Partial-likelihood recompute for SPR moves — fixed and shipped

## 2026-08-12 retry: root cause found and implementation restored

The unresolved corruption was not caused by reuse of the six scratch
buffers. The six-direction invalidation set itself was incomplete. A distant
SPR changes the subtree represented by one directed partial on every
persistent edge between the removal and insertion points, even though those
edges are not physically rewired. Leaving those path partials valid explains
both the plausible small discrepancies and the occasionally enormous false
improvements seen below.

The restored implementation in `tree/spr_topology_test.cpp` now swaps fresh
scratch buffers onto the six rewired directions, `prune_node->prune_dad`, and
the affected direction of every intermediate path edge. Candidate batches
are evaluated locally and restored exactly; rejected winners roll back
without disturbing the valid baseline cache. A locally improving winner is
still recomputed once from scratch before commit, providing a correctness
boundary and a fully populated baseline for the next step. Branch-length
reoptimization retains the full-reset path.

The corrected deferred verifier (which evaluates the entire local batch
before doing any full-reset reference work) reported zero mismatches in 160
fast-mode comparisons and zero mismatches in 38 exhaustive-mode comparisons
on `sim.treefile`, in addition to shorter follow-up runs after final cleanup.

### Speed and likelihood comparison on `sim.treefile`

The final implementation was benchmarked on the repository's demo dataset
(`sim.treefile` plus `sim.fa`: 100 DNA sequences, 10,000 sites). The hidden
`verifylocal` diagnostic scores each candidate locally first, then scores the
exact same move using the previous full `deleteAllPartialLh()` +
`initializeAllPartialLh()` method. Its timers cover candidate scoring itself
(including apply/rollback and local path discovery) and use IQ-TREE's
`getCPUTime()` clock; candidate enumeration, winner commit, and reporting are
outside both timers.

| Search command suffix | Candidates | Local CPU s | Full-reset CPU s | Speedup | Max abs. logL difference | Mean abs. difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `6 20 fast 20 quiet verifylocal` | 400 | 0.578 | 10.469 | 18.108x | 2.328306e-10 | 2.299203e-11 |
| `6 10 quiet verifylocal` (exhaustive) | 386 | 0.469 | 10.125 | 21.600x | 6.984919e-10 | 2.111159e-11 |

Both runs reported **zero mismatches** at the verifier's 1e-3 tolerance.
Across all 786 paired evaluations, the largest observed difference was
6.98e-10 log-likelihood units. That is many orders of magnitude below a
meaningful likelihood change and is consistent with floating-point operation
ordering. In particular, there were no false improvements like the
hundreds-of-thousands-of-units errors from the reverted six-direction
attempt. On this demo, partial recomputation is therefore about 18-22 times
faster per candidate without evidence of likelihood degradation. These
numbers apply to non-`reopt` scoring; `reopt` deliberately retains the
full-reset implementation.

Commands used to produce these numbers (from the repository root; the
`verifylocal` flag itself has since been removed from the tool now that the
fix is verified — these are kept for reference, not to be re-run as-is):

```text
build/spr_topology_test.exe --hillclimb sim.treefile 6 20 fast 20 quiet verifylocal
build/spr_topology_test.exe --hillclimb sim.treefile 6 10 quiet verifylocal
```

Date: 2026-08-12. This is the third documented attempt at replacing
`resetLikelihoodBuffers()` (full `deleteAllPartialLh()`+`initializeAllPartialLh()`,
O(whole tree) every call) with a LOCAL recompute that only touches the
`PhyloNeighbor` directions an SPR move actually changes. The first two
attempts are described in comments still in `tree/spr_topology_test.cpp`
(search for "abandoned" near `chooseGraft`/`scoreTrialSPRMove` — both
crashed on a `reorientPartialLh` assertion under `LM_PER_NODE`). This one
did NOT crash on its first pass — it failed silently, producing wrong (but
plausible-looking) log-likelihood values, which is worse than a crash, so
that first pass was reverted rather than shipped (its post-mortem is
preserved below, unmodified, since the internals it documents are still
accurate and the root-cause hunt is genuinely useful reading). **The same
day, a follow-up pass (see "2026-08-12 retry" above) found and fixed the
actual bug — an incomplete invalidation set, not buffer reuse — and that
fix is what's in `tree/spr_topology_test.cpp` today.** The `verifylocal`
diagnostic flag used to produce the benchmark numbers above has since been
removed from the tool now that the fix is verified; the local-recompute
path itself is unconditional (no flag needed) whenever branch-length
reoptimization isn't requested.

## Why this was attempted

Every `scoreTrialSPRMove` call (the hot path of `--hillclimb`) was paying
for 2-4 full-tree buffer resets (`resetLikelihoodBuffers`) plus a full
`computeLikelihood()`, per candidate, regardless of how local the SPR move
was. On `misofproteinalignment.nex` (144 taxa, 595k sites, 585,742 distinct
patterns — ~98% unique, almost no pattern compression) this made every
step take 20-30 seconds. The user asked to "try re-adding" local recompute,
using IQ-TREE's own technique as a basis.

## What was learned about IQ-TREE's internals (still true, still useful)

- `PhyloTree::applySPR`/`rollbackSPR` (`tree/phylotree.cpp`) do **not**
  touch `partial_lh_computed` at all — they only repoint `Neighbor::node`/
  `length`. Any partial-likelihood invalidation is entirely the caller's
  responsibility. `SPRRollback::clear_all_partial_lh` is set `true` by
  `applySPR` but is **never read anywhere** — a dead hint.
- `PhyloTree::computeLikelihood()` does not root-anchor at `tree.root` —
  it anchors at a lazily-initialized, then **permanently fixed**
  `current_it`/`current_it_back` pair (`phylotree.cpp:1264-1277`,
  `findFarthestLeaf()`), reused for the lifetime of the tree. A single
  anchored computation (`computeLikelihoodBranch`) only resolves **one**
  direction of each non-anchor edge (whichever is "away from" the anchor)
  — never both. This means a single `computeLikelihood()`/
  `computeLikelihoodBranch()` call after invalidating an SPR move's 6
  directions leaves **one of the six** (`regraft_node`'s own direction
  back toward `prune_dad`) still genuinely dirty. Fix used here: two
  explicit calls, `computeLikelihoodBranch(prune_dad->findNeighbor(dad2),
  prune_dad)` then `computeLikelihoodBranch(prune_dad->findNeighbor(node2),
  prune_dad)` — the second is cheap since it reuses almost everything the
  first call already computed. This part of the reasoning checked out in
  testing (the `ASSERT(get_partial_lh_computed())` guard never fired).
- `PhyloTree::reorientPartialLh` (`LM_PER_NODE` memory-saving mode) only
  tries to "steal" a buffer from a sibling direction when `partial_lh ==
  nullptr` (`if (dad_branch->partial_lh) return;`). The two earlier
  crashing attempts both, at some point, left `partial_lh` null for two
  directions at the same node simultaneously (the node's small buffer
  pool assumes at most one "away and uncomputed" direction at a time,
  matching NNI's own access pattern — SPR routinely needs two at
  `prune_dad`). **A direction whose `partial_lh` is swapped to a
  dedicated, independently-`aligned_alloc`'d buffer and never set to
  null never triggers this steal path at all**, regardless of what
  state neighboring directions are in. This is the part that actually
  worked: no crash was ever reproduced in this attempt, across dozens of
  runs including 60-step "fast" runs.
- IQ-TREE's own `PhyloTree::swapSPR`/`optimizeSPR` (`phylotree.cpp`,
  *not* `_old`) already do something similar — `newPartialLh()` +
  substituting a fresh buffer directly, bypassing the shared pool. But
  `phylotree.h` labels them explicitly `LEGACY PROTOTYPE: used only by
  the dormant optimizeSPR()/swapSPR() search` — they are not exercised
  by any live code path, and `swapSPR`'s own buffer-restore bookkeeping
  looked suspect on inspection (a `vector<double*>` pre-sized *and*
  `push_back`'d into, meaning the first N entries are never actually
  read back). Do not trust that code as a reference implementation
  without independently re-deriving it, as was done here.
- `getBestNNIForBran` (the actually-live NNI code path) uses a
  *different*, heavier mechanism: it swaps in brand-new `Neighbor`
  objects (`newNeighbor()`) with buffers carved out of a dedicated
  `nni_partial_lh`/`nni_scale_num` array sized for exactly `IT_NUM` (2 or
  6) directions, tracked via `mem_slots.addSpecialNei`/`eraseSpecialNei`.
  This attempt used a lighter-weight alternative — swapping the 4 fields
  of the *existing* `PhyloNeighbor` object in place, via a small new
  public method — rather than replacing the object. That part appears to
  have been sound (see below); the bug was elsewhere.

## The actual (unresolved) bug

Confirmed via direct instrumentation (temporarily comparing the "best
candidate" score recorded during the exhaustive-mode candidate loop
against an independent, freshly-computed reference for that exact same
move, immediately after the loop): **the winning candidate's own score,
as tracked during a sequence of several back-to-back local trial
evaluations reusing the same 6 scratch buffers, was already wrong before
any "commit" logic ran at all.** Example from one run: local-path
tracked score `-100090.1`, independent fresh recompute of the identical
move `-728309.0` — not a rounding-level discrepancy, a completely
different (implausibly high, i.e. "better than the true tree") value.

This rules out several hypotheses that were investigated and rejected in
the process:

- **Not** a commit/keep-path-specific bug — plain "fast" mode
  (`numCandidates=1`, i.e. exactly one trial per step, immediately
  followed by commit-or-reject) broke too, and a per-step check (reset
  to a clean baseline after *every* step, both accepted and rejected)
  showed the *first* accepted step could be fine while a *later*
  accepted step (after several rejects) was wrong, and vice versa in a
  different run — no fixed pattern tied to accept/reject.
- **Not** node-identity coincidences (e.g. `regraft_dad` happening to
  equal one of the pre-move siblings) — observed both with and without
  such overlaps, on both correct and incorrect steps; no correlation.
- **Not** `current_it`/`current_it_back` staleness — the diagnostic
  check itself always called `resetLikelihoodBuffers()` before its own
  reference `computeLikelihood()`, which forces `current_it`'s own
  flags to dirty regardless of history, so a stale anchor can't explain
  matches/mismatches in that specific check. (Whether `current_it`
  drifting to point at a since-repurposed `Neighbor` object is a
  *latent*, pre-existing correctness hazard in the tool more broadly —
  independent of this attempt — was not resolved; the OLD code sidesteps
  it by always fully resetting before every `computeLikelihood()` call,
  so it likely doesn't matter today, but flagging it here in case it
  ever does.)
- The `ASSERT(get_partial_lh_computed())` guard in the "keep" path never
  fired — so it is not simply "one of the six directions was never
  actually computed."

The leading remaining hypothesis, not yet tested, is that something
about **reusing the same 6 physical scratch buffers across a sequence of
several `beginLocalSPRInvalidation`/`...Discard` cycles within one step**
(exhaustive mode evaluates many candidates per step, all sharing one
`prune_dad`/sibling pair but different `regraft_dad`/`regraft_node`) —
possibly related to `computeLikelihoodBranch`'s internal traversal
scheduling (`computeTraversalInfo`'s neighbor-sorts-by-`size` heuristic,
or some other shared/scratch array touched by back-to-back
`computeLikelihoodBranch` calls, e.g. `_pattern_lh`/`theta_all`/
`buffer_partial_lh`-style tree-level scratch members that `phylotree.cpp`
declares alongside `central_partial_lh`) — leaves residual state that a
*later* candidate's computation silently reads. This was **not**
isolated before time/risk tradeoff called for reverting; whoever
attempts this next should instrument at the level of "does candidate
i's score change if it's evaluated as candidate 0 in its own fresh step
vs. candidate N after N prior local evaluations in the same step" to
localize it precisely, rather than re-deriving the buffer-swap mechanism
from scratch (that part seems fine).

## What was reverted (for reference, not present in the tree anymore)

Two files were touched; both are now byte-identical to commit `5db2d1fb`.

### `tree/phylonode.h`

Added three small public methods to `PhyloNeighbor` (all removed):

```cpp
UBYTE* getScaleNum() { return scale_num; }
double getLhScaleFactor() { return lh_scale_factor; }

void swapPartialLhState(double *&lh, UBYTE *&sn, double &sf, int &computed) {
    double *tmpLh = partial_lh; partial_lh = lh; lh = tmpLh;
    UBYTE *tmpSn = scale_num; scale_num = sn; sn = tmpSn;
    double tmpSf = lh_scale_factor; lh_scale_factor = sf; sf = tmpSf;
    int tmpComputed = partial_lh_computed; partial_lh_computed = computed; computed = tmpComputed;
}
```

`swapPartialLhState` is the one genuinely reusable piece here: a clean,
minimal, symmetric primitive letting non-friend code (this standalone
tool is not a friend of `PhyloNeighbor`) install/restore a dedicated
buffer on a direction without it ever going through `nullptr`. If this
optimization is retried, this method (or something like it) is very
likely still the right starting point — it was never implicated in the
bug.

### `tree/spr_topology_test.cpp`

New free functions added right after `resetLikelihoodBuffers` (all
removed — reconstruct from this description, not copied verbatim since
the exact surrounding line numbers have since shifted):

- `struct SPRLocalLhCache` — 6 persistent `(double* lh, UBYTE* scaleNum)`
  scratch buffer pairs, `lhCount`/`scaleCount` sizes from
  `tree.getPartialLhSize()`/`getScaleNumSize()`, allocated once via
  `allocateSPRLocalLhCache(tree, cache)` (uses the free template
  `aligned_alloc<T>`/`aligned_free<T>` already in `phylotree.h`), freed
  via `freeSPRLocalLhCache`.
- `struct SPRLocalInvalidationState` — per-call record of the 6
  `PhyloNeighbor*` touched and each one's saved `(partial_lh, scale_num,
  lh_scale_factor, partial_lh_computed)`.
- `beginLocalSPRInvalidation(cache, state, dad1, dad2, node2, sibling1,
  sibling2)` — for each of the 6 directions
  (`dad1<->dad2`, `dad1<->node2`, `sibling1<->sibling2`), passes a
  *copy* of `cache.lh[i]`/`cache.scaleNum[i]` into
  `swapPartialLhState(...)` (critically: not `cache.lh[i]` itself by
  reference, so the cache's own slots are never mutated by the swap —
  only the local copy is, and that local copy is what gets saved into
  `state.savedLh[i]` etc.). Marks each dirty (`computed=0`, `sf=0.0`
  passed in).
- `computeLocalSPRLikelihood(tree, dad1, dad2, node2)` — the two-call
  `computeLikelihoodBranch` pattern described above, returns the second
  call's score.
- `endLocalSPRInvalidationDiscard(state)` — swaps each direction back to
  its saved original state; a provably perfect round trip (this part
  was never observed to be wrong in isolation).
- `endLocalSPRInvalidationKeep(tree, cache, state)` — asserts all 6 are
  computed; if any `state.savedLh[i]` was null (original buffer had been
  stolen by unrelated shared-pool churn and never reclaimed), falls back
  to discard + full `resetLikelihoodBuffers` for that one occurrence;
  otherwise `memcpy`s each scratch buffer's freshly-computed content
  into the *original* buffer, marks it valid, and swaps the original
  pointer back onto the `PhyloNeighbor` — this is the step that (per the
  unresolved bug above) is suspected to interact badly with prior
  reuse of the same 6 slots within a step, though it was never proven
  to be the actual fault vs. `computeLocalSPRLikelihood` itself.

`scoreTrialSPRMove` gained an optional trailing `SPRLocalLhCache
*localCache = nullptr` parameter; when non-null and
`!reoptimizeBranchLengths`, used the local path instead of
`resetLikelihoodBuffers`; every other call site (branch-length-compare,
`computeSiblingCompatibilityScore`/"sweep" ranking, the "sweep"
post-processing block itself) was left passing the implicit `nullptr`
default, deliberately never touched.

`runHillClimb` allocated one `SPRLocalLhCache sprLocalCache;` right after
its own `tree.initializeAllPartialLh()`, freed it right before its single
`delete aln; return 2;` exit point. The two exhaustive/fast-mode
candidate-scoring call sites (inside the main step loop only — *not* the
"sweep" block) passed `&sprLocalCache`. The commit/reject bookkeeping
(`applySPRTracked` → decide → keep-or-`rollbackSPRTracked`) had its
`resetLikelihoodBuffers()` calls replaced with
`beginLocalSPRInvalidation`/`computeLocalSPRLikelihood`/
`endLocalSPRInvalidationKeep` (on accept) or `endLocalSPRInvalidationDiscard`
(on reject), gated the same way (`!reoptimizeBranchLengths`).

## If retrying

1. Keep `swapPartialLhState` — it's sound and is the actual fix for the
   two *previous* attempts' crash class.
2. Keep the two-call `computeLikelihoodBranch` pattern for full 6/6
   coverage — also verified sound.
3. Do **not** trust `endLocalSPRInvalidationKeep`/the exhaustive-mode
   multi-candidate reuse of one `SPRLocalLhCache` without first writing
   an isolated reproduction: apply a move, do N dummy
   begin/compute/discard cycles on *unrelated* directions elsewhere in
   the tree (simulating N prior candidates in the same step), *then*
   evaluate the real move and compare against a from-scratch reference.
   That isolates whether it's genuinely about *reuse count* rather than
   about "commit" specifically.
4. Verification harness worth keeping (recreate, don't have to guess):
   after each step, `resetLikelihoodBuffers(tree);
   tree.computeLikelihood();` and compare against the incrementally
   tracked `curScore` — cheap, catches drift immediately, was
   instrumental in finding this bug in under an hour once added (versus
   the much longer, inconclusive reasoning-only phase before it).
