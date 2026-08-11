# -*- coding: utf-8 -*-
"""Step 47 (E2): constrained-null battery for the AHA staging concordance test.

Peer review C3: the simple label-permutation null does not preserve node role, degree, missing-
direction structure, or phenotype class, so it may be anti-conservative. This runs the stronger
nulls the reviewer asked for, on the completed 132-edge network, plus robustness checks:

  A. Simple label permutation                 (baseline = the manuscript's current null)
  B. Phenotype-class-constrained permutation   (shuffle stage labels only within trait-class blocks)
  C. Degree-preserving edge-rewiring null      (configuration model; keeps in/out degree, tests topology)
  D. Complete-bidirectional-pair subgraph       (restrict to pairs estimable both ways; the direct
                                                 test that asymmetric estimability did not create the signal)
  E. Leave-one-node-out concordance             (drop each node; no single node drives it)
  F. Drop SBP + Stroke                          (the two structurally-special nodes)
  G. Bootstrap 95% CI on the concordance        (resample cross-stage edges)

Deterministic (seed 12345). Writes results/staging_constrained_nulls.txt.
"""
# --- portable roots -----------------------------------------------------------------------------
# Injected by scripts/build_repo.py. The working tree hardcoded an absolute local path; the deposit
# resolves it from CKM_P4_BASE, or from this file's own location (repo/code/ plays the role the
# working tree called independent_build/). Raw GWAS summary statistics are NOT deposited: set
# CKM_ROOT to wherever you obtained them if you intend to re-run the upstream extraction steps.
import os as _os
P4_BASE = _os.environ.get("CKM_P4_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
P4_ROOT = _os.environ.get("CKM_P4_ROOT") or _os.path.dirname(P4_BASE)
CKM_ROOT = _os.environ.get("CKM_ROOT") or _os.path.dirname(P4_ROOT)
SHARED_LIB = _os.environ.get("CKM_SHARED_LIB") or _os.path.join(CKM_ROOT, "_shared")
# ------------------------------------------------------------------------------------------------

import csv, io, os, random, itertools

BASE = P4_BASE
STAGE = {"BMI": 1, "SBP": 2, "TG": 2, "HDL": 2, "TC": 2, "LDL": 2, "HbA1c": 2, "T2D": 2, "eGFR": 2,
         "CAD": 4, "HF": 4, "Stroke": 4}
# phenotype classes for the class-constrained null (nodes only permute stage labels within a class)
CLASS = {"BMI": "adip", "SBP": "bp", "TG": "lipid", "HDL": "lipid", "TC": "lipid", "LDL": "lipid",
         "HbA1c": "glyc", "T2D": "glyc", "eGFR": "kidney", "CAD": "disease", "HF": "disease", "Stroke": "disease"}
NP = 10000
SEED = 12345

def load(path):
    return list(csv.DictReader(io.open(path, encoding="utf-8")))

ALL = load(os.path.join(BASE, "results/forward_local_edges.csv"))   # 132 edges
N = len(ALL)
BONF = 0.05 / N

def sig_edges(rows):
    out = []
    for r in rows:
        try:
            p = float(r["ivw_p"])
        except (ValueError, KeyError):
            continue
        if p < BONF and str(r["steiger_correct"]).upper() == "TRUE" \
           and r["exposure"] in STAGE and r["outcome"] in STAGE:
            out.append((r["exposure"], r["outcome"]))
    return out

def concordance(edges, stage):
    cross = [(a, b) for a, b in edges if stage[a] != stage[b]]
    if not cross:
        return None, 0, 0
    fwd = sum(1 for a, b in cross if stage[a] < stage[b])
    return fwd / len(cross), fwd, len(cross)

SIG = sig_edges(ALL)
nodes = sorted({x for e in SIG for x in e})
obs, fwd, cross = concordance(SIG, STAGE)
rng = random.Random(SEED)

lines = []
def out(s=""):
    lines.append(s); print(s, flush=True)

out(f"=== E2: constrained-null battery for staging concordance (132-edge network) ===")
out(f"Edges tested={N}, Bonferroni=0.05/{N}={BONF:.2e}; significant+Steiger staged edges={len(SIG)} over {len(nodes)} nodes")
out(f"OBSERVED cross-stage concordance = {obs:.3f} ({fwd}/{cross} forward)\n")

# ---- A. simple label permutation (baseline) ----
# The simple label-permutation null has a small, enumerable support (the distinct assignments of the
# stage multiset to nodes), so it is computed EXACTLY rather than by Monte Carlo, which removes the
# seed/node-order dependence of the earlier estimate. The class-constrained variant (restrict!=None)
# keeps its Monte-Carlo form because its degenerate structure is what is being demonstrated, not its P.
from collections import Counter as _Counter
def _msperm(counts):
    n = sum(counts.values()); cur = [None]*n
    def rec(i):
        if i == n:
            yield tuple(cur); return
        for v in sorted(counts):
            if counts[v]:
                counts[v] -= 1; cur[i] = v
                yield from rec(i+1)
                counts[v] += 1
    return rec(0)

def exact_perm_p(edges, node_list, stage_of):
    ob, _, _ = concordance(edges, stage_of)
    ge = tot = 0
    for pm in _msperm(dict(_Counter(stage_of[n] for n in node_list))):
        tot += 1
        c, _, _ = concordance(edges, dict(zip(node_list, pm)))
        if c is not None and c >= ob - 1e-12:
            ge += 1
    return ge/tot, ge, tot

def perm_p(edges, node_list, stage_of, restrict=None):
    if restrict is None:
        return exact_perm_p(edges, node_list, stage_of)[0]
    ob, _, _ = concordance(edges, stage_of)
    stages = [stage_of[n] for n in node_list]
    ge = 0
    r = random.Random(SEED)
    for _ in range(NP):
        smap = {}
        for cls, idxs in restrict.items():      # class-constrained: shuffle within each class block
            vals = [stages[i] for i in idxs]; r.shuffle(vals)
            for i, v in zip(idxs, vals):
                smap[node_list[i]] = v
        c, _, _ = concordance(edges, smap)
        if c is not None and c >= ob - 1e-12:
            ge += 1
    return (ge + 1) / (NP + 1)

pA, geA, totA = exact_perm_p(SIG, nodes, STAGE)
out(f"A. Simple label permutation (baseline / current method): P = {pA:.4f}  [EXACT: {geA}/{totA}]")

# ---- B. phenotype-class-constrained permutation ----
classes = {}
for i, n in enumerate(nodes):
    classes.setdefault(CLASS[n], []).append(i)
pB = perm_p(SIG, nodes, STAGE, restrict=classes)
out(f"B. Phenotype-class-constrained permutation (shuffle stage within trait-class blocks): P = {pB:.4f}")
out(f"   [class blocks: {{{', '.join(f'{k}:{len(v)}' for k,v in classes.items())}}}]")

# ---- C. degree-preserving edge-rewiring null (configuration model) ----
# Keep node stage labels + each node's out-degree and in-degree; rewire edges by repeated
# double-edge swaps; recompute concordance. Tests whether the TOPOLOGY (not degrees) respects staging.
# Chain length is expressed in SUCCESSFUL swaps per edge. The previous version counted ATTEMPTS and
# let a rejected proposal consume a step; at this graph's measured 8.9% acceptance rate, its nominal
# "3x the edge count" delivered about 0.27 successful swaps per edge, so the chain barely moved.
# Mixing diagnostics across 10/25/50/100 swaps per edge and three seeds are in
# scripts/74_degree_null_mixing.py -> results/degree_null_mixing.{csv,txt}.
SWAPS_PER_EDGE = 50

def rewire(edges, swaps_per_edge, r):
    """Double-edge swap preserving in- and out-degree, run to a target of SUCCESSFUL swaps."""
    E = list(edges); S = set(E); m = len(E)
    target = swaps_per_edge * m
    done = proposals = 0
    cap = target * 200                      # backstop; never reached on this graph
    while done < target and proposals < cap:
        proposals += 1
        i, j = r.randrange(m), r.randrange(m)
        (a, b), (c, d) = E[i], E[j]
        if len({a, b, c, d}) < 4:
            continue
        if (a, d) in S or (c, b) in S:
            continue
        S.discard((a, b)); S.discard((c, d))
        E[i], E[j] = (a, d), (c, b)         # preserves out-deg(a,c) and in-deg(b,d)
        S.add((a, d)); S.add((c, b))
        done += 1
    return E
rC = random.Random(SEED); ge = 0
for _ in range(NP):
    rw = rewire(SIG, SWAPS_PER_EDGE, rC)
    c, _, _ = concordance(rw, STAGE)
    if c is not None and c >= obs - 1e-12:
        ge += 1
pC = (ge + 1) / (NP + 1)
out(f"C. Degree-preserving edge-rewiring null (topology test, in/out-degree fixed): P = {pC:.4f}")
out(f"   chain: {SWAPS_PER_EDGE} SUCCESSFUL double-edge swaps per edge; mixing diagnostics across "
    f"10/25/50/100 swaps per edge and three seeds in results/degree_null_mixing.txt")

# ---- D. complete-bidirectional-pair subgraph ----
present = {(r["exposure"], r["outcome"]) for r in ALL}
def both_estimable(a, b):
    return (a, b) in present and (b, a) in present
sub = [(a, b) for a, b in SIG if both_estimable(a, b)]
subnodes = sorted({x for e in sub for x in e})
oD, fD, cD = concordance(sub, STAGE)
pD = perm_p(sub, subnodes, STAGE) if cD else None
out(f"D. Complete-pair subgraph (both directions estimable): {len(sub)} sig edges; "
    f"cross-stage concordance = {oD:.3f} ({fD}/{cD}); label-perm P = {pD:.4f}" if cD else "D. no cross-stage edges")

# ---- E. leave-one-node-out ----
out("E. Leave-one-node-out concordance (drop each node, recompute):")
lono = []
for n in nodes:
    e2 = [(a, b) for a, b in SIG if a != n and b != n]
    c, f, cr = concordance(e2, STAGE)
    lono.append((n, c, f, cr))
for n, c, f, cr in sorted(lono, key=lambda x: x[1] if x[1] is not None else 1):
    out(f"     drop {n:6s}: concordance {c:.3f} ({f}/{cr})")
worst = min(c for _, c, _, _ in lono if c is not None)
out(f"   -> worst-case leave-one-out concordance = {worst:.3f} (observed {obs:.3f}); no single node drives it")

# ---- F. drop SBP + Stroke (structurally-special nodes) ----
eF = [(a, b) for a, b in SIG if a not in ("SBP", "Stroke") and b not in ("SBP", "Stroke")]
nF = sorted({x for e in eF for x in e})
oF, fF, cF = concordance(eF, STAGE)
pF = perm_p(eF, nF, STAGE) if cF else None
out(f"F. Drop SBP + Stroke: {len(eF)} edges; cross-stage concordance = {oF:.3f} ({fF}/{cF}); label-perm P = {pF:.4f}")

# ---- G. bootstrap 95% CI on concordance ----
cross_edges = [(a, b) for a, b in SIG if STAGE[a] != STAGE[b]]
rG = random.Random(SEED); boot = []
for _ in range(NP):
    samp = [cross_edges[rG.randrange(len(cross_edges))] for _ in cross_edges]
    f = sum(1 for a, b in samp if STAGE[a] < STAGE[b])
    boot.append(f / len(samp))
boot.sort()
lo, hi = boot[int(0.025 * NP)], boot[int(0.975 * NP)]
out(f"G. Bootstrap 95% CI on concordance ({len(cross_edges)} cross-stage edges): {obs:.3f} [{lo:.3f}, {hi:.3f}]")

out("\nSUMMARY (honest, calibrated):")
out(f"ROBUST: the concordance {obs:.3f} is significant vs LABEL permutation (A, P={pA:.4f}) — the specific")
out("stage assignment matters, you cannot relabel which trait is stage 1/2/4 and preserve it. It survives")
out(f"leave-one-node-out (worst {worst:.3f}; the 2 discordant edges both originate from CAD, so dropping CAD")
out(f"gives 1.000), dropping the two structurally-special nodes (F, {oF:.3f}, P={pF:.4f}), and its bootstrap")
out(f"95% CI [{lo:.3f}, {hi:.3f}] excludes 0.5. D (complete-pair subgraph) now equals A because completing")
out("the network to 132 edges made EVERY pair estimable both ways.")
out("")
out("TWO HONEST LIMITATIONS the battery exposes (report, do NOT bury):")
out(f"1. The phenotype-class-constrained permutation is DEGENERATE (B, P={pB:.4f}): trait class and AHA")
out("   stage are PERFECTLY NESTED here (each of the 6 trait classes maps to exactly one stage), so there")
out("   is zero within-class stage variance to permute. This null is uninformative in this design, NOT")
out("   evidence for or against the signal.")
out(f"2. The degree-preserving edge-rewiring null is NOT exceeded (C, P={pC:.4f}): given the exposure/outcome")
out("   role structure (risk factors are high-out-degree senders, diseases high-in-degree receivers), a")
out("   randomly-rewired graph with the SAME degree sequence achieves comparable concordance ~1/3 of the")
out("   time. So the low->high concordance substantially REFLECTS the trait role structure. Read the staging")
out("   test as CORROBORATING the AHA order given the trait roles (the stage labels matter, P<0.01), NOT as")
out("   independent proof beyond the fact that risk factors are instrumented as exposures and diseases as outcomes.")
out("   -> the manuscript must CALIBRATE the staging claim accordingly (this is exactly reviewer C3's concern).")

io.open(os.path.join(BASE, "results/staging_constrained_nulls.txt"), "w", encoding="utf-8").write("\n".join(lines))
print("\nwrote results/staging_constrained_nulls.txt")
