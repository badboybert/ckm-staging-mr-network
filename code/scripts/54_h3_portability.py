# -*- coding: utf-8 -*-
"""Step 54 (H3): network-aware cross-ancestry portability with dependence-robust resampling.

Reviewer M4: the node net-flow / Kendall-tau portability statistic is power-dependent and the
binomial sign-concordance test treats correlated edges as independent. This recomputes portability
on the COMMON ESTIMABLE SUBGRAPH (edges estimated in both EUR and EAS-BBJ) with:
  - sign concordance + effect correlation (on comparable-scale edges only, per the scale dictionary),
  - EDGE-block bootstrap (resample edges) and NODE-block bootstrap (resample nodes with all their
    edges) -> CIs that respect edge dependence,
  - a sign-flip permutation null.
EAS = Biobank Japan (network_eas_edges.csv), the cohort the paper's cross-ancestry comparison uses.
Writes results/h3_portability.txt.
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

import csv, io, os, random, math

BASE = P4_BASE
# non-comparable-scale traits (per scale_dictionary.md): exclude edges touching these from the
# magnitude CORRELATION (sign is scale-invariant, so keep them for sign concordance). SBP is included
# here for the CORRELATION only, because EUR SBP is per-mmHg while EAS SBP is per-SD (needs ×19.3);
# leaving it in drags the correlation toward zero. Sign concordance keeps SBP (sign is scale-free).
NONCOMP = {"HbA1c", "eGFR", "SBP"}

def load(p):
    return {(r["exposure"], r["outcome"]): r for r in csv.DictReader(io.open(p, encoding="utf-8"))}

EUR = load(os.path.join(BASE, "results/forward_local_edges.csv"))
EAS = load(os.path.join(BASE, "results/network_eas_edges.csv"))
BONF = 0.05 / len(EUR)

common = sorted(set(EUR) & set(EAS))
rows = []
for k in common:
    be, pe = float(EUR[k]["ivw_b"]), float(EUR[k]["ivw_p"])
    ba = float(EAS[k]["ivw_b"])
    rows.append(dict(edge=k, be=be, ba=ba, pe=pe,
                     comparable=(k[0] not in NONCOMP and k[1] not in NONCOMP)))

lines = []
def out(s=""):
    lines.append(s); print(s, flush=True)

out("=== H3: network-aware cross-ancestry portability (EUR vs EAS-BBJ), dependence-robust ===")
out(f"Common estimable edges (both ancestries): {len(rows)}\n")

def sign_conc(rs):
    ss = [r for r in rs if r["be"] != 0 and r["ba"] != 0]
    return sum(1 for r in ss if (r["be"] > 0) == (r["ba"] > 0)) / len(ss), len(ss)

def pearson(rs):
    xs = [r["be"] for r in rs if r["comparable"]]
    ys = [r["ba"] for r in rs if r["comparable"]]
    n = len(xs)
    if n < 3:
        return None, n
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs); syy = sum((y - my) ** 2 for y in ys)
    return (sxy / math.sqrt(sxx * syy) if sxx and syy else None), n

# --- observed, over ALL common edges and over the Bonferroni-significant-in-EUR subset ---
for label, rs in [("all common edges", rows),
                  ("EUR-Bonferroni-significant subset", [r for r in rows if r["pe"] < BONF])]:
    sc, nsc = sign_conc(rs); r_, nr = pearson(rs)
    out(f"[{label}] n={len(rs)}")
    out(f"   sign concordance = {sc:.3f} ({int(round(sc*nsc))}/{nsc})")
    out(f"   effect correlation (comparable-scale, n={nr}) = " + (f"{r_:.3f}" if r_ is not None else "n/a"))

RS = [r for r in rows if r["pe"] < BONF]      # the load-bearing subset the paper claims ports
NP = 10000
rng = random.Random(12345)

# --- edge-block bootstrap ---
def boot(sampler):
    scs, corrs = [], []
    for _ in range(NP):
        samp = sampler()
        sc, n = sign_conc(samp)
        if n:
            scs.append(sc)
        c, nn = pearson(samp)
        if c is not None:
            corrs.append(c)
    scs.sort(); corrs.sort()
    def ci(v):
        return (v[int(0.025*len(v))], v[int(0.975*len(v))]) if v else (float("nan"),)*2
    return ci(scs), ci(corrs)

def edge_sampler():
    return [RS[rng.randrange(len(RS))] for _ in RS]

# --- node-block bootstrap: resample nodes, keep edges whose BOTH endpoints are sampled ---
nodes = sorted({x for r in RS for x in r["edge"]})
def node_sampler():
    keep = {n for n in nodes if rng.random() < 0.632}    # ~bootstrap inclusion prob
    return [r for r in RS if r["edge"][0] in keep and r["edge"][1] in keep] or RS

sc_ci_e, corr_ci_e = boot(edge_sampler)
sc_ci_n, corr_ci_n = boot(node_sampler)
out("")
out(f"EUR-significant subset ({len(RS)} edges), dependence-robust 95% CIs:")
out(f"   sign concordance : edge-block [{sc_ci_e[0]:.3f}, {sc_ci_e[1]:.3f}]  node-block [{sc_ci_n[0]:.3f}, {sc_ci_n[1]:.3f}]")
out(f"   effect correlation: edge-block [{corr_ci_e[0]:.3f}, {corr_ci_e[1]:.3f}]  node-block [{corr_ci_n[0]:.3f}, {corr_ci_n[1]:.3f}]")

# --- sign-flip permutation null: is observed sign concordance above chance? ---
obs_sc, nsc = sign_conc(RS)
ge = 0
for _ in range(NP):
    flips = sum(1 for r in RS if r["be"] != 0 and r["ba"] != 0 and (rng.random() < 0.5) == (r["be"] > 0))
    if flips / nsc >= obs_sc - 1e-12:
        ge += 1
out(f"   sign-flip permutation P (observed {obs_sc:.3f} vs random-sign) = {(ge+1)/(NP+1):.4f}")

out("")
out("READ: report portability on the COMMON ESTIMABLE SUBGRAPH with dependence-robust CIs (edge- and")
out("node-block bootstrap), NOT the binomial that assumes edge independence. Sign concordance is the")
out("robust portability signal (scale-invariant); the magnitude correlation is reported only on")
out("comparable-scale edges. State 'suggestive topological/directional concordance' with these CIs, not")
out("'the topology ports/replicates' (reviewer M4).")

io.open(os.path.join(BASE, "results/h3_portability.txt"), "w", encoding="utf-8").write("\n".join(lines))
print("\nwrote results/h3_portability.txt")
