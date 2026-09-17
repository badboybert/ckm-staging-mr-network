# -*- coding: utf-8 -*-
"""Step 77: EXACT staging-order permutation P for the two East-Asian meta networks.

WHY THIS EXISTS (round-7 item C057; PI decision Q5, 2026-09-16): the Results quoted the BBJ+TPMI
meta staging concordances as "0.909, P = 0.059" and "1.000, P = 0.131", but those P-values came from
27_eas_meta.py's MONTE-CARLO label permutation (n = 10,000) while the Methods say the null is
"enumerated exhaustively". This step re-scores the same two meta networks with the SAME exhaustive
enumerator the European, BBJ and TPMI staging tests use (07_staging.py / 59_staging_gwstrict.py), so
the Methods sentence is true of every staging P in the paper.

Edge inclusion is unchanged from 27_eas_meta.py: IVW P < 0.05/34 (the meta network's own Bonferroni
bar over the 34 edges that could be meta-analysed) and Steiger-direction-correct. The enumerated null
permutes AHA stage labels over the nodes that carry at least one retained edge - the same convention
as every other staging test here, and now stated in Supplementary Methods 4.

Writes: results/staging_eas_meta_exact.txt  and  results/staging_eas_meta_exact.csv
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

import csv, io, os
from collections import Counter

BASE = P4_BASE
RES = os.path.join(BASE, "results")
STAGE = {"BMI": 1, "SBP": 2, "TG": 2, "HDL": 2, "TC": 2, "LDL": 2, "HbA1c": 2, "T2D": 2, "eGFR": 2,
         "CKD": 2, "CAD": 4, "HF": 4, "Stroke": 4, "AF": 4}
VARIANTS = [("fixed-effect", "network_eas_meta_fixed.csv"),
            ("random-effect", "network_eas_meta_random.csv")]


def load(p):
    with io.open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def retained(rows, bonf):
    """Significant, Steiger-correct, staged edges - the selection 27_eas_meta.py's staging used."""
    out = []
    for r in rows:
        try:
            p = float(r["ivw_p"])
        except (TypeError, ValueError, KeyError):
            continue
        if p < bonf and str(r.get("steiger", r.get("steiger_correct", ""))).upper() == "TRUE" \
           and r["exposure"] in STAGE and r["outcome"] in STAGE:
            out.append((r["exposure"], r["outcome"]))
    return out


def concordance(edges, smap=None):
    s = smap or STAGE
    cross = [(a, b) for a, b in edges if s[a] != s[b]]
    fwd = [(a, b) for a, b in cross if s[a] < s[b]]
    back = [(a, b) for a, b in cross if s[a] > s[b]]
    return (len(fwd) / len(cross) if cross else float("nan")), fwd, back, cross


def _msperm(counts):
    n = sum(counts.values()); cur = [None] * n

    def rec(i):
        if i == n:
            yield tuple(cur); return
        for v in sorted(counts):
            if counts[v]:
                counts[v] -= 1; cur[i] = v
                yield from rec(i + 1)
                counts[v] += 1
    return rec(0)


def exact_perm_p(edges):
    """Enumerate every distinct assignment of the observed multiset of stage labels to the nodes
    that carry an edge (identical to 59_staging_gwstrict.exact_perm_p)."""
    nodes = sorted({x for e in edges for x in e})
    ob, _, _, _ = concordance(edges)
    ge = tot = 0
    for pm in _msperm(dict(Counter(STAGE[n] for n in nodes))):
        tot += 1
        smap = dict(zip(nodes, pm))
        cc = [(a, b) for a, b in edges if smap[a] != smap[b]]
        if cc and sum(1 for a, b in cc if smap[a] < smap[b]) / len(cc) >= ob - 1e-12:
            ge += 1
    return ge / tot, ge, tot, len(nodes)


def main():
    L = []
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("# Step 77 - EXACT staging permutation P for the BBJ+TPMI meta networks")
    out("# Enumerated null (no Monte Carlo); supersedes the n=10,000 empirical P in")
    out("# results/staging_eas_meta.txt for the two numbers quoted in the Results.")
    out("")
    recs = []
    for name, fn in VARIANTS:
        rows = load(os.path.join(RES, fn))
        bonf = 0.05 / len(rows)
        edges = retained(rows, bonf)
        conc, fwd, back, cross = concordance(edges)
        p, ge, tot, nnodes = exact_perm_p(edges)
        out(f"===== {name.upper()} META =====")
        out(f"  edges meta-analysed        : {len(rows)}   (Bonferroni bar 0.05/{len(rows)} = {bonf:.2e})")
        out(f"  significant Steiger-correct: {len(edges)} over {nnodes} nodes")
        out(f"  cross-stage edges          : {len(cross)}  forward {len(fwd)}  backward {len(back)}")
        out(f"  backward edges             : {back if back else 'none'}")
        out(f"  stage-order CONCORDANCE    : {conc:.3f}")
        out(f"  EXACT permutation P        : {ge}/{tot} = {p:.4f}   ({nnodes} nodes carrying an edge)")
        out(f"  => {'concordant beyond chance' if p < 0.05 else 'NOT beyond chance'}")
        out("")
        recs.append(dict(variant=name, edges_meta=len(rows), bonferroni=f"{bonf:.3e}",
                         retained_edges=len(edges), nodes=nnodes, cross=len(cross), forward=len(fwd),
                         backward=len(back), concordance=round(conc, 4), perm_ge=ge,
                         perm_total=tot, exact_p=round(p, 6)))
    out("NOTE: tests necessary POPULATION-LEVEL causal-ordering conditions, not within-person progression.")
    io.open(os.path.join(RES, "staging_eas_meta_exact.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    with io.open(os.path.join(RES, "staging_eas_meta_exact.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recs[0].keys())); w.writeheader(); w.writerows(recs)

    # Self-check: the enumerated P must agree with the Monte-Carlo run to within its own noise
    # (n = 10,000 -> SE ~ 0.003 at p ~ 0.06), otherwise the two tests are not scoring the same graph.
    mc = {"fixed-effect": 0.0592, "random-effect": 0.1314}
    print("\n--- SELF-CHECKS ---", flush=True)
    ok = True
    for r in recs:
        d = abs(r["exact_p"] - mc[r["variant"]])
        cond = d < 0.02
        print(f"  [{'PASS' if cond else 'FAIL'}] {r['variant']}: exact {r['exact_p']:.4f} vs "
              f"Monte-Carlo {mc[r['variant']]:.4f} (|d| = {d:.4f} < 0.02)", flush=True)
        ok &= cond
    for r, want in zip(recs, [0.909, 1.000]):
        cond = abs(r["concordance"] - want) < 5e-4
        print(f"  [{'PASS' if cond else 'FAIL'}] {r['variant']}: concordance {r['concordance']:.3f} "
              f"reproduces the published {want:.3f}", flush=True)
        ok &= cond
    print(f"\nALL TARGETS REPRODUCED: {ok}", flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
