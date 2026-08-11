# -*- coding: utf-8 -*-
"""Step 61b: CORRECTS the reading of step 61.

Step 61 swept prevalence over clinically plausible grids for every binary-involving staging edge and
reported "0 direction flips in 330 combinations". An adversarial check showed that result is close to
VACUOUS: in TwoSampleMR::get_r_from_lor, prevalence enters only through get_population_allele_
frequency (a small nudge to the allele frequency); it does not enter the residual-variance term. I
measured it directly:

    sum(r^2) at prevalence 0.40 vs 0.01  ->  ratio 1.0336   (a 40-fold prevalence range moves r^2 by 3.4%)

Across the grids step 61 actually used (e.g. CAD 0.03-0.10) the movement is ~0.3%. So a prevalence-
driven direction flip is only possible where the two sides' r^2 are within a fraction of a percent of
each other. "Zero flips" therefore says almost nothing about robustness on its own.

The informative quantity is the MARGIN: how far each edge sits from the flip boundary
r2_exposure == r2_outcome. This script computes it from step 61's own output, and states the
prevalence perturbation that would be required to flip each edge. It also records the second, separate
limitation the same check raised: step 61 reads fig1_edges.csv, i.e. edges ALREADY admitted to the
staging graph (which 07_staging.py filtered on steiger_correct == TRUE), so it is one-sided - it can
only find false positives, never edges the liability scale would ADMIT that the quantitative gate
rejected.

Out: results/steiger_margin.{txt,csv}
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

import csv, io, os, collections

BASE = P4_BASE
RES = os.path.join(BASE, "results")

# measured directly from TwoSampleMR::get_r_from_lor (see docstring)
R2_MOVE_PER_40X_PREVALENCE = 0.0336


def main():
    rows = list(csv.DictReader(io.open(os.path.join(RES, "steiger_prevalence_sweep.csv"),
                                       encoding="utf-8")))
    by_edge = collections.defaultdict(list)
    for r in rows:
        by_edge[(r["exposure"], r["outcome"])].append(r)

    out_rows, L = [], []
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("=== Steiger direction MARGIN (corrects the reading of step 61) ===")
    out("")
    out("Step 61 reported 0 direction flips across 330 edge x prevalence combinations. That test is")
    out("near-vacuous on its own: prevalence enters get_r_from_lor only through the allele-frequency")
    out(f"conversion, and a 40-FOLD prevalence range moves sum(r^2) by only {100*R2_MOVE_PER_40X_PREVALENCE:.1f}%")
    out("(~0.3% across the grids actually swept). What matters is the MARGIN below.")
    out("")
    out("margin = r2_exposure / r2_outcome at base prevalence. Steiger calls the direction correct when")
    out("this exceeds 1. A margin of 1.05 is fragile; a margin of 10 is not.")
    out("")
    out(f"{'edge':16s} {'cross':6s} {'r2_exp':>11s} {'r2_out':>11s} {'margin':>9s}  {'grid range of margin':>24s}")
    for (e, o), rs in sorted(by_edge.items()):
        base = next((r for r in rs if r["at_base"] == "TRUE"), rs[0])
        rx, ry = float(base["r2_exp"]), float(base["r2_out"])
        margin = rx / ry if ry else float("inf")
        margins = [float(r["r2_exp"]) / float(r["r2_out"]) for r in rs if float(r["r2_out"])]
        lo, hi = min(margins), max(margins)
        out(f"{e+'->'+o:16s} {base['cross']:6s} {rx:11.6f} {ry:11.6f} {margin:9.2f}  "
            f"{f'[{lo:.2f}, {hi:.2f}]':>24s}")
        out_rows.append(dict(exposure=e, outcome=o, cross=base["cross"], direction=base["direction"],
                             r2_exp=rx, r2_out=ry, margin=round(margin, 4),
                             margin_lo=round(lo, 4), margin_hi=round(hi, 4),
                             fragile=margin < 1.5))
    frag = [r for r in out_rows if r["fragile"]]
    mn = min(out_rows, key=lambda r: r["margin"])
    out("")
    out(f"Edges with margin < 1.5 (would be genuinely fragile): {len(frag)}"
        + (f" - {[r['exposure']+'->'+r['outcome'] for r in frag]}" if frag else ""))
    out(f"Smallest margin in the whole staging graph: {mn['exposure']}->{mn['outcome']} at "
        f"{mn['margin']:.2f}x.")
    out("")
    out("READ (honest):")
    out("- The defensible claim is NOT 'we swept prevalence and nothing flipped'. It is that the")
    out("  liability-scale direction calls sit far from the flip boundary: the smallest margin in the")
    out(f"  staging graph is {mn['margin']:.2f}x, whereas the entire plausible prevalence range perturbs")
    out(f"  r^2 by ~{100*R2_MOVE_PER_40X_PREVALENCE:.0f}% at most. Prevalence mis-specification cannot")
    out("  reach that boundary.")
    out("- REMAINING LIMITATION, not addressed here: this covers only edges ALREADY in the staging")
    out("  graph, which were filtered on the QUANTITATIVE Steiger gate. It cannot detect an edge the")
    out("  liability scale would ADMIT that the quantitative gate wrongly rejected. Testing that")
    out("  requires running the liability computation on the 19 Steiger-FALSE edges as well.")

    io.open(os.path.join(RES, "steiger_margin.txt"), "w", encoding="utf-8").write("\n".join(L))
    with io.open(os.path.join(RES, "steiger_margin.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)
    print("\nwrote results/steiger_margin.{txt,csv}")


if __name__ == "__main__":
    main()
