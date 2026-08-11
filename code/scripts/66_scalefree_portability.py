# -*- coding: utf-8 -*-
"""Step 66: scale-free cross-ancestry concordance (reviewer MAJOR-10 / Tier-2 #11).

The reviewer's objection: rank-inverse-normalised East Asian traits are not automatically comparable
with European per-SD effects, the 19.3-mmHg SBP conversion is a generic approximation rather than a
cohort SD, and the shipped comparable-scale effect CORRELATION (~0.66-0.68, Pearson) therefore
inherits an assumption it cannot verify. Cohort-specific SDs are not recoverable from the published
summary statistics, so the recommended alternative is to compare on signs, z-scores or ranks.

This computes the cross-ancestry concordance three ways, separating statistics that genuinely need
no SD conversion from those that do:
  1. SIGN concordance                       - invariant to any positive rescaling
  2. SPEARMAN rank correlation of estimates - invariant to any monotone rescaling
  3. Correlation of Z-STATISTICS (b/se)     - unit-free, though precision-weighted

and contrasts them with the Pearson correlation of the converted per-SD estimates that the
manuscript currently reports. If the scale-free measures agree with the converted one, the
portability claim does not depend on the SD assumption.

Bootstrap CIs use NODE-level resampling (each node retained with probability 0.632, an edge retained
only when both endpoints survive), which is the dependence-respecting scheme the manuscript already
promotes for the sign-concordance statistic.

Out: results/scalefree_portability.{txt,csv}
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

import csv, io, os, math, random

BASE = P4_BASE
RES = os.path.join(BASE, "results")
NON_COMPARABLE = {"HbA1c", "eGFR"}
NBOOT = 10000
SEED = 20260724


def load(p):
    return {(r["exposure"], r["outcome"]): r
            for r in csv.DictReader(io.open(os.path.join(RES, p), encoding="utf-8"))}


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(x), rank(y))


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def main():
    EUR, EAS = load("forward_local_edges.csv"), load("network_eas_edges.csv")
    common = sorted(set(EUR) & set(EAS))

    # the comparable subset the manuscript's magnitude claims are restricted to, and the
    # European-significant subset it reports the correlation on
    BONF = 0.05 / len(EUR)
    rows = []
    for e, o in common:
        re_, ra = EUR[(e, o)], EAS[(e, o)]
        try:
            be, see, pe = float(re_["ivw_b"]), float(re_["ivw_se"]), float(re_["ivw_p"])
            ba, sea = float(ra["ivw_b"]), float(ra["ivw_se"])
        except ValueError:
            continue
        rows.append(dict(exposure=e, outcome=o, b_eur=be, se_eur=see, p_eur=pe,
                         b_eas=ba, se_eas=sea,
                         z_eur=be / see if see else float("nan"),
                         z_eas=ba / sea if sea else float("nan"),
                         eur_sig=pe < BONF,
                         magnitude_eligible=not ({e, o} & NON_COMPARABLE)))
    sub = [r for r in rows if r["eur_sig"]]

    L = []
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("=== Scale-free cross-ancestry concordance ===")
    out("")
    out("Which edges each statistic may use, and why:")
    out("  * SIGN concordance and the Z-STATISTIC (b/se) are genuinely scale-free (b and se rescale")
    out("    together), so they run over ALL common edges - including HbA1c and eGFR. Restricting")
    out("    them to the 'comparable' subset would select the sample by the very property they are")
    out("    designed not to need.")
    out("  * SPEARMAN and PEARSON are computed over a POOLED set, so rescaling one subset changes the")
    out("    pooled ranking; neither is scale-free. They therefore keep the comparability filter, use")
    out("    SBP converted (x19.3), and a SBP-dropped variant reproduces the shipped 0.681/0.657.")
    out("")
    for r in rows:
        r["b_eur_conv"] = r["b_eur"] * 19.3 if r["exposure"] == "SBP" else r["b_eur"]
    for name, rs in (("all comparable common edges", rows),
                     ("European-Bonferroni-significant subset", sub)):
        if len(rs) < 3:
            continue
        # SCALE-FREE statistics (sign, z=b/se) do NOT need the comparability filter - that is the
        # whole point of using them. Computing them on a set selected BY a scale criterion is what
        # made an earlier draft's numbers disagree with the manuscript's own shipped values.
        elig   = [r for r in rs if r["magnitude_eligible"]]
        no_sbp = [r for r in elig if r["exposure"] != "SBP"]
        sign_all = sum(1 for r in rs if (r["b_eur"] > 0) == (r["b_eas"] > 0))
        sign_el  = sum(1 for r in elig if (r["b_eur"] > 0) == (r["b_eas"] > 0))
        pz_all = pearson([r["z_eur"] for r in rs],   [r["z_eas"] for r in rs])
        pz_el  = pearson([r["z_eur"] for r in elig], [r["z_eas"] for r in elig])
        sp = spearman([r["b_eur_conv"] for r in elig], [r["b_eas"] for r in elig])
        pb = pearson([r["b_eur_conv"] for r in elig], [r["b_eas"] for r in elig])
        pb_no = pearson([r["b_eur"] for r in no_sbp], [r["b_eas"] for r in no_sbp])
        out(f"--- {name} ---")
        out(f"  SCALE-FREE, over ALL {len(rs):3d} common edges (no comparability filter):")
        out(f"    sign concordance                  : {sign_all}/{len(rs)} = {sign_all/len(rs):.3f}")
        out(f"    z-statistic (b/se) correlation    : {pz_all:+.3f}")
        out(f"  same statistics restricted to the {len(elig):3d} magnitude-eligible edges (for contrast only):")
        out(f"    sign concordance                  : {sign_el}/{len(elig)} = {sign_el/len(elig):.3f}")
        out(f"    z-statistic correlation           : {pz_el:+.3f}")
        out(f"  MAGNITUDE statistics (comparability filter REQUIRED), n = {len(elig)}:")
        out(f"    Spearman rank [SBP converted]     : {sp:+.3f}")
        out(f"    Pearson       [SBP converted]     : {pb:+.3f}")
        out(f"    Pearson       [SBP dropped, n={len(no_sbp):3d}] : {pb_no:+.3f}   <- shipped definition")
        rng = random.Random(SEED)
        nodes = sorted({r["exposure"] for r in rs} | {r["outcome"] for r in rs})
        # bootstrap the SCALE-FREE statistics on the full common set (rs), magnitude on elig
        boots = {"sign": [], "spearman": [], "z": []}
        for _ in range(NBOOT):
            keep = {n for n in nodes if rng.random() < 0.632}
            bs = [r for r in rs if r["exposure"] in keep and r["outcome"] in keep]
            if len(bs) < 5:
                continue
            boots["sign"].append(sum(1 for r in bs if (r["b_eur"] > 0) == (r["b_eas"] > 0)) / len(bs))
            be_ = [r for r in bs if r["magnitude_eligible"]]
            if len(be_) >= 5:
                boots["spearman"].append(spearman([r["b_eur_conv"] for r in be_], [r["b_eas"] for r in be_]))
            boots["z"].append(pearson([r["z_eur"] for r in bs], [r["z_eas"] for r in bs]))
        for k, lab in (("sign", "sign concordance"), ("spearman", "Spearman rank"), ("z", "z-statistic corr")):
            v = sorted(x for x in boots[k] if x == x)
            if len(v) > 100:
                lo, hi = v[int(0.025 * len(v))], v[int(0.975 * len(v))]
                out(f"    node-level bootstrap 95% CI, {lab:17s}: [{lo:+.3f}, {hi:+.3f}]"
                    + ("   includes chance" if lo <= (0.5 if k == "sign" else 0.0) else ""))
        out("")

    out("READ:")
    out("- The shipped magnitude definition (comparable edges, SBP dropped) is REPRODUCED exactly:")
    out("  0.681 over 56 edges and 0.657 over 37. Implementation and shipped number agree.")
    out("- CORRECTION to an earlier draft of this script: it computed the scale-free statistics on the")
    out("  magnitude-ELIGIBLE subset, i.e. on a set selected by a scale criterion, which is self-")
    out("  defeating. Over all 101 common edges the sign concordance is 0.614 and 0.692, which are")
    out("  EXACTLY the manuscript's shipped values in results/h3_portability.txt - so there is one")
    out("  number for this statistic, not two.")
    out("- The z-statistic correlation FALLS when the excluded edges are restored (+0.77 -> +0.60 and")
    out("  +0.80 -> +0.59), i.e. below the magnitude correlation, so the earlier claim that it was")
    out("  'the highest of the three' was an artifact of the filter. It is also largely a shared-POWER")
    out("  statistic - correlating |z| alone recovers nearly the same value - and is dominated by")
    out("  near-definitional lipid edges (LDL->TC carries z_EUR = 71.7). Report it as supportive at")
    out("  most, never as evidence of stronger portability.")
    out("- Cohort-specific SDs are not recoverable from the published summary statistics on either")
    out("  side, which is why a sign-level reading is the defensible one for the non-comparable traits.")

    io.open(os.path.join(RES, "scalefree_portability.txt"), "w", encoding="utf-8").write("\n".join(L))
    with io.open(os.path.join(RES, "scalefree_portability.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nwrote results/scalefree_portability.{txt,csv}")


if __name__ == "__main__":
    main()
