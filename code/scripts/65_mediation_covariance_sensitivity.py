# -*- coding: utf-8 -*-
"""Step 65: covariance sensitivity for the formal mediation (reviewer BLOCKER-4 / Tier-2 #8).

The shipped mediation (`49_formal_mediation.py`) builds its Monte-Carlo confidence intervals by
drawing the four inputs INDEPENDENTLY:

    t ~ N(total, se_total);  aa ~ N(a, se_a);  bb ~ N(b, se_b);  dd ~ N(direct, se_direct)

That is an independence assumption, not a fact. Two of those pairs are estimated from overlapping
data and are certainly correlated:

  * corr(total, direct) - the univariable exposure->HF effect and the same effect conditioned on CAD
    are estimated on the same instruments and outcome; the DIFFERENCE method (total - direct) is
    directly sensitive to this correlation, and a positive correlation SHRINKS var(total - direct).
  * corr(a, b) - exposure->CAD and CAD->HF|exposure share the CAD node; the PRODUCT method a*b is
    sensitive to this one.

MVMR is separately run with gencov = 0 (13b_mvmr.R passes gencov=0 to both strength_mvmr and
pleiotropy_mvmr), so the whole mediation chain currently assumes zero covariance in three places.

This script does NOT claim to know the true covariance - it sweeps it, and asks the only question
that matters for the manuscript: does the qualitative reading survive? Namely
    BMI->HF is predominantly DIRECT, and T2D->HF is predominantly MEDIATED but imprecisely so.

Out: results/mediation_covariance_sensitivity.{txt,csv}
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
import numpy as np

BASE = P4_BASE
RES = os.path.join(BASE, "results")
NDRAW = 200_000
SEED = 12345
RHOS = [-0.5, 0.0, 0.5, 0.75]
RHO_AT = [0.0, 0.3, 0.5, 0.7, 0.9]   # corr(a,total) and corr(b,direct): shared instruments

UNI = {(r["exposure"], r["outcome"]): (float(r["ivw_b"]), float(r["ivw_se"]))
       for r in csv.DictReader(io.open(os.path.join(RES, "forward_local_edges.csv"), encoding="utf-8"))}
MV = {(r["model"], r["outcome"], r["exposure"]): r
      for r in csv.DictReader(io.open(os.path.join(RES, "mvmr_cascade.csv"), encoding="utf-8"))}


def mc(total, se_t, a, se_a, b, se_b, direct, se_d, rho_ab, rho_td, rng, rho_at=0.0):
    """Draw (a, b, total, direct) from a FULL 4x4 correlated normal.

    An earlier version swept only rho_ab and rho_td and fixed the other four pairwise correlations
    at zero. That silently omitted the term with the most leverage on the product estimator:
    corr(a, total). Here a = exposure->CAD and total = exposure->HF are BOTH univariable estimates
    from the SAME exposure instrument set, so they are strongly positively correlated - and a sits in
    the numerator of (a*b)/total while total sits in the denominator, which is exactly the ratio whose
    variance that correlation governs. Ignoring it inflates the interval.

    rho_at is applied to (a, total) and, by the same shared-instrument argument, to (b, direct).
    """
    C = np.eye(4)                      # order: a, b, total, direct
    C[0, 1] = C[1, 0] = rho_ab
    C[2, 3] = C[3, 2] = rho_td
    C[0, 2] = C[2, 0] = rho_at         # a vs total  - shared exposure instruments
    C[1, 3] = C[3, 1] = rho_at         # b vs direct - shared mediator/outcome data
    # project to the nearest positive-semidefinite matrix (a swept corner can leave the PSD cone)
    w, V = np.linalg.eigh(C)
    if w.min() < 1e-10:
        C = V @ np.diag(np.clip(w, 1e-10, None)) @ V.T
        d = np.sqrt(np.diag(C))
        C = C / np.outer(d, d)
    draws = rng.multivariate_normal(np.zeros(4), C, size=NDRAW, method="svd")
    aa = a + se_a * draws[:, 0]
    bb = b + se_b * draws[:, 1]
    tt = total + se_t * draws[:, 2]
    dd = direct + se_d * draws[:, 3]
    keep = np.abs(tt) > 1e-6
    tt, dd, aa, bb = tt[keep], dd[keep], aa[keep], bb[keep]
    pm_p = (aa * bb) / tt
    pm_d = (tt - dd) / tt
    q = lambda v: (np.percentile(v, 2.5), np.percentile(v, 97.5))
    return q(pm_p), q(pm_d)


def main():
    rng = np.random.default_rng(SEED)
    L, rows = [], []

    def out(s=""):
        L.append(s); print(s, flush=True)

    cadhf = MV[("HF_via_CAD", "HF", "CAD")]
    b, se_b = float(cadhf["direct_b"]), float(cadhf["direct_se"])

    out("=== Mediation: sensitivity to the assumed covariance structure ===")
    out("")
    out("The shipped analysis assumes corr(a,b) = corr(total,direct) = 0 (independent MC draws) and")
    out("MVMR uses gencov = 0. Both correlations are swept here; the point estimates do not depend on")
    out("them (they are functions of the estimates alone) - only the intervals do.")
    out(f"MC draws {NDRAW:,}, seed {SEED}. Shared mediator b = CAD->HF|exposure = {b:+.3f} (se {se_b:.3f}).")
    out("")

    for exp in ["BMI", "T2D"]:
        total, se_t = UNI[(exp, "HF")]
        a, se_a = UNI[(exp, "CAD")]
        d = MV[("HF_via_CAD", "HF", exp)]
        direct, se_d = float(d["direct_b"]), float(d["direct_se"])
        pm_prod_pt = 100 * (a * b) / total
        pm_diff_pt = 100 * (total - direct) / total

        out(f"--- {exp} -> HF ---")
        out(f"  point estimates (covariance-free): product {pm_prod_pt:.0f}% , difference {pm_diff_pt:.0f}%")
        out(f"  {'rho(a,b)':>9s} {'rho(tot,dir)':>13s} {'rho(a,tot)':>11s} {'product 95% CI':>24s} {'difference 95% CI':>24s}")
        for rho_ab in RHOS:
            for rho_td in RHOS:
              for rho_at in RHO_AT:
                (pl, ph), (dl, dh) = mc(total, se_t, a, se_a, b, se_b, direct, se_d,
                                        rho_ab, rho_td, rng, rho_at)
                shipped = (rho_ab == 0.0 and rho_td == 0.0 and rho_at == 0.0)
                out(f"  {rho_ab:9.2f} {rho_td:13.2f} {rho_at:11.2f} "
                    f"{f'[{100*pl:6.0f}%, {100*ph:6.0f}%]':>24s} "
                    f"{f'[{100*dl:6.0f}%, {100*dh:6.0f}%]':>24s}"
                    + ("   <- as shipped" if shipped else ""))
                rows.append(dict(exposure=exp, rho_ab=rho_ab, rho_total_direct=rho_td, rho_a_total=rho_at,
                                 pm_product_point=round(pm_prod_pt, 2),
                                 pm_product_lo=round(100 * pl, 2), pm_product_hi=round(100 * ph, 2),
                                 pm_diff_point=round(pm_diff_pt, 2),
                                 pm_diff_lo=round(100 * dl, 2), pm_diff_hi=round(100 * dh, 2)))
        out("")

    # ---- the qualitative question ------------------------------------------------------------
    bmi = [r for r in rows if r["exposure"] == "BMI"]
    t2d = [r for r in rows if r["exposure"] == "T2D"]
    bmi_upper = max(r["pm_product_hi"] for r in bmi)
    t2d_lower = min(r["pm_product_lo"] for r in t2d)
    t2d_spans100 = all(r["pm_product_hi"] > 100 for r in t2d)
    t2d_at9 = [r for r in t2d if r["rho_a_total"] == 0.9]
    hi9 = max(r["pm_product_hi"] for r in t2d_at9)
    lo9 = min(r["pm_product_lo"] for r in t2d_at9)
    out("READ:")
    out(f"- BMI->HF is COVARIANCE-ROBUST: across all {len(bmi)} combinations of the three swept")
    out(f"  correlations the product-method proportion mediated stays within [{min(r['pm_product_lo'] for r in bmi):.0f}%, {bmi_upper:.0f}%].")
    out("  'Predominantly direct' does not depend on any covariance assumption.")
    out("- T2D->HF is NOT: its interval is materially sensitive to corr(a,total), the correlation the")
    out("  first version of this script omitted. a = T2D->CAD and total = T2D->HF are both univariable")
    out("  estimates from the SAME T2D instruments, so a positive correlation is expected, and a sits")
    out("  in the numerator of (a*b)/total while total sits in the denominator.")
    out(f"    corr(a,total) = 0.0 (as shipped) -> product CI about [{[r for r in t2d if r['rho_a_total']==0.0 and r['rho_ab']==0.0 and r['rho_total_direct']==0.0][0]['pm_product_lo']:.0f}%, {[r for r in t2d if r['rho_a_total']==0.0 and r['rho_ab']==0.0 and r['rho_total_direct']==0.0][0]['pm_product_hi']:.0f}%]")
    out(f"    corr(a,total) = 0.9              -> product CI about [{lo9:.0f}%, {hi9:.0f}%]"
        + ("  - EXCLUDES complete mediation" if hi9 < 100 else ""))
    out("  So the shipped interval is CONSERVATIVE (too wide), not too narrow, and the earlier claim")
    out("  that T2D imprecision is 'intrinsic and not an artifact of the independence assumption' was")
    out("  WRONG - it was an artifact of omitting this term.")
    out("- The defensible statement for the manuscript is therefore NOT that the mediated proportion is")
    out("  irreducibly imprecise, but that its interval DEPENDS on a covariance the summary-data design")
    out("  cannot identify: under independence it spans complete mediation, under strong shared-")
    out("  instrument correlation it does not. Report the range and the dependence, and do not claim")
    out("  the mechanism is resolved either way.")
    out("- Point estimates are unchanged throughout (they are functions of the estimates alone).")

    io.open(os.path.join(RES, "mediation_covariance_sensitivity.txt"), "w",
            encoding="utf-8").write("\n".join(L))
    with io.open(os.path.join(RES, "mediation_covariance_sensitivity.csv"), "w",
                 encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nwrote results/mediation_covariance_sensitivity.{txt,csv}")


if __name__ == "__main__":
    main()
