# -*- coding: utf-8 -*-
"""Step 67: does the cross-ancestry interaction family survive absorbing a single global scale offset?

Adversarial verification of step 60 found that the z-test's null (b_EUR = b_EAS on a common scale) is
already FALSE FAMILY-WIDE, for a reason that has nothing to do with any individual edge: one global
multiplicative offset fits the whole 67-edge family. If East Asian estimates are systematically a
constant fraction of European ones - because Biobank Japan rank-inverse-normalises skewed traits, or
because the EAS arm is one-sample within-cohort MR against a two-sample EUR arm, or both - then a
per-edge "significant interaction" may be reading that shared offset rather than effect modification.

This script fits the offset and re-runs the family with it absorbed:

    global:   b_EAS = s * b_EUR      (inverse-variance weighted, through the origin)
    absorbed: z = (s*b_EUR - b_EAS) / sqrt((s*se_EUR)^2 + se_EAS^2)

and reports BOTH survivor sets plus their overlap. It also reports the per-exposure median
|b_EAS|/|b_EUR|, because a 3-fold spread among traits the scale dictionary marks "Yes - both per-SD"
is itself evidence that the comparability rows are asserted from GWAS documentation rather than
checked against the data.

This does NOT decide which analysis is right. A global offset could be a real shared biological or
design difference (in which case absorbing it is correct and the residual interactions are the honest
set) or an artifact of the fit (in which case the raw set is). The point is that the two disagree, so
neither can be reported as established without saying so.

Out: results/interaction_global_scale.{txt,csv}
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

import csv, io, os, math, statistics, collections

BASE = P4_BASE
RES = os.path.join(BASE, "results")


def z_p(be, se_e, ba, se_a):
    z = (be - ba) / math.sqrt(se_e ** 2 + se_a ** 2)
    return z, 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def main():
    rows = [r for r in csv.DictReader(io.open(os.path.join(RES, "ancestry_interaction_family.csv"),
                                              encoding="utf-8"))
            if r["magnitude_eligible"] == "True"]
    for r in rows:
        r["be"] = float(r["b_eur_used"]); r["see"] = float(r["se_eur_used"])
        r["ba"] = float(r["b_eas"]);      r["sea"] = float(r["se_eas"])
    K = len(rows)
    BONF = 0.05 / K

    # inverse-variance weighted regression through the origin: s = sum(w*be*ba)/sum(w*be^2)
    w = [1.0 / (r["sea"] ** 2) for r in rows]
    num = sum(wi * r["be"] * r["ba"] for wi, r in zip(w, rows))
    den = sum(wi * r["be"] ** 2 for wi, r in zip(w, rows))
    s = num / den
    se_s = math.sqrt(1.0 / den)
    z_s = (s - 1.0) / se_s

    # regression-dilution check: is the slope just measurement error in b_EUR?
    var_be = statistics.pvariance([r["be"] for r in rows])
    mean_see2 = statistics.fmean([r["see"] ** 2 for r in rows])
    dilution = var_be / (var_be - mean_see2) if var_be > mean_see2 else float("nan")

    L = []
    def out(t=""):
        L.append(t); print(t, flush=True)

    out("=== Does the interaction family survive a single global EUR->EAS scale offset? ===")
    out("")
    out(f"Family: {K} magnitude-eligible edges; Bonferroni 0.05/{K} = {BONF:.3e}")
    out(f"Global IV-weighted slope through the origin:  b_EAS = {s:.4f} x b_EUR   "
        f"(SE {se_s:.4f}, z vs 1 = {z_s:.2f})")
    out(f"Regression-dilution correction factor: {dilution:.4f}  -> dilution-corrected slope "
        f"{s*dilution:.4f}")
    out("  (a correction factor near 1 means the offset is NOT explained by measurement error in b_EUR)")
    out("")

    med = collections.defaultdict(list)
    for r in rows:
        if r["be"]:
            med[r["exposure"]].append(abs(r["ba"]) / abs(r["be"]))
    out("Per-exposure median |b_EAS| / |b_EUR| (all of these traits are marked comparable):")
    for e in sorted(med, key=lambda k: statistics.median(med[k])):
        out(f"    {e:8s} {statistics.median(med[e]):6.2f}   (n={len(med[e])})")
    out("")

    raw, absorbed = [], []
    for r in rows:
        zr, pr = z_p(r["be"], r["see"], r["ba"], r["sea"])
        za, pa = z_p(s * r["be"], s * r["see"], r["ba"], r["sea"])
        r["z_raw"], r["p_raw"], r["z_abs"], r["p_abs"] = zr, pr, za, pa
        if pr < BONF:
            raw.append(r)
        if pa < BONF:
            absorbed.append(r)
    key = lambda r: r["exposure"] + "->" + r["outcome"]
    sraw, sabs = {key(r) for r in raw}, {key(r) for r in absorbed}

    out(f"Bonferroni survivors, RAW (as step 60 reports)     : {len(sraw)}")
    out(f"Bonferroni survivors, GLOBAL SCALE ABSORBED        : {len(sabs)}")
    out(f"Overlap                                            : {len(sraw & sabs)}")
    out(f"  lost when the offset is absorbed  ({len(sraw - sabs)}): {sorted(sraw - sabs)}")
    out(f"  gained when the offset is absorbed ({len(sabs - sraw)}): {sorted(sabs - sraw)}")
    out("")
    flips_raw = [r for r in raw if (r["be"] > 0) != (r["ba"] > 0)]
    flips_abs = [r for r in absorbed if (s * r["be"] > 0) != (r["ba"] > 0)]
    out(f"Opposite-sign among survivors: raw {len(flips_raw)}, absorbed {len(flips_abs)}")
    out("")
    out("READ:")
    out(f"- The two analyses share only {len(sraw & sabs)} of {len(sraw)} / {len(sabs)} edges. The IDENTITY of the")
    out("  'significant ancestry interactions' is therefore NOT robust to a single nuisance parameter,")
    out("  and no per-edge interaction claim should be made without disclosing this.")
    out(f"- The offset is large ({s:.2f}) and precisely estimated (z = {z_s:.1f} against 1). It is not")
    out(f"  regression dilution (correction factor {dilution:.3f}).")
    out("- The per-exposure spread above is the diagnostic that matters for the scale dictionary: traits")
    out("  it marks equally comparable behave very differently, which is what rank-inverse-normalising")
    out("  skewed traits (TG, HDL) against per-SD European traits would produce.")
    out("- The manuscript already applies exactly this logic to SBP->CKD ('a power difference, not a")
    out("  proven divergence'). Consistency requires extending it to the rest of the family.")

    io.open(os.path.join(RES, "interaction_global_scale.txt"), "w", encoding="utf-8").write("\n".join(L))
    with io.open(os.path.join(RES, "interaction_global_scale.csv"), "w", encoding="utf-8", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["exposure", "outcome", "b_eur_used", "se_eur_used", "b_eas", "se_eas",
                     "z_raw", "p_raw", "z_absorbed", "p_absorbed", "sig_raw", "sig_absorbed",
                     "global_slope", "family_size", "bonferroni"])
        for r in sorted(rows, key=lambda x: x["p_raw"]):
            wr.writerow([r["exposure"], r["outcome"], f"{r['be']:.6f}", f"{r['see']:.6f}",
                         f"{r['ba']:.6f}", f"{r['sea']:.6f}", f"{r['z_raw']:.4f}",
                         f"{r['p_raw']:.6e}", f"{r['z_abs']:.4f}", f"{r['p_abs']:.6e}",
                         r["p_raw"] < BONF, r["p_abs"] < BONF, f"{s:.6f}", K, f"{BONF:.6e}"])
    print("\nwrote results/interaction_global_scale.{txt,csv}")


if __name__ == "__main__":
    main()
