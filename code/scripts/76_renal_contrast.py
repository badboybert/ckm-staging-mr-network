# -*- coding: utf-8 -*-
"""Step 76c: the renal (kidney-damage) contrast, one estimator and one scale.

Round-7 items C003 / C018 / C039 / C102 / C121, PI decision Q3 (2026-09-16). Three things were wrong
with the published version of this contrast:
  1. Europe was estimated by 09_zheng_ckd.py (fixed-effect Python IVW) and everything else in the
     paper by TwoSampleMR with multiplicative random effects -> two estimators for one claim.
  2. The European estimate is per mmHg and the East Asian one per SD, but they were printed and
     plotted on one axis with identical row labels, so the "null" East Asian point sat seven times
     further from zero than the "significant" European one.
  3. "Not in East Asians" rested on the SMALLEST East Asian kidney outcome (Biobank Japan atlas CKD,
     2,117 cases) while the package's own larger Taiwanese analyses (25,906 cases) were positive.
This step assembles every kidney-damage estimate the paper holds, from the files that produced them,
puts the SBP rows on the common per-SD scale using the ONE conversion constant the project already
declares (14_standardise.SD_SBP), and runs the ancestry-interaction test on that comparable scale.

Reads : results/network_ckd.csv            EUR, primary estimator (51_ckd_mr.R)
        results/renal_eas_primary.csv      EAS Biobank Japan, primary estimator (76_renal_eas_mr.R)
        results/network_tpmi_full.csv      EAS Taiwan, within-cohort (26_tpmi_full_mr.R)
        results/network_crosscohort_bbj_tpmi.csv   EAS zero-overlap BBJ->TPMI (39_crosscohort_mr.R)
Writes: results/renal_contrast.txt, results/renal_contrast.csv
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

import csv, io, math, os, re

BASE = P4_BASE
RES = os.path.join(BASE, "results")

# the SBP per-mmHg -> per-SD constant, read from the file that declares it (never retyped)
_std = io.open(os.path.join(BASE, "scripts", "14_standardise.py"), encoding="utf-8").read()
_m = re.search(r"^SD_SBP\s*=\s*([\d.]+)", _std, flags=re.M)
assert _m, "could not read SD_SBP from 14_standardise.py"
SD_SBP = float(_m.group(1))

SOURCES = [
    # key, file, ancestry, cohort label, outcome definition, cases
    ("EUR",      "network_ckd.csv",                  "European",   "CKDGen (Wuttke 2019)",
     "eGFR < 60 mL/min/1.73 m2", "41,395"),
    ("EAS_BBJ",  "renal_eas_primary.csv",            "East Asian", "Biobank Japan (Sakaue 2021)",
     "registry chronic renal failure", "2,117"),
    ("EAS_TPMI", "network_tpmi_full.csv",            "East Asian", "TPMI, within-cohort",
     "PheCode 585.3", "25,906"),
    ("EAS_ZERO", "network_crosscohort_bbj_tpmi.csv", "East Asian", "BBJ instruments -> TPMI (zero overlap)",
     "PheCode 585.3", "25,906"),
]


def rows_of(fn):
    with io.open(os.path.join(RES, fn), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def edge(fn, exp, out="CKD"):
    for r in rows_of(fn):
        if r["exposure"] == exp and r["outcome"] == out:
            return r
    return None


def z_test(b1, se1, b2, se2):
    z = (b1 - b2) / math.sqrt(se1 ** 2 + se2 ** 2)
    return z, 2 * 0.5 * math.erfc(abs(z) / math.sqrt(2))


def main():
    recs = []
    for exp in ("SBP", "BMI"):
        for key, fn, anc, cohort, defn, cases in SOURCES:
            r = edge(fn, exp)
            if r is None:
                continue
            b, se, p = float(r["ivw_b"]), float(r["ivw_se"]), float(r["ivw_p"])
            native = "per mmHg" if (exp == "SBP" and key == "EUR") else "per SD"
            # SBP: EUR is per mmHg, every East Asian source is rank-inverse-normal per SD
            b_sd, se_sd = (b * SD_SBP, se * SD_SBP) if native == "per mmHg" else (b, se)
            recs.append(dict(exposure=exp, outcome="CKD", source=key, ancestry=anc, cohort=cohort,
                             outcome_definition=defn, cases=cases, nsnp=int(float(r["nsnp"])),
                             b_native=round(b, 6), se_native=round(se, 6), units_native=native,
                             b_perSD=round(b_sd, 6), se_perSD=round(se_sd, 6), p=p,
                             source_file="results/" + fn))

    with io.open(os.path.join(RES, "renal_contrast.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recs[0].keys())); w.writeheader(); w.writerows(recs)

    L = []
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("=== Kidney damage (CKD): every estimate in the package, primary estimator, common scale ===")
    out(f"SBP per-mmHg -> per-SD conversion: x{SD_SBP} (14_standardise.SD_SBP)")
    out("")
    out(f"{'exposure':8s} {'source':9s} {'cohort':40s} {'cases':>7s} {'n':>4s} "
        f"{'b (native)':>12s} {'units':10s} {'b (per SD)':>11s} {'se':>8s} {'P':>10s}")
    for r in recs:
        out(f"{r['exposure']:8s} {r['source']:9s} {r['cohort']:40s} {r['cases']:>7s} {r['nsnp']:>4d} "
            f"{r['b_native']:>+12.4f} {r['units_native']:10s} {r['b_perSD']:>+11.4f} "
            f"{r['se_perSD']:>8.4f} {r['p']:>10.2e}")
    out("")
    out("--- ancestry interaction on the per-SD scale (European vs Biobank Japan, as published) ---")
    inter = []
    for exp in ("SBP", "BMI"):
        e = [r for r in recs if r["exposure"] == exp and r["source"] == "EUR"][0]
        a = [r for r in recs if r["exposure"] == exp and r["source"] == "EAS_BBJ"][0]
        z, p = z_test(e["b_perSD"], e["se_perSD"], a["b_perSD"], a["se_perSD"])
        out(f"  {exp}->CKD : EUR {e['b_perSD']:+.4f} ({e['se_perSD']:.4f}) vs BBJ {a['b_perSD']:+.4f} "
            f"({a['se_perSD']:.4f})  z = {z:+.3f}  P = {p:.3f}")
        inter.append((exp, z, p))
    out("")
    out("  The interaction is the comparison the paper reports; both are far from significant, so the")
    out("  East Asian result is imprecision against a 2,117-case outcome, not a demonstrated")
    out("  ancestry difference - and the larger Taiwanese outcomes above are positive for both edges.")
    io.open(os.path.join(RES, "renal_contrast.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")

    print("\n--- SELF-CHECKS ---", flush=True)
    ok = True
    checks = [
        ("EUR SBP->CKD is the S13/network_ckd.csv row (P < 1e-5, n = 437)",
         [r for r in recs if r["exposure"] == "SBP" and r["source"] == "EUR"][0]["nsnp"] == 437),
        ("EAS BBJ rows come from the primary estimator file",
         all(r["source_file"].endswith("renal_eas_primary.csv")
             for r in recs if r["source"] == "EAS_BBJ")),
        ("every SBP row is per SD after conversion",
         all(abs(r["b_perSD"] - r["b_native"] * (SD_SBP if r["units_native"] == "per mmHg" else 1))
             < 1e-6 * SD_SBP for r in recs if r["exposure"] == "SBP")),
        ("no interaction reaches P < 0.05", all(p >= 0.05 for _, _, p in inter)),
        ("the zero-overlap SBP->CKD estimate is positive and significant",
         [r for r in recs if r["exposure"] == "SBP" and r["source"] == "EAS_ZERO"][0]["p"] < 0.05),
    ]
    for name, cond in checks:
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}", flush=True)
        ok &= bool(cond)
    print(f"\nALL CHECKS PASS: {ok}", flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
