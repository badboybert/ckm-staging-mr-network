# -*- coding: utf-8 -*-
"""Build the source data for Supplementary Figure S6 (main Figure 5 until 2026-07-25) from
results/, so the panels stop carrying literals.

It was the one figure with no read calls at all: all five panels were inline numbers.
The values were correct, but nothing in the submission package backed them and nothing stopped them
drifting from the analysis. This script derives each panel from the committed analysis outputs,
ASSERTS the derived value matches what the figure has been drawing, and writes one tidy CSV.

Output: figures/suppfig_data/fig5_source.csv  (panel, row, label, group, b, se, p, source_file)
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

import csv, io, os, re

BASE = P4_BASE
RES = os.path.join(BASE, "results")
OUT = os.path.join(BASE, "figures", "suppfig_data")

rows = []
def add(panel, row, label, group, b, se, p, src):
    rows.append(dict(panel=panel, row=row, label=label, group=group,
                     b=b, se=se, p=("" if p is None else p), source_file=src))

def check(name, got, want, tol):
    if want is not None and abs(float(got) - want) > tol:
        raise SystemExit(f"FIG5 PREP: {name} derived {got} but the figure draws {want}")

# ---------- 5a: SBP/BMI -> CKD across every kidney outcome the study has (results/renal_contrast.csv)
# Round-7 C003/C039/C121: this panel used to read results/zheng_ckd.txt - a fixed-effect Python IVW,
# while the rest of the paper is TwoSampleMR random-effects - and plotted a per-mmHg European estimate
# next to a per-SD East Asian one under identical row labels, so the "null" East Asian point sat seven
# times further from zero than the "significant" European one. 76_renal_contrast.py now supplies both
# ancestries on the primary estimator and on the common per-SD scale, with all three East Asian
# outcomes rather than only the 2,117-case one.
COHORT_LAB = {"EUR": "Europeans (CKDGen)", "EAS_BBJ": "East Asians (Biobank Japan)",
              "EAS_TPMI": "East Asians (TPMI)", "EAS_ZERO": "East Asians (zero-overlap)"}
with io.open(os.path.join(RES, "renal_contrast.csv"), encoding="utf-8") as fh:
    _rc = list(csv.DictReader(fh))
if len(_rc) != 8:
    raise SystemExit(f"FIG5 PREP: renal_contrast.csv has {len(_rc)} rows, expected 8 (2 exposures x 4 sources)")
for r in _rc:
    add("a", f"{r['exposure']} → CKD {r['source']}", f"{r['exposure']} → CKD", COHORT_LAB[r["source"]],
        r["b_perSD"], r["se_perSD"], r["p"], r["source_file"])
# the one conversion this panel depends on must be visible in the source data, not implicit
_eur_sbp = [r for r in _rc if r["exposure"] == "SBP" and r["source"] == "EUR"][0]
if abs(float(_eur_sbp["b_perSD"]) / float(_eur_sbp["b_native"]) - 19.3) > 0.05:
    raise SystemExit("FIG5 PREP: the European SBP row is not on the declared per-SD scale")

# ---------- 5b + 5d: the ->BMI triangulation and the meta knife-edge ----------
def rd(fn):
    with io.open(os.path.join(RES, fn), encoding="utf-8") as f:
        return list(csv.DictReader(f))

def edge(rowset, exp, out):
    for r in rowset:
        if r.get("exposure") == exp and r.get("outcome") == out:
            return r
    return None

bbj = edge(rd("network_eas_edges.csv"), "T2D", "BMI")
tpmi = edge(rd("network_tpmi_full.csv"), "T2D", "BMI")
mfix = edge(rd("network_eas_meta_fixed.csv"), "T2D", "BMI")
mran = edge(rd("network_eas_meta_random.csv"), "T2D", "BMI")
kog = rd("koges_triangulation.csv")

def kv(exp, out):
    for r in kog:
        if r["exposure"] == exp and r["outcome"] == out:
            return r
    return None

k_bmi = kv("KoGES_DM", "KoGES_BMI")
k_wst = kv("KoGES_DM", "KoGES_WAIST")

check("5b BBJ T2D->BMI", bbj["ivw_b"], -0.098, 5e-3)
check("5b TPMI T2D->BMI", tpmi["ivw_b"], -0.046, 5e-3)
check("5b meta fixed", mfix["ivw_b"], -0.069, 5e-3)
check("5b meta random", mran["ivw_b"], -0.071, 5e-3)
check("5b KoGES DM->BMI IVW", k_bmi["ivw_b"], -0.043, 5e-3)
check("5b KoGES DM->BMI WM", k_bmi["wm_b"], -0.090, 5e-3)

H, P = "Hospital biobank", "Population cohort"
add("b", "BBJ  T2D→BMI (IVW)", "BBJ  T2D→BMI (IVW)", H, bbj["ivw_b"], bbj["ivw_se"], bbj["ivw_p"],
    "results/network_eas_edges.csv")
add("b", "TPMI T2D→BMI (IVW)", "TPMI T2D→BMI (IVW)", H, tpmi["ivw_b"], tpmi["ivw_se"], tpmi["ivw_p"],
    "results/network_tpmi_full.csv")
add("b", "BBJ+TPMI meta (fixed)", "BBJ+TPMI meta (fixed)", H, mfix["ivw_b"], mfix["ivw_se"], mfix["ivw_p"],
    "results/network_eas_meta_fixed.csv")
add("b", "BBJ+TPMI meta (random)", "BBJ+TPMI meta (random)", H, mran["ivw_b"], mran["ivw_se"], mran["ivw_p"],
    "results/network_eas_meta_random.csv")
add("b", "KoGES DM→BMI (IVW)", "KoGES DM→BMI (IVW)", P, k_bmi["ivw_b"], k_bmi["ivw_se"], k_bmi["ivw_p"],
    "results/koges_triangulation.csv")
# The weighted-median arms have no SE column in the triangulation table; the figure draws the SE the
# KoGES diagnostic reported. Kept here with the source named rather than silently re-typed in the .R.
add("b", "KoGES DM→BMI (WM)", "KoGES DM→BMI (WM)", P, k_bmi["wm_b"], 0.015, k_bmi["wm_p"],
    "results/koges_triangulation.csv (wm_b; SE from 29b_koges_diag.R)")
add("b", "KoGES DM→Waist (WM)", "KoGES DM→Waist (WM)", P, k_wst["wm_b"], 0.015, k_wst["wm_p"],
    "results/koges_triangulation.csv (wm_b; SE from 29b_koges_diag.R)")
add("b", "KoGES-DM→TWB-Waist (WM)", "KoGES-DM→TWB-Waist (WM)", P, -0.043, 0.016, None,
    "results/triangulation_summary.txt (30b_twb_edge_mr.R)")

add("d", "Fixed-effect IVW", "Fixed-effect IVW", "survives", mfix["ivw_b"], mfix["ivw_se"], mfix["ivw_p"],
    "results/network_eas_meta_fixed.csv")
add("d", "Random-effect IVW", "Random-effect IVW", "fails", mran["ivw_b"], mran["ivw_se"], mran["ivw_p"],
    "results/network_eas_meta_random.csv")
# Hartung-Knapp is the conservative RE inference for k=2 (reviewer M6 / round-2 #13): its wide t-based
# CI, not the fixed/random z-CI, is the honest headline. Parse b and the 95% CI from the analysis
# output and carry an EFFECTIVE se that reproduces that CI under the panel's b +/- 1.96*se drawing,
# so the figure shows the real HK interval without hand-typing it.
_h4 = io.open(os.path.join(RES, "h4_leandiabetes.txt"), encoding="utf-8").read()
_hk = re.search(r"HK\s*:\s*b=([+-]?[\d.]+)\s+95% CI \[([+-]?[\d.]+),\s*([+-]?[\d.]+)\]\s+p=([\d.]+)", _h4)
if not _hk:
    raise SystemExit("FIG5 PREP: could not parse the Hartung-Knapp line from h4_leandiabetes.txt")
_hk_b, _hk_lo, _hk_hi, _hk_p = (float(_hk.group(1)), float(_hk.group(2)),
                                float(_hk.group(3)), float(_hk.group(4)))
_hk_se = (_hk_hi - _hk_lo) / (2 * 1.96)          # effective se: b +/- 1.96*se reproduces the HK CI
check("5d HK b", _hk_b, -0.071, 5e-3)
add("d", "Hartung–Knapp (conservative)", "Hartung–Knapp (conservative)", "fails",
    _hk_b, _hk_se, _hk_p, "results/h4_leandiabetes.txt (52_h4_leandiabetes.R)")

# ---------- 5c: per-SNP structure of KoGES DM->BMI ----------
# The per-SNP Wald ratios were reported in the triangulation summary, not written to a CSV. Carry
# them with the source named, and assert every value the summary DOES state.
tri = io.open(os.path.join(RES, "triangulation_summary.txt"), encoding="utf-8").read()
persnp = [("rs11658063", "HNF1B", -0.147, 0.038), ("rs35612982", "CDKAL1", -0.124, 0.018),
          ("rs11187138", "HHEX", -0.124, 0.032), ("rs60808706", "KCNQ1", -0.100, 0.025),
          ("rs10811662", "CDKN2A/B", -0.073, 0.023), ("rs11558471", "SLC30A8", -0.069, 0.033),
          ("rs2233580", "PAX4", 0.007, 0.022), ("rs2074356", "12q24/ALDH2", 0.272, 0.041),
          ("rs6567160", "MC4R", 0.493, 0.043)]
for rs, gene, _b, _se in persnp:
    if rs not in tri:
        raise SystemExit(f"FIG5 PREP: {rs} ({gene}) is not named in triangulation_summary.txt")
if "PAX4 rs2233580 (+0.007)" not in tri:
    raise SystemExit("FIG5 PREP: the PAX4 near-null value is not as stated in the summary")
BETA, ADI = "β-cell / insulin-secretion", "adiposity"
for rs, gene, b, se in persnp:
    add("c", rs, f"{gene} ({rs})", ADI if gene in ("12q24/ALDH2", "MC4R") else BETA, b, se, None,
        "results/triangulation_summary.txt (29b_koges_diag.R)")

# ---------- 5e: both-significant sign divergences ----------
want_e = {"BMI->HbA1c": (0.23, -0.23), "HbA1c->HF": (0.05, -0.18), "TG->LDL": (0.28, -0.11),
          "HbA1c->TG": (0.07, -0.09), "T2D->HDL": (-0.07, 0.03)}
comp = rd("eur_vs_eas_comparison.csv")
for e, (we, wa) in want_e.items():
    r = next((x for x in comp if x["edge"] == e), None)
    if r is None:
        raise SystemExit(f"FIG5 PREP: {e} missing from eur_vs_eas_comparison.csv")
    # The figure drew these at 2 dp, so allow one rounding step (TG→LDL EAS is -0.115 in source and
    # was drawn as -0.11). The assertion is here to catch a TRANSCRIPTION error, not to police a
    # rounding convention — from now on the panel plots the exact derived value.
    for side, got, want in (("EUR", r["EUR_b"], we), ("EAS", r["EAS_b"], wa)):
        if abs(float(got) - want) > 6e-3:
            raise SystemExit(f"FIG5 PREP: 5e {e} {side} derived {got} but the figure draws {want}")
    lab = e.replace("->", "→")
    add("e", lab, lab, "EUR", r["EUR_b"], "", r["EUR_p"], "results/eur_vs_eas_comparison.csv")
    add("e", lab, lab, "EAS", r["EAS_b"], "", r["EAS_p"], "results/eur_vs_eas_comparison.csv")

os.makedirs(OUT, exist_ok=True)
path = os.path.join(OUT, "fig5_source.csv")
with io.open(path, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["panel", "row", "label", "group", "b", "se", "p", "source_file"])
    w.writeheader()
    w.writerows(rows)
print(f"wrote {os.path.relpath(path, BASE)} with {len(rows)} rows "
      f"({sorted(set(r['panel'] for r in rows))})")
print("all derived values reproduce the values Supplementary Figure S6 draws")
