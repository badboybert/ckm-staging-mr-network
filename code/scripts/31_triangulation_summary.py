# -*- coding: utf-8 -*-
"""Step 31 (ITEM 3): consolidate the population-vs-hospital ->BMI triangulation into one table.
Pulls the T2D/DM -> adiposity edge across all EAS cohorts + designs and writes the verdict."""
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

BASE = P4_BASE

# hospital (disease-ascertained) cohorts, within-cohort T2D->BMI (from prior milestones/results)
HOSP = [
    ("BBJ  (hospital, JP)", "T2D->BMI",  127, -0.0977, 0.0188, 2.1e-7, None,   "IVW sig neg"),
    ("TPMI (hospital, TW)", "T2D->BMI",   56, -0.0456, 0.0169, 7.1e-3, None,   "IVW sig neg"),
]
# meta of the two hospital cohorts (from results/eas_meta_summary.txt)
META = ("BBJ+TPMI meta",    "T2D->BMI", 183, -0.069, 0.013, 4.5e-8, -0.071, "fixed sig / random p6.5e-3")

def read_koges():
    rows = {}
    with io.open(os.path.join(BASE, "results/koges_triangulation.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows[(r["exposure"], r["outcome"])] = r
    return rows
K = read_koges()

def g(edge, k): return float(K[edge][k])

out = io.StringIO()
def P(*a): print(*a, file=out)
P("="*96)
P("ITEM 3 — POPULATION vs HOSPITAL  ->adiposity  triangulation  (the decisive selection-vs-biology test)")
P("="*96)
P()
P("QUESTION: the negative T2D/DM -> BMI edge appears in BOTH hospital EAS cohorts (BBJ, TPMI).")
P("  Selection/collider bias (disease-ascertained biobank) predicts it VANISHES in population cohorts.")
P("  Real EAS lean-diabetes biology predicts it PERSISTS in population cohorts.")
P()
P(f"{'cohort / design':<26}{'edge':<12}{'nSNP':>5}{'IVW b':>9}{'IVW p':>10}{'WM b':>9}{'WM p':>10}  note")
P("-"*96)
for name, edge, n, b, se, p, wm, note in HOSP:
    P(f"{name:<26}{edge:<12}{n:>5}{b:>9.3f}{p:>10.1e}{'':>9}{'':>10}  {note}")
name, edge, n, b, se, p, wm, note = META
P(f"{name:<26}{edge:<12}{n:>5}{b:>9.3f}{p:>10.1e}{wm:>9.3f}{'':>10}  {note}")
P("-"*96)
# KoGES / TWB population edges
def line(name, exp, outc, note):
    e = (exp, outc)
    P(f"{name:<26}{exp.replace('KoGES_','').replace('TWB_','')+'->'+outc.replace('KoGES_','').replace('TWB_',''):<12}"
      f"{int(float(K[e]['nsnp'])):>5}{g(e,'ivw_b'):>9.3f}{g(e,'ivw_p'):>10.1e}{g(e,'wm_b'):>9.3f}{g(e,'wm_p'):>10.1e}  {note}")
line("KoGES(pop,KR) within",  "KoGES_DM", "KoGES_BMI",   "IVW null(het Q p1e-52) / WM sig NEG")
line("KoGES(pop,KR) within",  "KoGES_DM", "KoGES_WAIST", "IVW null / WM sig NEG")
# TWB two-sample (hard-coded from step 30 run; no koges_triangulation row)
P(f"{'KoGES-DM->TWB-Waist 2smpl':<26}{'DM->Waist':<12}{6:>5}{-0.0189:>9.3f}{0.441:>10.1e}{-0.0432:>9.3f}{0.0065:>10.1e}  no overlap; MC4R absent; WM sig NEG")
P("-"*96)
P("POSITIVE CONTROL (forward adiposity->diabetes, must be strongly POSITIVE):")
line("KoGES(pop,KR) within",  "KoGES_BMI", "KoGES_DM",   "matches TPMI +0.99 / EUR +0.99  [pipeline OK]")
P()
P("PER-SNP STRUCTURE of KoGES DM->BMI (9 analysed after harmonise; rs11257657 dropped as ambiguous palindrome):")
P("  6/9 instruments NEGATIVE = insulin-secretion/beta-cell T2D loci: CDKAL1 rs35612982, CDKN2A/B rs10811662,")
P("     HHEX rs11187138, KCNQ1 rs60808706, SLC30A8 rs11558471, HNF1B rs11658063.")
P("  1 near-null = PAX4 rs2233580 (+0.007). 2 POSITIVE outliers = adiposity loci: MC4R rs6567160 (+0.49),")
P("     12q24/ALDH2(HECTD4) rs2074356 (+0.27). Sign test 6/9 negative = exact binomial p=0.51 (NOT significant).")
P("  => IVW-random null; only median/mode go negative, and the negative cluster is MECHANISTICALLY CORRELATED")
P("     (all beta-cell loci) = one shared insulin->adiposity pathway represented many times, NOT an independent")
P("     valid majority. Q=264 (I2=0.97, p1.5e-52). LOO: no sign flip; dropping MC4R -> IVW-random -0.067 (p0.052).")
P()
P("VERDICT  (CALIBRATED after adversarial verification, wf_b1e5d2ce stats-lens = PARTIALLY_REFUTED):")
P("  The population KoGES DM->BMI IVW POINT ESTIMATE (-0.043) MATCHES the hospital TPMI edge (-0.046) in sign")
P("  and magnitude => the negative ->BMI edge is NOT FALSIFIED in a population EAS cohort; it does not vanish as")
P("  a pure selection/collider would predict. **HOSPITAL SELECTION IS THEREFORE NOT REQUIRED** to explain the")
P("  ->BMI edges, and they cannot be dismissed as a BBJ/TPMI ascertainment artifact.")
P("  BUT a CAUSAL EAS T2D->lower-BMI effect is **NOT ESTABLISHED**: the population estimate is null/underpowered")
P("  (IVW-random -0.043, p0.40 under I2=0.97, 9 SNPs); the negative sign is carried by median/mode whose SEs do")
P("  NOT absorb the Q=264 overdispersion (WM p2.7e-9 is anti-conservative; heterogeneity-inflated WM p~0.24);")
P("  under BALANCED pleiotropy (Egger-int p0.49, no directional pleiotropy) IVW-random is the honest primary,")
P("  and the negative is carried by a mechanistically-correlated beta-cell cluster = plausible pleiotropy, not a")
P("  demonstrated exclusion-restriction-valid effect.")
P()
P("  HONEST CEILING for the write-up:")
P("  'The negative ->BMI edge is not falsified in a population EAS cohort (KoGES point estimate matches the")
P("   hospital cohorts in sign and magnitude), so hospital selection is not required to explain it and it cannot")
P("   be dismissed as a disease-ascertainment artifact. However, population inference is null/underpowered and")
P("   the sign rests on robust estimators plus a mechanistically-correlated insulin-secretion locus cluster with")
P("   plausible pleiotropy, so a causal EAS T2D->lower-BMI effect is not established. The pattern is CONSISTENT")
P("   WITH the lean-diabetes hypothesis (Ke 2022) but does not prove it.'")
P()
P("  UPDATES INTERPRETATION.md 3: reclassify the ->BMI edges from 'artifact-tier (BBJ selection)' to")
P("  'not a demonstrable selection artifact; an underpowered, mechanistically-plausible ancestry signal")
P("   consistent with lean-diabetes biology' (a hedged Tier-2 hypothesis, NOT an asserted Tier-1 finding).")
P("   The AHA stage-1(adiposity)->stage-2(dysglycemia) ordering is neither cleanly portable nor cleanly")
P("   falsified in EAS. Positive control BMI->DM +0.893 (p5e-7) confirms the FORWARD axis ports; only the")
P("   backward ->BMI edge is the ancestry-suggestive, underpowered signal.")
txt = out.getvalue()
with io.open(os.path.join(BASE, "results/triangulation_summary.txt"), "w", encoding="utf-8") as f:
    f.write(txt)
print(txt)
