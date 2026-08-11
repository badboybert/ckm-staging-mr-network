# -*- coding: utf-8 -*-
"""Step 33 (ITEM 1): summarise the ApoB atherogenic-axis MVMR -> resolves the HDL->CAD CAUSE false positive."""
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
rows = list(csv.DictReader(io.open(os.path.join(BASE, "results/mvmr_apob.csv"), encoding="utf-8")))
def g(model, exp, k):
    for r in rows:
        if r["model"] == model and r["exposure"] == exp: return r[k]
    return "NA"
out = io.StringIO()
def P(*a): print(*a, file=out)
P("="*92)
P("ITEM 1 — ApoB atherogenic-axis MVMR: resolving the HDL->CAD CAUSE false positive (Milestone 6.10)")
P("="*92)
P()
P("BACKGROUND: CAUSE returned HDL->CAD = CAUSAL (z=-2.46) = a FALSE POSITIVE at its documented boundary")
P("(proportional confounding via a correlated MARKER trait). The atherogenic axis (LDL/ApoB/TG) is the")
P("driver; HDL merely marks it. MVMR is the correct arbiter. ApoB (Sinnott-Armstrong 2021 UKB, N=435,744)")
P("added as an exposure; 124 common lipid instruments with complete ApoB/LDL/HDL/TG/CAD data (hg38->hg19).")
P()
P("Univariable (marginal) total effects on CAD:")
P("  ApoB +0.452 (p2e-19) | LDL +0.511 (p3e-53) | HDL -0.305 (p7e-30) | TG +0.393 (p6e-37)")
P()
P(f"{'model':<20}{'exposure':<7}{'direct b':>10}{'direct p':>11}{'condF':>8}   note")
P("-"*92)
def line(model, exp, note=""):
    P(f"{model:<20}{exp:<7}{float(g(model,exp,'direct_b')):>+10.3f}{float(g(model,exp,'direct_p')):>11.1e}"
      f"{float(g(model,exp,'cond_F')):>8.1f}   {note}")
P("# CLEAN atherogenic-axis test (ApoB represents the axis; NOT competing with its collinear twin LDL):")
line("ApoB_HDL_TG","ApoB","ApoB SURVIVES (strong condF 57.7)  <<< key")
line("ApoB_HDL_TG","HDL","124-SNP ApoB-subset: HDL -0.076 NS = POWER ARTIFACT (SE0.069); see 630-SNP below")
line("ApoB_HDL_TG","TG","attenuated, borderline")
P("# Reference axis-representative = LDL (same story, LDL as axis):")
line("LDL_HDL_TG","LDL","LDL SURVIVES (condF 132)")
line("LDL_HDL_TG","HDL","HDL attenuates to borderline-null")
line("LDL_HDL_TG","TG","attenuated, borderline")
P("# Why the joint 4-lipid model is uninformative for ApoB vs LDL (collinearity):")
line("ApoB_CAD","ApoB","ApoB condF collapses to 9.2 (weak) when LDL is co-modelled")
line("ApoB_CAD","LDL","LDL retains +0.51 but condF only 13.6")
line("ApoB_CAD","HDL","attenuates to -0.139 (p0.056)")
line("ApoB_LDL","ApoB","head-to-head: both weak; cannot separate ApoB from LDL")
line("ApoB_LDL","LDL","")
P("# ★ FULL-DENSITY reference (630 SNPs, NO ApoB-coverage restriction; verifier-prompted, independently reproduced):")
P(f"{'LDL_HDL_TG_630':<20}{'LDL':<7}{0.499:>+10.3f}{2.3e-47:>11.1e}{117.3:>8.1f}   LDL survives strongly")
P(f"{'LDL_HDL_TG_630':<20}{'HDL':<7}{-0.165:>+10.3f}{1.15e-5:>11.1e}{58.5:>8.1f}   ★ HDL RESIDUAL SURVIVES (not null)")
P(f"{'LDL_HDL_TG_630':<20}{'TG':<7}{0.140:>+10.3f}{9.0e-4:>11.1e}{47.6:>8.1f}   TG survives")
P("-"*92)
P("VERDICT (CALIBRATED after adversarial verification, wf_b1e5d2ce item1-lens = PARTIALLY_REFUTED):")
P("  1. The atherogenic lipoprotein axis is the causal lipid driver of CAD: ApoB direct +0.292 (p3.9e-5,")
P("     condF 57.7) representing the axis; LDL direct +0.499 (p2e-47, condF 117) at full density. Both survive")
P("     strongly. HONEST BOUND: ApoB & LDL are too collinear to separate (ApoB condF collapses to 9.2 in the")
P("     joint ApoB+LDL model) => claim is AXIS-level (Richardson 2020, Sniderman), NOT ApoB-over-LDL.")
P("  2. ★ HDL->CAD ATTENUATES ~50% under MVMR but does NOT dissolve: marginal -0.305 (p7e-30) -> conditional")
P("     -0.165 (p1.15e-5) at FULL instrument density (630-SNP LDL+HDL+TG). The near-null in the ApoB-restricted")
P("     124-SNP model (-0.076/-0.119) was a POWER ARTIFACT (HDL SE 0.069 vs 0.037 at full density; ApoB covers")
P("     only 143/817 union instruments). => MVMR shows HDL's marginal CAD signal is SUBSTANTIALLY (~46%)")
P("     atherogenic-axis confounding, but a RESIDUAL direct HDL->CAD association SURVIVES.")
P("  3. ★ CONSEQUENCE for the CAUSE HDL->CAD false-positive (6.10): MVMR does NOT by itself refute a causal HDL")
P("     contribution — a residual survives. The 'HDL is non-causal' negative-control status rests on the WIDER")
P("     evidence (Voight 2012 PMID 22607825 HDL genetic-score->MI null; failed HDL-raising trials [CETP inhib.];")
P("     Holmes 2015), which this MVMR is CONSISTENT with (large attenuation) but does not independently prove.")
P("     Honest read: 'MVMR attenuates HDL ~50%, consistent with major atherogenic-axis confounding of the CAUSE")
P("     verdict; full non-causality is established by external trial/MR evidence, not this analysis alone.'")
P("  4. TG retains a modest direct effect (+0.14, p9e-4 at full density). Does NOT impugn LDL->CAD (true")
P("     positive) or CAD->HF. No computational/harmonisation/liftover error (all reproduce; liftover is")
P("     hg19->hg38 to match ApoB then harmonise ApoB back to the hg19 panel allele; univariable ApoB uses 210).")
txt = out.getvalue()
with io.open(os.path.join(BASE, "results/apob_mvmr_summary.txt"), "w", encoding="utf-8") as f:
    f.write(txt)
print(txt)
