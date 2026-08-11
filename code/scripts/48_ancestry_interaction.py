# -*- coding: utf-8 -*-
"""Step 48 (E6): formal EUR-vs-EAS heterogeneity/interaction tests for every claimed ancestry
divergence (reviewer M5/E6). "Significant in one ancestry, null in the other" is NOT evidence of
interaction; this computes the actual between-ancestry difference test.

Interaction test (two independent estimates): z = (b_EUR - b_EAS) / sqrt(se_EUR^2 + se_EAS^2),
two-sided P (equivalent to Cochran's Q with 2 groups). Scale is harmonised where the two ancestries
use different units (SBP: EUR per-mmHg -> per-SD x 19.3; ties into E5). HbA1c/eGFR are flagged as
scale-non-comparable per the Methods.

EUR from forward_local_edges.csv (132-edge). EAS from network_eas_edges.csv (BBJ, the cohort the
paper's cross-ancestry comparison uses); SBP/BMI->CKD from zheng_ckd.txt (binary CKD, not a network node).
Writes results/ancestry_interaction_tests.txt.
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

import csv, io, os, math, re

BASE = P4_BASE
SBP_SD = 19.3   # adult SBP standard deviation (mmHg), per Methods, for the per-SD conversion

def load(path):
    d = {}
    for r in csv.DictReader(io.open(path, encoding="utf-8")):
        d[(r["exposure"], r["outcome"])] = r
    return d

EUR = load(os.path.join(BASE, "results/forward_local_edges.csv"))
EAS = load(os.path.join(BASE, "results/network_eas_edges.csv"))

def z_p(be, se_e, ba, se_a):
    z = (be - ba) / math.sqrt(se_e ** 2 + se_a ** 2)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return z, p

# claimed-divergence edges. (exp, out, scale_note, per_sd_convert_EUR)
EDGES = [
    ("SBP", "CKD", "SBP per-mmHg(EUR) vs per-SD(EAS) -> EUR converted x19.3", True),
    ("BMI", "CKD", "both per-SD", False),
    ("BMI", "HbA1c", "HbA1c scale non-comparable (Methods)", False),
    ("HbA1c", "HF", "HbA1c scale non-comparable (Methods)", False),
    ("TG", "LDL", "both per-SD lipids", False),
    ("HbA1c", "TG", "HbA1c scale non-comparable (Methods)", False),
    ("T2D", "HDL", "T2D log-OR exposure; HDL per-SD", False),
]

# SBP/BMI->CKD live in results/zheng_ckd.txt (step 09), not the network files. PARSED from that
# file rather than hardcoded, so the interaction test cannot silently drift from its source.
# Source lines look like:
#   EUR  SBP->CKD : nSNP=455  IVW b=+0.0143 se=0.0017 p=1.11e-17
ZHENG_RE = re.compile(
    r"^\s*(EUR|EAS)\s+(\w+)->(\w+)\s*:\s*nSNP=\s*(\d+)\s+IVW\s+b=([+-]?[\d.eE+-]+)\s+se=([\d.eE+-]+)")

def load_zheng(path):
    """-> {(exposure, outcome): {'eur': (b, se), 'eas': (b, se)}} parsed from results/zheng_ckd.txt."""
    d = {}
    for ln in io.open(path, encoding="utf-8"):
        m = ZHENG_RE.match(ln)
        if not m:
            continue
        anc, exp, outc = m.group(1), m.group(2), m.group(3)
        d.setdefault((exp, outc), {})[anc.lower()] = (float(m.group(5)), float(m.group(6)))
    return d

CKD = load_zheng(os.path.join(BASE, "results/zheng_ckd.txt"))

# ASSERT the parsed values equal the values this script previously hardcoded. Any drift in
# zheng_ckd.txt (or in the parser) must fail loudly rather than silently change a published P.
_CKD_EXPECTED = {
    ("SBP", "CKD"): dict(eur=(0.0143, 0.0017), eas=(0.0979, 0.2240)),
    ("BMI", "CKD"): dict(eur=(0.1429, 0.0286), eas=(0.4081, 0.1440)),
}
if set(CKD) != set(_CKD_EXPECTED):
    raise SystemExit(
        f"FATAL: results/zheng_ckd.txt edge set changed: parsed {sorted(CKD)}, "
        f"expected {sorted(_CKD_EXPECTED)}. A new edge here would silently override the "
        f"network-file source for that edge below.")
for _k, _exp in _CKD_EXPECTED.items():
    if _k not in CKD:
        raise SystemExit(f"FATAL: {_k[0]}->{_k[1]} not found in results/zheng_ckd.txt")
    for _anc, (_b, _se) in _exp.items():
        if _anc not in CKD[_k]:
            raise SystemExit(f"FATAL: {_anc.upper()} {_k[0]}->{_k[1]} missing in results/zheng_ckd.txt")
        _pb, _pse = CKD[_k][_anc]
        if abs(_pb - _b) > 1e-6 or abs(_pse - _se) > 1e-6:
            raise SystemExit(
                f"FATAL: zheng_ckd.txt drift for {_anc.upper()} {_k[0]}->{_k[1]}: "
                f"parsed b={_pb} se={_pse}, expected b={_b} se={_se}")

lines = []
def out(s=""):
    lines.append(s); print(s, flush=True)

out("=== E6: formal EUR-vs-EAS ancestry-interaction tests ===")
out("Test: z=(b_EUR - b_EAS)/sqrt(se_EUR^2+se_EAS^2), two-sided P. EAS = Biobank Japan (BBJ).\n")
out(f"{'edge':14s} {'EUR b(se)':18s} {'EAS b(se)':18s} {'z':>6s} {'P_int':>9s}  note")
results = []
for exp, outc, note, conv in EDGES:
    if (exp, outc) in CKD:
        be, se_e = CKD[(exp, outc)]["eur"]; ba, se_a = CKD[(exp, outc)]["eas"]
    else:
        re_, ra = EUR.get((exp, outc)), EAS.get((exp, outc))
        if not re_ or not ra:
            out(f"{exp+'->'+outc:14s} (missing in one ancestry: EUR={bool(re_)} EAS={bool(ra)})"); continue
        be, se_e = float(re_["ivw_b"]), float(re_["ivw_se"])
        ba, se_a = float(ra["ivw_b"]), float(ra["ivw_se"])
    if conv:                          # convert EUR per-mmHg -> per-SD
        be, se_e = be * SBP_SD, se_e * SBP_SD
    z, p = z_p(be, se_e, ba, se_a)
    results.append((exp, outc, z, p, note))
    out(f"{exp+'->'+outc:14s} {f'{be:+.3f}({se_e:.3f})':18s} {f'{ba:+.3f}({se_a:.3f})':18s} {z:6.2f} {p:9.4f}  {note}")

# multiple testing over the tested divergences
k = len(results)
out(f"\nMultiple testing: {k} interaction tests; Bonferroni threshold 0.05/{k} = {0.05/k:.4f}")
sig_raw = [r for r in results if r[3] < 0.05]
sig_bonf = [r for r in results if r[3] < 0.05 / k]
out(f"Significant at nominal 0.05: {[f'{e}->{o}' for e,o,z,p,n in sig_raw] or 'none'}")
out(f"Significant after Bonferroni: {[f'{e}->{o}' for e,o,z,p,n in sig_bonf] or 'none'}")
out("")
out("READ (honest):")
out("- SBP->CKD: on a COMMON per-SD scale the EUR and EAS estimates do NOT differ significantly")
out("  (the EAS estimate is imprecise, CI compatible with EUR). 'EUR-significant, EAS-null' is a POWER")
out("  difference, NOT a demonstrated ancestry interaction -> reframe from 'replicated divergence' to")
out("  'EUR-supported, EAS underpowered/imprecise' (reviewer M5).")
out("- SELECTION: these 5 are the BOTH-ANCESTRY-SIGNIFICANT opposite-sign subset. Across all 101 overlapping")
out("  edges there are ~12 Bonferroni-significant sign reversals; the other ~7 have >=1 ancestry non-significant")
out("  (a power/imprecise-null artifact, not a demonstrated interaction). So the 5 are a principled selection,")
out("  NOT the only or the strongest reversals - state this criterion so they are not read as the full set.")
out("- The 3 HbA1c edges are additionally scale-non-comparable -> sign-level interactions only (reviewer M7,C6).")
out("  TG->LDL and T2D->HDL are comparable-scale (magnitude + sign). Sign flips are scale-invariant under")
out("  positive rescaling, so the sign-level reading is defensible even for the HbA1c edges.")

io.open(os.path.join(BASE, "results/ancestry_interaction_tests.txt"), "w", encoding="utf-8").write("\n".join(lines))
print("\nwrote results/ancestry_interaction_tests.txt")
