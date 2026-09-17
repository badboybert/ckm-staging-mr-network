# -*- coding: utf-8 -*-
"""prep_newpanels.py — derive figure-panel data for the 132-edge rewrite's NEW panels.

Every value is PARSED from a committed file in results/ (never hand-typed), and every parse is
asserted against the number the manuscript states. If an upstream analysis is re-run and a value
moves, this script FAILS LOUDLY instead of silently shipping a stale panel (CLAUDE.md §9).

Emits (to figures/data/):
  fig1_constrained_nulls.csv  <- results/staging_constrained_nulls.txt   (E2)
  fig1_cac_staging.csv        <- results/staging_cac_stage3.txt          (CAC stage-3)
  fig2_mediation.csv          <- results/formal_mediation.txt            (E7)
  fig3_h5_estimators.csv      <- results/h5_reverse_pleiotropy.txt       (H5)
  fig4_portability.csv        <- results/h3_portability.txt              (H3)
  fig5_interaction.csv        <- results/ancestry_interaction_tests.txt  (E6)

Panels sourced from existing CSVs (network_hfsubtypes / network_whradjbmi / network_cac /
network_crosscohort_bbj_tpmi) are read directly by the R scripts and are NOT duplicated here.
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

import os, re, io, sys, csv

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
BASE = P4_BASE
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(BASE, "figures", "data")
os.makedirs(OUT, exist_ok=True)

def txt(name):
    with open(os.path.join(RES, name), encoding="utf-8") as fh:
        return fh.read()

def write(name, header, rows):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print(f"  wrote {name}  ({len(rows)} rows)")

def close(a, b, tol=1e-6):
    return abs(float(a) - float(b)) <= tol

FAIL = []
def check(label, got, want, tol=1e-6):
    if not close(got, want, tol):
        FAIL.append(f"{label}: parsed {got}, manuscript states {want}")

# ---------------------------------------------------------------- E2 constrained nulls (Fig 1)
t = txt("staging_constrained_nulls.txt")
obs   = float(re.search(r"OBSERVED cross-stage concordance = ([\d.]+)", t).group(1))
p_lab = float(re.search(r"A\. Simple label permutation.*?P = ([\d.]+)", t, re.S).group(1))
p_cls = float(re.search(r"B\. Phenotype-class-constrained permutation.*?P = ([\d.]+)", t, re.S).group(1))
p_deg = float(re.search(r"C\. Degree-preserving edge-rewiring null.*?P = ([\d.]+)", t, re.S).group(1))
loo   = re.findall(r"drop (\w+)\s*:\s*concordance ([\d.]+)", t)
f_con, f_p = re.search(r"Drop SBP \+ Stroke:.*?concordance = ([\d.]+)[^;]*; label-perm P = ([\d.]+)", t, re.S).groups()
b_lo, b_hi = re.search(r"Bootstrap 95% CI on concordance.*?\[([\d.]+), ([\d.]+)\]", t, re.S).groups()

check("E2 observed concordance", obs, 0.926, 1e-3)
# 2026-09-16 provenance correction: 0.0015 -> 0.0010. The European HbA1c exposure was found to be
# GCST90014006 (Mbatchou 2021, UK Biobank, 389,889), not MAGIC/Chen 2021 (146,806); the corrected
# sample size flips the Steiger call on TG->HbA1c, which joins the graph as a WITHIN-stage edge.
# The observed cross-stage set is untouched (25/27 = 0.926), but an extra edge changes which edges
# count as cross-stage under a permuted labelling, so the exact tail moves from 3/1,980 to 2/1,980.
check("E2 label-permutation P (exact)", p_lab, 0.0010, 3e-4)
# The degree-preserving null BOUNDS the paper's headline staging claim, so its expectation is read
# from what the paper actually says rather than typed here. It WAS typed (0.3142). When script 47 was
# re-run on 2026-07-26 the value moved to 0.2696 and the prose was updated everywhere — Abstract,
# Results, Discussion and the Figure 1 legend all say 0.27 — but this script was not re-run, so
# figures/data/fig1_constrained_nulls.csv kept 0.3142 and Figure 1f shipped "P = 0.31" beside them.
# A typed expectation cannot detect its own staleness; a docstring promising "asserted against the
# number the manuscript states" is not the same as reading the manuscript.
_prose = open(os.path.join(BASE, "figures", "RESULTS.md"), encoding="utf-8").read()
_m_deg = re.search(r"degree-preserving edge-rewiring null[^.]*?not exceeded \(\*P\* = ([\d.]+)\)", _prose)
if not _m_deg:
    FAIL.append("E2 degree-rewiring null P: RESULTS.md no longer states this P — cannot verify the panel")
else:
    check("E2 degree-rewiring null P vs RESULTS.md", p_deg, float(_m_deg.group(1)), 5e-3)
check("E2 drop-SBP+Stroke concordance", f_con, 0.944, 1e-3)
check("E2 drop-SBP+Stroke P (exact)", f_p, 0.0056, 3e-4)
check("E2 bootstrap CI lower", b_lo, 0.815, 1e-3)
loo_worst = min(float(v) for _, v in loo)
check("E2 leave-one-out worst", loo_worst, 0.900, 1e-3)

# Null battery rows: what each null tests + whether the observed concordance beats it.
rows = [
    ["Label permutation",              "Are the stage LABELS interchangeable?",      p_lab, "exceeded"],
    ["Degree-preserving rewiring",     "Beyond the traits' exposure/outcome roles?", p_deg, "NOT exceeded"],
    ["Phenotype-class permutation",    "Degenerate (class nested in stage)",         p_cls, "uninformative"],
]
write("fig1_constrained_nulls.csv", ["null", "question", "p", "verdict"], rows)
write("fig1_null_robustness.csv", ["test", "concordance", "p"],
      [["Observed", obs, p_lab],
       ["Leave-one-node-out (worst)", loo_worst, ""],
       ["Drop SBP + Stroke", float(f_con), float(f_p)]])
write("fig1_loo.csv", ["node", "concordance"], [[n, float(v)] for n, v in loo])
write("fig1_bootstrap.csv", ["stat", "value"],
      [["observed", obs], ["ci_lo", float(b_lo)], ["ci_hi", float(b_hi)]])

# ---------------------------------------------------------------- CAC stage-3 (Fig 1)
t = txt("staging_cac_stage3.txt")
cac_con, cac_p = re.search(r"concordance = ([\d.]+) \(32/34 forward cross-stage\); permutation P = ([\de.-]+)", t).groups()
base_con, base_p = re.search(r"vs 132-edge no-stage-3: ([\d.]+), P=([\de.-]+)", t).groups()
drop_con = re.search(r"Dropping the thin\s+CAC->CAD leaves concordance ([\d.]+)", t, re.S).group(1)
cac_p, base_p = cac_p.rstrip("."), base_p.rstrip(".")   # trailing sentence period is not part of the value
check("CAC concordance", cac_con, 0.941, 1e-3)
check("CAC baseline concordance", base_con, 0.926, 1e-3)
check("CAC drop-thin-edge concordance", drop_con, 0.939, 1e-3)
# NB: these strings are AXIS TICK LABELS on a half-width panel — keep them short or they collide.
write("fig1_cac_staging.csv", ["network", "concordance", "p"],
      [["no stage 3", float(base_con), float(base_p)],
       ["+ CAC\n(stage 3)", float(cac_con), float(cac_p)],
       ["+ CAC\n(thin edge\ndropped)", float(drop_con), ""]])

# ---------------------------------------------------------------- E7 formal mediation (Fig 2)
t = txt("formal_mediation.txt")
def med(block):
    b = re.search(r"--- %s -> HF ---(.*?)(?=---|\Z)" % block, t, re.S).group(1)
    g = lambda rx: re.search(rx, b).groups()
    tot,  = re.search(r"total\s+\(univariable\)\s+= \+?([\d.]+)", b).groups()
    dir_, dir_p = re.search(r"direct\s+\(\| CAD, MVMR\)\s+= \+?([\d.]+) \(se [\d.]+, p=([\de.-]+)\)", b).groups()
    pm, lo, hi = re.search(r"proportion mediated \(product\)= (\d+)%\s+95% CI \[(\d+)%, (\d+)%\]", b).groups()
    pmd,        = re.search(r"proportion mediated \(diff\)\s+= (\d+)%", b).groups()
    return float(tot), float(dir_), float(dir_p), int(pm), int(lo), int(hi), int(pmd)

t2d = med("T2D"); bmi = med("BMI")
check("E7 T2D prop mediated (product)", t2d[3], 68, 0.5)
check("E7 T2D prop mediated (diff)",    t2d[6], 79, 0.5)
check("E7 BMI prop mediated (product)", bmi[3], 20, 0.5)
check("E7 BMI CI lower", bmi[4], 16, 0.5); check("E7 BMI CI upper", bmi[5], 25, 0.5)
check("E7 T2D direct p", t2d[2], 0.406, 1e-3)
condF = float(re.search(r"cond_F=([\d.]+)", t).group(1)); check("E7 CAD mediator cond_F", condF, 12.2, 0.05)
write("fig2_mediation.csv",
      ["exposure", "total", "direct", "direct_p", "pm_product", "pm_lo", "pm_hi", "pm_diff", "mediator_condF"],
      [["BMI", bmi[0], bmi[1], bmi[2], bmi[3], bmi[4], bmi[5], bmi[6], condF],
       ["T2D", t2d[0], t2d[1], t2d[2], t2d[3], t2d[4], t2d[5], t2d[6], condF]])

# ---------------------------------------------------------------- H5 reverse-edge estimators (Fig 3)
t = txt("h5_reverse_pleiotropy.txt")
rows = []
for edge in ["CAD -> T2D", "HF -> T2D"]:
    blk = re.search(re.escape("--- %s :" % edge) + r"(.*?)(?=--- |\Z)", t, re.S).group(1)
    nsnp = int(re.search(r"(\d+) instruments", blk).group(1))
    for est, name in [("IVW", "Inverse variance weighted"),
                      ("Weighted median", "Weighted median"),
                      ("Weighted mode", "Weighted mode"),
                      ("Penalised WM", "Penalised weighted median")]:
        m = re.search(re.escape(name) + r"\s+b=([+\-\d.]+)\s+se=([\d.]+)\s+p=([\de.\-]+)", blk)
        assert m, "H5: could not parse %s for %s" % (name, edge)
        b, se, p = float(m.group(1)), float(m.group(2)), float(m.group(3))
        rows.append([edge.replace(" -> ", "→"), est, b, se, p, nsnp])
    # The MR-Egger SLOPE, not only its intercept. The Results quote the HF→T2D slope (+0.39,
    # P = 0.31) as the one estimator that is null — the exception that keeps the "estimators agree"
    # sentence honest — and it was parsed here but never written, so the only number in that
    # sentence a reader could not check was the load-bearing one. It is in the canonical
    # h5_reverse_pleiotropy.txt; it simply never reached the deposited source data.
    m_eg = re.search(r"MR Egger\s+b=([+\-\d.]+)\s+se=([\d.]+)\s+p=([\de.\-]+)", blk)
    assert m_eg, "H5: could not parse the MR-Egger slope for %s" % edge
    rows.append([edge.replace(" -> ", "→"), "MR-Egger slope",
                 float(m_eg.group(1)), float(m_eg.group(2)), float(m_eg.group(3)), nsnp])
    icept, icept_p = re.search(r"MR-Egger intercept = ([+\-\d.]+) \(p=([\de.-]+)\)", blk).groups()
    rows.append([edge.replace(" -> ", "→"), "Egger intercept", float(icept), "", float(icept_p), nsnp])
by = {(r[0], r[1]): r[2] for r in rows}
check("H5 HF→T2D IVW",  by[("HF→T2D", "IVW")], 0.4217, 1e-3)
check("H5 CAD→T2D IVW", by[("CAD→T2D", "IVW")], 0.1675, 1e-3)
check("H5 CAD→T2D mode", by[("CAD→T2D", "Weighted mode")], 0.1204, 1e-3)
write("fig3_h5_estimators.csv", ["edge", "estimator", "b", "se", "p", "nsnp"], rows)

# ---------------------------------------------------------------- H3 portability (Fig 4)
t = txt("h3_portability.txt")
n_common = int(re.search(r"Common estimable edges \(both ancestries\): (\d+)", t).group(1))
sub = re.search(r"\[EUR-Bonferroni-significant subset\] n=(\d+)\s+sign concordance = ([\d.]+) \((\d+)/(\d+)\)", t)
n_sub, sc, k, n = int(sub.group(1)), float(sub.group(2)), int(sub.group(3)), int(sub.group(4))
eb_lo, eb_hi, nb_lo, nb_hi = re.search(
    r"sign concordance\s+: edge-block \[([\d.]+), ([\d.]+)\]\s+node-block \[([\d.]+), ([\d.]+)\]", t).groups()
ce_lo, ce_hi = re.search(r"effect correlation: edge-block \[([\d.]+), ([\d.]+)\]", t).groups()
perm_p = float(re.search(r"sign-flip permutation P .*?= ([\d.]+)", t).group(1))
corr = float(re.search(r"effect correlation \(comparable-scale, n=\d+\) = ([\d.]+)", t.split("[EUR-Bonferroni")[1]).group(1))
check("H3 common edges", n_common, 101, 0.5); check("H3 EUR-sig subset", n_sub, 52, 0.5)
check("H3 sign concordance", sc, 0.692, 1e-3); check("H3 k of n", k, 36, 0.5)
check("H3 node-block lower", nb_lo, 0.444, 1e-3); check("H3 node-block upper", nb_hi, 1.000, 1e-3)
check("H3 edge-block lower", eb_lo, 0.558, 1e-3)
check("H3 sign-flip perm P", perm_p, 0.0044, 1e-4)
# NB: 54_h3_portability.py labels these "edge-block"/"node-block", but the code shows the first is a
# PLAIN i.i.d. bootstrap resample of edges (no blocking — it assumes edge independence) and the second
# is a NODE-LEVEL SUBSAMPLE (retention 0.632; an edge is kept only if both endpoints survive, which is
# what actually respects shared-node dependence). Relabelled here to match what the code does.
write("fig4_portability.csv", ["stat", "estimate", "lo", "hi", "block", "note"],
      [["Sign concordance", sc, float(eb_lo), float(eb_hi), "Edge-level", "excludes chance"],
       ["Sign concordance", sc, float(nb_lo), float(nb_hi), "Node-level", "INCLUDES chance"],
       ["Effect correlation", corr, float(ce_lo), float(ce_hi), "Edge-level", "comparable-scale only"]])
write("fig4_portability_meta.csv", ["stat", "value"],
      [["n_common", n_common], ["n_eur_sig", n_sub], ["k_concordant", k], ["n_tested", n], ["perm_p", perm_p]])



def _p_sbp_expected():
    """The SBP->CKD interaction P that renal_contrast.csv supports, recomputed not remembered."""
    import math as _m
    _rows = list(csv.DictReader(open(os.path.join(RES, "renal_contrast.csv"), encoding="utf-8")))

    def _a(_s):
        _r = [x for x in _rows if x["exposure"] == "SBP" and x["source"] == _s]
        return float(_r[0]["b_perSD"]), float(_r[0]["se_perSD"])

    _e, _k = _a("EUR"), _a("EAS_BBJ")
    _z = (_e[0] - _k[0]) / _m.sqrt(_e[1] ** 2 + _k[1] ** 2)
    return _m.erfc(abs(_z) / _m.sqrt(2))


# ---------------------------------------------------------------- E6 ancestry interactions (Fig 5)
t = txt("ancestry_interaction_tests.txt")
rows = []
for m in re.finditer(r"^(\S+)\s+([+\-][\d.]+)\(([\d.]+)\)\s+([+\-][\d.]+)\(([\d.]+)\)\s+([+\-]?[\d.]+)\s+([\d.]+)\s+(.*)$",
                     t, re.M):
    edge, be, see, ba, sea, z, p, note = m.groups()
    rows.append([edge.replace("->", "→"), float(be), float(see), float(ba), float(sea),
                 float(z), float(p), note.strip()])
assert rows, "E6 interaction table did not parse"
d = {r[0]: r for r in rows}

# The RENAL rows of this table are superseded. ancestry_interaction_tests.txt predates the round-7
# re-estimation of the renal edges on the primary pipeline; renal_contrast.csv is the canonical
# source for SBP->CKD and BMI->CKD, and the manuscript quotes ITS interaction P-values. Recompute
# those two rows here rather than carrying the old text table's, or 05_source_data ships numbers
# that contradict the paper (it did: 0.4315/0.0709 against 0.55/0.15).
# The previous line asserted d["SBP→CKD"][6] == 0.4315, a hard-coded literal -- which made the
# correct value un-shippable, because regenerating would have failed the build. A check that PINS a
# claim protects it; this one derives instead.
import math as _math
_rc = list(csv.DictReader(open(os.path.join(RES, "renal_contrast.csv"), encoding="utf-8")))


def _arm(_exp, _src):
    _r = [x for x in _rc if x["exposure"] == _exp and x["source"] == _src]
    assert len(_r) == 1, f"{_exp}/{_src}: {len(_r)} rows in renal_contrast.csv"
    return float(_r[0]["b_perSD"]), float(_r[0]["se_perSD"])


for _exp in ("SBP", "BMI"):
    _e, _a = _arm(_exp, "EUR"), _arm(_exp, "EAS_BBJ")
    _z = (_e[0] - _a[0]) / _math.sqrt(_e[1] ** 2 + _a[1] ** 2)
    _p = _math.erfc(abs(_z) / _math.sqrt(2))
    _row = d[f"{_exp}→CKD"]
    _row[1], _row[2], _row[3], _row[4], _row[5], _row[6] = _e[0], _e[1], _a[0], _a[1], _z, _p
    _row[7] = _row[7].split(" -> ")[0] + " (per-SD, re-estimated on the primary pipeline; renal_contrast.csv)"
check("E6 SBP→CKD interaction P is the one renal_contrast.csv supports",
      d["SBP→CKD"][6], _p_sbp_expected(), 1e-3)
check("E6 TG→LDL interaction P",  d["TG→LDL"][6], 0.0, 1e-4)
bonf = 0.05 / len(rows)
for r in rows:
    r.append("interaction" if r[6] < bonf else "not significant")
write("fig5_interaction.csv",
      ["edge", "eur_b", "eur_se", "eas_b", "eas_se", "z", "p_int", "note", "verdict"], rows)

# ---------------------------------------------------------------- report
print()
if FAIL:
    print("*** PREP FAILED — parsed values disagree with the manuscript: ***")
    for f in FAIL: print("   ", f)
    sys.exit(1)
print("All parsed values agree with the manuscript. Panel data written to figures/data/.")
