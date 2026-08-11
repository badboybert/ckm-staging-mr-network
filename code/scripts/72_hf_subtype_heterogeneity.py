# -*- coding: utf-8 -*-
"""Formal HFpEF-versus-HFrEF between-subtype test (round-3 blocker 5).

The manuscript said BMI reaches both subtypes "comparably" and T2D reaches HFrEF but not HFpEF. Both
are statements about a DIFFERENCE, and neither was tested: they were read off two significance
verdicts. A significant estimate beside a non-significant one is not a difference.

Two complications make a naive z-test wrong, and both are handled rather than assumed away:

1. SHARED CONTROLS. Enzan 2025 reports HFpEF (26,743 cases) and HFrEF (23,749 cases) against the
   SAME 456,520 controls, so the two estimates are positively correlated. The Lin-Sullivan overlap
   formula gives the induced correlation from the shared controls; cases are disjoint by definition.

2. SHARED INSTRUMENTS. Both ratio estimates divide by the same SNP-exposure coefficients, which
   correlates them further by an amount the summary data cannot identify.

Positive correlation REDUCES Var(b1 - b2), so assuming independence INFLATES the P-value: the
independence test is conservative for detecting a difference. That asymmetry matters for how each
claim may be worded, and it is reported per exposure rather than stated once:
  * a difference that is significant under independence is significant a fortiori;
  * a NON-significant difference under independence does NOT support "comparable", because the
    better-powered correlated test could still reject.

Out: results/hf_subtype_heterogeneity.{csv,txt}
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

import io, os, csv, math, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE
RES = os.path.join(BASE, "results")
SRC = os.path.join(RES, "network_hfsubtypes_allcause.csv")

# Enzan 2025 (PMID 41184235), read from hf_allcause_subtypes.txt provenance block.
HFPEF_CASES, HFPEF_N = 26_743, 483_263
HFREF_CASES, HFREF_N = 23_749, 480_269
HFPEF_CTRL = HFPEF_N - HFPEF_CASES
HFREF_CTRL = HFREF_N - HFREF_CASES
assert HFPEF_CTRL == HFREF_CTRL, "control counts differ; the shared-control formula below assumes one control set"
SHARED_CTRL = HFPEF_CTRL
SHARED_CASES = 0                      # HFpEF and HFrEF cases are disjoint by definition


def norm_sf(z):
    """Two-sided normal tail, no SciPy dependency."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def lin_sullivan_rho(n1a, n0a, n1b, n0b, ns1, ns0):
    """Correlation between two case-control effect estimates sharing samples (Lin & Sullivan 2009)."""
    num = ns0 * math.sqrt(n1a * n1b / (n0a * n0b)) + ns1 * math.sqrt(n0a * n0b / (n1a * n1b))
    return num / math.sqrt((n1a + n0a) * (n1b + n0b))


RHO_SHARED = lin_sullivan_rho(HFPEF_CASES, HFPEF_CTRL, HFREF_CASES, HFREF_CTRL,
                              SHARED_CASES, SHARED_CTRL)

rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
by = {}
for r in rows:
    by.setdefault(r["exposure"], {})[r["outcome"]] = r

RHO_GRID = [0.0, RHO_SHARED, 0.10, 0.25, 0.50]
out_rows, lines = [], []

lines += [
    "=== Formal HFpEF vs HFrEF between-subtype heterogeneity test ===",
    "",
    "H0: beta_HFpEF = beta_HFrEF for the same exposure.",
    "z = (b1 - b2) / sqrt(se1^2 + se2^2 - 2*rho*se1*se2);  Cochran Q on 1 df at rho = 0.",
    "",
    f"Enzan 2025 shares ONE control set between the two outcomes:",
    f"  HFpEF {HFPEF_CASES:,} cases / {HFPEF_CTRL:,} controls",
    f"  HFrEF {HFREF_CASES:,} cases / {HFREF_CTRL:,} controls  (same controls; cases disjoint)",
    f"  Lin-Sullivan induced correlation from shared controls: rho = {RHO_SHARED:.4f}",
    "",
    "Positive rho SHRINKS Var(b1-b2), so the rho = 0 column is the CONSERVATIVE test for a difference.",
    "A difference significant at rho = 0 is significant a fortiori; a NON-significant one at rho = 0",
    "does not establish equivalence, because the better-powered correlated test could still reject.",
    "",
]

hdr = f"{'exposure':10s} {'b_HFpEF':>10s} {'b_HFrEF':>10s} {'diff':>9s} " + \
      " ".join(f"{'P(rho='+format(r, '.2f')+')':>14s}" for r in RHO_GRID) + f" {'Q_p(rho=0)':>11s}  verdict"
lines.append(hdr)
lines.append("-" * len(hdr))

for exp in sorted(by):
    if "HFpEF" not in by[exp] or "HFrEF" not in by[exp]:
        continue
    a, b = by[exp]["HFpEF"], by[exp]["HFrEF"]
    b1, s1 = float(a["ivw_b"]), float(a["ivw_se"])
    b2, s2 = float(b["ivw_b"]), float(b["ivw_se"])
    diff = b1 - b2
    ps = []
    for rho in RHO_GRID:
        var = s1 * s1 + s2 * s2 - 2.0 * rho * s1 * s2
        ps.append(norm_sf(diff / math.sqrt(var)) if var > 0 else float("nan"))
    # Cochran Q on the two estimates, independence
    w1, w2 = 1.0 / (s1 * s1), 1.0 / (s2 * s2)
    bbar = (w1 * b1 + w2 * b2) / (w1 + w2)
    Q = w1 * (b1 - bbar) ** 2 + w2 * (b2 - bbar) ** 2
    q_p = norm_sf(math.sqrt(Q))            # chi2 on 1 df == squared standard normal
    p0 = ps[0]
    verdict = ("DIFFERENT (significant at the conservative rho = 0)" if p0 < 0.05
               else "no difference detected — NOT equivalence")
    lines.append(f"{exp:10s} {b1:+10.4f} {b2:+10.4f} {diff:+9.4f} " +
                 " ".join(f"{p:14.3g}" for p in ps) + f" {q_p:11.3g}  {verdict}")

# machine-readable rows for the workbook and the figure, built in one place
out_rows = []
for exp in sorted(by):
    if "HFpEF" not in by[exp] or "HFrEF" not in by[exp]:
        continue
    a, b = by[exp]["HFpEF"], by[exp]["HFrEF"]
    b1, s1 = float(a["ivw_b"]), float(a["ivw_se"])
    b2, s2 = float(b["ivw_b"]), float(b["ivw_se"])
    diff = b1 - b2
    rec = dict(exposure=exp, b_hfpef=round(b1, 6), se_hfpef=round(s1, 6),
               b_hfref=round(b2, 6), se_hfref=round(s2, 6), diff=round(diff, 6))
    for rho in RHO_GRID:
        var = s1 * s1 + s2 * s2 - 2.0 * rho * s1 * s2
        rec[f"p_rho_{rho:.2f}"] = norm_sf(diff / math.sqrt(var))
    w1, w2 = 1.0 / (s1 * s1), 1.0 / (s2 * s2)
    bbar = (w1 * b1 + w2 * b2) / (w1 + w2)
    rec["Q"] = w1 * (b1 - bbar) ** 2 + w2 * (b2 - bbar) ** 2
    rec["Q_p"] = norm_sf(math.sqrt(rec["Q"]))
    rec["significant_at_rho0"] = rec["p_rho_0.00"] < 0.05
    out_rows.append(rec)

BONF = 0.05 / len(out_rows)
lines += [
    "",
    f"Within-family Bonferroni over {len(out_rows)} exposures: {BONF:.4g}",
    "",
    "READ:",
]
for r in out_rows:
    if r["p_rho_0.00"] < BONF:
        lines.append(f"  {r['exposure']}: subtype difference survives Bonferroni "
                     f"(P = {r['p_rho_0.00']:.3g}); the larger estimate is "
                     f"{'HFrEF' if r['diff'] < 0 else 'HFpEF'}.")
for r in out_rows:
    if r["p_rho_0.00"] >= 0.05:
        lines.append(f"  {r['exposure']}: NO difference detected (P = {r['p_rho_0.00']:.2f}). This does not "
                     f"license 'comparable' or 'equivalent' — it is an absence of evidence, and the "
                     f"correlated test is better powered.")

with open(os.path.join(RES, "hf_subtype_heterogeneity.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)
open(os.path.join(RES, "hf_subtype_heterogeneity.txt"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nwrote results/hf_subtype_heterogeneity.{{csv,txt}}")
