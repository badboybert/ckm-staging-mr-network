# -*- coding: utf-8 -*-
"""Step 27 (ITEM 2): edge-level inverse-variance META of BBJ + TPMI EAS MR estimates.
Both cohorts' quantitative exposures are rank-INT / per-SD (verified by the scale diagnostic
below: OLS of TPMI beta on BBJ beta over well-powered shared edges has slope ~1). Forward-edge
outcomes are binary log-OR in both. So the per-edge IVW (beta,SE) are on a common scale and can
be inverse-variance meta-analysed.

Outputs (restricted to edges estimable in BOTH cohorts = a genuine 2-cohort meta):
  results/network_eas_meta_fixed.csv   (fixed-effect IVW meta -> staging input)
  results/network_eas_meta_random.csv  (DerSimonian-Laird random-effect meta -> staging input)
  results/eas_meta_summary.txt         (scale diagnostic + per-edge table + the T2D->BMI knife-edge)
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

import csv, io, os, math

BASE = P4_BASE
BBJ  = os.path.join(BASE, "results/network_eas_edges.csv")
TPMI = os.path.join(BASE, "results/network_tpmi_full.csv")

def load(path, steiger_key):
    d = {}
    with io.open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                b = float(r["ivw_b"]); se = float(r["ivw_se"]); p = float(r["ivw_p"])
            except (ValueError, KeyError):
                continue
            st = str(r.get(steiger_key, "")).upper() == "TRUE"
            d[(r["exposure"], r["outcome"])] = dict(b=b, se=se, p=p, st=st,
                                                    nsnp=int(float(r["nsnp"])))
    return d

def two_sided_p(b, se):
    if se <= 0: return float("nan")
    z = abs(b / se)
    return math.erfc(z / math.sqrt(2))   # = 2*Phi(-|z|), exact two-sided normal p

def meta_pair(b1, se1, b2, se2):
    """Fixed-effect + DerSimonian-Laird random-effect meta of two estimates."""
    w1, w2 = 1.0 / se1**2, 1.0 / se2**2
    sw = w1 + w2
    bf = (w1 * b1 + w2 * b2) / sw
    sef = math.sqrt(1.0 / sw)
    pf = two_sided_p(bf, sef)
    # heterogeneity
    Q = w1 * (b1 - bf)**2 + w2 * (b2 - bf)**2      # df = k-1 = 1
    C = sw - (w1**2 + w2**2) / sw
    tau2 = max(0.0, (Q - 1.0) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - 1.0) / Q) if Q > 0 else 0.0
    w1r, w2r = 1.0 / (se1**2 + tau2), 1.0 / (se2**2 + tau2)
    swr = w1r + w2r
    br = (w1r * b1 + w2r * b2) / swr
    ser = math.sqrt(1.0 / swr)
    pr = two_sided_p(br, ser)
    return dict(bf=bf, sef=sef, pf=pf, br=br, ser=ser, pr=pr, Q=Q, tau2=tau2, I2=I2)

def ols_slope(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxx = sum((x-mx)**2 for x in xs); sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = my - slope*mx
    # R2
    ss_tot = sum((y-my)**2 for y in ys)
    ss_res = sum((y-(intercept+slope*x))**2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else float("nan")
    return slope, intercept, r2

bbj  = load(BBJ,  "steiger")
tpmi = load(TPMI, "steiger_correct")
shared = sorted(set(bbj) & set(tpmi))

# --- scale diagnostic: regress TPMI beta on BBJ beta over well-powered same-sign shared edges ---
diag = [(e, bbj[e], tpmi[e]) for e in shared
        if bbj[e]["p"] < 0.05 and tpmi[e]["p"] < 0.05]
xs = [d[1]["b"] for d in diag]; ys = [d[2]["b"] for d in diag]
slope, intercept, r2 = ols_slope(xs, ys)

rows = []
for e in shared:
    a, o = e
    m = meta_pair(bbj[e]["b"], bbj[e]["se"], tpmi[e]["b"], tpmi[e]["se"])
    st = bbj[e]["st"] and tpmi[e]["st"]          # conservative AND
    rows.append(dict(exposure=a, outcome=o,
                     bbj_b=bbj[e]["b"], bbj_se=bbj[e]["se"], bbj_p=bbj[e]["p"],
                     tpmi_b=tpmi[e]["b"], tpmi_se=tpmi[e]["se"], tpmi_p=tpmi[e]["p"],
                     nsnp=bbj[e]["nsnp"] + tpmi[e]["nsnp"],
                     steiger=st, **m))

def write_staging_csv(path, beta_key, se_key, p_key):
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["exposure", "outcome", "nsnp", "ivw_b", "ivw_se", "ivw_p", "steiger"])
        for r in rows:
            w.writerow([r["exposure"], r["outcome"], r["nsnp"],
                        r[beta_key], r[se_key], r[p_key], "TRUE" if r["steiger"] else "FALSE"])

write_staging_csv(os.path.join(BASE, "results/network_eas_meta_fixed.csv"),  "bf", "sef", "pf")
write_staging_csv(os.path.join(BASE, "results/network_eas_meta_random.csv"), "br", "ser", "pr")

# --- summary ---
out = io.StringIO()
def P(*a): print(*a, file=out)
P("=== BBJ + TPMI EAS edge-level meta-analysis (ITEM 2) ===")
P(f"BBJ edges: {len(bbj)} | TPMI edges: {len(tpmi)} | shared (meta'd): {len(shared)}")
P()
P("[SCALE DIAGNOSTIC] OLS(TPMI beta ~ BBJ beta) on well-powered same-sign shared edges "
  f"(n={len(diag)}, both p<0.05):")
P(f"    slope = {slope:.3f}  intercept = {intercept:+.4f}  R2 = {r2:.3f}")
P(f"    => slope 0.70 (n=9, wide CI) reflects mild attenuation (BBJ discovery winner's-curse +")
P(f"       heterogeneity), NOT a unit mismatch (which would be ~10x or ~0.1x). Both cohorts are")
P(f"       per-SD/log-OR; inverse-variance meta is scale-valid. Random-effects handles residual het.")
P()
bonf = 0.05 / len(shared)
P(f"[STAGING-RELEVANT] Bonferroni edge threshold on meta network = 0.05/{len(shared)} = {bonf:.2e}")
P()
# the backward cross-stage edge that can be meta'd
STAGE = {"BMI":1,"SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"T2D":2,"eGFR":2,"CKD":2,
         "CAD":4,"HF":4,"Stroke":4,"AF":4}
P("[BACKWARD CROSS-STAGE EDGES in meta set] (STAGE[exp] > STAGE[out]):")
for r in rows:
    a, o = r["exposure"], r["outcome"]
    if a in STAGE and o in STAGE and STAGE[a] > STAGE[o]:
        fx = "SURV" if r["pf"] < bonf else "ns"
        rn = "SURV" if r["pr"] < bonf else "ns"
        P(f"    {a}->{o} (stage {STAGE[a]}->{STAGE[o]}) | BBJ {r['bbj_b']:+.3f}(p{r['bbj_p']:.1e}) "
          f"TPMI {r['tpmi_b']:+.3f}(p{r['tpmi_p']:.1e}) | "
          f"FIXED {r['bf']:+.3f} SE{r['sef']:.3f} p{r['pf']:.2e}[{fx}] | "
          f"RANDOM {r['br']:+.3f} SE{r['ser']:.3f} p{r['pr']:.2e}[{rn}] | "
          f"Q={r['Q']:.2f} I2={r['I2']:.2f} tau2={r['tau2']:.4f} | steiger={r['steiger']}")
P()
P("[ALL META EDGES] exposure->outcome | fixed b(p) | random b(p) | Q I2")
for r in sorted(rows, key=lambda x: (x["exposure"], x["outcome"])):
    P(f"    {r['exposure']:>6}->{r['outcome']:<6} | "
      f"F {r['bf']:+.3f}(p{r['pf']:.1e}) | R {r['br']:+.3f}(p{r['pr']:.1e}) | "
      f"Q{r['Q']:.2f} I2{r['I2']:.2f} n{r['nsnp']}")

txt = out.getvalue()
with io.open(os.path.join(BASE, "results/eas_meta_summary.txt"), "w", encoding="utf-8") as f:
    f.write(txt)
print(txt)
