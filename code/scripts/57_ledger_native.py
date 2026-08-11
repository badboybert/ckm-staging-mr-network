# -*- coding: utf-8 -*-
"""Step 57 (reviewer C4/E3): NATIVE-SCALE rebuild of the falsifiable AHA-staging ledger.

WHY THIS REPLACES THE POWER GATE IN 11_ledger.py
------------------------------------------------
The original ledger adjudicated each forward transition X->Y by testing the reverse edge Y->X
with (a) an 80%-power minimum-detectable-effect gate, `MDE_rev = (z_a + z_pwr) * SE_rev <= |b_fwd|`,
and (b) a two-one-sided-test equivalence check against a single fixed bound SESOI = 0.05
"standardised units". Both are unit-inconsistent:

  * the power gate compares a quantity in the REVERSE edge's units (MDE_rev, e.g. log-odds of CAD)
    against the FORWARD effect in the FORWARD edge's units (|b_fwd|, e.g. SD of BMI). For a
    continuous<->binary transition these are simply different scales, so the comparison has no
    defined meaning. 11_ledger.py acknowledged this with a `cross_type` flag but still applied
    the gate.
  * the TOST bound of 0.05 was applied to every reverse edge regardless of whether its outcome
    was measured per-SD, per-mmHg, or as a log-odds ratio.

This script adjudicates each reverse edge ON ITS OWN NATIVE SCALE: it reports the reverse estimate
with a 95% confidence interval in the units that edge was actually estimated in, and compares that
interval with a NEGLIGIBILITY BOUND defined for that outcome's scale. Verdicts become graded
evidence labels rather than a binary powered/underpowered call. The power gate and TOST are still
computed and carried as SECONDARY columns so the two formulations can be compared directly.

NEGLIGIBILITY BOUNDS (native units, pre-specified here, not tuned to the data)
  binary outcome (CAD/HF/Stroke/T2D/CKD), log-odds: log(1.05) = 0.0488  -> a 5% odds change
  SBP outcome, per mmHg                            : 1.0 mmHg
  other continuous outcomes, per SD                : 0.05 SD

GRADED LABELS (assigned to the REVERSE edge Y->X)
  reverse-excluded          : reverse 95% CI lies entirely inside +/- bound -> no practically
                              meaningful reverse effect is compatible with the data => CONCORDANT
  reverse-supported         : reverse CI excludes 0, |b_rev| > bound, Steiger-correct, clean Egger
                              intercept -> genuine reverse causation              => DISCORDANT
  reverse-pleiotropy-caveated: reverse CI excludes 0, |b_rev| > bound, Steiger-CORRECT, but the
                              Egger intercept is non-zero (P<0.05). The reverse signal is REAL;
                              only its magnitude is pleiotropy-inflated      => DISCORDANT_CAVEATED
  reverse-wrong-direction   : reverse CI excludes 0 but the edge is Steiger-WRONG-direction, i.e.
                              the variants explain more variance in the nominal outcome than in the
                              nominal exposure -> confounded/index-event, not reverse causation
                                                                                  => CONCORDANT
  reverse-inconclusive      : CI spans both 0 and the bound -> cannot separate null from a
                              meaningful reverse effect                           => INDETERMINATE

NB on the pleiotropy split: an earlier cut of this script collapsed the last three labels into a
single "reverse-artefactual" class. That was wrong and would have contradicted the project's frozen
calibration: CAD->T2D is Steiger-CORRECT with a non-zero Egger intercept, and the manuscript's
position (corroborated independently by the H5 pleiotropy-robust suite) is that it is genuine
feedback WITH a pleiotropy caveat, NOT an artifact. Only Steiger-wrong-direction edges are treated
as non-threatening to the forward order.

KNIFE-EDGE FLAG: a verdict that turns on a CI endpoint sitting within 10% of the negligibility bound
is flagged, because such a classification would flip under trivial re-rounding. These are reported,
never silently resolved.

Reads  results/forward_local_edges.csv (+ reverse_arms_edges.csv if present)
Writes results/staging_ledger_native.csv and results/staging_ledger_native.txt
Deterministic; no randomness.
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
Z_A   = 1.959963985   # two-sided alpha = 0.05
Z_PWR = 0.8416212     # 80% power
SESOI = 0.05          # legacy fixed bound, retained only for the secondary TOST column

STAGE = {"BMI":1, "SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"T2D":2,"eGFR":2,"CKD":2,
         "CAD":4,"HF":4,"Stroke":4}
BINARY = {"CAD","HF","Stroke","T2D","CKD"}

# Same pre-specified transition list as 11_ledger.py — the ledger's content is unchanged, only the
# adjudication rule is rebuilt.
PREREG = [
    ("BMI","T2D"), ("BMI","CAD"), ("BMI","HF"), ("BMI","Stroke"),
    ("SBP","CAD"), ("SBP","HF"), ("SBP","Stroke"),
    ("LDL","CAD"), ("LDL","HF"),
    ("HbA1c","T2D"), ("HbA1c","CAD"),
    ("T2D","CAD"), ("T2D","HF"), ("T2D","Stroke"),
    ("CAD","HF"),
]

def native_bound(outcome):
    """Negligibility bound for an effect ON `outcome`, in that outcome's own units."""
    if outcome in BINARY:
        return math.log(1.05)      # 0.0488 log-odds = a 5% odds change
    if outcome == "SBP":
        return 1.0                 # 1 mmHg
    return 0.05                    # 0.05 SD for per-SD continuous outcomes

def native_units(outcome):
    if outcome in BINARY: return "log-odds"
    if outcome == "SBP":  return "mmHg"
    return "per SD"

def load(path, into=None):
    d = into if into is not None else {}
    with io.open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                eip = r.get("egger_intercept_p", "")
                d.setdefault((r["exposure"], r["outcome"]), dict(
                    b=float(r["ivw_b"]), se=float(r["ivw_se"]), p=float(r["ivw_p"]),
                    n=int(float(r["nsnp"])),
                    steiger=str(r.get("steiger_correct", "")).upper() == "TRUE",
                    egger_int_p=(float(eip) if eip not in ("", "NA") else float("nan"))))
            except (ValueError, KeyError):
                continue
    return d

def norm_cdf(x): return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def tost(b, se, sesoi):
    if se <= 0: return None
    return max(1 - norm_cdf((b + sesoi) / se), norm_cdf((b - sesoi) / se))

def knife_edge(lo, hi, bound, frac=0.10):
    """True when a CI endpoint sits within `frac` of the negligibility bound, so the verdict would
    flip under trivial re-rounding. Reported, never silently resolved."""
    if lo is None: return False
    return min(abs(bound - hi), abs(bound + lo)) < frac * bound

def adjudicate(bf, pf, rev, X, Y):
    """Return (label, verdict) for transition X->Y given the reverse edge Y->X."""
    if pf >= 0.05:
        return "forward-not-established", "FWD_NOT_SIG"
    if rev is None:
        return "reverse-not-estimable", "NO_REVERSE_DATA"
    br, sr = rev["b"], rev["se"]
    lo, hi = br - Z_A * sr, br + Z_A * sr
    bound = native_bound(X)                      # reverse edge Y->X has OUTCOME X
    ci_excludes_zero = (lo > 0) or (hi < 0)
    ci_within_bound  = (lo > -bound) and (hi < bound)
    steiger_ok = rev["steiger"]
    eip = rev["egger_int_p"]
    # 2026-07-26 (round-4 Tier 1). This line used to read `(eip == eip) and (eip < 0.05)`, i.e. a
    # MISSING intercept was silently treated as "no directional pleiotropy". That is the wrong
    # default for a classifier: absence of a statistic is not evidence of its null. It bit exactly
    # once and it mattered — step 73 blanked the Egger columns for the X->SBP rows when the strict
    # extraction became primary, so CAD->SBP arrived here with eip = NaN, `NaN < 0.05` was False,
    # and the transition was recorded as reverse-SUPPORTED instead of reverse-PLEIOTROPY-CAVEATED.
    # The ledger then said 3 supported + 1 caveated while the Results and Figure 1 legend said 2 + 2,
    # which is the cross-file contradiction the round-4 review caught. The intercept is now computed
    # on the primary strict set (step 75: CAD->SBP intercept +0.062, P = 0.033).
    # A missing intercept is now a HARD ERROR rather than a silent downgrade.
    eip = rev["egger_int_p"]
    if eip != eip:                                      # NaN
        raise ValueError(
            f"reverse edge {rev.get('exposure','?')}->{rev.get('outcome','?')} has no Egger "
            f"intercept P. Refusing to classify: a missing pleiotropy statistic must not default to "
            f"'clean'. Compute it (see scripts/75_sbp_strict_pleiotropy.R) or drop the edge.")
    dir_pleiotropy = eip < 0.05
    if ci_excludes_zero and abs(br) > bound:
        if not steiger_ok:
            # Variants explain more variance in the nominal outcome than the nominal exposure:
            # a confounded / index-event signal, not reverse causation.
            return "reverse-wrong-direction", "CONCORDANT"
        if dir_pleiotropy:
            # Steiger-CORRECT but pleiotropy-inflated: the reverse effect is real, its magnitude is not.
            return "reverse-pleiotropy-caveated", "DISCORDANT_CAVEATED"
        return "reverse-supported", "DISCORDANT"
    if ci_within_bound:
        return "reverse-excluded", "CONCORDANT"
    return "reverse-inconclusive", "INDETERMINATE"

def main():
    d = load(os.path.join(BASE, "results/forward_local_edges.csv"))
    rev_path = os.path.join(BASE, "results/reverse_arms_edges.csv")
    if os.path.exists(rev_path):
        load(rev_path, into=d)

    rows = []
    for X, Y in PREREG:
        f = d.get((X, Y)); rev = d.get((Y, X))
        if f is None:
            continue
        bf, sef, pf, nf = f["b"], f["se"], f["p"], f["n"]
        label, verdict = adjudicate(bf, pf, rev, X, Y)
        if rev is None:
            br = sr = pr = nr = lo = hi = mde = tp = None
            rev_st = rev_eip = None
        else:
            br, sr, pr, nr = rev["b"], rev["se"], rev["p"], rev["n"]
            lo, hi = br - Z_A * sr, br + Z_A * sr
            mde = (Z_A + Z_PWR) * sr                 # secondary: legacy power gate
            tp = tost(br, sr, SESOI)                 # secondary: legacy TOST
            rev_st, rev_eip = rev["steiger"], rev["egger_int_p"]
        ke = knife_edge(lo, hi, native_bound(X)) if lo is not None else False
        rows.append(dict(
            transition=f"{X}->{Y}", X=X, Y=Y, stage_X=STAGE.get(X), stage_Y=STAGE.get(Y),
            b_fwd=bf, se_fwd=sef, p_fwd=pf, n_fwd=nf,
            rev_edge=f"{Y}->{X}", b_rev=br, se_rev=sr, p_rev=pr, n_rev=nr,
            rev_ci_lo=lo, rev_ci_hi=hi,
            rev_units=native_units(X), rev_bound=native_bound(X),
            rev_steiger_correct=rev_st, rev_egger_int_p=rev_eip,
            evidence_label=label, verdict=verdict, knife_edge=ke,
            legacy_mde=mde, legacy_tost_p=tp))

    out_csv = os.path.join(BASE, "results/staging_ledger_native.csv")
    cols = ["transition","X","Y","stage_X","stage_Y","b_fwd","se_fwd","p_fwd","n_fwd",
            "rev_edge","b_rev","se_rev","p_rev","n_rev","rev_ci_lo","rev_ci_hi",
            "rev_units","rev_bound","rev_steiger_correct","rev_egger_int_p",
            "evidence_label","verdict","knife_edge","legacy_mde","legacy_tost_p"]
    with io.open(out_csv, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow(r)

    # ---- summary ----
    from collections import Counter
    vc = Counter(r["verdict"] for r in rows)
    lc = Counter(r["evidence_label"] for r in rows)
    L = []
    L.append("=== NATIVE-SCALE falsifiable staging ledger (reviewer C4/E3) ===")
    L.append(f"Transitions adjudicated: {len(rows)}")
    L.append("")
    L.append("Each reverse edge Y->X is judged on ITS OWN native scale: the reverse 95% CI is compared")
    L.append("with a negligibility bound for that outcome's units (binary log-odds: log(1.05)=0.0488;")
    L.append("SBP: 1 mmHg; other continuous: 0.05 SD). This replaces the cross-unit power gate")
    L.append("(MDE_rev vs |b_fwd|) and the single fixed TOST bound, both of which compared quantities")
    L.append("measured on different scales. Legacy power/TOST values are retained as secondary columns.")
    L.append("")
    L.append("VERDICTS: " + ", ".join(f"{k}={v}" for k, v in sorted(vc.items())))
    L.append("LABELS  : " + ", ".join(f"{k}={v}" for k, v in sorted(lc.items())))
    L.append("")
    hdr = f"{'transition':<14}{'b_fwd':>9}{'reverse':<14}{'b_rev':>9}{'95% CI':>22}{'bound':>8}  {'label':<22}{'verdict'}"
    L.append(hdr); L.append("-" * len(hdr))
    for r in rows:
        ci = "        n/a" if r["b_rev"] is None else f"[{r['rev_ci_lo']:+.3f}, {r['rev_ci_hi']:+.3f}]"
        br = "      n/a" if r["b_rev"] is None else f"{r['b_rev']:+.4f}"
        L.append(f"{r['transition']:<14}{r['b_fwd']:+9.4f}{r['rev_edge']:<14}{br:>9}{ci:>22}"
                 f"{r['rev_bound']:>8.4f}  {r['evidence_label']:<22}{r['verdict']}")
    L.append("")
    conc = vc.get("CONCORDANT", 0); disc = vc.get("DISCORDANT", 0)
    dcav = vc.get("DISCORDANT_CAVEATED", 0); ind = vc.get("INDETERMINATE", 0)
    L.append(f"READ: {conc} concordant, {disc} discordant (genuine reverse effect), "
             f"{dcav} discordant-with-pleiotropy-caveat, {ind} indeterminate.")
    L.append("DISCORDANT = the reverse edge is significant on its native scale, exceeds the negligibility")
    L.append("bound, is Steiger-correct and carries no directional pleiotropy. Both such edges here point")
    L.append("INTO diabetes (HF->T2D, Stroke->T2D): they are reciprocal feedback into a stage-2 trait, not")
    L.append("reverse stage progression, and the Stroke->T2D arm rests on only 18 reverse instruments.")
    L.append("DISCORDANT_CAVEATED = Steiger-correct and significant, but the Egger intercept is non-zero, so")
    L.append("the reverse effect is real while its magnitude is pleiotropy-inflated (CAD->SBP index-event;")
    L.append("CAD->T2D feedback-with-caveat, matching the independent H5 pleiotropy-robust suite).")
    ke_rows = [r for r in rows if r["knife_edge"]]
    if ke_rows:
        L.append("")
        L.append("KNIFE-EDGE VERDICTS (a CI endpoint within 10% of the negligibility bound; these would flip")
        L.append("under trivial re-rounding and are reported as such, not resolved):")
        for r in ke_rows:
            L.append(f"  {r['transition']:<12} reverse {r['rev_edge']}: CI [{r['rev_ci_lo']:+.4f}, "
                     f"{r['rev_ci_hi']:+.4f}] vs bound {r['rev_bound']:.4f} -> {r['verdict']}")
    L.append("")
    L.append("COMPARISON WITH THE LEGACY POWER-GATE LEDGER (11_ledger.py, results/staging_ledger.csv):")
    L.append("the legacy rule returned 11 concordant / 2 discordant / 2 discordant-pleiotropic. The native-")
    L.append("scale rule is more conservative where the reverse arm is thin -- BMI->HF, BMI->Stroke and")
    L.append("SBP->Stroke move to INDETERMINATE because their reverse CIs span both zero and the")
    L.append("negligibility bound -- and it separates pleiotropy-caveated reverse effects from")
    L.append("Steiger-wrong-direction ones, which the legacy rule pooled.")
    txt = "\n".join(L) + "\n"
    out_txt = os.path.join(BASE, "results/staging_ledger_native.txt")
    io.open(out_txt, "w", encoding="utf-8").write(txt)
    print(txt)
    print("wrote", out_csv)
    print("wrote", out_txt)

if __name__ == "__main__":
    main()
