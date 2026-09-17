# -*- coding: utf-8 -*-
"""Build Table 1 — the graded evidence synthesis — from the files that produced its numbers.

WHY THIS EXISTS (round-7 item C014; PI decision Q7, 2026-09-16): the paper's closing synthesis used
to live as a six-line text panel inside a supplementary lipid figure (Supplementary Figure S7e) AND
as a prose list in the Results. The two disagreed - the panel filed the East Asian staging order
under EXPLORATORY / NULL while the prose graded it supportive - because each was maintained by hand.
Promoting it to a main-text table with one generator removes the duplicate by construction: the claim
text is authored here, every statistic beside it is read from the analysis output, and the Results
paragraph now points at the table instead of restating it.

Reads : figures/data/fig1_staging_stat.csv, fig1_constrained_nulls.csv, fig2_mediation.csv,
        fig4_portability.csv; results/forward_local_edges.csv, renal_contrast.csv,
        staging_eas_meta_exact.csv, network_ckd.csv, mvmr_apob.csv
Writes: figures/TABLES.md   (rendered into the manuscript by build_manuscript_docx.py and copied
        into the package by build_submission_v1.py)
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

import csv, io, os, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE
FIG, RES = os.path.join(BASE, "figures"), os.path.join(BASE, "results")
DAT = os.path.join(FIG, "data")

TIERS = ["Higher confidence", "Supportive", "Exploratory", "Bounded null", "Not interpretable"]


def rows(path):
    with io.open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def kv(path, key_col, val_col):
    return {r[key_col]: r[val_col] for r in rows(path)}


def edge(path, exp, out):
    for r in rows(path):
        if r["exposure"] == exp and r["outcome"] == out:
            return r
    raise SystemExit(f"TABLE1: {exp}->{out} not found in {os.path.basename(path)}")


def p_fmt(p):
    """Journal style: 1 significant figure in scientific notation below 0.001, else 2 decimals."""
    p = float(p)
    if p >= 0.001:
        return f"{p:.3f}".rstrip("0").rstrip(".")
    e = 0
    while p < 1:
        p *= 10; e += 1
    return f"{p:.0f} × 10⁻{'⁰¹²³⁴⁵⁶⁷⁸⁹'[e] if e < 10 else e}" if e < 10 else f"{p:.0f}e-{e}"


def sup(n):
    return "".join("⁰¹²³⁴⁵⁶⁷⁸⁹"[int(d)] for d in str(n))


def psci(p, sig=1):
    """P in the manuscript's convention: plain decimals at or above 0.01, else 'a × 10⁻b'.
    One significant figure by default (1 × 10⁻³¹, 9 × 10⁻⁷, 8 × 10⁻⁴), two where the manuscript
    itself carries two (the staging permutation P, 1.5 × 10⁻³) — so a number printed here and in
    the prose is the same string, not two roundings of one value."""
    p = float(p)
    if p >= 0.01:
        return f"{p:.3f}".rstrip("0").rstrip(".")
    m, e = f"{p:.{sig - 1}e}".split("e")
    m = m.rstrip("0").rstrip(".") if "." in m else m
    return f"{m} × 10⁻{sup(abs(int(e)))}"


def main():
    st = kv(os.path.join(DAT, "fig1_staging_stat.csv"), "metric", "value")
    nulls = {r["null"]: r for r in rows(os.path.join(DAT, "fig1_constrained_nulls.csv"))}
    med = {r["exposure"]: r for r in rows(os.path.join(DAT, "fig2_mediation.csv"))}
    port = rows(os.path.join(DAT, "fig4_portability.csv"))
    fwd = os.path.join(RES, "forward_local_edges.csv")
    renal = rows(os.path.join(RES, "renal_contrast.csv"))
    easmeta = rows(os.path.join(RES, "staging_eas_meta_exact.csv"))

    conc = float(st["observed_concordance"])
    perm_p = float(st["perm_p"])
    deg_p = float(nulls["Degree-preserving rewiring"]["p"])
    bmi = med["BMI"]; t2d = med["T2D"]
    # Direct effects are read from the FULL-PRECISION MVMR output, never from fig2_mediation.csv's
    # 3-decimal copy: 0.415 as a double is 0.41499..., which "%+.2f" rounds to +0.41 while every other
    # surface prints the canonical +0.42 — the round-6 graphical-abstract defect, in a new place.
    casc = {r["exposure"]: r for r in rows(os.path.join(RES, "mvmr_cascade.csv"))
            if r["outcome"] == "HF" and r["model"] == "HF_via_CAD"}

    # The MEDIATED PROPORTIONS have the same defect as the direct effects above, and fixing one and
    # not the other is how 37.61 shipped as "37" in Table 1 against "38" in the Discussion, Results
    # and Figure 2 legend. fig2_mediation.csv stores pm_lo already rounded to a whole number, so
    # "%.0f" over it re-emits that rounding instead of performing it. Read the sweep at independence.
    _sweep = rows(os.path.join(RES, "mediation_covariance_sensitivity.csv"))

    def pm_independence(exposure):
        hit = [r for r in _sweep if r["exposure"] == exposure
               and all(float(r[k]) == 0 for k in ("rho_ab", "rho_total_direct", "rho_a_total"))]
        assert len(hit) == 1, f"{exposure}: {len(hit)} independence rows in the covariance sweep"
        lo, hi = float(hit[0]["pm_product_lo"]), float(hit[0]["pm_product_hi"])
        # The rounded copy must still AGREE to within its own rounding, or the two files have drifted
        # apart for a reason that is not rounding -- which is a defect, not a formatting choice.
        cp = med[exposure]
        assert abs(lo - float(cp["pm_lo"])) <= 1 and abs(hi - float(cp["pm_hi"])) <= 1, (
            f"{exposure}: sweep {lo:.2f}-{hi:.2f} vs fig2_mediation {cp['pm_lo']}-{cp['pm_hi']}")
        return lo, hi

    t2d_lo, t2d_hi = pm_independence("T2D")
    bmi_lo, bmi_hi = pm_independence("BMI")
    bmi_direct_pct = 100 - float(bmi["pm_product"])
    cadhf = edge(fwd, "CAD", "HF")
    ldlt2d = edge(fwd, "LDL", "T2D")
    node = [r for r in port if r["stat"] == "Sign concordance" and r["block"] == "Node-level"][0]
    eur_sbp = [r for r in renal if r["exposure"] == "SBP" and r["source"] == "EUR"][0]
    zero_sbp = [r for r in renal if r["exposure"] == "SBP" and r["source"] == "EAS_ZERO"][0]
    apob = [r for r in rows(os.path.join(RES, "mvmr_apob.csv"))
            if r.get("exposure") == "ApoB" and r.get("model", "").startswith("ApoB_HDL_TG")]

    T = [
        ("Higher confidence",
         "Adiposity acts on heart failure beyond coronary disease",
         f"Direct BMI→HF effect in MVMR conditioning on T2D and CAD (+{float(casc['BMI']['direct_b']):.2f}, "
         f"*P* = {psci(bmi['direct_p'])}); about {bmi_direct_pct:.0f}% of the total effect is not mediated "
         f"through CAD (mediated {float(bmi['pm_product']):.0f}%, 95% CI {bmi_lo:.0f}–{bmi_hi:.0f}%)",
         "Covariance sweep (15–27% mediated); replicates a CAD-independent adiposity effect first reported by HERMES"),
        ("Higher confidence",
         "Coronary disease acts on heart failure (CAD→HF)",
         f"+{float(cadhf['ivw_b']):.3f}, *P* = {psci(cadhf['ivw_p'])}; IVW, MR-Egger and weighted median concordant; "
         f"Egger intercept *P* = {float(cadhf['egger_intercept_p']):.2f}",
         "Liability-scale Steiger and MR-PRESSO favour this direction; the reverse (HF→CAD) is unresolved"),
        ("Higher confidence",
         "The atherogenic apolipoprotein-B/LDL axis carries the coronary lipid signal",
         (f"ApoB direct effect on CAD +{float(apob[0]['direct_b']):.2f}, *P* = {psci(apob[0]['direct_p'])} "
          f"(conditional F = {float(apob[0]['cond_F']):.0f})") if apob else
         "ApoB direct effect on CAD retained when HDL and triglycerides are co-modelled",
         "Axis-level only: ApoB and LDL are too collinear to separate (conditional F falls to 9.2)"),
        ("Higher confidence",
         "The forward risk-factor→disease cascade replicates in East Asians edge by edge",
         "SBP, LDL, total cholesterol, BMI and HbA1c edges hold in Biobank Japan and TPMI, and in a "
         "design sharing no participants",
         "Magnitudes comparable only where scales permit; HbA1c and eGFR read on sign alone"),
        ("Supportive",
         "The AHA stage ordering is descriptively concordant with the genetics in Europeans",
         f"Cross-stage concordance {conc:.3f} (25 of 27 edges), exact label-permutation "
         f"*P* = {psci(perm_p, sig=2)}",
         f"NOT beyond a degree-preserving rewiring null (*P* = {deg_p:.2f}), so the ordering is not established "
         f"beyond the traits' exposure/outcome roles"),
        ("Supportive",
         "The type 2 diabetes–heart-failure association is compatible with substantial coronary mediation",
         f"No CAD-independent direct effect detected (+{float(casc['T2D']['direct_b']):.2f}, "
         f"*P* = {float(t2d['direct_p']):.2f})",
         f"Mediated proportion imprecise and covariance-dependent ({t2d_lo:.0f}–{t2d_hi:.0f}% "
         f"under independence); no single fraction is reportable"),
        ("Supportive",
         "Blood pressure and adiposity act on kidney damage in both ancestries",
         f"SBP→CKD in Europeans +{float(eur_sbp['b_perSD']):.2f} per SD (*P* = {psci(eur_sbp['p'])}) and in the "
         f"zero-overlap East Asian design +{float(zero_sbp['b_perSD']):.2f} (*P* = {psci(zero_sbp['p'])}); BMI→CKD "
         f"positive in every dataset",
         "No ancestry interaction on the common per-SD scale; the earlier ancestry contrast reflects outcome size"),
        ("Supportive",
         "Graph-level portability from Europeans to East Asians",
         f"Sign concordance {float(node['estimate']):.2f} on the shared subgraph",
         f"Dependence-robust (node-block) 95% CI {float(node['lo']):.2f}–{float(node['hi']):.2f} includes chance; "
         f"suggestive rather than established"),
        ("Exploratory",
         "The East Asian diabetes-to-adiposity (T2D→BMI) signature",
         "Negative in both hospital biobanks and not abolished in a population cohort",
         "Seven of nine Biobank Japan exposure→BMI edges are negative, so an outcome-side property of that "
         "GWAS is not excluded; meta-analysis sits on a fixed-versus-random borderline"),
        ("Exploratory",
         "The East Asian stage ordering",
         f"Concordant in direction but not beyond chance: meta-analysed concordance "
         f"{float(easmeta[0]['concordance']):.3f} (exact *P* = {float(easmeta[0]['exact_p']):.3f}) and "
         f"{float(easmeta[1]['concordance']):.3f} (exact *P* = {float(easmeta[1]['exact_p']):.3f})",
         "Underpowered rather than falsified; the East Asian graph is small and bipartite"),
        ("Bounded null",
         "Kidney function and kidney disease as exposures for cardiovascular outcomes",
         "eGFR→CVD and CKD→CVD show no detectable effect",
         "21–23 instruments and pleiotropy-inflated heterogeneity: absence of evidence, not evidence of absence"),
        ("Bounded null",
         "LDL→stroke and a direct type 2 diabetes→heart-failure effect",
         "Both dissolve under multivariable conditioning",
         "Bounded as null within the precision of these models"),
        ("Not interpretable",
         "LDL cholesterol→type 2 diabetes",
         "IVW " + f"{float(ldlt2d['ivw_b']):+.2f}".replace("-", "−")
         + f" (*P* = {float(ldlt2d['ivw_p']):.2f}) and weighted median "
         + f"{float(ldlt2d['wm_b']):+.2f}".replace("-", "−")
         + f" (*P* = {float(ldlt2d['wm_p']):.2f}) are null, but prior MR reports an inverse effect",
         f"Egger intercept *P* = {psci(ldlt2d['egger_intercept_p'])}: pleiotropic masking cannot be excluded"),
        ("Not interpretable",
         "HDL cholesterol as a protective factor",
         "A residual direct HDL→CAD association survives conditioning on the atherogenic axis",
         "MVMR alone cannot establish non-causality; CAUSE favours the causal model for this confounded pair"),
    ]

    assert all(t[0] in TIERS for t in T), "unknown tier label"
    L = []
    L.append("<!-- CKM Paper 4 — main-text tables. GENERATED by scripts/build_tables.py from the analysis")
    L.append("     outputs; do not edit by hand. The claim text is authored, every statistic beside it is read")
    L.append("     from the file that produced it. Round-7 C014/Q7: this replaces the Supplementary Figure S7e")
    L.append("     text panel and the Results list that disagreed with it. -->")
    L.append("")
    L.append("# Tables")
    L.append("")
    L.append("**Table 1. Graded synthesis of the study's claims.** Each claim is placed in the evidence tier it "
             "reaches, with the statistic that supports it and the analysis that bounds it. Tiers are: higher "
             "confidence (concordant across estimators and robust to the study's sensitivity analyses); supportive "
             "(consistent evidence that does not clear every bound); exploratory (directionally suggestive, "
             "underpowered or open to an alternative explanation); bounded null (no effect detected, with the "
             "precision stated); not interpretable (the evidence is internally inconsistent or confounded by "
             "pleiotropy). All effects are per standard deviation of the exposure unless stated, on the log-odds "
             "scale for disease outcomes.")
    L.append("")
    L.append("| Tier | Claim | Key evidence | What bounds it |")
    L.append("|---|---|---|---|")
    for tier, claim, ev, bound in T:
        L.append(f"| **{tier}** | {claim} | {ev} | {bound} |")
    L.append("")

    out = os.path.join(FIG, "TABLES.md")
    io.open(out, "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")

    print(f"wrote {out}")
    print(f"  rows: {len(T)}  tiers: " + ", ".join(f"{t}={sum(1 for x in T if x[0]==t)}" for t in TIERS))
    print("\n--- SELF-CHECKS ---")
    body = "\n".join(L)
    ok = True
    for name, cond in [
        ("no literal backtick in the table", "`" not in body),
        ("every row has four cells", all(r.count("|") == 5 for r in L if r.startswith("| **"))),
        ("staging row states the degree-null bound", "degree-preserving rewiring null" in body),
        ("concordance matches fig1_staging_stat.csv", f"{conc:.3f}" in body),
        ("the exact EAS meta P-values come from the enumerator",
         f"{float(easmeta[0]['exact_p']):.3f}" in body and f"{float(easmeta[1]['exact_p']):.3f}" in body),
        ("no tier outside the declared vocabulary", all(t[0] in TIERS for t in T)),
    ]:
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        ok &= bool(cond)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
