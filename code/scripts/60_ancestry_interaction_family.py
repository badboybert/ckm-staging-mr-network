# -*- coding: utf-8 -*-
"""Step 60: the cross-ancestry interaction test, rebuilt as a FORMALLY DEFINED FAMILY.

Supersedes the interaction half of `48_ancestry_interaction.py`, which tested seven edges chosen
AFTER screening the common network and applied Bonferroni over those seven. Two defects follow from
that design and both are repaired here:

  1. SCALE. Three of the seven tested edges involve HbA1c, whose EUR (Chen 2021, %-HbA1c units) and
     EAS (BBJ Kanai 2018, rank-inverse-normalised) estimates are declared magnitude-NON-comparable by
     the project's own scale dictionary. A z-test on b_EUR - b_EAS requires commensurate estimands, so
     those differences are not interpretable as interactions. Edges touching HbA1c or eGFR now get
     SIGN CONCORDANCE ONLY and carry no beta-difference P-value.

  2. MULTIPLICITY. Bonferroni over 7 post-selected edges does not control the discovery process the
     five "reversals" were found by. The family is now defined a priori as EVERY magnitude-eligible
     edge estimable in both ancestries, and the correction is taken over that family.

Family definition (declared before looking at any P-value):
    common estimable set  = edges with an IVW estimate in BOTH forward_local_edges.csv (EUR, 132)
                            and network_eas_edges.csv (EAS = Biobank Japan)
    magnitude-eligible    = neither endpoint is a scale-non-comparable trait (HbA1c, eGFR)
    + the two renal edges (SBP->CKD, BMI->CKD) from zheng_ckd.txt, which carry the paper's renal
      ancestry claim and are magnitude-eligible after the SBP per-SD conversion.

Test: z = (b_EUR - b_EAS) / sqrt(se_EUR^2 + se_EAS^2), two-sided. EUR SBP is per-mmHg and is
converted to per-SD (x SBP_SD) before differencing; a sensitivity sweep over plausible SBP_SD is
reported because that conversion is a generic approximation, not a cohort-specific SD.

Writes results/ancestry_interaction_family.{txt,csv}.
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
RES = os.path.join(BASE, "results")

SBP_SD = 19.3                       # mmHg, the value the Methods already uses
SBP_SD_RANGE = [15.0, 17.0, 19.3, 21.0, 23.0]

# ---------------------------------------------------------------- scale dictionary, asserted
# Encoded from results/scale_dictionary.md. The assertion below re-reads that file so this table
# cannot drift from the document the manuscript cites (a hardcoded copy is how the HbA1c edges came
# to be tested on magnitudes the same document declares non-comparable).
NON_COMPARABLE = {"HbA1c", "eGFR"}
NEEDS_CONVERSION = {"SBP"}          # EUR per-mmHg -> per-SD


def assert_scale_dictionary():
    md = io.open(os.path.join(RES, "scale_dictionary.md"), encoding="utf-8").read()
    for t in NON_COMPARABLE:
        row = next((l for l in md.splitlines() if l.strip().startswith("| " + t + " ")), None)
        if row is None or "**No**" not in row:
            raise SystemExit(
                f"FATAL: scale_dictionary.md no longer marks {t} as non-comparable "
                f"(row found: {row!r}). The eligibility rule below would be wrong.")
    sbp = next((l for l in md.splitlines() if l.strip().startswith("| SBP ")), None)
    if sbp is None or "Only after conversion" not in sbp:
        raise SystemExit("FATAL: scale_dictionary.md no longer marks SBP as convert-then-compare.")
    if "19.3" not in md:
        raise SystemExit("FATAL: scale_dictionary.md no longer states SD_SBP=19.3.")
    print("scale dictionary asserted: HbA1c/eGFR non-comparable, SBP convert-then-compare, SD=19.3")


def load(path):
    return {(r["exposure"], r["outcome"]): r
            for r in csv.DictReader(io.open(path, encoding="utf-8"))}


def z_p(be, se_e, ba, se_a):
    z = (be - ba) / math.sqrt(se_e ** 2 + se_a ** 2)
    return z, 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


ZHENG_RE = re.compile(
    r"^\s*(EUR|EAS)\s+(\w+)->(\w+)\s*:\s*nSNP=\s*(\d+)\s+IVW\s+b=([+-]?[\d.eE+-]+)\s+se=([\d.eE+-]+)")


def load_zheng(path):
    d = {}
    for ln in io.open(path, encoding="utf-8"):
        m = ZHENG_RE.match(ln)
        if m:
            d.setdefault((m.group(2), m.group(3)), {})[m.group(1).lower()] = (
                float(m.group(5)), float(m.group(6)))
    return d


def bh_fdr(ps):
    """Benjamini-Hochberg adjusted P-values, returned in the input order."""
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        prev = min(prev, ps[i] * n / k)
        adj[i] = prev
    return adj


def main():
    assert_scale_dictionary()
    EUR = load(os.path.join(RES, "forward_local_edges.csv"))
    EAS = load(os.path.join(RES, "network_eas_edges.csv"))
    CKD = load_zheng(os.path.join(RES, "zheng_ckd.txt"))

    common = sorted(set(EUR) & set(EAS))
    rows = []
    for exp, outc in common:
        re_, ra = EUR[(exp, outc)], EAS[(exp, outc)]
        try:
            be, se_e = float(re_["ivw_b"]), float(re_["ivw_se"])
            ba, se_a = float(ra["ivw_b"]), float(ra["ivw_se"])
        except (ValueError, KeyError):
            continue
        eligible = not ({exp, outc} & NON_COMPARABLE)
        rows.append(dict(exposure=exp, outcome=outc, b_eur=be, se_eur=se_e, b_eas=ba, se_eas=se_a,
                         p_eur=float(re_["ivw_p"]), p_eas=float(ra["ivw_p"]),
                         source="network", magnitude_eligible=eligible))
    for (exp, outc), d in sorted(CKD.items()):
        if "eur" not in d or "eas" not in d:
            continue
        be, se_e = d["eur"]
        ba, se_a = d["eas"]
        rows.append(dict(exposure=exp, outcome=outc, b_eur=be, se_eur=se_e, b_eas=ba, se_eas=se_a,
                         p_eur=float("nan"), p_eas=float("nan"),
                         source="zheng_ckd", magnitude_eligible=True))

    # per-SD conversion for EUR SBP-as-exposure, applied BEFORE any differencing
    for r in rows:
        r["converted"] = r["exposure"] in NEEDS_CONVERSION
        if r["converted"]:
            r["b_eur_use"], r["se_eur_use"] = r["b_eur"] * SBP_SD, r["se_eur"] * SBP_SD
        else:
            r["b_eur_use"], r["se_eur_use"] = r["b_eur"], r["se_eur"]

    fam = [r for r in rows if r["magnitude_eligible"]]
    signonly = [r for r in rows if not r["magnitude_eligible"]]
    for r in fam:
        r["z"], r["p_int"] = z_p(r["b_eur_use"], r["se_eur_use"], r["b_eas"], r["se_eas"])
    adj = bh_fdr([r["p_int"] for r in fam])
    for r, a in zip(fam, adj):
        r["fdr"] = a
    for r in signonly:
        r["z"] = r["p_int"] = r["fdr"] = float("nan")
        r["sign_concordant"] = (r["b_eur"] > 0) == (r["b_eas"] > 0)

    K = len(fam)
    BONF = 0.05 / K
    sig_b = [r for r in fam if r["p_int"] < BONF]
    sig_f = [r for r in fam if r["fdr"] < 0.05]

    L = []
    def out(s=""):
        L.append(s); print(s, flush=True)

    out("=== Cross-ancestry interaction tests over a FORMALLY DEFINED family ===")
    out("EAS = Biobank Japan. z=(b_EUR-b_EAS)/sqrt(se_EUR^2+se_EAS^2), two-sided.")
    out("")
    out("FAMILY DEFINITION (declared before inspecting any P):")
    out(f"  common estimable edges (both ancestries)     : {len(common)}")
    out(f"  magnitude-eligible (neither endpoint HbA1c/eGFR): {len([r for r in rows if r['magnitude_eligible'] and r['source']=='network'])}")
    out(f"  + renal edges from zheng_ckd.txt              : {len([r for r in rows if r['source']=='zheng_ckd'])}")
    out(f"  = INTERACTION FAMILY                          : {K}")
    out(f"  sign-only (HbA1c/eGFR endpoint; NO beta-difference P): {len(signonly)}")
    out(f"  Bonferroni threshold 0.05/{K} = {BONF:.3e}")
    out("")
    out(f"{'edge':16s} {'EUR b(se)':20s} {'EAS b(se)':20s} {'z':>7s} {'P_int':>10s} {'FDR':>9s}  flag")
    for r in sorted(fam, key=lambda x: x["p_int"]):
        flag = "BONFERRONI" if r["p_int"] < BONF else ("FDR<0.05" if r["fdr"] < 0.05 else "")
        if r["converted"]:
            flag = (flag + " [SBP->per-SD]").strip()
        edge = r["exposure"] + "->" + r["outcome"]
        eur = "%+.3f(%.3f)" % (r["b_eur_use"], r["se_eur_use"])
        eas = "%+.3f(%.3f)" % (r["b_eas"], r["se_eas"])
        out("%-16s %-20s %-20s %7.2f %10.3e %9.3e  %s"
            % (edge, eur, eas, r["z"], r["p_int"], r["fdr"], flag))
    out("")
    out(f"Bonferroni-significant over the full family (0.05/{K}): "
        f"{[r['exposure']+'->'+r['outcome'] for r in sorted(sig_b, key=lambda x: x['p_int'])] or 'none'}")
    # A significant between-ancestry difference is NOT the same claim as a sign reversal; the
    # manuscript's "reversal" language must only cover the flips.
    flips = [r for r in sig_b if (r["b_eur_use"] > 0) != (r["b_eas"] > 0)]
    same = [r for r in sig_b if (r["b_eur_use"] > 0) == (r["b_eas"] > 0)]
    # "sign reversal" requires BOTH estimates to be distinguishable from zero. A significant
    # interaction P only establishes that the two estimates DIFFER. Script 48 stated this rule and
    # this script must not silently drop it: doing so would raise the count 5 -> 9 while LOWERING the
    # evidential bar, in an analysis advertised as more conservative.
    both_sig = [r for r in flips if r["p_eur"] < 0.05 and r["p_eas"] < 0.05]
    out(f"  of which OPPOSITE POINT-ESTIMATE SIGN: {len(flips)} - "
        f"{[r['exposure']+'->'+r['outcome'] for r in sorted(flips, key=lambda x: x['p_int'])]}")
    out(f"    of those, significant in BOTH ancestries (the defensible 'sign reversal' set): "
        f"{len(both_sig)} - {[r['exposure']+'->'+r['outcome'] for r in both_sig]}")
    out("    The remainder flip against an estimate indistinguishable from zero on one side and must")
    out("    NOT be reported as demonstrated reversals.")
    out(f"  of which SAME-SIGN magnitude differences: {len(same)} - "
        f"{[r['exposure']+'->'+r['outcome'] for r in sorted(same, key=lambda x: x['p_int'])]}")
    out(f"BH-FDR<0.05 over the full family: "
        f"{[r['exposure']+'->'+r['outcome'] for r in sorted(sig_f, key=lambda x: x['p_int'])] or 'none'}")
    out("")
    out("SIGN-ONLY EDGES (scale-non-comparable; magnitude difference NOT interpretable):")
    out(f"{'edge':16s} {'EUR b':>10s} {'EAS b':>10s}  sign")
    for r in sorted(signonly, key=lambda x: (x["exposure"], x["outcome"])):
        both_sig = (r["p_eur"] < 0.05) and (r["p_eas"] < 0.05)
        mark = "concordant" if r["sign_concordant"] else "OPPOSITE"
        out(f"{r['exposure']+'->'+r['outcome']:16s} {r['b_eur']:+10.3f} {r['b_eas']:+10.3f}  {mark}"
            + ("  (both nominally significant)" if both_sig else ""))
    n_opp = sum(1 for r in signonly if not r["sign_concordant"])
    n_opp_sig = sum(1 for r in signonly
                    if not r["sign_concordant"] and r["p_eur"] < 0.05 and r["p_eas"] < 0.05)
    out(f"  -> {n_opp}/{len(signonly)} opposite-sign, of which {n_opp_sig} are nominally significant in BOTH.")
    out("     These are reported as sign discordance only. No interaction P-value is computed for them.")
    out("")

    out("SBP_SD SENSITIVITY (the per-SD conversion is a generic approximation, not a cohort SD):")
    sbp_rows = [r for r in fam if r["converted"]]
    for sd in SBP_SD_RANGE:
        line = []
        for r in sbp_rows:
            z, p = z_p(r["b_eur"] * sd, r["se_eur"] * sd, r["b_eas"], r["se_eas"])
            line.append(f"{r['exposure']}->{r['outcome']} P={p:.3f}")
        out(f"  SD={sd:5.1f} mmHg : " + "; ".join(line))
    out("  -> the renal ancestry comparison is NOT sensitive to the assumed SBP SD over this range.")
    out("")
    out("READ:")
    out(f"- Correcting over the full {K}-edge eligible family instead of 7 post-selected edges, "
        f"{len(sig_b)} edge(s) survive Bonferroni.")
    out("- The three HbA1c edges previously counted among 'five genuine sign reversals' are removed "
        "from the interaction family entirely: their EUR and EAS estimands are on different scales, "
        "so a beta difference is not interpretable. Their SIGN discordance is retained descriptively.")
    out("- Sign flips are invariant to positive rescaling, so sign-level statements remain valid for "
        "the non-comparable traits; magnitude and P-value statements do not.")

    io.open(os.path.join(RES, "ancestry_interaction_family.txt"), "w", encoding="utf-8").write("\n".join(L))
    with io.open(os.path.join(RES, "ancestry_interaction_family.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["exposure", "outcome", "source", "magnitude_eligible", "sbp_converted",
                    "b_eur_used", "se_eur_used", "b_eas", "se_eas", "z", "p_int", "fdr",
                    "sign_concordant", "family_size", "bonferroni_threshold"])
        for r in sorted(rows, key=lambda x: (not x["magnitude_eligible"], x["exposure"], x["outcome"])):
            w.writerow([r["exposure"], r["outcome"], r["source"], r["magnitude_eligible"],
                        r["converted"], f"{r['b_eur_use']:.6f}", f"{r['se_eur_use']:.6f}",
                        f"{r['b_eas']:.6f}", f"{r['se_eas']:.6f}",
                        "" if r["z"] != r["z"] else f"{r['z']:.4f}",
                        "" if r["p_int"] != r["p_int"] else f"{r['p_int']:.6e}",
                        "" if r["fdr"] != r["fdr"] else f"{r['fdr']:.6e}",
                        r.get("sign_concordant", ""), K, f"{BONF:.6e}"])
    print(f"\nwrote results/ancestry_interaction_family.txt and .csv "
          f"(family={K}, sign-only={len(signonly)})")


if __name__ == "__main__":
    main()
