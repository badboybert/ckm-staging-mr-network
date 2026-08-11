# -*- coding: utf-8 -*-
"""Task 1: annotate every headline edge with native UNITS and cross-ancestry COMPARABILITY,
and put SBP edges on a common per-SD scale so magnitudes are comparable.
SBP is the only mismatch that is RESCALED here: EUR-SBP (Evangelou, per-mmHg) vs EAS-BBJ-SBP
(Kanai, rank-INT per-SD) -> multiply EUR by SD_SBP to reach per-SD. Two OTHER cross-ancestry
unit mismatches exist and are handled by flagging comparable=False (not rescaled, because no
clean conversion): HbA1c (EUR per-unit MAGIC vs EAS rank-INT per-SD) and eGFR (EUR log-per-SD
vs EAS per-SD). All remaining continuous exposures (BMI, LDL, HDL, TC, TG) are per-SD in both
ancestries, and binary outcomes are log-OR in both. Writes results/edges_standardised.csv.
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

import csv, io, os
BASE = P4_BASE
SD_SBP = 19.3   # mmHg, adult SBP SD (EUR/UKB ~19-20; BBJ ~18.6). Approx; used to rescale EUR SBP->per-SD.

# native units per node when used as EXPOSURE (EUR / EAS)
EXP_UNITS = {
 "BMI":("per-SD","per-SD"), "SBP":("per-mmHg","per-SD (INT)"),
 "LDL":("per-SD","per-SD"), "HDL":("per-SD","per-SD"), "TC":("per-SD","per-SD"), "TG":("per-SD","per-SD"),
 "HbA1c":("per-unit (MAGIC)","per-SD (INT)"), "T2D":("per-log-OR","per-log-OR"),
 "eGFR":("per-SD (log)","per-SD"), "CAD":("per-log-OR","per-log-OR"), "HF":("per-log-OR","per-log-OR"),
}
BIN = {"CAD","HF","Stroke","T2D","CKD"}

def out_units(o): return "log-OR" if o in BIN else "per-SD"

def parse_edge(e):
    a,b = e.split("->"); return a,b

def _full_precision(path, key_exp="exposure", key_out="outcome", key_b="ivw_b"):
    """Map (exposure, outcome) -> full-precision IVW beta from a committed network table."""
    d = {}
    if not os.path.exists(path):
        return d
    with io.open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                d[(r[key_exp], r[key_out])] = float(r[key_b])
            except (KeyError, ValueError):
                continue
    return d


def main():
    src = os.path.join(BASE,"results/eur_vs_eas_comparison.csv")
    # eur_vs_eas_comparison.csv stores betas ROUNDED to 3 decimals. Rescaling a rounded value put
    # SBP->CAD on the per-SD scale at 19.3 x 0.034 = 0.6562 instead of 19.3 x 0.0335822 = 0.6481.
    # Take the beta from the network tables, at full precision, and fall back to the comparison file
    # only where an edge cannot be resolved there.
    EUR_FULL = _full_precision(os.path.join(BASE, "results/forward_local_edges.csv"))
    EAS_FULL = _full_precision(os.path.join(BASE, "results/network_eas_edges.csv"))
    n_full_eur = n_full_eas = 0
    rows=[]
    with io.open(src,encoding="utf-8") as f:
        for r in csv.DictReader(f):
            a,b = parse_edge(r["edge"])
            eu_units = EXP_UNITS.get(a,("per-SD","per-SD"))
            exp_u_eur, exp_u_eas = eu_units
            out_u = out_units(b)
            # comparability: True unless exposure units differ across ancestry (only SBP)
            comparable = (exp_u_eur == exp_u_eas)
            eur_b = float(r["EUR_b"]); eas_b = float(r["EAS_b"])
            if (a, b) in EUR_FULL:
                assert abs(EUR_FULL[(a, b)] - eur_b) < 5e-3, \
                    "EUR beta disagrees with the network table for %s->%s" % (a, b)
                eur_b = EUR_FULL[(a, b)]; n_full_eur += 1
            if (a, b) in EAS_FULL:
                assert abs(EAS_FULL[(a, b)] - eas_b) < 5e-3, \
                    "EAS beta disagrees with the network table for %s->%s" % (a, b)
                eas_b = EAS_FULL[(a, b)]; n_full_eas += 1
            eur_b_std, eas_b_std, note = eur_b, eas_b, ""
            if a=="SBP":
                # put EUR SBP on per-SD to match EAS: multiply per-mmHg beta by SD_SBP
                eur_b_std = eur_b * SD_SBP
                note = f"EUR rescaled per-mmHg->per-SD (x{SD_SBP})"
                comparable = True
            rows.append(dict(edge=r["edge"], exposure=a, outcome=b,
                exp_units_eur=exp_u_eur, exp_units_eas=exp_u_eas, out_units=out_u,
                EUR_b=eur_b, EAS_b=eas_b, EUR_b_perSD=round(eur_b_std,4), EAS_b_perSD=round(eas_b_std,4),
                EUR_p=r["EUR_p"], EAS_p=r["EAS_p"], dir_concordant=r["dir_concordant"],
                both_sig=r["both_sig"], comparable=comparable, scale_note=note))
    out=os.path.join(BASE,"results/edges_standardised.csv")
    with io.open(out,"w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows: w.writerow(r)
    # show the SBP edges rescaled
    print("full-precision betas taken from the network tables: EUR %d, EAS %d of %d edges"
          % (n_full_eur, n_full_eas, len(rows)))
    print("=== SBP edges on a common per-SD scale (the only cross-ancestry scale fix) ===")
    print(f"{'edge':14s} {'EUR/mmHg':>9} {'EUR/SD':>8} {'EAS/SD':>8}  ratio(EAS/EUR)")
    for r in rows:
        if r["exposure"]=="SBP" and r["both_sig"]=="True":
            ratio = r["EAS_b_perSD"]/r["EUR_b_perSD"] if r["EUR_b_perSD"] else float('nan')
            print(f"{r['edge']:14s} {r['EUR_b']:>9.3f} {r['EUR_b_perSD']:>8.3f} {r['EAS_b_perSD']:>8.3f}  {ratio:>6.2f}x")
    print(f"\n-> {out} ({len(rows)} edges annotated)")

if __name__=="__main__": main()
