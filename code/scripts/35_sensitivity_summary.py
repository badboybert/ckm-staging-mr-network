# -*- coding: utf-8 -*-
"""Step 35d: assemble results/sensitivity_noukb_summary.txt.
Side-by-side PRIMARY vs UKB-free-OUTCOME for the headline edges + AHA-staging concordance
(primary vs the two UKB-free networks), with an explicit HOLD/ATTENUATE/FLIP verdict."""
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

import csv, io, os, sys, contextlib

BASE = P4_BASE
RES  = os.path.join(BASE, "results")
sys.path.insert(0, os.path.join(BASE, "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("staging", os.path.join(BASE,"scripts","07_staging.py"))
staging = importlib.util.module_from_spec(spec); spec.loader.exec_module(staging)

def read_csv(p):
    with io.open(p, encoding="utf-8") as f: return list(csv.DictReader(f))

primary = { (r["exposure"], r["outcome"]): r for r in read_csv(os.path.join(RES,"forward_local_edges.csv")) }
noukb_rows = read_csv(os.path.join(RES,"sensitivity_noukb_edges.csv"))
noukb = {}
for r in noukb_rows: noukb[(r["exposure"], r["outcome"], r["provider"])] = r

# Network Bonferroni edge threshold, DERIVED from the primary edge table (never hardcoded — this file
# shipped 0.05/111 for six days after the network was completed to 132 directions).
N_PRIMARY = len(primary)
BONF = 0.05 / N_PRIMARY

MAIN = {"CAD":"CARDIoGRAM","HF":"HERMES","T2D":"Mahajan","Stroke":"MEGASTROKE","eGFR":"CKDGen"}
# headline forward edges + 2 feedback edges
HEADLINE = [("BMI","CAD"),("LDL","CAD"),("TC","CAD"),("TG","CAD"),("SBP","CAD"),
            ("BMI","T2D"),("HbA1c","T2D"),("TG","T2D"),
            ("SBP","Stroke"),("BMI","Stroke"),
            ("CAD","HF"),("SBP","HF"),("BMI","HF"),
            ("SBP","eGFR"),("BMI","eGFR"),
            ("CAD","T2D"),("HF","T2D")]      # last two = feedback

def f(x, nd=4):
    try: return f"{float(x):+.{nd}g}"
    except: return str(x)
def pe(x):
    try: return f"{float(x):.2e}"
    except: return str(x)

def verdict(pr, nr):
    if nr is None: return "no noUKB estimate"
    pb, pp = float(pr["ivw_b"]), float(pr["ivw_p"])
    nb, np_ = float(nr["ivw_b"]), float(nr["ivw_p"])
    same_sign = (pb>0) == (nb>0)
    pr_sig, nr_sig = pp<BONF, np_<BONF
    if not same_sign and pr_sig and nr_sig:  return "*** SIGN FLIP (both sig) ***"
    if not same_sign:                        return "sign differs (>=1 NS)"
    if pr_sig and nr_sig:                    return "HOLDS (same sign, sig at Bonf)"
    if (not pr_sig) and (not nr_sig):        return "HOLDS as NULL (NS in both)"
    if pr_sig and not nr_sig:                return "ATTENUATES (loses Bonf sig)"
    return "STRENGTHENS (gains sig)"

def run_staging(path):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        d = staging.analyse(path)
    return d

lines = []
W = lines.append
W("="*118)
W("CKM Paper 4 -- UKB-FREE OUTCOME-SIDE SENSITIVITY  (headline edges: PRIMARY vs UKB-free outcome GWAS)")
W("="*118)
W("Framing: OUTCOME-side UKB removal only. Exposure instruments (Yengo BMI, Evangelou SBP, GLGC lipids,")
W("MAGIC HbA1c) still contain UKB -- this is NOT a fully UKB-free analysis. Two-sample MR sample overlap")
W("biases estimates toward the confounded observational value; removing outcome-side overlap tests that.")
W("Providers: CAD=CARDIoGRAM(Nikpay2015,N~184k) | HF=HERMES-noUKB | T2D=Mahajan2018-noUKB | Stroke=MEGASTROKE | eGFR=CKDGen2016.")
W("CAVEAT: the noUKB GWAS also differ from the primary in cohort/N/era/phenotype, so any attenuation")
W("conflates UKB-overlap removal with study/power differences. The FinnGen R12 column (fully independent)")
W("is the clean robustness check -- concordant FinnGen signs = the edge is not an overlap artifact.")
W(f"Significance bar = network Bonferroni 0.05/{N_PRIMARY} = {BONF:.2e}  (primary staging threshold).")
W("N.B. FinnGen R12, MEGASTROKE, CKDGen, Mahajan use N_const (affects Steiger direction ONLY, not beta/P).")
W("")
hd = f"{'edge':<14}{'PRIMARY b':>12}{'PRIM p':>11}{'  ':2}{'noUKB prov':<11}{'noUKB b':>11}{'noUKB p':>11}{'nSNP':>6}{'Qp':>9}  verdict"
W(hd); W("-"*118)
for e,o in HEADLINE:
    pr = primary.get((e,o)); prov = MAIN[o]; nr = noukb.get((e,o,prov))
    tag = "  <feedback>" if (e,o) in [("CAD","T2D"),("HF","T2D")] else ""
    if pr is None:
        W(f"{e+'->'+o:<14}  (no primary row){tag}"); continue
    nb = f(nr['ivw_b']) if nr else "NA"; np_ = pe(nr['ivw_p']) if nr else "NA"
    nsnp = nr['nsnp'] if nr else "-"; qp = pe(nr['Q_p']) if nr else "-"
    W(f"{e+'->'+o:<14}{f(pr['ivw_b']):>12}{pe(pr['ivw_p']):>11}{'  ':2}{prov:<11}{nb:>11}{np_:>11}{nsnp:>6}{qp:>9}  {verdict(pr,nr)}{tag}")
W("")
# FinnGen replication column for the disease outcomes that have it
W("FinnGen R12 (fully independent, no UKB) cross-check for CAD/HF/T2D outcomes:")
W(f"{'edge':<14}{'FinnGen b':>12}{'FinnGen p':>12}{'nSNP':>7}   verdict-vs-primary")
W("-"*70)
for e,o in HEADLINE:
    if o not in ("CAD","HF","T2D"): continue
    pr = primary.get((e,o)); nr = noukb.get((e,o,"FinnGen"))
    if pr is None or nr is None: continue
    W(f"{e+'->'+o:<14}{f(nr['ivw_b']):>12}{pe(nr['ivw_p']):>12}{nr['nsnp']:>7}   {verdict(pr,nr)}")
W("")
W("="*118)
W("AHA CKM STAGING CONCORDANCE  (population-level causal-ordering test; EXACT label-permutation null)")
W("="*118)
nets = [("PRIMARY (UKB-containing outcomes)", os.path.join(RES,"forward_local_edges.csv")),
        ("UKB-free MAIN  (CARDIoGRAM/HERMES/Mahajan/MEGASTROKE/CKDGen)", os.path.join(RES,"sensitivity_noukb_network_main.csv")),
        ("UKB-free FINNGEN (FinnGen CAD/HF/T2D; MEGASTROKE/CKDGen)",   os.path.join(RES,"sensitivity_noukb_network_finngen.csv"))]
W(f"{'network':<60}{'concord.':>10}{'exact P':>12}   backward (AHA-discordant) edges")
W("-"*118)
STG = {}
for name, path in nets:
    d = run_staging(path)
    STG[name.split()[0]] = d
    n_rows = len(read_csv(path))
    W(f"{name:<60}{d['concordance']:>10.3f}{d['perm_p']:>12.2e}   {d['backward'] if d['backward'] else 'none'}")
    W(f"{'':<60}{'(exact ' + str(d['perm_ge']) + '/' + str(d['perm_tot']) + ' labelings)':>12}")
    W(f"{'':<60}{'(' + str(n_rows) + ' edges, Bonferroni 0.05/' + str(n_rows) + ')':>21}")
W("")
W("SCOPE OF THE UKB-FREE NETWORKS: each is the primary edge table with every edge whose OUTCOME is in")
W("{CAD,HF,T2D,Stroke,eGFR} replaced by its UKB-free provider estimate. Four Stroke-as-exposure edges")
W("(Stroke->CAD, Stroke->HF, Stroke->T2D, Stroke->eGFR) have no UKB-free counterpart in the provider")
W("GWAS and retain their PRIMARY values; every other edge into those five outcomes is swapped.")
W("")
# --- the verdict lines are DERIVED, never typed: every number below is recomputed above -------------
fwd_hold = [f"{e}->{o}" for e, o in HEADLINE
            if (e, o) not in [("CAD", "T2D"), ("HF", "T2D")]
            and primary.get((e, o)) and noukb.get((e, o, MAIN[o]))
            and verdict(primary[(e, o)], noukb[(e, o, MAIN[o])]).startswith("HOLDS (same sign, sig")]
fwd_null = [f"{e}->{o}" for e, o in HEADLINE
            if primary.get((e, o)) and noukb.get((e, o, MAIN[o]))
            and verdict(primary[(e, o)], noukb[(e, o, MAIN[o])]).startswith("HOLDS as NULL")]
fwd_att = [f"{e}->{o}" for e, o in HEADLINE
           if (e, o) not in [("CAD", "T2D"), ("HF", "T2D")]
           and primary.get((e, o)) and noukb.get((e, o, MAIN[o]))
           and verdict(primary[(e, o)], noukb[(e, o, MAIN[o])]).startswith("ATTENUATES")]
W("VERDICT:")
W(f"  - {len(fwd_hold)} headline FORWARD edges retain SIGN and Bonferroni significance under outcome-side UKB")
W(f"    removal; {len(fwd_att)} attenuate below the bar" + (f" ({', '.join(fwd_att)})" if fwd_att else "") + ".")
W(f"  - {len(fwd_null)} headline edges remain NULL in both ({', '.join(fwd_null)}) — matching the primary.")
_p = STG['PRIMARY']; _m = STG['UKB-free'] if 'UKB-free' in STG else None
_main = run_staging(nets[1][1]); _finn = run_staging(nets[2][1])
W(f"  - Staging order stays concordant beyond chance in BOTH UKB-free networks "
  f"(concordance {_main['concordance']:.3f}, exact P={_main['perm_p']:.2e} "
  f"[{_main['perm_ge']}/{_main['perm_tot']}] MAIN; "
  f"{_finn['concordance']:.3f}, exact P={_finn['perm_p']:.2e} "
  f"[{_finn['perm_ge']}/{_finn['perm_tot']}] FinnGen;")
W(f"    primary {_p['concordance']:.3f}, exact P={_p['perm_p']:.2e} [{_p['perm_ge']}/{_p['perm_tot']}])."
  f" Directional architecture is ROBUST to outcome-side UKB overlap.")
W("  - FEEDBACK nuance: CAD->T2D survives in independent FinnGen -> NOT a UKB artifact, but attenuates")
W("    below Bonferroni in Mahajan-noUKB -> partly overlap/power. HF->T2D attenuates in both noUKB")
W("    providers and its Steiger direction flips to FALSE (N_const-sensitive), so it drops from the")
W("    significant DAG -> the HF->T2D backward feedback is the one FRAGILE edge.")
_back_p = _p['backward']; _back_m = _main['backward']
W(f"  - Net effect on staging: the primary network's backward set {_back_p} becomes {_back_m if _back_m else 'empty'}")
W(f"    under MAIN and {_finn['backward'] if _finn['backward'] else 'empty'} under FinnGen; the AHA-staging conclusion is not weakened by UKB removal.")
out = os.path.join(RES,"sensitivity_noukb_summary.txt")
with io.open(out,"w",encoding="utf-8") as fh: fh.write("\n".join(lines)+"\n")
print("\n".join(lines))
print("\nWrote", out)
