# -*- coding: utf-8 -*-
"""Step 24: TPMI (Taiwan) vs BBJ (Japan) vs EUR — cross-cohort/cross-ancestry replication of the
forward network. TPMI = independent 2nd EAS cohort. Writes results/tpmi_bbj_eur_comparison.csv + txt."""
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

def load(path, bcol, pcol):
    d={}
    with io.open(os.path.join(BASE,path),encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try: d[(r["exposure"],r["outcome"])]=(float(r[bcol]),float(r[pcol]))
            except (ValueError,KeyError): continue
    return d

tpmi = load("results/network_tpmi_edges.csv","ivw_b","ivw_p")
bbj  = load("results/network_eas_edges.csv","ivw_b","ivw_p")
eur  = load("results/forward_local_edges.csv","ivw_b","ivw_p")

# forward edges present in TPMI
edges = sorted(tpmi.keys())
rows=[]; conc=0; both_sig_conc=0; n_bbj=0
for (e,o) in edges:
    t=tpmi[(e,o)]; b=bbj.get((e,o)); u=eur.get((e,o))
    # concordance TPMI vs BBJ (sign), among edges present in both
    dir_tb = None
    if b is not None:
        n_bbj+=1
        dir_tb = (t[0]>0)==(b[0]>0)
        if dir_tb: conc+=1
        if dir_tb and t[1]<0.05 and b[1]<0.05: both_sig_conc+=1
    rows.append(dict(edge=f"{e}->{o}",
        TPMI_b=round(t[0],3),TPMI_p=f"{t[1]:.1e}",
        BBJ_b=(round(b[0],3) if b else ""),BBJ_p=(f"{b[1]:.1e}" if b else ""),
        EUR_b=(round(u[0],3) if u else ""),EUR_p=(f"{u[1]:.1e}" if u else ""),
        TPMI_BBJ_concordant=("" if b is None else dir_tb)))

out=os.path.join(BASE,"results/tpmi_bbj_eur_comparison.csv")
with io.open(out,"w",encoding="utf-8",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)

lines=["=== TPMI (Taiwan) vs BBJ (Japan) vs EUR — forward network replication ===",
       f"Edges in TPMI: {len(edges)}; shared with BBJ: {n_bbj}",
       f"TPMI-BBJ sign-concordance: {conc}/{n_bbj} = {conc/n_bbj:.2f}" if n_bbj else "",
       f"TPMI-BBJ both-significant & concordant: {both_sig_conc}",
       "",
       f"{'edge':14s} {'TPMI_b':>8} {'TPMI_p':>9} {'BBJ_b':>7} {'BBJ_p':>9} {'EUR_b':>7} {'EUR_p':>9}"]
for r in rows:
    lines.append(f"{r['edge']:14s} {str(r['TPMI_b']):>8} {r['TPMI_p']:>9} {str(r['BBJ_b']):>7} {str(r['BBJ_p']):>9} {str(r['EUR_b']):>7} {str(r['EUR_p']):>9}")
# highlight the power-rescue edges (sig in TPMI, was NS/absent in BBJ)
lines.append("\n★ POWER-RESCUE candidates (TPMI-significant, BBJ NS or absent):")
for (e,o) in edges:
    t=tpmi[(e,o)]; b=bbj.get((e,o))
    if t[1]<0.05 and (b is None or b[1]>=0.05):
        lines.append(f"   {e}->{o}: TPMI b={t[0]:+.3f} p={t[1]:.1e} | BBJ {'absent' if b is None else f'b={b[0]:+.3f} p={b[1]:.1e} (NS)'}")
txt="\n".join(lines)
io.open(os.path.join(BASE,"results/tpmi_bbj_eur_comparison.txt"),"w",encoding="utf-8").write(txt)
print(txt)

if __name__=="__main__": pass
