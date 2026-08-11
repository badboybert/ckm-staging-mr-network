# -*- coding: utf-8 -*-
"""Step 35c: build UKB-free staging networks (schema == forward_local_edges.csv).
Hybrid = primary network with every edge whose OUTCOME is in {CAD,HF,T2D,Stroke,eGFR}
REPLACED by its UKB-free provider estimate (outcome-side UKB removal). Risk-factor-outcome
edges (into BMI/HDL/TC/TG/LDL/HbA1c/SBP) keep primary values -- no noUKB replacement exists
and they are same-stage or sub-threshold, so they do not drive cross-stage concordance.
Two variants: MAIN (CARDIoGRAM/HERMES/Mahajan/MEGASTROKE/CKDGen) and FINNGEN (FinnGen for
CAD/HF/T2D; MEGASTROKE/CKDGen for Stroke/eGFR which have no FinnGen equivalent)."""
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
RES  = os.path.join(BASE, "results")
PRIMARY = os.path.join(RES, "forward_local_edges.csv")
NOUKB   = os.path.join(RES, "sensitivity_noukb_edges.csv")

HDR = ["exposure","outcome","nsnp","ivw_b","ivw_se","ivw_p","egger_b","egger_p",
       "egger_intercept","egger_intercept_p","wm_b","wm_p","Q","Q_p","steiger_correct","steiger_p"]

SWAP = {"CAD","HF","T2D","Stroke","eGFR"}
MAIN = {"CAD":"CARDIoGRAM","HF":"HERMES","T2D":"Mahajan","Stroke":"MEGASTROKE","eGFR":"CKDGen"}
FINN = {"CAD":"FinnGen","HF":"FinnGen","T2D":"FinnGen","Stroke":"MEGASTROKE","eGFR":"CKDGen"}

def read_csv(p):
    with io.open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))

primary = read_csv(PRIMARY)
noukb   = read_csv(NOUKB)
# index noUKB by (exposure,outcome,provider)
nidx = {(r["exposure"], r["outcome"], r["provider"]): r for r in noukb}

def build(provmap, outpath):
    n_swapped = 0; missing = []
    with io.open(outpath, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f); w.writerow(HDR)
        for r in primary:
            e, o = r["exposure"], r["outcome"]
            if o in SWAP:
                prov = provmap[o]
                nr = nidx.get((e, o, prov))
                if nr is None:
                    missing.append((e, o, prov))
                    row = [r[c] for c in HDR]           # keep primary if no noUKB estimate
                else:
                    row = [nr.get(c, r.get(c, "")) for c in HDR]
                    n_swapped += 1
            else:
                row = [r[c] for c in HDR]
            w.writerow(row)
    print(f"[net] {os.path.basename(outpath)}: {len(primary)} edges, {n_swapped} swapped to UKB-free"
          + (f", MISSING noUKB for {missing}" if missing else ""))
    return outpath

build(MAIN, os.path.join(RES, "sensitivity_noukb_network_main.csv"))
build(FINN, os.path.join(RES, "sensitivity_noukb_network_finngen.csv"))
print("done")
