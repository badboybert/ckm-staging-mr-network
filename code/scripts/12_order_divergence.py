# -*- coding: utf-8 -*-
"""Component D: cross-ancestry ORDER-divergence statistic (graph-level; scale-free).
Zheng2021 already owns EDGE-level cross-ancestry divergence; the unclaimed leg is a
GRAPH-level node-ordering divergence. We derive, per ancestry, a scale-FREE node causal
'depth' score from the SIGN + SIGNIFICANCE of edges (net downstream flow), rank the shared
nodes, and compute Kendall-tau between the EUR and EAS orderings with a permutation null.
Scale-free because it uses only edge direction/significance, never magnitudes (this sidesteps
the BBJ-per-SD vs EUR-per-mmHg SBP scale problem). We also report the simpler, robust
edge-level SIGN-concordance rate on well-powered shared edges (binomial null), and run a
sensitivity excluding the flagged ->BMI artifact edges and SBP-origin edges.

Inputs: results/network_forward_20260711.csv (EUR bidirectional) OR forward_local_edges.csv,
        results/network_eas_edges.csv (EAS).
Output: results/order_divergence.csv + .txt
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

import csv, io, os, math, itertools
BASE = P4_BASE
SHARED = ["BMI","SBP","LDL","HDL","TC","TG","HbA1c","T2D","eGFR","CAD","HF","Stroke"]
SIG = 0.05
MIN_SNP = 10

def load(path, bcol, pcol, ncol, scol):
    d={}
    with io.open(path,encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                b=float(r[bcol]); p=float(r[pcol]); n=int(float(r[ncol]))
                st=str(r.get(scol,"")).upper()=="TRUE"
                d[(r["exposure"],r["outcome"])]=dict(b=b,p=p,n=n,st=st)
            except (ValueError,KeyError): continue
    return d

def netflow(edges, nodes, drop_to_bmi=False, drop_sbp=False):
    """scale-free node score: (#sig outgoing) - (#sig incoming) among Steiger-correct edges."""
    score={k:0 for k in nodes}
    for (a,b),e in edges.items():
        if a not in nodes or b not in nodes: continue
        if drop_to_bmi and b=="BMI": continue
        if drop_sbp and a=="SBP": continue
        if e["p"]<SIG and e["st"] and e["n"]>=MIN_SNP:
            score[a]+=1; score[b]-=1
    return score

def rankvec(score, nodes):
    # higher netflow = more upstream. rank ascending by score -> ordinal position
    order=sorted(nodes, key=lambda k:score[k])
    return {k:i for i,k in enumerate(order)}, score

def kendall_tau(x, y):
    keys=list(x); nc=nd=0
    for i,j in itertools.combinations(keys,2):
        dx=x[i]-x[j]; dy=y[i]-y[j]
        s=dx*dy
        if s>0: nc+=1
        elif s<0: nd+=1
    tot=nc+nd
    return (nc-nd)/tot if tot else float("nan"), nc, nd

def perm_p(x, y, tau_obs, nperm=20000):
    import random
    rng=random.Random(12345)
    keys=list(x); yv=[y[k] for k in keys]
    ge=0
    for _ in range(nperm):
        rng.shuffle(yv)
        yp={k:yv[i] for i,k in enumerate(keys)}
        t,_,_=kendall_tau(x,yp)
        if abs(t)>=abs(tau_obs)-1e-12: ge+=1
    return (ge+1)/(nperm+1)

def sign_concordance(eur, eas, drop_to_bmi=False, drop_sbp=False):
    pairs=[]
    for (a,b) in set(eur)&set(eas):
        if a not in SHARED or b not in SHARED: continue
        if drop_to_bmi and b=="BMI": continue
        if drop_sbp and a=="SBP": continue
        e1,e2=eur[(a,b)],eas[(a,b)]
        # well-powered anchor: EUR nominally sig + both have instruments
        if e1["p"]<SIG and e1["n"]>=MIN_SNP and e2["n"]>=MIN_SNP:
            pairs.append(((a,b), (e1["b"]>0)==(e2["b"]>0), e2["p"]<SIG))
    k=sum(1 for _,c,_ in pairs if c); m=len(pairs)
    # binomial two-sided p vs 0.5
    from math import comb
    def binom_p(k,m):
        if m==0: return float("nan")
        tail=sum(comb(m,i) for i in range(m+1) if abs(i-m/2)>=abs(k-m/2)-1e-9)
        return min(1.0, tail/(2**m))
    return k,m,binom_p(k,m),pairs

def main():
    eur = load(os.path.join(BASE,"results/forward_local_edges.csv"),
               "ivw_b","ivw_p","nsnp","steiger_correct")
    eas = load(os.path.join(BASE,"results/network_eas_edges.csv"),
               "ivw_b","ivw_p","nsnp","steiger")
    nodes=[n for n in SHARED if any((n==a or n==b) for (a,b) in set(eur)&set(eas))]

    out=["=== Component D: cross-ancestry ORDER-divergence (graph-level, scale-free) ===",
         f"shared nodes: {nodes}\n"]

    # --- node-order Kendall-tau, main + sensitivities ---
    # Compute tau-a on the RAW netflow scores (ties excluded), NOT on tie-broken ordinal
    # ranks: an arbitrary tie-break lands near the min of admissible values and manufactures
    # a spurious 'NS until edges dropped' story. drop_SBP is reported only as a sensitivity
    # (dropping SBP is incoherent with a magnitude-invariant score; kept for transparency).
    out.append("[1] NODE-ORDER Kendall-tau_a on raw netflow scores (ties excluded; higher score=upstream):")
    rows=[]
    for label,dt,ds in [("main",False,False),("drop_->BMI",True,False),
                        ("drop_SBP(sens)",False,True),("drop_both(sens)",True,True)]:
        s_eur=netflow(eur,nodes,dt,ds); s_eas=netflow(eas,nodes,dt,ds)
        tau,nc,nd=kendall_tau(s_eur,s_eas)            # tie-aware
        pp=perm_p(s_eur,s_eas,tau)
        rows.append(("node_order_tau",label,tau,pp,len(nodes)))
        out.append(f"  {label:15s} tau_a={tau:+.3f}  divergence(1-tau)={1-tau:.3f}  perm_p={pp:.4f}  "
                   f"(concordant_pairs={nc}, discordant={nd})")
        if label=="main":
            out.append("     EUR netflow: "+", ".join(f"{k}:{s_eur[k]:+d}" for k in sorted(nodes,key=lambda k:-s_eur[k])))
            out.append("     EAS netflow: "+", ".join(f"{k}:{s_eas[k]:+d}" for k in sorted(nodes,key=lambda k:-s_eas[k])))
    out.append("     PRIMARY = 'main' (no drops); 'drop_->BMI' = pre-specified BBJ-selection artifact"
               " sensitivity; 'drop_SBP' shown only for transparency (incoherent with scale-free score).")

    # --- edge sign-concordance, main + sensitivities ---
    out.append("\n[2] EDGE SIGN-concordance on well-powered shared edges (EUR-sig anchor, both nSNP>=10):")
    for label,dt,ds in [("main",False,False),("drop_->BMI",True,False),
                        ("drop_SBP",False,True),("drop_both",True,True)]:
        k,m,bp,pairs=sign_concordance(eur,eas,dt,ds)
        both_sig=sum(1 for _,c,bs in pairs if c and bs)
        rows.append(("edge_sign_conc",label,(k/m if m else float('nan')),bp,m))
        out.append(f"  {label:12s} sign-concordant {k}/{m} = {(k/m if m else float('nan')):.3f}  "
                   f"binom_p={bp:.4g}  (both-sig&concordant={both_sig})")

    # discordant well-powered edges: split GENUINE (both-sig) from DARK (EAS underpowered/NS)
    k,m,bp,pairs=sign_concordance(eur,eas)
    genuine=[f"{a}->{b}(EUR{eur[(a,b)]['b']:+.2f}/EAS{eas[(a,b)]['b']:+.2f})"
             for (a,b),c,bs in pairs if not c and bs]
    dark   =[f"{a}->{b}" for (a,b),c,bs in pairs if not c and not bs]
    out.append(f"\n[3] Sign-DISCORDANT well-powered edges, split by EAS significance:")
    out.append(f"    GENUINE divergences (BOTH ancestries significant, opposite sign): {genuine if genuine else 'none'}")
    out.append(f"    DARK (EAS NS/underpowered — sign flip is noise, not divergence): {dark}")

    txt="\n".join(out)
    io.open(os.path.join(BASE,"results/order_divergence.txt"),"w",encoding="utf-8").write(txt)
    with io.open(os.path.join(BASE,"results/order_divergence.csv"),"w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh); w.writerow(["statistic","variant","value","p","k_or_nodes"])
        for r in rows: w.writerow([r[0],r[1],f"{r[2]:.4f}",f"{r[3]:.4g}",r[4]])
    print(txt); print("\n-> results/order_divergence.csv + .txt")

if __name__=="__main__": main()
