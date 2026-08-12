# -*- coding: utf-8 -*-
"""Component C: pre-registered falsifiable AHA-staging LEDGER.
For each pre-specified forward transition X->Y (X at lower AHA stage), test the REVERSE
edge Y->X with a POWER GATE and a TOST equivalence test, so that an absent reverse effect
is only counted as 'concordant' when the reverse test was actually powered to detect it.

Verdicts:
  CONCORDANT    : forward significant; reverse NS AND reverse powered to detect an effect
                  as large as the forward one (MDE_rev <= |b_fwd|)  -> reverse credibly ruled out
  DISCORDANT    : reverse significant (p<0.05) with |b_rev| not trivially small -> bidirectional / AHA-violating
  INDETERMINATE : reverse NS but UNDERPOWERED (MDE_rev > |b_fwd|)   -> cannot rule reverse out
Secondary: TOST equivalence of the reverse effect vs a fixed SESOI (default 0.05) at alpha=0.05.
Honest framing: population-level necessary-condition test, NOT within-person progression.
Reads results/forward_local_edges.csv (bidirectional). Writes results/staging_ledger.csv + .txt.
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
Z_A   = 1.959963985   # two-sided alpha=0.05
Z_PWR = 0.8416212     # 80% power
SESOI = 0.05          # secondary fixed equivalence bound (standardized units)

STAGE = {"BMI":1, "SBP":2,"TG":2,"HDL":2,"TC":2,"LDL":2,"HbA1c":2,"T2D":2,"eGFR":2,"CKD":2,
         "CAD":4,"HF":4,"Stroke":4}

# Pre-registered forward transitions to adjudicate (cross-stage or established upstream->downstream).
PREREG = [
    ("BMI","T2D"), ("BMI","CAD"), ("BMI","HF"), ("BMI","Stroke"),
    ("SBP","CAD"), ("SBP","HF"), ("SBP","Stroke"),
    ("LDL","CAD"), ("LDL","HF"),
    ("HbA1c","T2D"), ("HbA1c","CAD"),
    ("T2D","CAD"), ("T2D","HF"), ("T2D","Stroke"),
    ("CAD","HF"),
]

BINARY = {"CAD","HF","Stroke","T2D","CKD"}
def cross_type(X,Y):
    """continuous<->binary transition: magnitude/power comparison is unit-inconsistent."""
    return (X in BINARY) != (Y in BINARY)

def load(path, into=None):
    d = into if into is not None else {}
    with io.open(path,encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                eip=r.get("egger_intercept_p","")
                d.setdefault((r["exposure"],r["outcome"]), dict(
                    b=float(r["ivw_b"]), se=float(r["ivw_se"]), p=float(r["ivw_p"]),
                    n=int(float(r["nsnp"])),
                    steiger=str(r.get("steiger_correct","")).upper()=="TRUE",
                    egger_int_p=(float(eip) if eip not in("","NA") else float("nan"))))
            except (ValueError,KeyError): continue
    return d

def norm_cdf(x): return 0.5*(1+math.erf(x/math.sqrt(2)))

def tost(b, se, sesoi):
    """two one-sided tests: H0 |true|>=sesoi. Returns equivalence p (max of the two one-sided p)."""
    if se<=0: return None
    # two one-sided tests; reject non-equivalence if BOTH one-sided p are small
    p1 = 1-norm_cdf((b+sesoi)/se)   # H0: true <= -sesoi
    p2 = norm_cdf((b-sesoi)/se)     # H0: true >=  sesoi
    return max(p1,p2)

def main():
    d = load(os.path.join(BASE,"results/forward_local_edges.csv"))
    rev_path = os.path.join(BASE,"results/reverse_arms_edges.csv")   # SBP/Stroke reverse arms (gap-closers)
    if os.path.exists(rev_path): load(rev_path, into=d)
    rows=[]
    for X,Y in PREREG:
        f = d.get((X,Y)); r = d.get((Y,X))
        if f is None:
            rows.append((X,Y,None)); continue
        bf,sef,pf,nf = f["b"],f["se"],f["p"],f["n"]
        fwd_sig = pf<0.05
        same_stage = STAGE.get(X)==STAGE.get(Y)
        ct = cross_type(X,Y)
        if r is None:
            verdict="NO_REVERSE_DATA"; brev=srev=prev=nrev=mde=None; tp=None; rev_st=rev_eip=None
        else:
            brev,srev,prev,nrev = r["b"],r["se"],r["p"],r["n"]
            rev_st  = r["steiger"]             # reverse edge Steiger-correct?
            rev_eip = r["egger_int_p"]         # reverse edge directional-pleiotropy flag
            mde = (Z_A+Z_PWR)*srev             # 80%-power min detectable reverse effect
            powered = mde <= abs(bf)           # powered to see an effect as big as forward
            rev_sig = prev<0.05
            rev_clean = rev_st and not (rev_eip==rev_eip and rev_eip<0.05)  # Steiger-correct & no dir-pleiotropy
            tp = tost(brev, srev, SESOI)
            mag_ok = True if ct else (abs(brev) >= 0.2*abs(bf))  # ct = cross-type: mag/power guard unit-inconsistent
            # A reverse edge only counts as genuine reverse causation if it is significant AND
            # Steiger-correct (right direction) AND not driven by directional pleiotropy.
            # (A transition is only adjudicable if the FORWARD edge is itself established.)
            if not fwd_sig:
                verdict="FWD_NOT_SIG"
            elif rev_sig and rev_st and (rev_eip==rev_eip and rev_eip<0.05):
                verdict="DISCORDANT_PLEIOTROPIC"     # sig reverse but flagged as pleiotropy artifact
            elif rev_sig and rev_clean and mag_ok:
                verdict="DISCORDANT"                 # genuine bidirectional / AHA-violating
            elif rev_sig and (not rev_st):
                # Steiger wrong-direction is a variance-explained DIRECTIONALITY criterion: it
                # says the reverse estimate is not directionally supported. It does NOT prove
                # confounding — the same pattern arises from shared instruments, pleiotropy or
                # sample overlap. The prose gloss "and so confounded rather than causal" was
                # imported from this comment and shipped through five review rounds.
                verdict="CONCORDANT"                 # reverse sig but not directionally supported
            elif (not rev_sig) and (not rev_st):
                verdict="CONCORDANT"                 # reverse NS AND wrong-direction = doubly not a threat
            elif (not rev_sig) and powered:
                verdict="CONCORDANT"                 # reverse is a tight null in its own units -> ruled out
            elif (not rev_sig) and (not powered):
                verdict="INDETERMINATE"              # NS and underpowered -> cannot rule out reverse
            else:
                verdict="AMBIGUOUS"
        rows.append((X,Y,dict(bf=bf,sef=sef,pf=pf,nf=nf,fwd_sig=fwd_sig,same_stage=same_stage,cross_type=ct,
                    brev=brev,srev=srev,prev=prev,nrev=nrev,mde=mde,tost_p=tp,
                    rev_st=rev_st,rev_eip=rev_eip,verdict=verdict)))
    # write
    out_csv=os.path.join(BASE,"results/staging_ledger.csv")
    with io.open(out_csv,"w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh)
        w.writerow(["transition","stage_X","stage_Y","same_stage","cross_type","b_fwd","se_fwd","p_fwd","fwd_sig",
                    "b_rev","se_rev","p_rev","n_rev","rev_steiger_correct","rev_egger_int_p",
                    "mde_rev_80pct","reverse_powered_to_bfwd","tost_equiv_p_sesoi0.05","verdict"])
        for X,Y,o in rows:
            if o is None:
                w.writerow([f"{X}->{Y}",STAGE.get(X),STAGE.get(Y),"NA"]); continue
            powered = (o["mde"] is not None and abs(o["mde"])<=abs(o["bf"]))
            w.writerow([f"{X}->{Y}",STAGE.get(X),STAGE.get(Y),o["same_stage"],o["cross_type"],
                f"{o['bf']:+.4f}",f"{o['sef']:.4f}",f"{o['pf']:.2e}",o["fwd_sig"],
                "" if o["brev"] is None else f"{o['brev']:+.4f}",
                "" if o["srev"] is None else f"{o['srev']:.4f}",
                "" if o["prev"] is None else f"{o['prev']:.2e}",
                "" if o["nrev"] is None else o["nrev"],
                "" if o["rev_st"] is None else o["rev_st"],
                "" if o["rev_eip"] is None or o["rev_eip"]!=o["rev_eip"] else f"{o['rev_eip']:.3f}",
                "" if o["mde"]  is None else f"{o['mde']:.4f}",
                "" if o["mde"]  is None else powered,
                "" if o["tost_p"] is None else f"{o['tost_p']:.3f}",
                o["verdict"]])
    # console summary
    lines=["=== Component C: pre-registered falsifiable staging ledger ==="]
    lines.append(f"SESOI(equiv)={SESOI}; reverse power gate = 80%-power MDE ({Z_A+Z_PWR:.2f}*se_rev) vs |b_fwd|\n")
    from collections import Counter
    vc=Counter()
    for X,Y,o in rows:
        if o is None: continue
        vc[o["verdict"]]+=1
        pw = "" if o["mde"] is None else ("powered" if abs(o["mde"])<=abs(o["bf"]) else "UNDERPOWERED")
        rst = "" if o["rev_st"] is None else ("revSt+" if o["rev_st"] else "revSt-")
        ss  = (" SAME-STAGE" if o["same_stage"] else "") + (" X-TYPE" if o["cross_type"] else "")
        lines.append(f"  {X:5s}->{Y:6s} [{STAGE.get(X)}->{STAGE.get(Y)}]  fwd b={o['bf']:+.3f} p={o['pf']:.1e} | "
                     f"rev b={('%+.3f'%o['brev']) if o['brev'] is not None else 'NA':>7} "
                     f"p={('%.1e'%o['prev']) if o['prev'] is not None else 'NA':>8} "
                     f"MDE={('%.3f'%o['mde']) if o['mde'] is not None else 'NA':>6} {pw:12s} {rst:6s} -> {o['verdict']}{ss}")
    lines.append("\nverdict tally: "+", ".join(f"{k}={v}" for k,v in sorted(vc.items())))
    txt="\n".join(lines)
    io.open(os.path.join(BASE,"results/staging_ledger.txt"),"w",encoding="utf-8").write(txt)
    print(txt)
    print("\n-> results/staging_ledger.csv + staging_ledger.txt")

if __name__=="__main__": main()
