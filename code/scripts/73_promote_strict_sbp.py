# -*- coding: utf-8 -*-
"""Promote the strict allele-compatible X->SBP estimates to PRIMARY (round-3 blocker 2).

The round-3 review's objection is procedural and correct: the Methods described a position-only
primary extraction for the disease->SBP edges, with allele compatibility enforced only later at
harmonisation. Harmonisation cannot repair a row that was selected wrongly at extraction, so the
strict set has to be the primary analysis rather than a sensitivity that sits beside it.

What this script does NOT do is re-derive the estimates: script 62 already audited the matching and
script 63 already re-estimated every X->SBP edge on the strict instrument set. This step swaps those
estimates into the canonical network table, keeps the superseded values beside them for the record,
and prints EVERY consequence -- edge counts, significance flips, staging membership -- so the
propagation is visible rather than asserted.

Backup: results/forward_local_edges.pre_strict_sbp.csv (written once, never overwritten).
Out:    results/forward_local_edges.csv (X->SBP rows replaced)
        results/strict_sbp_promotion.txt (the change log)
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

import io, os, csv, sys, shutil

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE
RES = os.path.join(BASE, "results")
NET = os.path.join(RES, "forward_local_edges.csv")
BAK = os.path.join(RES, "forward_local_edges.pre_strict_sbp.csv")
STRICT = os.path.join(RES, "sbp_strict_mr.csv")

if not os.path.exists(BAK):
    shutil.copy2(NET, BAK)
    print(f"backup written: {os.path.basename(BAK)}")
else:
    print(f"backup already exists ({os.path.basename(BAK)}); re-promoting from it")
    shutil.copy2(BAK, NET)          # idempotent: always promote from the pristine pre-strict table

rows = list(csv.DictReader(open(NET, encoding="utf-8")))
fields = list(rows[0].keys())
strict = {(r["exposure"], r["outcome"]): r for r in csv.DictReader(open(STRICT, encoding="utf-8"))}

BONF = 0.05 / len(rows)


def is_sig(p):
    return p not in (None, "") and float(p) < BONF


def steiger_ok(r):
    return str(r.get("steiger_correct", "")).strip().upper() in ("TRUE", "1", "T")


before_sig = {(r["exposure"], r["outcome"]) for r in rows if is_sig(r["ivw_p"])}
before_graph = {k for k in before_sig if steiger_ok(next(r for r in rows if (r["exposure"], r["outcome"]) == k))}

log = ["=== X->SBP promoted to strict allele-compatible matching (PRIMARY) ===", "",
       "Round-3 blocker 2. Estimates come from scripts/63_sbp_strict_mr.R on the step-62 strict set",
       "('match' + 'strand_flip' only; every strand-ambiguous variant dropped). Nothing is re-derived",
       f"here. Network Bonferroni = 0.05/{len(rows)} = {BONF:.3g}.", "",
       f"{'edge':16s} {'nSNP was':>9s} {'nSNP now':>9s} {'b was':>10s} {'b now':>10s} "
       f"{'P was':>11s} {'P now':>11s}  change"]

# Strict-set pleiotropy statistics (step 75). Loaded here so the promoted row describes ONE
# instrument set end to end; without it the Egger columns stay blank and a NaN intercept silently
# reclassifies a reverse edge downstream.
_PL = os.path.join(RES, "sbp_strict_pleiotropy.csv")
if not os.path.exists(_PL):
    sys.exit("ABORT: results/sbp_strict_pleiotropy.csv missing — run "
             "`Rscript scripts/75_sbp_strict_pleiotropy.R` first. Promoting without it leaves the "
             "Egger/WM/Q columns blank, which is how CAD->SBP was silently misclassified.")
STRICT_PLEIO = {r["exposure"]: r for r in csv.DictReader(open(_PL, encoding="utf-8"))}
missing_pleio = []

changed = []
for r in rows:
    key = (r["exposure"], r["outcome"])
    if r["outcome"] != "SBP" or key not in strict:
        continue
    s = strict[key]
    was_n, was_b, was_p = r["nsnp"], float(r["ivw_b"]), float(r["ivw_p"])
    now_n, now_b = s["strict_nsnp"], float(s["strict_b"])
    now_se, now_p = float(s["strict_se"]), float(s["strict_p"])
    flip = ("SIGNIFICANCE LOST" if is_sig(was_p) and not is_sig(now_p) else
            "SIGNIFICANCE GAINED" if not is_sig(was_p) and is_sig(now_p) else
            "sign flip" if was_b * now_b < 0 else "status unchanged")
    if flip != "status unchanged":
        changed.append((f"{key[0]}->SBP", flip))
    log.append(f"{key[0]+'->SBP':16s} {was_n:>9s} {now_n:>9s} {was_b:+10.4f} {now_b:+10.4f} "
               f"{was_p:11.3g} {now_p:11.3g}  {flip}")
    # replace the primary estimate; the superseded value stays in the backup and in Table S-strict
    r["nsnp"], r["ivw_b"], r["ivw_se"], r["ivw_p"] = now_n, f"{now_b}", f"{now_se}", f"{now_p}"
    # The Egger/WM/Q columns were computed on the base set and no longer describe the primary
    # instrument set. Blanking them was honest but INCOMPLETE: nothing refilled them, and a
    # downstream classifier read the resulting NaN intercept as "no pleiotropy" and silently
    # downgraded CAD->SBP (round-4 Tier 1). They are now REPLACED with the strict-set values
    # computed by scripts/75_sbp_strict_pleiotropy.R, so the whole row describes one instrument set.
    s = STRICT_PLEIO.get(key[0])
    if s is None:
        for col in ("egger_b", "egger_p", "egger_intercept", "egger_intercept_p",
                    "wm_b", "wm_p", "Q", "Q_p"):
            if col in r:
                r[col] = ""
        missing_pleio.append(key[0])
    else:
        for col, src in (("egger_b", "egger_b"), ("egger_p", "egger_p"),
                         ("egger_intercept", "egger_intercept"),
                         ("egger_intercept_p", "egger_intercept_p"),
                         ("wm_b", "wm_b"), ("wm_p", "wm_p"), ("Q", "Q"), ("Q_p", "Q_p")):
            if col in r:
                r[col] = s[src]

with open(NET, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, quoting=csv.QUOTE_NONNUMERIC)
    w.writeheader()
    w.writerows(rows)

after_sig = {(r["exposure"], r["outcome"]) for r in rows if is_sig(r["ivw_p"])}
after_graph = {k for k in after_sig
               if steiger_ok(next(r for r in rows if (r["exposure"], r["outcome"]) == k))}

log += ["",
        f"Bonferroni-significant edges : {len(before_sig)} -> {len(after_sig)}",
        f"Steiger-directed graph edges : {len(before_graph)} -> {len(after_graph)}",
        f"  lost   : {sorted('->'.join(k) for k in before_graph - after_graph) or 'none'}",
        f"  gained : {sorted('->'.join(k) for k in after_graph - before_graph) or 'none'}",
        ""]
if missing_pleio:
    log.append("")
    log.append("X->SBP rows with NO strict-set pleiotropy statistics (Egger/WM/Q left blank): "
               + ", ".join(sorted(missing_pleio)))
if changed:
    log.append("Edges whose significance status changed:")
    for e, f_ in changed:
        log.append(f"  {e}: {f_}")
else:
    log.append("No edge changed significance status.")
log += ["",
        "The Egger, weighted-median and Cochran-Q columns for the X->SBP rows are REPLACED with the",
        "strict-set values from scripts/75_sbp_strict_pleiotropy.R, not carried over from the base",
        "set and not left blank. Carrying them would mix two instrument sets in one row; leaving",
        "them blank let a NaN Egger intercept read downstream as \"no pleiotropy detected\", which",
        "silently downgraded CAD->SBP in the native-scale ledger (round-4 Tier 1)."]

open(os.path.join(RES, "strict_sbp_promotion.txt"), "w", encoding="utf-8").write("\n".join(log) + "\n")
print("\n".join(log))
print(f"\nwrote {os.path.relpath(NET, BASE)} and results/strict_sbp_promotion.txt")
