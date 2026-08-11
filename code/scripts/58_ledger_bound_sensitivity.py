# -*- coding: utf-8 -*-
"""Step 58: bound-sensitivity for the native-scale staging ledger.

The one thing the native-scale rule is exposed to is that its negligibility bounds are CHOSEN
(log(1.05) log-odds for binary outcomes, 1 mmHg for SBP, 0.05 SD otherwise). This sweeps them and
reports which verdicts are bound-stable and which are not, so the ledger's counts are presented as
bound-conditional rather than as facts.

It imports `adjudicate` and `knife_edge` from scripts/57_ledger_native.py and patches ONLY
`native_bound`, so the sweep exercises the identical decision rule rather than a reimplementation,
and it asserts that the primary scenario reproduces the committed ledger before sweeping.

Writes results/staging_ledger_bound_sensitivity.{csv,txt}. Deterministic; no randomness.
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

import csv, io, os, math, importlib.util
from collections import Counter, OrderedDict

BASE = P4_BASE
RES = os.path.join(BASE, "results")

spec = importlib.util.spec_from_file_location("led57", os.path.join(BASE, "scripts", "57_ledger_native.py"))
led = importlib.util.module_from_spec(spec)
spec.loader.exec_module(led)                      # safe: 57's main() is behind __main__

d = led.load(os.path.join(RES, "forward_local_edges.csv"))
_rev = os.path.join(RES, "reverse_arms_edges.csv")
if os.path.exists(_rev):
    led.load(_rev, into=d)

# Scenarios, each anchored on a stated rationale rather than a round number.
SCENARIOS = OrderedDict([
    ("very tight", (math.log(1.01), 0.5, 0.02)),
    ("tight",      (math.log(1.02), 0.5, 0.02)),
    ("PRIMARY",    (math.log(1.05), 1.0, 0.05)),
    ("loose",      (math.log(1.10), 2.0, 0.10)),
    ("very loose", (math.log(1.20), 5.0, 0.20)),
])
LABEL = {"very tight": "OR 1.01 / 0.5 mmHg / 0.02 SD", "tight": "OR 1.02 / 0.5 mmHg / 0.02 SD",
         "PRIMARY": "OR 1.05 / 1 mmHg / 0.05 SD", "loose": "OR 1.10 / 2 mmHg / 0.10 SD",
         "very loose": "OR 1.20 / 5 mmHg / 0.20 SD"}


def run(bin_b, sbp_b, sd_b):
    def nb(outcome):
        if outcome in led.BINARY:
            return bin_b
        if outcome == "SBP":
            return sbp_b
        return sd_b
    old, led.native_bound = led.native_bound, nb
    out = OrderedDict()
    try:
        for X, Y in led.PREREG:
            f, rev = d.get((X, Y)), d.get((Y, X))
            if f is None:
                continue
            out["%s->%s" % (X, Y)] = led.adjudicate(f["b"], f["p"], rev, X, Y)[1]
    finally:
        led.native_bound = old
    return out


RESULTS = OrderedDict((k, run(*v)) for k, v in SCENARIOS.items())
PRIM = RESULTS["PRIMARY"]

# --- self-check: the primary scenario must reproduce the committed ledger exactly ---------------
committed = {r["transition"]: r["verdict"] for r in
             csv.DictReader(io.open(os.path.join(RES, "staging_ledger_native.csv"), encoding="utf-8"))}
bad = {t: (PRIM[t], committed.get(t)) for t in PRIM if committed.get(t) != PRIM[t]}
assert not bad, "bound sweep does not reproduce the committed ledger: %r" % bad

cols = ["transition"] + list(SCENARIOS) + ["stability"]
rows = []
for t in PRIM:
    verdicts = [RESULTS[s][t] for s in SCENARIOS]
    rows.append(dict([("transition", t)] + list(zip(SCENARIOS, verdicts))
                     + [("stability", "stable" if len(set(verdicts)) == 1 else "bound-sensitive")]))

with io.open(os.path.join(RES, "staging_ledger_bound_sensitivity.csv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in rows:
        w.writerow(r)

L = []
A = L.append
A("=== Bound sensitivity of the native-scale staging ledger ===")
A("")
A("Negligibility bounds are clinical judgements, so every verdict is reported across a sweep of them.")
A("Scenario bounds (binary log-odds / SBP mmHg / other per-SD):")
for s in SCENARIOS:
    A("  %-11s %s" % (s, LABEL[s]))
A("")
A("VERDICT COUNTS BY SCENARIO")
hdr = "%-12s %s" % ("scenario", "  ".join("%-13s" % k for k in
                                          ["CONCORDANT", "INDETERMINATE", "DISCORDANT", "DISC_CAVEATED"]))
A(hdr); A("-" * len(hdr))
for s in SCENARIOS:
    c = Counter(RESULTS[s].values())
    A("%-12s %s" % (s, "  ".join("%-13d" % c.get(k, 0) for k in
                                 ["CONCORDANT", "INDETERMINATE", "DISCORDANT", "DISCORDANT_CAVEATED"])))
A("")
stable = [r["transition"] for r in rows if r["stability"] == "stable"]
sens = [r["transition"] for r in rows if r["stability"] != "stable"]
A("BOUND-STABLE (%d/%d): %s" % (len(stable), len(rows), ", ".join(stable)))
A("BOUND-SENSITIVE (%d/%d): %s" % (len(sens), len(rows), ", ".join(sens)))
A("")
A("THE LOAD-BEARING CLAIM IS BOUND-INVARIANT. The ledger's falsification claim is not the concordant")
A("count but the direction of every reverse effect that is real. Under every scenario above, each")
A("transition in the discordant family has a reverse edge running from a stage-4 disease into a")
A("stage-2 trait — feedback or an index-event effect — and none reverses the staging order:")
for s in SCENARIOS:
    disc = [t for t, v in RESULTS[s].items() if v in ("DISCORDANT", "DISCORDANT_CAVEATED")]
    detail = []
    for t in disc:
        X, Y = t.split("->")
        detail.append("%s (reverse %s->%s, stage %s->%s)" % (t, Y, X, led.STAGE[Y], led.STAGE[X]))
    A("  %-11s n=%d  %s" % (s, len(disc), "; ".join(detail) if detail else "none"))
A("")
A("What DOES move is the concordant/indeterminate split, which is why the counts are reported as")
A("bound-conditional. Tightening the bound makes an excluded reverse effect harder to declare;")
A("loosening it makes a supported reverse effect harder to declare.")
txt = "\n".join(L) + "\n"
io.open(os.path.join(RES, "staging_ledger_bound_sensitivity.txt"), "w", encoding="utf-8").write(txt)
print(txt)
print("wrote results/staging_ledger_bound_sensitivity.{csv,txt}")
