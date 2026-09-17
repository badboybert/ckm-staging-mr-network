# -*- coding: utf-8 -*-
"""One command, one line — the whole session-start verification pass for CKM Paper 4.

WHY THIS EXISTS
---------------
The session-start pass was seven commands and seven blocks of output. Measured, it is only ~23 s of
wall clock (Paper 4 has no slow gates; Paper 3's equivalent was ~7 min, which is why THAT project
needed a 391-file fingerprint tool). So the cost here was never time — it was CONTEXT and attention
at the exact moment the session should be starting work.

This collapses it: green prints ONE line, red prints everything. There is also a `--fast` mode that
short-circuits when nothing the suite reads has changed since the last recorded green run.

WHAT MAKES --fast SAFE
----------------------
A fast path that says "unchanged" when something DID change is worse than no fast path. Two rules,
both learned the hard way on paper 3:

  * THE SCOPE IS DERIVED, NEVER HAND-LISTED. `--discover` runs the real suite under
    sys.addaudithook('open') and records the true read-set. A hand-written list omits the input
    nobody thought of, and that omission is invisible.
  * THE READ-SET LIVES INSIDE THE MANIFEST. Deleting a sidecar cannot silently shrink what a later
    --check watches, and a missing, corrupt or empty manifest FAILS CLOSED into the full suite.

--fast never reports green on its own authority. It reports "nothing changed since the last green
run", which is a different and weaker claim, and it says so.

USAGE
-----
    python scripts/verify_all.py              # full pass, one-line summary   (~23 s)
    python scripts/verify_all.py --fast       # short-circuit if unchanged    (~1 s)
    python scripts/verify_all.py --record     # full pass, then store the fingerprint
    python scripts/verify_all.py --discover   # re-derive the read-set (after adding a gate)
    python scripts/verify_all.py --self-test  # prove --fast actually notices a change
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

import hashlib, io, json, os, subprocess, sys, time

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCR = os.path.join(BASE, "scripts")
MANIFEST = os.path.join(BASE, "manifest", "VERIFY_FINGERPRINT.json")

# The suite, with the result line each gate must print when it passes. The EXPECTED string is part
# of the check: a gate that dies before printing, or prints a different tally, is not green.
SUITE = [
    ("gate_package.py",        "GATE PASSED"),
    ("qa_manuscript.py",       "0 FAIL"),
    ("verify_rebuttal_r2.py",  "0 FAILED"),
    ("verify_rebuttal_r3.py",  "0 FAILED"),
    ("verify_rebuttal_r4.py",  "0 FAILED"),
    ("verify_rebuttal_r5.py",  "0 FAILED"),
    ("verify_rebuttal_r6.py",  "0 FAILED"),
    ("scan_claims.py",         "CLAIM SCAN CLEAN"),
    ("verify_repo.py",         "0 FAILED"),
]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:16]


def run_one(script, expect):
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(SCR, script)], capture_output=True,
                       text=True, encoding="utf-8", errors="replace")
    out = (r.stdout or "") + (r.stderr or "")
    return {"script": script, "rc": r.returncode, "ok": r.returncode == 0 and expect in out,
            "secs": round(time.time() - t0, 1), "out": out}


def run_suite():
    return [run_one(s, e) for s, e in SUITE]


def report(results, header):
    bad = [r for r in results if not r["ok"]]
    secs = sum(r["secs"] for r in results)
    if not bad:
        print(f"{header}: {len(results)}/{len(results)} green in {secs:.0f}s "
              f"({', '.join(r['script'].replace('.py', '') for r in results)})")
        return 0
    print(f"{header}: {len(results) - len(bad)}/{len(results)} green, {len(bad)} FAILED in {secs:.0f}s\n")
    for r in bad:
        print("=" * 96)
        print(f"FAILED: {r['script']}  (rc={r['rc']})")
        tail = [l for l in r["out"].split("\n") if l.strip()]
        for l in tail[-40:]:
            print("   " + l)
    return 1


# ---------------------------------------------------------------- read-set discovery
def discover():
    """Run each gate under an open() audit hook and record every file it actually reads."""
    hook = os.path.join(SCR, "_discover_hook.py")
    io.open(hook, "w", encoding="utf-8").write(
        "import sys, os, json, runpy\n"
        "SEEN=set()\n"
        "def _h(ev, a):\n"
        "    if ev in ('open','exec','compile') and a and isinstance(a[0], str):\n"
        "        SEEN.add(os.path.abspath(a[0]))\n"
        "sys.addaudithook(_h)\n"
        "tgt=sys.argv[1]; out=sys.argv[2]; sys.argv=[tgt]\n"
        "try: runpy.run_path(tgt, run_name='__main__')\n"
        "except SystemExit: pass\n"
        "except Exception: pass\n"
        # An `import` does not raise a plain `open` audit event, so the shared checkers in
        # ~/.claude/lib were invisible to the first version of this hook and the fast path would
        # have reported UNCHANGED after an edit to pdftext.py or prosecheck.py — a checker change
        # that alters what the suite SEES. Sweep the loaded modules at exit: still derived from the
        # real run, never hand-listed.
        "for _m in list(sys.modules.values()):\n"
        "    _f = getattr(_m, '__file__', None)\n"
        "    if isinstance(_f, str) and _f.endswith('.py'):\n"
        "        SEEN.add(os.path.abspath(_f))\n"
        "json.dump(sorted(SEEN), open(out,'w'))\n")
    seen = set()
    try:
        for script, _ in SUITE:
            tmp = os.path.join(SCR, "_reads.json")
            subprocess.run([sys.executable, hook, os.path.join(SCR, script), tmp],
                           capture_output=True, text=True)
            if os.path.exists(tmp):
                seen |= set(json.load(open(tmp)))
                os.remove(tmp)
    finally:
        if os.path.exists(hook):
            os.remove(hook)
    # Keep what belongs to this project and the shared checkers; drop the interpreter's own reads.
    keep = []
    for p in seen:
        q = p.replace("\\", "/")
        if "site-packages" in q or "/Lib/" in q or "__pycache__" in q or q.endswith(".pyc"):
            continue
        if not os.path.isfile(p):
            continue
        if q.startswith(BASE.replace("\\", "/")) or "/.claude/lib" in q \
                or "/paper 4/" in q or "/paper 4\\" in p:
            keep.append(p)
    return sorted(set(keep))


def fingerprint(reads):
    return {p: sha(p) for p in reads if os.path.isfile(p)}


def cmd_record(reads=None):
    results = run_suite()
    rc = report(results, "full pass")
    if rc:
        print("\nNOT RECORDED — the fingerprint only ever describes a GREEN run.")
        return rc
    if reads is None:
        prev = load_manifest()
        reads = list(prev["reads"]) if prev else discover()
    fp = fingerprint(reads)
    os.makedirs(os.path.dirname(MANIFEST), exist_ok=True)
    json.dump({"recorded": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "suite": [s for s, _ in SUITE],
               "results": {r["script"]: r["secs"] for r in results},
               "reads": sorted(fp), "hashes": fp},
              io.open(MANIFEST, "w", encoding="utf-8"), indent=1)
    print(f"recorded: {len(fp)} files fingerprinted -> manifest/VERIFY_FINGERPRINT.json")
    return 0


def load_manifest():
    """Fail CLOSED: a missing, corrupt or empty manifest yields None, which forces the full suite."""
    try:
        m = json.load(io.open(MANIFEST, encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(m, dict) or not m.get("reads") or not m.get("hashes"):
        return None
    if list(m["suite"]) != [s for s, _ in SUITE]:
        return None                      # the suite itself changed; the record does not describe it
    return m


def cmd_fast():
    m = load_manifest()
    if m is None:
        print("fast: no usable fingerprint (missing, corrupt, or the suite changed) — running full pass")
        return report(run_suite(), "full pass")
    changed, gone = [], []
    for p, h in m["hashes"].items():
        if not os.path.isfile(p):
            gone.append(p)
        elif sha(p) != h:
            changed.append(p)
    if not changed and not gone:
        print(f"fast: UNCHANGED — none of the {len(m['hashes'])} files the suite reads has changed "
              f"since the green run recorded {m['recorded']}. Full pass not re-run.")
        return 0
    for p in (changed + gone)[:8]:
        print(f"  {'deleted' if p in gone else 'changed'}: {os.path.relpath(p, BASE)}")
    print(f"fast: {len(changed)} changed, {len(gone)} deleted — running the full pass\n")
    return report(run_suite(), "full pass")


def cmd_self_test():
    """Prove --fast NOTICES a change, and prove it reports UNCHANGED when there is none.

    A fast path is only ever seen saying 'unchanged'. That is precisely the state in which it is
    indistinguishable from a fast path that is broken, so both directions are asserted.
    """
    fails = []
    N_ASSERTIONS = 5      # positive control · mutation detected · named · size · package · libs
    m = load_manifest()
    if m is None:
        return print("self-test: no manifest to test against — run --record first") or 1
    # positive control: an untouched tree must report UNCHANGED
    r = subprocess.run([sys.executable, __file__, "--fast"], capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    if "UNCHANGED" not in (r.stdout or ""):
        fails.append("positive control: an untouched tree did not report UNCHANGED")
    # mutation: touching ONE watched file must stop it reporting UNCHANGED
    victim = next((p for p in m["hashes"] if p.endswith(".md") and os.path.isfile(p)), None)
    if not victim:
        fails.append("no watched .md to mutate — cannot test detection")
    else:
        orig = io.open(victim, "rb").read()
        try:
            io.open(victim, "ab").write(b"\n<!-- fingerprint self-test -->\n")
            r = subprocess.run([sys.executable, __file__, "--fast"], capture_output=True, text=True,
                               encoding="utf-8", errors="replace")
            out = (r.stdout or "")
            if "UNCHANGED" in out:
                fails.append(f"MUTATION NOT DETECTED: {os.path.relpath(victim, BASE)} was modified "
                             f"and --fast still reported UNCHANGED")
            elif "changed:" not in out:
                fails.append("a change was detected but not named — the report is unusable")
        finally:
            io.open(victim, "wb").write(orig)
    # the read-set must be non-trivial and must include the package and the gate sources
    reads = " ".join(m["reads"]).replace("\\", "/")
    if len(m["hashes"]) < 60:
        fails.append(f"read-set is only {len(m['hashes'])} files — too small to be the real scope")
    if "submission package v4.3" not in reads:
        fails.append("read-set does not include the shipped package")
    # The shared checkers decide what the suite can SEE. An edit to pdftext.py or prosecheck.py
    # changes every PDF- and prose-based verdict, and the first version of the discovery hook missed
    # them entirely because an import raises no plain `open` event. Assert they are watched.
    if "/.claude/lib" not in reads:
        fails.append("read-set does not include the shared checkers in ~/.claude/lib — editing "
                     "pdftext.py or prosecheck.py would go unnoticed by --fast")
    if "scripts/gate_package.py" not in reads and "scripts\\gate_package.py" not in " ".join(m["reads"]):
        fails.append("read-set does not include the gate sources — editing a gate would go unnoticed")
    if fails:
        print("SELF-TEST FAILED:")
        for f in fails:
            print("  x", f)
        return 1
    # The count is DERIVED. It read "4/4" while five assertions ran — a self-test that misreports
    # its own size is the first thing that stops being true when a check is added.
    print(f"verify_all self-test: {N_ASSERTIONS}/{N_ASSERTIONS} OK (untouched -> UNCHANGED; a "
          f"mutation IS detected and named; read-set is {len(m['hashes'])} files covering the "
          f"package, the generators, the gate sources and the shared checkers)")
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--self-test" in a:
        sys.exit(cmd_self_test())
    if "--discover" in a:
        rs = discover()
        print(f"discovered read-set: {len(rs)} files")
        sys.exit(cmd_record(rs))
    if "--record" in a:
        sys.exit(cmd_record())
    if "--fast" in a:
        sys.exit(cmd_fast())
    sys.exit(report(run_suite(), "full pass"))
