# -*- coding: utf-8 -*-
"""Zip the frozen submission package and verify the archive against its own checksum manifest.

A zip is only useful if what comes out equals what went in, so this does not just compress: it
re-reads every entry back out of the finished archive and md5-compares it to MANIFEST_checksums.txt.
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

import hashlib, os, re, sys, zipfile

ROOT = P4_ROOT
# Version-parameterised so the script cannot zip a different version than the one just cut.
# Must match PKG_VERSION in build_submission_v1.py. Override with:  python zip_submission_v1.py v1
_bsv = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_submission_v1.py"),
            encoding="utf-8").read()
_mv = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _mv, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _mv.group(1)   # DERIVED, never typed: a gate pointed at a stale package is worse than none
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")
TOP = f"CKM_Paper4_submission_package_{PKG_VERSION}"   # no spaces: safer across upload portals
ZIP = os.path.join(ROOT, TOP + ".zip")
if not os.path.isdir(PKG):
    print("ABORT: no such package directory:", PKG); sys.exit(1)

def md5_bytes(b):
    h = hashlib.md5(); h.update(b); return h.hexdigest()

def md5_file(p, chunk=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

# ---- enumerate the source of truth ----
files = []
for r, _, fs in os.walk(PKG):
    for f in fs:
        p = os.path.join(r, f)
        files.append((os.path.relpath(p, PKG).replace("\\", "/"), p))
files.sort()
print(f"packaging {len(files)} files from {PKG}")

# ---- write ----
if os.path.exists(ZIP):
    os.remove(ZIP)
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for rel, p in files:
        z.write(p, arcname=f"{TOP}/{rel}")

# ---- verify by reading the archive back ----
problems = []
with zipfile.ZipFile(ZIP) as z:
    bad = z.testzip()
    if bad:
        problems.append(f"corrupt entry: {bad}")
    names = [n for n in z.namelist() if not n.endswith("/")]
    if len(names) != len(files):
        problems.append(f"archive holds {len(names)} entries, package has {len(files)} files")
    on_disk = {rel: md5_file(p) for rel, p in files}
    for rel in on_disk:
        try:
            got = md5_bytes(z.read(f"{TOP}/{rel}"))
        except KeyError:
            problems.append(f"missing from archive: {rel}"); continue
        if got != on_disk[rel]:
            problems.append(f"md5 differs after round-trip: {rel}")

    # cross-check against the package's own manifest, so the zip is verified against the
    # checksums the recipient will use, not only against the bytes we just wrote.
    man = z.read(f"{TOP}/MANIFEST_checksums.txt").decode("utf-8")
    listed = {m.group(3): m.group(1) for m in
              (re.match(r"^([0-9a-f]{32})\s+(\d+)\s+(.+)$", ln.strip()) for ln in man.split("\n"))
              if m}
    for rel, h in listed.items():
        if rel not in on_disk:
            problems.append(f"manifest lists a file absent from the archive: {rel}")
        elif on_disk[rel] != h:
            problems.append(f"manifest md5 disagrees with archived bytes: {rel}")
    unlisted = set(on_disk) - set(listed) - {"MANIFEST_checksums.txt"}
    if unlisted:
        problems.append(f"archived but not in manifest: {sorted(unlisted)}")

size = os.path.getsize(ZIP)
raw = sum(os.path.getsize(p) for _, p in files)
print(f"wrote {os.path.relpath(ZIP, ROOT)}")
print(f"  {len(files)} files, {size/1e6:.2f} MB compressed from {raw/1e6:.2f} MB "
      f"({100*size/raw:.0f}%)")
print(f"  manifest cross-checked: {len(listed)} entries")

if problems:
    print(f"\nZIP VERIFICATION FAILED — {len(problems)} problem(s):")
    for p_ in problems:
        print("  x", p_)
    sys.exit(1)
print("\nverified: every file round-trips md5-identical and matches MANIFEST_checksums.txt")
