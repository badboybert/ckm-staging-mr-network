# -*- coding: utf-8 -*-
"""Verify the public deposit at `paper 4/repo/` — independently of the script that built it.

The build script's own assertions prove that IT thinks the deposit is sound. This one reads the
deposit as an outsider would: nothing is imported from build_repo.py, every expectation is
re-derived from the shipped submission package or from the manuscript's canonical facts, and the
load-bearing check is not a string scan but an actual RENDER — a figure is regenerated from the
deposit alone, with the working tree's absolute paths gone, and its output compared to the shipped
figure. A deposit that passes a grep and cannot run is the failure mode this exists to catch.

Run:  PYTHONIOENCODING=utf-8 python scripts/verify_repo.py
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

import ast, csv, glob, hashlib, io, json, os, re, shutil, subprocess, sys, tempfile

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
REPO = os.path.join(ROOT, "repo")
_bsv = io.open(os.path.join(BASE, "scripts", "build_submission_v1.py"), encoding="utf-8").read()
PKG = os.path.join(ROOT, "submission package " +
                   re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M).group(1))
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))

OK, BAD = [], []
def chk(label, cond, detail=""):
    (OK if cond else BAD).append(label)
    print(("  ok   " if cond else "  ✗ FAIL "), label, "" if cond else f" [{detail}]")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

ALL = [p for p in glob.glob(os.path.join(REPO, "**", "*"), recursive=True)
       if os.path.isfile(p) and ".git" + os.sep not in p]
rel = lambda p: os.path.relpath(p, REPO).replace("\\", "/")
print(f"verify_repo — {len(ALL)} files under {REPO}\n")

# ---- 1. nothing private, nothing third-party ----------------------------------------------------
leaks = [rel(p) for p in ALL if not p.lower().endswith((".xlsx", ".png", ".pdf", ".gz"))
         and re.search(r"[A-Za-z]:[\\/]Users[\\/]", io.open(p, encoding="utf-8", errors="replace").read())]
chk("no absolute local path anywhere in the deposit", not leaks, leaks[:5])
sums = [rel(p) for p in ALL if re.search(r"\.(tsv|vcf|bgen|bed|bim|fam)(\.gz)?$", p, re.I)]
chk("no raw summary-statistic file reached the deposit", not sums, sums[:5])
chk("no raw-data directory was copied", not os.path.isdir(os.path.join(REPO, "code", "data"))
    and not os.path.isdir(os.path.join(REPO, "data")))
CTRL = re.compile(r"tpmi|taiwan|koges|biobank[ _]japan", re.I)
PER_SNP = re.compile(r"^(snp|rsid|variant|marker)$", re.I)
risky = []
for p in [x for x in ALL if x.lower().endswith(".csv")]:
    with io.open(p, encoding="utf-8", errors="replace", newline="") as fh:
        hdr = next(csv.reader(fh), [])
        body = fh.read(40000)
    if any(PER_SNP.match(h.strip()) for h in hdr) and CTRL.search(body):
        risky.append(rel(p))
chk("no per-variant file from a controlled-access cohort", not risky, risky[:5])

# ---- 2. it is complete, and identical to what the journal gets ----------------------------------
for d, n_min in (("code/scripts", 100), ("code/figures", 14), ("code/results", 90),
                 ("source_data", 45)):
    n = len(glob.glob(os.path.join(REPO, d, "*")))
    chk(f"{d} holds {n} files (>= {n_min})", n >= n_min, n)
diffs = []
for p in glob.glob(os.path.join(PKG, "05_source_data", "*")):
    if os.path.basename(p) == "SOURCE_DATA_MANIFEST.md":
        continue
    q = os.path.join(REPO, "source_data", os.path.basename(p))
    if not os.path.exists(q) or md5(p) != md5(q):
        diffs.append(os.path.basename(p))
chk("every shipped source-data file is md5-identical in the deposit", not diffs, diffs[:5])
chk("the Supplementary Tables workbook is md5-identical to the shipped one",
    md5(os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"))
    == md5(os.path.join(REPO, "supplementary_tables", "Supplementary_Tables.xlsx")))

# ---- 3. it still parses ---------------------------------------------------------------------------
bad = []
for p in glob.glob(os.path.join(REPO, "code", "**", "*.py"), recursive=True):
    try:
        ast.parse(io.open(p, encoding="utf-8").read())
    except SyntaxError as e:
        bad.append(f"{os.path.basename(p)}:{e.lineno}")
chk("every Python file in the deposit parses", not bad, bad[:5])
rs = glob.glob(os.path.join(REPO, "code", "**", "*.R"), recursive=True)
if shutil.which("Rscript"):
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as fh:
        fh.write('for (f in commandArgs(TRUE)) tryCatch(invisible(parse(f)), '
                 'error=function(e) cat("FAIL", basename(f), "\\n"))\n')
        t = fh.name
    out = subprocess.run(["Rscript", "--vanilla", t] + rs, capture_output=True, text=True,
                         encoding="utf-8", errors="replace")
    os.remove(t)
    bad_r = [l for l in (out.stdout or "").split("\n") if l.startswith("FAIL")]
    chk(f"every one of the {len(rs)} R files in the deposit parses", not bad_r, bad_r[:5])
else:
    chk("R parse check ran", False, "Rscript not on PATH — this check did NOT run")

# ---- 4. the metadata says what the manuscript says -----------------------------------------------
_front = io.open(os.path.join(BASE, "figures", "FRONT_MATTER.md"), encoding="utf-8").read()
TITLE = " ".join(re.search(r"##\s*Title\s*\n\*\*(.+?)\*\*", _front, flags=re.S).group(1).split())
zen = json.load(open(os.path.join(REPO, ".zenodo.json"), encoding="utf-8"))
readme = io.open(os.path.join(REPO, "README.md"), encoding="utf-8").read()
chk("the Zenodo title is the manuscript's canonical title", TITLE in zen["title"])
chk("the README title is the manuscript's canonical title", readme.startswith("# " + TITLE))
chk("the Zenodo record names the corresponding author with an ORCID",
    any(c.get("orcid") == "0000-0002-2218-7115" for c in zen["creators"]))
chk("the licence pair is declared (MIT code + CC BY 4.0 data)",
    zen["license"] == "cc-by-4.0"
    and "MIT" in io.open(os.path.join(REPO, "LICENSE"), encoding="utf-8").read()
    and "CC BY 4.0" in io.open(os.path.join(REPO, "LICENSE"), encoding="utf-8").read())
chk("the README states that raw summary statistics are NOT redistributed",
    "not redistributed" in readme and "controlled-access" in readme)
# Co-authors are now supplied in .zenodo.json (2026-08-14), so the only pre-release item the deposit
# still discloses is the Zenodo DOI, which cannot exist until the release is archived; it lives in the
# README. The old check also required AUTHOR-SUPPLIED inside .zenodo.json, which was satisfied only by
# the co-author placeholder that has legitimately been filled.
# 2026-09-17: v1.2.0 is archived, so "the DOI cannot exist yet" is no longer true and the old check
# (README must still say AUTHOR-SUPPLIED) would now ENFORCE a stale placeholder. Repointed to the
# truth: the README cites the CONCEPT DOI from manifest/DEPOSIT.json, and never the version DOI.
_DEP = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                   "manifest", "DEPOSIT.json"), encoding="utf-8"))
chk("the deposit README cites the Zenodo CONCEPT DOI, not the version DOI",
    _DEP["concept_doi"] in readme and _DEP["version_doi"] not in readme
    and "AUTHOR-SUPPLIED" not in readme)
for pat, want in ((r"\*\*(\d+)\*\* non-self directed", F["network_edges_total"]),
                  (r"\*\*(\d+)\*\*\s*\n?edges pass", F["edges_passing_bonferroni"])):
    m = re.search(pat, readme)
    chk(f"the README figure {want} is derived, not typed", m and int(m.group(1)) == want,
        f"README={m and m.group(1)}")

# ---- 5. THE LOAD-BEARING CHECK: it actually runs -------------------------------------------------
# A grep proves the strings are gone. Only a render proves the rewrite works. suppfig2.R is chosen
# because it reads results/ and writes a figure, needs no TwoSampleMR and no raw data.
if shutil.which("Rscript"):
    env = dict(os.environ, CKM_P4_BASE=os.path.join(REPO, "code").replace("\\", "/"))
    out = subprocess.run(["Rscript", "--vanilla", os.path.join(REPO, "code", "figures", "suppfig2.R")],
                         capture_output=True, text=True, encoding="utf-8", errors="replace", env=env,
                         cwd=os.path.join(REPO, "code", "figures"))
    made = os.path.join(REPO, "code", "figures", "SupplFig2.png")
    ran = out.returncode == 0 and os.path.exists(made)
    chk("a figure RENDERS from the deposit alone, with CKM_P4_BASE set", ran,
        (out.stderr or "")[-300:])
    if ran:
        shipped = os.path.join(PKG, "04_supplementary", "SupplFig2.png")
        chk("the deposit's render is the same size as the shipped figure (within 2%)",
            abs(os.path.getsize(made) - os.path.getsize(shipped)) / os.path.getsize(shipped) < 0.02,
            f"{os.path.getsize(made)} vs {os.path.getsize(shipped)}")
        for f in glob.glob(os.path.join(REPO, "code", "figures", "SupplFig2.*")):
            os.remove(f)                       # rendered outputs are not tracked
else:
    chk("render smoke test ran", False, "Rscript not on PATH — this check did NOT run")

print("\n" + "=" * 92)
print(f"RESULT: {len(OK)} verified, {len(BAD)} FAILED")
for b in BAD:
    print("  x", b)
print("=" * 92)
sys.exit(1 if BAD else 0)
