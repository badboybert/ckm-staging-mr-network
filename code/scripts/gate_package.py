# -*- coding: utf-8 -*-
"""Pre-submission gate for the CKM Paper 4 frozen package. Hard-fails.

Three FORWARD-LOOKING gates. They assert what should be TRUE rather than sweeping for past
mistakes, because a regression sweep cannot catch a number that is simply wrong:

  Gate A  no internal provenance or unresolved placeholders in any shipped file
  Gate B  cross-file identity — every copy is md5-identical to its source, and facts that appear
          in more than one file are mechanically diffed
  Gate C  numbers against source — every count quoted anywhere in the package is recomputed from
          the section files, legends, workbook and reference builder

SCOPE IS THE GATE. The package is enumerated first; every file must be claimed by at least one
gate, and unclaimed files are a FAILURE. A gate that silently sweeps a subset manufactures false
confidence, which is worse than no gate.
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

import os, re, sys, csv, json, glob, zlib, hashlib, subprocess, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_vancouver as VAN

BASE = P4_BASE
ROOT = P4_ROOT
FIG = os.path.join(BASE, "figures")
MAN = os.path.join(BASE, "manuscript")
RES = os.path.join(BASE, "results")
# Must match PKG_VERSION in build_submission_v1.py. Pointing this at a version you did not just cut
# reports that OTHER package's staleness as if it were the new one -- which is exactly what happened
# when v2 was first cut while this still read "v1".
_bsv = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_submission_v1.py"),
            encoding="utf-8").read()
_mv = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _mv, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _mv.group(1)   # DERIVED, never typed: a gate pointed at a stale package is worse than none
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")

FAILS, NOTES = [], []
def fail(m): FAILS.append(m)
def note(m): NOTES.append(m)

def rd(p):
    return open(p, encoding="utf-8").read() if os.path.exists(p) else ""

def md5(p, chunk=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

# ---------------- text extraction, so .docx and .xlsx are not blind spots -------------
def docx_text(p):
    from docx import Document
    d = Document(p)
    parts = [para.text for para in d.paragraphs]
    for t in d.tables:
        for row in t.rows:
            for c in row.cells:
                parts.append(c.text)
    return "\n".join(parts)

def xlsx_text(p):
    import openpyxl
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    parts = list(wb.sheetnames)
    for ws in wb.worksheets:
        for row in ws.iter_rows(values_only=True):
            for v in row:
                if isinstance(v, str):
                    parts.append(v)
    wb.close()
    return "\n".join(parts)

def text_of(p):
    e = os.path.splitext(p)[1].lower()
    if e in (".md", ".txt", ".csv", ".tsv", ".json"):
        return open(p, encoding="utf-8", errors="replace").read()
    if e == ".docx":
        return docx_text(p)
    if e == ".xlsx":
        return xlsx_text(p)
    return None                     # binary (png/pdf) — no text surface

# ================= 0. ENUMERATE THE ARTIFACT FIRST =================
ALL = []
for root, _, files in os.walk(PKG):
    for f in files:
        ALL.append(os.path.join(root, f))
ALL.sort()
if not ALL:
    print("GATE ABORT: package not found at", PKG); sys.exit(1)
rel = lambda p: os.path.relpath(p, PKG).replace("\\", "/")
claimed = set()

# facts, recomputed now (never read from a stored copy)
subprocess.run([sys.executable, os.path.join(BASE, "scripts", "derive_facts.py")],
               capture_output=True, text=True, check=True)
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))
VAN_MASTER, VAN_ORDER, VAN_AMAP = VAN.numbering()

# ================= GATE A — no internal provenance =================
# Split by audience: everything is checked, but the journal-facing folders are held to the
# stricter vocabulary. _author_todo and README are ours and may legitimately discuss next steps.
# Case matters. The working-tree markers are ALL-CAPS tokens ("REMAINING =", "TODO:"); the same
# letters in ordinary prose ("the remaining ambiguity", "Remaining author items") are legitimate
# English and must not fire. Each pattern therefore carries its own flags.
IC = re.I | re.M
CS = re.M
STRICT_ALL = [
    (r"wf_[0-9a-f]{6,}", "workflow ID", IC),
    (r"\badversarial\w*", "verification jargon", IC),
    (r"\bverifier\b", "verification jargon", IC),
    (r"(?:PASS|SOUND|CONFIRMED)_WITH_(?:FIXES|CAVEATS)", "verifier verdict", CS),
    (r"[A-Za-z]:[\\/]Users[\\/]", "local absolute path", IC),
    (r"\bMILESTONE \d", "working-tree milestone marker", CS),
]
STRICT_UPLOAD = [
    (r"\bTODO\b", "todo marker", CS),
    (r"\bREMAINING\b", "status marker", CS),
    (r"\[confirm\]|\[REF\]|\[citation to confirm\]", "unresolved citation", IC),
    (r"^_?DRAFT\b", "draft marker", CS),
    (r"\bBioRender\b", "production note", IC),
]
UPLOAD_DIRS = ("01_manuscript", "02_cover_declarations", "03_main_figures",
               "04_supplementary", "05_source_data")
INTERNAL = ("_author_todo", "README_SUBMISSION.md", "MANIFEST_checksums.txt")

placeholders = []
binary_files = 0
for p in ALL:
    r = rel(p)
    t = text_of(p)
    claimed.add(r)
    if t is None:
        binary_files += 1
        continue
    is_upload = r.split("/")[0] in UPLOAD_DIRS
    pats = STRICT_ALL + (STRICT_UPLOAD if is_upload else [])
    for pat, why, flags in pats:
        for m in re.finditer(pat, t, flags=flags):
            frag = " ".join(t[max(0, m.start() - 45):m.start() + 55].split())
            fail(f"GATE A [{why}] {r}: ...{frag}...")
    for m in re.finditer(r"\[AUTHOR-SUPPLIED[^\]]*\]|\[Author list to be supplied\]|"
                         r"\[GitHub URL[^\]]*\]|\[AUTHOR-SUPPLIED\]", t):
        placeholders.append((r, " ".join(m.group(0).split())[:70]))

# ================= GATE B — cross-file identity =================
# B1: every copied file must be md5-identical to the file it came from.
COPY_MAP = {}
for i in range(1, F["main_figures"] + 1):
    for e in ("png", "pdf"):
        COPY_MAP[f"03_main_figures/Figure{i}.{e}"] = os.path.join(FIG, f"Figure{i}.{e}")
for i in range(1, F["supp_figures"] + 1):
    for e in ("png", "pdf"):
        COPY_MAP[f"04_supplementary/SupplFig{i}.{e}"] = os.path.join(FIG, f"SupplFig{i}.{e}")
COPY_MAP["04_supplementary/Supplementary_Tables.xlsx"] = os.path.join(MAN, "Supplementary_Tables.xlsx")
COPY_MAP["01_manuscript/CKM_Paper4_manuscript.docx"] = os.path.join(MAN, "CKM_Paper4_manuscript.docx")

for r, src in COPY_MAP.items():
    dst = os.path.join(PKG, r.replace("/", os.sep))
    if not os.path.exists(dst):
        fail(f"GATE B [missing copy] {r}"); continue
    if not os.path.exists(src):
        fail(f"GATE B [missing source] {src}"); continue
    if md5(dst) != md5(src):
        fail(f"GATE B [copy differs from source] {r}")

# B1b: source-data copies must match their origin
sd_dir = os.path.join(PKG, "05_source_data")
for f in os.listdir(sd_dir):
    if f == "SOURCE_DATA_MANIFEST.md":
        continue
    for d_ in (RES, os.path.join(FIG, "suppfig_data"), os.path.join(FIG, "data")):
        cand = os.path.join(d_, f)
        if os.path.exists(cand):
            if md5(os.path.join(sd_dir, f)) != md5(cand):
                fail(f"GATE B [source-data copy differs] 05_source_data/{f}")
            break
    else:
        fail(f"GATE B [source-data file has no origin in the analysis tree] {f}")

# B2: the title must be byte-identical everywhere it appears.
front = rd(os.path.join(FIG, "FRONT_MATTER.md"))
mt = re.search(r"##\s*Title\s*\n\*\*(.+?)\*\*", front, flags=re.S)
TITLE = " ".join(mt.group(1).split()) if mt else None
if not TITLE:
    fail("GATE B [no title in FRONT_MATTER.md]")
else:
    tp = rd(os.path.join(PKG, "02_cover_declarations", "TITLE_PAGE.md"))
    if TITLE not in " ".join(tp.split()):
        fail("GATE B [title differs] TITLE_PAGE.md does not carry the front-matter title")
    man_txt = " ".join(docx_text(os.path.join(PKG, "01_manuscript",
                                              "CKM_Paper4_manuscript.docx")).split())
    if TITLE not in man_txt:
        fail("GATE B [title differs] manuscript .docx does not carry the front-matter title")
    rm = " ".join(rd(os.path.join(PKG, "README_SUBMISSION.md")).split())
    if TITLE not in rm:
        fail("GATE B [title differs] README_SUBMISSION.md does not carry the front-matter title")
    # B2b: the title is ALSO hand-typed in the abstract and the cover letter (both shipped, both
    # read by the editor). Those two surfaces were outside every gate until 2026-07-25, which is
    # exactly the hand-duplicated-fact hole that let a retired title survive a rename before.
    for _sub, _f in (("01_manuscript", "ABSTRACT.md"),
                     ("02_cover_declarations", "COVER_LETTER.md"),
                     ("02_cover_declarations", "COVER_LETTER.docx")):
        _p = os.path.join(PKG, _sub, _f)
        _t = docx_text(_p) if _f.endswith(".docx") else rd(_p)
        if TITLE not in " ".join(_t.split()):
            fail(f"GATE B [title differs] {_sub}/{_f} does not carry the front-matter title")

# B3: shipped section .md must equal the cleaned source (no drift, no leftover provenance)
def clean_md(t):
    t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    t = "\n".join(ln for ln in t.split("\n") if not re.match(r"^_.+_$", ln.rstrip()))
    return re.sub(r"\n{3,}", "\n\n", t).strip() + "\n"

# The shipped .md is the cleaned source AFTER the build-time author-date -> Vancouver transform.
# The gate applies the same transform from the same module rather than trusting the builder, so a
# change in either one shows up here as drift instead of shipping silently.
for f in ["INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md", "CONCLUSIONS.md",
          "ABBREVIATIONS.md", "ABSTRACT.md", "FIGURE_LEGENDS.md"]:
    shipped = rd(os.path.join(PKG, "01_manuscript", f))
    src = rd(os.path.join(FIG, f))
    # ABSTRACT.md is the one section the builder TRUNCATES: only the FULL abstract ships, because
    # the source also holds a stricter-cap variant whose wording predates the round-3 recalibrations.
    # This gate must apply the same truncation, or it either fails on every build or - worse - is
    # weakened to a substring test that would stop noticing real drift.
    if f == "ABSTRACT.md":
        _cut = src.split("## CAP variant")
        if len(_cut) != 2:
            fail("GATE B [abstract] source ABSTRACT.md must contain exactly one '## CAP variant'")
        src = _cut[0].rstrip().rstrip("-").rstrip() + "\n"
    if f in VAN.REWRITE:
        src = VAN.to_vancouver(src, VAN_AMAP, VAN_MASTER)
    expect = clean_md(src)
    if shipped != expect:
        fail(f"GATE B [section drift] 01_manuscript/{f} is not the cleaned, Vancouver-converted source")

# B3b: no shipped text surface may still carry an author-date citation, and the in-text numbers must
# be exactly 1..N with no gaps. A gap means a reference is listed but never cited, which is the
# defect Vancouver conversions produce most often.
_cited = set()
for _f in sorted(glob.glob(os.path.join(PKG, "01_manuscript", "*.md"))):
    # REFERENCES_NUMBERED.md is the BIBLIOGRAPHY plus its author-date -> number map. Author names
    # and years are its content, not unconverted citations, so scanning it for author-date tokens
    # reports the file doing its job as a defect.
    if os.path.basename(_f) == "REFERENCES_NUMBERED.md":
        continue
    _t = rd(_f)
    _left = VAN.unconverted(_t, VAN_MASTER, VAN_AMAP)
    if _left:
        fail(f"GATE B [author-date survives Vancouver conversion] "
             f"01_manuscript/{os.path.basename(_f)}: {_left[:3]}")
    for _m in re.finditer(r"\[(\d+(?:,\d+)*)\]", _t):
        _cited |= {int(x) for x in _m.group(1).split(",")}
_expected = set(range(1, len(VAN_ORDER) + 1))
if _cited and _cited != _expected:
    fail(f"GATE B [Vancouver numbering] cited numbers are not 1..{len(VAN_ORDER)} "
         f"(missing {sorted(_expected - _cited)}, extra {sorted(_cited - _expected)})")

# B4: checksum manifest must be self-consistent and complete
mf = rd(os.path.join(PKG, "MANIFEST_checksums.txt"))
listed = {}
for line in mf.split("\n"):
    m = re.match(r"^([0-9a-f]{32})\s+(\d+)\s+(.+)$", line.strip())
    if m:
        listed[m.group(3)] = (m.group(1), int(m.group(2)))
for p in ALL:
    r = rel(p)
    if r == "MANIFEST_checksums.txt":
        continue
    if r not in listed:
        fail(f"GATE B [not in manifest] {r}")
    else:
        h, sz = listed[r]
        if md5(p) != h:
            fail(f"GATE B [manifest md5 mismatch] {r}")
        if os.path.getsize(p) != sz:
            fail(f"GATE B [manifest size mismatch] {r}")
for r in listed:
    if not os.path.exists(os.path.join(PKG, r.replace("/", os.sep))):
        fail(f"GATE B [manifest lists a missing file] {r}")

# ================= GATE C — numbers against source =================
w = F["word_counts"]
EXPECT = {
    "abstract full": F["abstract_full_words"],
    "abstract cap": F["abstract_cap_words"],
    "introduction": w["INTRODUCTION"],
    "methods": w["METHODS"],
    "results": w["RESULTS"],
    "discussion": w["DISCUSSION"],
    "conclusions": w["CONCLUSIONS"],
    "main text": F["main_text_words"],
    "panels": F["total_panels"],
    "main figures": F["main_figures"],
    "supp figures": F["supp_figures"],
    "supp tables": F["supp_tables"],
    "worksheets": F["supp_table_sheets"],
    "references": F["references_cited"],
}

def num(s):
    return int(str(s).replace(",", ""))

# Round 5 (M7) moved the derived metrics OFF the journal-facing title page — a reviewer should not
# be shown internal word counts, an unshipped "short variant" or a pointer to an editable source
# file. The numbers must still be derived and still be checked, so this block now reads the README,
# which is internal by Gate A's audience split and already carried the same table. The gate was not
# stale in INTENT, only in SCOPE: it was asserting the right facts against the wrong surface. The
# title page is separately asserted to be free of them, below.
tp_txt = rd(os.path.join(PKG, "README_SUBMISSION.md"))
checks = [
    (r"\| Abstract \| ([\d,]+) words full / ([\d,]+) words cap variant \|",
     ("abstract full", "abstract cap")),
    (r"\| Introduction \| ([\d,]+) words \|", ("introduction",)),
    (r"\| Methods \| ([\d,]+) words \|", ("methods",)),
    (r"\| Results \| ([\d,]+) words \|", ("results",)),
    (r"\| Discussion \| ([\d,]+) words \|", ("discussion",)),
    (r"\| Conclusions \| ([\d,]+) words \|", ("conclusions",)),
    (r"\| Main text total \| ([\d,]+) words \|", ("main text",)),
    (r"\| Main figures \| (\d+) figures, (\d+) panels \|", ("main figures", "panels")),
    (r"\| Supplementary tables \| S1–S(\d+) in (\d+) worksheets \|", ("supp tables", "worksheets")),
    (r"\| Supplementary figures \| (\d+) \|", ("supp figures",)),
    (r"\| References \| (\d+) cited", ("references",)),
]
for pat, keys in checks:
    m = re.search(pat, tp_txt)
    if not m:
        fail(f"GATE C [title page missing a metric] /{pat}/")
        continue
    for i, k in enumerate(keys):
        got = num(m.group(i + 1))
        if got != EXPECT[k]:
            fail(f"GATE C [title page] {k}: says {got}, source says {EXPECT[k]}")

readme = rd(os.path.join(PKG, "README_SUBMISSION.md"))
for pat, k in [
    (r"\|\s*Introduction\s*\|\s*([\d,]+) words", "introduction"),
    (r"\|\s*Methods\s*\|\s*([\d,]+) words", "methods"),
    (r"\|\s*Results\s*\|\s*([\d,]+) words", "results"),
    (r"\|\s*Discussion\s*\|\s*([\d,]+) words", "discussion"),
    (r"\|\s*Conclusions\s*\|\s*([\d,]+) words", "conclusions"),
    (r"\|\s*Main text total\s*\|\s*([\d,]+) words", "main text"),
    (r"\|\s*Main figures\s*\|\s*(\d+) figures", "main figures"),
    (r"\|\s*Main figures\s*\|\s*\d+ figures,\s*(\d+) panels", "panels"),
    (r"\|\s*Supplementary figures\s*\|\s*(\d+)", "supp figures"),
    (r"\|\s*References\s*\|\s*(\d+) cited", "references"),
]:
    m = re.search(pat, readme)
    if not m:
        fail(f"GATE C [README missing a metric] /{pat}/")
    elif num(m.group(1)) != EXPECT[k]:
        fail(f"GATE C [README] {k}: says {num(m.group(1))}, source says {EXPECT[k]}")

# C2: references must be internally consistent (0 uncited / 0 unmatched), and the docx must
# actually render that many.
if F["references_uncited"] != 0:
    fail(f"GATE C [references] {F['references_uncited']} uncited")
if F["references_unmatched"] != 0:
    fail(f"GATE C [references] {F['references_unmatched']} unmatched")

# C3: display items on disk must match the counts claimed
n_fig_files = len(glob.glob(os.path.join(PKG, "03_main_figures", "Figure*.png")))
n_pdf_files = len(glob.glob(os.path.join(PKG, "03_main_figures", "Figure*.pdf")))
if n_fig_files != F["main_figures"] or n_pdf_files != F["main_figures"]:
    fail(f"GATE C [main figures] {n_fig_files} png / {n_pdf_files} pdf vs {F['main_figures']} claimed")
n_sf_png = len(glob.glob(os.path.join(PKG, "04_supplementary", "SupplFig*.png")))
n_sf_pdf = len(glob.glob(os.path.join(PKG, "04_supplementary", "SupplFig*.pdf")))
if n_sf_png != F["supp_figures"] or n_sf_pdf != F["supp_figures"]:
    fail(f"GATE C [supp figures] {n_sf_png} png / {n_sf_pdf} pdf vs {F['supp_figures']} claimed")

# C3b: the network size, Bonferroni threshold and significant-edge count are quoted in the Methods,
# the Results and the Supplementary Figure S1 legend. Recompute all three from the edge table.
n_edges = F["network_edges_total"]
n_sig = F["edges_passing_bonferroni"]
methods_txt = rd(os.path.join(PKG, "01_manuscript", "METHODS.md"))
m = re.search(r"0\.05\s*/\s*(\d+)\s*=\s*([\d.]+)\s*×\s*10⁻⁴", methods_txt)
if not m:
    fail("GATE C [Methods does not state the Bonferroni threshold in the expected form]")
else:
    if int(m.group(1)) != n_edges:
        fail(f"GATE C [Bonferroni denominator] Methods says 0.05/{m.group(1)}, "
             f"edge table has {n_edges} edges")
    stated = float(m.group(2)) * 1e-4
    if abs(stated - F["bonferroni_threshold"]) > 0.06e-4:
        fail(f"GATE C [Bonferroni threshold] Methods says {stated:.2e}, "
             f"recomputed {F['bonferroni_threshold']:.2e}")
sf_txt = rd(os.path.join(FIG, "SUPPL_FIGURES.md"))
m = re.search(r"(\d+) of (\d+) edges are Bonferroni-significant", sf_txt)
if m:
    if int(m.group(1)) != n_sig or int(m.group(2)) != n_edges:
        fail(f"GATE C [significant-edge count] legend says {m.group(1)} of {m.group(2)}, "
             f"recomputed {n_sig} of {n_edges}")

# C4: every source path cited in a legend must resolve to a shipped file
have = set(os.listdir(sd_dir))
legend_txt = rd(os.path.join(FIG, "FIGURE_LEGENDS.md")) + rd(os.path.join(FIG, "SUPPL_FIGURES.md"))
for m in re.finditer(r"`([^`]+\.(?:csv|tsv))`", legend_txt):
    path = m.group(1)
    b = os.path.basename(path)
    if "<" in path or "*" in path:          # a documented pipeline-input template, not one file
        continue
    if b not in have:
        fail(f"GATE C [legend cites unshipped source data] {path}")

# C5: docx author metadata (global rule 8)
from docx import Document
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    a = Document(p).core_properties.author
    if a != "Bertrand Chin-Ming Tan":
        fail(f"GATE C [docx creator] {rel(p)}: author is {a!r}")

# C6: every shipped .docx must actually CONTAIN something. A 36 KB .docx with zero paragraphs
# looks fine on disk and in a checksum manifest; DECLARATIONS.docx shipped empty because a
# DOTALL regex in the builder swallowed its source text. Assert content, not just existence.
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    d = Document(p)
    paras = [x.text for x in d.paragraphs if x.text.strip()]
    ncells = sum(1 for t in d.tables for row in t.rows for c in row.cells if c.text.strip())
    if len(paras) + ncells < 5:
        fail(f"GATE C [empty deliverable] {rel(p)}: {len(paras)} paragraphs, {ncells} table cells")

# C7: no raw markdown may survive into a rendered .docx. Inline spans that wrap across lines
# (the two-line bold title) render as literal asterisks unless the builder joins them first.
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    t = docx_text(p)
    # The underscore pattern is boundary-anchored so filenames (Supplementary_Tables.xlsx,
    # forward_local_edges.csv) are not mistaken for unrendered emphasis.
    for pat, why in [(r"\*\*", "literal bold markers"),
                     (r"(?<![\w`])_[^_`\n]+_(?![\w])", "literal italic markers"),
                     (r"<!--", "unstripped HTML comment"),
                     (r"\]\(http", "unrendered markdown link")]:
        if re.search(pat, t):
            m = re.search(pat, t)
            frag = " ".join(t[max(0, m.start() - 40):m.start() + 60].split())
            fail(f"GATE C [{why}] {rel(p)}: ...{frag}...")

# C-extra: a legend must not quote a different value for a statistic than the panel prints.
# Figure 1 reports the label-permutation P in panels a, b, d, f and g and in its legend; the battery
# script and the main staging engine are two Monte-Carlo runs of the SAME test, and the legend went
# on quoting the battery's last digit after the panel was aligned to the canonical value.
_legends = rd(os.path.join(PKG, "01_manuscript", "FIGURE_LEGENDS.md"))
_canon_p = float(open(os.path.join(FIG, "data", "fig1_staging_stat.csv"), encoding="utf-8")
                 .read().split("perm_p,")[1].splitlines()[0])
def _as_float(s):
    s = s.replace("×", "x").replace("−", "-")
    for a, b in zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-"):
        s = s.replace(a, b)
    m = re.match(r"([\d.]+)\s*x\s*10(-?\d+)", s.replace(" ", " ").strip())
    return float(m.group(1)) * 10 ** int(m.group(2)) if m else float(s)
def _sig_digits(s):
    """Significant digits in the quoted mantissa, so the canonical value can be compared at the
    legend's OWN precision. The store now holds full precision (3/1980 = 0.001515...), so an exact
    equality test would reject the correct 2-s.f. rendering '1.5 × 10⁻³' while a loose tolerance
    would accept a genuinely wrong one. Rounding the canon to the quoted precision does neither."""
    d = re.sub(r"[^\d]", "", re.split(r"\s*[x×]\s*10", s.replace("−", "-"))[0].lstrip("0.").lstrip("0"))
    return max(1, len(d.rstrip("0")) or 1)
def _round_sig(x, n):
    return 0.0 if x == 0 else round(x, -int(math.floor(math.log10(abs(x)))) + (n - 1))
for _m in re.finditer(r"label permutation \(\*P\* = ([^\s—)]+(?:\s*×\s*10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)?)", _legends):
    try:
        _v = _as_float(_m.group(1))
    except (ValueError, AttributeError):
        continue
    if abs(_v - _round_sig(_canon_p, _sig_digits(_m.group(1)))) > 1e-12:
        fail(f"GATE C [legend/panel disagree] FIGURE_LEGENDS quotes label-permutation P = {_m.group(1)} "
             f"but the panels print {_canon_p:.4g}")

# ================= GATE D — defect classes found by the 2026-07-18 audit =================
# Each of these shipped in v1's first cut. They are asserted here so the class cannot return.

# D0: SUPERSEDED NUMBERS ANYWHERE IN THE PACKAGE.
# Added 2026-07-23. Gate C checks the main sections; it does NOT read the cover letter, the STROBE
# checklist or the supplementary-information document. The v2 cut shipped the retired staging P in
# COVER_LETTER.md and the retired 111-edge Bonferroni in STROBE_MR_CHECKLIST.md, and every existing
# gate passed. A gate's SCOPE is the gate: this one walks EVERY shipped text file.
SUPERSEDED_PKG = {
    "111-edge Bonferroni":      r"0\.05/111|4\.5 × 10⁻⁴",
    "111 as current network":   r"111[- ]edge|111 directed|the 111\b",
    "superseded staging P":     r"0\.926[^.]{0,80}1\.8 × 10⁻³|0\.926[^.]{0,80}1\.2 × 10⁻³|permutation[^.]{0,40}1\.[28] × 10⁻³",
    "50-edge graph (now 54)":   r"\b50-edge\b",
    "retired Kendall tau":      r"Kendall|τ = 0\.49",
    "abandoned exclusive T2D":  r"only via CAD|reaches heart failure only",
    "pre-registered":           r"[Pp]re-registered",
}
for p in ALL:
    t = text_of(p)
    if t is None:
        continue
    for label, pat in SUPERSEDED_PKG.items():
        m = re.search(pat, t)
        if m:
            frag = " ".join(t[max(0, m.start()-50):m.start()+60].split())
            fail(f"GATE D [superseded: {label}] {rel(p)}: ...{frag}...")

# D1: no pointer to a supplementary figure beyond the number actually shipped (S6 was retired).
n_sf = F["supp_figures"]
for p in ALL:
    t = text_of(p)
    if t is None:
        continue
    for m in re.finditer(r"Supplementary Figures? S(\d+)(?:\s*[–-]\s*S(\d+))?", t):
        hi = int(m.group(2) or m.group(1))
        if hi > n_sf:
            fail(f"GATE D [pointer to unshipped Supplementary Figure S{hi}] {rel(p)}: "
                 f"{' '.join(m.group(0).split())} (package ships S1–S{n_sf})")

# D2: a document may not name a companion .md file that is not in the package.
for p in ALL:
    r = rel(p)
    if r.split("/")[0] not in UPLOAD_DIRS:
        continue
    t = text_of(p)
    if t is None:
        continue
    for m in re.finditer(r"\b([A-Z][A-Z0-9_]{3,})\.md\b", t):
        name = m.group(0)
        if not any(os.path.basename(x) == name for x in ALL):
            fail(f"GATE D [names a file absent from the package] {r}: {name}")

# D3: bare causal assertions breach calibration invariant 1, including inside rendered figures.
for p in ALL:
    t = text_of(p)
    if t is None:
        continue
    for m in re.finditer(r"\b(?:is|are|was|were)\s+causal\b", t, flags=re.I):
        frag = " ".join(t[max(0, m.start() - 60):m.start() + 60].split())
        fail(f"GATE D [bare causal assertion] {rel(p)}: ...{frag}...")

# D4: no document may ship python-docx's template creation stamp.
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    c = Document(p).core_properties.created
    if c is not None and c.year == 2013:
        fail(f"GATE D [python-docx template timestamp] {rel(p)}: created={c.isoformat()}")

# D5: the canonical title must also hold inside the workbook, which Gate B's .docx/.md sweep missed.
if TITLE:
    import openpyxl as _ox
    _wb = _ox.load_workbook(os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"),
                            read_only=True)
    _a1 = str(_wb["Contents"]["A1"].value or "")
    if TITLE not in " ".join(_a1.split()):
        fail(f"GATE D [workbook title differs] Contents!A1 = {_a1[:90]!r}")
    _desc = str(_wb.properties.description or "")
    _ids = sorted({int(m.group(1)) for s in _wb.sheetnames
                   for m in [re.match(r"S(\d+)", s)] if m})
    if _ids and f"S{_ids[0]}–S{_ids[-1]}" not in _desc:
        fail(f"GATE D [workbook description range] says {_desc[:70]!r}, "
             f"sheets span S{_ids[0]}–S{_ids[-1]}")
    _wb.close()

# D6: every supplementary table and figure should be cited somewhere in the main text. Reported as
# a NOTE, not a failure: where each callout belongs is an editorial decision, not a mechanical one.
body = "".join(rd(os.path.join(PKG, "01_manuscript", f)) for f in
               ["ABSTRACT.md", "INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md",
                "FIGURE_LEGENDS.md"])
tab_ids = sorted({m.group(1) for s in
                  __import__("openpyxl").load_workbook(
                      os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"),
                      read_only=True).sheetnames
                  for m in [re.match(r"(S\d+b?)", s)] if m})
# Callouts appear in plural and list form ("Supplementary Tables S4 and S4b"), so collect every
# identifier that follows a "Supplementary Table(s)/Figure(s)" lead-in rather than requiring the
# singular. Matching only "Table S4" reported four cited tables as uncited.
def cited_ids(kind):
    seen = set()
    for m in re.finditer(rf"Supplementary\s+{kind}s?\.?\s+((?:S\d+[a-z]?(?:\s*(?:,|and|–|-)\s*)?)+)",
                         body, flags=re.I):
        seen |= set(re.findall(r"S\d+[a-z]?", m.group(1)))
    return seen

cited_t, cited_f = cited_ids("Table"), cited_ids("Fig(?:ure)?")
uncited_t = [t_ for t_ in tab_ids if t_ not in cited_t]
uncited_f = [f"S{i}" for i in range(1, n_sf + 1) if f"S{i}" not in cited_f]
if uncited_t:
    note(f"supplementary tables with no main-text callout: {', '.join(uncited_t)}")
if uncited_f:
    note(f"supplementary figures with no main-text callout: {', '.join(uncited_f)}")

# ================= GATE E — claims and rendered numbers =================
# Gate D0 tests a fixed list of strings found by an earlier audit: it is a REGRESSION gate and cannot
# see a claim that is simply no longer true. Gate E adds the two checks that would have caught this
# round's blockers — a claim scan with a stated scope over every shipped text file, and a scan of the
# figures' PDF TEXT LAYER, where a superseded threshold had been rendered into the axis of a
# supplementary figure while every text gate passed.

# E1: claim scan (revision history / retracted claims / superseded numbers) over the whole package.
_scan = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_claims.py")
if os.path.exists(_scan):
    # encoding must be forced: the default console codec (cp950 here) raised inside subprocess's
    # reader thread on the scan's UTF-8 output, so this gate appeared to pass without running.
    _r = subprocess.run([sys.executable, _scan, PKG], capture_output=True, text=True,
                        encoding="utf-8", errors="replace")
    if "CLAIM SCAN" not in (_r.stdout or ""):
        fail("GATE E [claim scan] produced no verdict line — treat as NOT RUN. stderr: "
             + (_r.stderr or "")[:300])
    if _r.returncode != 0:
        for _ln in (_r.stdout or "").splitlines():
            if _ln.strip().startswith("["):
                fail("GATE E [claim scan] " + _ln.strip())
        if not any(l.strip().startswith("[") for l in (_r.stdout or "").splitlines()):
            fail("GATE E [claim scan] scan_claims.py exited %d without reportable hits — "
                 "treat as NOT RUN. stderr: %s" % (_r.returncode, (_r.stderr or "")[:300]))
else:
    fail("GATE E [claim scan] scripts/scan_claims.py is missing — the claim scan cannot run")

# E2: the figures' text layer. A number drawn INTO a panel is invisible to every .md/.docx check.
_FIG_BAD = {
    "superseded Bonferroni": r"0\.05/111|4\.5e-0?4",
    "111 as the current network": r"111[- ]edge",
    "retired staging P": r"(?<![\d.])0\.0018(?![\d])|1\.8e-0?3",
    "retired portability statistic": r"Kendall|44/65",
    "pre-registered": r"[Pp]re-registered",
    # Added 2026-07-24. Figure 6 shipped the pre-enumeration staging P and this list could not see it:
    # every entry above came from an earlier audit, so the list only ever caught yesterday's mistake.
    "pre-enumeration staging P": r"(?<![\d.])0\.0012(?![\d])|1\.2e-0?3",
    "pre-enumeration CAC staging P": r"(?<![\d.])0\.0001(?![\d])",
    # A retired METHOD DESCRIPTION is as wrong as a retired number, and nothing was checking for one.
    # The staging null is an exhaustive enumeration; no figure legitimately prints "10,000" (verified
    # across all 11 shipped figures before this rule was added).
    "Monte-Carlo null description": r"10,?000\s*(shuffle|iteration|draw|permutation|replicate)",
}
try:
    import pypdf
    _pdfs = sorted(glob.glob(os.path.join(PKG, "03_main_figures", "*.pdf")) +
                   glob.glob(os.path.join(PKG, "04_supplementary", "*.pdf")))
    if not _pdfs:
        fail("GATE E [figure text layer] no figure PDFs found to scan")
    for _p in _pdfs:
        _txt = pypdf.PdfReader(_p).pages[0].extract_text() or ""
        for _label, _pat in _FIG_BAD.items():
            _m = re.search(_pat, _txt)
            if _m:
                _frag = " ".join(_txt[max(0, _m.start() - 45):_m.end() + 45].split())
                fail(f"GATE E [rendered into a figure: {_label}] {rel(_p)}: ...{_frag}...")
except ImportError:
    fail("GATE E [figure text layer] pypdf is not installed — the figure scan cannot run")

# ================= SCOPE ASSERTION =================
unclaimed = [rel(p) for p in ALL if rel(p) not in claimed]
if unclaimed:
    fail(f"SCOPE [files no gate opened] {unclaimed}")

# ================= REPORT =================
print(f"Package : {PKG}")
print(f"Scope   : {len(ALL)} files enumerated, {len(claimed)} opened by Gate A "
      f"({binary_files} binary — no text surface), {len(COPY_MAP)} md5-compared to source")
print(f"Facts   : intro {w['INTRODUCTION']:,} / methods {w['METHODS']:,} / results {w['RESULTS']:,} / "
      f"discussion {w['DISCUSSION']:,} = {F['main_text_words']:,} words; "
      f"{F['main_figures']} figs / {F['total_panels']} panels; {F['supp_figures']} supp figs; "
      f"S1–S{F['supp_tables']} in {F['supp_table_sheets']} sheets; {F['references_cited']} refs")

if placeholders:
    print(f"\nAuthor-supplied placeholders ({len(placeholders)}) — intentional, tracked in AUTHOR_TODO.md:")
    seen = set()
    for r, frag in placeholders:
        if (r, frag) in seen:
            continue
        seen.add((r, frag))
        print(f"  - {r}: {frag}")

# The target journal is now fixed (Cardiovascular Diabetology), so its abstract limit is a HARD gate,
# not an author decision. CVD Research articles: "The Abstract should not exceed 350 words", structured
# Background / Methods / Results / Conclusions. Two counts are asserted because journals differ on
# whether standalone math operators count: the project convention (operators excluded) and the
# worst-case whitespace count (operators included). Both must clear 350.
CVD_ABSTRACT_LIMIT = 350
_ab_src = rd(os.path.join(FIG, "ABSTRACT.md"))
_ab_src = re.sub(r"<!--.*?-->", "", _ab_src, flags=re.S)
_m_full = re.search(r"##\s*FULL abstract[^\n]*\n(.*?)\n---", _ab_src, flags=re.S)
if not _m_full:
    fail("GATE C [abstract] FULL abstract block not found in ABSTRACT.md")
else:
    _body = re.sub(r"^#.*$", "", _m_full.group(1), flags=re.M).replace("**", "").replace("---", " ")
    _worst = len(_body.split())
    if F["abstract_full_words"] > CVD_ABSTRACT_LIMIT:
        fail(f"GATE C [abstract over journal limit] {F['abstract_full_words']} words "
             f"> {CVD_ABSTRACT_LIMIT} (project convention)")
    if _worst > CVD_ABSTRACT_LIMIT:
        fail(f"GATE C [abstract over journal limit] {_worst} words > {CVD_ABSTRACT_LIMIT} "
             f"(worst-case count, math operators included)")
    _heads = re.findall(r"\*\*(\w+)\.\*\*", _m_full.group(1))
    if _heads != ["Background", "Methods", "Results", "Conclusions"]:
        fail(f"GATE C [abstract structure] CVD requires Background/Methods/Results/Conclusions, got {_heads}")
if F["abstract_cap_words"] > 250:
    note(f"abstract cap variant is {F['abstract_cap_words']} words, above a strict 250-word limit "
         f"(kept in the build source for a journal stricter than Cardiovascular Diabetology; "
         f"asserted below as NOT shipped)")

# ================= GATE F — the rendered panels must not contradict the text =================
# Added 2026-07-26. Every check below exists because the artifact ALREADY FAILED it while every
# earlier gate passed, and in each case only the prose half of a fix had landed. Gates A-E read
# .md/.docx/.xlsx; the figure PDFs were a blind spot, so a panel could assert what the text denied.
# These assert TRUTH about the shipped artifact, not the absence of a past mistake.

# The PDF reader is the SHARED one (~/.claude/lib/pdftext.py, on the Python path via a .pth), not a
# local copy. This file used to carry its own scanner, and on 2026-07-27 `pdftext.py --audit` found
# FIVE such copies across this project. Measuring them against the shared reader showed the local
# scanner was not blind but NOISY: it returned ~70,000 characters for a one-page cairo_pdf figure
# against ~1,000 from PyMuPDF, the extra 69,000 being CMap and font-encoding junk (CIDInit,
# beginbfchar, Adobe-Identity-UCS, glyph names). Every real panel string was recovered by BOTH, so
# no absence check was ever wrong — but a PRESENCE probe over 70,000 characters of junk is the
# inverted control CLAUDE.md warns about: the noisier the extractor, the likelier a probe passes by
# chance. Import raises loudly if the shared library is missing; a checker that cannot read must
# never report clean.
sys.path.insert(0, SHARED_LIB)
from pdftext import pdf_text

_PDFS = {rel(p): pdf_text(p) for p in ALL if p.lower().endswith(".pdf")}
if len(_PDFS) < 12:
    fail(f"GATE F [scope] only {len(_PDFS)} figure PDFs opened; expected at least 12")
_TXT = {rel(p): (text_of(p) or "") for p in ALL if text_of(p) is not None}
_SURF = dict(_TXT); _SURF.update(_PDFS)

# F1. The retired T2D mediated range. The covariance sweep replaced 68-79% with an envelope, and no
#     surface may quote the old pair. The PDF text layer encodes the en-dash as a literal \226, so a
#     "68.79%" regex could not match it - that is exactly how it survived in SupplFig7 undetected.
for k, t in _SURF.items():
    if re.search(r"68.{0,6}79\s*%", t):
        fail(f"GATE F [retired mediated range 68-79%] {k}")
# F2. Unqualified "DIRECT" as a panel label. MVMR identifies a component conditional on the modelled
#     exposures; the prose says "not mediated through CAD" and the panels must say the same thing.
for k, t in _PDFS.items():
    if re.search(r"%\s*DIRECT", t):
        fail(f"GATE F [unqualified '% DIRECT' panel label] {k}")
# F3. The synthesis panel's evidence tiers. Staging is bounded by a non-significant degree-preserving
#     null, so it cannot sit in the highest tier; only the BMI row may. Assert the COUNT, because the
#     row text had been corrected while the tier label beside it had not.
_s7 = _PDFS.get("04_supplementary/SupplFig7.pdf", "")
if not _s7:
    fail("GATE F [scope] SupplFig7.pdf text layer not readable")
else:
    # The expected count is DERIVED from the renderer's own tier vector, never typed here: a literal
    # would have to be edited by hand every time the synthesis strip is re-graded, and a hand-edited
    # expectation is how a gate stops asserting anything. fig6.R itself carries the one assertion
    # that is a claim about the evidence (the AHA-staging row may never be top-tier).
    _fig6 = rd(os.path.join(FIG, "fig6.R"))
    _mt = re.search(r"tier=c\(([^)]*)\)", _fig6)
    if not _mt:
        fail("GATE F [evidence tiers] could not read the tier vector from fig6.R")
    else:
        _tiers = re.findall(r'"([^"]+)"', _mt.group(1))
        _expect_hi = sum(1 for t in _tiers if t.startswith("HIGHER CONF"))
        _n_hi = len(re.findall(r"HIGHER CONF", _s7))
        if _n_hi != _expect_hi:
            fail(f"GATE F [evidence tiers] SupplFig7.pdf shows {_n_hi} HIGHER CONF. rows but fig6.R "
                 f"assigns {_expect_hi} — the rendered panel is stale or the renderer changed")
        # Do NOT re-implement "which row is staging" by position here — that assumption would go
        # stale silently the next time a row is added. fig6.R asserts it at render time and ABORTS
        # the render if it fails; this gate asserts only that the guard is still in the renderer.
        if not re.search(r'stopifnot\(tier\$tier\[grep\("AHA staging"', _fig6):
            fail("GATE F [evidence tiers] fig6.R has lost its render-time assertion that the "
                 "AHA-staging row is not in the top tier")
    if "not beyond a degree-preserving null" not in _s7:
        fail("GATE F [evidence tiers] SupplFig7 staging row must state the degree-preserving bound")
# F4. Panel labels that assert more than the test behind them supports.
for _pat, _why in [(r"resolved in Fig", "Figure 3 grades these edges, it does not resolve them"),
                   (r"T2D reaches HFrEF,\s*not HFpEF",
                    "asserts a subtype difference the formal test supports only nominally")]:
    for k, t in _SURF.items():
        if re.search(_pat, t):
            fail(f"GATE F [over-claiming label: {_why}] {k}")
# F5. One abstract ships, not two. The CAP variant's wording predates the round-3 recalibrations.
_ab_shipped = _TXT.get("01_manuscript/ABSTRACT.md", "")
if "CAP variant" in _ab_shipped or _ab_shipped.count("**Background.**") != 1:
    fail("GATE F [two abstracts] shipped ABSTRACT.md must carry exactly one structured abstract")
# F9. A bare "negative control". HDL->CAD is a CONFOUNDED correlated marker, not a null control, and
#     round 2 replaced the label in the prose, on Figure 3 and in the workbook. It survived on
#     SupplFig2 until 2026-07-26 because the project's PDF reader saw only `(...) Tj` operators and
#     that legend is drawn with `[...] TJ` arrays. The check is only worth having now that the
#     extractor reads both — which is the point: a gate is only as wide as what it can see.
for k, t in _SURF.items():
    for m in re.finditer(r"negative[\s-]control", t, re.I):
        _ctx = t[max(0, m.start() - 40):m.start()]
        if not re.search(r"correlated.marker\s*$", _ctx, re.I):
            fail(f"GATE F [bare 'negative control'] {k}: ...{' '.join(_ctx.split())[-45:]}"
                 f"{t[m.start():m.end()]}...")

# F7. The size of the Steiger-directed graph, DERIVED from the file that draws it. This gate exists
#     because "54-edge directed causal graph" shipped in v4 Results AND in the STROBE checklist long
#     after the strict-SBP promotion took the graph to 53. Every gate passed, because no check ever
#     asked the source file how many edges it contains. A number in prose is a claim about a file.
_edges_csv = os.path.join(PKG, "05_source_data", "fig1_edges.csv")
if not os.path.exists(_edges_csv):
    fail("GATE F [graph size] 05_source_data/fig1_edges.csv missing; cannot derive the graph size")
else:
    import csv as _csv
    with open(_edges_csv, encoding="utf-8") as _fh:
        _n_graph = sum(1 for _ in _csv.DictReader(_fh))
    # Search FORWARD for the phrase, then take the NEAREST preceding "<n>-edge". A forward regex
    # from the number is wrong and was caught in testing: "the 132-edge network, the 53-edge
    # Bonferroni-significant Steiger-directed causal graph" made the gate read 132 as the graph size
    # and fail a sentence that is correct. Verify the verifier.
    _seen = False
    for k, t in _SURF.items():
        for m in re.finditer(r"(?:Steiger-directed|directed causal graph)", t):
            _pre = t[max(0, m.start() - 80): m.start()]
            _nums = re.findall(r"(\d+)[- ]edge", _pre)
            if not _nums:
                continue
            _seen = True
            if int(_nums[-1]) != _n_graph:
                fail(f"GATE F [graph size] {k} says a {_nums[-1]}-edge Steiger-directed graph; "
                     f"fig1_edges.csv contains {_n_graph}")
    if not _seen:
        fail(f"GATE F [graph size] no surface states the size of the Steiger-directed graph "
             f"({_n_graph} edges) — the claim was removed rather than corrected")
    # and the Bonferroni count, derived the same way from derive_facts rather than typed
    _n_bonf = F["edges_passing_bonferroni"]
    if str(_n_bonf) not in _TXT.get("01_manuscript/RESULTS.md", ""):
        fail(f"GATE F [Bonferroni count] RESULTS must state the derived count {_n_bonf}")

# F8. Supplementary Methods must reach BOTH shipped twins, and must not be duplicated in the main
#     Methods. The move was deferred for exactly this reason: SUPPLEMENTARY_INFORMATION.md is written
#     from source while the .docx is assembled block by block, so text can land in one and not the
#     other. Counts are DERIVED from the source, never typed.
_sm_src = rd(os.path.join(FIG, "SUPPL_METHODS.md"))
_n_src = len(re.findall(r"^## Supplementary Methods ", _sm_src, flags=re.M))
if _n_src < 15:
    fail(f"GATE F [suppl methods] source has only {_n_src} sections")
else:
    for _surf in ("04_supplementary/SUPPLEMENTARY_INFORMATION.md",
                  "04_supplementary/SUPPLEMENTARY_INFORMATION.docx"):
        _n = len(re.findall(r"Supplementary Methods \d+ \|", _TXT.get(_surf, "")))
        if _n != _n_src:
            fail(f"GATE F [suppl methods] {_surf} carries {_n} sections; source has {_n_src} — "
                 f"the two shipped twins disagree")
    # nothing moved may still stand in the main Methods: assert on the first clause of each block
    _meth_txt = _TXT.get("01_manuscript/METHODS.md", "")
    for _m in re.finditer(r"^## Supplementary Methods [^\n]*\n\n(.{40,110}?)[,.;]", _sm_src,
                          flags=re.M | re.S):
        _frag = re.sub(r"\s+", " ", _m.group(1)).strip()
        if _frag and _frag in re.sub(r"\s+", " ", _meth_txt):
            fail(f"GATE F [suppl methods duplicated in Methods] '{_frag[:60]}...'")
    # and every Supplementary Methods number the Methods points at must exist
    for _n in sorted({int(x) for x in re.findall(r"Supplementary Methods (\d+)", _meth_txt)}):
        if not re.search(rf"^## Supplementary Methods {_n} \|", _sm_src, flags=re.M):
            fail(f"GATE F [dangling pointer] Methods cites Supplementary Methods {_n}, which does not exist")

# F10. Round-4 Tier 1: staging places traits in earlier clinical STAGES; it does not establish that
#      they are causally upstream. The phrase was corrected once in RESEARCH_INSIGHTS.md and came
#      back on the next build, because that file is GENERATED by build_research_insights.py and the
#      edit had been made to the artifact rather than to its generator. This gate catches that class.
for k, t in _SURF.items():
    if re.search(r"causally upstream", t, re.I):
        fail(f"GATE F [staging overclaim] {k}: 'causally upstream' — staging places traits in earlier "
             f"clinical stages; it does not prove causality")

# ---- F11-F15: round-4 Tier 2. Each is a CLASS the reviewer found one instance of; every one of the
#      five turned out to have instances the review did not name, in files the fix list did not list.
#      Reader-facing = what the journal receives. README_SUBMISSION.md, _author_todo/ and the
#      checksum/source-data manifests are ours; they may name build paths.
_READER = {k: t for k, t in _SURF.items()
           if k.startswith(("01_manuscript/", "02_cover_declarations/", "04_supplementary/",
                            "03_main_figures/"))}
if len(_READER) < 30:
    fail(f"GATE F [scope] only {len(_READER)} reader-facing surfaces opened for F11-F15; expected 30+")

# F11. Internal typography on a shipped panel. Round 4 asked for "beta_iv" and the ASCII arrows on
#      four SUPPLEMENTARY figures; the PDF text layers showed the same class on MAIN Figures 1 and 4,
#      which the review had not opened. Edge keys stay ASCII in the CSVs and are typeset via
#      fig_setup.R::disp_edge at plot time, so this asserts the DISPLAY, not the data.
#      NB: the arrow glyph itself is not recoverable from these subsetted fonts, so this is an
#      absence check with no presence twin - hence the char-count floor below, which is the only
#      thing standing between "clean" and "the reader went blind".
for k, t in _PDFS.items():
    # A blind reader passes every absence check, and a FLOODED one passes every presence probe, so
    # the plausible band is asserted at BOTH ends. Measured on the shipped set with the shared
    # reader: 389-1,336 characters for a one-page composite figure. The floor catches an unreadable
    # panel; the ceiling catches the noise mode that made the retired local scanner return ~70,000.
    if not 200 <= len(t) <= 6000:
        fail(f"GATE F [scope] {k} yielded {len(t)} characters, outside the plausible 200-6,000 band "
             f"for a one-page figure; this reader is either blind or flooded, so no absence or "
             f"presence claim over it can be trusted")
    for _m in re.finditer(r"\b[a-z]+_[a-z]+\b", t):
        fail(f"GATE F [variable name on a panel] {k}: '{_m.group(0)}'")
# The ASCII arrow is banned on EVERY reader-facing surface, not only the panels: 259 cells of the
# Supplementary Tables workbook (Additional file 3, which reviewers read as tables) carried it while
# the review had named only four supplementary figures. The machine-readable copies keep the ASCII
# key and ship separately as Additional file 5, so 05_source_data is deliberately outside _READER.
for k, t in _READER.items():
    for _m in re.finditer(r"[A-Za-z0-9]+->[A-Za-z0-9]+", t):
        fail(f"GATE F [ASCII arrow in reader-facing text] {k}: '{_m.group(0)}' — typeset it")
# PRESENCE TWIN. The rule above is an absence check, and the shared reader decodes U+2192, so the
# panels that carry edge labels must be shown to actually carry the TYPESET arrow. Without this, a
# renderer that dropped its edge labels entirely would report cleaner than one that kept them.
for _k in ("03_main_figures/Figure1.pdf", "03_main_figures/Figure4.pdf",
           "04_supplementary/SupplFig2.pdf", "04_supplementary/SupplFig3.pdf",
           "04_supplementary/SupplFig4.pdf", "04_supplementary/SupplFig5.pdf"):
    _t = _PDFS.get(_k)
    if _t is None:
        fail(f"GATE F [scope] {_k} was not opened; the arrow presence twin cannot run")
    elif "→" not in _t:
        fail(f"GATE F [typeset arrow missing] {_k} carries no '→'; its edge labels are gone, not fixed")

# F12. Revision-history language. A first submission must not narrate its own build history. The
#      review named one clause in METHODS; two more were shipping, one of them in a supplementary
#      figure legend, one row away from the sentence that had been fixed.
for k, t in _READER.items():
    for _m in re.finditer(r"(?:an?|the)\s+(?:earlier|previous|prior)\s+(?:build|draft|version)"
                          r"|\bearlier build\b", t, re.I):
        fail(f"GATE F [revision-history language] {k}: '{_m.group(0)}'")

# F13. Figure 4c must PRINT the fitted global offset, and it must be the number the fit produced.
#      This is a PRESENCE check on purpose: F11 above is pure absence, and an absence-only gate over
#      a PDF cannot tell a clean panel from an unreadable one.
_gs_rows = list(csv.DictReader(open(os.path.join(RES, "interaction_global_scale.csv"),
                                    encoding="utf-8")))
_slopes = {r["global_slope"] for r in _gs_rows}
if len(_slopes) != 1:
    fail(f"GATE F [global offset] interaction_global_scale.csv carries {len(_slopes)} distinct slopes")
else:
    _slope_s = f"{float(_slopes.pop()):.3f}"
    _f4 = _PDFS.get("03_main_figures/Figure4.pdf", "")
    if _slope_s not in _f4:
        fail(f"GATE F [global offset] Figure4.pdf does not print the fitted slope {_slope_s}; the "
             f"panel's solid line is unlabelled or the label was typed rather than derived")
    _s6leg = _TXT.get("04_supplementary/SUPPLEMENTARY_INFORMATION.md", "")
    if _slope_s not in _s6leg:
        fail(f"GATE F [global offset] the Supplementary Figure S6 legend must report {_slope_s}, the "
             f"same fit Figure 4c draws")

# F14. Every figure legend, main and supplementary, within the 300-word limit. The main legends were
#      trimmed in v4.3 and S6 was left at 419 words.
def _legend_words(txt, pat):
    out, hits = {}, [(m.start(), m.group(1)) for m in re.finditer(pat, txt)]
    for i, (s, name) in enumerate(hits):
        e = hits[i + 1][0] if i + 1 < len(hits) else len(txt)
        out[name] = len(txt[s:e].split())
    return out
_LEGEND_CAP = 300
_sfw = _legend_words(_TXT.get("04_supplementary/SUPPLEMENTARY_INFORMATION.md", ""),
                     r"\*\*(Supplementary Figure S\d+) \|")
if len(_sfw) < 7:
    fail(f"GATE F [scope] only {len(_sfw)} supplementary legends found; expected 7")
_mfw = _legend_words(_TXT.get("01_manuscript/FIGURE_LEGENDS.md", ""), r"\*\*(Figure \d+) \|")
if len(_mfw) < 4:
    fail(f"GATE F [scope] only {len(_mfw)} main legends found; expected 4")
for _n, _w in sorted({**_mfw, **_sfw}.items()):
    if _w > _LEGEND_CAP:
        fail(f"GATE F [legend length] {_n} legend is {_w} words, over the {_LEGEND_CAP}-word limit")

# F15. Build-tree paths in reader-facing text. The criterion is REACHABILITY: no surface the journal
#      receives may name a location the reader cannot open. `05_source_data/forward_local_edges.csv`
#      is the package's own internal layout (the reader is given "Additional file 5"), and
#      `results/` is our build tree. A BARE filename is allowed — those name artifacts in the
#      deposited repository that Code availability points at — so the pattern requires a directory
#      component. Written this way the gate also caught the Supplementary Tables Contents column,
#      which was headed "Source (results/)" and which no review round had opened.
for k, t in _READER.items():
    for _m in re.finditer(r"0[1-5]_[a-z_]+/[\w./-]*|(?<![\w-])(?:scripts|figures|results|panels)/[\w./-]+"
                          r"|independent_build|suppfig_data/", t):
        fail(f"GATE F [build-tree path in reader-facing text] {k}: '{_m.group(0)}'")

# F16. Raw markdown syntax reaching a shipped .docx. The Additional-files table shipped in the
#      submitted manuscript as seven literal paragraphs, separator row included, because the
#      Declarations loop in build_manuscript_docx.py had no table branch while the packager's
#      renderer did. Nothing had ever looked at the .docx as a READER sees it.
for _r in sorted(k for k in _READER if k.endswith(".docx")):
    _t = _TXT.get(_r, "")
    for _ln in _t.split("\n"):
        _s = _ln.strip()
        if _s.startswith("|") or "|---" in _s:
            fail(f"GATE F [raw markdown in a .docx] {_r}: '{_s[:60]}' — render it as a Word table")
    for _m in re.finditer(r"(?<!\*)\*\*[^*\n]{1,60}\*\*(?!\*)", _t):
        fail(f"GATE F [raw markdown in a .docx] {_r}: literal bold markers '{_m.group(0)[:40]}'")

# F17. THE CEILING STATISTIC: the degree-preserving edge-rewiring null bounds the paper's headline
#      staging claim, and Figure 1f printed P = 0.31 from a two-day-stale CSV while the Abstract,
#      Results, Discussion and the figure's own legend all said 0.27 — the value the analysis
#      produced. Four review rounds and every gate passed over it, because every check was an
#      ABSENCE check over prose and none compared a number on a PANEL to its source. Assert the
#      chain: canonical output == prose == what the figure prints.
_cn_txt = rd(os.path.join(RES, "staging_constrained_nulls.txt"))
_m_can = re.search(r"C\. Degree-preserving edge-rewiring null[^\n]*?P = ([\d.]+)", _cn_txt)
if not _m_can:
    fail("GATE F [ceiling statistic] cannot read the degree-preserving null P from "
         "staging_constrained_nulls.txt — the gate has no source to check against")
else:
    _p_can = float(_m_can.group(1))
    _shown = f"{_p_can:.2f}"                      # the panel and the prose both print 2 dp
    for _surf in ("01_manuscript/ABSTRACT.md", "01_manuscript/RESULTS.md",
                  "01_manuscript/DISCUSSION.md", "01_manuscript/FIGURE_LEGENDS.md"):
        _t = _TXT.get(_surf, "")
        _ps = re.findall(r"degree-preserving[^.]{0,120}?\(\*?P\*? = ([\d.]+)\)", _t)
        for _p in _ps:
            if abs(float(_p) - _p_can) > 5e-3:
                fail(f"GATE F [ceiling statistic] {_surf} states degree-preserving P = {_p}, but "
                     f"staging_constrained_nulls.txt gives {_p_can}")
    _f1 = _PDFS.get("03_main_figures/Figure1.pdf", "")
    if not re.search(rf"P\s*=\s*{re.escape(_shown)}\b", _f1):
        fail(f"GATE F [ceiling statistic] Figure1.pdf does not print P = {_shown} for the "
             f"degree-preserving null; the panel and the prose disagree")

# F18. Truncated Supplementary Methods. The 20-section move relocated blocks verbatim, and because
#      SUPPL_METHODS.md is not in lib_vancouver.REWRITE a citation could not travel with them — so
#      two sentences were cut at their citation and left dangling. Section 5 ended at "A
#      degree-preserving", which is the DEFINITION of the null the Abstract, Results, Discussion and
#      Figure 1f all rest on; section 13 ended at "the conventional threshold of 10". Both shipped in
#      the .md and the .docx. Assert every section ends in a terminal stop.
_si = _TXT.get("04_supplementary/SUPPLEMENTARY_INFORMATION.md", "")
_sm_secs = re.findall(r"^## (Supplementary Methods \d+)[^\n]*\n\n([^\n]+)", _si, flags=re.M)
if len(_sm_secs) < 20:
    fail(f"GATE F [scope] only {len(_sm_secs)} Supplementary Methods sections read for the "
         f"truncation check; expected at least 20")
for _name, _body in _sm_secs:
    if _body.rstrip()[-1:] not in ".!?":
        fail(f"GATE F [truncated section] {_name} ends mid-sentence: "
             f"'...{_body.rstrip()[-55:]}'")

# F19. DIALECT CONSISTENCY. Never checked on this paper until 2026-07-27, because both local prose
#      scanners were forks of one another and neither carried the spelling lens the shared
#      ~/.claude/lib/prosecheck.py has. The manuscript is British (231 tokens vs 32); the minority
#      forms turned out to be proper nouns, CRediT terms and verbatim reference titles — all correct
#      — except "labeling/labelings", which sat beside "labelled" x4 in the same package.
#      Only TRUE variant pairs are listed, so "size", "arise", "imprecise" and the plural "analyses"
#      cannot pollute the count, and REFERENCES_NUMBERED.md is excluded because a citation title is
#      quoted verbatim and must never be re-spelled.
_DIALECT = [("harmonis", "harmoniz"), ("normalis", "normaliz"), ("summaris", "summariz"),
            ("prioritis", "prioritiz"), ("formalis", "formaliz"), ("generalis", "generaliz"),
            ("standardis", "standardiz"), ("recognis", "recogniz"), ("organis", "organiz"),
            ("minimis", "minimiz"), ("maximis", "maximiz"), ("categoris", "categoriz"),
            ("analys", "analyz"), ("behaviour", "behavior"), ("colour", "color"),
            ("favour", "favor"), ("modelling", "modeling"), ("modelled", "modeled"),
            ("labelling", "labeling"), ("labelled", "labeled"), ("signalling", "signaling"),
            ("centred", "centered"), ("grey", "gray"), ("haemoglobin", "hemoglobin")]
_SUF = r"(?:e|es|ed|ing|ation|ations|s)?"
# The ae/oe families, added 2026-07-27 after an adversarial pass found ischaemic x1 against
# ischemic x5 — in the .docx the two spellings sat in the SAME sentence stem in adjacent sections.
# The earlier sweep missed it for one reason: ischaem/ischem was not in the list. Complete FILE
# coverage and a passing mutation test cannot rescue a PATTERN nobody wrote.
# These stems take adjectival/nominal endings the verb suffix group above does not allow
# ("ischaemic", "aetiology", "paediatric"), so they carry their own.
_DIALECT_AE = [("ischaem", "ischem"), ("aetiolog", "etiolog"), ("haemodynam", "hemodynam"),
               ("oedema", "edema"), ("oestrog", "estrog"), ("paediatr", "pediatr"),
               ("foetal", "fetal"), ("anaemi", "anemi"), ("haemostas", "hemostas")]
_SUF_AE = r"(?:ic|ics|ical|ically|y|ies|en|ens|a|as|is|ous|ed|es|s)?"
# NOT a generic -aemia/-emia rule: glycemia, dysglycemia and dyslipidemia are used consistently in
# the American form throughout as terms of art, and a blanket pattern would report them as defects.
# fixed field terms and proper nouns that are NOT dialect choices
_DIA_OK = re.compile(r"Mendelian[- ]randomi[sz]ation|GWAS Catalog|Conceptualization|"
                     r"Medical Cent(?:er|re)|Center for Bioenergetics|harmonise_data", re.I)
_br = _am = 0
_am_hits = []
for k, t in _READER.items():
    if "REFERENCES" in k:
        continue
    # The manuscript .docx CONTAINS the rendered reference list, so sweeping it whole re-imports
    # every citation title. Cut at the References heading: titles are quoted verbatim and a gate
    # that asked them to change dialect would be asking for a corrupted citation.
    if k.endswith("CKM_Paper4_manuscript.docx"):
        _cut = re.search(r"^References\s*$", t, flags=re.M)
        if not _cut:
            fail("GATE F [scope] manuscript .docx has no 'References' heading; the dialect sweep "
                 "cannot exclude the reference titles and would report false positives")
            continue
        t = t[:_cut.start()]
    for _tbl, _suf in ((_DIALECT, _SUF), (_DIALECT_AE, _SUF_AE)):
        for _b, _a in _tbl:
            for _m in re.finditer(rf"\b{_b}{_suf}\b", t, re.I):
                if not _DIA_OK.search(t[max(0, _m.start() - 30):_m.end() + 30]):
                    _br += 1
            for _m in re.finditer(rf"\b{_a}{_suf}\b", t, re.I):
                if not _DIA_OK.search(t[max(0, _m.start() - 30):_m.end() + 30]):
                    _am += 1; _am_hits.append((k, _m.group(0)))
if _br + _am < 50:
    fail(f"GATE F [scope] dialect sweep saw only {_br + _am} variant tokens; the manuscript has "
         f"hundreds, so this sweep is not reading the package")
elif _br and _am:
    fail(f"GATE F [mixed dialect] {_br} British vs {_am} American variant tokens; the minority set "
         f"is {sorted({h for _, h in _am_hits})[:8]}")

# F20. CROSS-REPRESENTATION. An artifact has representations, and a gate over one says nothing about
#      the others. Every .md/.docx twin must carry the same numbers, the same figure/table callouts
#      and the same citations — that is the boundary Paper A's reference-numbering blocker lived on,
#      and it is the boundary on which SUPPLEMENTARY_INFORMATION's twins were found to disagree
#      (the .docx told the reader where the tables were and the .md did not). Uses the SHARED
#      prosecheck.guard, not a local copy.
sys.path.insert(0, SHARED_LIB)
try:
    import prosecheck as _PC
except Exception as _e:                      # a missing checker must be loud, never silently skipped
    fail(f"GATE F [cross-representation] shared prosecheck unavailable ({_e}); the twin check "
         f"cannot run and must not be reported as clean")
    _PC = None
if _PC is not None:
    import io as _io, tempfile as _tf
    _twins = [(k, k[:-5] + ".md") for k in _SURF if k.endswith(".docx") and k[:-5] + ".md" in _TXT]
    if len(_twins) < 5:
        fail(f"GATE F [scope] only {len(_twins)} .md/.docx twin pairs found; expected at least 5")
    for _dx, _md in sorted(_twins):
        _a = os.path.join(_tf.gettempdir(), "gate_" + os.path.basename(_md))
        _b = os.path.join(_tf.gettempdir(), "gate_" + os.path.basename(_dx) + ".txt")
        _io.open(_a, "w", encoding="utf-8").write(_TXT[_md])
        _io.open(_b, "w", encoding="utf-8").write(_SURF[_dx])
        _sink = _io.StringIO()
        if not _PC.guard(_a, _b, out=_sink):
            _detail = " | ".join(l.strip() for l in _sink.getvalue().split("\n")
                                 if "DROPPED" in l or "ADDED" in l)
            fail(f"GATE F [cross-representation] {_md} and its .docx twin disagree: {_detail[:220]}")

# F21. DANGLING MAIN-FIGURE CALLOUT. The STROBE checklist cited "Figures 2d–e, 3, 6b–d" while the
#      package ships four main figures: migrate_fig56_to_supp.py maps main 5/6 to Supplementary
#      S6/S7, but its XREF regex requires the digit to follow "Figure(s)" directly, and here the
#      "Figures" prefix had already been consumed by "2d–e" in the same comma list. Its survivor
#      check used the SAME regex, so it could not see what the regex could not match. Assert the
#      truth instead: every main-figure number a reader-facing surface cites must exist.
_n_main = F["main_figures"]
for k, t in _READER.items():
    # PANELS ARE IN SCOPE. This rule previously skipped .pdf, and that is exactly how Supplementary
    # Figure S2 shipped an on-panel "(attenuates in MVMR, Fig 6)" through five rounds: the migration
    # that renamed main Figure 6 rewrote .md sources only, and the gate that catches dangling
    # main-figure callouts declined to read the one surface the stale reference was on.
    for _m in re.finditer(r"\bFigs?\.?\s+((?:\d+[a-h]?(?:\s*[,–-]\s*[a-h])?(?:\s*(?:,|and)\s*)?)+)"
                          r"|\bFigures?\s+((?:\d+[a-h]?(?:\s*[,–-]\s*[a-h])?(?:\s*(?:,|and)\s*)?)+)", t):
        for _n in re.findall(r"\b(\d+)", _m.group(1) or _m.group(2) or ""):
            if int(_n) > _n_main:
                fail(f"GATE F [dangling figure callout] {k} cites main Figure {_n}, but the package "
                     f"ships {_n_main} main figures: '{' '.join(_m.group(0).split())[:60]}'")

# F23. The journal-facing TITLE PAGE carries no internal production metrics and no pointer to an
#      editable source file. It was shipping a word-count inventory, a "287 words (short variant)"
#      naming an abstract the package does not ship, and "see ABSTRACT.md".
_tp = _TXT.get("02_cover_declarations/TITLE_PAGE.md", "")
if not _tp:
    fail("GATE F [scope] TITLE_PAGE.md was not read; the title-page rule cannot run")
for _bad in ("Manuscript metrics", "short variant", "ABSTRACT.md", "Main text total",
             "Main display items"):
    if _bad in _tp:
        fail(f"GATE F [title page] carries internal production detail: '{_bad}'")

# F24. Each declared Additional file must be ONE physical upload, sequentially named. The inventory
#      declared "Additional file 4 | PDF" for seven SupplFig PDFs and "Additional file 5 | CSV" for
#      forty-seven CSVs — many files described as one, which the journal's instructions do not allow.
_decl = _TXT.get("02_cover_declarations/DECLARATIONS.md", "")
_af = re.findall(r"\| (Additional file (\d+)) \| (\w+) \|([^|]*)\|", _decl)
if len(_af) < 3:
    fail(f"GATE F [scope] only {len(_af)} Additional-file rows parsed from DECLARATIONS.md")
else:
    _nums = [int(n) for _, n, _f, _d in _af]
    if _nums != list(range(1, len(_nums) + 1)):
        fail(f"GATE F [additional files] not numbered sequentially from 1: {_nums}")
    _EXPECT = {1: "02_cover_declarations/STROBE_MR_CHECKLIST.docx",
               2: "04_supplementary/SUPPLEMENTARY_INFORMATION.docx",
               3: "04_supplementary/Supplementary_Tables.xlsx",
               4: "04_supplementary/Source_Data.zip"}
    for _, _n, _fmt, _ in _af:
        _want = _EXPECT.get(int(_n))
        if _want is None:
            fail(f"GATE F [additional files] Additional file {_n} is declared but has no known "
                 f"physical upload")
        elif not os.path.exists(os.path.join(PKG, _want.replace("/", os.sep))):
            fail(f"GATE F [additional files] Additional file {_n} is declared but {_want} is not "
                 f"in the package")
        elif _want.rsplit(".", 1)[1].upper() != _fmt.upper():
            fail(f"GATE F [additional files] Additional file {_n} is declared {_fmt} but the "
                 f"physical upload is {_want.rsplit('.', 1)[1].upper()}")
    _cited = {int(x) for t in _READER.values() for x in re.findall(r"Additional files? (\d+)", t)}
    if _cited - set(_nums):
        fail(f"GATE F [additional files] the text cites Additional file(s) {sorted(_cited - set(_nums))} "
             f"which the inventory does not declare")

# F22. THE GRAPHICAL ABSTRACT MUST AGREE WITH THE MANUSCRIPT. It printed the BMI→HF direct effect as
#      +0.41 while the Abstract, Results, Discussion and Figure 2 all printed +0.42, through every
#      review round. Cause: it read data/fig2_mediation.csv, which prep_newpanels.py writes at the
#      3 decimals results/formal_mediation.txt prints, and 0.415 as a double is 0.41499999999999998,
#      so "%+.2f" rounded down. "Read not typed" was true and still not enough — the value came from
#      a lossier COPY of the same quantity. Assert the printed label against the canonical source.
_ga = _PDFS.get("03_main_figures/GraphicalAbstract.pdf", "")
_mvc = [r for r in csv.DictReader(open(os.path.join(RES, "mvmr_cascade.csv"), encoding="utf-8"))
        if r.get("model") == "HF_via_CAD" and r.get("exposure") == "BMI" and r.get("outcome") == "HF"]
if len(_mvc) != 1:
    fail(f"GATE F [graphical abstract] expected exactly one BMI→HF row in mvmr_cascade.csv, "
         f"found {len(_mvc)}")
else:
    _want = f"{float(_mvc[0]['direct_b']):+.2f}"
    if _want not in _ga:
        fail(f"GATE F [graphical abstract] does not print the BMI→HF direct effect {_want} that "
             f"mvmr_cascade.csv gives; it is showing a differently-rounded copy")
    if _want not in _TXT.get("01_manuscript/ABSTRACT.md", ""):
        fail(f"GATE F [graphical abstract] the Abstract does not state {_want}, so the figure and "
             f"the manuscript cannot both be right")

# F6. Cardiovascular Diabetology requires LLM use to be documented in Methods. Assert the section is
#     present and names both the assistance and the author's responsibility - not merely the word.
_meth = _TXT.get("01_manuscript/METHODS.md", "")
if not re.search(r"large language model", _meth, re.I):
    fail("GATE F [LLM disclosure] METHODS must document any LLM use (journal requirement)")
elif not re.search(r"full responsibility", _meth, re.I):
    fail("GATE F [LLM disclosure] the disclosure must affirm author responsibility")
if NOTES:
    print("\nLength notes (author decision, not gate failures):")
    for n_ in NOTES:
        print("  !", n_)

if FAILS:
    print(f"\n===== GATE FAILED — {len(FAILS)} problem(s) =====")
    for f_ in FAILS:
        print("  x", f_)
    sys.exit(1)
print(f"\n===== GATE PASSED — Gates A, B, C, D, E, F clean over all {len(ALL)} files "
      f"({len(_PDFS)} of them figure PDFs opened by Gate F) =====")
