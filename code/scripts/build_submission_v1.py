# -*- coding: utf-8 -*-
"""Cut the CKM Paper 4 FROZEN submission package v1.

Differs from the working-tree builder (build_submission_package.py) in four ways that matter:

  1. TWO-FOLDER SEPARATION. Writes to "paper 4/submission package v1/", outside independent_build,
     so the upload set is never the working tree.
  2. NOTHING IS HAND-TYPED. Every count in the title page, README and manifest comes from
     derive_facts.py. The AUTO:METRICS block in FRONT_MATTER.md is filled here, not maintained by hand.
  3. DELIVERABLES CARRY NO INTERNAL PROVENANCE. Shipped .md are stripped of HTML comments and
     whole-line italic meta notes (the same rule the docx renderer applies), so workflow IDs and
     verification notes cannot reach a reviewer.
  4. IT REFUSES TO SHIP STALE OR INCOMPLETE CONTENT. Generated artifacts must postdate their inputs,
     and the source-data set is derived from what the figure scripts actually read.

Author/creator metadata = Bertrand Chin-Ming Tan (global rule 8).
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

import os, re, sys, csv, json, glob, shutil, hashlib, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_vancouver as VAN
import lib_docx
from datetime import date, datetime
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

BASE = P4_BASE
ROOT = P4_ROOT
FIG = os.path.join(BASE, "figures")
MAN = os.path.join(BASE, "manuscript")
RES = os.path.join(BASE, "results")
SD_FIG = os.path.join(FIG, "suppfig_data")
SD_F1 = os.path.join(FIG, "data")
# Package version. The script name is historical (it cut v1); the version it CUTS is set here.
# v1 = the 111-edge pre-peer-review package, now SUPERSEDED and deliberately left untouched on disk
# as a historical record. v2 = the 132-edge rewrite. Because PKG is derived from this constant, the
# clear-and-rebuild block below can only ever empty the directory of the version being cut.
PKG_VERSION = "v4.3"
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")

PROBLEMS = []
def require(cond, msg):
    if not cond:
        PROBLEMS.append(msg)

def rd(p):
    return open(p, encoding="utf-8").read() if os.path.exists(p) else ""

def md5(p, chunk=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

# ============================ 0. canonical facts ============================
subprocess.run([sys.executable, os.path.join(BASE, "scripts", "derive_facts.py")],
               capture_output=True, text=True, check=True)
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))
VAN_MASTER, VAN_ORDER, VAN_AMAP = VAN.numbering()

# ============================ 1. staleness gate =============================
# copy2 preserves mtime, so mtime tells us nothing about COPIES; it is still the right signal for
# GENERATED artifacts, which must not predate the sources they were generated from.
SECTION_MD = ["INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION", "CONCLUSIONS", "ABBREVIATIONS",
              "ABSTRACT", "FIGURE_LEGENDS", "FRONT_MATTER", "REFERENCES_MASTER"]
docx_path = os.path.join(MAN, "CKM_Paper4_manuscript.docx")
newest_src, newest_name = 0, ""
for n in SECTION_MD:
    p = os.path.join(FIG, f"{n}.md")
    if os.path.exists(p) and os.path.getmtime(p) > newest_src:
        newest_src, newest_name = os.path.getmtime(p), f"{n}.md"
require(os.path.exists(docx_path), "manuscript docx missing")
if os.path.exists(docx_path):
    require(os.path.getmtime(docx_path) >= newest_src,
            f"STALE: CKM_Paper4_manuscript.docx predates {newest_name} — "
            f"re-run build_manuscript_docx.py before packaging")

xlsx_path = os.path.join(MAN, "Supplementary_Tables.xlsx")
require(os.path.exists(xlsx_path), "Supplementary_Tables.xlsx missing")

# ---- 1b. FIGURES vs the data they are drawn from -------------------------------------------
# The gate above covered exactly one generated artifact. Every FIGURE is generated too, and on
# 2026-07-24 Figure 6 shipped a retired staging P because it was rendered at 10:47 from a CSV that
# was rewritten at 13:09 — no gate compared the two. A rendered panel is a claim; if its input moved
# after it was drawn, the claim is unverified. Each fig*.R is parsed for the data files it reads and
# every render must post-date all of them.
_fig_stale = []
for _rs in sorted(glob.glob(os.path.join(FIG, "fig*.R")) + glob.glob(os.path.join(FIG, "suppfig*.R"))):
    _base = os.path.splitext(os.path.basename(_rs))[0]
    if _base == "fig_setup":
        continue
    _src = rd(_rs)
    # 2026-07-26: output names are DERIVED from each script's own save_fig() call. The old rule
    # inferred them from the FILENAME (fig<N>.R -> Figure<N>), which silently stopped matching when
    # the v3.2 split made fig5.R render SupplFig6 and fig6.R render SupplFig7: `Figure5.pdf` and
    # `Figure6.pdf` no longer exist, so `_outs` came out empty and BOTH renderers were skipped by
    # this gate entirely — including the synthesis panel that changed twice today.
    _names = re.findall(r'save_fig\([^,]+,\s*"([^"]+)"', _src)
    _outs = [os.path.join(FIG, f"{n}.{e}") for n in _names for e in ("pdf", "png")]
    _outs = [o for o in _outs if os.path.exists(o)]
    if not _outs:
        continue                      # script with no shipped render (e.g. a retired supp figure)
    _omt = min(os.path.getmtime(o) for o in _outs)
    # A render must post-date the CODE that draws it, not only the data it reads. Editing a label in
    # a .R and forgetting to re-render ships the old label with no trace in the copy.
    if os.path.getmtime(_rs) > _omt + 1:
        _fig_stale.append(f"{os.path.basename(_outs[0])} predates its renderer {os.path.basename(_rs)}")
    for _f in sorted(set(re.findall(r"[\"']([A-Za-z0-9_./-]+\.(?:csv|txt|tsv))[\"']", _src))):
        for _cand in (os.path.join(FIG, "data", os.path.basename(_f)),
                      os.path.join(FIG, "suppfig_data", os.path.basename(_f)),
                      os.path.join(BASE, "results", os.path.basename(_f))):
            if os.path.exists(_cand):
                if os.path.getmtime(_cand) > _omt + 1:      # 1 s slack for same-batch writes
                    _fig_stale.append(f"{os.path.basename(_outs[0])} predates {os.path.basename(_f)}")
                break
require(not _fig_stale,
        "STALE FIGURE(S) — a render is older than data it reads; re-run the figure script(s):\n    "
        + "\n    ".join(_fig_stale))

# ============================ 2. shipping-clean markdown ====================
def clean_md(t):
    """Strip internal provenance so it cannot reach a reviewer: HTML comments and whole-line
    italic meta notes. Mirrors build_manuscript_docx.add_body, so the shipped .md and the .docx
    are rendered from the same content."""
    t = re.sub(r"<!--.*?-->", "", t, flags=re.S)
    t = "\n".join(ln for ln in t.split("\n") if not re.match(r"^_.+_$", ln.rstrip()))
    return re.sub(r"\n{3,}", "\n\n", t).strip() + "\n"

# ---- fill the AUTO:METRICS block in the front matter (derived, never hand-typed) ----
def metrics_block():
    w = F["word_counts"]
    return "\n".join([
        f"- Abstract: {F['abstract_full_words']} words (shipped, against the Cardiovascular "
        f"Diabetology 350-word limit) / {F['abstract_cap_words']} words (short variant) — see ABSTRACT.md",
        f"- Introduction: {w['INTRODUCTION']:,} words",
        f"- Methods: {w['METHODS']:,} words",
        f"- Results: {w['RESULTS']:,} words",
        f"- Discussion: {w['DISCUSSION']:,} words",
        f"- Conclusions: {w['CONCLUSIONS']:,} words",
        f"- Main text total: {F['main_text_words']:,} words",
        f"- Main display items: {F['main_figures']} figures ({F['total_panels']} panels)",
        f"- Supplementary: {F['supp_tables']} tables (S1–S{F['supp_tables']}) in "
        f"{F['supp_table_sheets']} worksheets; {F['supp_figures']} supplementary figures",
        f"- References: {F['references_cited']}",
    ])

front_raw = rd(os.path.join(FIG, "FRONT_MATTER.md"))
require("[AUTO]" in front_raw or "AUTO:METRICS" in front_raw,
        "FRONT_MATTER.md has no AUTO:METRICS block — metrics would be hand-typed")

def fill_metrics(text):
    """Replace the metrics bullet list line-wise. Deliberately NOT a DOTALL regex: `(?:- .*\\n?)+`
    under re.S matches across newlines and eats the rest of the file, which silently emptied the
    Declarations section and shipped a 0-paragraph DECLARATIONS.docx."""
    lines = text.split("\n")
    try:
        i = next(k for k, l in enumerate(lines) if l.startswith("## Manuscript metrics"))
    except StopIteration:
        return text, False
    j = i + 1
    if j < len(lines) and lines[j].lstrip().startswith("<!--"):      # skip the guard comment
        while j < len(lines) and "-->" not in lines[j]:
            j += 1
        j += 1
    start = j
    while j < len(lines) and lines[j].startswith("- "):              # the bullet list only
        j += 1
    if j == start:
        return text, False
    return "\n".join(lines[:start] + metrics_block().split("\n") + lines[j:]), True

front_filled, ok_metrics = fill_metrics(front_raw)
require(ok_metrics, "AUTO:METRICS block not found or empty in FRONT_MATTER.md")
require("[AUTO]" not in front_filled, "AUTO:METRICS substitution failed — [AUTO] placeholders remain")
require("# Declarations" in front_filled,
        "metrics substitution destroyed the Declarations section")

# ---- title, taken from ONE source and asserted against the docx builder ----
mt = re.search(r"##\s*Title\s*\n\*\*(.+?)\*\*", front_raw, flags=re.S)
require(bool(mt), "could not read the title from FRONT_MATTER.md")
TITLE = " ".join(mt.group(1).split()) if mt else ""
builder_src = rd(os.path.join(BASE, "scripts", "build_manuscript_docx.py"))
mb = re.search(r'title = \("(.+?)"\s*\n\s*"(.+?)"\)', builder_src, flags=re.S)
if mb:
    docx_title = " ".join((mb.group(1) + mb.group(2)).split())
    require(docx_title == TITLE,
            f"TITLE MISMATCH between FRONT_MATTER.md and build_manuscript_docx.py:\n"
            f"    front matter: {TITLE}\n    docx builder: {docx_title}")

# ============================ 3. folder tree ================================
SUB = ["01_manuscript", "02_cover_declarations", "03_main_figures",
       "04_supplementary", "05_source_data", "_author_todo"]
# A frozen cut is rebuilt whole, never patched in place. Clear the CONTENTS rather than the
# directory itself: on Windows an open Explorer window or a shell sitting in the folder holds a
# handle on the directory node, and rmtree(PKG) then dies at the final rmdir having already
# deleted everything inside it.
if os.path.isdir(PKG):
    for entry in os.listdir(PKG):
        p = os.path.join(PKG, entry)
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
for s in SUB:
    os.makedirs(os.path.join(PKG, s), exist_ok=True)

# ---------- docx helpers ----------
# Underscore italics need word boundaries, or filenames (forward_local_edges.csv,
# Supplementary_Tables.xlsx) would be parsed as emphasis and lose their underscores.
INLINE = re.compile(r"\*\*(.+?)\*\*|\*(.+?)\*|(?<![\w`])_([^_`\n]+)_(?![\w])")
def add_inline(p, text):
    pos = 0
    for m in INLINE.finditer(text):
        if m.start() > pos:
            p.add_run(text[pos:m.start()])
        if m.group(1) is not None:
            p.add_run(m.group(1)).bold = True
        elif m.group(2) is not None:
            p.add_run(m.group(2)).italic = True
        else:
            p.add_run(m.group(3)).italic = True
        pos = m.end()
    if pos < len(text):
        p.add_run(text[pos:])

def render_table(doc, rows):
    # ONE implementation, shared with build_manuscript_docx.py (scripts/lib_docx.py). This renderer
    # worked; the manuscript builder had no table branch at all and shipped raw pipe syntax.
    lib_docx.render_md_table(doc, rows, add_inline)

def render_md(doc, md):
    """Render markdown, joining hard-wrapped continuation lines into one paragraph.

    Rendering line-by-line breaks inline spans that wrap: the title is written as
    `**Adiposity ...\\nMendelian-randomization ...**`, and each half alone has an unmatched `**`,
    so the asterisks shipped as literal text on the title page.
    """
    lines = clean_md(md).split("\n"); i = 0
    BREAKS = ("# ", "## ", "- ", "* ", "|", "---")
    def is_break(s):
        return (not s) or s.startswith(("#", "|", "- ", "* ")) or set(s.strip()) <= set("-*_ ")
    while i < len(lines):
        s = lines[i].rstrip()
        if not s:
            i += 1; continue
        if s.lstrip().startswith("|"):
            tbl = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                tbl.append(lines[i]); i += 1
            render_table(doc, tbl); continue
        if s.startswith("## "):
            doc.add_heading(s[3:].strip(), 2); i += 1; continue
        if s.startswith("# "):
            doc.add_heading(s[2:].strip(), 1); i += 1; continue
        if set(s.strip()) <= set("-") and len(s.strip()) >= 3:      # horizontal rule
            i += 1; continue
        buf = [s]; i += 1
        if not s.startswith(("- ", "* ")):                          # bullets stay one per line
            while i < len(lines) and not is_break(lines[i].rstrip()):
                buf.append(lines[i].strip()); i += 1
        p = doc.add_paragraph(); add_inline(p, " ".join(buf))

def new_doc():
    d = Document(); n = d.styles["Normal"]
    n.font.name = "Calibri"; n.font.size = Pt(10.5)
    n.paragraph_format.space_after = Pt(6)
    return d

def finalize(d, path, title):
    cp = d.core_properties
    cp.author = "Bertrand Chin-Ming Tan"
    cp.last_modified_by = "Bertrand Chin-Ming Tan"
    cp.title = title
    cp.comments = ""
    # python-docx ships a template whose created/modified stamps are 2013-12-23T23:15:00Z. Left
    # alone, every document in the package advertises that date and a generator default.
    cp.created = cp.modified = datetime.now()
    d.save(path)

# ============================ 4. 01_manuscript ==============================
shutil.copy2(docx_path, os.path.join(PKG, "01_manuscript", "CKM_Paper4_manuscript.docx"))
for f in ["INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md",
          # Cardiovascular Diabetology requires both of these as their own sections.
          "CONCLUSIONS.md", "ABBREVIATIONS.md",
          # Cardiovascular Diabetology asks for this table below the abstract, uploaded as a table.
          "RESEARCH_INSIGHTS.md",
          "ABSTRACT.md", "FIGURE_LEGENDS.md",
          # The manuscript's References preamble names this file, so it must travel with the
          # manuscript: a reviewer told a numbered list "accompanies this manuscript" must find it.
          "REFERENCES_NUMBERED.md"]:
    _src = rd(os.path.join(FIG, f))
    # 2026-07-26: ABSTRACT.md used to ship BOTH abstracts. The .docx carries only the FULL one, but
    # the shipped .md also carried the 286-word CAP variant, whose wording predates the round-3
    # recalibrations ("corroborating the order ... rather than proving it", "~80% direct") — so the
    # package contained two abstracts making differently-calibrated claims about the same result,
    # and an editor opening the .md would read the retired one. One journal is chosen; only its
    # abstract ships. The variant stays in the SOURCE for a future stricter-cap submission.
    if f == "ABSTRACT.md":
        _cut = _src.split("## CAP variant")
        require(len(_cut) == 2, "ABSTRACT.md: expected exactly one '## CAP variant' header")
        _src = _cut[0].rstrip().rstrip("-").rstrip() + "\n"
        require("Background." in _src and _src.count("**Background.**") == 1,
                "ABSTRACT.md: shipped copy must contain exactly one structured abstract")
    # Same build-time author-date -> Vancouver transform the .docx builder applies, from the same
    # module, so the shipped .md and the shipped .docx cannot cite differently.
    if f in VAN.REWRITE:
        _out = VAN.to_vancouver(_src, VAN_AMAP, VAN_MASTER)
        _left = VAN.unconverted(_out, VAN_MASTER, VAN_AMAP)
        require(not _left, f"{f}: citations left unconverted -> {_left[:5]}")
        _bd = VAN.bracket_delta(_src, _out)
        require(_bd["("] == _bd[")"] and _bd["["] == _bd["]"],
                f"{f}: Vancouver rewrite unbalanced brackets {_bd}")
        _src = _out
    open(os.path.join(PKG, "01_manuscript", f), "w", encoding="utf-8").write(clean_md(_src))

# ============================ 5. 02_cover_declarations ======================
cd = os.path.join(PKG, "02_cover_declarations")
for src, out, title in [("COVER_LETTER.md", "COVER_LETTER", "CKM Paper 4 cover letter"),
                        ("STROBE_MR_CHECKLIST.md", "STROBE_MR_CHECKLIST", "CKM Paper 4 STROBE-MR checklist")]:
    txt = rd(os.path.join(FIG, src))
    d = new_doc(); render_md(d, txt); finalize(d, os.path.join(cd, out + ".docx"), title)
    open(os.path.join(cd, out + ".md"), "w", encoding="utf-8").write(clean_md(txt))

title_md = front_filled.split("# Declarations")[0]
decl_md = "# Declarations" + front_filled.split("# Declarations")[1] if "# Declarations" in front_filled else ""
d = new_doc(); render_md(d, title_md); finalize(d, os.path.join(cd, "TITLE_PAGE.docx"), "CKM Paper 4 title page")
d = new_doc(); render_md(d, decl_md); finalize(d, os.path.join(cd, "DECLARATIONS.docx"), "CKM Paper 4 declarations")
open(os.path.join(cd, "TITLE_PAGE.md"), "w", encoding="utf-8").write(clean_md(title_md))
open(os.path.join(cd, "DECLARATIONS.md"), "w", encoding="utf-8").write(clean_md(decl_md))

# ============================ 6. 03_main_figures ============================
mf = os.path.join(PKG, "03_main_figures")
n_main = F["main_figures"]
for i in range(1, n_main + 1):
    for ext in ("png", "pdf"):
        src = os.path.join(FIG, f"Figure{i}.{ext}")
        require(os.path.exists(src), f"missing main figure file Figure{i}.{ext}")
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(mf, f"Figure{i}.{ext}"))
# Graphical abstract: "strongly recommended" by Cardiovascular Diabetology, ~920x300 px, uploaded
# under "figures". Its caption is baked into the image, as the journal requires.
for ext in ("png", "pdf"):
    _ga = os.path.join(FIG, f"GraphicalAbstract.{ext}")
    require(os.path.exists(_ga), f"missing GraphicalAbstract.{ext} — run figures/graphical_abstract.R")
    if os.path.exists(_ga):
        shutil.copy2(_ga, os.path.join(mf, f"GraphicalAbstract.{ext}"))

# ============================ 7. 04_supplementary ===========================
sup = os.path.join(PKG, "04_supplementary")
n_supp = F["supp_figures"]
si = new_doc(); si.add_heading("Supplementary Information", 1)

# ---- Supplementary Methods (added 2026-07-26, v4.2) -------------------------------------------
# The C1 length trim finished by MOVING mechanics out of the main Methods. That is only safe if the
# two shipped twins agree: SUPPLEMENTARY_INFORMATION.md is written wholesale from source, while the
# .docx used to pull ONLY "**Supplementary Figure Sn |" blocks by regex — so text appended to the
# source would have reached the .md and silently NOT the .docx. Both are rendered here from the same
# string, and gate_package asserts the section count matches across source, .md and .docx.
methtext = clean_md(rd(os.path.join(FIG, "SUPPL_METHODS.md")))
require(bool(methtext.strip()), "SUPPL_METHODS.md is empty or missing")
_n_meth = len(re.findall(r"^## Supplementary Methods ", methtext, flags=re.M))
require(_n_meth >= 15, f"only {_n_meth} Supplementary Methods sections found; expected the full set")
si.add_heading("Supplementary Methods", 2)
render_md(si, re.sub(r"^#\s+Supplementary Methods\s*$", "", methtext, flags=re.M))

si.add_heading("Supplementary Figures", 2)
# The source file's own H1 is "# Paper 4 — Supplementary Figure Legends", which the .docx replaces
# with the heading above while the .md kept it — so the twins differed by the "4" in "Paper 4".
# Strip it from the .md too; the heading is emitted once, here, for both representations.
legtext = clean_md(rd(os.path.join(FIG, "SUPPL_FIGURES.md")))
legtext = re.sub(r"^#\s+Paper \d+\s+—\s+Supplementary Figure Legends\s*$", "",
                 legtext, flags=re.M).strip() + "\n"
embedded = 0
for i in range(1, n_supp + 1):
    img = os.path.join(FIG, f"SupplFig{i}.png")
    if os.path.exists(img):
        si.add_picture(img, width=Inches(6.5))
        si.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        embedded += 1
    m = re.search(rf"\*\*Supplementary Figure S{i} \|.*?(?=\n\n|\*\*Supplementary Figure|\Z)",
                  legtext, flags=re.S)
    require(bool(m), f"no legend found for Supplementary Figure S{i}")
    if m:
        p = si.add_paragraph(); add_inline(p, " ".join(m.group(0).split()))
        for r_ in p.runs:
            r_.font.size = Pt(9)
require(embedded == n_supp, f"embedded {embedded} supplementary figures, expected {n_supp}")
# The two shipped twins of Additional file 2 must be the SAME document. This paragraph reached the
# .docx and not the .md, so a numbers-multiset guard across the representation boundary showed
# 245 -> 247 and the .md never told the reader where the tables were. Built once, emitted twice.
_tables_note = (f"Supplementary Tables S1–S{F['supp_tables']} are provided as a separate Excel "
                f"workbook (Supplementary_Tables.xlsx) containing {F['supp_table_sheets']} "
                f"worksheets; the contents sheet lists all tables.")
si.add_heading("Supplementary Tables", 2)
si.add_paragraph(_tables_note)
finalize(si, os.path.join(sup, "SUPPLEMENTARY_INFORMATION.docx"), "CKM Paper 4 Supplementary Information")
open(os.path.join(sup, "SUPPLEMENTARY_INFORMATION.md"), "w", encoding="utf-8").write(
    methtext.rstrip() + "\n\n" + legtext.rstrip()
    + "\n\n## Supplementary Tables\n\n" + _tables_note + "\n")

for i in range(1, n_supp + 1):
    for ext in ("png", "pdf"):
        src = os.path.join(FIG, f"SupplFig{i}.{ext}")
        require(os.path.exists(src), f"missing supplementary figure file SupplFig{i}.{ext}")
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(sup, f"SupplFig{i}.{ext}"))
shutil.copy2(xlsx_path, os.path.join(sup, "Supplementary_Tables.xlsx"))

# ============================ 8. 05_source_data =============================
# The set is DERIVED from what the figure scripts actually read, so a new panel cannot silently
# ship without its data. Extra files backing supplementary TABLES are listed explicitly.
def scripts_reading():
    """Panel data also comes from the Python preparation scripts, not only the .R renderers:
    Figure 4b's permutation nulls are generated by suppfig6_staging_null.py from the EAS and TPMI
    network files. Scanning only *.R silently dropped those two source files."""
    m = {}
    cand = (sorted(glob.glob(os.path.join(FIG, "fig[1-9].R"))) +
            sorted(glob.glob(os.path.join(FIG, "suppfig[1-9].R"))) +
            sorted(glob.glob(os.path.join(FIG, "*.py"))))
    for p in cand:
        b = os.path.basename(p)
        # A renderer is RETIRED if the file it writes is not on disk. suppfig6.R is one: its panel a
        # duplicated main Figure 1b and its 3-cohort panel became Figure 4b, so including it minted a
        # manifest row for a figure the package does not ship.
        #
        # This used to be the arithmetic test `suppfig<N>.R where N > supp_figures`. That rule
        # SILENTLY CHANGED MEANING on 2026-07-25 when the supplementary count rose from 5 to 7 and
        # 6 > 7 became false — the retired renderer would have re-entered the manifest. Deriving the
        # output name from the script's own save_fig() call cannot drift with a count.
        t = open(p, encoding="utf-8", errors="replace").read()
        out = re.search(r'save_fig\([^,]+,\s*"([A-Za-z0-9]+)"', t)
        if out and not os.path.exists(os.path.join(FIG, out.group(1) + ".pdf")):
            continue
        RENDER_OUT[b] = out.group(1) if out else None
        m[b] = sorted(set(re.findall(r"[\"']([A-Za-z0-9_\-]+\.csv)[\"']", t)))
    return m

RENDER_OUT = {}

READS = scripts_reading()
SEARCH_DIRS = [RES, SD_FIG, SD_F1]
TABLE_ONLY = ["edges_standardised.csv", "koges_triangulation.csv", "steiger_lor.csv",
              "network_eas_meta_fixed.csv", "network_eas_meta_random.csv"]

wanted = sorted({f for fs in READS.values() for f in fs} | set(TABLE_ONLY))
sd = os.path.join(PKG, "05_source_data")
resolved, unresolved = {}, []
for f in wanted:
    for d_ in SEARCH_DIRS:
        cand = os.path.join(d_, f)
        if os.path.exists(cand):
            shutil.copy2(cand, os.path.join(sd, f))
            resolved[f] = os.path.relpath(cand, BASE).replace("\\", "/")
            break
    else:
        unresolved.append(f)
require(not unresolved, f"source-data files cited by figure scripts but not found on disk: {unresolved}")

# ---- source-data manifest: figure -> file map, derived ----
lines = ["# Source data manifest — CKM Paper 4", "",
         "Every file in `05_source_data/` with the display item it backs. This map is generated by",
         "`scripts/build_submission_v1.py` from the files each figure script actually reads, so it",
         "cannot fall out of step with the figures.", "",
         "| Display item | Source-data file(s) | Origin in the analysis tree |", "|---|---|---|"]
NAMED = {"prep_fig1.py": "Figure 1 (panel data preparation)",
         "suppfig6_staging_null.py": "Figure 4b (staging permutation nulls)"}
def label(script):
    """Label a renderer by the file it WRITES, not by its own name. fig5.R and fig6.R render
    SupplFig6 and SupplFig7 since the 2026-07-25 main/supplementary split, so keying off the script
    name would have labelled two supplementary figures 'Figure 5' and 'Figure 6' in the manifest."""
    if script in NAMED:
        return NAMED[script]
    out = RENDER_OUT.get(script)
    if out:
        m = re.match(r"(SupplFig|Figure)(\d+)$", out)
        if m:
            return (f"Supplementary Figure S{m.group(2)}" if m.group(1) == "SupplFig"
                    else f"Figure {m.group(2)}")
    m = re.match(r"(supp)?fig(\d)\.R", script)
    return (f"Supplementary Figure S{m.group(2)}" if m.group(1) else f"Figure {m.group(2)}") if m else script
for s in sorted(READS, key=lambda x: (label(x).startswith("Supplementary"), label(x), x)):
    if not READS[s]:
        continue
    files = ", ".join(f"`{f}`" for f in READS[s])
    origins = ", ".join(sorted({os.path.dirname(resolved[f]) for f in READS[s] if f in resolved}))
    lines.append(f"| {label(s)} | {files} | {origins} |")
lines.append(f"| Supplementary Tables S1–S{F['supp_tables']} | " +
             ", ".join(f"`{f}`" for f in TABLE_ONLY) +
             " (plus the figure files above) | results |")
lines += ["",
          "## Not shipped as files",
          "Per-edge pipeline intermediates (`data/instruments/<EXPOSURE>.clumped.tsv`,",
          "`data/harmonised/<EXPOSURE>__<OUTCOME>.outcome.tsv`) are inputs to the harmonisation step",
          "described in the Methods and in the Supplementary Figure S3/S4 legends. They are",
          "regenerated deterministically by the analysis code rather than shipped, because they are",
          "per-trait GWAS extracts; the code repository named in the Declarations reproduces them.", ""]
open(os.path.join(sd, "SOURCE_DATA_MANIFEST.md"), "w", encoding="utf-8").write("\n".join(lines))

# ============================ 9. _author_todo ===============================
over_full = F["abstract_full_words"] > 350
over_cap = F["abstract_cap_words"] > 250
todo = f"""# Author to-do before submission — CKM Paper 4

Everything below needs information only the authors hold, or a decision only they can make. The
scientific content, figures, tables, references and package structure are complete and verified.

## Required
- [ ] Co-author list, affiliations, ORCIDs and CRediT contributions (TITLE_PAGE.docx; DECLARATIONS.docx).
      The corresponding author is filled in; co-authors are deliberately left as a placeholder rather
      than invented.
- [ ] Funding sources and grant numbers (DECLARATIONS.docx; STROBE-MR item 19).
- [ ] Competing-interests statement (DECLARATIONS.docx).
- [ ] Mint the GitHub repository and Zenodo DOI for the analysis code, then insert both
      (DECLARATIONS.docx; STROBE-MR item 20).
- [ ] Confirm TPMI and KoGES data-access wording (DECLARATIONS.docx; Supplementary Table S1).
      TPMI has no PubMed cohort-descriptor paper, so it is cited through its data-access route.
- [ ] Cover-letter date, and any suggested or excluded reviewers (COVER_LETTER.docx).
- [ ] Confirm the large-language-model disclosure in Methods, "Software and reproducibility"
      (METHODS.md; CKM_Paper4_manuscript.docx). Cardiovascular Diabetology requires that LLM use be
      documented in the Methods section. A draft statement is in place naming the assistance given
      and affirming author verification and responsibility; only the authors can confirm that it
      describes the assistance actually used, so amend it if it does not.

## Done, and no longer author decisions
- [x] Target journal chosen: **Cardiovascular Diabetology**. Citations are converted to numbered
      Vancouver at build time by `scripts/lib_vancouver.py` ({F['references_cited']} references, ordered by
      first appearance); the shipped .docx and .md carry the numbers, and the gate re-derives them.
- [x] Abstract is {F['abstract_full_words']} words against the journal's 350-word limit{'' if not over_full else ' — STILL OVER, FIX BEFORE SUBMITTING'},
      structured Background / Methods / Results / Conclusions as the journal requires. A shorter
      {F['abstract_cap_words']}-word variant is kept in the build source for a journal with a
      stricter cap; it is deliberately NOT shipped, so the package contains one abstract only.
- [x] Main text is {F['main_text_words']:,} words (Introduction {F['word_counts']['INTRODUCTION']:,},
      Methods {F['word_counts']['METHODS']:,}, Results {F['word_counts']['RESULTS']:,},
      Discussion {F['word_counts']['DISCUSSION']:,}, Conclusions {F['word_counts']['CONCLUSIONS']:,}).
      Cardiovascular Diabetology sets no main-text limit for Research articles, so no trim is needed.
- [x] Graphical abstract (920 x 300 px, caption baked in) and the Research Insights table are built
      and included; both are journal-specific items, not author decisions.

## Optional
- [ ] Redraw the two conceptual panels (Figure 1a staged DAG, Figure 2a common-driver-versus-cascade
      cartoon) in BioRender or Illustrator per figures/BIORENDER_SPEC.md. The data panels are already
      at reference caliber; these two are ggplot-drawn.
"""
open(os.path.join(PKG, "_author_todo", "AUTHOR_TODO.md"), "w", encoding="utf-8").write(todo)

# ---- R environment lockfile (round-2 item B-2 #22) ----
_lock = os.path.join(BASE, "manifest", "renv.lock")
require(os.path.exists(_lock), "manifest/renv.lock missing — run scripts/build_renv_lock.R")
if os.path.exists(_lock):
    shutil.copy2(_lock, os.path.join(PKG, "renv.lock"))

# ============================ 10. README ====================================
# renv.lock facts are READ from the lockfile, so the README cannot claim an R version or a package
# count the deposit does not actually record.
_lockj = json.load(open(os.path.join(PKG, "renv.lock"), encoding="utf-8"))
r_ver = _lockj["R"]["Version"]
n_pkgs = len(_lockj["Packages"])
w = F["word_counts"]
# Counts quoted in the README are DERIVED here, never typed: the previous cut claimed a "23-item"
# checklist that has 22 rows and a code deposit of "scripts 01-35" that runs to 58.
_strobe = rd(os.path.join(FIG, "STROBE_MR_CHECKLIST.md"))
n_strobe_rows = len([ln for ln in _strobe.splitlines()
                     if ln.startswith("| **") and "| STROBE-MR item |" not in ln])
_nums = sorted({int(re.match(r"^(\d+)", f).group(1))
                for f in os.listdir(os.path.join(BASE, "scripts")) if re.match(r"^\d", f)})
script_lo, script_hi = "%02d" % _nums[0], "%02d" % _nums[-1]

readme = f"""# Submission package {PKG_VERSION} — CKM Paper 4

**Title:** {TITLE}
**Corresponding author:** Bertrand Chin-Ming Tan, PhD (btan@mail.cgu.edu.tw), ORCID 0000-0002-2218-7115.
**Target journal (primary):** Cardiovascular Diabetology; thematic alternatives JAHA and
Circulation: Genomic and Precision Medicine.
**Assembled:** {date.today().isoformat()} by `independent_build/scripts/build_submission_v1.py`.
All files are COPIES — the working tree under `paper 4/independent_build/` is untouched.

## Contents
```
submission package {PKG_VERSION}/
├── 01_manuscript/
│   ├── CKM_Paper4_manuscript.docx      main manuscript, numbered Vancouver citations
│   └── *.md                            editable sections (Introduction, Methods, Results,
│                                       Discussion, Conclusions, List of abbreviations,
│                                       Research Insights, Abstract, Figure legends)
├── 02_cover_declarations/
│   ├── COVER_LETTER.docx / .md
│   ├── TITLE_PAGE.docx / .md           derived word counts and display-item counts
│   ├── DECLARATIONS.docx / .md
│   └── STROBE_MR_CHECKLIST.docx / .md  {n_strobe_rows}-row map (Skrivankova 2021, PMID 34698778)
├── 03_main_figures/                    Figure1–Figure{n_main} + GraphicalAbstract, PNG + vector PDF
├── 04_supplementary/
│   ├── SUPPLEMENTARY_INFORMATION.docx  {n_supp} supplementary figures embedded with legends
│   ├── SupplFig1–SupplFig{n_supp}             PNG + vector PDF
│   └── Supplementary_Tables.xlsx       S1–S{F['supp_tables']} across {F['supp_table_sheets']} worksheets
├── 05_source_data/                     per-figure source values + SOURCE_DATA_MANIFEST.md
├── _author_todo/AUTHOR_TODO.md         the only outstanding work
├── renv.lock                           R environment the analysis ran in (R {r_ver}, {n_pkgs} packages)
├── README_SUBMISSION.md
└── MANIFEST_checksums.txt              md5 + size for every file
```

## The manuscript in numbers
Every figure below is derived from the source files by `scripts/derive_facts.py` at build time and
asserted by `scripts/gate_package.py`; none is hand-maintained.

| | |
|---|---|
| Abstract | {F['abstract_full_words']} words full / {F['abstract_cap_words']} words cap variant |
| Introduction | {w['INTRODUCTION']:,} words |
| Methods | {w['METHODS']:,} words |
| Results | {w['RESULTS']:,} words |
| Discussion | {w['DISCUSSION']:,} words |
| Conclusions | {w['CONCLUSIONS']:,} words |
| Main text total | {F['main_text_words']:,} words |
| Main figures | {n_main} figures, {F['total_panels']} panels |
| Supplementary figures | {n_supp} |
| Supplementary tables | S1–S{F['supp_tables']} in {F['supp_table_sheets']} worksheets |
| References | {F['references_cited']} cited, {F['references_uncited']} uncited, {F['references_unmatched']} unmatched |
| Causal network | {F['network_edges_total']} directed edges tested; {F['edges_passing_bonferroni']} pass Bonferroni (0.05/{F['network_edges_total']} = {F['bonferroni_threshold']:.1e}) |

## Readiness
Present and gate-checked: the manuscript and all sections, {n_main} main figures ({F['total_panels']} panels),
{n_supp} supplementary figures, Supplementary Tables S1–S{F['supp_tables']}, {F['references_cited']} PubMed-verified
references with no uncited or unmatched entries, the STROBE-MR checklist, the cover letter, and
per-figure source data. The gates below check provenance, cross-file identity, counts against source
and claims against the current analysis; they are not a substitute for scientific review.

Outstanding work is author-supplied only and listed in `_author_todo/AUTHOR_TODO.md`: co-authors and
ORCIDs, funding, competing interests, the code repository DOI, the journal choice and its citation
style, and a short abstract trim.

## How this package is verified
`scripts/gate_package.py` runs four gate families over **every file in this package**, not a subset:

- **Gate A — no internal provenance.** No workflow IDs, verification notes, draft markers, local
  paths or unresolved placeholders in any shipped file, including inside .docx and .xlsx text.
- **Gate B — cross-file identity.** Every copied file is md5-compared against its source, and the
  title is compared across the front matter and the manuscript builder.
- **Gate C — numbers against source.** Every count quoted in this README, the title page and the
  supplementary information is recomputed from the section files, legends, workbook and reference
  builder, and must match exactly.
- **Gate D — claims against the current analysis.** Every shipped text file, including the cover
  letter, the STROBE checklist and the workbook, is scanned for superseded numbers, retracted
  wording and revision-history language, and every figure PDF is scanned in its text layer for the
  same. `scripts/scan_claims.py` runs that scan standalone.

## Reproducing the package
From `paper 4/independent_build/`, in order:
`build_abbreviations.py` → `build_research_insights.py` → `derive_facts.py` →
`build_numbered_refs.py` → `qa_manuscript.py` → `build_manuscript_docx.py` →
`build_supp_tables.py` → `build_submission_v1.py` → `gate_package.py`.
Figures are rendered by `figures/fig1.R`–`fig6.R`, `figures/suppfig1.R`–`suppfig5.R` and
`figures/graphical_abstract.R`; the R environment is recorded in `renv.lock`
(`Rscript scripts/build_renv_lock.R`).
Analysis code that produced `results/` is in `independent_build/scripts/` (numbered {script_lo}–{script_hi},
plus the build, gate and figure scripts); the deposit named in the Declarations is the whole directory.
"""
open(os.path.join(PKG, "README_SUBMISSION.md"), "w", encoding="utf-8").write(readme)

# ============================ 11. checksums =================================
rows = []
for root, _, files in os.walk(PKG):
    for f in sorted(files):
        if f == "MANIFEST_checksums.txt":
            continue
        p = os.path.join(root, f)
        rows.append((os.path.relpath(p, PKG).replace("\\", "/"), md5(p), os.path.getsize(p)))
rows.sort()
with open(os.path.join(PKG, "MANIFEST_checksums.txt"), "w", encoding="utf-8") as fh:
    fh.write(f"# CKM Paper 4 — submission package {PKG_VERSION} — {len(rows)} files\n")
    fh.write(f"# md5  size(bytes)  path   (generated {date.today().isoformat()})\n")
    for rel, h, sz in rows:
        fh.write(f"{h}  {sz:>10}  {rel}\n")

# ============================ 12. verdict ===================================
print(f"package: {PKG}")
print(f"files:   {len(rows) + 1} across {len(SUB)} folders")
print(f"figures: {n_main} main ({F['total_panels']} panels), {n_supp} supplementary")
print(f"source data: {len(resolved)} files")
if PROBLEMS:
    print(f"\nBUILD REFUSED — {len(PROBLEMS)} problem(s):")
    for p_ in PROBLEMS:
        print("  x", p_)
    sys.exit(1)
print("\nbuild assertions: all passed")
