# -*- coding: utf-8 -*-
"""Assemble the CKM Paper 4 submission package into a numbered folder tree (mirroring the
4-layer-null layout). Builds .docx for cover letter, title page, declarations and the STROBE-MR
checklist, and a Supplementary Information .docx with the 5 supplementary figures embedded; copies
the manuscript, main + supplementary figures, the supplementary-tables workbook, and source data.
Author/creator metadata = Bertrand Chin-Ming Tan (global rule 8). Run AFTER the manuscript docx is final."""
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

import os, re, shutil, sys
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

# DEPRECATED 2026-07-18. Superseded by build_submission_v1.py, which writes the frozen package to
# "paper 4/submission package v1/". This builder shipped internal provenance in the section .md,
# hand-typed title-page metrics that had drifted, and an incomplete source-data set; its output
# folder was deleted. Kept only as build history. Running it would recreate that folder and invite
# someone to upload the wrong one.
if __name__ == "__main__" and "--i-know-this-is-deprecated" not in sys.argv:
    sys.exit("build_submission_package.py is DEPRECATED — run scripts/build_submission_v1.py "
             "instead (then scripts/gate_package.py).")

BASE = P4_BASE
FIG  = os.path.join(BASE, "figures")
MAN  = os.path.join(BASE, "manuscript")
RES  = os.path.join(BASE, "results")
PKG  = os.path.join(BASE, "submission_package")

SUB = ["01_manuscript","02_cover_declarations","03_main_figures","04_supplementary","05_source_data","_author_todo"]
for s in SUB:
    os.makedirs(os.path.join(PKG, s), exist_ok=True)

def rd(p): return open(p, encoding="utf-8").read() if os.path.exists(p) else ""

INLINE = re.compile(r"\*\*(.+?)\*\*|\*(.+?)\*")
def add_inline(p, text):
    pos = 0
    for m in INLINE.finditer(text):
        if m.start() > pos: p.add_run(text[pos:m.start()])
        if m.group(1) is not None: p.add_run(m.group(1)).bold = True
        else: p.add_run(m.group(2)).italic = True
        pos = m.end()
    if pos < len(text): p.add_run(text[pos:])

def render_table(doc, rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    cells = [c for c in cells if not all(set(x) <= set("-: ") for x in c)]  # drop |---| separator
    if not cells: return
    t = doc.add_table(rows=len(cells), cols=len(cells[0])); t.style = "Table Grid"
    for ri, row in enumerate(cells):
        for ci, val in enumerate(row):
            cell = t.cell(ri, ci); cell.paragraphs[0].text = ""
            add_inline(cell.paragraphs[0], val)
            for run in cell.paragraphs[0].runs:
                run.font.size = Pt(8)
                if ri == 0: run.font.bold = True

def render_md(doc, md):
    md = re.sub(r"<!--.*?-->", "", md, flags=re.S)
    lines = md.split("\n"); i = 0
    while i < len(lines):
        s = lines[i].rstrip()
        if not s: i += 1; continue
        if re.match(r"^_.+_$", s): i += 1; continue
        if s.lstrip().startswith("|"):
            tbl = []
            while i < len(lines) and lines[i].lstrip().startswith("|"): tbl.append(lines[i]); i += 1
            render_table(doc, tbl); continue
        if s.startswith("## "): doc.add_heading(s[3:].strip(), 2); i += 1; continue
        if s.startswith("# "): doc.add_heading(s[2:].strip(), 1); i += 1; continue
        p = doc.add_paragraph(); add_inline(p, s); i += 1

def new_doc():
    d = Document(); n = d.styles["Normal"]; n.font.name = "Calibri"; n.font.size = Pt(10.5)
    n.paragraph_format.space_after = Pt(6); return d

def finalize(d, path, title):
    d.core_properties.author = "Bertrand Chin-Ming Tan"
    d.core_properties.last_modified_by = "Bertrand Chin-Ming Tan"
    d.core_properties.title = title
    d.save(path); print("  wrote", os.path.relpath(path, BASE))

# ---------- 01_manuscript ----------
for f in ["CKM_Paper4_manuscript.docx"]:
    shutil.copy2(os.path.join(MAN, f), os.path.join(PKG, "01_manuscript", f))
for f in ["INTRODUCTION.md","METHODS.md","RESULTS.md","DISCUSSION.md","ABSTRACT.md","FIGURE_LEGENDS.md"]:
    shutil.copy2(os.path.join(FIG, f), os.path.join(PKG, "01_manuscript", f))

# ---------- 02_cover_declarations ----------
cd = os.path.join(PKG, "02_cover_declarations")
d = new_doc(); render_md(d, rd(os.path.join(FIG, "COVER_LETTER.md"))); finalize(d, os.path.join(cd, "COVER_LETTER.docx"), "CKM Paper 4 cover letter")
d = new_doc(); render_md(d, rd(os.path.join(FIG, "STROBE_MR_CHECKLIST.md"))); finalize(d, os.path.join(cd, "STROBE_MR_CHECKLIST.docx"), "CKM Paper 4 STROBE-MR checklist")
# title page + declarations from FRONT_MATTER.md (split at "# Declarations")
fm = rd(os.path.join(FIG, "FRONT_MATTER.md"))
title_md = fm.split("# Declarations")[0]
decl_md = "# Declarations" + fm.split("# Declarations")[1] if "# Declarations" in fm else ""
d = new_doc(); render_md(d, title_md); finalize(d, os.path.join(cd, "TITLE_PAGE.docx"), "CKM Paper 4 title page")
d = new_doc(); render_md(d, decl_md); finalize(d, os.path.join(cd, "DECLARATIONS.docx"), "CKM Paper 4 declarations")

# ---------- 03_main_figures ----------
mf = os.path.join(PKG, "03_main_figures")
for i in range(1, 7):   # Figure1-6 (main figures)
    for ext in ("png", "pdf"):
        src = os.path.join(FIG, f"Figure{i}.{ext}")
        if os.path.exists(src): shutil.copy2(src, os.path.join(mf, f"Figure{i}.{ext}"))

# ---------- 04_supplementary: SI docx (figs embedded) + supp figs + workbook ----------
sup = os.path.join(PKG, "04_supplementary")
# Supplementary Information document = intro + each supp figure image with its legend
si = new_doc(); si.add_heading("Supplementary Information", 1)
si.add_heading("Supplementary Figures", 2)
legtext = rd(os.path.join(FIG, "SUPPL_FIGURES.md"))
# intro paragraph (the leading _italic_ meta is skipped by render_md; add a plain note)
for i in range(1, 6):   # SupplFig1-5 (SF6 retired -> promoted to Fig 4b)
    img = os.path.join(FIG, f"SupplFig{i}.png")
    if os.path.exists(img):
        si.add_picture(img, width=Inches(6.5))
        si.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    m = re.search(rf"\*\*Supplementary Figure S{i} \|.*?(?=\n\n|\*\*Supplementary Figure|\Z)", legtext, flags=re.S)
    if m:
        p = si.add_paragraph(); add_inline(p, " ".join(m.group(0).split()))
        for r in p.runs: r.font.size = Pt(9)
si.add_heading("Supplementary Tables", 2)
si.add_paragraph("Supplementary Tables S1–S12 are provided as a separate Excel workbook "
                 "(Supplementary_Tables.xlsx); the contents sheet lists all tables.")
finalize(si, os.path.join(sup, "SUPPLEMENTARY_INFORMATION.docx"), "CKM Paper 4 Supplementary Information")
# copy supp figure files + the STROBE checklist md + the workbook
for i in range(1, 6):   # SupplFig1-5 (SF6 retired -> promoted to Fig 4b)
    for ext in ("png", "pdf"):
        src = os.path.join(FIG, f"SupplFig{i}.{ext}")
        if os.path.exists(src): shutil.copy2(src, os.path.join(sup, f"SupplFig{i}.{ext}"))
shutil.copy2(os.path.join(MAN, "Supplementary_Tables.xlsx"), os.path.join(sup, "Supplementary_Tables.xlsx"))

# ---------- 05_source_data ----------
sd = os.path.join(PKG, "05_source_data")
for f in ["forward_local_edges.csv","mvmr_cascade.csv","mvmr_apob.csv","cause_headline.csv",
          "cause_negcontrol_pairs.csv","staging_ledger.csv","eur_vs_eas_comparison.csv",
          "network_eas_edges.csv","network_tpmi_full.csv","network_eas_meta_fixed.csv",
          "network_eas_meta_random.csv","koges_triangulation.csv","edges_standardised.csv",
          "steiger_lor.csv","mrpresso_headline.csv","sensitivity_noukb_edges.csv","leaveoneout_headline.csv"]:
    src = os.path.join(RES, f)
    if os.path.exists(src): shutil.copy2(src, os.path.join(sd, f))

# ---------- _author_todo ----------
todo = """# Author to-do before submission (CKM Paper 4)

- [ ] Co-author list, affiliations, ORCIDs, and CRediT contributions (TITLE_PAGE.docx; DECLARATIONS.docx).
- [ ] Funding sources + grant numbers (DECLARATIONS.docx; STROBE-MR item 19).
- [ ] Competing-interests statement (DECLARATIONS.docx).
- [ ] Mint GitHub repository + Zenodo DOI for the analysis code; insert both (DECLARATIONS.docx; STROBE-MR item 20).
- [ ] TPMI / KoGES data-access wording confirmation (DECLARATIONS.docx; Supplementary Table S1).
- [ ] Confirm the target journal, then convert author-date in-text citations to numbered Vancouver
      (REFERENCES_NUMBERED.md is ready, ordered by appearance; 45 references).
- [ ] Cover-letter date and any suggested/excluded reviewers (COVER_LETTER.docx).
- [ ] Confirm title-page word counts against the final manuscript.
"""
open(os.path.join(PKG, "_author_todo", "AUTHOR_TODO.md"), "w", encoding="utf-8").write(todo)

# ---------- manifest ----------
lines = ["# CKM Paper 4 — submission package manifest\n"]
for root, _, files in os.walk(PKG):
    for f in sorted(files):
        p = os.path.join(root, f); rel = os.path.relpath(p, PKG)
        lines.append(f"- {rel}  ({os.path.getsize(p):,} bytes)")
open(os.path.join(PKG, "MANIFEST.md"), "w", encoding="utf-8").write("\n".join(lines))
print("\nPackage assembled at", os.path.relpath(PKG, BASE))
print(f"{sum(len(fs) for _,_,fs in os.walk(PKG))} files across {len(SUB)} folders")
