# -*- coding: utf-8 -*-
"""Render RESPONSE_TO_INTERNAL_REVIEW_R1.md to a .docx, reusing the package builder's renderer so
formatting and author metadata match the rest of the deliverables. Author/creator = Bertrand
Chin-Ming Tan (global rule 8). Deterministic apart from the modified-date stamp.

Out: RESPONSE_TO_INTERNAL_REVIEW_R1.docx (alongside the .md)
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

import os, sys, importlib.util

BASE = P4_BASE
# Which round to render. Pass "R1".."R4"; defaults to the latest.
ROUND = sys.argv[1] if len(sys.argv) > 1 else "R4"
SRC = os.path.join(BASE, f"RESPONSE_TO_INTERNAL_REVIEW_{ROUND}.md")
OUT = os.path.join(BASE, f"RESPONSE_TO_INTERNAL_REVIEW_{ROUND}.docx")

# borrow render_md / new_doc / finalize from the package builder without running its top-level cut
spec = importlib.util.spec_from_file_location("bsv1", os.path.join(BASE, "scripts", "build_submission_v1.py"))
# build_submission_v1 runs a full package cut at import; instead, lift just the helpers we need by
# re-implementing the tiny glue here against python-docx, matching new_doc/finalize exactly.
import io
from docx import Document
from docx.shared import Pt
from datetime import datetime


def read(p):
    return io.open(p, encoding="utf-8").read()


# Minimal, self-contained markdown renderer (headings, bold/italic inline, paragraphs, hrules).
import re
def add_inline(p, text):
    # split on ** ... ** and * ... * / _ ... _
    tokens = re.split(r"(\*\*.+?\*\*|\*[^*]+?\*|_[^_]+?_|`[^`]+?`)", text)
    for tok in tokens:
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**"):
            r = p.add_run(tok[2:-2]); r.bold = True
        elif (tok.startswith("*") and tok.endswith("*")) or (tok.startswith("_") and tok.endswith("_")):
            r = p.add_run(tok[1:-1]); r.italic = True
        elif tok.startswith("`") and tok.endswith("`"):
            r = p.add_run(tok[1:-1]); r.font.name = "Consolas"
        else:
            p.add_run(tok)


def render(md, doc):
    lines = md.split("\n"); i = 0
    def is_break(s):
        return (not s.strip()) or s.startswith(("#", "- ", "* ")) or set(s.strip()) <= set("-*_ ")
    while i < len(lines):
        s = lines[i].rstrip()
        if not s.strip():
            i += 1; continue
        if s.startswith("### "):
            doc.add_heading(s[4:].strip(), 3); i += 1; continue
        if s.startswith("## "):
            doc.add_heading(s[3:].strip(), 2); i += 1; continue
        if s.startswith("# "):
            doc.add_heading(s[2:].strip(), 1); i += 1; continue
        if set(s.strip()) <= set("-") and len(s.strip()) >= 3:
            i += 1; continue
        if s.lstrip().startswith(("- ", "* ")):
            buf = [s.lstrip()[2:]]; i += 1                     # a bullet may wrap across lines
            while i < len(lines) and lines[i].strip() and not is_break(lines[i]):
                buf.append(lines[i].strip()); i += 1
            p = doc.add_paragraph(style="List Bullet"); add_inline(p, " ".join(buf)); continue
        buf = [s]; i += 1
        while i < len(lines) and not is_break(lines[i]):
            buf.append(lines[i].strip()); i += 1
        p = doc.add_paragraph(); add_inline(p, " ".join(buf))


d = Document()
n = d.styles["Normal"]; n.font.name = "Calibri"; n.font.size = Pt(10.5)
n.paragraph_format.space_after = Pt(6)
render(read(SRC), d)
cp = d.core_properties
cp.author = "Bertrand Chin-Ming Tan"
cp.last_modified_by = "Bertrand Chin-Ming Tan"
cp.title = f"CKM Paper 4 — response to internal review ({ROUND})"
cp.comments = ""
cp.created = cp.modified = datetime.now()
d.save(OUT)
print("wrote", OUT, "|", len(d.paragraphs), "paragraphs | author =", cp.author)
