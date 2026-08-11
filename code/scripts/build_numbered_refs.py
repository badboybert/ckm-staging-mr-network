# -*- coding: utf-8 -*-
"""Order the author-date citations by first appearance across the assembled manuscript
(Introduction -> Methods -> Results -> Discussion) and emit a Vancouver numbered reference
list (REFERENCES_NUMBERED.md) + an author-date -> number map. Vancouver strings are read
from REFERENCES_MASTER.md so numbering stays in sync with the verified bibliography."""
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

import os, re, unicodedata

FIG = os.path.join(P4_BASE, "figures")
SECTIONS = ["INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md"]  # assembly order

# --- parse master bibliography table: | **Key** | Citation | PMID | DOI | Status | ---
master = {}   # "Zheng 2021" -> dict(cite, pmid, doi)
for line in open(os.path.join(FIG, "REFERENCES_MASTER.md"), encoding="utf-8"):
    m = re.match(r"\|\s*\*\*(.+?)\*\*\s*\|\s*(.+?)\s*\|\s*(.+?)\s*\|\s*(.+?)\s*\|", line)
    if m:
        key, cite, pmid, doi = (x.strip() for x in m.groups())
        pmid = pmid.replace("**", "").strip()
        master[key] = dict(cite=cite, pmid=pmid, doi=doi)

# --- citation token -> master key ---
# regex captures (lastname, year) from: "X et al. 2021", "X et al. (2021)", "X and Y 2017",
# and bare "X 2021" for the short single-author conceptual refs.
# name char class includes accented letters (e.g. Munafo with diacritic) so surnames are not truncated
def norm(s):
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))
master_norm = {norm(k): k for k in master}   # normalized "Lastname YYYY" -> canonical master key

# 3 alternatives: (1,2) "X et al. YYYY"  (3,4) "X and Y YYYY"  (5,6) bare "X YYYY".
# The bare form only counts when it resolves to a real master key, so noise ("AHA 2023") is ignored.
CITE_RE = re.compile(
    r"([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+et\s+al\.?\s*\(?(\d{4})\)?"
    r"|([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+and\s+[A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+\s+(\d{4})"
    r"|([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+\(?(\d{4})\)?"
)

def resolve(last, year):
    return master_norm.get(norm(f"{last} {year}"))

order = []      # first-appearance list of canonical master keys
seen = set()
unmatched = set()
for sec in SECTIONS:
    txt = re.sub(r"<!--.*?-->", "", open(os.path.join(FIG, sec), encoding="utf-8").read(), flags=re.S)
    for m in CITE_RE.finditer(txt):
        last = m.group(1) or m.group(3) or m.group(5)
        year = m.group(2) or m.group(4) or m.group(6)
        key = resolve(last, year)
        if key:
            if key not in seen:
                seen.add(key); order.append(key)
        elif m.group(1) or m.group(3):   # explicit "et al."/"and" citation that did NOT resolve -> flag
            unmatched.add(f"{last} {year}")
unmatched = sorted(unmatched)

# master keys never cited (should be empty for a clean list; data-source refs all cited in Methods)
uncited = sorted(k for k in master if k not in seen)

# --- emit numbered list ---
out = [f"# Paper 4 — Numbered reference list (Vancouver, by order of appearance)",
       f"_Auto-generated from REFERENCES_MASTER.md + section files. {len(order)} references cited._\n"]
amap = []
for i, key in enumerate(order, 1):
    r = master[key]
    out.append(f"{i}. {r['cite']} PMID: {r['pmid']}. https://doi.org/{r['doi']}")
    amap.append(f"{key} -> [{i}]")

out.append("\n---\n## Author-date -> number map (for in-text replacement)\n")
out += amap
if uncited:
    out.append("\n## Master refs NOT cited in the 4 body sections (check if needed elsewhere):")
    out += [f"- {k}" for k in uncited]
if unmatched:
    out.append("\n## In-text tokens with NO master entry (RESOLVE THESE):")
    out += [f"- {k}" for k in unmatched]

open(os.path.join(FIG, "REFERENCES_NUMBERED.md"), "w", encoding="utf-8").write("\n".join(out))
print(f"wrote REFERENCES_NUMBERED.md: {len(order)} cited, {len(uncited)} uncited, {len(unmatched)} unmatched")
print("ORDER:", ", ".join(f"{i}.{k}" for i, k in enumerate(order, 1)))
if uncited: print("UNCITED:", uncited)
if unmatched: print("UNMATCHED:", unmatched)
