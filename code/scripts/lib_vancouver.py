# -*- coding: utf-8 -*-
"""Author-date -> numbered Vancouver, applied at BUILD TIME.

Cardiovascular Diabetology: "References must follow Vancouver style" (numbered, by order of first
appearance). The manuscript sources stay in author-date because that is the form a human can edit
without renumbering anything; the numbered form is DERIVED here, so the two cannot disagree and
switching journals is one flag, not a re-edit of every section.

`build_manuscript_docx.py`, `build_submission_v1.py` and `gate_package.py` all call the same
functions, so the shipped .md, the shipped .docx and the gate's expectation are the same transform
of the same source by construction.

Ordering is `build_numbered_refs.py`'s ordering, imported rather than re-implemented.
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

import os, re, unicodedata

FIG = os.path.join(P4_BASE, "figures")

# Assembly order for first-appearance numbering. CONCLUSIONS follows DISCUSSION in the manuscript,
# so it must be here too or a reference cited only in the Conclusions would get no number.
SECTIONS = ["INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md", "CONCLUSIONS.md"]

# Sections whose citations are rewritten. The abstract carries none (CVD forbids them there) and the
# figure legends do cite, so they are included.
REWRITE = SECTIONS + ["FIGURE_LEGENDS.md"]


def norm(s):
    return "".join(c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c))


def _master():
    m = {}
    for line in open(os.path.join(FIG, "REFERENCES_MASTER.md"), encoding="utf-8"):
        mm = re.match(r"\|\s*\*\*(.+?)\*\*\s*\|\s*(.+?)\s*\|\s*(.+?)\s*\|\s*(.+?)\s*\|", line)
        if mm:
            key, cite, pmid, doi = (x.strip() for x in mm.groups())
            m[key] = dict(cite=cite, pmid=pmid.replace("**", "").strip(), doi=doi)
    return m


# DETECTION regex. The optional parentheses let it recognise "Zheng et al. (2021)" as one token.
CITE_RE = re.compile(
    r"([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+et\s+al\.?\s*\(?(\d{4})\)?"
    r"|([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+and\s+[A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+\s+(\d{4})"
    r"|([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+\(?(\d{4})\)?")

# SUBSTITUTION regex: identical except that it NEVER consumes a parenthesis.
# Using the detection regex to substitute silently deleted the closing bracket of any citation that
# ended a parenthesis -- "(... Graham et al. 2021), and apolipoprotein B" became "... Graham et al.
# [12], and apolipoprotein B", losing the ")" and running two clauses together. The token check was
# clean throughout, because nothing author-date-shaped was left behind; only a bracket-balance check
# catches it, and that check is now part of unconverted().
SUB_RE = re.compile(
    r"([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+et\s+al\.?\s*(\d{4})"
    r"|([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+and\s+([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+(\d{4})"
    r"|([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+(\d{4})")


def _strip(p):
    return re.sub(r"<!--.*?-->", "", open(p, encoding="utf-8").read(), flags=re.S)


def numbering():
    """-> (master, ordered_keys, {master_key: number}). Order = first appearance across SECTIONS."""
    master = _master()
    mnorm = {norm(k): k for k in master}
    order, seen = [], set()
    for sec in SECTIONS:
        p = os.path.join(FIG, sec)
        if not os.path.exists(p):
            continue
        for m in CITE_RE.finditer(_strip(p)):
            last = m.group(1) or m.group(3) or m.group(5)
            year = m.group(2) or m.group(4) or m.group(6)
            key = mnorm.get(norm(f"{last} {year}"))
            if key and key not in seen:
                seen.add(key)
                order.append(key)
    return master, order, {k: i for i, k in enumerate(order, 1)}


def to_vancouver(text, amap=None, master=None):
    """Replace every resolvable author-date token with its bracketed number.

    Three in-text shapes occur and all three must survive the rewrite as grammatical English:
      "(Zheng et al. 2021)"        -> "[4]"                  parenthetical, alone
      "(GLGC 2021; Graham et al. 2021)" -> "[12,13]"         parenthetical, several
      "Zheng et al. (2021) mapped" -> "Zheng et al. [4] mapped"   narrative -- the NAME stays
    """
    if amap is None or master is None:
        master, _order, amap = numbering()
    mnorm = {norm(k): k for k in master}

    def number_of(last, year):
        key = mnorm.get(norm(f"{last} {year}"))
        return amap.get(key) if key else None

    # ---- 0. same-author multi-year: "Ndumele et al. 2023, 2026" is TWO references sharing one
    # author string. Without this rule the later years are orphaned and silently dropped -- the
    # Discussion's "(Ndumele et al. 2023, 2026)" became "(Ndumele et al. [1], 2026)", losing a
    # citation entirely. Must run before the single-token rules consume the first year.
    def multiyear(m):
        years = re.findall(r"\d{4}", m.group(2))
        nums = [number_of(m.group(1), y) for y in years]
        if any(n is None for n in nums):
            return m.group(0)
        return f"{m.group(1)} et al. [" + ",".join(str(n) for n in sorted(set(nums))) + "]"
    text = re.sub(r"([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+et\s+al\.\s*(\d{4}(?:\s*,\s*\d{4})+)", multiyear, text)

    # ---- 1. narrative form: "Name et al. (2021)" keeps the name, the year becomes the number.
    def narrative(m):
        n = number_of(m.group(1), m.group(2))
        return f"{m.group(0)[:m.group(0).index('(')]}[{n}]" if n else m.group(0)
    text = re.sub(r"([A-ZÀ-ɏ][A-Za-zÀ-ɏ\-]+)\s+et\s+al\.\s+\((\d{4})\)", narrative, text)

    # ---- 2. parenthetical groups: rewrite the WHOLE parenthesis when every element resolves.
    def parenthetical(m):
        inner = m.group(1)
        parts = re.split(r"\s*;\s*", inner)
        nums, ok = [], True
        for part in parts:
            mm = CITE_RE.fullmatch(part.strip())
            if not mm:
                ok = False
                break
            n = number_of(mm.group(1) or mm.group(3) or mm.group(5),
                          mm.group(2) or mm.group(4) or mm.group(6))
            if n is None:
                ok = False
                break
            nums.append(n)
        if ok and nums:
            return "[" + ",".join(str(n) for n in sorted(set(nums))) + "]"
        return m.group(0)
    text = re.sub(r"\(([^()]*\d{4}[^()]*)\)", parenthetical, text)

    # ---- 3. any remaining bare token inside prose or a mixed parenthesis. SUB_RE, not CITE_RE:
    # see the comment on SUB_RE for the closing-bracket bug this distinction fixes.
    def bare(m):
        if m.group(1):                                  # "Name et al. 2021"
            n = number_of(m.group(1), m.group(2))
            return f"{m.group(1)} et al. [{n}]" if n else m.group(0)
        if m.group(3):                                  # "Name and Other 2021"
            n = number_of(m.group(3), m.group(5))
            return f"{m.group(3)} and {m.group(4)} [{n}]" if n else m.group(0)
        n = number_of(m.group(6), m.group(7))           # bare "Name 2021"
        return f"{m.group(6)} [{n}]" if n else m.group(0)
    text = SUB_RE.sub(bare, text)
    return text


def unconverted(text, master=None, amap=None):
    """Everything the rewrite may have missed. Must be empty.

    Two failure modes, not one. The obvious one is a surviving author-date token. The dangerous one
    is an ORPHANED YEAR -- "Ndumele et al. [1], 2026" -- where the first citation converted and the
    second was silently dropped. That leaves no author-date token to find, so a check that only
    looks for tokens reports success on a manuscript that has lost a reference."""
    if master is None or amap is None:
        master, _order, amap = numbering()
    mnorm = {norm(k): k for k in master}
    left = []
    for m in CITE_RE.finditer(text):
        last = m.group(1) or m.group(3) or m.group(5)
        year = m.group(2) or m.group(4) or m.group(6)
        if mnorm.get(norm(f"{last} {year}")):
            left.append("token: " + " ".join(m.group(0).split()))
    for m in re.finditer(r"\[\d+(?:,\d+)*\]\s*,\s*((?:19|20)\d{2})\b", text):
        left.append("orphaned year after a citation: " +
                    " ".join(text[max(0, m.start() - 40):m.end() + 10].split()))
    return left


def bracket_delta(src, out):
    """Punctuation the rewrite must not touch. Returns the (open, close, square-open) differences;
    all three must be zero apart from the square brackets the citations themselves introduce."""
    return {"(": out.count("(") - src.count("("),
            ")": out.count(")") - src.count(")"),
            "[": out.count("[") - src.count("["),
            "]": out.count("]") - src.count("]")}


def reference_list(master, order):
    """Vancouver reference strings, numbered by first appearance, in CVD's stated example format
    (authors. Title. Journal abbrev. Year;Volume(Issue):Pages.) with the DOI appended -- CVD's own
    "Article within a journal by DOI" example carries a DOI, and the PMID aids the editor."""
    return [f"{i}. {master[k]['cite']} PMID: {master[k]['pmid']}. https://doi.org/{master[k]['doi']}"
            for i, k in enumerate(order, 1)]
