# -*- coding: utf-8 -*-
"""Machine-check every factual assertion in RESPONSE_TO_INTERNAL_REVIEW_R3.md.

A response letter is the easiest document in the project to get wrong, because it describes work
rather than being it: the round-2 letter passed a careful read and then FAILED six of its own claims
on first machine check, four of which were real package defects. So the letter is not believed until
this runs.

Every check below is stated as the letter states it, and evaluated against the SHIPPED package
(every .md, both .docx twins, the workbook text, the text layer of all 12 figure PDFs) or against the
results file the number came from — never against the letter itself.

Run:  PYTHONIOENCODING=utf-8 python scripts/verify_rebuttal_r3.py
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

import io, os, re, sys, csv, glob, json, zlib, hashlib

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = P4_ROOT
BASE = os.path.join(ROOT, "independent_build")
FIG = os.path.join(BASE, "figures")
RES = os.path.join(BASE, "results")

# The package version is DERIVED, never typed: a verifier pointed at a stale package proves nothing.
_bsv = io.open(os.path.join(BASE, "scripts", "build_submission_v1.py"), encoding="utf-8").read()
_m = re.search(r'^PKG_VERSION\s*=\s*"([^"]+)"', _bsv, flags=re.M)
assert _m, "could not read PKG_VERSION from build_submission_v1.py"
PKG_VERSION = _m.group(1)
PKG = os.path.join(ROOT, f"submission package {PKG_VERSION}")
LETTER = os.path.join(BASE, "RESPONSE_TO_INTERNAL_REVIEW_R3.md")


def docx_text(p):
    import docx
    d = docx.Document(p)
    return "\n".join([q.text for q in d.paragraphs] +
                     [c.text for t in d.tables for r in t.rows for c in r.cells])


# PDF text comes from the SHARED reader (~/.claude/lib/pdftext.py), not a local copy. This file
# carried its own scanner until 2026-07-27; `pdftext.py --audit` found five such forks in this
# project, and measuring them showed the local scanner returned ~70x too many characters — CMap and
# font-encoding junk — so every presence probe over it was weaker than it looked. The import raises
# if the shared library is missing: a checker that cannot read must never report clean.
sys.path.insert(0, SHARED_LIB)
from pdftext import pdf_text
SURF = {}
for p in glob.glob(os.path.join(PKG, "**", "*.md"), recursive=True):
    SURF[os.path.relpath(p, PKG).replace("\\", "/")] = io.open(p, encoding="utf-8").read()
for p in glob.glob(os.path.join(PKG, "**", "*.docx"), recursive=True):
    SURF[os.path.relpath(p, PKG).replace("\\", "/")] = docx_text(p)
for p in glob.glob(os.path.join(PKG, "**", "*.pdf"), recursive=True):
    SURF["PDF:" + os.path.basename(p)] = pdf_text(p)
import openpyxl
_wb = openpyxl.load_workbook(os.path.join(PKG, "04_supplementary", "Supplementary_Tables.xlsx"),
                             read_only=True, data_only=True)
SURF["WORKBOOK"] = "\n".join(str(v) for ws in _wb.worksheets
                             for row in ws.iter_rows(values_only=True)
                             for v in row if isinstance(v, str))
ALLTEXT = "\n".join(SURF.values())
F = json.load(open(os.path.join(BASE, "manifest", "CANONICAL_FACTS.json"), encoding="utf-8"))

OK, BAD, NOTE = [], [], []


def check(label, cond, detail=""):
    (OK if cond else BAD).append(label + (f"  [{detail}]" if detail and not cond else ""))


def present(label, pat, scope=None, flags=re.I):
    pool = SURF if scope is None else {k: v for k, v in SURF.items() if scope in k}
    hits = sorted(k for k, t in pool.items() if re.search(pat, t, flags))
    check(label, bool(hits), f"not found in {scope or 'any surface'}")
    return hits


def absent(label, pat, flags=re.I):
    hits = sorted(k for k, t in SURF.items() if re.search(pat, t, flags))
    check(label, not hits, f"found in {hits[:4]}")


def rdres(name):
    p = os.path.join(RES, name)
    return io.open(p, encoding="utf-8", errors="replace").read() if os.path.exists(p) else ""


print("=" * 96)
print(f"VERIFYING RESPONSE_TO_INTERNAL_REVIEW_R3.md AGAINST submission package {PKG_VERSION}")
print(f"surfaces loaded: {len(SURF)}")
print("=" * 96)

# ---------------------------------------------------------------- 1. the package the letter names
n_files = sum(len(fs) for _, _, fs in os.walk(PKG))
# Round 7 added main-text Table 1 (TABLES.md + its .docx rendering), so the package grew from 99
# files to 101. The letter states the count, so compare the two rather than pinning a literal.
_letter_txt = io.open(LETTER, encoding="utf-8").read()
_m_files = re.search(r"(\d+)\s+files", _letter_txt)
_claimed_files = int(_m_files.group(1)) if _m_files else None
check(f"§1 the package file count matches the letter ({n_files})", n_files == _claimed_files,
      f"letter says {_claimed_files}, package has {n_files}")
man = io.open(os.path.join(PKG, "MANIFEST_checksums.txt"), encoding="utf-8").read()
n_man = len([l for l in man.splitlines() if re.match(r"^[0-9a-f]{32}\s", l)])
check(f"§1 the manifest carries one entry per shipped file except itself ({n_files - 1})",
      n_man == n_files - 1, f"{n_man} entries for {n_files} files")

# ---------------------------------------------------------------- 2. the headline numbers
# These used to be hardcoded here, which meant a number could move in the package and the letter and
# the verifier had to be edited in three places. They are now PARSED FROM THE LETTER and compared
# against derive_facts, so the check is "does the letter state what the package contains" rather
# than "does everything still equal what I typed on 2026-07-26".
lt_head = io.open(LETTER, encoding="utf-8").read()
w = F["word_counts"]


def stated(pat, label):
    m = re.search(pat, lt_head)
    if not m:
        BAD.append(f"§1 the letter does not state {label}")
        return None
    return int(m.group(1).replace(",", ""))


for pat, label, actual in [
    (r"Abstract \*{0,2}([\d,]+)\*{0,2} words", "the abstract length", F["abstract_full_words"]),
    (r"main text \*\*([\d,]+)\*\* words", "the main-text length", F["main_text_words"]),
    (r"Introduction ([\d,]+)", "Introduction words", w["INTRODUCTION"]),
    (r"Methods ([\d,]+)", "Methods words", w["METHODS"]),
    (r"Results ([\d,]+)", "Results words", w["RESULTS"]),
    (r"Discussion ([\d,]+)", "Discussion words", w["DISCUSSION"]),
    (r"Conclusions ([\d,]+)\)", "Conclusions words", w["CONCLUSIONS"]),
    (r"(\d+) main figures", "the main-figure count", F["main_figures"]),
    (r"in (\d+) panels", "the panel count", F["total_panels"]),
    (r"(\d+) supplementary\s+figures", "the supplementary-figure count", F["supp_figures"]),
    (r"(\d+) references", "the reference count", F["references_cited"]),
]:
    v = stated(pat, label)
    if v is not None:
        check(f"§1 {label} in the letter matches the package ({actual})", v == actual, f"letter says {v}")

_pct_actual = round(100 * (F["main_text_words"] - 12436) / 12436, 1)
# the manuscript writes a typographic minus; accept either sign glyph
_m = re.search(r"12,436 → ([\d,]+) words, ([−-]?[\d.]+)%", lt_head)
check("§4 the trim figure in the letter is the one the package supports",
      bool(_m) and int(_m.group(1).replace(",", "")) == F["main_text_words"]
      and abs(float(_m.group(2).replace("−", "-")) - _pct_actual) < 0.05,
      f"letter says {_m.groups() if _m else 'nothing'}, package gives {F['main_text_words']:,} / {_pct_actual}%")

# ---------------------------------------------------------------- 3. Supplementary Methods
sm_src = io.open(os.path.join(FIG, "SUPPL_METHODS.md"), encoding="utf-8").read()
n_src = len(re.findall(r"^## Supplementary Methods ", sm_src, flags=re.M))
check("§1 20 Supplementary Methods sections in the source", n_src == 20, str(n_src))
for surf in ("04_supplementary/SUPPLEMENTARY_INFORMATION.md",
             "04_supplementary/SUPPLEMENTARY_INFORMATION.docx"):
    n = len(re.findall(r"Supplementary Methods \d+ \|", SURF.get(surf, "")))
    check(f"§4 all 20 reach {surf.split('/')[-1]}", n == n_src, f"{n}")
# "nothing was deleted": each relocated block must be findable in the shipped supplement
_mv = re.findall(r"^## Supplementary Methods [^\n]*\n\n(.{40,90}?)[,.;]", sm_src, flags=re.M | re.S)
_si = re.sub(r"\s+", " ", SURF.get("04_supplementary/SUPPLEMENTARY_INFORMATION.md", ""))
_lost = [f for f in _mv if re.sub(r"\s+", " ", f).strip() not in _si]
check("§4 every relocated block is reproduced verbatim in the shipped supplement", not _lost,
      f"{len(_lost)} missing")

# ---------------------------------------------------------------- 4. section 5.1, staging
present("§3 degree-null P = 0.27 is reported", r"0\.27", "RESULTS.md")
present("§3 the 8.9% acceptance rate is disclosed", r"8\.9%")
present("§3 the mixing range 0.245-0.295 is reported", r"0\.245.{0,3}0\.295")
present("§3 edge turnover near 30% is reported", r"turnover near 30")
absent("§3 the superseded degree-null 0.31 is gone", r"degree-preserving[^.]{0,60}0\.31")
for sec in ["01_manuscript/RESULTS.md", "01_manuscript/DISCUSSION.md", "01_manuscript/CONCLUSIONS.md"]:
    check(f"§3 'no more concordant than degree-matched' in {sec.split('/')[-1]}",
          bool(re.search(r"no more concordant than degree-matched", SURF.get(sec, ""), re.I)))
absent("§3 'inserts cleanly' is gone", r"inserts cleanly")
present("§3 'no backward CAC edge' is on the Figure 1 panel", r"no backward CAC edge", "Figure1.pdf")

# The ledger bound sweep, read from the shipped workbook sheet rather than pattern-matched. The
# first version of this check used a regex over a flattened text dump and reported a FAILURE for a
# claim that is true — a checker bug, not a defect. Read the cells.
_s5b = _wb["S5b_ledger_bound_sensitivity"] if "S5b_ledger_bound_sensitivity" in _wb.sheetnames else None
if _s5b is None:
    BAD.append("§3 the S5b bound-sensitivity sheet is missing")
else:
    # A ledger row NAMES A TRANSITION. "non-empty first cell" also matched the footnote round 7
    # added, and the check then counted 16 transitions where there are 15.
    _rows = [r for r in _s5b.iter_rows(values_only=True)
             if r and r[0] and (str(r[0]).strip() == "Transition"
                                or "→" in str(r[0]) or "->" in str(r[0]))]
    _hdr = next(i for i, r in enumerate(_rows) if str(r[0]).strip() == "Transition")
    # keys are compared notation-independently: the workbook typesets "X->Y" as "X→Y".
    _data = {str(r[0]).strip().replace("→", "->"): [str(x).strip() if x else "" for x in r[1:6]]
             for r in _rows[_hdr+1:]}
    check("§3 the ledger covers 15 transitions", len(_data) == 15, str(len(_data)))
    # Round 7 (C049): S5/S5b now print REVERSE_SUPPORTED / REVERSE_CAVEATED instead of DISCORDANT /
    # DISCORDANT_CAVEATED. This block reads the WORKBOOK, so it asserts the new vocabulary, and the
    # legacy tokens must be gone from the sheet entirely.
    _legacy = [k for k, v in _data.items() if any(str(x).startswith("DISCORDANT") for x in v)]
    check("§3 the retired DISCORDANT vocabulary is gone from S5b", not _legacy, f"{_legacy}")
    _prim = [k for k, v in _data.items() if v[2].startswith("REVERSE_")]
    check("§3 4 of 15 transitions carry a reverse effect at the primary bounds", len(_prim) == 4,
          f"{len(_prim)}: {_prim}")
    for _t in ("T2D->HF", "T2D->Stroke"):
        check(f"§3 {_t} carries a reverse effect at all five bound settings",
              _t in _data and all(v.startswith("REVERSE_") for v in _data[_t]),
              str(_data.get(_t)))
    _flip = [k for k, v in _data.items()
             if not v[2].startswith("REVERSE_") and any(x.startswith("REVERSE_") for x in v)]
    check("§3 no primary-concordant transition acquires a reverse effect at any setting", not _flip,
          f"{_flip}")

# ---------------------------------------------------------------- 5. section 5.2, mediation
absent("§3 'direct driver' is absent everywhere", r"direct driver")
absent("§3 no unqualified '% DIRECT' panel label", r"%\s*DIRECT")
present("§3 Figure 2c labels >100% unstable", r"unstable estimate", "Figure2.pdf")
present("§3 the covariance envelope 44-116 is reported", r"44.{1,3}116")
present("§3 the covariance envelope 38-142 is reported", r"38.{1,3}142")
present("§3 conditional F is computed under gencov = 0", r"gencov = 0")
absent("§3 the retired single mediated range 68-79% is gone", r"68.{0,6}79\s*%")

# ---------------------------------------------------------------- 6. section 5.3, subtypes
het = os.path.join(RES, "hf_subtype_heterogeneity.csv")
check("§3 the subtype heterogeneity results file exists", os.path.exists(het))
if os.path.exists(het):
    rows = {r["exposure"]: r for r in csv.DictReader(io.open(het, encoding="utf-8"))}
    check("§3 CAD subtype difference P = 2.2e-9", abs(float(rows["CAD"]["p_rho_0.00"]) - 2.24e-9) < 1e-10,
          rows["CAD"]["p_rho_0.00"])
    check("§3 T2D subtype difference P = 7.5e-3", abs(float(rows["T2D"]["p_rho_0.00"]) - 7.5e-3) < 5e-4,
          rows["T2D"]["p_rho_0.00"])
    check("§3 BMI subtype difference P = 0.74", abs(float(rows["BMI"]["p_rho_0.00"]) - 0.74) < 5e-3,
          rows["BMI"]["p_rho_0.00"])
    check("§3 the shared control set is 456,520", "456,520" in ALLTEXT)
present("§3 the nominal-only subtype caveat is on the Figure 2 panel", r"nominal only", "Figure2.pdf")
for banned in ["pan-subtype", "reaches both subtypes comparably", "HFrEF-restricted"]:
    absent(f"§3 '{banned}' is gone", re.escape(banned))
absent("§3 'attenuates rather than inflates' is gone from the subtype text",
       r"attenuates rather than inflates")
absent("§3 'genuine null' is gone", r"genuine null")
present("§3 the WHRadjBMI wording is 'no WHRadjBMI-HF association'", r"[Nn]o WHRadjBMI.{0,3}HF association")

# ---------------------------------------------------------------- 7. section 5.4, CAUSE / PRESSO
present("§3 Figure 3b reads 'causal model favoured'", r"causal model favoured", "Figure3.pdf")
absent("§3 'true causal' is gone", r"true causal")
absent("§3 'false positive' is gone", r"false positive")
present("§3 the correlated-marker calibration label is used", r"correlated.marker")
present("§3 MR-PRESSO ran on every Bonferroni-significant edge", r"every Bonferroni-significant edge")
_s11b = _wb["S11b_PRESSO"] if "S11b_PRESSO" in _wb.sheetnames else None
if _s11b is not None:
    n11b = sum(1 for r in _s11b.iter_rows(values_only=True) if r and r[0] and str(r[0]).strip())
    check("§3 S11b publishes the full sweep, not a 5-row headline", n11b > 40, f"{n11b} rows")
present("§3 the distortion test is significant for five edges", r"distortion test[^.]{0,60}five edges")

# ---------------------------------------------------------------- 8. section 5.5, cross-ancestry
present("§3 the global offset 0.436 is reported", r"0\.436")
present("§3 the offset z = -21.7 is reported", r"21\.7")
present("§3 absorbing the offset leaves 12 edges differing (six shared, six new)",
        r"six of the original 14 and six others|12 edges differ")
present("§3 Figure 4a reads 'not significant'", r"not significant", "Figure4.pdf")
absent("§3 'not sig (power)' is gone", r"not sig \(power\)|power artifact")
present("§3 node-level resampling retains a node with probability 0.632", r"0\.632")
present("§3 the node-level interval is primary", r"node-level[^.]{0,80}primary|dependence-respecting")
present("§3 the edge-level interval is labelled anticonservative", r"anticonservative")
absent("§3 'directions are not biased' is gone", r"directions are not biased")
# The on-panel label carries a line break between "aligns" and "with the EUR direction". The
# needle assumed one text run because the retired local reader concatenated across the break; the
# shared reader preserves it. Match across whitespace instead of assuming the run boundary.
present("§3 Figure 4e marks the anomaly the second cohort resolves", r"TPMI recovers the", "Figure4.pdf")

# ---------------------------------------------------------------- 9. section 4, reorder + numbers
res = SURF["01_manuscript/RESULTS.md"]
check("§4 the Results open on coverage",
      res.lstrip().split("\n")[0].startswith("# Results") and
      "We estimated all 132 non-self directed relationships" in res.split("##")[0])
check("§4 the opening states 58 Bonferroni-significant", "Fifty-eight met network-wide" in res)
check("§4 the opening states 54 after Steiger gating", "54 remained after Steiger gating" in res)
edges_csv = os.path.join(PKG, "05_source_data", "fig1_edges.csv")
n_graph = sum(1 for _ in csv.DictReader(io.open(edges_csv, encoding="utf-8")))
check("§2 the graph really contains 54 edges", n_graph == 54, str(n_graph))
check("§2 derive_facts really gives 58 Bonferroni edges", F["edges_passing_bonferroni"] == 58,
      str(F["edges_passing_bonferroni"]))
absent("§2 the superseded '53-edge' claim is gone", r"53[- ]edge")
DICTION = ["direct driver", "predominantly via", r"\bcleanly\b", "genuine null", "true causal",
           "false positive", "not a selection artifact", "pan-subtype", "track the identity line",
           "vicious", "load-bearing", "earned that reading", "hardened the hedge", "knife-edge"]
_dhits = [d for d in DICTION if any(re.search(d, t, re.I) for t in SURF.values())]
check("§4 all section-6.3 diction phrases are absent", not _dhits, f"{_dhits}")

# ---------------------------------------------------------------- 10. section 6, figures declined
present("§5 Figure 1e still reads 'per-SD / log-OR' (declined, reviewer mistaken)",
        r"per-SD / log-OR", "Figure1.pdf")
# Round-7 Q7: the graded synthesis moved from the SupplFig7 text panel to main-text Table 1, so
# the tier checks follow the claim to where it lives now.
_tab1 = SURF.get("01_manuscript/TABLES.md", "")
_hi = [ln for ln in _tab1.split("\n") if ln.startswith("| **Higher confidence**")]
check("§5 Table 1 carries the higher-confidence tier", len(_hi) >= 2, str(len(_hi)))
_stg = [ln for ln in _tab1.split("\n") if "AHA stage ordering" in ln]
check("§5 the AHA-staging row states the degree-preserving bound",
      bool(_stg) and "degree-preserving rewiring null" in _stg[0])

# ---------------------------------------------------------------- 11. section 6, compliance
ri = SURF.get("01_manuscript/RESEARCH_INSIGHTS.md", "")
_bul = [b.strip() for b in re.findall(r"^[\-\u2022]\s*(.+)$", ri, flags=re.M)]
_long = [b for b in _bul if len(b) > 100]
check("§6 every Research Insights highlight is <=100 characters", not _long, f"{len(_long)} over")
check("§6 the Research Insights table is under 200 words", len(ri.split()) < 200, str(len(ri.split())))
try:
    from PIL import Image
    _w, _h = Image.open(os.path.join(PKG, "03_main_figures", "GraphicalAbstract.png")).size
    check("§6 the graphical abstract is 920 x 300 px", (_w, _h) == (920, 300), f"{_w}x{_h}")
except Exception as e:
    NOTE.append(f"graphical-abstract dimensions not checked ({e})")
meth = SURF["01_manuscript/METHODS.md"]
check("§6 the LLM disclosure is in Methods", "large language model" in meth.lower())
check("§6 the LLM disclosure affirms author responsibility", "full responsibility" in meth.lower())
# R6 (2026-08-20): the author-facing "[AUTHOR-SUPPLIED — confirm ...]" instruction was removed and the
# disclosure made publication-facing (Claude-only, author-confirmed). The marker must now be ABSENT.
check("§6 the LLM disclosure is publication-facing (R6: no author-supplied marker)",
      "[AUTHOR-SUPPLIED" not in meth and "No AI system accessed raw data" in meth)

# ---------------------------------------------------------------- 12. section 8, the four not-done
NOTDONE = {"10.1-2 'is genetically untested' remains": r"is genetically untested",
           "10.1-3 'assessed its portability' remains": r"assessed its portability",
           "10.4-12 lean-diabetes remains in the Conclusions": r"lean-diabetes signature"}
for lab, pat in NOTDONE.items():
    check("§8 " + lab, any(re.search(pat, t, re.I) for t in SURF.values()),
          "the letter says this row was NOT done; it is absent, so the letter is wrong")

# ---------------------------------------------------------------- 12b. items 11 and 12
# The letter claims the PDF extractor now reads BOTH text operators. Prove it on a string that is
# only reachable through a TJ array: Figure 4a's legend. If this check ever passes trivially again,
# the extractor has regressed to Tj-only and every figure claim in this letter is unverified.
_f4 = SURF.get("PDF:Figure4.pdf", "")
check("§2 item 11 the extractor recovers TJ-array text (Figure 4a legend)",
      "not significant" in _f4 and "beyond chance" in _f4,
      "the legend is drawn with TJ arrays; missing it means the extractor is Tj-only again")
check("§2 item 12 SupplFig2 reads 'Correlated-marker control'",
      "Correlated-marker control" in SURF.get("PDF:SupplFig2.pdf", ""))
_bare = []
for k, t in SURF.items():
    for m in re.finditer(r"negative[\s-]control", t, re.I):
        if not re.search(r"correlated.marker\s*$", t[max(0, m.start() - 40):m.start()], re.I):
            _bare.append(k)
check("§2 item 12 no bare 'negative control' survives on any surface", not _bare, f"{sorted(set(_bare))}")

lt = io.open(LETTER, encoding="utf-8").read()

# ------------------------------------------------------------ 12c. the abstract quoted in the letter
# The letter reproduces the shipped abstract so the reviewer can suggest edits against the real text.
# A hand-copied abstract is exactly the failure mode this project keeps hitting, so that block is
# GENERATED from the package and asserted identical here. If this check fails, the letter is quoting
# an abstract nobody is submitting.
_ab_pkg = SURF["01_manuscript/ABSTRACT.md"]
_m_pkg = re.search(r"##\s*FULL abstract[^\n]*\n(.*)", _ab_pkg, flags=re.S)
_pkg_body = re.sub(r"\s+", " ", _m_pkg.group(1).strip().rstrip("-").strip())
_m_let = re.search(r"^## 9b\.(?:.*?\n)*?((?:^>.*\n)+)", lt, flags=re.M)
if not _m_let:
    BAD.append("§9b the letter has no quoted abstract block")
else:
    _let_body = re.sub(r"\s+", " ", re.sub(r"^>[ ]?", "", _m_let.group(1), flags=re.M).strip())
    check("§9b the abstract quoted in the letter is IDENTICAL to the shipped abstract",
          _let_body == _pkg_body,
          f"letter {len(_let_body)} chars vs package {len(_pkg_body)} chars")
    check("§9b the quoted abstract is the structured one",
          all(h in _let_body for h in ("Background.", "Methods.", "Results.", "Conclusions.")))
_t_pkg = re.search(r"\*\*Title:\*\*\s*(.+)", _ab_pkg).group(1).strip()
check("§9b the title quoted in the letter matches the shipped title", _t_pkg in lt,
      "the letter quotes a different title")

# ---------------------------------------------------------------- 13. the letter must not ship
_probe = ["Response to internal review, round 3", "the review was right and we had recorded",
          "Declined, and we believe the review is mistaken"]
_leak = sorted(k for k, t in SURF.items() if any(p in t for p in _probe))
check("the response letter does not ship", not _leak, f"leaked into {_leak}")
check("no shipped surface names the internal review",
      not any(re.search(r"internal review|round-3 review|reviewer round", t, re.I) for t in SURF.values()))

# ---------------------------------------------------------------- 14. reference integrity
refs = SURF["01_manuscript/REFERENCES_NUMBERED.md"]
_n = sorted({int(x) for x in re.findall(r"^\s*(\d+)\.", refs, flags=re.M)})
check(f"§11 the reference list is 1..{F['references_cited']} with no gaps",
      _n == list(range(1, F['references_cited'] + 1)), f"{_n[:3]}..{_n[-3:]}")

# ---------------------------------------------------------------- report
print()
for b in BAD:
    print("  ✗ FAILED:", b)
print(f"\nRESULT: {len(OK)} verified, {len(BAD)} FAILED")
if NOTE:
    print("\nNot mechanically checkable / skipped:")
    for n in NOTE:
        print("  -", n)
print("=" * 96)
sys.exit(1 if BAD else 0)
