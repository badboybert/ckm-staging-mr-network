# -*- coding: utf-8 -*-
"""Derive the List of abbreviations required by Cardiovascular Diabetology.

CVD: "If abbreviations are used in the text they should be defined in the text at first use, and a
list of abbreviations should be provided."

Design, so the list cannot rot:

  * GLOSSES is a *superset* dictionary of known expansions. Only entries the manuscript actually
    uses are emitted, so "listed but never used" is impossible by construction rather than by
    somebody remembering to prune.
  * Expansions are also harvested from the manuscript's own first-use "Full name (ABBR)"
    constructions, so where the text defines a term the list agrees with it automatically.
  * The check that matters, and that a hand-written list cannot do: any acronym-shaped token used
    twice or more in the body sections that has no expansion and is not explicitly classified as a
    non-abbreviation (gene symbol, software name, cross-reference, ancestry code already glossed)
    is a REPORTING DEFECT and fails the build.

Run:  python scripts/build_abbreviations.py        -> writes figures/ABBREVIATIONS.md
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

import io, os, re, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE
FIG = os.path.join(BASE, "figures")

# The body sections an editor reads. FIGURE_LEGENDS is included because its abbreviations are read
# alongside the figures; ABSTRACT because CVD asks that abbreviations be minimised there.
SECTIONS = ["ABSTRACT", "INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION",
            "CONCLUSIONS", "FIGURE_LEGENDS"]

GLOSSES = {
    # ---- traits and phenotypes in the network ----
    "CKM":       "cardiovascular–kidney–metabolic",
    "CVD":       "cardiovascular disease",
    "BMI":       "body-mass index",
    "T2D":       "type 2 diabetes",
    "DM":        "diabetes mellitus",
    "CAD":       "coronary artery disease",
    "HF":        "heart failure",
    "HFpEF":     "heart failure with preserved ejection fraction",
    "HFrEF":     "heart failure with reduced ejection fraction",
    "SBP":       "systolic blood pressure",
    "LDL":       "low-density-lipoprotein cholesterol",
    "HDL":       "high-density-lipoprotein cholesterol",
    "TG":        "triglycerides",
    "TC":        "total cholesterol",
    "HbA1c":     "glycated haemoglobin",
    "eGFR":      "estimated glomerular filtration rate",
    "CKD":       "chronic kidney disease",
    "CAC":       "coronary artery calcium",
    "ApoB":      "apolipoprotein B",
    "WHR":       "waist-to-hip ratio",
    "WHRadjBMI": "waist-to-hip ratio adjusted for body-mass index",
    # ---- methods and statistics ----
    "MR":        "Mendelian randomization",
    "MVMR":      "multivariable Mendelian randomization",
    "IVW":       "inverse-variance weighted",
    "MR-PRESSO": "MR Pleiotropy RESidual Sum and Outlier",
    "MR-Egger":  "Mendelian-randomization Egger regression",
    "rsID":      "reference SNP cluster identifier",
    "CAUSE":     "Causal Analysis Using Summary Effect estimates",
    "ELPD":      "expected log pointwise predictive density",
    "LD":        "linkage disequilibrium",
    "SNP":       "single-nucleotide polymorphism",
    "GWAS":      "genome-wide association study",
    "QC":        "quality control",
    "SD":        "standard deviation",
    "CI":        "confidence interval",
    "CrI":       "credible interval",
    "FDR":       "false discovery rate",
    "DAG":       "directed acyclic graph",
    "STROBE-MR": "Strengthening the Reporting of Observational Studies in Epidemiology — "
                 "Mendelian Randomization",
    # ---- ancestry codes ----
    "EUR":       "European ancestry",
    "EAS":       "East Asian ancestry",
    # ---- organisations, cohorts and consortia ----
    "AHA":       "American Heart Association",
    "ACC":       "American College of Cardiology",
    "ADA":       "American Diabetes Association",
    "ASN":       "American Society of Nephrology",
    "BBJ":       "Biobank Japan",
    "TPMI":      "Taiwan Precision Medicine Initiative",
    "KoGES":     "Korean Genome and Epidemiology Study",
    "UKB":       "UK Biobank",
    "GIANT":     "Genetic Investigation of Anthropometric Traits consortium",
    "GLGC":      "Global Lipids Genetics Consortium",
    "ICBP":      "International Consortium of Blood Pressure",
    "MAGIC":     "Meta-Analyses of Glucose and Insulin-related traits Consortium",
    "CKDGen":    "Chronic Kidney Disease Genetics consortium",
    "HERMES":    "Heart Failure Molecular Epidemiology for Therapeutic Targets consortium",
    "DIAMANTE":  "Diabetes Meta-Analysis of Trans-Ethnic association studies",
    "DIAGRAM":   "Diabetes Genetics Replication and Meta-analysis consortium",
    "AGEN":      "Asian Genetic Epidemiology Network",
    "GIGASTROKE": "GIGASTROKE stroke genetics consortium",
    "MEGASTROKE": "MEGASTROKE stroke genetics consortium",
    "CARDIoGRAMplusC4D": "Coronary ARtery DIsease Genome-wide Replication and Meta-analysis plus "
                         "The Coronary Artery Disease Genetics consortium",
}

# Acronym-SHAPED tokens that are not abbreviations of anything a reader needs expanded.
SOFTWARE = {"TwoSampleMR", "RadialMR", "MRPRESSO", "PLINK", "SAIGE", "REGENIE", "BOLT", "R"}
GENES = {"CDKAL1", "CDKN2A", "HHEX", "KCNQ1", "SLC30A8", "HNF1B", "MC4R", "CETP", "FTO",
         "APOE", "LPA", "PCSK9", "ITIH3", "ITIH4", "FGFR4", "LARP4B", "ARG1"}
NOT_ABBREV = {
    "UK", "Biobank",                          # parts of "UK Biobank", glossed as UKB
    "GRCh37", "GRCh38",                       # genome builds
    "INFO",                                   # imputation quality field name
    "NBDC", "FTP",                            # data-portal names appearing in Methods paths
    "FinnGen",                                # cohort proper name, not an abbreviation
    "DerSimonian", "Hartung", "Knapp", "Laird", "Steiger", "Cochran",  # surnames in method names
}
# Structural tokens the regex sees but which carry no expansion: units, roman numerals, statistics,
# supplementary cross-references (S1, S5b), and hyphenated compounds built from already-glossed parts.
IGNORE_RE = [
    re.compile(r"^S\d+[a-z]?$"),              # Supplementary Table/Figure references
    re.compile(r"^(II|III|IV|VI|VII|IX|XI)$"),
    re.compile(r"^(P|Q|F|R2|Z|N|SE|OR|HR|RR|NA|ID|DOI|PMID|URL|CSV|PDF|PNG|XLSX)$"),
    re.compile(r"-(versus|vs|to|and|adjusted|free|level|wide|strict|matched|specific|restricted)-",
               re.I),                          # hyphenated compounds, not acronyms
    re.compile(r"^Per-"), re.compile(r"^Non-"), re.compile(r"^Pre-"), re.compile(r"^Anti-"),
]

def rd(p):
    return open(p, encoding="utf-8").read() if os.path.exists(p) else ""

def body_text():
    parts = []
    for s in SECTIONS:
        t = rd(os.path.join(FIG, f"{s}.md"))
        t = re.sub(r"<!--.*?-->", "", t, flags=re.S)   # provenance comments are not manuscript text
        parts.append(t)
    return "\n\n".join(parts)

TEXT = body_text()

# ---- 1. harvest the manuscript's own first-use definitions: "Full name (ABBR)" -------------------
DEF_RE = re.compile(
    r"((?:[A-Za-z][A-Za-z0-9'\u2013-]*(?:\s+(?:of|and|the|to|for|in|per))?\s+){0,5}"
    r"[A-Za-z][A-Za-z0-9'\u2013-]*)\s+\(([A-Z][A-Za-z0-9-]{1,11})\)")
CONNECT = {"of", "and", "the", "to", "for", "in", "per"}

harvested = {}
for m in DEF_RE.finditer(TEXT):
    phrase, abbr = m.group(1).strip(), m.group(2).strip()
    if abbr in harvested:
        continue
    words = phrase.split()
    for k in range(1, len(words) + 1):
        tail = words[-k:]
        initials = "".join(w[0] for w in tail if w.lower() not in CONNECT)
        if initials.upper() == "".join(c for c in abbr if c.isalpha()).upper():
            harvested[abbr] = " ".join(tail)
            break

# ---- 2. every acronym-shaped token the text actually uses ---------------------------------------
USED_RE = re.compile(r"\b([A-Za-z]*[A-Z][A-Za-z0-9]*[A-Z0-9][A-Za-z0-9]*(?:-[A-Za-z0-9]+)*)\b")
used = {}
for m in USED_RE.finditer(TEXT):
    tok = m.group(1)
    if len(tok) < 2 or any(rx.search(tok) for rx in IGNORE_RE):
        continue
    used[tok] = used.get(tok, 0) + 1

# ---- 3. build the list from what is USED, not from what is known --------------------------------
table = {}
for tok in used:
    if tok in GLOSSES:
        table[tok] = GLOSSES[tok]
    elif tok in harvested:
        table[tok] = harvested[tok]

# ---- 4. the check a hand-written list cannot do -------------------------------------------------
KNOWN = set(table) | set(GLOSSES) | SOFTWARE | GENES | NOT_ABBREV

def derived_from_known(tok):
    """A hyphenated compound built on an already-defined term ("CAD-independent", "HFpEF-specific",
    "ApoB-over-LDL") is not itself an abbreviation needing its own entry. Entries that ARE genuine
    hyphenated abbreviations (MR-PRESSO, STROBE-MR, MR-Egger) never reach here because they are in
    GLOSSES and therefore already in `table`. Plurals of a known term are treated the same way."""
    if tok.endswith("s") and tok[:-1] in KNOWN:
        return True
    return "-" in tok and any(seg in KNOWN for seg in tok.split("-"))

problems = []
for tok, n in sorted(used.items(), key=lambda kv: (-kv[1], kv[0])):
    if n < 2 or tok in table or tok in SOFTWARE or tok in GENES or tok in NOT_ABBREV:
        continue
    if derived_from_known(tok):
        continue
    problems.append(f"used {n}x in the body sections but has no expansion anywhere: {tok}")

CORE = ["CKM", "MR", "MVMR", "BMI", "T2D", "CAD", "HF", "SBP", "LDL", "HDL", "eGFR", "CKD"]
for abbr in CORE:
    if abbr in used and abbr not in table:
        problems.append(f"core abbreviation used but not listed: {abbr}")

# ---- 5. write ------------------------------------------------------------------------------------
rows = sorted(table.items(), key=lambda kv: kv[0].lower())
out = [
    "<!-- CKM Paper 4 — List of abbreviations. REQUIRED by Cardiovascular Diabetology. GENERATED by",
    "scripts/build_abbreviations.py from the abbreviations the manuscript actually uses, with",
    "expansions taken from the manuscript's own first-use definitions where it defines them. Do not",
    "edit by hand — edit the manuscript, or the GLOSSES map in the script. The generator FAILS if an",
    "acronym is used twice or more with no expansion anywhere, and it can never list one the text",
    "does not use, because the list is built from the used set. -->",
    "",
    "# List of abbreviations",
    "",
    "; ".join(f"**{a}**, {e}" for a, e in rows) + ".",
]
open(os.path.join(FIG, "ABBREVIATIONS.md"), "w", encoding="utf-8").write("\n".join(out) + "\n")

print(f"acronym-shaped tokens seen : {len(used)}")
print(f"expansions harvested from text : {len(harvested)}")
print(f"abbreviations listed : {len(rows)}")
if problems:
    print(f"\n===== ABBREVIATION CHECK FAILED — {len(problems)} problem(s) =====")
    for p in problems:
        print("  x", p)
    sys.exit(1)
print("\n===== ABBREVIATION CHECK CLEAN =====")
