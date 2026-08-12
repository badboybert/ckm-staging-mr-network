# -*- coding: utf-8 -*-
"""Assemble the CKM Paper 4 Supplementary Tables workbook from the verified results/ CSVs.
Each sheet is a faithful, human-readable rendering of one results file (renamed columns +
rounding only; no re-computation). Author metadata = Bertrand Chin-Ming Tan (global rule 8).

Run:  python scripts/build_supp_tables.py
Out:  manuscript/Supplementary_Tables.xlsx

S1 (data sources) and S12 (UKB-free sensitivity) are appended by later blocks once the
reference-bibliography and sensitivity-pass agents deliver; the DATA-DRIVEN sheets below
(S2-S11) are complete now."""
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

import os, csv, io, re, math
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

BASE = P4_BASE
RES  = os.path.join(BASE, "results")
OUT  = os.path.join(BASE, "manuscript")
os.makedirs(OUT, exist_ok=True)

# ---------- styling ----------
CAP_FILL  = PatternFill("solid", fgColor="1F3864")     # dark navy
HDR_FILL  = PatternFill("solid", fgColor="D9E1F2")     # light blue
CAP_FONT  = Font(bold=True, color="FFFFFF", size=11)
HDR_FONT  = Font(bold=True, size=10)
THIN      = Side(style="thin", color="BFBFBF")
BORDER    = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

def to_num(x):
    if x is None: return None
    s = str(x).strip()
    if s in ("", "NA", "NaN", "nan", "None"): return None
    try:
        return float(s)
    except ValueError:
        return s  # keep strings like "<0.001", "TRUE"

def read_csv(fname):
    with open(os.path.join(RES, fname), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def read_csv_at(path):
    """Read a CSV by absolute path (figure-data files live outside results/)."""
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
def _sci(v, sig=2):
    """Format a P-value as 'm × 10⁻ⁿ' for a caption, from the stored value — never hand-typed."""
    x = float(v)
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    m = round(x / (10 ** e), sig - 1)
    if m >= 10:                      # rounding can push 9.97 -> 10.0
        m, e = m / 10, e + 1
    ms = ("%g" % m)
    return (ms + " × 10" + str(e).translate(_SUP)) if e else ms

def typeset(v):
    """Typeset an ASCII edge key for display: "BMI->CAD" -> "BMI→CAD".

    Round 4 (T2-14) named the ASCII arrows on four supplementary FIGURES. The same notation was
    also in 259 cells of this workbook, which is Additional file 3 and is read as a table by
    reviewers. The machine-readable copies of these data keep the ASCII key: they ship separately
    in Additional file 5, and nothing joins on a workbook cell. One implementation, applied at the
    single point where every sheet writes a cell, so a new sheet cannot reintroduce the class.
    """
    return _re_disp.sub("→", v) if isinstance(v, str) else v

import re as _re_mod
_re_disp = _re_mod.compile(r"(?<=[A-Za-z0-9])->(?=[A-Za-z0-9])")


COL_W = 13.5          # the width every data column is set to, below
COL_W_A = 16          # column A is wider
CHARS_PER_UNIT = 1.0  # openpyxl column width is ~1 character per unit at the default font
PT_PER_LINE = 12.6    # 11 pt caption text at single spacing


def cap_height(text, ncol, pt=PT_PER_LINE, floor=30.0, usable=None):
    """Row height that actually fits a wrapped caption across a merged row.

    Round 5 (M4). Every caption row was pinned at 30 pt, which fits two lines. Captions here run to
    1,000+ characters, so the long ones were visibly CLIPPED in the shipped workbook — S14, S17, S18
    and S21 among them — and the workbook is Additional file 3, a published record. A fixed height is
    a hand-maintained constant standing in for a computed one: it was right for the captions that
    existed when it was typed and wrong for every one written since. Derive it instead.

    `usable` is the merged width in openpyxl column units; pass it when the sheet does not use the
    default widths (the Contents sheet is 8/62/46, not 16/13.5/13.5, and assuming the default
    under-counted its width by two thirds and left it the only clipped row in the book).
    """
    if usable is None:
        usable = COL_W_A + COL_W * max(0, ncol - 1)
    per_line = max(20.0, usable * CHARS_PER_UNIT)
    lines = max(1, math.ceil(len(str(text)) / per_line))
    return max(floor, lines * pt + 6.0)


def write_sheet(wb, name, caption, cols, rows, note=None):
    """cols = list of (src_key, display_name, kind) where kind in {b,se,p,int,str}.
    number formats: b=0.0000, se=0.0000, p=0.00E+00, int=0, str=general."""
    caption = typeset(caption)
    note = typeset(note) if note else note
    ws = wb.create_sheet(name)
    ncol = len(cols)
    # caption
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=ncol)
    c = ws.cell(1, 1, caption); c.fill = CAP_FILL; c.font = CAP_FONT
    c.alignment = Alignment(wrap_text=True, vertical="center")
    ws.row_dimensions[1].height = cap_height(caption, ncol)
    # header
    for j, (_, disp, _) in enumerate(cols, 1):
        h = ws.cell(2, j, disp); h.fill = HDR_FILL; h.font = HDR_FONT
        h.alignment = Alignment(wrap_text=True, vertical="center", horizontal="center")
        h.border = BORDER
    ws.row_dimensions[2].height = 28
    # data
    fmt = {"b": "0.0000", "se": "0.0000", "p": "0.00E+00", "int": "0", "str": None}
    for i, r in enumerate(rows, 3):
        for j, (key, _, kind) in enumerate(cols, 1):
            v = typeset(to_num(r.get(key)))
            cell = ws.cell(i, j, v)
            cell.border = BORDER
            cell.font = Font(size=9)
            if isinstance(v, float) and fmt[kind]:
                cell.number_format = fmt[kind]
            if kind in ("b", "se", "p", "int"):
                cell.alignment = Alignment(horizontal="center")
    # note row
    if note:
        nr = 3 + len(rows) + 1
        ws.merge_cells(start_row=nr, start_column=1, end_row=nr, end_column=ncol)
        nc = ws.cell(nr, 1, note); nc.font = Font(italic=True, size=9, color="595959")
        nc.alignment = Alignment(wrap_text=True, vertical="top")
        # The bottom note is merged across the same width and clips for exactly the same reason the
        # caption did; it was never given a height at all, so it inherited the single-line default.
        ws.row_dimensions[nr].height = cap_height(note, ncol, pt=11.4, floor=14.0)
    # widths
    for j in range(1, ncol + 1):
        ws.column_dimensions[get_column_letter(j)].width = 13.5
    ws.column_dimensions["A"].width = 16
    ws.freeze_panes = "A3"
    return ws

wb = openpyxl.Workbook()
wb.remove(wb.active)

# ---------- Contents sheet ----------
toc = wb.create_sheet("Contents")
toc.merge_cells("A1:C1")
# Read the title from the ONE place it is authored. Retyping it here is how this cell came to ship
# an older wording than the seven other places that carry the canonical title.
import re as _re
_front = open(os.path.join(BASE, "figures", "FRONT_MATTER.md"), encoding="utf-8").read()
_mt = _re.search(r"##\s*Title\s*\n\*\*(.+?)\*\*", _front, flags=_re.S)
assert _mt, "could not read the canonical title from FRONT_MATTER.md"
CANON_TITLE = " ".join(_mt.group(1).split())
tc = toc.cell(1, 1, f"{CANON_TITLE} — Supplementary Tables")
tc.fill = CAP_FILL; tc.font = CAP_FONT; tc.alignment = Alignment(wrap_text=True, vertical="center")
# The Contents title is built here rather than through write_sheet, so it needs the same derived
# height — with this sheet's own column widths (8 + 62 + 46), not the default ones.
toc.row_dimensions[1].height = cap_height(f"{CANON_TITLE} — Supplementary Tables", 3,
                                          usable=8 + 62 + 46, floor=44.0)
# Round 4 (T2-10): this column used to be headed "Source (results/)" — a directory in OUR build tree
# that no reader has. The filenames themselves are fine and useful: they name artifacts in the
# deposited analysis repository, which the Code availability statement points at. Only the location
# was a production note. Gate F15 bans any directory-qualified build path in reader-facing text.
for j, h in enumerate(["Table", "Title", "Source file (deposited analysis repository)"], 1):
    x = toc.cell(2, j, h); x.fill = HDR_FILL; x.font = HDR_FONT
TOC = [
    ("S1", "GWAS data sources (accessions, ancestry, sample size, reference)", "hardcoded from accessions"),
    ("S2", "Full European bidirectional MR network (132 directed edges)", "forward_local_edges.csv"),
    ("S3", "Multivariable MR: cascade-vs-common-driver direct effects", "mvmr_cascade.csv"),
    ("S3b","Multivariable MR robustness (orientation-based multivariable Egger, Q-heterogeneity)", "mvmr_robust.csv"),
    ("S4", "CAUSE: correlated-pleiotropy vs causation (10 headline edges)", "cause_headline.csv"),
    ("S4b","CAUSE correlated-marker negative control (confounded lipid pairs)", "cause_negcontrol_pairs.csv"),
    ("S5", "Falsifiable staging ledger, native-scale (15 AHA transitions)", "staging_ledger_native.csv"),
    ("S5b", "Staging-ledger bound sensitivity (verdicts across negligibility bounds)", "staging_ledger_bound_sensitivity.csv"),
    ("S6", "Cross-ancestry EUR-vs-EAS edge comparison", "eur_vs_eas_comparison.csv"),
    ("S7", "East Asian networks: Biobank Japan and TPMI", "network_eas_edges.csv / network_tpmi_full.csv"),
    ("S7b","East Asian BBJ+TPMI inverse-variance meta-analysis", "network_eas_meta_fixed.csv / _random.csv"),
    ("S8", "Population-vs-hospital →BMI triangulation (KoGES)", "koges_triangulation.csv"),
    ("S9", "Atherogenic-axis / ApoB multivariable MR", "mvmr_apob.csv"),
    ("S10","Cross-ancestry scale standardisation", "edges_standardised.csv"),
    ("S11","Liability-scale Steiger and MR-PRESSO detail", "steiger_lor.csv / mrpresso_all_edges.csv"),
    ("S12","UK-Biobank-free outcome sensitivity analysis", "sensitivity_noukb_edges.csv"),
    ("S13","Binary chronic-kidney-disease node (Wuttke 2019)", "network_ckd.csv"),
    ("S14","Sample-overlap structure (EUR UK Biobank; EAS designs)", "overlap_matrix.md"),
    ("S15","Cross-ancestry scale dictionary", "scale_dictionary.md"),
    ("S16","Coronary-artery-calcium node (Kavousi 2023, stage 3)", "network_cac.csv"),
    ("S17","All-aetiology heart-failure subtypes (Enzan 2025)", "network_hfsubtypes_allcause.csv"),
    ("S18","Adiposity distribution: WHRadjBMI (Pulit 2019)", "network_whradjbmi.csv"),
    ("S19","Zero-overlap cross-cohort MR (BBJ instruments → TPMI outcomes)", "network_crosscohort_bbj_tpmi.csv"),
    ("S20","Steiger direction margins for every staging-graph edge", "steiger_margin.csv"),
    ("S21","Strict allele-compatible X→SBP primary analysis with position-only sensitivity comparison", "sbp_strict_mr.csv"),
    ("S22","Degree-preserving null: mixing diagnostics across chain lengths and seeds", "degree_null_mixing.csv"),
]
for i, (t, title, src) in enumerate(TOC, 3):
    toc.cell(i, 1, t).font = Font(bold=True, size=9)
    toc.cell(i, 2, title).font = Font(size=9)
    toc.cell(i, 3, src).font = Font(italic=True, size=9, color="595959")
toc.column_dimensions["A"].width = 8
toc.column_dimensions["B"].width = 62
toc.column_dimensions["C"].width = 46
toc.freeze_panes = "A3"

# ---------- S1: data sources (hardcoded from accessions/references; N filled only where confident) ----------
S1_COLS = [("trait","Trait","str"),("role","Role","str"),("anc","Ancestry","str"),
           ("cohort","Cohort / consortium","str"),("acc","Accession","str"),
           ("n","Sample size","str"),("build","Build","str"),("ref","Reference","str")]
S1_ROWS = [
 # EUR exposures
 dict(trait="BMI", role="Exposure", anc="European", cohort="GIANT", acc="—", n="~700,000", build="hg19", ref="Yengo 2018"),
 dict(trait="SBP", role="Exposure", anc="European", cohort="ICBP + UK Biobank", acc="—", n=">1,000,000", build="hg19", ref="Evangelou 2018"),
 dict(trait="HbA1c", role="Exposure", anc="European", cohort="MAGIC", acc="—", n="146,806", build="hg19", ref="Chen 2021"),
 dict(trait="HDL, LDL, TC, TG", role="Exposure", anc="European", cohort="Global Lipids Genetics Consortium", acc="—", n="~1,650,000", build="hg19", ref="Graham 2021"),
 dict(trait="ApoB", role="Exposure", anc="European", cohort="UK Biobank", acc="GCST90025952", n="435,744", build="hg19/hg38", ref="Barton 2021"),
 # EUR outcomes
 dict(trait="Coronary artery disease", role="Outcome", anc="European", cohort="CARDIoGRAMplusC4D + UKB (Aragam)", acc="GCST90132314", n="181,522 cases / 1,165,690", build="hg19", ref="Aragam 2022"),
 dict(trait="Heart failure", role="Outcome", anc="European", cohort="HERMES (Henry)", acc="GCST90728695", n="139,533 cases", build="hg19", ref="Henry 2025"),
 dict(trait="Stroke", role="Outcome", anc="European", cohort="GIGASTROKE", acc="GCST90104539", n="—", build="hg19", ref="Mishra 2022"),
 dict(trait="Type 2 diabetes", role="Outcome", anc="European", cohort="DIAGRAM (Xue)", acc="GCST006867", n="62,892 cases / 596,424", build="hg19", ref="Xue 2018"),
 dict(trait="eGFR", role="Outcome", anc="European", cohort="CKDGen + UKB (Stanzick)", acc="—", n="~1,200,000", build="hg19", ref="Stanzick 2021"),
 dict(trait="Chronic kidney disease", role="Outcome", anc="European", cohort="CKDGen (Wuttke)", acc="GCST008065", n="41,395 cases / 439,303", build="hg19", ref="Wuttke 2019"),
 # EAS exposures / cohorts
 dict(trait="BMI", role="Exposure", anc="East Asian", cohort="Biobank Japan (Akiyama)", acc="GCST004904", n="173,430", build="hg19", ref="Akiyama 2017"),
 dict(trait="SBP, HbA1c, HDL, LDL, TC, TG", role="Exposure", anc="East Asian", cohort="Biobank Japan (Kanai)", acc="NBDC hum0014", n="~160,000", build="hg19", ref="Kanai 2018"),
 dict(trait="Type 2 diabetes", role="Outcome", anc="East Asian", cohort="AGEN (Spracklen)", acc="GCST010118", n="77,418 cases / 433,540", build="hg19", ref="Spracklen 2020"),
 dict(trait="Coronary artery disease", role="Outcome", anc="East Asian", cohort="Biobank Japan (Ishigaki)", acc="—", n="~29,000 cases", build="hg19", ref="Ishigaki 2020"),
 dict(trait="HF", role="Outcome", anc="East Asian", cohort="Biobank Japan", acc="GCST90668009", n="16,251 cases / 197,577 controls", build="hg38", ref="Enzan 2025"),
 dict(trait="Stroke", role="Outcome", anc="East Asian", cohort="GIGASTROKE EAS stratum", acc="GCST90104545", n="19,032 cases / 237,242 controls", build="hg38", ref="Mishra 2022"),
 dict(trait="CKD", role="Outcome", anc="East Asian", cohort="Biobank Japan (cross-population atlas)", acc="GCST90018602", n="2,117 cases / 174,345 controls", build="hg38", ref="Sakaue 2021"),
 dict(trait="BMI, SBP, lipids, HbA1c, T2D, CAD, HF, stroke, CKD", role="Exposure / Outcome", anc="East Asian (Taiwan)", cohort="Taiwan Precision Medicine Initiative (TPMI)", acc="controlled access", n="hospital-based", build="hg38", ref="TPMI data-access statement"),
 dict(trait="BMI, waist, DM, lipids, SBP, HbA1c", role="Exposure / Outcome", anc="East Asian (Korea)", cohort="KoGES (population)", acc="controlled access", n="~211,000", build="hg19", ref="Kim & Han 2017"),
 # UKB-free sensitivity outcomes
 dict(trait="Coronary artery disease", role="Sensitivity outcome (UKB-free)", anc="European", cohort="CARDIoGRAMplusC4D (1000G)", acc="ieu-a-7", n="60,801 cases / 184,305", build="hg19", ref="Nikpay 2015"),
 dict(trait="Heart failure", role="Sensitivity outcome (UKB-free)", anc="European", cohort="HERMES (no UKB)", acc="—", n="47,309 cases", build="hg19", ref="Shah 2020"),
 dict(trait="Stroke", role="Sensitivity outcome (UKB-free)", anc="European", cohort="MEGASTROKE (no UKB)", acc="GCST005838", n="40,585 cases / 446,696", build="hg19", ref="Malik 2018"),
 dict(trait="eGFR", role="Sensitivity outcome (UKB-free)", anc="European", cohort="CKDGen 2016", acc="GCST003372", n="—", build="hg19", ref="Pattaro 2016"),
 dict(trait="Type 2 diabetes", role="Sensitivity outcome (UKB-free)", anc="European", cohort="DIAMANTE (no UKB)", acc="—", n="74,124 cases", build="hg19", ref="Mahajan 2018"),
 dict(trait="CAD, HF, T2D", role="Sensitivity outcome (independent)", anc="European", cohort="FinnGen R12", acc="—", n="—", build="hg19", ref="Kurki 2023"),
]
ws1 = write_sheet(wb, "S1_data_sources",
    "Table S1 | GWAS data sources. Exposure and outcome summary statistics used to construct the European and East "
    "Asian CKM networks and the UK-Biobank-free sensitivity analysis. Sample sizes are shown where unambiguous; full "
    "details are in the cited reference. '—' indicates not separately tabulated here (see reference/accession).",
    S1_COLS, S1_ROWS,
    note="TPMI and KoGES are available under controlled access. Reference keys map to the reference list in the manuscript References section.")
ws1.column_dimensions["A"].width = 26; ws1.column_dimensions["D"].width = 30
ws1.column_dimensions["E"].width = 20; ws1.column_dimensions["F"].width = 20; ws1.column_dimensions["H"].width = 22

# ---------- S2: full EUR network (with derived I²) ----------
def add_i2(rows):
    for r in rows:
        try:
            Q = float(r.get("Q")); n = float(r.get("nsnp")); df = n - 1
            r["I2"] = round(max(0.0, (Q - df) / Q) * 100) if Q and Q > 0 else None
        except (TypeError, ValueError):
            r["I2"] = None
    return rows
write_sheet(wb, "S2_EUR_network",
    "Table S2 | Full European-ancestry bidirectional MR network: all 132 directed exposure–outcome edges. "
    "IVW is the primary estimate; MR-Egger, weighted median, Egger intercept, Cochran Q, the derived "
    "heterogeneity I² = (Q − df)/Q, and liability-scale Steiger direction are reported for each edge. Effect "
    "sizes: per-SD (continuous exposure) or per-log-OR (binary).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
     ("egger_b","Egger β","b"),("egger_p","Egger P","p"),
     ("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
     ("wm_b","WM β","b"),("wm_p","WM P","p"),
     ("Q","Cochran Q","b"),("Q_p","Q P","p"),("I2","I² (%)","int"),
     ("steiger_correct","Steiger correct","str"),("steiger_p","Steiger P","p")],
    add_i2(read_csv("forward_local_edges.csv")),
    note="Bonferroni threshold for the network = 0.05/132 = 3.8×10⁻⁴. Steiger direction for the network table uses the quantitative-trait approximation; the liability-scale recomputation for the edges carrying directional claims is Table S11. I², heterogeneity index; WM, weighted median.")

# ---------- S3: MVMR cascade ----------
write_sheet(wb, "S3_MVMR_cascade",
    "Table S3 | Multivariable MR adjudicating cascade versus common-driver mechanism. Each disease outcome is "
    "conditioned jointly on its candidate upstream exposures; the direct effect is the MVMR-IVW estimate and the "
    "total effect is the univariable estimate. Conditional F > 12 indicates adequate conditional instrument strength.",
    [("model","Model","str"),("outcome","Outcome","str"),("exposure","Exposure","str"),("nsnp","nSNP","int"),
     ("direct_b","Direct β","b"),("direct_se","Direct SE","se"),("direct_p","Direct P","p"),
     ("total_b","Total β","b"),("total_p","Total P","p"),("cond_F","Conditional F","b")],
    read_csv("mvmr_cascade.csv"),
    note="Key contrasts: BMI retains a direct effect on CAD and HF; no CAD-independent T2D→HF effect is detected, and the mediated component is substantial but imprecise and covariance-dependent, so no single proportion is reported (Table S3/mediation); LDL→stroke collapses.")

# ---------- S3b: MVMR robust ----------
write_sheet(wb, "S3b_MVMR_robust",
    "Table S3b | Multivariable MR robustness: an orientation-based multivariable Egger regression (variants oriented "
    "on the primary exposure, the intercept serving as the directional-pleiotropy test) and Q-minimisation (qhet) "
    "direct effects for the same models as Table S3. Direct slopes survive across estimators; the multivariable-Egger "
    "intercept is orientation-dependent for the CAD-outcome models (InSIDE bounded, not absent).",
    [("model","Model","str"),("outcome","Outcome","str"),("exposure","Exposure","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_p","IVW P","p"),
     ("egger_b","mv-Egger β","b"),("egger_p","mv-Egger P","p"),
     ("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
     ("qhet_b","qhet β","b")],
    read_csv("mvmr_robust.csv"))

# ---------- S4: CAUSE ----------
write_sheet(wb, "S4_CAUSE",
    "Table S4 | CAUSE analysis of 10 headline edges, separating correlated pleiotropy from causation. A negative "
    "Δ-ELPD z (sharing-vs-causal) with P < 0.05 favours the causal model; γ is the causal effect (posterior "
    "median, 95% credible interval); q is the proportion of variants showing correlated pleiotropy.",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("n_merged","N variants","int"),
     ("z_sharing_vs_causal","Δ-ELPD z","b"),("p","P","p"),
     ("gamma_med","γ (median)","b"),("gamma_lo","γ low","b"),("gamma_hi","γ high","b"),
     ("eta_med","η (median)","b"),("q_med","q (median)","b"),("verdict","Verdict","str")],
    read_csv("cause_headline.csv"),
    note="CAD→HF strongly causal; CAD→T2D and HF→T2D causal (reciprocal feedback); T2D→HF and HbA1c→T2D return as sharing.")

# ---------- S4b: CAUSE neg control ----------
write_sheet(wb, "S4b_CAUSE_negcontrol",
    "Table S4b | CAUSE correlated-marker negative control on confounded-not-causal lipid pairs. CAUSE misclassifies HDL→CAD as "
    "causal — a documented failure mode for correlated-marker exposures — so every lipid CAUSE verdict is read "
    "alongside multivariable MR (Table S9).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("n_merged","N variants","int"),
     ("z_sharing_vs_causal","Δ-ELPD z","b"),("p","P","p"),
     ("gamma_med","γ (median)","b"),("gamma_lo","γ low","b"),("gamma_hi","γ high","b"),
     ("q_med","q (median)","b"),("expected","Expected","str"),("verdict","CAUSE verdict","str")],
    read_csv("cause_negcontrol_pairs.csv"))

# ---------- S5: staging ledger ----------
write_sheet(wb, "S5_staging_ledger",
    "Table S5 | Falsifiable staging ledger, adjudicated on NATIVE effect-size scales. "
    "Each of 15 AHA stage transitions is tested against its reverse edge by comparing the reverse 95% confidence "
    "interval, in that edge's own units, with a negligibility bound set for the outcome's scale (log(1.05) = 0.049 "
    "log-odds for binary outcomes; 1 mmHg for SBP; 0.05 SD otherwise). Cross-direction power comparisons are avoided "
    "because forward and reverse estimands are measured on different scales; a power-based minimum detectable "
    "effect and a fixed-bound equivalence statistic are retained in the last two columns for reference only. "
    "The bounds are clinical judgements, and Table S5b reports every verdict across a sweep of them.",
    [("transition","Transition","str"),("stage_X","Stage X","int"),("stage_Y","Stage Y","int"),
     ("b_fwd","Forward β","b"),("p_fwd","Forward P","p"),
     ("rev_edge","Reverse edge","str"),("b_rev","Reverse β","b"),("p_rev","Reverse P","p"),("n_rev","Reverse nSNP","int"),
     ("rev_ci_lo","Rev. CI low","b"),("rev_ci_hi","Rev. CI high","b"),
     ("rev_units","Rev. units","str"),("rev_bound","Negligibility bound","b"),
     ("rev_steiger_correct","Rev. Steiger","str"),("rev_egger_int_p","Rev. Egger int. P","b"),
     ("evidence_label","Evidence label","str"),("verdict","Verdict","str"),("knife_edge","Borderline","str"),
     ("legacy_mde","Legacy 80% MDE","b"),("legacy_tost_p","Legacy TOST P","b")],
    read_csv("staging_ledger_native.csv"),
    note="Verdicts: 7 concordant, 4 indeterminate, 2 with a genuine reverse effect (HF→T2D, Stroke→T2D) and 2 pleiotropy-caveated (CAD→SBP index-event; CAD→T2D feedback-with-caveat, matching the independent H5 robust suite). All four reverse effects point into a lower-stage trait as feedback or an index-event effect, not reverse stage progression. BMI→Stroke is flagged borderline: its reverse CI limit (+0.0505) sits 0.0005 from the 0.05 SD bound, so that verdict would flip under trivial re-rounding.")

# ---------- S5b: staging-ledger bound sensitivity ----------
# The negligibility bounds in S5 are clinical judgements. This sheet is the answer to "what if you
# had chosen differently": every verdict across a sweep, produced by scripts/58_ledger_bound_sensitivity.py.
_bs = read_csv("staging_ledger_bound_sensitivity.csv")
write_sheet(wb, "S5b_ledger_bound_sensitivity",
    "Table S5b | Bound sensitivity of the staging ledger. Each transition's verdict under five negligibility-bound "
    "scenarios, from very tight (odds ratio 1.01 / 0.5 mmHg / 0.02 SD) to very loose (odds ratio 1.20 / 5 mmHg / "
    "0.20 SD); the PRIMARY column is the rule reported in Table S5. The concordant/indeterminate split is "
    "bound-conditional, but the falsification claim is not: under every scenario each reverse effect in the "
    "discordant family runs from a stage-4 disease into a stage-2 trait, as feedback or an index-event effect, and "
    "none reverses the staging order.",
    [("transition","Transition","str"),
     ("very tight","Very tight","str"),("tight","Tight","str"),("PRIMARY","PRIMARY","str"),
     ("loose","Loose","str"),("very loose","Very loose","str"),
     ("stability","Stability","str")],
    _bs)

# ---------- S6: cross-ancestry comparison ----------
write_sheet(wb, "S6_crossanc_edges",
    "Table S6 | Cross-ancestry comparison of every shared edge between European and East Asian (Biobank Japan) "
    "networks: effect sizes, P-values, direction concordance and joint significance.",
    [("edge","Edge","str"),("EUR_b","EUR β","b"),("EUR_p","EUR P","p"),
     ("EAS_b","EAS β","b"),("EAS_p","EAS P","p"),
     ("dir_concordant","Dir. concordant","str"),("both_sig","Both significant","str")],
    read_csv("eur_vs_eas_comparison.csv"))

# ---------- S7: EAS networks ----------
write_sheet(wb, "S7_EAS_BBJ",
    "Table S7 | East Asian MR network, Biobank Japan (BBJ). IVW primary estimate with weighted-median and MR-Egger "
    "sensitivity and Steiger direction.",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
     ("wm_b","WM β","b"),("egger_b","Egger β","b"),("steiger","Steiger correct","str")],
    read_csv("network_eas_edges.csv"))

write_sheet(wb, "S7_EAS_TPMI",
    "Table S7 (cont.) | East Asian MR network, Taiwan Precision Medicine Initiative (TPMI). Full battery.",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
     ("egger_b","Egger β","b"),("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
     ("wm_b","WM β","b"),("Q","Cochran Q","b"),("Q_p","Q P","p"),
     ("steiger_correct","Steiger correct","str"),("steiger_p","Steiger P","p")],
    read_csv("network_tpmi_full.csv"))

# ---------- S7b: EAS meta ----------
def merge_meta():
    fx = {(r["exposure"], r["outcome"]): r for r in read_csv("network_eas_meta_fixed.csv")}
    rd = {(r["exposure"], r["outcome"]): r for r in read_csv("network_eas_meta_random.csv")}
    out = []
    for k, r in fx.items():
        rr = rd.get(k, {})
        out.append({"exposure": k[0], "outcome": k[1], "nsnp": r.get("nsnp"),
                    "fx_b": r.get("ivw_b"), "fx_se": r.get("ivw_se"), "fx_p": r.get("ivw_p"),
                    "rd_b": rr.get("ivw_b"), "rd_se": rr.get("ivw_se"), "rd_p": rr.get("ivw_p"),
                    "steiger": r.get("steiger")})
    return out
write_sheet(wb, "S7b_EAS_meta",
    "Table S7b | Inverse-variance meta-analysis of the two East Asian hospital cohorts (BBJ + TPMI), fixed- and "
    "random-effects. Fixed-vs-random divergence is reported because high-heterogeneity edges are borderline at "
    "Bonferroni (e.g. T2D→BMI).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("fx_b","Fixed β","b"),("fx_se","Fixed SE","se"),("fx_p","Fixed P","p"),
     ("rd_b","Random β","b"),("rd_se","Random SE","se"),("rd_p","Random P","p"),
     ("steiger","Steiger correct","str")],
    merge_meta())

# ---------- S8: triangulation ----------
write_sheet(wb, "S8_triangulation",
    "Table S8 | Population-vs-hospital triangulation of the diabetes→adiposity (→BMI) edge. Within the "
    "population-based KoGES cohort, DM→BMI/WAIST (population) is compared with the hospital BBJ/TPMI estimates; "
    "the forward BMI→DM positive control is intact.",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
     ("wm_b","WM β","b"),("wm_p","WM P","p"),
     ("egger_b","Egger β","b"),("egger_int","Egger int.","b"),("egger_int_p","Egger int. P","p"),
     ("Q_p","Q P","p"),("meanF","Mean F","b"),("steiger","Steiger correct","str")],
    read_csv("koges_triangulation.csv"),
    note="DM→BMI IVW point estimate (−0.043) matches the hospital TPMI edge (−0.046): the edge does not vanish in a population cohort, so it is not a demonstrable ascertainment artifact; a causal effect is not established (IVW-random null, balanced pleiotropy). Robust suite: pleiotropy-robust estimators (weighted median −0.090, weighted mode −0.102, penalised weighted median, radial IVW) stay negative, while the conservative two-cohort Hartung–Knapp meta-analysis is non-significant (95% CI [−0.40, +0.26], P = 0.22; k = 2, I² = 0.76). The forward BMI→DM positive control is intact.")

# ---------- S9: ApoB MVMR ----------
write_sheet(wb, "S9_ApoB_MVMR",
    "Table S9 | Atherogenic-axis multivariable MR. ApoB retains a direct CAD effect representing the atherogenic "
    "lipoprotein axis; ApoB and LDL are too collinear to separate (conditional F collapses when co-modelled), so the "
    "claim is axis-level. The full-instrument-density models are included in the same sheet: the HDL→CAD direct "
    "effect attenuates by about half at full density but a residual survives, which is why HDL non-causality is "
    "not declared from MVMR alone.",
    [("model","Model","str"),("outcome","Outcome","str"),("exposure","Exposure","str"),("nsnp","nSNP","int"),
     ("direct_b","Direct β","b"),("direct_se","Direct SE","se"),("direct_p","Direct P","p"),
     ("total_b","Total β","b"),("total_p","Total P","p"),("cond_F","Conditional F","b")],
    read_csv("mvmr_apob.csv") + read_csv("mvmr_lipid_fulldensity.csv"))

# ---------- S10: standardisation ----------
write_sheet(wb, "S10_standardised",
    "Table S10 | Cross-ancestry scale standardisation. Effect sizes on native and per-SD scales; only SBP required "
    "rescaling (EUR per-mmHg → per-SD, SD = 19.3 mmHg), computed from the full-precision IVW estimate rather than "
    "the rounded native value shown. HbA1c and eGFR were flagged non-comparable rather than rescaled.",
    [("edge","Edge","str"),("exp_units_eur","EUR units","str"),("exp_units_eas","EAS units","str"),
     ("out_units","Outcome units","str"),("EUR_b","EUR β (native)","b"),("EAS_b","EAS β (native)","b"),
     ("EUR_b_perSD","EUR β (per-SD)","b"),("EAS_b_perSD","EAS β (per-SD)","b"),
     ("comparable","Comparable","str"),("scale_note","Note","str")],
    read_csv("edges_standardised.csv"))

# ---------- S11: Steiger + PRESSO ----------
write_sheet(wb, "S11_Steiger",
    "Table S11 | Liability-scale Steiger directionality (get_r_from_lor) for headline edges. HF→CAD is "
    "decisively wrong-direction (variance explained in CAD exceeds that in HF liability).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("manual_r2x","r² exposure","b"),("manual_r2y","r² outcome","b"),
     ("manual_correct","Correct direction","str"),("steiger_p","Steiger P","p"),
     ("exp_binary","Exp. binary","str"),("out_binary","Out. binary","str")],
    read_csv("steiger_lor.csv"))

write_sheet(wb, "S11b_PRESSO",
    "Table S11b | MR-PRESSO across EVERY Bonferroni-significant edge (tiered: full global, outlier-corrected and "
    "distortion battery for edges with at most 130 harmonised instruments; global test alone above that, since the "
    "outlier and distortion tests scale with resampling count times squared instrument count). The global test is "
    "significant for every edge, which is expected at these instrument counts and indicates heterogeneity rather "
    "than bias; the interpretable quantity is the distortion test. Only HF→CAD is materially distorted by outlier "
    "correction (raw 0.807 → corrected 0.509; 13 outliers; distortion P < 0.001).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("tier","Tier","str"),
     ("nsnp_harmonised","nSNP harmonised","int"),
     ("raw_b","Raw β","b"),("raw_p","Raw P","p"),
     ("corrected_b","Corrected β","b"),("corrected_p","Corrected P","p"),
     ("n_outliers","N outliers","int"),("global_p","Global P","str"),("distortion_p","Distortion P","str")],
    read_csv("mrpresso_all_edges.csv"),
    note="Tier column records which edges received the full battery. Five edges are distortion-significant: HF→CAD, HF→Stroke, T2D→HbA1c, T2D→HDL and T2D→SBP; HF→CAD is the reverse edge the manuscript already flags.")

# ---------- S12: UKB-free sensitivity ----------
# The staging numbers in this legend are PARSED from the analysis output, never typed: the sheet
# shipped 1.000 / 1e-3 for six days after the network was completed, because they were literals.
_noukb_txt = io.open(os.path.join(RES, "sensitivity_noukb_summary.txt"), encoding="utf-8").read()
# The staging P is now an EXACT enumerated tail printed in scientific notation. The old pattern
# required `0.dddd`, so it silently stopped matching the moment the summary switched format — a
# caption regex that cannot match is a caption that ships whatever it last hardcoded. Both groups
# are asserted below instead.
_m = re.search(r"^UKB-free MAIN.*?([01]\.\d{3})\s+(\d\.\d{2}e-\d{2})", _noukb_txt, re.M)
_f = re.search(r"^UKB-free FINNGEN.*?([01]\.\d{3})\s+(\d\.\d{2}e-\d{2})", _noukb_txt, re.M)
_h = re.search(r"^\s*-\s*(\d+) headline FORWARD edges retain", _noukb_txt, re.M)
assert _m and _f and _h, (
    "S12 caption: a source pattern did not match sensitivity_noukb_summary.txt. Refusing to build a "
    "caption from stale literals — regenerate the summary (scripts/35_sensitivity_summary.py) or fix "
    "the pattern. matched: main=%s finngen=%s headline=%s" % (bool(_m), bool(_f), bool(_h)))
assert _m and _f and _h, "could not parse the UKB-free staging summary — refusing to hand-type S12"
write_sheet(wb, "S12_UKBfree_sensitivity",
    "Table S12 | UK-Biobank-free outcome-side sensitivity analysis. Headline European edges re-estimated against "
    "UK-Biobank-free outcome GWAS (CARDIoGRAM CAD, HERMES-noUKB HF, Mahajan-noUKB T2D, MEGASTROKE, CKDGen) and the "
    "fully independent FinnGen R12. All %s headline forward edges retain sign and Bonferroni significance at "
    "0.05/132; the staging order rebuilt on these outcomes stays beyond chance (concordance %s, exact P = %s main; "
    "%s, exact P = %s FinnGen). Four Stroke-as-exposure directions have no UK-Biobank-free counterpart and retain their "
    "primary estimates. This removes outcome-side overlap only — the BMI and SBP exposure GWAS still include "
    "UK Biobank." % (_h.group(1), _m.group(1), _sci(_m.group(2)), _f.group(1), _sci(_f.group(2))),
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("provider","Outcome GWAS","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
     ("egger_b","Egger β","b"),("egger_p","Egger P","p"),
     ("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
     ("wm_b","WM β","b"),("wm_p","WM P","p"),
     ("Q","Cochran Q","b"),("Q_p","Q P","p"),
     ("steiger_correct","Steiger correct","str"),("steiger_p","Steiger P","p")],
    read_csv("sensitivity_noukb_edges.csv"),
    note="Feedback edges CAD→T2D and HF→T2D attenuate on the outcome-overlap-removed Mahajan T2D, but independent FinnGen keeps CAD→T2D strongly (P ≈ 9 × 10⁻⁹); HF→T2D is the one fragile edge. Steiger for FinnGen/MEGASTROKE/CKDGen/Mahajan uses a constant N (affects Steiger direction only, not β/P).")

# ============================================================================================
# New nodes and analyses added in the peer-review round (132-edge rewrite). Each is a faithful
# render of a committed results/ file (S13/S16-S19) or of the analysis text (S14/S15).
# ============================================================================================
# 16-column network schema shared by the full network, CKD and CAC nodes:
NET16_COLS = [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
    ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
    ("egger_b","Egger β","b"),("egger_p","Egger P","p"),
    ("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
    ("wm_b","WM β","b"),("wm_p","WM P","p"),
    ("Q","Cochran Q","b"),("Q_p","Q P","p"),("steiger_correct","Steiger correct","str"),("steiger_p","Steiger P","p")]
# 12-column schema (no wm_p / no Steiger) for the subtype, WHRadjBMI and cross-cohort tables:
NET12_COLS = [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
    ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
    ("egger_b","Egger β","b"),("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
    ("wm_b","WM β","b"),("Q","Cochran Q","b"),("Q_p","Q P","p")]

# ---------- S13: binary CKD node ----------
write_sheet(wb, "S13_CKD_node",
    "Table S13 | Binary chronic-kidney-disease (CKD) node (Wuttke 2019, GCST008065), a kidney-damage endpoint "
    "distinct from the quantitative eGFR node. Risk-factor→CKD edges fall in graded tiers: robust (BMI, SBP, HF: "
    "IVW and weighted median significant, clean Egger intercept) versus pleiotropy-suspect (HDL, TG: weighted median "
    "null, significant Egger intercept). CKD→cardiovascular edges are null but underpowered (~24 instruments).",
    NET16_COLS, read_csv("network_ckd.csv"),
    # Round 5. The fourth surface of the stale main-figure class: old main Figure 5 became
    # Supplementary Figure S6 in the 2026-07-25 migration, whose rewrite touched .md sources only.
    # This is a Python string inside a workbook note — a surface no rewrite and, until Gate F21 was
    # extended to panels and cells, no gate had ever opened.
    note="This node is distinct from the Zheng-replication SBP→CKD / BMI→CKD comparison shown for the ancestry-divergence analysis (Supplementary Figure S6a), which uses a separate CKD analysis.")

# ---------- S14: sample-overlap structure ----------
S14_COLS = [("item","GWAS / design","str"),("role","Role","str"),("anc","Ancestry","str"),
            ("overlap","UK Biobank / overlap","str"),("note","Note","str")]
S14_ROWS = [
 dict(item="BMI (Yengo 2018)", role="Exposure", anc="EUR", overlap="Yes (GIANT+UKB)", note="Exposure-side UKB overlap remains after the outcome-side sensitivity."),
 dict(item="SBP (Evangelou 2018)", role="Exposure", anc="EUR", overlap="Yes (UKB+ICBP)", note="Exposure-side UKB overlap remains."),
 dict(item="Lipids (GLGC 2021)", role="Exposure", anc="EUR", overlap="Yes (meta includes UK Biobank)", note="Exposure-side; the outcome-side UKB-free sensitivity (Table S12) is independent of the lipid exposure GWAS."),
 dict(item="HbA1c (Chen 2021, MAGIC)", role="Exposure", anc="EUR", overlap="Exposure-side", note="HbA1c enters only as an exposure; its UK Biobank status does not affect the outcome-side UKB-free sensitivity (Table S12)."),
 dict(item="T2D (Xue 2018, GCST006867)", role="Exposure/Outcome", anc="EUR", overlap="Yes (DIAGRAM+UKB)", note="Both roles include UKB."),
 dict(item="CAD (Aragam 2022)", role="Outcome", anc="EUR", overlap="Yes", note="UKB-free provider = CARDIoGRAM (Nikpay 2015) + FinnGen R12 (S12)."),
 dict(item="HF (HERMES/Henry 2025)", role="Outcome", anc="EUR", overlap="Yes", note="UKB-free provider = HERMES-noUKB (Shah 2020) + FinnGen (S12)."),
 dict(item="Stroke (GIGASTROKE)", role="Outcome", anc="EUR", overlap="Yes", note="UKB-free provider = MEGASTROKE (S12)."),
 dict(item="eGFR (Stanzick 2021)", role="Outcome", anc="EUR", overlap="Yes (~36% UKB)", note="UKB-free eGFR = Pattaro 2016 (S1)."),
 dict(item="CKD (Wuttke 2019)", role="Outcome", anc="EUR", overlap="Mostly non-UKB", note="Closest to non-overlapping of the EUR outcomes."),
 dict(item="BBJ exposure → BBJ outcome", role="Design", anc="EAS", overlap="Complete (one-sample)", note="Within-cohort MR; flagged as overlapping-sample."),
 dict(item="TPMI exposure → TPMI outcome", role="Design", anc="EAS", overlap="Complete (one-sample)", note="Within-cohort MR."),
 dict(item="KoGES DM → KoGES BMI", role="Design", anc="EAS", overlap="Complete (within-cohort)", note="Explicitly one-sample triangulation edge."),
 dict(item="BBJ exposure → TPMI outcome", role="Design", anc="EAS", overlap="None (zero-overlap two-sample)", note="Japanese instruments vs Taiwanese outcomes; no shared participants (Table S19)."),
]
write_sheet(wb, "S14_overlap",
    "Table S14 | Sample-overlap structure. Two-sample MR assumes little exposure–outcome sample overlap; overlap "
    "biases estimates toward the confounded observational association. In Europeans overlap arises mainly through "
    "shared UK Biobank participants; the outcome-side UK-Biobank-free sensitivity analysis (Table S12) removes it on "
    "the outcome side, and the exposure-side UKB overlap (BMI, SBP) is disclosed as a limitation. In East Asians the "
    "within-cohort networks are one-sample; the zero-overlap BBJ→TPMI design (Table S19) provides the clean check.",
    S14_COLS, S14_ROWS)

# ---------- S15: cross-ancestry scale dictionary ----------
S15_COLS = [("trait","Trait","str"),("eur","EUR source / scale","str"),("eas","EAS source / scale","str"),
            ("cmp","Comparable?","str")]
S15_ROWS = [
 dict(trait="BMI", eur="Yengo 2018, per-SD", eas="BBJ Akiyama 2017, rank-inverse-normal (≈per-SD)", cmp="Yes"),
 dict(trait="SBP", eur="Evangelou 2018, per-mmHg", eas="BBJ Kanai 2018, rank-inverse-normal (≈per-SD)", cmp="Only after conversion (EUR ×19.3 mmHg → per-SD)"),
 dict(trait="LDL / HDL / TC / TG", eur="GLGC 2021, per-SD", eas="BBJ Kanai 2018, rank-inverse-normal (≈per-SD)", cmp="Yes"),
 dict(trait="HbA1c", eur="Chen 2021 (MAGIC), %-HbA1c units", eas="BBJ Kanai 2018, rank-inverse-normal", cmp="No — sign only"),
 dict(trait="eGFR", eur="Stanzick 2021, log(eGFR) per-SD", eas="BBJ Kanai 2018, rank-inverse-normal", cmp="No — sign only"),
 dict(trait="T2D (as exposure)", eur="Xue 2018, log-odds ratio", eas="AGEN Spracklen 2020, log-odds ratio", cmp="Yes"),
 dict(trait="CAD / HF / Stroke / CKD (outcomes)", eur="observed-scale log-odds ratio", eas="observed-scale log-odds ratio", cmp="Yes"),
]
ws15 = write_sheet(wb, "S15_scale_dictionary",
    "Table S15 | Cross-ancestry scale dictionary. Per-trait exposure scale on each ancestry side and whether raw "
    "magnitudes are directly comparable. Magnitude or ratio comparisons are made only for the 'comparable' rows on a "
    "common scale; SBP is compared only after converting the European per-mmHg estimate to per-SD (SD = 19.3 mmHg); "
    "HbA1c and eGFR are read on sign concordance only, never magnitude.",
    S15_COLS, S15_ROWS)
for col, w in {"A": 30, "B": 34, "C": 40, "D": 34}.items(): ws15.column_dimensions[col].width = w
ws14 = wb["S14_overlap"]
for col, w in {"A": 30, "B": 16, "C": 8, "D": 26, "E": 46}.items(): ws14.column_dimensions[col].width = w

# ---------- S16: CAC stage-3 node ----------
# The three staging numbers are READ from the panel's own source file. This caption shipped
# "P = 1 × 10⁻⁴" — the pre-enumeration value — while Results and Figure 1g both printed 1.2 × 10⁻⁴.
_cac = read_csv_at(os.path.join(BASE, "figures", "data", "fig1_cac_staging.csv"))
_cac_base = next(r for r in _cac if r["network"].startswith("no stage 3"))
_cac_full = next(r for r in _cac if "CAC" in r["network"] and "thin" not in r["network"])
_cac_thin = next(r for r in _cac if "thin" in r["network"])
assert _cac_full["p"], "S16 caption: CAC staging P missing from fig1_cac_staging.csv"
write_sheet(wb, "S16_CAC_node",
    "Table S16 | Coronary-artery-calcium (CAC) node (Kavousi 2023, GCST90278456; European stratum, 26,909), the "
    "guideline-named subclinical stage-3 marker. Stage-1/2→CAC edges are robust (hundreds of instruments); the "
    "stage-3→4 bridge (CAC→CAD) rests on only 5–6 instruments and is pleiotropy-prone. Adding this all-forward tier "
    "raises the staging concordance from %s to %s (exact permutation P = %s); the informative result is that "
    "no backward CAC edge appeared among the directions tested — nine CAC-to-risk-factor directions were "
    "not estimated, so a backward edge could not have been detected." % (
        _cac_base["concordance"], _cac_full["concordance"], _sci(_cac_full["p"])),
    NET16_COLS, read_csv("network_cac.csv"),
    note="The 9p21 lead rs4977575 is dropped as an ambiguous palindrome, leaving CAC→CAD on 5 instruments; dropping this thin edge leaves concordance %s." % _cac_thin["concordance"])

# ---------- S17: non-ischaemic HF subtypes ----------
write_sheet(wb, "S17_HF_subtypes",
    "Table S17 | All-aetiology heart-failure subtypes (Enzan 2025, PMID 41184235; HFpEF GCST90654629, HFrEF "
    "GCST90654628). All-aetiology subtype outcomes are used because definitions excluding ischaemic "
    "cases condition directly on CAD, which would make coronary disease null against non-ischaemic HF "
    "by case definition. With ischaemic cases retained, CAD reaches both subtypes "
    "(HFpEF P = 6e-13; HFrEF P = 6e-42), adiposity reaches both, and T2D reaches HFrEF but not HFpEF (P = 0.55). "
    "Between-subtype differences are TESTED, not inferred from significance patterns: CAD differs strongly (HFrEF > HFpEF, "
    "P = 2.2e-09, surviving Bonferroni over 8 exposures), T2D only nominally (P = 7.5e-03), and BMI shows no detected "
    "difference (P = 0.74) — an absence of evidence, not equivalence. These outcomes are cross-ancestry (EUR+EAS), so "
    "predominantly European instruments act on a ~42%-East-Asian sample, which limits interpretation of magnitude; "
    "edges are read on direction, not magnitude, and the EAS component overlaps the Biobank Japan HF outcome used elsewhere.",
    NET12_COLS, read_csv_at(os.path.join(BASE, "results", "network_hfsubtypes_allcause.csv")),
    note="T2D->HFrEF is IVW-positive but attenuates under weighted-median (~0.046) and MR-Egger (~0.022); report the direction and range, not a firm magnitude. Full between-subtype heterogeneity tests across five assumed correlations are given in the deposited hf_subtype_heterogeneity.csv; rho = 0 is the conservative column because the shared control set induces positive correlation.")

# ---------- S18: WHRadjBMI (adiposity distribution) ----------
write_sheet(wb, "S18_WHRadjBMI",
    "Table S18 | Adiposity distribution independent of overall mass (WHRadjBMI; Pulit 2019, 694,649 European, "
    "combined-sex). Central fat distribution is consistent with a causal effect on CAD and T2D but not on heart "
    "failure: no WHRadjBMI–HF association was detected. Because the exposure is conditioned on BMI and the "
    "estimates are heterogeneous, this does not establish whether overall mass or fat distribution carries the "
    "risk. Steiger filtering was not applied to these edges.",
    NET12_COLS, read_csv("network_whradjbmi.csv"),
    note="Two caveats: WHRadjBMI is BMI-adjusted (conditioning on a heritable covariate can induce collider bias, so the null-to-inverse HF association is not read as protective); and WHRadjBMI shares UK Biobank with the CAD, HF and T2D outcomes, so direction is firmer than magnitude.")

# ---------- S19: zero-overlap cross-cohort ----------
_XC19 = read_csv("network_crosscohort_bbj_tpmi.csv")
_cc = [r for r in _XC19 if r["exposure"] == "CAD" and r["outcome"] == "CAD"]
assert len(_cc) == 1, f"expected exactly one CAD->CAD positive-control row, found {len(_cc)}"
_XC_CADCAD, _XC_CADCAD_P = float(_cc[0]["ivw_b"]), float(_cc[0]["ivw_p"])
write_sheet(wb, "S19_crosscohort",
    "Table S19 | Zero-overlap cross-cohort East Asian MR: Biobank Japan instruments against Taiwan (TPMI) outcomes, "
    "sharing no participants (a clean two-sample design addressing within-cohort overlap). The core cascade "
    "(BMI→CAD/HF, LDL/TC→CAD, HbA1c→T2D, SBP→stroke) replicates with no overlap; BMI→T2D does not (null cross-cohort "
    "and within-BBJ, significant only within TPMI), a discrepancy compatible with overlap and/or "
    "cohort differences rather than attributable to either alone.",
    NET12_COLS, _XC19,
    # Round 4 (T2-8): the CAD->CAD row can be misread as a causal self-edge. It is the design's
    # positive control - Biobank Japan coronary instruments carried into the INDEPENDENT Taiwan
    # coronary outcome - so it is labelled as instrument validation, with the estimate READ from the
    # row rather than typed, so the sentence cannot drift from the table above it.
    note="Liftover hg19→hg38; allele-compatibility checked at extraction. SBP→CAD/HF attenuate to null "
         "cross-cohort (SBP→stroke survives), consistent with reduced power of the smaller blood-pressure "
         f"instrument set. The CAD→CAD row is not a causal self-edge: it applies the Biobank Japan coronary "
         f"instruments to the independent Taiwan coronary outcome and is reported as instrument / "
         f"positive-control validation of the cross-cohort design (β = {_XC_CADCAD:+.3f}, "
         f"*P* = {_sci(_XC_CADCAD_P)}).")

# ---------- S20: Steiger direction margins (round-3 item A8) ----------
# The review noted that no complete margin table was evident in the workbook, only the default gate
# and selected liability-scale analyses. The margin is the informative quantity -- a prevalence sweep
# is near-vacuous because prevalence enters get_r_from_lor only through the allele-frequency
# conversion -- so the full per-edge table is surfaced here rather than left in results/.
STEIGER_MARGIN_COLS = [
    ("exposure", "Exposure", "str"), ("outcome", "Outcome", "str"),
    ("cross", "Cross-stage", "str"), ("direction", "Direction", "str"),
    ("r2_exp", "r² exposure", "b"), ("r2_out", "r² outcome", "b"),
    ("margin", "Margin (r²exp / r²out)", "b"),
    ("margin_lo", "Margin, low prevalence", "b"), ("margin_hi", "Margin, high prevalence", "b"),
    ("fragile", "Fragile (margin < 1.5)", "str"),
]
_marg = read_csv("steiger_margin.csv")
_mvals = sorted((float(r["margin"]), f"{r['exposure']}→{r['outcome']}") for r in _marg)
write_sheet(wb, "S20_Steiger_margins",
    "Table S20 | Steiger direction margins for every edge in the staging graph. The margin is "
    "r²exposure / r²outcome at the base prevalence; Steiger calls the direction correct when it "
    "exceeds 1, so the margin — not a prevalence sweep — is what shows whether a direction call is "
    f"fragile. The smallest margin anywhere in the graph is {_mvals[0][1]} at {_mvals[0][0]:.2f}×, and "
    f"the largest is {_mvals[-1][1]} at {_mvals[-1][0]:.0f}×. A 40-fold prevalence range moves the "
    "summed r² by about 3%, far less than any of these margins, so the direction calls are not "
    "conditional on the assumed prevalence.",
    STEIGER_MARGIN_COLS, _marg,
    note="Margin low/high are the endpoints across the prevalence grid. This analysis is one-sided: it "
         "covers edges the quantitative Steiger gate already admitted and cannot recover an edge the "
         "gate wrongly rejected.")

# ---------- S21: strict allele-compatible X->SBP promoted to primary (round-3 blocker 2) ----------
SBP_STRICT_COLS = [
    ("exposure", "Exposure", "str"), ("outcome", "Outcome", "str"),
    ("base_nsnp", "nSNP, position-only", "int"), ("base_b", "β, position-only", "b"),
    ("base_p", "P, position-only", "p"),
    ("strict_nsnp", "nSNP, strict (PRIMARY)", "int"), ("strict_b", "β, strict (PRIMARY)", "b"),
    ("strict_se", "SE, strict", "se"), ("strict_p", "P, strict (PRIMARY)", "p"),
    ("sign_same", "Sign unchanged", "str"), ("both_bonf", "Bonferroni status unchanged", "str"),
]
write_sheet(wb, "S21_SBP_strict_primary",
    "Table S21 | Disease- and trait-to-SBP edges re-estimated on strict allele-compatible instruments, "
    "which are the primary analysis. The SBP outcome statistics are keyed by position with no rsID, "
    "so extraction matched on position alone and allele compatibility was enforced only at "
    "harmonisation; harmonisation cannot repair a row selected wrongly at extraction. The strict set "
    "keeps only allele-matched and strand-flip variants, dropping every strand-ambiguous one. Point "
    "estimates differ between the two extractions — CAD→SBP is +1.48 on the strict set against +1.82 on "
    "the position-only set — and exactly one edge differs in significance status: HDL→SBP sits below the "
    "network Bonferroni threshold on the strict set, giving 58 Bonferroni-significant edges and a "
    "53-edge Steiger-directed graph. HDL→SBP is intra-stage, "
    "so cross-stage concordance is unchanged at 0.926 (25/27).",
    SBP_STRICT_COLS, read_csv("sbp_strict_mr.csv"),
    note="Position-only estimates are retained as a transparency sensitivity comparator, not as an alternative "
         "primary analysis. At the genome-wide-strict bar CAD→SBP does not clear 5e-8 (P = 1.11e-07), so that "
         "robustness variant reads 1.000 (19/19) — higher than the primary only because the "
         "discordant edge falls below a threshold.")

# ---------- S22: degree-preserving null mixing diagnostics (round-3 blocker 3) ----------
MIXING_COLS = [
    ("swaps_per_edge", "Successful swaps per edge", "int"), ("seed", "Seed", "int"),
    ("p", "Null P", "b"), ("accept_pct", "Acceptance (%)", "b"),
    ("turnover_pct", "Edge turnover (%)", "b"), ("orig_kept_pct", "Original edges kept (%)", "b"),
    ("mean_successful_swaps", "Mean successful swaps", "b"),
    ("target_swaps", "Target swaps", "int"), ("replicates_capped", "Replicates hitting the cap", "int"),
]
write_sheet(wb, "S22_degree_null_mixing",
    "Table S22 | Mixing diagnostics for the degree-preserving rewiring null, the analysis that bounds "
    "the staging conclusion. Chain length is measured in SUCCESSFUL double-edge swaps per edge — "
    "invalid proposals are retried rather than skipped, so the realised chain is the reported one — "
    "and the null P is reported at 10×, 25×, 50× and 100× across three independent seeds, alongside "
    "acceptance rate and edge turnover, so that the conclusion rests on a measured chain.",
    MIXING_COLS, read_csv("degree_null_mixing.csv"),
    note="Turnover is the fraction of the observed graph's edges no longer present after rewiring; it "
         "shows the chain leaves the neighbourhood of the observed graph rather than resampling it.")

# ---------- metadata (global rule 8) ----------
wb.properties.creator = "Bertrand Chin-Ming Tan"
wb.properties.lastModifiedBy = "Bertrand Chin-Ming Tan"
wb.properties.title = "CKM Paper 4 — Supplementary Tables"
# Derive the range from the sheets actually written; the hardcoded "S2–S11" understated a
# workbook that holds S1–S12.
_ids = sorted({int(m.group(1)) for s in wb.sheetnames
               for m in [_re.match(r"S(\d+)", s)] if m})
# The description used to paraphrase the title as a literal, and duly went stale when the title
# changed on 2026-07-25 -- still advertising the retired "cross-ancestry MR map of selected CKM
# traits" in the workbook's file properties, where no gate reads it. Reuse CANON_TITLE, which is
# already derived from FRONT_MATTER.md for the Contents sheet, rather than deriving it twice.
wb.properties.description = (
    f"Supplementary Tables S{_ids[0]}–S{_ids[-1]} ({len(wb.sheetnames)} worksheets) for: "
    f"{CANON_TITLE}")

path = os.path.join(OUT, "Supplementary_Tables.xlsx")
wb.save(path)
print("wrote", path, "with", len(wb.sheetnames), "sheets:", wb.sheetnames)
