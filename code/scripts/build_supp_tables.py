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


def bh_q(pvals):
    """Benjamini-Hochberg adjusted P (q) for a list of P-values, preserving input order.

    Methods and Supplementary Methods 3 both promise "a Benjamini-Hochberg FDR computed over the
    same 132-edge table is provided as a supplement"; nothing in the workbook carried one, so the
    supplement the Methods point at did not exist. Computed here from the sheet's OWN IVW P column,
    so the family is exactly the table the reader is looking at. scipy is the reference
    implementation; the closed-form BH step-up below is the fallback and is asserted to agree with
    it when both are available."""
    idx = [i for i, p in enumerate(pvals) if isinstance(p, float)]
    ps = [pvals[i] for i in idx]
    m = len(ps)
    out = [None] * len(pvals)
    if m == 0:
        return out, 0
    order = sorted(range(m), key=lambda k: ps[k])
    q = [0.0] * m
    prev = 1.0
    for rank, k in enumerate(reversed(order), 1):        # step-up from the largest P
        val = min(prev, ps[k] * m / (m - rank + 1))
        q[k] = val
        prev = val
    try:
        from scipy.stats import false_discovery_control
        ref = list(false_discovery_control(ps, method="bh"))
        assert max(abs(a - b) for a, b in zip(q, ref)) < 1e-12, \
            "BH fallback disagrees with scipy.stats.false_discovery_control"
        q = ref
    except ImportError:
        pass
    for i, k in zip(idx, range(m)):
        out[i] = q[k]
    return out, m


def write_sheet(wb, name, caption, cols, rows, note=None, blocks=None):
    """cols = list of (src_key, display_name, kind) where kind in {b,se,p,int,str}.
    number formats: b=0.0000, se=0.0000, p=0.00E+00, int=0, str=general.

    `blocks` appends one or more clearly separated sub-tables BELOW the main data and ABOVE the
    note, each as (title, [(display_name, kind), ...], [[value, ...], ...]). Three sheets cite an
    analysis the sheet did not contain (S8's robust-estimator suite, S17's between-subtype
    heterogeneity test, S7b's East-Asian staging result); one implementation serves all three so a
    fourth cannot invent a fourth layout."""
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
    # appended sub-tables (see `blocks` above)
    r = 3 + len(rows) - 1
    for btitle, bcols, brows in (blocks or []):
        r += 2                                   # one blank spacer row separates the blocks
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=ncol)
        bt = ws.cell(r, 1, typeset(btitle)); bt.fill = CAP_FILL; bt.font = CAP_FONT
        bt.alignment = Alignment(wrap_text=True, vertical="center")
        ws.row_dimensions[r].height = cap_height(btitle, ncol, floor=18.0)
        r += 1
        for j, (disp, _) in enumerate(bcols, 1):
            h = ws.cell(r, j, disp); h.fill = HDR_FILL; h.font = HDR_FONT
            h.alignment = Alignment(wrap_text=True, vertical="center", horizontal="center")
            h.border = BORDER
        ws.row_dimensions[r].height = 28
        for brow in brows:
            r += 1
            for j, (v, (_, kind)) in enumerate(zip(brow, bcols), 1):
                v = typeset(to_num(v))
                cell = ws.cell(r, j, v); cell.border = BORDER; cell.font = Font(size=9)
                if isinstance(v, float) and fmt[kind]:
                    cell.number_format = fmt[kind]
                if kind in ("b", "se", "p", "int"):
                    cell.alignment = Alignment(horizontal="center")
    # note row
    if note:
        nr = r + 2
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
# Round 7 (C103): this list is keyed by the WORKSHEET it describes, and the assertion at the
# bottom of this script refuses to save a workbook whose data sheets and Contents rows disagree.
# The hand-maintained version shipped 26 entries for 28 data sheets — S11b had no line at all and
# the two S7 sheets shared one — because nothing compared the list to the book it indexes.
# Column C names bare artifact filenames only (Gate F15 bans directory-qualified build paths).
TOC = [
    ("S1_data_sources", "S1", "GWAS data sources (accessions, ancestry, sample size, reference)",
     "compiled from the source GWAS releases; loader provenance in 02_download.sh / lib_formats.py"),
    ("S2_EUR_network", "S2", "Full European bidirectional MR network (132 directed edges), with Benjamini-Hochberg FDR", "forward_local_edges.csv"),
    ("S3_MVMR_cascade", "S3", "Multivariable MR: cascade-vs-common-driver direct effects", "mvmr_cascade.csv"),
    ("S3b_MVMR_robust", "S3b", "Multivariable MR robustness (orientation-based multivariable Egger, Q-heterogeneity)", "mvmr_robust.csv"),
    ("S4_CAUSE", "S4", "CAUSE: correlated-pleiotropy vs causation (10 headline edges)", "cause_headline.csv"),
    ("S4b_CAUSE_negcontrol", "S4b", "CAUSE correlated-marker negative control (confounded lipid pairs)", "cause_negcontrol_pairs.csv"),
    ("S5_staging_ledger", "S5", "Falsifiable staging ledger, native-scale (15 AHA transitions)", "staging_ledger_native.csv"),
    ("S5b_ledger_bound_sensitivity", "S5b", "Staging-ledger bound sensitivity (verdicts across negligibility bounds)", "staging_ledger_bound_sensitivity.csv"),
    ("S6_crossanc_edges", "S6", "Cross-ancestry EUR-vs-EAS edge comparison", "eur_vs_eas_comparison.csv"),
    ("S7_EAS_BBJ", "S7", "East Asian MR network, Biobank Japan", "network_eas_edges.csv"),
    ("S7_EAS_TPMI", "S7 (cont.)", "East Asian MR network, Taiwan Precision Medicine Initiative", "network_tpmi_full.csv"),
    ("S7b_EAS_meta", "S7b", "East Asian BBJ+TPMI inverse-variance meta-analysis, with the enumerated staging result", "network_eas_meta_fixed.csv / _random.csv / staging_eas_meta_exact.csv"),
    ("S8_triangulation", "S8", "Population-vs-hospital →BMI triangulation (KoGES), with the robust-estimator suite", "koges_triangulation.csv / h4_leandiabetes.txt / triangulation_summary.txt"),
    ("S9_ApoB_MVMR", "S9", "Atherogenic-axis / ApoB multivariable MR", "mvmr_apob.csv"),
    ("S10_standardised", "S10", "Cross-ancestry scale standardisation", "edges_standardised.csv"),
    ("S11_Steiger", "S11", "Liability-scale Steiger directionality for headline edges", "steiger_lor.csv"),
    ("S11b_PRESSO", "S11b", "MR-PRESSO across every Bonferroni-significant edge (tiered battery)", "mrpresso_all_edges.csv"),
    ("S12_UKBfree_sensitivity", "S12", "UK-Biobank-free outcome sensitivity analysis", "sensitivity_noukb_edges.csv"),
    ("S13_CKD_node", "S13", "Binary chronic-kidney-disease node (Wuttke 2019)", "network_ckd.csv"),
    ("S14_overlap", "S14", "Sample-overlap structure (EUR UK Biobank; EAS designs)", "overlap_matrix.md"),
    ("S15_scale_dictionary", "S15", "Cross-ancestry scale dictionary", "scale_dictionary.md"),
    ("S16_CAC_node", "S16", "Coronary-artery-calcium node (Kavousi 2023, stage 3)", "network_cac.csv"),
    ("S17_HF_subtypes", "S17", "All-aetiology heart-failure subtypes (Enzan 2025), with the between-subtype heterogeneity test", "network_hfsubtypes_allcause.csv / hf_subtype_heterogeneity.csv"),
    ("S18_WHRadjBMI", "S18", "Adiposity distribution: WHRadjBMI (Pulit 2019)", "network_whradjbmi.csv"),
    ("S19_crosscohort", "S19", "Zero-overlap cross-cohort MR (BBJ instruments → TPMI outcomes)", "network_crosscohort_bbj_tpmi.csv"),
    ("S20_Steiger_margins", "S20", "Steiger direction margins for every staging-graph edge", "steiger_margin.csv"),
    ("S21_SBP_strict_primary", "S21", "Strict allele-compatible X→SBP primary analysis with position-only sensitivity comparison", "sbp_strict_mr.csv"),
    ("S22_degree_null_mixing", "S22", "Degree-preserving null: mixing diagnostics across chain lengths and seeds", "degree_null_mixing.csv"),
]
for i, (_sheet, label, title, src) in enumerate(TOC, 3):
    toc.cell(i, 1, label).font = Font(bold=True, size=9)
    toc.cell(i, 2, typeset(title)).font = Font(size=9)
    toc.cell(i, 3, src).font = Font(italic=True, size=9, color="595959")
toc.column_dimensions["A"].width = 8
toc.column_dimensions["B"].width = 62
toc.column_dimensions["C"].width = 46
toc.freeze_panes = "A3"

# ---------- S1: data sources ----------
# Round 7 (C045/C048). Every "—" accession and every UK-Biobank-INCLUSIVE case count in the
# UK-Biobank-FREE rows is resolved here, and the sample sizes for the sensitivity providers are
# DERIVED from the shipped files rather than copied from the source papers: the papers' counts
# describe the full release, the shipped files are the UKB-free subsets actually analysed.
def _noukb_n(tag):
    """(min, median, max, n_distinct) of the per-SNP N column across every harmonised outcome file
    written for one UKB-free provider. This is the sample size the analysis actually saw."""
    import glob, statistics
    vals = []
    for f in glob.glob(os.path.join(BASE, "data/harmonised_noukb", f"*_{tag}.outcome.tsv")):
        with io.open(f, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if "N" not in hdr:
                continue
            i = hdr.index("N")
            for ln in fh:
                try:
                    vals.append(float(ln.rstrip("\n").split("\t")[i]))
                except (ValueError, IndexError):
                    pass
    assert vals, f"S1: no per-SNP N recovered for the {tag} UKB-free provider"
    return min(vals), statistics.median(vals), max(vals), len(set(vals))


def _vcf_sample_meta(fname):
    """TotalCases / TotalControls from a GWAS-VCF ##SAMPLE header line."""
    import gzip
    path = os.path.join(BASE, "data/uploads/eur_noukb", fname)
    with gzip.open(path, "rt", errors="replace") as f:
        for ln in f:
            if ln.startswith("##SAMPLE="):
                d = dict(kv.split("=", 1) for kv in ln.strip().rstrip(">").split("<", 1)[1].split(",")
                         if "=" in kv)
                return d
            if not ln.startswith("##"):
                break
    raise AssertionError(f"S1: no ##SAMPLE header in {fname}")


_HERMES_N = _noukb_n("HF_HERMES")
_MEGA_N   = _noukb_n("Stroke_MEGASTROKE")
_MAHAJAN_N = _noukb_n("T2D_Mahajan")
_FINNGEN_N = _noukb_n("CAD_FinnGen")
_CKDGEN_N  = _noukb_n("eGFR_CKDGen")
_MEGA_META = _vcf_sample_meta("Stroke_MEGASTROKE_ebi-a-GCST005838.vcf.gz")
_CKDGEN_META = _vcf_sample_meta("eGFR_CKDGen_ebi-a-GCST003372.vcf.gz")

_CARDIO_N = _noukb_n("CAD_CARDIoGRAM")


def _instrument_n(fname):
    """The single per-SNP N carried by a clumped instrument file (asserted to be constant)."""
    vals = set()
    with io.open(os.path.join(BASE, "data/instruments", fname), encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); i = hdr.index("N")
        for ln in fh:
            try:
                vals.add(int(float(ln.rstrip("\n").split("\t")[i])))
            except (ValueError, IndexError):
                pass
    assert len(vals) == 1, f"S1: {fname} does not carry a single N: {sorted(vals)}"
    return vals.pop()


def _koges_n():
    """KoGES sample size, read from the KoGES instrument files the analysis clumped."""
    import glob
    vals = set()
    for f in glob.glob(os.path.join(BASE, "data/instruments", "KoGES_*.eas.clumped.tsv")):
        with io.open(f, encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t"); i = hdr.index("N")
            for ln in fh:
                try:
                    vals.add(int(float(ln.rstrip("\n").split("\t")[i])))
                except (ValueError, IndexError):
                    pass
    assert len(vals) == 1, f"S1: KoGES instrument files disagree on N: {sorted(vals)}"
    return vals.pop()


_KOGES_N = _koges_n()


def _eas_cad_n():
    """Per-SNP N of the Biobank Japan coronary release, read from its clumped instrument file."""
    vals = set()
    with io.open(os.path.join(BASE, "data/instruments", "CAD.eas.clumped.tsv"), encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); i = hdr.index("N")
        for ln in fh:
            try:
                vals.add(int(float(ln.rstrip("\n").split("\t")[i])))
            except (ValueError, IndexError):
                pass
    assert len(vals) == 1, f"S1: BBJ CAD instruments disagree on N: {sorted(vals)}"
    return vals.pop()


_EAS_CAD_N = _eas_cad_n()

S1_COLS = [("trait","Trait","str"),("role","Role","str"),("anc","Ancestry","str"),
           ("cohort","Cohort / consortium","str"),("acc","Accession","str"),
           ("n","Sample size","str"),("build","Build","str"),("ref","Reference","str")]

# ---------------------------------------------------------------------------------------------
# Round-7 provenance: exposure sample sizes are READ from the instrument files, never typed. Each
# instrument row carries the N its loader assigned when it read the source GWAS, so this is the N
# the analysis used - which is the number Table S1 is for.
def _instr_n_stats(stem):
    """(modal N, max N, n instruments) from data/instruments/<stem>.clumped.tsv.

    Distinct from _instrument_n() above, which takes a FILE NAME and returns a single count."""
    _f = os.path.join(BASE, "data", "instruments", f"{stem}.clumped.tsv")
    if not os.path.exists(_f):
        raise SystemExit(f"S1: no instrument file for {stem}")
    from collections import Counter as _C
    with io.open(_f, encoding="utf-8") as _fh:
        _h = _fh.readline().rstrip("\n").split("\t")
        _j = _h.index("N")
        _ns = []
        for _ln in _fh:
            _t = _ln.rstrip("\n").split("\t")
            if len(_t) > _j:
                try:
                    _ns.append(int(float(_t[_j])))
                except ValueError:
                    pass
    if not _ns:
        raise SystemExit(f"S1: instrument file for {stem} carries no sample size")
    return _C(_ns).most_common(1)[0][0], max(_ns), len(_ns)

def _n_cell(stem, note=""):
    _mode, _max, _k = _instr_n_stats(stem)
    _s = f"{_max:,}" if _mode == _max else f"{_mode:,} (median); up to {_max:,}"
    return _s + (f" {note}" if note else "")

def _n_multi(stems):
    """One cell for a row covering several traits from one release, each with its own N.

    Accepts "stem" or ("label", "stem") so the cell names the TRAIT, not the instrument file."""
    _out = []
    for _s in stems:
        _lab, _stem = _s if isinstance(_s, tuple) else (_s, _s)
        _out.append(f"{_lab} {_instr_n_stats(_stem)[1]:,}")
    return "; ".join(_out)

S1_ROWS = [
 # EUR exposures
 dict(trait="BMI", role="Exposure", anc="European", cohort="GIANT + UK Biobank", acc="GCST006900", n=_n_cell("BMI"), build="hg19", ref="Yengo 2018"),
 dict(trait="SBP", role="Exposure", anc="European", cohort="ICBP + UK Biobank", acc="GCST006624", n=">1,000,000", build="hg19", ref="Evangelou 2018"),
 # PROVENANCE CORRECTION (round 7). This row read "MAGIC / Chen 2021 / 146,806" and carried no
 # accession. The file the pipeline actually loads is data/raw/HbA1c_MAGIC_EUR.h.tsv.gz, downloaded
 # by 02_download.sh from the GWAS Catalog harmonised path 34017140-GCST90014006-EFO_0004541 —
 # i.e. GCST90014006, a UK Biobank HbA1c GWAS (Mbatchou 2021, PMID 34017140, 389,889 individuals),
 # not a MAGIC release. The filename had become the fact. Accession, cohort and N are corrected to
 # the shipped file; the Methods, Supplementary Table S14 and the reference list still name
 # MAGIC/Chen 2021 and must be reconciled by the manuscript owner (see the S1 footnote).
 dict(trait="HbA1c", role="Exposure", anc="European", cohort="UK Biobank", acc="GCST90014006", n="389,889", build="hg19", ref="Mbatchou 2021"),
 dict(trait="HDL, LDL, TC, TG", role="Exposure", anc="European", cohort="Global Lipids Genetics Consortium",
      acc="no accession; GLGC 2021 ancestry-specific (European) release, csg.sph.umich.edu/willer/public/glgc-lipids2021",
      n=_n_multi(["HDL", "LDL", "TC", "TG"]), build="hg19", ref="Graham 2021"),
 dict(trait="ApoB", role="Exposure", anc="European", cohort="UK Biobank", acc="GCST90025952", n="435,744", build="hg19/hg38", ref="Barton 2021"),
 # EUR outcomes
 dict(trait="Coronary artery disease", role="Outcome", anc="European", cohort="CARDIoGRAMplusC4D + UKB (Aragam)", acc="GCST90132314", n="181,522 cases / 1,165,690", build="hg19", ref="Aragam 2022"),
 dict(trait="Heart failure", role="Outcome", anc="European", cohort="HERMES (Henry)", acc="GCST90728695", n="139,533 cases", build="hg19", ref="Henry 2025"),
 dict(trait="Stroke", role="Outcome", anc="European", cohort="GIGASTROKE", acc="GCST90104539", n="73,652 cases / 1,234,808 controls", build="hg19", ref="Mishra 2022"),
 dict(trait="Type 2 diabetes", role="Outcome", anc="European", cohort="DIAGRAM (Xue)", acc="GCST006867", n="62,892 cases / 596,424", build="hg19", ref="Xue 2018"),
 dict(trait="eGFR", role="Outcome", anc="European", cohort="CKDGen + UKB (Stanzick)", acc="GCST90103634", n=_n_cell("eGFR"), build="hg19", ref="Stanzick 2021"),
 dict(trait="Chronic kidney disease", role="Outcome", anc="European", cohort="CKDGen (Wuttke)", acc="GCST008065", n="41,395 cases / 439,303", build="hg19", ref="Wuttke 2019"),
 # EAS exposures / cohorts
 dict(trait="BMI", role="Exposure", anc="East Asian", cohort="Biobank Japan (Akiyama)", acc="GCST004904", n=_n_cell("BMI.eas"), build="hg19", ref="Akiyama 2017"),
 dict(trait="SBP, HbA1c, HDL, LDL, TC, TG", role="Exposure", anc="East Asian", cohort="Biobank Japan (Kanai)", acc="NBDC hum0014",
      n=_n_multi([("SBP", "SBP.eas"), ("HbA1c", "HbA1c.eas"), ("HDL", "HDL.eas"),
                  ("LDL", "LDL.eas"), ("TC", "TC.eas"), ("TG", "TG.eas")]), build="hg19", ref="Kanai 2018"),
 dict(trait="Type 2 diabetes", role="Outcome", anc="East Asian", cohort="AGEN (Spracklen)", acc="GCST010118", n="77,418 cases / 433,540", build="hg19", ref="Spracklen 2020"),
 dict(trait="Coronary artery disease", role="Outcome", anc="East Asian", cohort="Biobank Japan (Ishigaki)",
      acc="no accession; Biobank Japan SAIGE release (file BBJ_CAD.txt.gz). The GWAS Catalog entry for this "
          "publication, GCST90013687, is a larger release (29,319 cases / 183,134 controls) and is NOT the file analysed",
      n=f"{_EAS_CAD_N:,} (per-SNP N in the shipped file; case count not carried)", build="hg19", ref="Ishigaki 2020"),
 dict(trait="HF", role="Outcome", anc="East Asian", cohort="Biobank Japan", acc="GCST90668009", n="16,251 cases / 197,577 controls", build="hg38", ref="Enzan 2025"),
 dict(trait="Stroke", role="Outcome", anc="East Asian", cohort="GIGASTROKE EAS stratum", acc="GCST90104545", n="19,032 cases / 237,242 controls", build="hg38", ref="Mishra 2022"),
 dict(trait="CKD", role="Outcome", anc="East Asian", cohort="Biobank Japan (cross-population atlas)", acc="GCST90018602", n="2,117 cases / 174,345 controls", build="hg38", ref="Sakaue 2021"),
 dict(trait="BMI, SBP, lipids, HbA1c, T2D, CAD, HF, stroke, CKD", role="Exposure / Outcome", anc="East Asian (Taiwan)", cohort="Taiwan Precision Medicine Initiative (TPMI)", acc="controlled access", n="hospital-based", build="hg38", ref="TPMI data-access statement"),
 # C045(b): "~211,000" was the KoGES cohort-profile total from Kim & Han 2017. The release analysed
 # here is the KoGES PheWeb GWAS (phenocode-KoGES_*.tsv.gz), whose N is read back from the clumped
 # instrument files the analysis used.
 dict(trait="BMI, waist, DM, lipids, SBP, HbA1c", role="Exposure / Outcome", anc="East Asian (Korea)",
      cohort="KoGES (population); PheWeb GWAS release analysed", acc="controlled access",
      n=f"{_KOGES_N:,}", build="hg19", ref="Kim & Han 2017"),
 dict(trait="Waist circumference", role="Outcome (no-overlap population comparator)", anc="East Asian (Taiwan)", cohort="Taiwan Biobank (population)", acc="controlled access", n="80,465", build="hg38", ref="Taiwan Biobank data-access statement"),
 # UKB-free sensitivity outcomes
 dict(trait="Coronary artery disease", role="Sensitivity outcome (UKB-free)", anc="European", cohort="CARDIoGRAMplusC4D (1000G)",
      acc="ieu-a-7",
      n=f"60,801 cases / 123,504 controls ({_CARDIO_N[2]:,.0f} total)",
      build="hg19", ref="Nikpay 2015"),
 # C048(c). The three rows below printed the case counts of the FULL (UK-Biobank-inclusive) releases.
 # What the UKB-free files carry is a sample size per SNP, and for HERMES it varies; the case counts
 # are not in these releases at all, so they are not asserted here.
 dict(trait="Heart failure", role="Sensitivity outcome (UKB-free)", anc="European", cohort="HERMES (no UKB)",
      acc="no accession; HERMES UK-Biobank-free release (file HF_HERMES_noUKB.tsv.gz)",
      n=f"per-SNP N: median {_HERMES_N[1]:,.0f} (range {_HERMES_N[0]:,.0f}–{_HERMES_N[2]:,.0f}); "
        f"case count not carried in the release file",
      build="hg19", ref="Shah 2020"),
 dict(trait="Stroke", role="Sensitivity outcome (UKB-free)", anc="European", cohort="MEGASTROKE (no UKB)", acc="GCST005838",
      n=f"{float(_MEGA_META['TotalCases']):,.0f} cases / {float(_MEGA_META['TotalControls']):,.0f} controls "
        f"({_MEGA_N[2]:,.0f} total)",
      build="hg19", ref="Malik 2018"),
 dict(trait="eGFR", role="Sensitivity outcome (UKB-free)", anc="European", cohort="CKDGen 2016", acc="GCST003372",
      n=f"{float(_CKDGEN_META['TotalControls']):,.0f} (continuous trait; source-file total)",
      build="hg19", ref="Pattaro 2016"),
 dict(trait="Type 2 diabetes", role="Sensitivity outcome (UKB-free)", anc="European", cohort="DIAMANTE (no UKB)",
      acc="no accession; DIAMANTE UK-Biobank-free release (file T2D_Mahajan_noUKB_rsid.txt)",
      n=f"{_MAHAJAN_N[2]:,.0f} total (constant N assigned at extraction; the release file carries no "
        f"sample-size column, and no case count)",
      build="hg19", ref="Mahajan 2018"),
 dict(trait="CAD, HF, T2D", role="Sensitivity outcome (independent)", anc="European", cohort="FinnGen R12 (endpoints I9_CHD, I9_HEARTFAIL, T2D)",
      acc="no accession; FinnGen R12 public release, r12.finngen.fi",
      n=f"{_FINNGEN_N[2]:,.0f} (R12 data-freeze total; endpoint-specific effective N not used)",
      build="hg19", ref="Kurki 2023"),
]
# The MEGASTROKE cases+controls in the VCF header must reconcile with the per-SNP N the pipeline
# read, and the GIGASTROKE cell must reconcile with the constant lib_formats assigns to that node:
# a hand-typed count that does not add up is exactly the class this table shipped.
assert float(_MEGA_META["TotalCases"]) + float(_MEGA_META["TotalControls"]) == _MEGA_N[2], (
    "S1: MEGASTROKE VCF cases+controls (%s+%s) != the per-SNP N used (%s)"
    % (_MEGA_META["TotalCases"], _MEGA_META["TotalControls"], _MEGA_N[2]))
# The CARDIoGRAM GWAS-VCF carries no TotalCases/TotalControls (the IEU record types it "Continuous"),
# so the published split is reconciled against the per-SNP N the pipeline read instead.
assert 60801 + 123504 == _CARDIO_N[2], (
    "S1: the CARDIoGRAM case/control split no longer sums to the per-SNP N in the shipped file (%s)"
    % _CARDIO_N[2])
assert 73652 + 1234808 == _instrument_n("Stroke.clumped.tsv"), (
    "S1: the GIGASTROKE case/control counts (GWAS Catalog GCST90104539) no longer sum to the N the "
    "pipeline assigns to the Stroke node (%s)" % _instrument_n("Stroke.clumped.tsv"))
ws1 = write_sheet(wb, "S1_data_sources",
    "Table S1 | GWAS data sources. Exposure and outcome summary statistics used to construct the European and East "
    "Asian CKM networks and the UK-Biobank-free sensitivity analysis. Every source carries either a repository "
    "accession or an explicit statement that no accession exists, with the distribution route named. Sample sizes "
    "for the UK-Biobank-free sensitivity providers are the per-SNP sample sizes carried by the files actually "
    "analysed, which are the UK-Biobank-free subsets and therefore smaller than the case counts reported for the "
    "full releases in the cited papers.",
    S1_COLS, S1_ROWS,
    note="TPMI, Taiwan Biobank and KoGES are available under controlled access. Reference keys map to the reference "
         "list in the manuscript References section. Two provenance notes. (i) The European HbA1c exposure file is "
         "GWAS Catalog GCST90014006, a UK Biobank HbA1c GWAS (Mbatchou 2021), not a MAGIC release: the accession, "
         "cohort and sample size in this row are taken from the file the pipeline loads, and the HbA1c exposure "
         "therefore shares UK Biobank participants with the UK-Biobank-containing outcomes. (ii) The East Asian "
         "coronary file is a Biobank Japan SAIGE release whose per-SNP N differs from the GWAS Catalog entry for the "
         "same publication, so no accession is claimed for it.")
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


def add_bh_q(rows, pkey="ivw_p", qkey="bh_q"):
    """Attach the Benjamini-Hochberg adjusted P to each row, over this table's own P column."""
    qs, m = bh_q([to_num(r.get(pkey)) if isinstance(to_num(r.get(pkey)), float) else None
                  for r in rows])
    for r, q in zip(rows, qs):
        r[qkey] = q
    return rows, m


_S2_ROWS, _S2_M = add_bh_q(add_i2(read_csv("forward_local_edges.csv")))
assert _S2_M == len(_S2_ROWS) == 132, \
    "S2: the BH family must be the full 132-edge table (got %d P-values over %d rows)" % (_S2_M, len(_S2_ROWS))
_S2_BONF = 0.05 / _S2_M
_S2_N_BONF = sum(1 for r in _S2_ROWS if isinstance(to_num(r.get("ivw_p")), float)
                 and to_num(r["ivw_p"]) < _S2_BONF)
_S2_N_Q05 = sum(1 for r in _S2_ROWS if isinstance(r.get("bh_q"), float) and r["bh_q"] < 0.05)
write_sheet(wb, "S2_EUR_network",
    "Table S2 | Full European-ancestry bidirectional MR network: all 132 directed exposure–outcome edges. "
    "IVW is the primary estimate; MR-Egger, weighted median, Egger intercept, Cochran Q, the derived "
    "heterogeneity I² = (Q − df)/Q, the Benjamini–Hochberg adjusted P (q) over this table, and liability-scale "
    "Steiger direction are reported for each edge. Effect "
    "sizes: per-SD (continuous exposure) or per-log-OR (binary).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("ivw_b","IVW β","b"),("ivw_se","IVW SE","se"),("ivw_p","IVW P","p"),
     ("bh_q","BH FDR (q)","p"),
     ("egger_b","Egger β","b"),("egger_p","Egger P","p"),
     ("egger_intercept","Egger int.","b"),("egger_intercept_p","Egger int. P","p"),
     ("wm_b","WM β","b"),("wm_p","WM P","p"),
     ("Q","Cochran Q","b"),("Q_p","Q P","p"),("I2","I² (%)","int"),
     ("steiger_correct","Steiger correct","str"),("steiger_p","Steiger P","p")],
    _S2_ROWS,
    note="Bonferroni threshold for the network = 0.05/%d = %s. BH FDR (q) is the Benjamini–Hochberg "
         "step-up adjusted P (scipy.stats.false_discovery_control, method='bh') computed over the IVW P "
         "column of this table — family size m = %d, the full edge set, so the family is exactly the table "
         "shown. It is reported because the Methods state that a Benjamini–Hochberg FDR over this table is "
         "provided as a supplement; no claim in the paper rests on it. %d edges clear the Bonferroni "
         "threshold and %d have q < 0.05; network membership and the staging graph use the Bonferroni "
         "threshold alone, and q informs only which findings are labelled exploratory. Steiger direction for "
         "the network table uses the quantitative-trait approximation; the liability-scale recomputation for "
         "the edges carrying directional claims is Table S11. I², heterogeneity index; WM, weighted median."
         % (_S2_M, _sci(_S2_BONF), _S2_M, _S2_N_BONF, _S2_N_Q05))

# ---------- S3: MVMR cascade ----------
# C019. The caption stated a criterion of "F > 12". The Methods threshold is the conventional 10
# (Sanderson 2019); 12.2 is the smallest value OBSERVED across these models, not a rule. The two are
# now separated and the observed minimum is read from the table rather than typed.
_CASC = read_csv("mvmr_cascade.csv")
_CASC_F = sorted((float(r["cond_F"]), f"{r['exposure']}->{r['outcome']}")
                 for r in _CASC if to_num(r.get("cond_F")) is not None
                 and isinstance(to_num(r.get("cond_F")), float))
assert _CASC_F, "S3 caption: no conditional F values in mvmr_cascade.csv"
assert _CASC_F[0][0] > 10, (
    "S3 caption claims every cascade model clears conditional F = 10, but the minimum is %.2f (%s)"
    % _CASC_F[0])
write_sheet(wb, "S3_MVMR_cascade",
    "Table S3 | Multivariable MR adjudicating cascade versus common-driver mechanism. Each disease outcome is "
    "conditioned jointly on its candidate upstream exposures; the direct effect is the MVMR-IVW estimate and the "
    "total effect is the univariable estimate. Conditional instrument strength is judged against the conventional "
    "threshold of 10 (Sanderson et al. 2019); every model here clears it, the smallest observed value being "
    "%.1f (%s)." % (_CASC_F[0][0], _CASC_F[0][1]),
    [("model","Model","str"),("outcome","Outcome","str"),("exposure","Exposure","str"),("nsnp","nSNP","int"),
     ("direct_b","Direct β","b"),("direct_se","Direct SE","se"),("direct_p","Direct P","p"),
     ("total_b","Total β","b"),("total_p","Total P","p"),("cond_F","Conditional F","b")],
    _CASC,
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
# C049. The ledger's analysis vocabulary calls the two reverse-effect verdicts DISCORDANT and
# DISCORDANT_CAVEATED. The manuscript calls the same two things "directionally supported reverse
# association" and "statistically supported but pleiotropy-caveated", and "discordant" invites the
# reading that the staging order was contradicted, which is the opposite of what these verdicts say
# (every reverse effect runs from a stage-4 disease back into a stage-2 trait). The workbook is the
# reader-facing surface, so it is renamed here, once, and applied to BOTH sheets that print a
# verdict. The analysis CSVs keep the original tokens — main Figure 1 maps them to colours and
# legend labels — so scripts/verify_rebuttal_r{2,3,4}.py were updated to assert this vocabulary on
# the workbook and to map the CSV tokens through the same table.
VERDICT_RENAME = {"DISCORDANT": "REVERSE_SUPPORTED",
                  "DISCORDANT_CAVEATED": "REVERSE_CAVEATED",
                  "DISC_CAVEATED": "REVERSE_CAVEATED"}


def rename_verdicts(rows, keys):
    seen = set()
    for r in rows:
        for k in keys:
            if k in r and r[k] in VERDICT_RENAME:
                r[k] = VERDICT_RENAME[r[k]]
            if k in r and r[k]:
                seen.add(r[k])
    assert not (seen & set(VERDICT_RENAME)), \
        "verdict rename missed a value: %s" % sorted(seen & set(VERDICT_RENAME))
    return rows


_LEDGER = rename_verdicts(read_csv("staging_ledger_native.csv"), ["verdict"])
_LEDGER_VC = {}
for _r in _LEDGER:
    _LEDGER_VC[_r["verdict"]] = _LEDGER_VC.get(_r["verdict"], 0) + 1
assert set(_LEDGER_VC) <= {"CONCORDANT", "INDETERMINATE", "REVERSE_SUPPORTED", "REVERSE_CAVEATED"}, \
    "S5 verdict vocabulary is not the documented one: %s" % sorted(_LEDGER_VC)
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
     ("legacy_mde","80% MDE (reference only)","b"),("legacy_tost_p","Fixed-bound TOST P (reference only)","b")],
    _LEDGER,
    note="Verdict vocabulary. CONCORDANT = the reverse 95%% confidence interval lies entirely inside the "
         "negligibility bound, or the reverse edge is Steiger-wrong-direction (a confounded or index-event "
         "signal, not reverse causation). INDETERMINATE = the interval spans both zero and the bound, so a "
         "null and a material reverse effect cannot be separated. REVERSE_SUPPORTED = a directionally "
         "supported reverse association: the reverse edge excludes zero, exceeds the bound, is "
         "Steiger-correct and carries no directional pleiotropy. REVERSE_CAVEATED = statistically supported "
         "but pleiotropy-caveated: as REVERSE_SUPPORTED, except that the reverse edge's MR-Egger intercept "
         "is non-zero (P < 0.05), so the reverse signal is real while its magnitude is pleiotropy-inflated. "
         "Counts here: %d concordant, %d indeterminate, %d reverse-supported (HF→T2D, Stroke→T2D) and %d "
         "reverse-caveated (CAD→SBP index-event; CAD→T2D feedback-with-caveat, consistent with the "
         "pleiotropy-robust re-estimation of the same reverse edges — a different estimator set on the same "
         "data, not an independent replication). All four reverse effects point into a lower-stage trait as "
         "feedback or an index-event effect, not reverse stage progression. BMI→Stroke is flagged "
         "borderline: its reverse CI limit (+0.0505) sits 0.0005 from the 0.05 SD bound, so that verdict "
         "would flip under trivial re-rounding. The last two columns are retained for reference only and "
         "are not used to adjudicate any verdict."
         % (_LEDGER_VC.get("CONCORDANT", 0), _LEDGER_VC.get("INDETERMINATE", 0),
            _LEDGER_VC.get("REVERSE_SUPPORTED", 0), _LEDGER_VC.get("REVERSE_CAVEATED", 0)))

# ---------- S5b: staging-ledger bound sensitivity ----------
# The negligibility bounds in S5 are clinical judgements. This sheet is the answer to "what if you
# had chosen differently": every verdict across a sweep, produced by scripts/58_ledger_bound_sensitivity.py.
_bs = rename_verdicts(read_csv("staging_ledger_bound_sensitivity.csv"),
                      ["very tight", "tight", "PRIMARY", "loose", "very loose"])
write_sheet(wb, "S5b_ledger_bound_sensitivity",
    "Table S5b | Bound sensitivity of the staging ledger. Each transition's verdict under five negligibility-bound "
    "scenarios, from very tight (odds ratio 1.01 / 0.5 mmHg / 0.02 SD) to very loose (odds ratio 1.20 / 5 mmHg / "
    "0.20 SD); the PRIMARY column is the rule reported in Table S5. The concordant/indeterminate split is "
    "bound-conditional, but the falsification claim is not: under every scenario each reverse effect in the "
    "REVERSE_SUPPORTED and REVERSE_CAVEATED family runs from a stage-4 disease into a stage-2 trait, as feedback or "
    "an index-event effect, and none reverses the staging order.",
    [("transition","Transition","str"),
     ("very tight","Very tight","str"),("tight","Tight","str"),("PRIMARY","PRIMARY","str"),
     ("loose","Loose","str"),("very loose","Very loose","str"),
     ("stability","Stability","str")],
    _bs,
    note="Verdict vocabulary is defined in the footnote to Table S5; the same four labels are used here.")

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
# C057. The Results now cite Table S7b for the East-Asian meta staging numbers, which live in
# staging_eas_meta_exact.csv (written by 77_eas_meta_staging_exact.py). They are appended as a
# separate block rather than retyped into the legend.
_EASSTAGE = read_csv("staging_eas_meta_exact.csv")
assert _EASSTAGE, "S7b: staging_eas_meta_exact.csv is empty"
_EASSTAGE_COLS = [("Meta model", "str"), ("Edges tested", "int"), ("Bonferroni threshold", "p"),
                  ("Retained edges", "int"), ("Nodes carrying a retained edge", "int"),
                  ("Cross-stage edges", "int"), ("Forward", "int"), ("Backward", "int"),
                  ("Concordance", "b"), ("Orderings ≥ observed", "int")]
_EASSTAGE_ROWS = [[r["variant"], r["edges_meta"], r["bonferroni"], r["retained_edges"], r["nodes"],
                   r["cross"], r["forward"], r["backward"], r["concordance"], r["perm_ge"]]
                  for r in _EASSTAGE]
_EASSTAGE_ROWS2 = [[r["variant"], r["perm_total"], r["exact_p"]] for r in _EASSTAGE]
write_sheet(wb, "S7b_EAS_meta",
    "Table S7b | Inverse-variance meta-analysis of the two East Asian hospital cohorts (BBJ + TPMI), fixed- and "
    "random-effects, and the staging result the meta networks support. Fixed-vs-random divergence is reported "
    "because high-heterogeneity edges are borderline at Bonferroni (e.g. T2D→BMI).",
    [("exposure","Exposure","str"),("outcome","Outcome","str"),("nsnp","nSNP","int"),
     ("fx_b","Fixed β","b"),("fx_se","Fixed SE","se"),("fx_p","Fixed P","p"),
     ("rd_b","Random β","b"),("rd_se","Random SE","se"),("rd_p","Random P","p"),
     ("steiger","Steiger correct","str")],
    merge_meta(),
    blocks=[("Staging result for the two East Asian meta networks (exact enumerated null)",
             _EASSTAGE_COLS, _EASSTAGE_ROWS),
            ("Exact permutation null for the same two networks",
             [("Meta model", "str"), ("Distinct stage orderings enumerated", "int"),
              ("Exact P", "b")], _EASSTAGE_ROWS2)],
    note="Concordance is the fraction of cross-stage retained edges pointing forward in stage order. The null is "
         "ENUMERATED, not sampled: every distinct assignment of stage labels to the nodes that carry a retained "
         "edge is scored, so the exact P is the proportion of those orderings whose concordance is at least the "
         "observed one. The two meta networks retain different edge sets, hence different node counts and "
         "different enumeration sizes.")

# ---------- S8: triangulation ----------
# C044. The footnote and Supplementary Methods 9 both promise weighted median, weighted mode,
# penalised weighted median, radial IVW and leave-one-locus-out for the population KoGES edge, plus a
# no-overlap Taiwan Biobank waist comparator. The sheet shipped IVW/WM/Egger only. Every value below
# is PARSED from the analysis output that produced it; nothing is typed.
_H4 = io.open(os.path.join(RES, "h4_leandiabetes.txt"), encoding="utf-8").read()
_TRI = io.open(os.path.join(RES, "triangulation_summary.txt"), encoding="utf-8").read()

_m = re.search(r"Population KoGES DM->BMI:\s*(\d+)\s+harmonised instruments", _H4)
assert _m, "S8: cannot read the KoGES DM->BMI instrument count from h4_leandiabetes.txt"
_KOGES_NSNP = int(_m.group(1))

def _h4_est(label):
    m = re.search(r"^\s*" + re.escape(label) + r"\s+b=([-+\d.eE]+)\s+se=([-+\d.eE]+)\s+p=([-+\d.eE]+)",
                  _H4, re.M)
    assert m, f"S8: estimator '{label}' not found in h4_leandiabetes.txt"
    return m.group(1), m.group(2), m.group(3)

_m = re.search(r"^Effect \(Mod\.2nd\)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+[-+\d.eE]+\s+([-+\d.eE]+)",
               _H4, re.M)
assert _m, "S8: radial-IVW (Mod.2nd) line not found in h4_leandiabetes.txt"
_RADIAL = _m.groups()

_m = re.search(r"KoGES-DM->TWB-Waist\s+\S+\s+\S+\s+(\d+)\s+([-+\d.]+)\s+([\d.eE+-]+)\s+"
               r"([-+\d.]+)\s+([\d.eE+-]+)", _TRI)
assert _m, "S8: the KoGES-DM -> Taiwan-Biobank-waist comparator row is not in triangulation_summary.txt"
_TWB = _m.groups()          # nsnp, ivw_b, ivw_p, wm_b, wm_p

_S8_EST_COLS = [("Estimator", "str"), ("Edge", "str"), ("nSNP", "int"),
                ("Estimate (β)", "b"), ("SE", "se"), ("P", "p"), ("Note", "str")]
_KEDGE = "KoGES_DM→KoGES_BMI"
_S8_EST_ROWS = []
for _lbl, _disp in [("Inverse variance weighted", "Inverse-variance weighted (primary)"),
                    ("MR Egger", "MR-Egger"),
                    ("Weighted median", "Weighted median"),
                    ("Weighted mode", "Weighted mode"),
                    ("Penalised weighted median", "Penalised weighted median")]:
    _b, _se, _p = _h4_est(_lbl)
    _S8_EST_ROWS.append([_disp, _KEDGE, _KOGES_NSNP, _b, _se, _p, ""])
_S8_EST_ROWS.append(["Radial IVW (modified 2nd-order weights)", _KEDGE, _KOGES_NSNP,
                     _RADIAL[0], _RADIAL[1], _RADIAL[2], "7 radial Q-outliers flagged"])
# h4_leandiabetes.txt prints P to four decimals, so two estimators come back as "0.0000". Printing a
# P of exactly zero in a published table is false. The weighted-median P exists at full precision in
# the canonical triangulation table, so read it from there; the penalised weighted median has no
# full-precision source anywhere in results/, so it is reported as a bound rather than invented.
_KTRI = next(r for r in read_csv("koges_triangulation.csv")
             if r["exposure"] == "KoGES_DM" and r["outcome"] == "KoGES_BMI")
for _row in _S8_EST_ROWS:
    if _row[0] == "Weighted median" and _row[1] == _KEDGE:
        assert abs(float(_row[3]) - float(_KTRI["wm_b"])) < 5e-4, \
            "S8: the weighted-median estimate disagrees between h4_leandiabetes.txt and koges_triangulation.csv"
        _row[3], _row[5] = _KTRI["wm_b"], _KTRI["wm_p"]
    elif float(_row[5]) == 0.0:
        _row[5] = "< 0.0001"
        _row[6] = (_row[6] + "; " if _row[6] else "") + "P below the precision of the deposited output"
_S8_EST_ROWS.append(["Inverse-variance weighted", "KoGES_DM→TWB_WAIST", _TWB[0], _TWB[1], "", _TWB[2],
                     "Taiwan Biobank population comparator; no sample overlap with KoGES"])
_S8_EST_ROWS.append(["Weighted median", "KoGES_DM→TWB_WAIST", _TWB[0], _TWB[3], "", _TWB[4],
                     "Taiwan Biobank population comparator; MC4R instrument absent"])

_LOO = re.findall(r"^\s*drop (rs\d+)\s+b=([-+\d.eE]+)\s+p=([-+\d.eE]+)", _H4, re.M)
assert len(_LOO) == _KOGES_NSNP, \
    "S8: %d leave-one-out rows for %d instruments" % (len(_LOO), _KOGES_NSNP)
_m = re.search(r"LOO range \[([-+\d.]+),\s*([-+\d.]+)\]", _H4)
assert _m, "S8: the leave-one-out range line is not in h4_leandiabetes.txt"
_LOO_RANGE = _m.groups()
_S8_LOO_COLS = [("Instrument dropped", "str"), ("Edge", "str"), ("nSNP", "int"),
                ("IVW β", "b"), ("P", "p")]
_S8_LOO_ROWS = [[s, _KEDGE, _KOGES_NSNP - 1, b, p] for s, b, p in _LOO]

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
    blocks=[("Robust-estimator suite for the population KoGES DM→BMI edge, and the no-overlap "
             "Taiwan Biobank waist comparator (Supplementary Methods 9)", _S8_EST_COLS, _S8_EST_ROWS),
            ("Leave-one-locus-out: IVW re-estimated dropping each KoGES DM instrument in turn",
             _S8_LOO_COLS, _S8_LOO_ROWS)],
    note="DM→BMI IVW point estimate (−0.043) matches the hospital TPMI edge (−0.046): the edge does not vanish "
         "in a population cohort, so it is not a demonstrable ascertainment artifact; a causal effect is not "
         "established (IVW-random null, balanced pleiotropy). The pleiotropy-robust estimators in the block above "
         "stay negative, while the conservative two-cohort Hartung–Knapp meta-analysis is non-significant "
         "(95%% CI [−0.40, +0.26], P = 0.22; k = 2, I² = 0.76). Leave-one-locus-out does not flip the sign: the "
         "IVW estimate ranges over [%s, %s] across the nine single-instrument drops, and the two adiposity loci "
         "(MC4R rs6567160, 12q24 rs2074356) are the drops that move it most negative. Estimator standard errors "
         "are not reported for the Taiwan Biobank comparator because the deposited summary does not carry them. "
         "The forward BMI→DM positive control is intact."
         % (_LOO_RANGE[0], _LOO_RANGE[1]))

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
# C046. "Comparable" was computed from the EXPOSURE units alone, so 16 edges whose OUTCOME is HbA1c
# or eGFR were flagged comparable here while Table S15 and the Methods declare those two traits
# readable on sign only. The rule now requires BOTH sides (scripts/14_standardise.py) and the column
# says which quantity it is about.
_STD = read_csv("edges_standardised.csv")
_STD_BAD = [r["edge"] for r in _STD
            if str(r.get("comparable")) == "True" and (r["exposure"] in ("HbA1c", "eGFR")
                                                       or r["outcome"] in ("HbA1c", "eGFR"))]
assert not _STD_BAD, (
    "S10: %d edges are still flagged magnitude-comparable although HbA1c/eGFR is on one side "
    "(Table S15 calls those sign-only). Re-run 14_standardise.py. Offenders: %s"
    % (len(_STD_BAD), _STD_BAD[:8]))
write_sheet(wb, "S10_standardised",
    "Table S10 | Cross-ancestry scale standardisation. Effect sizes on native and per-SD scales; only SBP required "
    "rescaling (EUR per-mmHg → per-SD, SD = 19.3 mmHg), computed from the full-precision IVW estimate rather than "
    "the rounded native value shown. Magnitudes are comparable only where BOTH the exposure and the outcome are on "
    "a common scale in the two ancestries: HbA1c and eGFR are read on sign concordance alone (Table S15), so every "
    "edge with either trait on either side is flagged not magnitude-comparable rather than rescaled.",
    [("edge","Edge","str"),("exp_units_eur","EUR units","str"),("exp_units_eas","EAS units","str"),
     ("out_units","Outcome units","str"),("EUR_b","EUR β (native)","b"),("EAS_b","EAS β (native)","b"),
     ("EUR_b_perSD","EUR β (per-SD)","b"),("EAS_b_perSD","EAS β (per-SD)","b"),
     ("comparable","Magnitude-comparable (exposure and outcome)","str"),("scale_note","Note","str")],
    _STD,
    note="A 'False' in the magnitude-comparable column does not mean the edge is uninformative: direction "
         "concordance is still read for every edge (Table S6). It means only that the two effect sizes are on "
         "different scales, so their ratio or difference carries no interpretation.")

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
# C017. The caption put HF→CKD in the same "robust" tier as BMI→CKD and SBP→CKD. Its IVW P does not
# clear the network Bonferroni threshold (0.05/132), so it is demoted to nominal and the P is read
# from the sheet's own source file rather than described.
_CKDN = read_csv("network_ckd.csv")
_HFCKD = next(r for r in _CKDN if r["exposure"] == "HF" and r["outcome"] == "CKD")
assert float(_HFCKD["ivw_p"]) > _S2_BONF, \
    "S13 caption calls HF→CKD nominal, but its IVW P now clears the network Bonferroni threshold"
write_sheet(wb, "S13_CKD_node",
    "Table S13 | Binary chronic-kidney-disease (CKD) node (Wuttke 2019, GCST008065), a kidney-damage endpoint "
    "distinct from the quantitative eGFR node. Risk-factor→CKD edges fall in graded tiers: robust (BMI, SBP: "
    "IVW and weighted median significant past the network Bonferroni threshold, clean Egger intercept); nominal "
    "(HF: IVW P = %s, above the network Bonferroni threshold of %s, with a concordant weighted median but only "
    "%s instruments); and pleiotropy-suspect (HDL, TG: weighted median null, significant Egger intercept). "
    "CKD→cardiovascular edges are null but underpowered (~24 instruments)."
    % ("%.3f" % float(_HFCKD["ivw_p"]), _sci(_S2_BONF), _HFCKD["nsnp"]),
    NET16_COLS, _CKDN,
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
 dict(item="HbA1c (Mbatchou 2021, UK Biobank, GCST90014006)", role="Exposure", anc="EUR", overlap="Exposure-side", note="The European HbA1c exposure is a UK Biobank GWAS, so HbA1c edges into UK-Biobank-containing outcomes carry exposure-side overlap, as BMI and SBP do. The outcome-side UKB-free sensitivity analysis (Table S12) is unaffected; exposure-side overlap remains a stated limitation."),
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
 dict(trait="HbA1c", eur="Mbatchou 2021 (UK Biobank, GCST90014006)", eas="BBJ Kanai 2018, rank-inverse-normal", cmp="No — sign only (the two releases are not on a declared common scale)"),
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
# C043. Methods and Supplementary Methods 2 both cite Table S17 for the FORMAL between-subtype
# heterogeneity test (difference, Q, and P across five assumed correlations). The sheet held only
# the per-edge IVW/Egger/WM rows; the test itself sat unrendered in hf_subtype_heterogeneity.csv.
_HET = read_csv("hf_subtype_heterogeneity.csv")
assert _HET, "S17: hf_subtype_heterogeneity.csv is empty"
_RHOS = [k for k in _HET[0] if k.startswith("p_rho_")]
assert _RHOS == ["p_rho_0.00", "p_rho_0.05", "p_rho_0.10", "p_rho_0.25", "p_rho_0.50"], \
    "S17: the rho grid in hf_subtype_heterogeneity.csv changed: %s" % _RHOS
_HET_COLS = ([("Exposure", "str"), ("HFpEF β", "b"), ("HFpEF SE", "se"),
              ("HFrEF β", "b"), ("HFrEF SE", "se"), ("Difference (HFpEF − HFrEF)", "b")]
             + [("P, ρ = %s" % k.split("_")[-1], "p") for k in _RHOS]
             + [("Cochran Q", "b")])
_HET_ROWS = [[r["exposure"], r["b_hfpef"], r["se_hfpef"], r["b_hfref"], r["se_hfref"], r["diff"]]
             + [r[k] for k in _RHOS] + [r["Q"]] for r in _HET]
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
    blocks=[("Formal between-subtype heterogeneity test (HFpEF versus HFrEF), per exposure: the "
             "difference in effects, Cochran's Q, and the two-sample P across five assumed "
             "correlations between the two subtype estimates", _HET_COLS, _HET_ROWS)],
    note="T2D->HFrEF is IVW-positive but attenuates under weighted-median (~0.046) and MR-Egger (~0.022); report "
         "the direction and range, not a firm magnitude. The heterogeneity block above is the test the Methods and "
         "Supplementary Methods 2 cite this table for. The two subtype outcomes share one control set (cases are "
         "disjoint), which induces a positive correlation between their estimates; because a positive correlation "
         "SHRINKS the variance of the difference, the rho = 0 column is the conservative test and rho = 0.05 is the "
         "value implied by the shared controls under the Lin-Sullivan formula. The rho grid is 0, 0.05, 0.10, 0.25 "
         "and 0.50; Q and its P are computed at rho = 0. A difference significant at rho = 0 therefore holds a "
         "fortiori, and a non-significant difference does not establish equivalence. Bonferroni over the %d "
         "exposures tested is 0.05/%d = %s." % (len(_HET), len(_HET), _sci(0.05 / len(_HET))))

# ---------- S18: WHRadjBMI (adiposity distribution) ----------
# C047. The sheet carried WHRadjBMI->niHF / niHFpEF / niHFrEF rows. Those non-ischaemic heart-failure
# subtype outcomes are defined nowhere in the package and were RETIRED as circular in an earlier
# review round (excluding ischaemic cases conditions directly on CAD) — Table S17 uses the
# all-aetiology subtypes for exactly that reason. Dropped here and from the source CSV. The filter is
# applied at the generator so a regenerated CSV cannot reintroduce them silently.
RETIRED_HF_OUTCOMES = {"niHF", "niHFpEF", "niHFrEF"}
_WHR_ALL = read_csv("network_whradjbmi.csv")
_WHR = [r for r in _WHR_ALL if r["outcome"] not in RETIRED_HF_OUTCOMES]
assert _WHR, "S18: every WHRadjBMI row was filtered out"
write_sheet(wb, "S18_WHRadjBMI",
    "Table S18 | Adiposity distribution independent of overall mass (WHRadjBMI; Pulit 2019, 694,649 European, "
    "combined-sex). Central fat distribution is consistent with a causal effect on CAD and T2D but not on heart "
    "failure: no WHRadjBMI–HF association was detected. Because the exposure is conditioned on BMI and the "
    "estimates are heterogeneous, this does not establish whether overall mass or fat distribution carries the "
    "risk. Steiger filtering was not applied to these edges.",
    NET12_COLS, _WHR,
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
    note="Margin low/high are the endpoints across the prevalence grid; for the 20 continuous-to-continuous edges prevalence cannot enter the calculation, so low = margin = high and they are flagged prevalence-invariant (step 61c). This analysis is one-sided: it "
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
    "54-edge Steiger-directed graph. HDL→SBP is intra-stage, "
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

# C103. The Contents sheet is the workbook's index; a hand-maintained index drifts from the book it
# indexes and nothing noticed. Refuse to save a workbook whose data sheets and Contents rows disagree.
_DATA_SHEETS = [s for s in wb.sheetnames if s != "Contents"]
_TOC_SHEETS = [t[0] for t in TOC]
assert len(_TOC_SHEETS) == len(set(_TOC_SHEETS)), "Contents lists a sheet twice"
_missing = [s for s in _DATA_SHEETS if s not in _TOC_SHEETS]
_extra = [s for s in _TOC_SHEETS if s not in _DATA_SHEETS]
assert not _missing and not _extra, (
    "Contents does not index the workbook: %d data sheets, %d Contents rows; missing %s; stale %s"
    % (len(_DATA_SHEETS), len(_TOC_SHEETS), _missing, _extra))

path = os.path.join(OUT, "Supplementary_Tables.xlsx")
wb.save(path)

# Global rule 8. wb.properties covers dc:creator / cp:lastModifiedBy, but openpyxl also writes
# docProps/app.xml advertising its own name and version, which is third-party tool provenance in a
# shipped file. Rewrite that one part in place; everything else in the package is byte-identical.
import zipfile, shutil, tempfile
_APP = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
        '<Application>Microsoft Excel</Application><AppVersion>16.0300</AppVersion></Properties>')
_tmp = tempfile.mktemp(suffix=".xlsx")
with zipfile.ZipFile(path) as _zin, zipfile.ZipFile(_tmp, "w", zipfile.ZIP_DEFLATED) as _zout:
    for _it in _zin.infolist():
        _data = _APP.encode("utf-8") if _it.filename == "docProps/app.xml" else _zin.read(_it.filename)
        _zout.writestr(_it, _data)
shutil.move(_tmp, path)
assert wb.properties.creator == wb.properties.lastModifiedBy == "Bertrand Chin-Ming Tan"

print("wrote", path, "with", len(wb.sheetnames), "sheets:", wb.sheetnames)
print("Contents indexes all %d data sheets; BH FDR family m=%d; S18 dropped %d retired niHF rows"
      % (len(_DATA_SHEETS), _S2_M, len(_WHR_ALL) - len(_WHR)))
