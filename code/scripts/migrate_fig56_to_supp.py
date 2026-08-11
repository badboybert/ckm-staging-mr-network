# -*- coding: utf-8 -*-
"""ONE-TIME structural migration: main Figures 5 and 6 become Supplementary Figures S6 and S7.

Round-2 review item B-2 #15, author decision 2026-07-25: keep four main figures and move Figures 5
and 6 wholesale to the supplement. The legibility argument (36 panels at 85-170 mm) is real and
Cardiovascular Diabetology sets no figure cap for Research articles, so this is presentation only --
no estimate, panel or number changes.

Everything here is mechanical and asserted:
  * every "Figure 5x" / "Figure 6x" cross-reference in the authored sources is rewritten by ONE
    regex pass, and the script FAILS if any survives;
  * the two legend blocks are MOVED, not retyped, so no legend text exists in two places;
  * the render scripts are pointed at the new file names, so the panels cannot be regenerated under
    the old names.

Idempotent: re-running after a successful migration is a no-op that still passes the assertions.
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

import io, os, re, sys, shutil

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

BASE = P4_BASE
FIG = os.path.join(BASE, "figures")

MAP = {"5": "6", "6": "7"}          # main Figure N -> Supplementary Figure S(MAP[N])

def rd(p):
    return open(p, encoding="utf-8").read()

def wr(p, t):
    open(p, "w", encoding="utf-8").write(t)

# ---------------------------------------------------------------- 1. cross-references
# Panels may be written 6b, 6b,d or 5a-c. "main Figure 6" collapses to the supplementary form.
XREF = re.compile(r"\b(?:main\s+)?Figures?\s*([56])((?:[a-h])(?:\s*[,\u2013-]\s*[a-h])*)?\b")

def sub_xref(m):
    panels = (m.group(2) or "").replace(" ", "")
    return f"Supplementary Figure S{MAP[m.group(1)]}{panels}"

SECTION_FILES = ["ABSTRACT.md", "INTRODUCTION.md", "METHODS.md", "RESULTS.md", "DISCUSSION.md",
                 "CONCLUSIONS.md", "COVER_LETTER.md", "STROBE_MR_CHECKLIST.md",
                 "FIGURE_LEGENDS.md", "SUPPL_FIGURES.md", "FRONT_MATTER.md"]

# ---------------------------------------------------------------- 2. move the legend blocks first
leg_p = os.path.join(FIG, "FIGURE_LEGENDS.md")
leg = rd(leg_p)
blocks = {}
for n in ("5", "6"):
    m = re.search(rf"^\*\*Figure {n} \|.*?(?=^\*\*Figure \d \||\Z)", leg, flags=re.S | re.M)
    if m:
        blocks[n] = m.group(0).rstrip() + "\n"
        leg = leg[:m.start()] + leg[m.end():]

if blocks:
    # header sentence in FIGURE_LEGENDS names the count and the build scripts
    leg = leg.replace("_Six main figures;", "_Four main figures;")
    leg = re.sub(r"`fig1\.R`[^\n`]*`fig6\.R`", "`fig1.R`–`fig4.R`", leg)
    leg = re.sub(r"\n{3,}", "\n\n", leg).rstrip() + "\n"
    wr(leg_p, leg)

sup_p = os.path.join(FIG, "SUPPL_FIGURES.md")
sup = rd(sup_p)
if blocks:
    # The retired-S6 note must lose its number: an S6 now exists, and a note calling S6 "retired"
    # beside a live S6 is precisely the cross-file contradiction the gates cannot see.
    sup = re.sub(
        r"_\(The former \*\*Supplementary Figure S6\*\*.*?\)_\n?",
        "_(The staging permutation-null data (`05_source_data/suppfig6_*.csv`, from "
        "`figures/suppfig6_staging_null.py`) is plotted in main **Figure 4b**, beside the "
        "concordance bars it explains; it is not a separate supplementary figure.)_\n",
        sup, flags=re.S)
    for n, block in blocks.items():
        newn = MAP[n]
        block = re.sub(rf"^\*\*Figure {n} \|", f"**Supplementary Figure S{newn} |", block, flags=re.M)
        sup = sup.rstrip() + "\n\n" + block
    # header sentence in SUPPL_FIGURES names counts and script ranges
    sup = sup.replace("Five supplementary figures rendered in the same house style as the six main figures",
                      "Seven supplementary figures rendered in the same house style as the four main figures")
    sup = sup.replace("(`SupplFig1`–`SupplFig5`)", "(`SupplFig1`–`SupplFig7`)")
    sup = sup.replace("Build scripts `suppfig1.R`–`suppfig5.R`;",
                      "Build scripts `suppfig1.R`–`suppfig5.R` and, for S6 and S7, `fig5.R` and `fig6.R`;")
    wr(sup_p, sup)

# ---------------------------------------------------------------- 3. rewrite cross-references
changed = {}
for f in SECTION_FILES:
    p = os.path.join(FIG, f)
    if not os.path.exists(p):
        continue
    t = rd(p)
    new, n = XREF.subn(sub_xref, t)
    # collapse the double prefix a nested match could produce
    new = new.replace("Supplementary Supplementary Figure", "Supplementary Figure")
    if new != t:
        wr(p, new)
        changed[f] = n

# ---------------------------------------------------------------- 4. point the renderers at the new names
for src, out in (("fig5.R", "SupplFig6"), ("fig6.R", "SupplFig7")):
    p = os.path.join(FIG, src)
    t = rd(p)
    old = re.search(r'save_fig\((\w+), "(Figure[56])"', t)
    if old:
        t = t.replace(f'"{old.group(2)}"', f'"{out}"')
        t = ("# 2026-07-25: this figure MOVED to the supplement (round-2 item B-2 #15). Panels and\n"
             f"# numbers are unchanged; only the output name is, so it renders as {out}.\n") + t
        wr(p, t)

# remove the now-orphaned main-figure files so nothing can ship under the retired names
removed = []
for n in ("5", "6"):
    for e in ("png", "pdf"):
        q = os.path.join(FIG, f"Figure{n}.{e}")
        if os.path.exists(q):
            os.remove(q); removed.append(os.path.basename(q))

# ---------------------------------------------------------------- 5. assert the migration is complete
problems = []
for f in SECTION_FILES:
    p = os.path.join(FIG, f)
    if not os.path.exists(p):
        continue
    t = re.sub(r"<!--.*?-->", "", rd(p), flags=re.S)
    for m in XREF.finditer(t):
        frag = " ".join(t[max(0, m.start() - 50):m.start() + 60].split())
        problems.append(f"{f}: surviving main-figure reference ...{frag}...")

leg = rd(leg_p)
main_ids = re.findall(r"^\*\*Figure (\d+) \|", leg, flags=re.M)
if main_ids != ["1", "2", "3", "4"]:
    problems.append(f"FIGURE_LEGENDS.md now defines figures {main_ids}, expected ['1','2','3','4']")
sup = rd(sup_p)
sup_ids = re.findall(r"^\*\*Supplementary Figure S(\d+) \|", sup, flags=re.M)
if sup_ids != ["1", "2", "3", "4", "5", "6", "7"]:
    problems.append(f"SUPPL_FIGURES.md now defines S{sup_ids}, expected S1-S7")

print(f"cross-references rewritten : {changed}")
print(f"legend blocks moved        : {sorted(blocks)}  ->  S{[MAP[k] for k in sorted(blocks)]}")
print(f"orphaned figure files removed: {removed}")
print(f"main figures defined       : {main_ids}")
print(f"supplementary figures      : S{sup_ids}")
if problems:
    print(f"\n===== MIGRATION INCOMPLETE — {len(problems)} problem(s) =====")
    for p_ in problems:
        print("  x", p_)
    sys.exit(1)
print("\n===== MIGRATION COMPLETE AND ASSERTED =====")
