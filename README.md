# Adiposity-centred causal architecture of cardiovascular–kidney–metabolic traits: a cross-ancestry Mendelian randomization study of coronary disease and heart failure

Analysis code and frozen derived data for the manuscript above (submitted to *Cardiovascular
Diabetology*).

The study estimates all **132** non-self directed relationships among 12
cardiovascular–kidney–metabolic (CKM) traits by two-sample Mendelian randomization,
tests whether the resulting graph is concordant with the American Heart Association CKM staging
order, adjudicates cascade versus common-driver mechanism by multivariable MR and formal mediation,
and evaluates directional portability in two East Asian hospital biobanks and one population cohort.
**58** edges pass the network-wide Bonferroni threshold and
**54**
survive Steiger gating.

## What is and is not here

**Here:** every analysis, figure and table script, and every *derived* file they read — the edge
tables, the null distributions, the mediation and portability outputs, the per-figure source data
shipped as Additional file 5, and the Supplementary Tables workbook.

**Not here: raw GWAS summary statistics.** They are third-party data and are not redistributed.
European statistics are available from their original consortia and the GWAS Catalog under the
accessions listed in the manuscript Methods; the TPMI, Taiwan Biobank and KoGES statistics are
available under their respective controlled-access procedures. The upstream extraction and clumping
steps therefore cannot be re-run from this deposit alone — everything downstream of them can.

## Layout

```
code/scripts/    the numbered analysis pipeline (01–75) plus the build, QA and gate scripts
code/figures/    the figure renderers (fig1–fig6, suppfig1–suppfig6) and the locked shared theme
code/results/    derived analysis outputs the figures and tables read
source_data/     the per-panel source data shipped with the manuscript as Additional file 5
supplementary_tables/  Supplementary Tables S1–S22 (29 worksheets)
renv.lock        pinned R package versions
```

## Running it

Paths resolve from the environment, falling back to each script's own location, so no configuration
is needed to regenerate the display items:

```bash
export CKM_P4_BASE="$(pwd)/code"      # Windows: set CKM_P4_BASE=%CD%\code
pip install -r requirements.txt
Rscript -e 'renv::restore()'          # or install the packages in R_packages.txt
Rscript code/figures/fig1.R           # …fig2.R, fig3.R, fig4.R, fig5.R, fig6.R, suppfig1–6.R
python code/scripts/qa_manuscript.py  # the manuscript number-vs-source checks
```

Set `CKM_ROOT` as well only if you have obtained the raw summary statistics and intend to re-run the
upstream extraction steps.

## Provenance and checks

Numbers in the manuscript are not typed; they are derived from these files and asserted at build
time. `code/scripts/gate_package.py` runs Gates A–F over the submission package — cross-file
identity, numbers against source, the rendered figure text layers, cross-representation agreement
between every `.md`/`.docx` twin, and the spelling, typography and callout rules. Each gate has been
mutation-tested: an injected defect must make it fail with the message it names.

## Citation

Please cite the manuscript. This deposit is archived on Zenodo under the concept DOI
[10.5281/zenodo.22813026](https://doi.org/10.5281/zenodo.22813026), which always resolves to the latest
version. Source: https://github.com/badboybert/ckm-staging-mr-network

## Licence

MIT for the code; CC BY 4.0 for the frozen data files. See `LICENSE`.
