# E4 — Sample-overlap matrices, EUR and EAS (2026-07-19)

Reviewer C5/E4/M11: two-sample MR assumes little exposure–outcome sample overlap; overlap biases
toward the confounded observational association. This records the overlap structure and what has
already been done to address it.

## EUR — UK Biobank is the shared cohort to track
Overlap between an exposure GWAS and an outcome GWAS arises mainly through shared UK Biobank (UKB)
participants. Per-GWAS UKB membership (⚠ = author-confirm against the source's cohort description):

| GWAS (role) | Includes UKB? |
|---|---|
| BMI Yengo 2018 (exposure) | Yes (GIANT+UKB) |
| SBP Evangelou 2018 (exposure) | Yes (UKB+ICBP) |
| Lipids GLGC 2021 (exposure) | ⚠ partial (meta includes UKB) |
| HbA1c Chen 2021 MAGIC (exposure) | ⚠ confirm |
| T2D Xue 2018 GCST006867 (exp/outcome) | **Yes** (DIAGRAM+UKB) |
| CAD Aragam 2022 GCST90132314 (outcome) | Yes |
| HF Henry/HERMES GCST90728695 (outcome) | Yes |
| Stroke GIGASTROKE Mishra 2022 (outcome) | Yes |
| eGFR Stanzick 2021 CKDGen (outcome) | **Yes** (CKDGen+UKB, ~36% UKB; the project keeps a UKB-free eGFR = Pattaro 2016 for this reason — Supp Table S1) |
| CKD Wuttke 2019 CKDGen (outcome) | ⚠ mostly non-UKB |

**Overlap classification:** any exposure×outcome pair where BOTH sides include UKB is a
partial-overlap MR (most risk-factor→CAD/HF/T2D/stroke headline edges, and the eGFR edges via
Stanzick). Only the **CKD (Wuttke)** kidney outcome is closer to non-overlapping.

**Already mitigated:** the outcome-side **UKB-free sensitivity** (Milestone 8; `results/sensitivity_noukb_*`)
re-estimated the headline edges against UKB-free outcome GWAS (CARDIoGRAM CAD, Mahajan-2018 T2D,
MEGASTROKE stroke, etc.); all 13 forward edges retained sign and Bonferroni significance against the
UKB-free main providers. FinnGen R12 (an additional fully independent check) supplies CAD/HF/T2D
outcomes only — the two stroke edges rest on MEGASTROKE. Direction is robust to outcome-side overlap. Exposure-side UKB (BMI, SBP) remains and
is disclosed as a limitation; a fully UKB-free exposure re-run is the residual follow-up. WHRadjBMI
(Path B) shares UKB with the CAD/HF/T2D outcomes → its magnitudes carry the same caveat.

## EAS — overlap is by cohort, and a zero-overlap design now exists
| Design | Overlap |
|---|---|
| BBJ exposure → BBJ outcome (`network_eas_edges`) | **Complete** (one-sample / within-cohort MR) |
| TPMI exposure → TPMI outcome (`network_tpmi_*`) | **Complete** (one-sample / within-cohort MR) |
| KoGES DM → KoGES BMI (triangulation) | **Complete** (within-cohort, explicitly one-sample) |
| **BBJ exposure → TPMI outcome** (`network_crosscohort_bbj_tpmi`, Path B) | **None** — Japanese vs Taiwanese, no shared participants = clean two-sample |

**Already mitigated:** the within-cohort EAS networks are now flagged as overlapping-sample MR, and
the **cross-cohort BBJ→TPMI** analysis (E6/H6) provides the zero-overlap check: the headline cascade
(BMI→CAD/HF, LDL/TC→CAD, HbA1c→T2D, SBP→Stroke) replicates zero-overlap, while BMI→T2D — significant
only within TPMI — does NOT, exactly the edge where within-cohort overlap inflated the estimate.

## Net
Direction of the headline edges is robust to overlap on both ancestries (UKB-free outcomes in EUR;
zero-overlap cross-cohort in EAS). Magnitudes remain overlap-sensitive and are reported as such.
Residual: a fully UKB-free EUR exposure re-run (optional follow-up).
