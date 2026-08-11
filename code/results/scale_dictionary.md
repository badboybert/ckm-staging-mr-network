# E5 — Scale dictionary for cross-ancestry comparison (2026-07-19)

Reviewer C6/E5: raw effect-size similarity is uninterpretable when exposure scaling differs. This
table records, per trait, the exposure scale on each ancestry side and whether magnitudes are
directly comparable. It governs which ancestry comparisons may be shown as magnitudes vs sign/rank.

| Trait | EUR source / scale | EAS source / scale | Comparable? |
|---|---|---|---|
| BMI | Yengo 2018 (GIANT+UKB), per-SD | BBJ Akiyama 2017, RINT (≈per-SD) | **Yes** (both per-SD) |
| SBP | Evangelou 2018, **per mmHg** | BBJ Kanai 2018, RINT (**per-SD**) | **Only after conversion** — multiply EUR by SD_SBP≈19.3 mmHg to reach per-SD; do NOT compare raw |
| LDL / HDL / TC / TG | GLGC 2021, per-SD | BBJ Kanai 2018, RINT (≈per-SD) | **Yes** (both per-SD) |
| HbA1c | Chen 2021 (MAGIC), %-HbA1c units | BBJ Kanai 2018, RINT | **No** — different units/transform; flagged non-comparable in the Methods. Sign only. |
| eGFR | Stanzick 2021 (CKDGen), log(eGFR) per-SD | BBJ Kanai 2018, RINT | **No** — flagged non-comparable in the Methods. Sign only. |
| T2D (as exposure) | Xue 2018, log-OR | AGEN Spracklen 2020, log-OR | **Yes** (both log-OR) |
| CAD / HF / Stroke / CKD (outcomes) | observed-scale log-OR | observed-scale log-OR | **Yes** (both log-OR) |

## Consequences for the figures/claims (apply in the rewrite)
- **Figure 4c** (EUR-vs-EAS raw identity line): remove the raw identity line for any panel mixing
  per-mmHg SBP with per-SD traits, or restrict the identity line to the genuinely per-SD/log-OR
  subset. Use **sign/rank concordance** for HbA1c and eGFR, never magnitude.
- **Figure 5a** (SBP→CKD EUR vs EAS): the two SBP points are on different scales (per-mmHg EUR vs
  per-SD EAS) and must not be juxtaposed as magnitudes. Convert EUR to per-SD (×19.3) before any
  magnitude statement; the formal interaction test (E6) uses the converted per-SD estimate.
- **Any ancestry "magnitude" or "ratio" claim** is valid only for the "Comparable = Yes" rows on a
  common scale; otherwise report sign concordance.
- **HbA1c-involving sign reversals** (BMI→HbA1c, HbA1c→HF, HbA1c→TG) pass the E6 interaction test but
  carry the HbA1c non-comparability caveat; report as sign-level interactions, not magnitude.

⚠ Cohort composition (which consortia each GWAS includes, e.g. UKB) is author-confirmable against
each source's cohort table; see the overlap matrix (E4). SD_SBP=19.3 mmHg is the value already used
in the Methods for the per-SD conversion.
