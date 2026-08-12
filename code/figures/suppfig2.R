# --- portable roots -----------------------------------------------------------------------------
# Injected by scripts/build_repo.py; see the note in the Python header. Set CKM_P4_BASE to override.
P4_BASE <- Sys.getenv("CKM_P4_BASE")
if (!nzchar(P4_BASE)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
  .here <- if (length(.self)) dirname(.self[1]) else "."
  P4_BASE <- normalizePath(file.path(.here, ".."), winslash = "/", mustWork = FALSE)
}
P4_ROOT <- Sys.getenv("CKM_P4_ROOT"); if (!nzchar(P4_ROOT)) P4_ROOT <- dirname(P4_BASE)
CKM_ROOT <- Sys.getenv("CKM_ROOT");   if (!nzchar(CKM_ROOT)) CKM_ROOT <- dirname(P4_ROOT)
SHARED_LIB <- Sys.getenv("CKM_SHARED_LIB")
if (!nzchar(SHARED_LIB)) SHARED_LIB <- file.path(CKM_ROOT, "_shared")
# ------------------------------------------------------------------------------------------------
# SupplFig2 — Canonical-recovery forest: the pipeline reproduces established causal estimates.
# Source: results/forward_local_edges.csv (IVW beta, 95% CI = beta +/- 1.96*SE).
# Positive controls (LDL->CAD ~+0.5, BMI->CAD, BMI->T2D, HbA1c->T2D, SBP->Stroke),
# supporting canonical edges, and the negative control (HDL->CAD; marginal inverse = confounded).
source(file.path(P4_BASE, "figures", "fig_setup.R"))

raw <- rd("forward_local_edges.csv")

spec <- data.frame(
  exp  = c("LDL","BMI","HbA1c","BMI","SBP","TC","TG","BMI","SBP","HDL"),
  out  = c("CAD","T2D","T2D","CAD","Stroke","CAD","CAD","HF","CAD","CAD"),
# 2026-07-26: "Negative control" -> "Correlated-marker control". HDL->CAD is not a null
# control; it is a CONFOUNDED correlated marker whose marginal inverse estimate is expected,
# and round 2 replaced the bare label on Figure 3 and in the workbook. This panel kept it,
# because every PDF checker in the project read only `(...) Tj` and this legend is drawn with
# `[...] TJ` arrays - it was invisible, not overlooked.
  role = c("Positive control","Positive control","Positive control","Positive control","Positive control",
           "Canonical","Canonical","Canonical","Canonical","Correlated-marker control"),
  expect = c("LDL causal (statin / PCSK9 RCTs)","Adiposity strongly causal","Glycaemia definitional",
             "Adiposity causal","BP causal (per mmHg)","Atherogenic cholesterol causal",
             "TG-rich lipoprotein causal","Adiposity causal","BP causal (per mmHg)",
             # Round 5 (M2). Two defects in one label. "Non-causal … confounded" asserts more than
             # the design gives — MVMR shows attenuation under conditioning, not non-causality, and
             # the manuscript says so. And "Fig 6" is a stale main-figure reference: the 2026-07-25
             # migration moved main Figure 6 to Supplementary Figure S7, but the rewrite only ever
             # touched .md sources, so this literal inside a RENDERER survived — and Gate F21, which
             # catches dangling main-figure callouts, reads text surfaces and not PDF text layers.
             "Marginal inverse association attenuates\nafter multivariable conditioning\n(Supplementary Figure S7)"),
  stringsAsFactors = FALSE)

d <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i){
  r <- raw[raw$exposure==spec$exp[i] & raw$outcome==spec$out[i], ]
  data.frame(edge=sprintf("%s → %s", spec$exp[i], spec$out[i]), role=spec$role[i], expect=spec$expect[i],
             b=r$ivw_b[1], lo=r$ivw_b[1]-1.96*r$ivw_se[1], hi=r$ivw_b[1]+1.96*r$ivw_se[1],
             p=r$ivw_p[1], nsnp=r$nsnp[1], stringsAsFactors=FALSE)
}))
# top-to-bottom order: positive controls, canonical, negative control (numeric y for clean annotation)
d$role  <- factor(d$role, levels=c("Positive control","Canonical","Correlated-marker control"))
d$ypos  <- rev(seq_len(nrow(d)))   # first row at top
ROLE_COL <- c("Positive control"="#2E7D32", "Canonical"="#0072B2", "Correlated-marker control"="#D55E00")

# DERIVE from the committed edge table (was hardcoded 0.05/111, the superseded denominator).
BONF <- 0.05/nrow(read.csv(file.path(RES, "forward_local_edges.csv")))
XT <- 1.52   # left edge of the expectation-text column
YTOP <- max(d$ypos)

sf2 <- ggplot(d, aes(b, ypos, colour=role)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  annotate("segment", x=XT-0.08, xend=XT-0.08, y=0.4, yend=YTOP+0.6, linewidth=0.3, colour="grey85") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.55, orientation="y") +
  geom_point(aes(shape=ifelse(p<BONF,"sig","ns")), size=2.3, fill="white", stroke=0.7) +
  geom_text(aes(label=sprintf("%.2f", b)), vjust=-0.95, size=FS/.pt-1.1, colour="black") +
  geom_text(aes(x=XT, label=expect), hjust=0, vjust=0.5, size=FS/.pt-1.2, colour="grey20", lineheight=0.9) +
  annotate("text", x=XT, y=YTOP+0.85, label="Canonical expectation", hjust=0, fontface=2,
           size=FS/.pt-0.9, colour="grey30") +
  scale_colour_manual(values=ROLE_COL, name=NULL) +
  scale_shape_manual(values=c(sig=16, ns=1), guide="none") +
  scale_y_continuous(breaks=d$ypos, labels=d$edge, expand=expansion(add=c(0.6, 1.1))) +
  labs(x="IVW causal estimate (per-SD / per-unit exposure; log-OR disease outcome)", y=NULL,
       title="Pipeline recovers canonical causal estimates and controls") +
  coord_cartesian(xlim=c(-0.45, 3.15), clip="off") +
  theme_ckm(legend="top") +
  theme(legend.position="top", plot.margin=margin(11,8,4,4))

save_fig(sf2, "SupplFig2", 178, 108)
cat("SF2 done:", nrow(d), "edges (", sum(d$p<BONF), "Bonferroni-sig )\n")
