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
# Figure 2 — Adiposity reaches heart failure directly; T2D reaches it predominantly via CAD.
# a mechanism schematic · b MVMR direct-effect forest · c formal mediation (proportion mediated, E7)
# d robustness across estimators · e conditional-F instrument strength · f HF-subtype specificity (H1)
#
# 2026-07-22 REVISION (132-edge rewrite):
#   - panel c: the old total->direct slopegraph held HARDCODED literals and has been REPLACED by the
#     formal network-MR mediation result, read from data/fig2_mediation.csv (derived + self-asserted
#     by prep_newpanels.py from results/formal_mediation.txt). No hand-typed effect sizes remain.
#   - panel a: the annotation "T2D->HF direct = 0 (fully mediated)" was FALSE under E7 and is now
#     "predominantly via CAD" (proportion mediated ~68-79%, wide CI; direct undetectable, not zero).
#   - panel f: NEW HF-subtype panel (Henry 2024 HERMES non-ischemic), the specificity control.
source(file.path(P4_BASE, "figures", "fig_setup.R"))
D   <- file.path(FIG, "data")
mv  <- rd("mvmr_cascade.csv"); rb <- rd("mvmr_robust.csv")
med <- read.csv(file.path(D, "fig2_mediation.csv"), stringsAsFactors = FALSE)
hs  <- rd("network_hfsubtypes_allcause.csv")   # Enzan 2025 all-aetiology HFpEF/HFrEF (round-2 #5)
mcov <- read.csv(file.path(FIG, "..", "results", "mediation_covariance_sensitivity.csv"), stringsAsFactors = FALSE)

OUTMAP <- c(CAD_full="CAD", HF_via_CAD="HF", Stroke_full="Stroke")
m <- mv[mv$model %in% names(OUTMAP), ]
m$outcome_f <- factor(OUTMAP[m$model], levels=c("CAD","HF","Stroke"))
m$lo <- m$direct_b - 1.96*m$direct_se; m$hi <- m$direct_b + 1.96*m$direct_se
m$surv <- ifelse(m$direct_p < 0.05, "survives", "attenuated to null")
m <- m[order(m$outcome_f, m$direct_b), ]
m$row <- factor(paste(m$outcome_f, m$exposure), levels=paste(m$outcome_f, m$exposure))

# ---------------------------------------------------------------- 2a. mechanism schematic
nd <- data.frame(node=c("Adiposity\n(BMI)","T2D","LDL / SBP","CAD","HF"),
                 x=c(0,1,1,2,3), y=c(0,1.1,-1.1,0,0),
                 col=c(STAGE_COL["1"],STAGE_COL["2"],STAGE_COL["2"],STAGE_COL["4"],STAGE_COL["4"]))
seg <- data.frame(
  x =c(0,   0,   0,   1,   1,   2,   0,   0),
  y =c(0,   0,   0,   1.1, -1.1,0,   0,   0),
  xe=c(1,   1,   2,   2,   2,   3,   2,   3),
  ye=c(1.1, -1.1,0,   0,   0,   0,   0,   0),
  kind=c("cascade","cascade","cascade","cascade","cascade","cascade","direct","direct"))
pm_t2d <- med$pm_product[med$exposure=="T2D"]; pm_t2d_d <- med$pm_diff[med$exposure=="T2D"]
pm_bmi <- 100 - med$pm_product[med$exposure=="BMI"]
pa <- ggplot() +
  geom_curve(data=seg[seg$kind=="cascade",], aes(x=x,y=y,xend=xe,yend=ye), curvature=0.12,
             linewidth=0.4, colour="grey55", arrow=arrow(length=unit(4,"pt"),type="closed")) +
  geom_curve(data=seg[seg$kind=="direct",], aes(x=x,y=y,xend=xe,yend=ye), curvature=-0.34,
             linewidth=1.1, colour=POS_COL, arrow=arrow(length=unit(5,"pt"),type="closed")) +
  annotate("segment", x=1.1, y=1.0, xend=2.9, yend=0.12, linewidth=0.5, colour="grey72", linetype="21") +
  geom_point(data=nd, aes(x,y), shape=21, size=15, fill=nd$col, colour="white", stroke=0.5) +
  geom_text(data=nd, aes(x,y,label=node), size=FS/.pt-1.3, fontface=2, lineheight=0.85) +
  # NB: this label must be drawn AFTER the nodes (they are opaque and would paint over it) and placed
  # in the empty band BELOW the CAD-HF axis — x~1.3 is clipped by T2D, x~1.6 is clipped by CAD.
  annotate("label", x=2.05, y=-0.85,
           label="no CAD-independent T2D→HF effect\n(compatible with CAD mediation)",
           colour="grey35", size=FS/.pt-1.7, lineheight=0.9, fill="white") +
  # 2026-07-26: "~80% DIRECT" -> "~80% not via CAD". MVMR estimates a component conditional on the
  # included exposures; it does not exclude unmeasured mediators, so an unqualified "DIRECT" on the
  # panel over-reads what the model identifies. The prose already says "not mediated through CAD".
  annotate("label", x=1.55, y=2.45,
           label=sprintf("adiposity → HF: ~%d%% not via CAD (common driver)", pm_bmi),
           colour=POS_COL, size=FS/.pt-1, fontface=2, fill="white", label.size=0) +
  coord_cartesian(xlim=c(-0.4,3.5), ylim=c(-1.7,2.65), clip="off") +
  theme_void(base_family=FONT) + theme(plot.margin=margin(4,6,4,6))

# ---------------------------------------------------------------- 2b. direct-effect forest
pb <- ggplot(m, aes(direct_b, row, colour=surv)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.1) +
  geom_text(aes(label=sprintf("%.2f", direct_b)), vjust=-0.9, size=FS/.pt-1, colour="black") +
  geom_text(aes(x=Inf, label=sprintf("F=%.0f", cond_F)), hjust=1.05, size=FS/.pt-1.3, colour="grey45") +
  # Upper expansion must clear the value label drawn above the topmost point of EACH facet. At the
  # default 0.6 the top label in the HF facet ("0.42") and in the Stroke facet ("0.16") were sliced
  # by the facet boundary, in both the PNG and the PDF.
  scale_y_discrete(labels=setNames(m$exposure, m$row),
                   expand=ggplot2::expansion(add=c(0.6, 1.0))) +
  scale_colour_manual(values=c("survives"=POS_COL, "attenuated to null"=NULL_COL), name=NULL) +
  facet_grid(outcome_f ~ ., scales="free_y", space="free_y", switch="y") +
  labs(x="MVMR direct effect (conditional)", y=NULL) +
  coord_cartesian(xlim=c(-0.05, 0.62)) +
  theme_ckm(legend="bottom") +
  theme(strip.placement="outside", strip.text.y.left=element_text(angle=0, face=2),
        panel.spacing=unit(4,"pt"))

# ---------------------------------------------------------------- 2c. formal mediation via CAD (E7)
# Proportion of the exposure->HF effect carried by CAD, with Monte-Carlo 95% CI.
# The honest contrast: BMI tightly ~20% mediated (=> ~80% DIRECT); T2D ~68-79% mediated but with a
# CI spanning 100% (the total effect is small) => "predominantly", never "only".
# BMI is tightly ~20% mediated (=> ~80% DIRECT) and covariance-robust; T2D's mediated proportion is
# imprecise AND covariance-DEPENDENT (round-2 #1): the exposure->CAD (a) and total exposure->HF share
# the T2D instruments, so corr(a,total) shrinks the product interval. Two intervals are drawn for T2D
# from results/mediation_covariance_sensitivity.csv - path independence and strong shared-instrument
# correlation - both spanning full mediation. No single fraction is reported. All values READ.
.bmi  <- med[med$exposure=="BMI", ]
.indep <- mcov[mcov$exposure=="T2D" & mcov$rho_ab==0 & mcov$rho_total_direct==0 & mcov$rho_a_total==0, ]
.strg  <- mcov[mcov$exposure=="T2D" & mcov$rho_a_total==0.9, ]
mc <- data.frame(
  exp_f = factor(c("BMI","T2D (independence)","T2D (strong corr.)"),
                 levels=c("T2D (strong corr.)","T2D (independence)","BMI")),
  # Round-7 C032: the legend and the Results both say no single T2D point estimate is reportable,
  # but the panel drew one (the product-method point, twice). The two T2D rows are interval-only.
  pt    = c(.bmi$pm_product[1], NA, NA),
  lo    = c(.bmi$pm_lo[1], .indep$pm_product_lo[1], min(.strg$pm_product_lo)),
  hi    = c(.bmi$pm_hi[1], .indep$pm_product_hi[1], max(.strg$pm_product_hi)),
  grp   = c("BMI","T2D","T2D"))
mc$lab <- sprintf("%.0f–%.0f%%", mc$lo, mc$hi)
# The two T2D rows carry no point estimate (C032), so the interval label is anchored at the middle of
# the interval instead of at the (missing) point; otherwise the label disappears with the point.
mc$labx <- ifelse(is.na(mc$pt), (mc$lo + mc$hi)/2, mc$pt)
pc <- ggplot(mc, aes(pt, exp_f)) +
  annotate("rect", xmin=100, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey94") +
  geom_vline(xintercept=100, linetype="dashed", linewidth=0.4, colour="grey55") +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi, colour=grp), width=0, linewidth=0.6, orientation="y") +
  geom_point(aes(colour=grp), size=2.6) +
  geom_text(aes(x=labx, label=lab), vjust=-1.1, size=FS/.pt-1.3, colour="black") +
  annotate("text", x=100, y=3.5, label=">100%: unstable estimate", size=FS/.pt-1.6, colour="grey45",
           vjust=0, hjust=0.5) +
  # moved below the BMI row: at y=3 it sat ON the BMI point and the marker overprinted the text.
  annotate("text", x=8, y=2.66, label=sprintf("~%.0f%% not via CAD", 100 - .bmi$pm_product[1]),
           size=FS/.pt-1.6, colour="grey35", hjust=0) +
  scale_colour_manual(values=c(BMI=POS_COL, T2D=NEG_COL), guide="none") +
  scale_x_continuous(breaks=c(0,25,50,75,100,150)) +
  coord_cartesian(xlim=c(0,155), ylim=c(0.55,3.7), clip="off") +
  # 2026-07-26 (round-4 Tier 1): this axis title was CLIPPED after "(95%" at the panel edge, so
  # the panel read as unfinished. Shortened rather than widened, because the composite width is
  # fixed by the journal column. A text-layer check could not see it: the string was complete in
  # the PDF, it was the DEVICE that cut it.
  # Still clipped after one shortening: the axis title is centred under a narrow panel in the
  # composite, so LENGTH is the only lever. "(95% CI)" is dropped here because the legend already
  # states that the bars are Monte-Carlo 95% CIs - the information is not lost, only relocated.
  labs(x="% of HF effect mediated by CAD", y=NULL) +
  theme_ckm() + theme(axis.text.y=element_text(face=2, size=8))

# ---------------------------------------------------------------- 2d. robustness across estimators
rbk <- rb[(rb$model=="CAD_full" & rb$exposure=="BMI") | (rb$model=="HF_via_CAD" & rb$exposure %in% c("BMI","T2D")) |
          (rb$model=="Stroke_full" & rb$exposure=="LDL"), ]
rbk$edge <- paste0(rbk$exposure, "→", OUTMAP[rbk$model])
rl <- reshape(rbk[,c("edge","ivw_b","egger_b","qhet_b")], direction="long",
              varying=c("ivw_b","egger_b","qhet_b"), v.names="b", times=c("IVW","Egger","Q-minimisation"), timevar="est")
rl$edge <- factor(rl$edge, levels=rev(unique(rbk$edge)))
pd <- ggplot(rl, aes(b, edge, shape=est, colour=est)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  scale_shape_manual(values=c(IVW=16, Egger=17, "Q-minimisation"=15), name=NULL) +
  scale_colour_manual(values=c(IVW="#333333", Egger="#0072B2", "Q-minimisation"="#009E73"), name=NULL) +
  labs(x="Direct effect", y=NULL) +
  coord_cartesian(xlim=c(-0.05,0.5)) + theme_ckm(legend="bottom")

# ---------------------------------------------------------------- 2e. conditional F
m$edge <- paste0(m$exposure, "→", m$outcome_f)
pe <- ggplot(m, aes(reorder(edge, cond_F), cond_F)) +
  geom_col(fill="#8DA0CB", width=0.7) +
  geom_hline(yintercept=10, linetype="dashed", colour=POS_COL, linewidth=0.4) +
  annotate("text", x=10.3, y=10.5, label="F = 10", colour=POS_COL, size=FS/.pt-1, hjust=-0.05, vjust=0.4) +
  geom_text(aes(label=sprintf("%.0f", cond_F)), hjust=-0.2, size=FS/.pt-1.2) +
  coord_flip(ylim=c(0,62)) +
  labs(x=NULL, y="Conditional F-statistic") +
  theme_ckm() + theme(axis.text.y=element_text(size=8))

# ---------------------------------------------------------------- 2f. HF-subtype reach (all-aetiology)
# Enzan 2025 all-aetiology HFpEF/HFrEF (replaces the circular non-ischemic panel; round-2 #5). CAD
# reaches BOTH subtypes (ischemia not excluded); adiposity reaches both with no between-subtype
# difference DETECTED (P = 0.74 — absence of evidence, not equivalence); T2D reaches HFrEF but not
# HFpEF. Cross-ancestry outcome -> read on direction, not magnitude (legend). Filled = P<0.05.
# 2026-07-26: the on-panel label used to read "T2D reaches HFrEF, not HFpEF", i.e. it asserted a
# SUBTYPE DIFFERENCE. The formal between-subtype test (script 72) supports that contrast only
# NOMINALLY (P = 7.5e-3, below 0.05 but not Bonferroni over 8 exposures); the Results say so and the
# panel did not. The label now states the same thing the prose does. Both P-values are READ from the
# heterogeneity results file so the panel cannot drift from the test.
.het <- rd("hf_subtype_heterogeneity.csv")
stopifnot(nrow(.het[.het$exposure == "T2D", ]) == 1)
sub <- hs[hs$exposure %in% c("BMI","T2D","CAD","SBP") & hs$outcome %in% c("HFpEF","HFrEF"), ]
sub$lo <- sub$ivw_b - 1.96*sub$ivw_se; sub$hi <- sub$ivw_b + 1.96*sub$ivw_se
sub$sig <- ifelse(sub$ivw_p < 0.05, "P < 0.05", "n.s.")
sub$outcome_f <- factor(sub$outcome, levels=c("HFpEF","HFrEF"),
                        labels=c("HFpEF (all-aetiology)","HFrEF (all-aetiology)"))
sub$exp_f <- factor(sub$exposure, levels=c("SBP","T2D","CAD","BMI"))
pf <- ggplot(sub, aes(ivw_b, exp_f, colour=outcome_f, shape=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5,
                position=position_dodge(width=0.55), orientation="y") +
  geom_point(size=2.3, position=position_dodge(width=0.55)) +
  scale_colour_manual(values=c("HFpEF (all-aetiology)"="#B2182B","HFrEF (all-aetiology)"="#2166AC"), name=NULL) +
  scale_shape_manual(values=c("P < 0.05"=16, "n.s."=1), name=NULL) +
  annotate("text", x=0.30, y=1.6, hjust=0, vjust=0.5, size=FS/.pt-1.6, colour="grey35",
           label=sprintf("T2D: HFrEF detected, HFpEF not\n(between-subtype difference\nnominal only, P = %.1g)",
                         .het$p_rho_0.00[.het$exposure=="T2D"])) +
  # 2026-07-26 (round-4 Tier 1): a single "log-OR per SD" label is inaccurate here — the exposures
  # mix disease liabilities with SBP in mmHg, so no one unit describes every row.
  labs(x="Effect on HF subtype (exposure units differ by row; read direction, not magnitude)", y=NULL) +
  theme_ckm(legend="bottom") + theme(axis.text.y=element_text(face=2))

# ---------------------------------------------------------------- assemble
# tag order = a schematic, b forest, c mediation, d robustness, e cond-F, f subtypes
fig2 <- pa / (pb | pc) / (pd | pe) / pf +
  plot_layout(heights=c(0.58, 1.12, 0.88, 0.86)) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(fig2, "Figure2", 170, 222)
