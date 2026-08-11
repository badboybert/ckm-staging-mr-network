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
# Figure 3 — Directionality and reciprocal feedback in the CKM network.
# a CAD<->HF asymmetry (multi-method) · b CAUSE across 10 edges · c MR-PRESSO raw->corrected
# d triangulation matrix (Steiger / PRESSO / CAUSE) · e CAUSE specificity, correlated-marker negative control
source(file.path(P4_BASE, "figures", "fig_setup.R"))
ca <- rd("cause_headline.csv"); pr <- rd("mrpresso_headline.csv"); nc <- rd("cause_negcontrol_pairs.csv")

# ---------- 3a. CAD<->HF asymmetry (IVW / Egger / WM) ----------
# Every estimate, intercept P and the "diverge N-fold" factor below is DERIVED from the committed
# edge table and asserted, never typed. These six betas and two intercept P-values were previously
# hand-entered literals: correct at the time, but silently detached from the analysis.
.net  <- rd("forward_local_edges.csv")
.row  <- function(e, o) {
  r <- .net[.net$exposure == e & .net$outcome == o, ]
  stopifnot(nrow(r) == 1)
  r
}
.ch <- .row("CAD", "HF"); .hc <- .row("HF", "CAD")
asym <- data.frame(
  edge=rep(c("CAD → HF","HF → CAD"), each=3),
  method=rep(c("IVW","MR-Egger","Weighted median"), 2),
  b=c(.ch$ivw_b, .ch$egger_b, .ch$wm_b, .hc$ivw_b, .hc$egger_b, .hc$wm_b),
  note=c("","","", "","",""))
stopifnot(all(is.finite(asym$b)))
asym$edge <- factor(asym$edge, levels=c("HF → CAD","CAD → HF"))
asym$method <- factor(asym$method, levels=c("IVW","MR-Egger","Weighted median"))
.fold <- max(c(.hc$ivw_b, .hc$egger_b, .hc$wm_b)) / min(c(.hc$ivw_b, .hc$egger_b, .hc$wm_b))
.pfmt <- function(p) if (p < 0.01) sprintf("%.3f", p) else sprintf("%.2f", p)
lab3a <- data.frame(edge=factor(c("CAD → HF","HF → CAD"), levels=c("HF → CAD","CAD → HF")),
                    # two SHORT lines: a single long line overruns the panel and is clipped at both ends
                    x=c(1.0,1.0),
                    txt=c(sprintf("methods concordant\nSteiger ✓ · Egger-int P=%s", .pfmt(.ch$egger_intercept_p)),
                          # %.0f rounded 5.61 up to "6×" while the Results said "fivefold"; one decimal
                          # makes the panel and the prose state the same number.
                          sprintf("methods diverge %.1f×\nSteiger ✗ · Egger-int P=%s", .fold, .pfmt(.hc$egger_intercept_p))))
pa <- ggplot(asym, aes(b, edge, colour=method, shape=method)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_point(size=2.4, position=position_dodge(width=0.55)) +
  geom_text(aes(label=sprintf("%.2f", b)), position=position_dodge(width=0.55), vjust=-1, size=FS/.pt-1.4, show.legend=FALSE) +
  geom_text(data=lab3a, aes(x=x, y=edge, label=txt), inherit.aes=FALSE, hjust=0.5, vjust=2.6,
            size=FS/.pt-1.7, colour="grey30", lineheight=0.9) +
  scale_colour_manual(values=c(IVW="#333333","MR-Egger"="#0072B2","Weighted median"="#009E73"),
                      labels=c("IVW","MR-Egger","W. median"), name=NULL) +
  scale_shape_manual(values=c(IVW=16,"MR-Egger"=17,"Weighted median"=15),
                     labels=c("IVW","MR-Egger","W. median"), name=NULL) +
  labs(x="Causal estimate (log-OR)", y=NULL, title="CAD↔HF estimates are directionally asymmetric") +
  coord_cartesian(xlim=c(-0.1,2.0)) + theme_ckm(legend="bottom")

# ---------- 3b. CAUSE across 10 edges ----------
ca$edge <- paste0(ca$exposure,"→",ca$outcome)
ca$verd <- ifelse(ca$verdict=="CAUSAL","causal model favoured","sharing model favoured")
ca <- ca[order(ca$z_sharing_vs_causal), ]
ca$edge <- factor(ca$edge, levels=ca$edge)
pb <- ggplot(ca, aes(z_sharing_vs_causal, edge, colour=verd)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  annotate("rect", xmin=-Inf, xmax=-1.96, ymin=-Inf, ymax=Inf, fill=CAUSAL_COL, alpha=0.05) +
  geom_point(size=2.2) +
# 2026-07-26 (round-4 Tier 1): the bottom legends of the SIDE-BY-SIDE panels ran into each
# other ("Weighted medi[ar]causal model favoured", "sharing model favou[red]VW"). The legend
# labels are shortened and the canvas is taller; the full wording stays in the legend text.
  scale_colour_manual(values=c("causal model favoured"=CAUSAL_COL, "sharing model favoured"=SHARING_COL),
                      labels=c("causal favoured","sharing favoured"), name=NULL) +
  labs(x="CAUSE Δ-ELPD z (negative → causal model favoured)", y=NULL,
       title="CAUSE model preference, on a continuous scale") +
  annotate("text", x=-8.5, y=9.4, label="causal model favoured", colour=CAUSAL_COL, size=FS/.pt-1.5, hjust=0) +
  theme_ckm(legend="bottom")

# ---------- 3c. MR-PRESSO raw -> corrected ----------
pr$edge <- paste0(pr$exposure,"→",pr$outcome)
pr$distort <- ifelse(pr$distortion_p=="<0.001" | suppressWarnings(as.numeric(pr$distortion_p))<0.05, "distorted","robust")
pr$distort[is.na(pr$distort)] <- "robust"
pc_ <- reshape(pr[,c("edge","raw_b","corrected_b","distort","n_outliers")], direction="long",
               varying=c("raw_b","corrected_b"), v.names="b", times=c("raw","outlier-\ncorrected"), timevar="stage")
pc_$stage <- factor(pc_$stage, levels=c("raw","outlier-\ncorrected"))
pc <- ggplot(pc_, aes(stage, b, group=edge, colour=distort)) +
  geom_line(linewidth=0.6) + geom_point(size=1.9) +
  ggrepel::geom_text_repel(data=pc_[pc_$stage=="outlier-\ncorrected",], aes(label=edge),
     size=FS/.pt-1.3, direction="y", hjust=0, nudge_x=0.1, segment.size=0.2, min.segment.length=0) +
  scale_colour_manual(values=c("distorted"=POS_COL, "robust"=NULL_COL), name=NULL) +
  labs(x=NULL, y="Causal estimate", title="MR-PRESSO: only HF→CAD is outlier-distorted") +
  coord_cartesian(xlim=c(1,2.7)) + theme_ckm(legend="bottom")

# ---------- 3d. the two feedback limbs differ in robustness (H5, CAUSE-independent) ----------
# REPLACES the former hardcoded Steiger/PRESSO/CAUSE triangulation tile-matrix, which asserted
# CAD→T2D as cleanly "robust / causal". H5 shows the two limbs are NOT equivalent: HF→T2D is
# concordant across every estimator with a clean Egger intercept, whereas CAD→T2D loses the
# weighted MODE and carries a significant directional-pleiotropy intercept. Data derived + asserted
# by prep_newpanels.py from results/h5_reverse_pleiotropy.txt.
h5 <- read.csv(file.path(FIG, "data", "fig3_h5_estimators.csv"), stringsAsFactors=FALSE)
# The panel shows the FOUR pleiotropy-robust point estimators. Two rows are deliberately excluded:
# the Egger intercept (a pleiotropy test, annotated separately below) and the MR-Egger SLOPE, added
# to the source data on 2026-07-27 so the reader can check the one estimator the Results report as
# null (+0.39, P = 0.31 for HF→T2D). Its CI is [-0.35, +1.13] and would triple this half-width
# panel's x-range. The exclusion is now EXPLICIT and asserted: the old filter dropped only "Egger
# intercept", so a new estimator row would have fallen out of the factor as NA and vanished from the
# panel with no error.
PANEL_EST <- c("IVW", "Weighted median", "Penalised WM", "Weighted mode")
EXCLUDED  <- c("Egger intercept", "MR-Egger slope")
stopifnot(setdiff(unique(h5$estimator), c(PANEL_EST, EXCLUDED)) |> length() == 0)
sl5 <- h5[h5$estimator %in% PANEL_EST, ]
sl5$lo <- sl5$b - 1.96*sl5$se; sl5$hi <- sl5$b + 1.96*sl5$se
sl5$estimator <- factor(sl5$estimator, levels=rev(PANEL_EST))
stopifnot(!any(is.na(sl5$estimator)))
sl5$sig <- ifelse(sl5$p < 0.05, "P < 0.05", "n.s.")
ic5 <- h5[h5$estimator == "Egger intercept", ]
ic5$lab <- sprintf("Egger intercept P = %s",
                   ifelse(ic5$p < 0.01, sprintf("%.3f", ic5$p), sprintf("%.2f", ic5$p)))
pd <- ggplot(sl5, aes(b, estimator, colour=edge, shape=sig)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.2) +
  geom_text(data=ic5, aes(x=0.62, y=0.62, label=lab), inherit.aes=FALSE,
            hjust=1, size=FS/.pt-1.6, colour="grey35") +
  facet_wrap(~edge, ncol=1, strip.position="right") +
  scale_colour_manual(values=c("CAD→T2D"="#E69F00","HF→T2D"="#2E7D32"), guide="none") +
  scale_shape_manual(values=c("P < 0.05"=16, "n.s."=1), name=NULL) +
  labs(x="Reverse-edge effect on T2D (95% CI)", y=NULL,
       title="Feedback limbs differ: HF→T2D robust, CAD→T2D pleiotropy-caveated") +
  coord_cartesian(xlim=c(-0.25, 0.68)) +
  theme_ckm(legend="bottom") +
  theme(axis.text.y=element_text(size=8), strip.text.y=element_text(angle=0, face=2))

# ---------- 3e. CAUSE specificity, correlated-marker negative control ----------
# The MVMR annotations were hand-typed and drifted from source: "−0.17" is not any model's estimate —
# the full-density 630-instrument conditional HDL→CAD effect is −0.1647, which Supplementary Figure S7d (deriving the
# same number from the same file) correctly prints as −0.16. Both are now read from the CSVs.
.fd  <- rd("mvmr_lipid_fulldensity.csv")
.ap  <- rd("mvmr_apob.csv")
.ldl_dir  <- .fd$direct_b[.fd$exposure=="LDL"][1]
.hdl_dir  <- .fd$direct_b[.fd$exposure=="HDL"][1]
.hdl_marg <- .ap$total_b[.ap$exposure=="HDL"][1]
stopifnot(is.finite(.ldl_dir), is.finite(.hdl_dir), is.finite(.hdl_marg))
sp <- data.frame(edge=c("LDL→CAD\n(positive-control\nassociation)","HDL→CAD\n(correlated-marker\ncalibration association)"),
                 z=c(-5.87,-2.46),
                 mvmr=c(sprintf("direct %+.2f (robust)", .ldl_dir),
                        sprintf("attenuates %.2f→%.2f", .hdl_marg, .hdl_dir)),
                 kind=c("causal model favoured","sharing model favoured"))
sp$hj <- ifelse(sp$z < -4, 0, 1)      # left point grows right, right point grows left
sp$edge <- factor(sp$edge, levels=rev(sp$edge))
pe <- ggplot(sp, aes(z, edge, colour=kind)) +
  geom_vline(xintercept=-1.96, linetype="dashed", colour="grey60", linewidth=0.3) +
  geom_point(size=2.6) +
  geom_text(aes(label=sprintf("CAUSE z=%.2f", z)), vjust=-1.1, size=FS/.pt-1.3, colour="black") +
  # 2026-07-26 (round-4 Tier 1): centred on the point, the LDL label ran off the LEFT edge and
  # rendered as "irect +0.50 (robust)"; left-aligning both then pushed the HDL label off the RIGHT.
  # The anchor is therefore chosen per point - grow away from the nearer edge.
  geom_text(aes(label=mvmr, hjust=hj), vjust=2, size=FS/.pt-1.4, colour="grey35") +
  scale_colour_manual(values=c("causal model favoured"=CAUSAL_COL,"sharing model favoured"="#E69F00"),
                      labels=c("causal favoured","sharing favoured"), name=NULL) +
  annotate("text", x=-1.96, y=2.5, label="CAUSE causal-model\npreference threshold", size=FS/.pt-1.6, colour="grey55", lineheight=0.9) +
  labs(x="CAUSE Δ-ELPD z", y=NULL,
       title="CAUSE favours the causal model for a confounded marker (HDL)") +
  coord_cartesian(xlim=c(-7,0.3)) + theme_ckm(legend="bottom")

# ---------- 3f. per-SNP scatter for the headline edge (instrument-level evidence for CAD->HF) ----------
# Source data exported by suppfig3.R (re-harmonised via TwoSampleMR, IVW reproduces the committed table to <1e-15).
# The full 3-edge diagnostic set (scatter + funnel) remains in Supplementary Fig. S3; here we show the headline edge.
DD3   <- file.path(FIG, "suppfig_data")
sc3   <- read.csv(file.path(DD3, "suppfig3_scatter.csv"), stringsAsFactors=FALSE)
sl3   <- read.csv(file.path(DD3, "suppfig3_slopes.csv"),  stringsAsFactors=FALSE)
sc3   <- sc3[sc3$edge == "CAD -> HF", ]; sl3 <- sl3[sl3$edge == "CAD -> HF", ]
MCOL3 <- c("IVW"="#B2182B", "MR-Egger"="#4393C3", "Weighted median"="#E69F00")
MLAB3 <- c("IVW"="IVW", "MR-Egger"="MR-Egger", "Weighted median"="W. median")
pf <- ggplot(sc3, aes(bx, by)) +
  geom_hline(yintercept=0, linewidth=0.2, colour="grey80") +
  geom_vline(xintercept=0, linewidth=0.2, colour="grey80") +
  geom_errorbar(aes(ymin=by-sy, ymax=by+sy), width=0, linewidth=0.2, colour="grey75") +
  geom_errorbarh(aes(xmin=bx-sx, xmax=bx+sx), height=0, linewidth=0.2, colour="grey75") +
  geom_point(size=0.8, colour="grey30", alpha=0.7) +
  geom_abline(data=sl3, aes(slope=slope, intercept=intercept, colour=method), linewidth=0.6) +
  scale_colour_manual(values=MCOL3, labels=MLAB3, name=NULL) +
  labs(x="SNP effect on CAD (log-OR)", y="SNP effect on HF (log-OR)",
       title="CAD→HF per-SNP scatter: the three estimators coincide") +
  theme_ckm(legend="bottom")

fig3 <- (pa | pb) / (pc | pd) / (pe | pf) +
  plot_layout(heights=c(1, 1, 0.9)) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(fig3, "Figure3", 183, 232)
