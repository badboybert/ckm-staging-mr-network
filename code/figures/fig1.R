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
# Figure 1 — The European CKM causal network and its concordance with AHA-2023 staging.
# Panels: a DAG (stage-lane flow) · b canonical/control edges · c staging permutation null
#         d pre-registered falsifiable ledger (Component C) · e concordance robustness
source(file.path(P4_BASE, "figures", "fig_setup.R"))
suppressMessages(library(dplyr))
D <- file.path(FIG, "data")
edges <- read.csv(file.path(D, "fig1_edges.csv"), stringsAsFactors = FALSE)
stat  <- read.csv(file.path(D, "fig1_staging_stat.csv"), stringsAsFactors = FALSE)
getv  <- function(k) stat$value[stat$metric == k][1]
perm  <- read.csv(file.path(D, "fig1_permnull.csv"), stringsAsFactors = FALSE)
raw   <- rd("forward_local_edges.csv")
# NATIVE-SCALE ledger (scripts/57_ledger_native.py, reviewer C4/E3). Supersedes staging_ledger.csv,
# whose power-gate/TOST rule compared quantities measured on different scales. The legacy file is
# still written by 11_ledger.py and its statistics survive as secondary columns here.
ledg  <- rd("staging_ledger_native.csv")
# Read the permutation P from the same file that supplies the plotted null, then ASSERT the two
# agree. Hardcoding it here printed P = 1.8e-3 over a histogram whose own 10,000 draws yielded
# 1.1e-3, because prep_fig1.py permuted over a differently ordered node list. The assertion makes
# that class of drift impossible: the panel can only print what its null vector supports.
PERM_P     <- as.numeric(getv("perm_p"))

# Round-7 C089: one P formatter for every panel in this figure, matching the manuscript's convention
# (three decimals at or above 0.01, otherwise a two-significant-figure mantissa with a real
# superscript exponent), so a panel and the prose print the SAME string.
.sup <- function(e) {
  d <- c("⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")
  paste0("⁻", paste(d[as.integer(strsplit(as.character(abs(e)), "")[[1]]) + 1], collapse=""))
}
.pf1 <- function(p, dec = 3) {
  if (!is.finite(p)) return("NA")
  if (p >= 0.01) return(sprintf(paste0("%.", dec, "f"), p))
  e <- floor(log10(p)); m <- p / 10^e
  paste0(sub("\\.0$", "", sprintf("%.1f", m)), " × 10", .sup(e))
}
.pf <- function(p, dec = 3) vapply(p, .pf1, character(1), dec = dec)

.obs_exact <- as.numeric(getv("n_forward")) / as.numeric(getv("n_cross"))
# The permutation null is now enumerated EXACTLY (fig1_permnull.csv holds every distinct labeling's
# concordance), so the printed P is the exact tail fraction, not an add-one Monte-Carlo estimate. The
# panel can only print what its null vector supports: the assertion below recomputes the exact
# fraction from the plotted null and must equal it.
.recomp    <- sum(perm$concordance >= .obs_exact - 1e-9) / nrow(perm)
if (abs(.recomp - PERM_P) > 1e-9)   # both sides exact now; MC-era 5e-4 would hide a real 1.5-vs-1.1e-3 drift
  stop(sprintf("Figure 1b: stated perm P = %.5f but the plotted null yields %.5f", PERM_P, .recomp))

# ---------- 1a. staged causal-network DAG (cross-stage edges) ----------
node_xy <- tibble::tribble(
  ~node,   ~x, ~y,
  "BMI",    0,  0,
  "SBP",    1,  3.2, "HbA1c", 1, 1.7, "T2D", 1, 0, "LDL", 1, -1.4, "TC", 1, -2.8, "TG", 1, -4.2, "HDL", 1, 4.6, "eGFR", 1, -5.6,
  "CAD",    2,  0.3, "HF", 2, -2.2, "Stroke", 2, 2.6)
node_xy$stage <- STAGE[node_xy$node]
ce <- edges[edges$cross == "True", ]
ce <- merge(ce, node_xy[, c("node","x","y")], by.x="exp", by.y="node")
ce <- merge(ce, node_xy[, c("node","x","y")], by.x="out", by.y="node", suffixes=c("",".e"))
ce$dir <- ifelse(ce$direction == "backward", "backward", "forward")

# Round-7 C031: pull each edge back from the centre of its target node by the disc radius, in
# data units (the x lane spacing is 1 unit; a size-8.5 disc is about 0.17 x-units and 0.42
# y-units on this layout).
.rx <- 0.11; .ry <- 0.80
.dx <- ce$x.e - ce$x; .dy <- ce$y.e - ce$y
.len <- pmax(sqrt((.dx/.rx)^2 + (.dy/.ry)^2), 1e-9)
ce$x.t <- ce$x.e - .dx/.len
ce$y.t <- ce$y.e - .dy/.len

pa <- ggplot() +
  # stage lane guides
  annotate("rect", xmin=-0.45, xmax=0.45, ymin=-6.4, ymax=5.4, fill=STAGE_COL["1"], alpha=0.06) +
  annotate("rect", xmin= 0.55, xmax=1.45, ymin=-6.4, ymax=5.4, fill=STAGE_COL["2"], alpha=0.06) +
  annotate("rect", xmin= 1.55, xmax=2.45, ymin=-6.4, ymax=5.4, fill=STAGE_COL["4"], alpha=0.06) +
  annotate("text", x=c(0,1,2), y=5.9, label=c("Stage 1","Stage 2","Stage 4"), size=FS_T/.pt, fontface=2, colour="grey25") +
  # Round-7 C031: every arrowhead used to land under a node disc (the discs are drawn after the
  # edges and are 8.5 pt across), so a figure whose whole point is DIRECTION showed none. Each
  # edge now stops short of the target node - the shortening is computed in the data units of
  # the layout, so it holds at any canvas size - and the discs no longer cover the heads.
  geom_curve(data=ce[ce$dir=="forward",], aes(x=x, y=y, xend=x.t, yend=y.t),
             curvature=0.16, linewidth=0.26, colour="grey65", alpha=0.55,
             arrow=arrow(length=unit(3.4,"pt"), type="closed")) +
  geom_curve(data=ce[ce$dir=="backward",], aes(x=x, y=y, xend=x.t, yend=y.t),
             curvature=-0.28, linewidth=0.7, colour=POS_COL, linetype="21",
             arrow=arrow(length=unit(4.4,"pt"), type="closed")) +
  geom_point(data=node_xy, aes(x=x, y=y, fill=stage), shape=21, size=9.6, stroke=0.4, colour="white") +
  geom_text(data=node_xy, aes(x=x, y=y, label=node), size=FS/.pt-1.5, fontface=2, colour="black") +
  scale_fill_manual(values=STAGE_COL, guide="none") +
  annotate("label", x=-0.5, y=-3.6,
           label=sprintf("AHA CKM staging\nconcordance = %.3f\nperm. P = %s",
                         as.numeric(getv("observed_concordance")), .pf(PERM_P)),
           hjust=0, vjust=1, size=FS/.pt-0.7, colour="grey15", label.size=0, fill=NA, lineheight=0.95) +
  # 2026-07-26: "resolved in Fig 3" -> "evaluated in Figure 3". Figure 3 grades these edges; it does
  # not resolve them (HF->CAD is reported as unresolved), so the old label promised more than the
  # figure it points at delivers.
  annotate("text", x=1.0, y=-7.15, label="red dashed = reverse/feedback association (evaluated in Figure 3)", size=FS/.pt-1.2, colour=POS_COL, hjust=0.5) +
  coord_cartesian(xlim=c(-0.55,2.55), ylim=c(-7.4,6.4), clip="off") +
  theme_void(base_family=FONT) + theme(plot.margin=margin(4,4,2,4))

# ---------- 1b. canonical / control edges (validity) ----------
ctrl <- data.frame(edge=c("LDL->CAD","BMI->CAD","BMI->T2D","HbA1c->T2D","TG->CAD","HDL->CAD"),
                   exp=c("LDL","BMI","BMI","HbA1c","TG","HDL"), out=c("CAD","CAD","T2D","T2D","CAD","CAD"),
                   role=c("+ control","+ control","+ control","+ control","canonical","- control"))
ci <- do.call(rbind, lapply(seq_len(nrow(ctrl)), function(i){
  r <- raw[raw$exposure==ctrl$exp[i] & raw$outcome==ctrl$out[i], ]
  data.frame(edge=ctrl$edge[i], role=ctrl$role[i], b=r$ivw_b[1], lo=r$ivw_b[1]-1.96*r$ivw_se[1],
             hi=r$ivw_b[1]+1.96*r$ivw_se[1], p=r$ivw_p[1]) }))
# The ASCII "X->Y" key matches rows in raw; only the printed tick label is typeset (see disp_edge).
ci$edge <- factor(ci$edge, levels=rev(ctrl$edge), labels=rev(disp_edge(ctrl$edge)))
ci$sign <- ifelse(ci$b>0, "pos", "neg")
pb <- ggplot(ci, aes(b, edge, colour=sign)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0, linewidth=0.5, orientation="y") +
  geom_point(size=2.1) +
  geom_text(aes(label=sprintf("%.2f", b)), vjust=-0.9, size=FS/.pt-0.8, colour="black") +
  scale_colour_manual(values=c(pos=POS_COL, neg=NEG_COL), guide="none") +
  labs(x="Causal effect (per-SD / log-OR)", y=NULL,
       title="Canonical edges recovered (+/- controls)") +
  coord_cartesian(xlim=c(-0.6, 1.25)) + theme_ckm()

# ---------- 1c. staging permutation null ----------
obs <- as.numeric(getv("observed_concordance")); pp <- PERM_P   # derived + asserted at the top
pc <- ggplot(perm, aes(concordance)) +
  geom_histogram(bins=40, fill="grey80", colour="white", linewidth=0.15) +
  geom_vline(xintercept=obs, colour=POS_COL, linewidth=0.7) +
  annotate("text", x=obs-0.02, y=Inf, label=sprintf("observed = %.3f\nP = %s", obs, .pf(pp)), colour=POS_COL,
           hjust=1, vjust=1.4, size=FS/.pt-0.7, lineheight=0.95) +
  # vjust must clear the "observed" annotation above it (at vjust=1.4 the two collide mid-panel), and
  # the label must be LEFT-ANCHORED inside xlim — centred at x=0.30 it overran the 0.15 limit and the
  # leading "la" was clipped.
  # The null is the EXHAUSTIVE enumeration of every distinct stage labeling, not a Monte-Carlo sample.
  # The count is read off the plotted vector so this caption cannot drift from what the panel shows —
  # it previously read "(10,000 shuffles)" and survived the switch to exact enumeration untouched,
  # because a hardcoded caption is invisible to both mtime staleness and number-vs-source checks.
  annotate("text", x=0.17, y=Inf,
           label=sprintf("label-permutation null\n(all %s labellings)", format(nrow(perm), big.mark=",")),
           colour="grey35", hjust=0, vjust=3.4, size=FS/.pt-0.8, lineheight=0.95) +
  labs(x="Cross-stage concordance", y="Permutations", title="AHA staging-order concordance vs chance") +
  coord_cartesian(xlim=c(0.15,1.02)) + theme_ckm()

# ---------- 1d. pre-registered falsifiable ledger (Component C) ----------
vc <- c(CONCORDANT="#2E7D32", INDETERMINATE="grey60", DISCORDANT="#E69F00", DISCORDANT_CAVEATED=POS_COL)
# tick labels are typeset (X→Y); nothing downstream keys on ledg$transition after this point.
ledg$transition <- factor(ledg$transition, levels=rev(ledg$transition),
                          labels=rev(disp_edge(ledg$transition)))
ledg$bnum <- as.numeric(sub("\\+","",ledg$b_fwd))
ledg$verdict <- factor(ledg$verdict, levels=c("CONCORDANT","INDETERMINATE","DISCORDANT","DISCORDANT_CAVEATED"))
# counts are DERIVED from the data, never typed — the legend previously hardcoded "(11)" and would
# have silently kept saying 11 after the ledger rule was rebuilt.
.vn <- table(ledg$verdict)
pd <- ggplot(ledg, aes(bnum, transition, colour=verdict)) +
  geom_point(size=2) +
  geom_text(aes(label=sprintf("%.2f", bnum)), hjust=-0.35, size=FS/.pt-1.1, colour="grey30") +
  scale_colour_manual(values=vc, name=NULL,
     labels=c(CONCORDANT=sprintf("concordant (%d)", .vn[["CONCORDANT"]]),
              INDETERMINATE=sprintf("indeterminate (%d)", .vn[["INDETERMINATE"]]),
              DISCORDANT=sprintf("reverse effect (%d)", .vn[["DISCORDANT"]]),
              # 2026-09-17: at 8 pt this label overran the legend clip and the PAGE showed
              # "reverse, pleiotropy-caveated (2" -- the closing bracket was cut. Shortened to the
              # caption's own term ("reverse-caveated") rather than shrinking all four entries.
              DISCORDANT_CAVEATED=sprintf("reverse, caveated (%d)", .vn[["DISCORDANT_CAVEATED"]]))) +
  guides(colour=guide_legend(nrow=2, byrow=TRUE)) +
  labs(x="Forward-transition causal effect", y=NULL,
       # invariant #10: "pre-specified", never "pre-registered" (there is no time-stamped registration)
       title=sprintf("Falsifiable staging ledger (%d transitions)", nrow(ledg))) +
  coord_cartesian(xlim=c(0, 1.25)) + theme_ckm(legend="bottom") +
  theme(legend.text=element_text(size=8), axis.text.y=element_text(size=6.4))

# ---------- 1e. concordance robustness across thresholds ----------
# Both rows are DERIVED, never typed. The genome-wide-strict pair is parsed out of the analysis
# output that produced it (results/staging_forward.txt); the Bonferroni pair comes from the same
# fig1_staging_stat.csv that panel c plots. Hardcoding 1.000/0.926/0.0038 here meant the panel could
# not notice if either analysis moved.
.gw_txt  <- readLines(file.path(RES, "staging_forward.txt"), warn = FALSE)
.gw_conc <- as.numeric(sub(".*CONCORDANCE\\s*=\\s*([0-9.]+).*", "\\1",
                           grep("Stage-order CONCORDANCE", .gw_txt, value = TRUE)[1]))
.gw_p    <- as.numeric(sub(".*empirical p\\s*=\\s*([0-9.eE+-]+).*", "\\1",
                           grep("empirical p", .gw_txt, value = TRUE)[1]))
stopifnot(is.finite(.gw_conc), is.finite(.gw_p), .gw_conc >= 0, .gw_conc <= 1)
rob <- data.frame(thr=c("Genome-wide-\nstrict threshold","Bonferroni\n(bidirectional)"),
                  conc=c(.gw_conc, obs), p=c(.gw_p, pp), x=c(1,2))
pe <- ggplot(rob, aes(x, conc)) +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey65", linewidth=0.3) +
  annotate("text", x=1.5, y=0.515, label="chance (0.5)", size=FS/.pt-1, colour="grey55") +
  geom_point(size=2.6, colour=POS_COL) +
  geom_text(aes(label=sprintf("%.3f\nP=%s", conc, .pf(p))), vjust=-0.5, size=FS/.pt-1, lineheight=0.9) +
  scale_x_continuous(breaks=c(1,2), labels=rob$thr) +
  labs(x=NULL, y="Cross-stage concordance", title="Robust to edge-inclusion threshold") +
  # ylim upper must clear the TWO-LINE label drawn above a point at y = 1.000. At 1.08 the top line
  # ("1.000") was sliced by the panel edge in both the PNG and the PDF.
  coord_cartesian(xlim=c(0.6,2.4), ylim=c(0.45,1.20)) + theme_ckm() +
  theme(axis.text.x=element_text(size=8))

# ---------- 1f. constrained-null battery: how far the staging claim goes (E2) ----------
# THE CEILING PANEL. The concordance beats a label permutation (the stage labels are not
# interchangeable) but does NOT beat a degree-preserving edge-rewiring null: given that risk factors
# are instrumented as senders and diseases as receivers, a same-degree graph reaches comparable
# concordance ~1/3 of the time. Reported, not buried.
cn  <- read.csv(file.path(D, "fig1_constrained_nulls.csv"), stringsAsFactors=FALSE)
rob2 <- read.csv(file.path(D, "fig1_null_robustness.csv"), stringsAsFactors=FALSE)
bs  <- read.csv(file.path(D, "fig1_bootstrap.csv"), stringsAsFactors=FALSE)
gb  <- function(k) bs$value[bs$stat==k]
# The "Label permutation" row of the battery is the SAME statistic panel c plots. Script 47 runs it
# under an alphabetical node enumeration and script 07 under first-appearance order, so the two
# 10,000-draw Monte-Carlo estimates differ in the last digit (0.0011 vs 0.0012) and Figure 1 was
# printing both. Display the canonical value, and assert the battery's own estimate agrees with it
# to within Monte-Carlo error for a P of this size (~4e-4 at n = 10,000).
.lab_i <- which(cn$null == "Label permutation")
stopifnot(length(.lab_i) == 1)
if (abs(cn$p[.lab_i] - pp) > 5e-4)
  stop(sprintf("Figure 1f: battery label-permutation P = %.5f disagrees with the canonical %.5f by more than Monte-Carlo error",
               cn$p[.lab_i], pp))
cn$p[.lab_i] <- pp
# The DEGREE-PRESERVING row had no such assertion, and it is the row that drifted: this panel shipped
# "P = 0.31" from a 2026-07-24 CSV while the Abstract, Results, Discussion and this figure's own
# legend all said 0.27, the value the 2026-07-26 analysis produced. Guarding one of three rows is
# what let it through, so the ceiling statistic is now checked against its canonical output directly.
.deg_i <- which(cn$null == "Degree-preserving rewiring")
stopifnot(length(.deg_i) == 1)
.deg_canon <- as.numeric(sub(".*P = ([0-9.]+).*", "\\1",
  grep("^C\\. Degree-preserving edge-rewiring null", readLines(file.path(RES, "staging_constrained_nulls.txt"),
                                                               warn=FALSE), value=TRUE)[1]))
stopifnot(is.finite(.deg_canon))
if (abs(cn$p[.deg_i] - .deg_canon) > 1e-6)
  stop(sprintf("Figure 1f: panel degree-preserving P = %.4f but staging_constrained_nulls.txt says %.4f — re-run prep_newpanels.py",
               cn$p[.deg_i], .deg_canon))
cn$null <- factor(cn$null, levels=rev(cn$null))
cn$verdict <- factor(cn$verdict, levels=c("exceeded","NOT exceeded","uninformative"))
pf <- ggplot(cn, aes(p, null, colour=verdict)) +
  annotate("rect", xmin=0.05, xmax=Inf, ymin=-Inf, ymax=Inf, fill="grey94") +
  geom_vline(xintercept=0.05, linetype="dashed", linewidth=0.4, colour="grey45") +
  geom_segment(aes(x=1e-3, xend=p, yend=null), linewidth=0.5) +
  geom_point(size=2.8) +
  # hjust must flip with position: at this (inset, ~40%-width) panel a centred label overflows the
  # axis at both ends — the leftmost ran into the y-axis text, the rightmost off the panel.
  geom_text(aes(label=paste("P =", .pf(p, dec = 2)),   # the prose states the ceiling P to two decimals (0.27)
                hjust=ifelse(p < 0.01, -0.12, 1.12)),
            vjust=-1.2, size=FS/.pt-1.3, colour="black") +
  scale_colour_manual(values=c("exceeded"="#2E7D32","NOT exceeded"=POS_COL,"uninformative"=NULL_COL),
                      name=NULL) +
  # 0.05 is NOT a tick: at this width its label overprints the 0.1 label. The dashed line and its
  # "P = 0.05" annotation mark the threshold instead.
  scale_x_log10(breaks=c(0.001,0.01,0.1,1), labels=c("0.001","0.01","0.1","1")) +
  # two-line y labels: the single-line forms protrude into the plot gutter and make this row's
  # plotting area visibly narrower than the rows above it.
  scale_y_discrete(labels=c("Label permutation"="Label\npermutation",
                            "Degree-preserving rewiring"="Degree-preserving\nrewiring",
                            "Phenotype-class permutation"="Phenotype-class\npermutation")) +
  annotate("text", x=0.05, y=0.62, label="P = 0.05", size=FS/.pt-1.6, colour="grey45", hjust=-0.1, vjust=0) +
  # The robustness figures (LOO / drop-SBP+Stroke / bootstrap CI) are stated in the FIGURE LEGEND,
  # not on-panel: at this width the two-line annotation printed straight over the grey lollipop.
  # They are still derived here so the legend and the panel cannot drift apart:
  #   min LOO = rob2$concordance, drop SBP+Stroke, bootstrap CI = gb("ci_lo")/gb("ci_hi").
  coord_cartesian(ylim=c(0.5,3.9), clip="off") +
  guides(colour=guide_legend(nrow=1)) +
  labs(x="P vs each null (log scale)", y=NULL) +
  theme_ckm(legend="bottom") + theme(axis.text.y=element_text(size=6.6, lineheight=0.85))

# ---------- 1g. inserting the missing subclinical stage-3 tier (CAC) ----------
cac <- read.csv(file.path(D, "fig1_cac_staging.csv"), stringsAsFactors=FALSE)
cac$network <- factor(cac$network, levels=cac$network)
cac$lab <- ifelse(is.na(cac$p) | cac$p=="", sprintf("%.3f", cac$concordance),
                  sprintf("%.3f\nP=%s", cac$concordance, .pf(as.numeric(cac$p))))
pg <- ggplot(cac, aes(network, concordance)) +
  geom_col(width=0.6, fill="#8DA0CB") +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey60", linewidth=0.3) +
  geom_text(aes(label=lab), vjust=-0.25, size=FS/.pt-1.4, lineheight=0.9) +
  # Round 4 (T2-6) asked for this caveat to be enlarged or shortened. Raised from FS/.pt-1.6 to the
  # size the bar-value labels use, and re-checked in the PNG for a collision with the bars.
  annotate("text", x=2, y=0.30, size=FS/.pt-1.4, colour="grey30", lineheight=0.9,
           label="no backward CAC edge\namong the directions tested") +
  coord_cartesian(ylim=c(0,1.16)) +
  labs(x=NULL, y="Cross-stage concordance") +
  theme_ckm() + theme(axis.text.x=element_text(size=8))

# ---------- compose ----------
# panel order: a DAG · b permutation-null · c falsifiable ledger · d threshold-robustness ·
#              e control-strip · f constrained-null battery (the ceiling) · g CAC stage-3 tier
# The DAG needs a tall row: squeezing it to 1.02 to fit the 4th row crushed the Stage-2 node column
# and pushed the red-dashed caption onto the eGFR node. Give it back its height and grow the canvas.
# f|g are inset to ~80% width via a nested spacer row so the composite reads with an even margin.
# plot_spacer() does NOT consume a tag letter, so f and g keep their letters.
# Do NOT use a patchwork `design` grid for this: it forces ONE shared column grid across every row.
row_fg <- (plot_spacer() | pf | pg | plot_spacer()) + plot_layout(widths=c(0.25, 1, 1, 0.25))
fig1 <- (pa) / (pc | pd) / (pe | pb) / row_fg +
  plot_layout(heights=c(1.16, 1.06, 0.84, 0.74)) +
  plot_annotation(tag_levels="a", theme=theme(plot.tag=element_text(size=FS_TAG, face="bold")))
save_fig(fig1, "Figure1", 170, 224)
