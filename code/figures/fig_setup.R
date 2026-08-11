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
# =============================================================================
# fig_setup.R — shared theme, palettes, helpers, and data loaders for Paper 4 figures.
# Target band: Cardiovascular Diabetology / JAHA / Circ Genom Precis Med.
# Output: 300-dpi PNG (ragg, native Arial on Windows) + vector PDF. source() this from each figXX.R.
# =============================================================================
suppressMessages({library(ggplot2); library(patchwork); library(scales); library(grid)})
BASE <- P4_BASE
RES  <- file.path(BASE, "results")
FIG  <- file.path(BASE, "figures")
PAN  <- file.path(FIG, "panels")

# Font sizes bumped +1 body / +2-3 annotations (user, 2026-07-15: intra-figure labels were too small).
# FS/FS_T anchor the geom_text/annotation sizes (size=FS/.pt-X); they were pinned to the OLD 7 pt base and
# rendered at ~4-6 pt. Now FS=10 -> annotations land at ~7-9 pt. Point/node sizes are literals (unaffected).
FONT <- "Arial"; FS <- 10; FS_T <- 11; FS_TAG <- 12   # annotation anchor 10, small-annotation 11, panel label 12 bold
FLAT_BASE <- 10   # uniform body font (axis text / titles / legend / strip = 10 pt)

# ---- minimal / Tufte L-axis style (fallback; former default) ----
theme_ckm_minimal <- function(base = FS, legend = "right") {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = "black"),
      plot.title = element_text(size = FS_T, face = "plain", hjust = 0, margin = margin(b = 3)),
      axis.title = element_text(size = FS_T, colour = "black"),
      axis.text = element_text(size = base, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      legend.position = legend,
      legend.text = element_text(size = base, colour = "black"),
      legend.title = element_text(size = base, colour = "black"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(8, "pt"), legend.margin = margin(1, 1, 1, 1),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(size = FS_T, face = "plain", margin = margin(2, 2, 2, 2)),
      plot.margin = margin(6, 6, 4, 4), plot.tag = element_text(size = FS_TAG, face = "bold")
    )
}

# ---- FLAT / framed style (theme_bw base, uniform font, panel border + major gridlines) ----
# plot.title is BLANKED: per-panel takeaways live only in the figure legends.
theme_ckm_flat <- function(base = FLAT_BASE, legend = "right") {
  theme_bw(base_size = base, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = "black"),
      plot.title = element_blank(), plot.subtitle = element_blank(),   # titles -> legends
      axis.title = element_text(size = base, colour = "black"),
      axis.title.x = element_text(margin = margin(t = 2)),
      axis.title.y = element_text(margin = margin(r = 2)),
      axis.text = element_text(size = base, colour = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      strip.text = element_text(size = base, face = "plain", margin = margin(2, 2, 2, 2)),
      legend.position = legend,
      legend.text = element_text(size = base, colour = "black"),
      legend.title = element_text(size = base, colour = "black"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(9, "pt"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.margin = margin(1, 1, 1, 1),
      plot.margin = margin(11, 6, 4, 4),   # extra top so the bold a/b/c tag clears the framed panel
      plot.tag = element_text(size = FS_TAG, face = "bold")
    )
}

# ---- PUBLICATION / Nature-minimal style (theme_classic base, clean L-axis, NO border, NO gridlines) ----
# Mirrors the reversibility-paper theme (paper 3.eas theme_publication_ckm): body 9 pt everywhere, panel
# label 11 pt bold, white strips (no box), grey30 subtitle. In-panel titles blanked -> live in the caption.
# This is the ACTIVE style (2026-07-15, user request to match the reversibility figures). Heatmap/tile panels
# blank their own axis lines per-panel (already done in fig3/suppfig1), so the L-axis default is safe.
theme_ckm_pub <- function(base = FLAT_BASE, legend = "right") {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = "black"),
      plot.title = element_blank(),                              # titles -> caption (clean panels)
      plot.subtitle = element_text(size = base - 1, colour = "grey30", hjust = 0, margin = margin(b = 4)),
      axis.title = element_text(size = base, colour = "black"),
      axis.title.x = element_text(margin = margin(t = 3)),
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.text = element_text(size = base, colour = "black"),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),                              # NO gridlines (Nature-minimal)
      panel.border = element_blank(),                            # NO frame
      strip.background = element_rect(fill = "white", colour = NA),   # clean strip, no box
      strip.text = element_text(size = base, face = "plain", margin = margin(2, 2, 2, 2)),
      legend.position = legend,
      legend.text = element_text(size = base, colour = "black"),
      legend.title = element_text(size = base, colour = "black"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(9, "pt"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.margin = margin(1, 1, 1, 1),
      plot.margin = margin(10, 6, 4, 4),
      plot.tag = element_text(size = FS_TAG, face = "bold")
    )
}

theme_ckm <- theme_ckm_pub   # active style (Nature-minimal; theme_ckm_flat kept as the framed fallback)

# ---- palettes (color-blind-conscious; redundant shape/pattern used alongside) ----
# AHA CKM stage colour (nodes): stage 1 adiposity, stage 2 metabolic/kidney, stage 4 clinical CVD
STAGE_COL <- c("1" = "#E69F00", "2" = "#56B4E9", "4" = "#D55E00")
STAGE_LAB <- c("1" = "Stage 1 (adiposity)", "2" = "Stage 2 (metabolic / kidney)", "4" = "Stage 4 (clinical CVD)")
# trait -> stage
STAGE <- c(BMI="1", SBP="2",TG="2",HDL="2",TC="2",LDL="2",HbA1c="2",FI="2",T2D="2",eGFR="2",CKD="2",
           ApoB="2", CAD="4",HF="4",Stroke="4",AF="4")
# cohort colours
COHORT_COL <- c(EUR = "#0072B2", BBJ = "#009E73", TPMI = "#D55E00", KoGES = "#CC79A7", TWB = "#E69F00",
                meta = "#333333", Hospital = "#009E73", Population = "#CC79A7")
# verdict / direction colours
POS_COL <- "#B2182B"; NEG_COL <- "#2166AC"; NULL_COL <- "grey60"
CAUSAL_COL <- "#B2182B"; SHARING_COL <- "#4393C3"

# ---- reader-facing typography ----
# Edge keys arrive from results/*.csv as ASCII "X->Y" and are used BOTH to match rows and to label
# panels. Matching must stay on the raw key; only the DISPLAYED string is typeset. Round 4 found
# "LDL->CAD" and "beta_iv" printed on shipped panels — internal notation reaching the reader. Keep ONE
# implementation here: the same defect was fixed independently in four supplementary renderers in an
# earlier pass and came back on the two main figures nobody had looked at.
# Gate F11 asserts no ASCII arrow survives in any shipped figure's PDF text layer.
disp_edge <- function(x) gsub("->", "→", as.character(x), fixed = TRUE)   # X->Y  =>  X→Y

# significance -> shape (filled = sig, open = ns), consistent across figures
sig_shape <- function(p, thr = 0.05) ifelse(p < thr, 16, 1)
fmt_p <- function(p) ifelse(p < 0.01, sprintf("%.0e", p), sprintf("%.2f", p))
fmt_b <- function(b, lo, hi) sprintf("%.2f (%.2f, %.2f)", b, lo, hi)

# ---- savers ----
save_panel <- function(p, name, w, h) {  # w,h in mm ; writes PNG(300) + PDF to panels/
  ragg::agg_png(file.path(PAN, paste0(name, ".png")), width = w, height = h, units = "mm", res = 300, background = "white")
  print(p); invisible(dev.off())
}
save_fig <- function(p, name, w, h) {     # composite -> figures/ PNG(400) + PDF
  ragg::agg_png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "mm", res = 400, background = "white")
  print(p); invisible(dev.off())
  grDevices::cairo_pdf(file.path(FIG, paste0(name, ".pdf")), width = w/25.4, height = h/25.4)
  print(p); invisible(dev.off())
  cat("wrote", name, ".png/.pdf  (", w, "x", h, "mm )\n")
}

# ---- data loaders (tidy frames from results/) ----
rd <- function(f) read.csv(file.path(RES, f), stringsAsFactors = FALSE, check.names = FALSE)
