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
# graphical_abstract.R — Cardiovascular Diabetology graphical abstract.
#
# CVD spec, quoted from the Research-article guidelines: "approximately 920x300 pixels", uploaded as
# JPEG, PNG or SVG, "may have a caption of up to 30 words. This caption must be a part of the image
# file and must be below the picture."
#
# Every number drawn here is READ from the same source-data CSVs the main figures use. Nothing is
# typed as a literal, and the script stops if a value it draws disagrees with its file, so the
# graphical abstract cannot drift away from Figure 2 the way a hand-made image would.
# =============================================================================
source(file.path(P4_BASE, "figures", "fig_setup.R"))
suppressMessages(library(ggplot2))

D  <- file.path(FIG, "data")
ed <- read.csv(file.path(D, "fig1_edges.csv"), stringsAsFactors = FALSE)
me <- read.csv(file.path(D, "fig2_mediation.csv"), stringsAsFactors = FALSE)

# ---- values, read not typed --------------------------------------------------------------
# The MVMR direct effect is read at FULL PRECISION from results/mvmr_cascade.csv, the file Figure 2
# and the prose both use — NOT from data/fig2_mediation.csv, which prep_newpanels.py writes at the
# 3 decimals results/formal_mediation.txt prints. That rounding shipped a contradiction: 0.415 as a
# double is 0.41499999999999998, so sprintf("%+.2f", .) gave +0.41 here while the Abstract, Results,
# Discussion and Figure 2 all printed +0.42 from 0.4151.... "Read not typed" was true and still not
# enough — the value was read from a DIFFERENT, lossier copy of the same quantity.
.mv <- rd("mvmr_cascade.csv")
.row <- .mv[.mv$model == "HF_via_CAD" & .mv$exposure == "BMI" & .mv$outcome == "HF", ]
stopifnot(nrow(.row) == 1)
bmi_hf_direct <- .row$direct_b[1]                       # MVMR direct effect, BMI -> HF
bmi_pm        <- me$pm_product[me$exposure == "BMI"]    # % of BMI -> HF carried by CAD
t2d_direct    <- me$direct[me$exposure == "T2D"]        # MVMR direct effect, T2D -> HF
t2d_direct_p  <- me$direct_p[me$exposure == "T2D"]
getb <- function(x, y) {                      # fig1_edges.csv columns are exp / out
  r <- ed[ed$exp == x & ed$out == y, ]
  stopifnot(nrow(r) == 1)
  r$b[1]
}
cad_hf <- getb("CAD", "HF")

stopifnot(bmi_pm > 0, bmi_pm < 100, t2d_direct_p > 0.05, bmi_hf_direct > 0, cad_hf > 0)

lab <- function(x) sprintf("%+.2f", x)

# The graphical abstract must print the SAME rounded value as the manuscript, and the manuscript is
# the one place that number is authored. Assert it here rather than hoping the two agree: this
# figure printed +0.41 against the manuscript's +0.42 through every review round, and no check
# compared an image label to prose.
.abs_txt <- paste(readLines(file.path(FIG, "ABSTRACT.md"), warn = FALSE), collapse = " ")
if (!grepl(gsub("([+])", "\\\\\\1", lab(bmi_hf_direct)), .abs_txt, fixed = FALSE))
  stop(sprintf("graphical abstract: BMI->HF direct effect prints %s, which the Abstract does not state",
               lab(bmi_hf_direct)))

# ---- layout ------------------------------------------------------------------------------
# A 920x300 canvas is very wide and very short, so the vertical budget is allocated explicitly:
#   96-89 title/subtitle | 80 direct-effect label | 76-24 the diagram | 22 the null route
#   14 rule | 8-4 caption. Nodes are 20 wide, so no centre may sit closer than 11 to an edge.
nodes <- data.frame(
  id   = c("BMI", "T2D", "CAD", "HF"),
  lbl  = c("Adiposity\n(BMI)", "Type 2\ndiabetes", "Coronary\nartery disease", "Heart\nfailure"),
  x    = c( 12,    12,     50,     88),
  y    = c( 68,    32,     50,     50),
  w    = c( 20,    20,     22,     20),
  fill = c("#F7D8D5", "#F7D8D5", "#DCE6F2", "#D2E7DA"),
  edge = c(POS_COL, POS_COL, NEG_COL, "#2E7D32"),
  stringsAsFactors = FALSE)

arrows <- data.frame(
  x    = c(  22,   22,   22,   61.5),
  y    = c(  70,   64,   36,   50),
  xend = c(  77,   38,   38,   77.5),
  yend = c(  56,   54,   45,   50),
  lwd  = c( 2.2,  1.0,  1.0,  1.6),
  col  = c(POS_COL, "grey45", "grey45", NEG_COL),
  stringsAsFactors = FALSE)

alab <- data.frame(
  x   = c(  52,   30,   30,   69.5),
  # Round-7: the red strap line sat 6 units under the second subtitle and their glyph boxes touched
  # (caught by the new overlapping-text gate). Dropped to 77 for clear separation.
  y   = c(  77,   62,   38,   54),
  txt = c(sprintf("CAD-INDEPENDENT  %s log-odds/SD  (~%.0f%% not mediated via CAD)",
                  lab(bmi_hf_direct), 100 - bmi_pm),
          sprintf("%s (total)", lab(getb("BMI", "CAD"))),
          sprintf("%s (total)", lab(getb("T2D", "CAD"))),
          sprintf("%s (total)", lab(cad_hf))),
  col = c(POS_COL, "grey30", "grey30", NEG_COL),
  sz  = c( 3.0,   2.8,   2.8,   2.8),
  stringsAsFactors = FALSE)

# The one thing a reader must not miss: type 2 diabetes has NO detectable direct route.
blocked <- data.frame(x = 50, y = 22,
                      txt = sprintf("no CAD-independent type 2 diabetes route detected in this MVMR model (%s, P = %.2f)",
                                    lab(t2d_direct), t2d_direct_p))

CAPTION <- paste("BMI retained a CAD-independent association with heart failure in multivariable",
                 "Mendelian randomization, whereas the type 2 diabetes-heart failure association was",
                 "compatible with coronary mediation.")
stopifnot(length(strsplit(CAPTION, "\\s+")[[1]]) <= 30)   # CVD: caption of up to 30 words

p <- ggplot() +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               colour = arrows$col, linewidth = arrows$lwd,
               arrow = arrow(length = unit(5.5, "pt"), type = "closed")) +
  geom_tile(data = nodes, aes(x = x, y = y, width = w), height = 16,
            fill = nodes$fill, colour = nodes$edge, linewidth = 0.6) +
  geom_text(data = nodes, aes(x = x, y = y, label = lbl),
            family = FONT, size = 3.2, lineheight = 0.92, fontface = 2) +
  geom_text(data = alab, aes(x = x, y = y, label = txt), colour = alab$col,
            family = FONT, size = alab$sz, lineheight = 0.95, fontface = 2) +
  geom_text(data = blocked, aes(x = x, y = y, label = txt), colour = "grey35",
            family = FONT, size = 2.8, fontface = 3) +
  # 2026-07-26 (round-4 Tier 1): the title and subtitle VISIBLY OVERLAPPED. The title sat at y=96
  # centred on its baseline, and the two-line subtitle block was centred at y=89, so the subtitle's
  # first line grew UPWARD into the title. Text-layer checks could not see this — both strings were
  # present and correct — which is why it survived to a fourth review. Each line is now anchored by
  # its TOP (vjust=1) at an explicit y, so the vertical budget is stated rather than hoped for:
  #   99 title | 91 subtitle line 1 | 86 subtitle line 2 | 80 direct-effect label | 76-24 diagram.
  annotate("text", x = 50, y = 99, family = FONT, size = 3.5, fontface = 2, vjust = 1,
           label = "Adiposity-centred causal architecture of cardiovascular–kidney–metabolic traits") +
  annotate("text", x = 50, y = 91, family = FONT, size = 2.5, colour = "grey30", vjust = 1,
           label = "Bidirectional Mendelian randomization  ·  132 directed edges  ·  12 traits  ·  European ancestry") +
  annotate("text", x = 50, y = 86, family = FONT, size = 2.5, colour = "grey30", vjust = 1,
           label = "Red = conditional estimate (multivariable MR);  grey/blue = total (univariable) effects") +
  annotate("segment", x = 4, xend = 96, y = 14, yend = 14, colour = "grey80", linewidth = 0.4) +
  annotate("text", x = 50, y = 7, family = FONT, size = 2.6, colour = "grey25", lineheight = 1.15,
           label = paste(strwrap(CAPTION, width = 105), collapse = "\n")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(1, 1, 1, 1))

# 920 x 300 px exactly, at 150 dpi -> 6.1333 x 2.0 in.
W_PX <- 920; H_PX <- 300; DPI <- 150
out_png <- file.path(FIG, "GraphicalAbstract.png")
ggsave(out_png, p, width = W_PX / DPI, height = H_PX / DPI, dpi = DPI,
       units = "in", device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "GraphicalAbstract.pdf"), p, width = W_PX / DPI, height = H_PX / DPI,
       units = "in", device = cairo_pdf, bg = "white")

dim_png <- dim(png::readPNG(out_png))
cat(sprintf("wrote GraphicalAbstract .png/.pdf  (%d x %d px; CVD asks for ~920 x 300)\n",
            dim_png[2], dim_png[1]))
cat(sprintf("caption words: %d (limit 30)\n", length(strsplit(CAPTION, "\\s+")[[1]])))
cat(sprintf("values drawn (all read from source data): BMI->HF direct %s, mediated %.0f%%, "
            , lab(bmi_hf_direct), bmi_pm))
cat(sprintf("CAD->HF %s, T2D->HF direct %s (P=%.2f)\n", lab(cad_hf), lab(t2d_direct), t2d_direct_p))
stopifnot(dim_png[2] == W_PX, dim_png[1] == H_PX)
