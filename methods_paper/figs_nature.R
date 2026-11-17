# =====================================================================
# Nature-style reproduction figures -- panel-for-panel composites
# Falcomata et al., Nature 658, 1330-1340 (2026)
#
# Every panel of the source paper's Fig. 1-5 and Extended Data Fig. 1-10
# gets a slot, lettered as in the paper.
#
#   * panels whose numbers are in the deposited Source Data (MOESM9-17) are
#     redrawn from those numbers;
#   * panels that can be recomputed from GSE289138 are redrawn from the
#     per-library tables exported by audit_mouse_confound.py;
#   * panels that are schematics, micrographs, or analyses of external
#     cohorts get an explicit grey placeholder stating why they cannot be
#     reproduced from public deposits.
#
# No p-value is invented.  Where a p-value is shown it is either recomputed
# from the deposited per-mouse values (Welch t-test, as in the paper's
# legend) and labelled with the paper's claimed value, or taken from the
# paper and quoted as "claim".
#
# Output: methods_paper/figs_nature/*.pdf
# Run:    Rscript methods_paper/figs_nature.R
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(gtable)
})

ROOT   <- "E:/pancreatic_cancer/submit"
SD_DIR <- file.path(ROOT, "serpin_repro", "data", "Nature_SourceData")
OUT    <- file.path(ROOT, "methods_paper", "figs_nature")
PT     <- file.path(ROOT, "methods_paper", "figs_reproduction", "plot_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

COL_CTRL <- "#4d4d4d"
COL_B    <- "#1f77b4"
COL_E    <- "#d62728"
COL_MIX  <- "#f2a541"

# geometry of a composite panel: 3 columns on a 7.2 in canvas
BS <- 6.6   # theme base size
AX <- 5.2   # axis text
TI <- 7.8   # panel letter
ROW_H <- 1.62

SD_FILES <- list(
  fig2 = "41586_2026_11002_MOESM9_ESM.xlsx",
  fig3 = "41586_2026_11002_MOESM10_ESM.xlsx",
  fig4 = "41586_2026_11002_MOESM11_ESM.xlsx",
  fig5 = "41586_2026_11002_MOESM12_ESM.xlsx",
  ed2  = "41586_2026_11002_MOESM13_ESM.xlsx",
  ed4  = "41586_2026_11002_MOESM14_ESM.xlsx",
  ed5  = "41586_2026_11002_MOESM15_ESM.xlsx",
  ed7  = "41586_2026_11002_MOESM16_ESM.xlsx",
  ed8  = "41586_2026_11002_MOESM17_ESM.xlsx"
)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------
# Source Data readers
# ---------------------------------------------------------------------

# Read a bar-panel sheet. Handles both layouts found in the deposit:
#   (a) one header row  -> colnames are the group names, data from row 2
#   (b) two header rows -> row 1 = super-group, row 2 = subgroup
sd_bar <- function(fn_key, sheet) {
  raw <- as.data.frame(readxl::read_excel(
    file.path(SD_DIR, SD_FILES[[fn_key]]), sheet = sheet, col_names = FALSE))
  clean <- function(x) gsub("\u00a0", " ", trimws(as.character(x)))

  r2  <- clean(unlist(raw[2, ]))
  num <- suppressWarnings(as.numeric(r2))
  hdr <- if (any(!is.na(r2) & is.na(num))) 2L else 1L

  top <- clean(unlist(raw[1, ]))
  sub <- if (hdr == 2L) r2 else rep(NA_character_, ncol(raw))
  if (hdr == 2L) {
    for (j in seq_along(top)) {
      if ((is.na(top[j]) || top[j] %in% c("", "NA")) && j > 1)
        top[j] <- top[j - 1]
    }
  }
  body <- raw[(hdr + 1L):nrow(raw), , drop = FALSE]

  groups <- list()
  for (j in seq_len(ncol(raw))) {
    t <- top[j]; s <- sub[j]
    if (hdr == 2L && (is.na(s) || s %in% c("", "NA"))) next
    nm <- if (hdr == 2L) {
      if (!is.na(t) && nzchar(t) && t != "NA" && t != s) paste0(s, " (", t, ")") else s
    } else t
    if (is.na(nm) || nm %in% c("", "NA")) next
    v <- suppressWarnings(as.numeric(as.character(body[[j]])))
    v <- v[!is.na(v)]
    if (length(v) > 0) groups[[nm]] <- v
  }
  groups
}

# Read a survival sheet: column 1 = day, one flag column per arm
sd_surv <- function(fn_key, sheet) {
  raw <- readxl::read_excel(file.path(SD_DIR, SD_FILES[[fn_key]]), sheet = sheet)
  days <- raw[[1]]
  out <- list()
  for (g in colnames(raw)[-1]) {
    g2 <- gsub("\u00a0", " ", trimws(g))
    flag <- suppressWarnings(as.numeric(raw[[g]]))
    ok <- !is.na(days) & !is.na(flag)
    out[[g2]] <- days[ok][flag[ok] == 1]
  }
  out
}

# n per arm: one row per mouse, flag column is 1 for every row of that arm
sd_surv_n <- function(fn_key, sheet) {
  raw <- readxl::read_excel(file.path(SD_DIR, SD_FILES[[fn_key]]), sheet = sheet)
  out <- numeric(0)
  for (g in colnames(raw)[-1]) {
    g2 <- gsub("\u00a0", " ", trimws(g))
    flag <- suppressWarnings(as.numeric(raw[[g]]))
    out[g2] <- sum(!is.na(flag))
  }
  out
}

welch_p <- function(a, b) t.test(a, b, var.equal = FALSE)$p.value

# Lighten a colour towards white (second level of a repeated genotype)
lighten <- function(col, amt = 0.45) {
  v <- grDevices::col2rgb(col) / 255
  v <- v + (1 - v) * amt
  grDevices::rgb(v[1], v[2], v[3])
}

# Genotype-driven palette; repeated genotypes get a lightened second level
pal_for <- function(nms) {
  base <- function(g) {
    if (grepl("Control|Vehicle|IgG/PBS|^F8", g)) COL_CTRL
    else if (grepl("Serpinb2", g))              COL_B
    else if (grepl("Serpine1", g))              COL_E
    else                                        COL_MIX
  }
  cols <- vapply(nms, base, character(1))
  seen <- character(0)
  for (i in seq_along(cols)) {
    if (cols[i] %in% seen) cols[i] <- lighten(cols[i])
    seen <- c(seen, cols[i])
  }
  stats::setNames(unname(cols), nms)
}

# Short labels for group names inside p-value captions
short_lab <- function(g) {
  g <- sub("^95% Se1 KO \\+ 5% Control$", "mix", g)
  g <- gsub("Anti-CSF1/Clodronate-Lipo", "aCSF1", g)
  g <- gsub("IgG/PBS-Lipo", "IgG", g)
  g <- gsub("Serpinb2 KO", "b2 KO", g)
  g <- gsub("Serpine1 KO", "e1 KO", g)
  g <- gsub("Control", "Ctrl", g)
  g <- gsub("Vehicle control", "Veh", g)
  g <- gsub(" +\\+ +", "+", g)
  g <- gsub("\\(([^)]*)\\)", "[\\1]", g)
  g
}

wrap_txt <- function(x, width = 36) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------------
# Panel builders
# ---------------------------------------------------------------------

# Grey placeholder for a panel that cannot be reproduced from public data.
ph <- function(letter, reason, width = 34) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "#f4f4f4", colour = "#8c8c8c",
             linetype = "dashed", linewidth = 0.28) +
    annotate("text", x = 0.025, y = 0.965, label = letter, hjust = 0, vjust = 1,
             fontface = "bold", size = 2.9) +
    annotate("text", x = 0.5, y = 0.50, label = wrap_txt(reason, width),
             size = 1.85, colour = "#5a5a5a", lineheight = 0.95) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5))
}

# Nature-style bar + dot panel.  `letter` is the panel label from the paper.
nature_bar <- function(groups, ylab, letter, pairs = NULL, claimed = NULL,
                       logy = FALSE, colors = NULL, point_size = 1.2,
                       footnote = NULL) {
  n <- names(groups)
  if (is.null(colors))
    colors <- stats::setNames(
      c(COL_CTRL, COL_B, COL_E, COL_MIX, "#a6a6a6", "#9ecae1", "#f2a3a3")[seq_along(n)], n)

  mn  <- data.frame(mean = sapply(groups, mean),
                    group = factor(n, levels = n))
  pts <- do.call(rbind, lapply(n, function(g)
    data.frame(value = groups[[g]], group = g)))
  pts$group <- factor(pts$group, levels = n)

  p <- ggplot(mn, aes(group, mean, fill = group)) +
    geom_col(width = 0.58, colour = "black", linewidth = 0.28) +
    scale_fill_manual(values = colors, guide = "none") +
    geom_point(data = pts, aes(group, value), inherit.aes = FALSE,
               shape = 21, fill = "white", colour = "black", stroke = 0.28,
               size = point_size,
               position = position_jitter(width = 0.085, height = 0, seed = 1)) +
    theme_minimal(base_size = BS) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#ededed", linewidth = 0.18),
      axis.line = element_line(linewidth = 0.24),
      axis.ticks = element_line(linewidth = 0.24),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = AX - 0.3, angle = 32, hjust = 1,
                                 vjust = 1, lineheight = 0.92),
      axis.text.y = element_text(size = AX),
      axis.title.y = element_text(size = AX + 0.1),
      plot.title = element_text(size = TI, face = "bold", hjust = 0),
      plot.title.position = "plot",
      plot.margin = margin(2, 3, 2, 2)
    ) +
    labs(y = ylab, title = letter)

  if (logy) p <- p + scale_y_log10(expand = expansion(mult = c(0.03, 0.12)))

  cap_lines <- character(0)
  if (!is.null(pairs)) {
    labs_out <- character(0)
    for (i in seq_along(pairs)) {
      a <- pairs[[i]][1]; b <- pairs[[i]][2]
      if (is.null(groups[[a]]) || is.null(groups[[b]])) next
      pv  <- welch_p(groups[[a]], groups[[b]])
      lab <- sprintf("%s vs %s: p = %.4g", short_lab(b), short_lab(a), pv)
      if (!is.null(claimed) && length(claimed) >= i &&
          !is.na(claimed[[i]]) && nzchar(claimed[[i]]))
        lab <- paste0(lab, " (claim ", claimed[[i]], ")")
      labs_out <- c(labs_out, lab)
    }
    if (length(labs_out) > 0)
      cap_lines <- c(cap_lines,
                     strwrap(paste(labs_out, collapse = "; "), width = 54))
  }
  if (!is.null(footnote)) cap_lines <- c(cap_lines, strwrap(footnote, width = 64))
  if (length(cap_lines) > 0) {
    p <- p + labs(caption = paste(cap_lines, collapse = "\n")) +
      theme(plot.caption = element_text(size = 4.2, hjust = 0, face = "italic",
                                        lineheight = 0.95, colour = "black"),
            plot.caption.position = "plot")
  }
  p
}

# Kaplan-Meier panel: one axes, one curve per arm (as in the paper)
km_panel <- function(surv, n, letter, cols, ltypes = NULL, note = NULL,
                     ylab = "Survival (%)") {
  km <- do.call(rbind, lapply(names(surv), function(gn) {
    t  <- sort(surv[[gn]])
    nt <- n[[gn]]
    data.frame(x = c(0, t), y = c(100, 100 * (1 - seq_along(t) / nt)),
               label = gn)
  }))
  km$label <- factor(km$label, levels = names(surv))

  p <- ggplot(km, aes(x, y, colour = label, linetype = label, group = label)) +
    geom_step(linewidth = 0.5) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_linetype_manual(
      values = ltypes %||% stats::setNames(rep("solid", length(surv)), names(surv)),
      name = NULL) +
    scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 101),
                       expand = expansion(mult = c(0.01, 0.03))) +
    theme_minimal(base_size = BS) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#ededed", linewidth = 0.18),
      axis.line = element_line(linewidth = 0.24),
      axis.ticks = element_line(linewidth = 0.24),
      axis.text = element_text(size = AX - 0.2),
      axis.title = element_text(size = AX + 0.1),
      legend.position = "bottom",
      legend.key.size = unit(0.26, "cm"),
      legend.text = element_text(size = 4.5),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-2, 0, 0, 0),
      plot.title = element_text(size = TI, face = "bold", hjust = 0),
      plot.title.position = "panel",
      plot.margin = margin(6, 3, 2, 2)
    ) +
    guides(colour = guide_legend(nrow = 3, byrow = TRUE),
           linetype = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(x = "Days", y = ylab, title = letter)

  if (!is.null(note))
    p <- p + labs(caption = wrap_txt(note, 52)) +
      theme(plot.caption = element_text(size = 4.2, hjust = 0, face = "italic",
                                        lineheight = 0.95),
            plot.caption.position = "plot")
  p
}

# Grouped bars: one bar per genotype within each subset
GENO_LAB <- c(control_KO = "Control", Serpinb2_KO = "Serpinb2 KO",
              Serpine1_KO = "Serpine1 KO")
GENO_COL <- c(control_KO = COL_CTRL, Serpinb2_KO = COL_B, Serpine1_KO = COL_E)

subset_bars <- function(df, subsets, letter, ylab = "% of cells",
                        legend = TRUE, legend_rows = 1) {
  long <- do.call(rbind, lapply(subsets, function(s)
    data.frame(subset = s, genotype = df$genotype, value = df[[s]])))
  long$subset   <- factor(long$subset, levels = subsets)
  long$genotype <- factor(long$genotype, levels = names(GENO_LAB))

  p <- ggplot(long, aes(subset, value, fill = genotype)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72,
             colour = "black", linewidth = 0.24) +
    scale_fill_manual(values = GENO_COL, labels = GENO_LAB, name = NULL) +
    theme_minimal(base_size = BS) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#ededed", linewidth = 0.18),
          axis.line = element_line(linewidth = 0.24),
          axis.ticks = element_line(linewidth = 0.24),
          axis.text.x = element_text(size = AX - 0.4, angle = 32, hjust = 1,
                                     vjust = 1, lineheight = 0.92),
          axis.text.y = element_text(size = AX),
          axis.title.y = element_text(size = AX + 0.1),
          axis.title.x = element_blank(),
          legend.position = if (legend) "bottom" else "none",
          legend.key.size = unit(0.24, "cm"),
          legend.text = element_text(size = 4.8),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-3, 0, 0, 0),
          plot.title = element_text(size = TI, face = "bold", hjust = 0),
          plot.title.position = "plot",
          plot.margin = margin(2, 3, 2, 2)) +
    guides(fill = guide_legend(nrow = legend_rows)) +
    labs(x = NULL, y = ylab, title = letter)
  p
}

# UMAP panel for an arbitrary coordinate table
umap_panel <- function(df, col, letter, pal = NULL, legend_rows = 2,
                       pt = 0.12, leg_cex = 4.4) {
  p <- ggplot(df, aes(UMAP1, UMAP2, colour = .data[[col]])) +
    geom_point(size = pt, stroke = 0) +
    theme_void(base_size = BS) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.18, "cm"),
          legend.text = element_text(size = leg_cex),
          legend.title = element_blank(),
          legend.spacing.x = unit(0.04, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-4, 0, 0, 0),
          plot.title = element_text(size = TI, face = "bold", hjust = 0),
          plot.title.position = "plot",
          plot.margin = margin(1, 2, 1, 1)) +
    guides(colour = guide_legend(nrow = legend_rows,
                                 override.aes = list(size = 1.1))) +
    labs(title = letter)
  if (!is.null(pal)) p <- p + scale_colour_manual(values = pal)
  p
}

# ---------------------------------------------------------------------
# Composite assembly
# ---------------------------------------------------------------------

grid_layout <- function(n, ncol = 3) {
  rows <- ceiling(n / ncol)
  m <- matrix(seq_len(rows * ncol), nrow = rows, ncol = ncol, byrow = TRUE)
  m[m > n] <- NA_integer_
  m
}

compose <- function(plots, file, ncol = 3, height = NULL, width = 7.2) {
  stopifnot(length(plots) > 0)
  lay <- grid_layout(length(plots), ncol)
  if (is.null(height)) height <- nrow(lay) * ROW_H + 0.15
  pdf(file.path(OUT, file), width = width, height = height, useDingbats = FALSE)
  grid.arrange(grobs = plots, layout_matrix = lay)
  dev.off()
  message("wrote ", file, "  (", length(plots), " panels, ",
          nrow(lay), " x ", ncol, ")")
}

have <- function(f) file.exists(file.path(PT, f))
rd   <- function(f) read.delim(file.path(PT, f), check.names = FALSE,
                               stringsAsFactors = FALSE)

umap_like <- function(df, col, letter, pal = NULL, legend_rows = 2,
                      leg_cex = 4.4, pt = 0.12) {
  # drop unused factor levels so the legend does not show empty entries
  df[[col]] <- factor(df[[col]])
  umap_panel(df, col, letter, pal = pal, legend_rows = legend_rows,
             leg_cex = leg_cex, pt = pt)
}

# =====================================================================
# Fig. 1 -- Perturb-map screen (no Source Data deposited for any panel)
# =====================================================================
message("=== Fig. 1 ===")

fig1 <- list(
  ph("a", "Perturb-map workflow schematic -- drawing, no underlying data", 24),
  ph("b", "Representative MICSSS images, days 7 / 14 / 21 -- raw microscopy not deposited", 24),
  ph("c", "Spatial heterogeneity maps (non-matching neighbours) -- per-cell coordinates not deposited", 24),
  ph("d", "Pro-Code composition of focal regions -- per-cell coordinates not deposited", 24),
  ph("e", "Fitness heat map, in vivo / in vitro log2FC -- per-gene values in no Source Data file", 24),
  ph("f", "Neighbourhood enrichment z-scores (Squidpy) -- output not deposited", 24),
  ph("g", "Rag2-/- experimental design schematic -- drawing, no underlying data", 24),
  ph("h", "Knockout abundance across Rag2-/- tumours -- per-mouse counts not deposited", 24)
)
compose(fig1, "fig1.pdf", ncol = 4, height = 2 * ROW_H + 0.15)

# =====================================================================
# Fig. 2 -- loss of Serpinb2 / Serpine1 promotes T-cell immunity
# =====================================================================
message("=== Fig. 2 ===")

g2d <- sd_bar("fig2", "Figure_2d")
g2e <- sd_bar("fig2", "Figure_2e")
g2h <- sd_bar("fig2", "Figure_2h")
g2i <- sd_bar("fig2", "Figure_2i")
g2k <- sd_bar("fig2", "Figure_2k")
s2l <- sd_surv("fig2", "Figure_2l")
s2n <- sd_surv("fig2", "Figure_2n")

n2l <- c(`Control IgG` = 5, `Serpinb2 KO IgG` = 5, `Serpine1 KO IgG` = 5,
         `Control aPD1` = 6, `Serpinb2 KO aPD1` = 6, `Serpine1 KO aPD1` = 7)
n2n <- c(`Vehicle control` = 6, `PAI-039` = 6, `aPD1` = 6, `PAI-039+aPD1` = 7)

km2l_cols <- c(`Control IgG` = COL_CTRL, `Serpinb2 KO IgG` = COL_B,
               `Serpine1 KO IgG` = COL_E, `Control aPD1` = lighten(COL_CTRL, 0.55),
               `Serpinb2 KO aPD1` = lighten(COL_B, 0.55),
               `Serpine1 KO aPD1` = lighten(COL_E, 0.55))
km2l_lty <- c(`Control IgG` = "solid", `Serpinb2 KO IgG` = "solid",
              `Serpine1 KO IgG` = "solid", `Control aPD1` = "dashed",
              `Serpinb2 KO aPD1` = "dashed", `Serpine1 KO aPD1` = "dashed")
km2n_cols <- c(`Vehicle control` = COL_CTRL, `PAI-039` = COL_B,
               `aPD1` = COL_MIX, `PAI-039+aPD1` = COL_E)

p2f <- if (have("tcell_subsets.tsv")) {
  tc <- rd("tcell_subsets.tsv"); tc <- tc[tc$chemistry == "classic", ]
  subset_bars(tc, c("CD8 Tex term", "CD8 Tex prog", "CD8 Tex int", "CD8 Teff",
                    "CD8 Tn-like", "CD4 Teff Gzmk+", "CD4 Tn-like", "CD4 Treg"),
              "f", ylab = "T cells (%)")
} else ph("f", "T-cell subsets -- tcell_subsets.tsv not exported")

fig2 <- list(
  ph("a", "TCGA/GTEx differential expression -- raw matrix not deposited (UCSC Xena)"),
  ph("b", "Forest plot, 8 PDAC cohorts (n = 1,176) -- per-study HRs not deposited"),
  ph("c", "In vivo validation design schematic -- drawing"),
  nature_bar(g2d, "Tumour mass (g)", "d",
             pairs = list(c(names(g2d)[1], names(g2d)[2]),
                          c(names(g2d)[1], names(g2d)[3])),
             claimed = c("0.001", "0.0003")),
  nature_bar(g2e, "CD8+ T cells / mg", "e",
             pairs = list(c(names(g2e)[1], names(g2e)[2]),
                          c(names(g2e)[1], names(g2e)[3])),
             claimed = c("0.0436", "")),
  p2f,
  ph("g", "Representative IHC (CK19 / CD8a / GZMB) -- raw images not deposited"),
  nature_bar(g2h, "CD8+ cells / mm2", "h",
             pairs = list(c(names(g2h)[1], names(g2h)[2]),
                          c(names(g2h)[1], names(g2h)[3])),
             claimed = c("0.0064", "0.0064")),
  nature_bar(g2i, "GZMB+CD8+ cells / mm2", "i",
             pairs = list(c(names(g2i)[1], names(g2i)[2]),
                          c(names(g2i)[1], names(g2i)[3])),
             claimed = c("0.0018", "0.041")),
  ph("j", "Rag2-/- and NSG design schematic -- drawing"),
  nature_bar(g2k, "Tumour mass (g), NSG", "k",
             pairs = list(c(names(g2k)[1], names(g2k)[2]),
                          c(names(g2k)[1], names(g2k)[3])),
             claimed = c("0.0043", "0.0275")),
  km_panel(s2l, n2l, "l", km2l_cols, km2l_lty,
           note = "log-rank claim: 4.7e-4 (IgG) / 1.8e-4 (anti-PD-1); solid = IgG, dashed = anti-PD-1"),
  ph("m", "PAI-039 +/- anti-PD-1 design schematic -- drawing"),
  km_panel(s2n, n2n, "n", km2n_cols,
           note = "log-rank claims 3.0e-4 / 2.4e-3 / 1.4e-2 (vehicle vs combo / PAI-039 vs combo / vehicle vs PAI-039)")
)
compose(fig2, "fig2.pdf", height = 5 * 1.24 + 0.35)

# =====================================================================
# Fig. 3 -- transcriptional programmes and ECM deposition
# =====================================================================
message("=== Fig. 3 ===")

g3e <- sd_bar("fig3", "Figure_3e")
g3f <- sd_bar("fig3", "Figure_3f")
g3j <- sd_bar("fig3", "Figure_3j")
g3l <- sd_bar("fig3", "Figure_3l")

p3g <- if (have("caf_umap.tsv")) {
  cu <- rd("caf_umap.tsv")
  umap_like(cu, "subset", "g", legend_rows = 2, leg_cex = 4.2, pt = 0.25)
} else ph("g", "CAF sub-cluster UMAP -- caf_umap.tsv not exported")

p3h <- if (have("caf_subsets.tsv")) {
  cf <- rd("caf_subsets.tsv"); cf <- cf[cf$chemistry == "Flex", ]
  subset_bars(cf, c("Matrix CAF", "apCAF", "MyCAF", "LRCC15+ MyCAF", "iCAF"),
              "h", ylab = "CAFs (%)", legend_rows = 1)
} else ph("h", "CAF subsets -- caf_subsets.tsv not exported")

fig3 <- list(
  ph("a", "Perturb-map Multi-modal workflow -- schematic"),
  ph("b", "Representative Xenium + t-CyCIF region -- raw images not deposited"),
  ph("c", "Pathway activity heat map (Xenium) -- on Zenodo, not in Source Data"),
  ph("d", "Volcano plots, Xenium DGE -- on Zenodo, not in Source Data"),
  nature_bar(g3e, "PDPN+ (%)", "e",
             pairs = list(c(names(g3e)[1], names(g3e)[2]),
                          c(names(g3e)[1], names(g3e)[3])),
             claimed = c("0.2218", "0.002")),
  nature_bar(g3f, "aSMA+PDPN+ (%)", "f",
             pairs = list(c(names(g3f)[1], names(g3f)[2]),
                          c(names(g3f)[1], names(g3f)[3])),
             claimed = c("0.7267", "0.0359")),
  p3g, p3h,
  ph("i", "H&E and Masson trichrome images -- raw images not deposited"),
  nature_bar(g3j, "Trichrome+ (%)", "j",
             pairs = list(c(names(g3j)[1], names(g3j)[2]),
                          c(names(g3j)[1], names(g3j)[3])),
             claimed = c("0.0564", "0.0282")),
  ph("k", "H&E and fibrin(ogen) IHC images -- raw images not deposited"),
  nature_bar(g3l, "Fibrin(ogen)+ (%)", "l",
             pairs = list(c(names(g3l)[1], names(g3l)[2]),
                          c(names(g3l)[1], names(g3l)[3])),
             claimed = c("0.0482", "0.0145"))
)
compose(fig3, "fig3.pdf")

# =====================================================================
# Fig. 4 -- fibrin-dependent macrophage retention and polarisation
# =====================================================================
message("=== Fig. 4 ===")

g4a <- sd_bar("fig4", "Figure_4a")
g4d <- sd_bar("fig4", "Figure_4d")
g4e <- sd_bar("fig4", "Figure_4e")
g4f <- sd_bar("fig4", "Figure_4f")
g4j <- sd_bar("fig4", "Figure_4j")

message("  4e groups: ", paste(names(g4e), collapse = " | "))
message("  4f groups: ", paste(names(g4f), collapse = " | "))
message("  4j groups: ", paste(names(g4j), collapse = " | "))

pal6 <- c(COL_CTRL, COL_B, COL_E, "#a6a6a6", "#9ecae1", "#f2a3a3")

p4g <- if (have("myeloid_subsets.tsv")) {
  my <- rd("myeloid_subsets.tsv"); my <- my[my$chemistry == "classic", ]
  subset_bars(my, c("Spp1/Marco+", "Cxcl9+", "Arg1+", "Proliferating",
                    "C1q+ homeostatic", "Folr2+", "Il1b+", "DC-like",
                    "Monocyte-like"), "g", ylab = "Myeloid (%)")
} else ph("g", "Myeloid subsets -- myeloid_subsets.tsv not exported")

fig4 <- list(
  nature_bar(g4a, "Macrophages / mg", "a",
             pairs = list(c(names(g4a)[1], names(g4a)[2]),
                          c(names(g4a)[1], names(g4a)[3])),
             claimed = c("", "0.0219")),
  ph("b", "Human PDAC MICSSS (CK19/CD11b/CD68/fibrinogen/PAI1/PAI2) -- raw images not deposited"),
  ph("c", "Spatial proximity density plots -- per-cell coordinates not deposited"),
  nature_bar(g4d, "Distance (\u00b5m)", "d",
             pairs = list(c(names(g4d)[1], names(g4d)[2]),
                          c(names(g4d)[1], names(g4d)[3])),
             claimed = c("", "")),
  nature_bar(g4e, "Tumour mass (g)", "e",
             pairs = list(c(names(g4e)[1], names(g4e)[4]),
                          c(names(g4e)[3], names(g4e)[6])),
             colors = stats::setNames(pal6[seq_along(g4e)], names(g4e)),
             footnote = "deposited p 0.0021; contrast not stated"),
  nature_bar(g4f, "Tumour mass (g)", "f",
             pairs = list(c(names(g4f)[2], names(g4f)[5]),
                          c(names(g4f)[3], names(g4f)[6])),
             colors = stats::setNames(pal6[seq_along(g4f)], names(g4f)),
             footnote = "deposited p 0.0014; contrast not stated"),
  p4g,
  ph("h", "Per-cell Arg1 / fibrinogen-binding scores -- per-cell expression not exported"),
  ph("i", "BMDM fibrin culture schematic -- drawing"),
  nature_bar(g4j, "MFI", "j",
             colors = stats::setNames(pal6[seq_along(g4j)], names(g4j))),
  ph("k", "Human spatial enrichment effect sizes (3 cohorts) -- external Visium not deposited"),
  ph("l", "Human spatial enrichment heat map -- external Visium not deposited")
)
compose(fig4, "fig4.pdf")

# =====================================================================
# Fig. 5 -- rare serpin+ cells dominate local immunity
# =====================================================================
message("=== Fig. 5 ===")

g5e <- sd_bar("fig5", "Figure_5e")
g5f <- sd_bar("fig5", "Figure_5f")
g5g <- sd_bar("fig5", "Figure_5g")
g5h <- sd_bar("fig5", "Figure_5h")

n5 <- names(g5e)
cols5 <- stats::setNames(c(COL_CTRL, COL_E, COL_MIX)[seq_along(n5)], n5)

# 5k: ROI-level CD8 fractions in three ROI categories (raw values, not means).
# readxl drops the empty leading column and the empty leading row, so the
# category header is row 1 and the ROI values start at row 2.
raw5k <- as.data.frame(readxl::read_excel(
  file.path(SD_DIR, SD_FILES$fig5), sheet = "Figure_5k", col_names = FALSE))
cat5k <- gsub("\u00a0", " ", trimws(as.character(unlist(raw5k[1, ]))))
d5k <- do.call(rbind, lapply(seq_along(cat5k), function(j) {
  v <- suppressWarnings(as.numeric(as.character(unlist(raw5k[2:nrow(raw5k), j]))))
  data.frame(cat = cat5k[j], val = v)
}))
d5k <- d5k[!is.na(d5k$val), ]
d5k$cat <- factor(d5k$cat, levels = cat5k)

p5k <- ggplot(d5k, aes(cat, val)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#ededed",
               colour = "black", linewidth = 0.28) +
  geom_point(aes(fill = cat), shape = 21, colour = "black", stroke = 0.22,
             size = 1.0, position = position_jitter(width = 0.13, height = 0,
                                                    seed = 1)) +
  scale_fill_manual(values = stats::setNames(
    c(COL_CTRL, COL_E, COL_MIX)[seq_along(cat5k)], cat5k), guide = "none") +
  theme_minimal(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#ededed", linewidth = 0.18),
        axis.line = element_line(linewidth = 0.24),
        axis.ticks = element_line(linewidth = 0.24),
        axis.text.x = element_text(size = AX - 0.4, angle = 32, hjust = 1,
                                   vjust = 1, lineheight = 0.92),
        axis.text.y = element_text(size = AX),
        axis.title = element_text(size = AX + 0.4),
        axis.title.x = element_blank(),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 3, 2, 2)) +
  labs(y = "CD8+ fraction (ROI)", title = "k")

# 5l: ROI-level control fraction vs CD8 abundance
raw5l <- as.data.frame(readxl::read_excel(
  file.path(SD_DIR, SD_FILES$fig5), sheet = "Figure_5l", col_names = FALSE))
d5l <- data.frame(
  f8  = suppressWarnings(as.numeric(as.character(unlist(raw5l[2:nrow(raw5l), 1])))),
  cd8 = suppressWarnings(as.numeric(as.character(unlist(raw5l[2:nrow(raw5l), 2])))))
d5l <- d5l[stats::complete.cases(d5l), ]

p5l <- ggplot(d5l, aes(f8, cd8)) +
  geom_point(shape = 21, fill = "#9ecae1", colour = "black", stroke = 0.25,
             size = 1.1) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.4,
              colour = COL_E, fill = "#f2a3a3", alpha = 0.35) +
  theme_minimal(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#ededed", linewidth = 0.18),
        axis.line = element_line(linewidth = 0.24),
        axis.ticks = element_line(linewidth = 0.24),
        axis.text = element_text(size = AX - 0.2),
        axis.title = element_text(size = AX + 0.2),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 3, 2, 2)) +
  labs(x = "Control (F8) tumour-cell fraction", y = "CD8+ T-cell fraction",
       title = "l")

fig5 <- list(
  ph("a", "SERPINB2+/SERPINE1+ fraction per patient (n = 71) -- human compendium not deposited"),
  ph("b", "Venn of shared upregulated genes -- DEG tables not deposited"),
  ph("c", "Pathway enrichment of serpin+ malignant cells -- not deposited"),
  ph("d", "Clonal mixing design schematic -- drawing"),
  nature_bar(g5e, "Tumour mass (g)", "e",
             pairs = list(c(n5[2], n5[3])), claimed = c("0.0578"),
             colors = cols5),
  nature_bar(g5f, "Macrophages (%)", "f",
             pairs = list(c(n5[1], n5[3])), claimed = c("0.0005"),
             colors = cols5),
  nature_bar(g5g, "CD8+ T cells (%)", "g",
             pairs = list(c(n5[1], n5[3])), claimed = c("0.9146"),
             colors = cols5),
  nature_bar(g5h, "GZMB+CD8+ (%)", "h",
             pairs = list(c(n5[1], n5[3])), claimed = c("0.0064"),
             colors = cols5),
  ph("i", "Injected vs recovered fraction -- per-mouse recovery not in Source Data"),
  ph("j", "Representative multiplex IF images -- raw images not deposited"),
  p5k, p5l
)
compose(fig5, "fig5.pdf")

# =====================================================================
# Extended Data Fig. 1 -- gene prioritisation
# =====================================================================
message("=== ED 1 ===")
ed1 <- list(
  ph("a", "Prioritisation workflow -- schematic"),
  ph("b", "DepMap dependency distribution (46 cell lines) -- per-gene scores not deposited"),
  ph("c", "Reactome enrichment of the candidate pool -- results only as gene lists"),
  ph("d", "Pfam / domain enrichment -- results only as gene lists"),
  ph("e", "Circular heat map -- requires the human and mouse scRNA-seq inputs")
)
compose(ed1, "ed1.pdf")

# =====================================================================
# Extended Data Fig. 2 -- Pro-Code barcoding library
# =====================================================================
message("=== ED 2 ===")

# 2e: relative in vivo/in vitro frequency of each PC, split by epitope tag
e2 <- as.data.frame(readxl::read_excel(
  file.path(SD_DIR, SD_FILES$ed2), sheet = "Ext_Data_Figure_2e", col_names = FALSE))
pc <- data.frame(name  = as.character(unlist(e2[, 1])),
                 value = suppressWarnings(as.numeric(as.character(unlist(e2[, 2])))))
pc <- pc[stats::complete.cases(pc), ]
tags2 <- c("AU1", "C", "HA", "HSV", "NWS", "S", "V5")
d2e <- do.call(rbind, lapply(tags2, function(t) {
  has <- vapply(strsplit(pc$name, ";"), function(z) t %in% z, logical(1))
  data.frame(tag = t, grp = ifelse(has, "with tag", "without"),
             value = pc$value)
}))
d2e$tag <- factor(d2e$tag, levels = tags2)
d2e$grp <- factor(d2e$grp, levels = c("with tag", "without"))
p2e <- ggplot(d2e, aes(tag, value, fill = grp)) +
  geom_point(shape = 21, colour = "black", stroke = 0.22, size = 1.0,
             position = position_jitterdodge(jitter.width = 0.08,
                                             dodge.width = 0.62, seed = 1)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, linewidth = 0.25,
               position = position_dodge(width = 0.62), show.legend = FALSE) +
  scale_fill_manual(values = c(`with tag` = COL_E, `without` = "#b3b3b3"),
                    name = NULL) +
  theme_minimal(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#ededed", linewidth = 0.18),
        axis.line = element_line(linewidth = 0.24),
        axis.ticks = element_line(linewidth = 0.24),
        axis.text = element_text(size = AX - 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = AX + 0.2),
        legend.position = "bottom",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 4.8),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-3, 0, 0, 0),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 3, 2, 2)) +
  labs(y = "In vivo/in vitro", title = "e",
       caption = wrap_txt(paste0(
         "Welch p per tag: ",
         paste(vapply(tags2, function(t) {
           s <- d2e[d2e$tag == t, ]
           sprintf("%s %.3g", t, welch_p(s$value[s$grp == "with tag"],
                                         s$value[s$grp == "without"]))
         }, character(1)), collapse = ", ")), 54)) +
  theme(plot.caption = element_text(size = 4.6, hjust = 0, face = "italic",
                                    lineheight = 0.95),
        plot.caption.position = "plot")

ed2 <- list(
  ph("a", "Pro-Code vector schematic -- drawing"),
  ph("b", "Pro-Code validation workflow -- schematic"),
  ph("c", "Representative section, H&E + single-marker stains -- raw images not deposited"),
  ph("d", "PC frequencies in vitro vs in vivo (ring plot) -- per-PC counts not deposited"),
  p2e,
  ph("f", "PC-PC neighbourhood enrichment heat map -- Squidpy output not deposited")
)
compose(ed2, "ed2.pdf")

# =====================================================================
# Extended Data Fig. 3 -- in vivo phenotypes vs in vitro fitness
# =====================================================================
message("=== ED 3 ===")
ed3 <- list(
  ph("a", "Perturb-map workflow and vector design -- schematic"),
  ph("b", "PC abundance, Cas9+ vs Cas9- after 2 weeks in culture -- counts not deposited"),
  ph("c", "DepMap dependency scores (46 cell lines) -- not deposited"),
  ph("d", "Image analysis / Pro-Code deconvolution workflow -- schematic")
)
compose(ed3, "ed3.pdf", ncol = 2)

# =====================================================================
# Extended Data Fig. 4 -- clonal and immune-context evolution
# =====================================================================
message("=== ED 4 ===")

e4 <- as.data.frame(readxl::read_excel(
  file.path(SD_DIR, SD_FILES$ed4), sheet = "Ext_Data_Figure_4c", col_names = FALSE))
# row 1 = metric name, row 2 = day label, data from row 3; column 4 is a spacer
days <- gsub("\u00a0", " ", trimws(as.character(unlist(e4[2, 1:3]))))
body <- e4[3:nrow(e4), , drop = FALSE]
sh <- body[, 1:3, drop = FALSE]; names(sh) <- days
ev <- body[, 5:7, drop = FALSE]; names(ev) <- days
to_long <- function(m, metric) {
  out <- do.call(rbind, lapply(names(m), function(d)
    data.frame(day = d, value = suppressWarnings(as.numeric(as.character(m[[d]]))),
               metric = metric)))
  out[stats::complete.cases(out), ]
}
d4c <- rbind(to_long(sh, "Shannon index"), to_long(ev, "Evenness index"))
d4c$day <- factor(d4c$day, levels = days)
d4c$metric <- factor(d4c$metric, levels = c("Shannon index", "Evenness index"))

p4c <- ggplot(d4c, aes(day, value)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#ededed",
               colour = "black", linewidth = 0.28) +
  geom_point(shape = 21, fill = COL_B, colour = "black", stroke = 0.22,
             size = 0.9, position = position_jitter(width = 0.11, height = 0,
                                                    seed = 1)) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  theme_minimal(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#ededed", linewidth = 0.18),
        strip.text = element_text(size = AX + 0.2, face = "bold"),
        axis.line = element_line(linewidth = 0.24),
        axis.ticks = element_line(linewidth = 0.24),
        axis.text = element_text(size = AX - 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 3, 2, 2)) +
  labs(title = "c", caption = "n = 5 (d7) / 7 (d14) / 8 (d21) mice") +
  theme(plot.caption = element_text(size = 4.0, hjust = 0, face = "italic"),
        plot.caption.position = "plot")

ed4 <- list(
  ph("a", "Digital reconstruction at day 14 -- per-cell coordinates not deposited"),
  ph("b", "PC frequencies within focal areas -- not deposited"),
  p4c,
  ph("d", "Radial neighbourhood-enrichment plots -- Squidpy output not deposited"),
  ph("e", "Longitudinal CD8 infiltration per knockout -- not deposited"),
  ph("f", "Rag2-/- fitness volcano -- not deposited"),
  ph("g", "Immunocompetent vs Rag2-/- volcano -- not deposited")
)
compose(ed4, "ed4.pdf")

# =====================================================================
# Extended Data Fig. 5 -- functional validation and immunotherapy
# =====================================================================
message("=== ED 5 ===")

g5c <- sd_bar("ed5", "Ext_Data_Figure_5c")
g5d <- sd_bar("ed5", "Ext_Data_Figure_5d")
g5j <- sd_bar("ed5", "Ext_Data_Figure_5j")
s5e <- sd_surv("ed5", "Ext_Data_Figure_5e")
n5e <- sd_surv_n("ed5", "Ext_Data_Figure_5e")
s5l <- sd_surv("ed5", "Ext_Data_Figure_5l")
n5l <- sd_surv_n("ed5", "Ext_Data_Figure_5l")

p5f <- if (have("umap_coords.tsv")) {
  um <- rd("umap_coords.tsv")
  set.seed(1)
  if (nrow(um) > 12000) um <- um[sample(nrow(um), 12000), ]
  anno_lv  <- sort(unique(um$anno))
  anno_pal <- stats::setNames(scales::hue_pal()(length(anno_lv)), anno_lv)
  umap_panel(um, "anno", "f", pal = anno_pal, legend_rows = 4, leg_cex = 3.9)
} else ph("f", "UMAP of major populations -- umap_coords.tsv not exported")

p5g <- if (have("tcell_umap.tsv")) {
  tu <- rd("tcell_umap.tsv")
  tcols <- stats::setNames(
    c("#7f7f7f", "#bcbd22", "#17becf", "#d62728", "#ff7f0e", "#2ca02c",
      "#9467bd", "#8c564b", "#1f77b4")[seq_len(length(unique(tu$subset)))],
    sort(unique(tu$subset)))
  umap_panel(tu, "subset", "g", pal = tcols, legend_rows = 3, leg_cex = 4.0, pt = 0.25)
} else ph("g", "T-cell subset UMAP -- tcell_umap.tsv not exported")

p5i <- if (have("tcell_subsets.tsv")) {
  tc <- rd("tcell_subsets.tsv"); tc <- tc[tc$chemistry == "classic", ]
  subset_bars(tc, c("CD8 Tex prog", "CD8 Tex int", "CD8 Tex term"),
              "i", ylab = "T cells (%)")
} else ph("i", "Terminally exhausted subsets -- tcell_subsets.tsv not exported")

km5e_cols <- stats::setNames(
  c(COL_CTRL, COL_B, COL_E)[seq_along(s5e)], names(s5e))
km5l_cols <- stats::setNames(
  c(COL_CTRL, COL_B, COL_MIX, COL_E)[seq_along(s5l)], names(s5l))

ed5 <- list(
  ph("a", "TCGA hazard-ratio lollipop -- per-cancer-type HRs not deposited"),
  ph("b", "Western blot (PAI1 / PAI2, two sgRNAs) -- raw blots in Supplementary Data 3"),
  nature_bar(g5c, "Doubling time", "c", colors = pal_for(names(g5c))),
  nature_bar(g5d, "Tumour mass (g)", "d",
             pairs = list(c(names(g5d)[1], names(g5d)[2])), claimed = c(""),
             colors = pal_for(names(g5d))),
  km_panel(s5e, n5e, "e", km5e_cols,
           note = "log-rank p quoted in the paper; n = 5 per group"),
  p5f, p5g,
  ph("h", "Hallmark gene heat map across T-cell subsets -- per-cell z-scores not exported"),
  p5i,
  nature_bar(g5j, "Tumour mass (g)", "j",
             colors = pal_for(names(g5j))),
  ph("k", "Serpin knockout + anti-PD-1 design schematic -- drawing"),
  km_panel(s5l, n5l, "l", km5l_cols,
           note = "KPC model; vehicle n = 6, PAI-039 n = 7, anti-PD-1 n = 6, combination n = 6"),
  ph("m", "JAVELIN Renal 101 progression-free survival -- trial data not deposited")
)
compose(ed5, "ed5.pdf", height = 5 * 1.24 + 0.35)

# =====================================================================
# Extended Data Fig. 6 -- Perturb-map Multi-modal
# =====================================================================
message("=== ED 6 ===")
ed6 <- list(
  ph("a", "Multi-modal workflow -- schematic"),
  ph("b", "Xenium -> t-CyCIF registration -- schematic"),
  ph("c", "Coagulation gene heat map in malignant cells -- not deposited"),
  ph("d", "Xenium UMAP by cell type -- raw Xenium object on Zenodo"),
  ph("e", "Xenium UMAP by knockout -- raw Xenium object on Zenodo"),
  ph("f", "Tumour-cell UMAP by knockout -- raw Xenium object on Zenodo")
)
compose(ed6, "ed6.pdf")

# =====================================================================
# Extended Data Fig. 7 -- stromal and endothelial composition
# =====================================================================
message("=== ED 7 ===")

g7b <- sd_bar("ed7", "Ext_Data_Figure_7b")
g7e <- sd_bar("ed7", "Ext_Data_Figure_7e")

p7c <- if (have("umap_coords.tsv")) {
  um <- rd("umap_coords.tsv"); um <- um[um$condition == "Flex", ]
  anno_lv  <- sort(unique(um$anno))
  anno_pal <- stats::setNames(scales::hue_pal()(length(anno_lv)), anno_lv)
  umap_panel(um, "anno", "c", pal = anno_pal, legend_rows = 4, leg_cex = 3.8)
} else ph("c", "Flex-library UMAP -- umap_coords.tsv not exported")

p7f <- if (have("endo_umap.tsv")) {
  eu <- rd("endo_umap.tsv")
  ecols <- stats::setNames(
    c(COL_CTRL, COL_B, COL_E, COL_MIX)[seq_along(unique(eu$subset))],
    sort(unique(eu$subset)))
  umap_panel(eu, "subset", "f", pal = ecols, legend_rows = 2, leg_cex = 4.2, pt = 0.25)
} else ph("f", "Endothelial subset UMAP -- endo_umap.tsv not exported")

p7g <- if (have("endo_subsets.tsv")) {
  en <- rd("endo_subsets.tsv")
  subs <- setdiff(names(en), c("filename", "genotype", "chemistry", "n"))
  subset_bars(en, subs, "g", ylab = "Endothelial (%)")
} else ph("g", "Endothelial subsets -- endo_subsets.tsv not exported")

ed7 <- list(
  ph("a", "t-CyCIF CAF images (PDPN / aSMA / PDGFRa) -- raw images not deposited"),
  nature_bar(g7b, "PDPN+/aSMA+ (%)", "b",
             colors = pal_for(names(g7b))),
  p7c,
  ph("d", "t-CyCIF CD31 images -- raw images not deposited"),
  nature_bar(g7e, "CD31+ cells (%)", "e", colors = pal_for(names(g7e))),
  p7f, p7g
)
compose(ed7, "ed7.pdf")

# =====================================================================
# Extended Data Fig. 8 -- PAI1/PAI2 and suppressive macrophages
# =====================================================================
message("=== ED 8 ===")

ed8_sheets <- c(a = "Ext_Data_Figure_8a",  h = "Ext_Data_Figure_8h",
                j = "Ext_Data_Figure_8j",  k = "Ext_Data_Figure_8k",
                l = "Ext_Data_Figure_8l",  m = "Ext_Data_Figure_8m",
                o = "Ext_Data_Figure_8o",  r = "Ext_Data_Figure_8r")
ed8_ylab <- c(a = "Neutrophils/mg", h = "Macrophages/mg",
              j = "GZMB+CD8+ (%)", k = "CD8+ (%)",
              l = "Macrophages (%)", m = "ARG1+ (%)",
              o = "Tumour mass (g)", r = "MFI")

ed8_bar <- function(letter) {
  g <- sd_bar("ed8", ed8_sheets[[letter]])
  cols <- if (letter == "r") {
    stats::setNames(c("#4d4d4d", "#1f77b4", "#f2a541",
                      "#b3b3b3", "#9ecae1", "#f2a3a3")[seq_along(g)], names(g))
  } else pal_for(names(g))
  nature_bar(g, ed8_ylab[[letter]], letter, colors = cols)
}

p8p <- if (have("myeloid_umap.tsv")) {
  mu <- rd("myeloid_umap.tsv")
  mcols <- stats::setNames(
    c("#7f7f7f", COL_B, COL_E, "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b",
      "#17becf", "#bcbd22")[seq_along(unique(mu$subset))],
    sort(unique(mu$subset)))
  umap_panel(mu, "subset", "p", pal = mcols, legend_rows = 3, leg_cex = 3.9, pt = 0.25)
} else ph("p", "Macrophage subset UMAP -- myeloid_umap.tsv not exported")

ed8 <- list(
  ed8_bar("a"),
  ph("b", "Mouse MICSSS (CK19/CD11b/F4/80/fibrinogen/PAI1/PAI2) -- raw images not deposited"),
  ph("c", "Tumour centre / periphery reconstruction -- per-cell coordinates not deposited"),
  ph("d", "Human spatial analysis workflow -- schematic"),
  ph("e", "Distance density plots, mouse and human -- per-cell distances not deposited"),
  ph("f", "Human macrophage distances (n = 30) -- not deposited"),
  ph("g", "Macrophage-depletion design schematic -- drawing"),
  ed8_bar("h"),
  ph("i", "CD18-blockade design schematic -- drawing"),
  ed8_bar("j"), ed8_bar("k"), ed8_bar("l"), ed8_bar("m"),
  ph("n", "Anti-PD-1 + anti-CD18 design schematic -- drawing"),
  ed8_bar("o"),
  p8p,
  ph("q", "GO enrichment in knockout macrophages -- requires re-running DE + Enrichr"),
  ed8_bar("r")
)
compose(ed8, "ed8.pdf", height = 6 * 1.20 + 0.30)

# =====================================================================
# Extended Data Fig. 9 / 10 -- human spatial and single-cell cohorts
# =====================================================================
message("=== ED 9 / ED 10 ===")
ed9 <- list(
  ph("a", "SERPINB2 compartment fractions (57 Visium sections) -- external data not deposited"),
  ph("b", "SERPINE1 compartment fractions (57 Visium sections) -- external data not deposited")
)
compose(ed9, "ed9.pdf", ncol = 2, height = ROW_H * 1.6)

ed10 <- list(
  ph("a", "Mean expression in serpin+ cells per patient (n = 71) -- compendium not deposited"),
  ph("b", "Fraction of serpin+ and double+ cells per patient -- not deposited"),
  ph("c", "Patient-level composition (stacked bars) -- not deposited"),
  ph("d", "Fraction by PDAC subtype (violin) -- not deposited"),
  ph("e", "Fraction across molecular subtypes -- not deposited"),
  ph("f", "Relative composition within subtype -- not deposited"),
  ph("g", "Fraction by treatment status -- not deposited"),
  ph("h", "Mean expression by treatment status -- not deposited"),
  ph("i", "Proposed model schematic -- drawing")
)
compose(ed10, "ed10.pdf")

# =====================================================================
# Supplementary panels -- our own readouts, not panels of the paper
# =====================================================================
message("=== Supplementary ===")

# S1: the reconstructed object
if (have("umap_coords.tsv")) {
  um <- rd("umap_coords.tsv")
  set.seed(1)
  if (nrow(um) > 10000) um <- um[sample(nrow(um), 10000), ]
  um$genotype <- factor(GENO_LAB[um$genotype], levels = GENO_LAB)
  um$condition <- factor(ifelse(um$condition == "classic", "classic (3')", "Flex"),
                         levels = c("classic (3')", "Flex"))
  anno_lv  <- sort(unique(um$anno))
  anno_pal <- stats::setNames(scales::hue_pal()(length(anno_lv)), anno_lv)
  pdf(file.path(OUT, "figS1_umap.pdf"), width = 7.2, height = 3.1,
      useDingbats = FALSE)
  grid.arrange(grobs = list(
    umap_panel(um, "anno", "Cell type", pal = anno_pal, legend_rows = 3,
               leg_cex = 5.0),
    umap_panel(um, "genotype", "Genotype", legend_rows = 1, leg_cex = 5.0),
    umap_panel(um, "condition", "Chemistry", legend_rows = 1, leg_cex = 5.0)),
    nrow = 1, ncol = 3)
  dev.off()
  message("wrote figS1_umap.pdf")
}

# S2: serpin detection in malignant epithelium (mosaic knockdown check)
if (have("serpin_detection.tsv")) {
  sdet <- rd("serpin_detection.tsv")
  sdet$chemistry <- factor(sdet$chemistry, levels = c("classic", "Flex"))
  sdet <- sdet[order(sdet$chemistry, match(sdet$genotype, names(GENO_LAB))), ]
  sdet$label <- paste0(GENO_LAB[sdet$genotype], "\n", sdet$chemistry)
  sdet$label <- factor(sdet$label, levels = sdet$label)

  serpin_panel <- function(col, letter, ylab) {
    ggplot(sdet, aes(label, .data[[col]], fill = genotype)) +
      geom_col(width = 0.66, colour = "black", linewidth = 0.25) +
      geom_text(aes(label = sprintf("%.2f", .data[[col]])),
                vjust = -0.35, size = 1.9) +
      scale_fill_manual(values = GENO_COL, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      theme_minimal(base_size = 8.2) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "#ededed", linewidth = 0.2),
            axis.text.x = element_text(size = 6.0, angle = 35, hjust = 1),
            axis.text.y = element_text(size = 6.4),
            axis.title.y = element_text(size = 7.0),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 8.0, face = "bold", hjust = 0),
            plot.title.position = "plot",
            plot.margin = margin(2, 3, 2, 2)) +
      labs(y = ylab, title = letter)
  }
  pdf(file.path(OUT, "figS2_serpin.pdf"), width = 7.2, height = 3.1,
      useDingbats = FALSE)
  grid.arrange(grobs = list(
    serpin_panel("mal_serpine1_pct", "a  SERPINE1+ malignant cells", "Detection (%)"),
    serpin_panel("mal_serpinb2_pct", "b  SERPINB2+ malignant cells", "Detection (%)")),
    nrow = 1, ncol = 2)
  dev.off()
  message("wrote figS2_serpin.pdf")
}

# S3: genotype effect vs chemistry effect per readout
if (have("effect_sizes.tsv")) {
  READOUT_LAB <- c(
    pct_TNK            = "T/NK cells (%)",
    pct_Macrophage     = "Macrophages (%)",
    pct_CAF            = "CAFs (%)",
    tnk_exhaust        = "T/NK exhaustion score",
    tnk_effector       = "T/NK effector score",
    tnk_naive          = "T/NK naive score",
    mac_spp1_marco_pct = "SPP1+MARCO+ macrophages",
    mac_cxcl9_pct      = "CXCL9+ macrophages",
    mal_serpine1_pct   = "SERPINE1+ malignant",
    mal_serpinb2_pct   = "SERPINB2+ malignant")
  ef <- rd("effect_sizes.tsv")
  ef$ratio   <- ef$genotype_eff / ef$chemistry_eff
  ef$readout <- factor(READOUT_LAB[ef$readout],
                       levels = READOUT_LAB[ef$readout][order(ef$ratio)])

  pF <- ggplot(ef, aes(ratio, readout)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3,
               colour = "grey35") +
    geom_col(width = 0.60, fill = "#9ecae1", colour = "black",
             linewidth = 0.25) +
    scale_x_continuous(limits = c(0, max(ef$ratio) * 1.12),
                       expand = expansion(mult = c(0, 0))) +
    theme_minimal(base_size = 8.2) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "#ededed", linewidth = 0.2),
          axis.text = element_text(size = 6.6),
          axis.title = element_text(size = 7.2),
          plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
          plot.title.position = "plot",
          plot.caption = element_text(size = 5.8, hjust = 0, face = "italic"),
          plot.caption.position = "plot",
          plot.margin = margin(2, 3, 2, 2)) +
    labs(x = "Genotype effect / chemistry effect", y = NULL,
         title = "Effect-size ratio per readout",
         caption = "ratio < 1: the genotype effect is smaller than the technical\nchemistry effect; dashed line = equal size")

  pdf(file.path(OUT, "figS3_effects.pdf"), width = 7.2, height = 3.2,
      useDingbats = FALSE)
  print(pF)
  dev.off()
  message("wrote figS3_effects.pdf")
}

message("All PDFs written to: ", OUT)
