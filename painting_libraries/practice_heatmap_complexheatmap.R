# =====================================================================
# practice_heatmap_complexheatmap.R
#
# Redraws panel F of practice_gse202051.py in R with ComplexHeatmap, so
# the Python (seaborn) and R (ComplexHeatmap) versions can be compared on
# exactly the same numbers.
#
# Input  : painting_libraries/practice_figs_gse202051/heatmap_*.tsv
#          (written by practice_gse202051.py -- no re-computation here)
# Output : painting_libraries/practice_figs_gse202051/
#          gse202051_H_heatmap_complexheatmap.png
#
# What to look at when comparing with the Python version:
#   1. row_split draws visible lineage blocks + per-block gene order
#   2. column annotation bar shows n nuclei per lineage, so a "hotter"
#      column can be checked against how much data it is based on
#   3. no per-cell statistics are involved -- this is a descriptive
#      mean-expression z-score, NOT a differential test
# =====================================================================

suppressPackageStartupMessages({
  has_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
  has_cr <- requireNamespace("circlize", quietly = TRUE)
  if (has_ch && has_cr) {
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  } else {
    message("ComplexHeatmap/circlize missing -- falling back to pheatmap")
    library(pheatmap)
  }
})

fig_dir <- "painting_libraries/practice_figs_gse202051"
f_mat <- file.path(fig_dir, "heatmap_zscore_matrix.tsv")
f_row <- file.path(fig_dir, "heatmap_row_annotation.tsv")
f_col <- file.path(fig_dir, "heatmap_col_annotation.tsv")

stopifnot(file.exists(f_mat), file.exists(f_row), file.exists(f_col))

mat <- as.matrix(read.delim(f_mat, row.names = 1, check.names = FALSE))
storage.mode(mat) <- "double"

row_ann <- read.delim(f_row, stringsAsFactors = FALSE)
col_ann <- read.delim(f_col, stringsAsFactors = FALSE)

# ---- hard assertions: annotations must cover the matrix exactly -------
stopifnot(!anyDuplicated(row_ann$gene_id))
stopifnot(all(rownames(mat) %in% row_ann$gene_id))
stopifnot(setequal(colnames(mat), col_ann$celltype))
stopifnot(all(is.finite(mat)))
stopifnot(all(abs(mat) < 1e3))
stopifnot(nrow(mat) >= 10, ncol(mat) >= 3)

row_ann <- row_ann[match(rownames(mat), row_ann$gene_id), , drop = FALSE]
rownames(row_ann) <- row_ann$gene_id
col_ann <- col_ann[match(colnames(mat), col_ann$celltype), , drop = FALSE]

# order lineage blocks the same way as the Python figure
lineage_levels <- colnames(mat)
row_ann$lineage <- factor(row_ann$lineage, levels = lineage_levels)
stopifnot(!any(is.na(row_ann$lineage)))

# ---- colours ---------------------------------------------------------
pal_lineage <- c(
  Epithelial  = "#d62728",
  Fibroblast  = "#ff7f0e",
  Immune      = "#2ca02c",
  Endothelial = "#1f77b4",
  Endocrine   = "#9467bd",
  Schwann     = "#8c564b",
  unknown     = "#bcbd22"
)
pal_lineage <- pal_lineage[names(pal_lineage) %in% lineage_levels]

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

out_png <- file.path(fig_dir, "gse202051_H_heatmap_complexheatmap.png")

if (exists("Heatmap")) {
  ha_top <- HeatmapAnnotation(
    `n nuclei` = anno_barplot(
      col_ann$n_cells,
      gp = gpar(fill = pal_lineage[col_ann$celltype], col = NA),
      border = FALSE,
      axis_param = list(gp = gpar(fontsize = 6)),
      height = unit(1.1, "cm")
    ),
    annotation_name_gp = gpar(fontsize = 7),
    annotation_name_side = "left",
    show_legend = FALSE
  )

  ra <- rowAnnotation(
    lineage = row_ann$lineage,
    col = list(lineage = pal_lineage),
    show_annotation_name = FALSE,
    simple_anno_size = unit(2.4, "mm"),
    annotation_legend_param = list(
      lineage = list(title = "Lineage", title_gp = gpar(fontsize = 7),
                     labels_gp = gpar(fontsize = 6),
                     grid_width = unit(3, "mm"))
    )
  )

  ht <- Heatmap(
    mat,
    name = "Z-score",
    col = col_fun,
    row_split = row_ann$lineage,
    cluster_row_slices = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = ha_top,
    left_annotation = ra,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 6, fontface = "italic"),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 35,
    column_names_side = "bottom",
    row_names_side = "right",
    column_title = "Top marker genes (mean expression z-score across lineages)",
    column_title_gp = gpar(fontsize = 7.5),
    heatmap_legend_param = list(
      title = "Z-score", title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6),
      legend_height = unit(2.3, "cm"),
      grid_width = unit(3, "mm")
    ),
    rect_gp = gpar(col = "white", lwd = 0.5)
  )

  png(out_png, width = 7.4, height = 6.6, units = "in", res = 300)
  draw(ht,
       column_title = "H",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"),
       ht_gap = unit(3, "mm"))
  grid.text("ComplexHeatmap 2.x redraw of the Python/seaborn panel F",
            x = unit(1, "npc") - unit(2, "mm"),
            y = unit(1, "npc") - unit(2, "mm"),
            just = c("right", "top"),
            gp = gpar(fontsize = 5.5, col = "grey40"))
  dev.off()

} else {
  # ---- fallback: pheatmap --------------------------------------------
  ann_col <- data.frame(lineage = col_ann$celltype, row.names = col_ann$celltype)
  ann_row <- data.frame(lineage = row_ann$lineage, row.names = rownames(mat))
  gaps <- cumsum(table(row_ann$lineage))[seq_len(nlevels(row_ann$lineage) - 1)]
  pheatmap(
    mat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    breaks = seq(-2, 2, length.out = 102),
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_row = ann_row, annotation_col = ann_col,
    annotation_colors = list(lineage = pal_lineage),
    gaps_row = gaps,
    fontsize = 7, fontsize_row = 6,
    filename = out_png, width = 7.4, height = 6.6, res = 300
  )
}

cat("WROTE", out_png, "\n")
cat("MATRIX", nrow(mat), "genes x", ncol(mat), "lineages\n")
cat("R", R.version.string, "\n")
cat(if (has_ch) paste("ComplexHeatmap",
                      as.character(utils::packageVersion("ComplexHeatmap")))
    else "pheatmap fallback", "\n")
