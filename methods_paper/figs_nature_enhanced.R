# =====================================================================
# Enhanced panels — supplemental figures with more visual variety.
# Place alongside figs_nature/ for comparison.
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

ROOT   <- "E:/pancreatic_cancer/submit"
PT     <- file.path(ROOT, "methods_paper", "figs_reproduction", "plot_tables")
OUT    <- file.path(ROOT, "methods_paper", "figs_nature")

BS <- 6.6; AX <- 5.2; TI <- 7.8; ROW_H <- 1.62
`%||%` <- function(a, b) if (is.null(a)) b else a

read_tsv <- function(f) {
  read.delim(file.path(PT, f), check.names = FALSE, stringsAsFactors = FALSE)
}

GENO_LAB <- c(control_KO = "Control", Serpinb2_KO = "Serpinb2 KO", Serpine1_KO = "Serpine1 KO")
GENO_COL <- c(control_KO = "#4d4d4d", Serpinb2_KO = "#1f77b4", Serpine1_KO = "#d62728")

# ------------------------------------------------------------------
# Panel f redone: T-cell composition as horizontal dot-blot
# ------------------------------------------------------------------
tc <- read_tsv("tcell_subsets.tsv")
tc <- tc[tc$chemistry == "classic", ]
subset_cols <- c("CD8 Teff", "CD8 Tex term", "CD8 Tex prog", "CD8 Tex int",
                 "CD4 Teff Gzmk+", "CD4 Tn-like", "CD4 Treg", "CD8 Tn-like")
tc_long <- do.call(rbind, lapply(subset_cols, function(s) {
  data.frame(genotype = tc$genotype, subset = s, value = tc[[s]], n = tc$n)
}))
tc_long$genotype <- factor(tc_long$genotype, levels = names(GENO_LAB))
tc_long$subset   <- factor(tc_long$subset, levels = subset_cols)

tc_pal <- stats::setNames(
  c("#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f"),
  subset_cols)

p_f_dot <- ggplot(tc_long, aes(x = value, y = genotype, colour = subset, size = n)) +
  geom_point(shape = 16, stroke = 0.3) +
  geom_segment(aes(xend = 0, yend = genotype), colour = "#d0d0d0", linewidth = 0.15) +
  scale_colour_manual(values = tc_pal) +
  scale_size(range = c(2.2, 5.5), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, 62)) +
  theme_minimal(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#ededed", linewidth = 0.15),
        axis.text.y = element_text(size = AX - 0.2),
        axis.text.x = element_text(size = AX),
        axis.title = element_text(size = AX + 0.1),
        legend.position = "right",
        legend.key.size = unit(0.2, "cm"),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 4, 2, 2)) +
  labs(x = "Fraction of total T cells (%)", y = NULL, title = "f*",
       colour = "Subset", size = "n cells") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 4)))

pdf(file.path(OUT, "fig2f_dotplot.pdf"), width = 7.2, height = ROW_H * 2.2, useDingbats = FALSE)
print(p_f_dot)
dev.off()
message("wrote fig2f_dotplot.pdf")

# ------------------------------------------------------------------
# Panel g redone: Myeloid heatmap (6×9 matrix)
# ------------------------------------------------------------------
my <- read_tsv("myeloid_subsets.tsv")
my_cols <- c("Arg1+","Cxcl9+","Spp1/Marco+","Proliferating",
             "C1q+ homeostatic","DC-like","Folr2+","Il1b+","Monocyte-like")
tmp <- data.frame(genotype = my$genotype, chemistry = my$chemistry)
for (col in my_cols) tmp <- cbind(tmp, as.numeric(my[[col]]))
rownames(tmp) <- paste0(tmp$genotype, "/", tmp$chemistry)
my_mat <- as.matrix(round(tmp[, -c(1,2)], 1))
my_long <- as.data.frame(as.table(my_mat))
names(my_long) <- c("Gen.Chem", "Subset", "Value")
my_long$Subset <- factor(my_long$Subset, levels = my_cols)
my_long$Gen.Chem <- factor(my_long$Gen.Chem, levels = rownames(my_mat))

p_g_hm <- ggplot(my_long, aes(Subset, Gen.Chem, fill = Value)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", Value)), size = 2.8, colour = "grey20") +
  scale_fill_gradientn(colours = c("#f0f0f0","#9ecae1","#f2a541","#d62728"),
                       limits = c(0, max(my_long$Value)), na.value = "white") +
  theme_minimal(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = AX - 0.4, angle = 25, hjust = 1),
        axis.text.y = element_text(size = AX - 0.3),
        axis.title = element_text(size = AX + 0.1),
        legend.key.size = unit(0.2, "cm"),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 3, 2, 2)) +
  labs(x = NULL, y = NULL, title = "g*", fill = "Fraction (%)")

pdf(file.path(OUT, "fig4g_heatmap.pdf"), width = 7.2, height = ROW_H * 1.6, useDingbats = FALSE)
print(p_g_hm)
dev.off()
message("wrote fig4g_heatmap.pdf")

# ------------------------------------------------------------------
# SERPIN detection heat-strip (redone figS2 style)
# ------------------------------------------------------------------
sd <- read_tsv("serpin_detection.tsv")
sd <- sd[order(sd$chemistry, match(sd$genotype, names(GENO_LAB))), ]
sd$label <- paste0(GENO_LAB[sd$genotype], "\n[", sd$chemistry, "]")
sd$label <- factor(sd$label, levels = sd$label)

p_s2a <- ggplot(sd, aes(x = label, y = "", fill = mal_serpine1_pct)) +
  geom_tile(colour = "white", linewidth = 0.4, height = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", mal_serpine1_pct)),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_gradient(low = "#e8e8e8", high = "#d62728", na.value = "#f4f4f4",
                      name = "SERPINE1+ mal epi (%)") +
  geom_vline(xintercept = c(1.5, 3.5), colour = "#aaa", linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = 1.5, y = 1.25, label = "classic (3')", size = 2.8, colour = "grey40") +
  annotate("text", x = 4.5, y = 1.25, label = "Flex", size = 2.8, colour = "grey40") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(size = 6, angle = 30, hjust = 1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 9, face = "bold", hjust = 0),
        plot.margin = margin(10, 10, 6, 10)) +
  labs(title = "a  SERPINE1+ malignant cells  (mosaic KO check)")

p_s2b <- ggplot(sd, aes(x = label, y = "", fill = mal_serpinb2_pct)) +
  geom_tile(colour = "white", linewidth = 0.4, height = 0.65) +
  geom_text(aes(label = sprintf("%.2f%%", mal_serpinb2_pct)),
            size = 3.0, colour = "white", fontface = "bold") +
  scale_fill_gradient(low = "#e8e8e8", high = "#1f77b4", na.value = "#f4f4f4",
                      name = "SERPINB2+\nmal epi (%)") +
  geom_vline(xintercept = c(1.5, 3.5), colour = "#aaa", linewidth = 0.3, linetype = "dashed") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title = element_blank(),
        plot.margin = margin(10, 10, 6, 10)) +
  labs(title = "b  SERPINB2+ malignant cells")

pdf(file.path(OUT, "figS2_heatstrip.pdf"), width = 7.2, height = 3.2, useDingbats = FALSE)
grid.arrange(p_s2a, p_s2b, nrow = 2)
dev.off()
message("wrote figS2_heatstrip.pdf")

# ------------------------------------------------------------------
# Combined T-cell + Myeloid comparison dot-blot (new panel)
# Key readouts across both lineages
# ------------------------------------------------------------------
key_markers <- c("CD8 Teff","CD8 Tex term","Cxcl9+","Spp1/Marco+","Proliferating","Arg1+")
# build combined df from classic chemistry only
tc_classic <- tc[tc$chemistry == "classic", ]
my_all <- my
my_classic <- my_all[my_all$chemistry == "classic", ]

comb_df <- rbind(
  data.frame(lineage = "T cell", marker = "CD8 Teff", value = tc_classic[["CD8 Teff"]], genotype = tc_classic$genotype),
  data.frame(lineage = "T cell", marker = "CD8 Tex term", value = tc_classic[["CD8 Tex term"]], genotype = tc_classic$genotype),
  data.frame(lineage = "Myeloid", marker = "Cxcl9+", value = my_classic[["Cxcl9+"]], genotype = my_classic$genotype),
  data.frame(lineage = "Myeloid", marker = "Spp1/Marco+", value = my_classic[["Spp1/Marco+"]], genotype = my_classic$genotype),
  data.frame(lineage = "Myeloid", marker = "Proliferating", value = my_classic[["Proliferating"]], genotype = my_classic$genotype),
  data.frame(lineage = "Myeloid", marker = "Arg1+", value = my_classic[["Arg1+"]], genotype = my_classic$genotype)
)
comb_df$genotype <- factor(comb_df$genotype, levels = names(GENO_LAB))
comb_df$marker <- factor(comb_df$marker, levels = key_markers)

pal_lineage <- c("T cell" = "#1f77b4", "Myeloid" = "#d62728")

p_combo <- ggplot(comb_df, aes(x = value, y = marker, colour = genotype, shape = lineage)) +
  geom_point(size = 3.5, stroke = 0.4) +
  geom_hline(yintercept = seq_along(key_markers) - 0.5, colour = "#e0e0e0", linewidth = 0.2) +
  scale_colour_manual(values = GENO_COL) +
  scale_shape_manual(values = c("T cell" = 16, "Myeloid" = 17)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(-5, 60)) +
  theme_minimal(base_size = BS) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.15),
        axis.text.y = element_text(size = AX - 0.2),
        axis.text.x = element_text(size = AX),
        axis.title = element_text(size = AX + 0.1),
        legend.position = "right",
        legend.key.size = unit(0.22, "cm"),
        plot.title = element_text(size = TI, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 4, 2, 2)) +
  labs(x = "Fraction (%)", y = NULL, title = "T cell vs Myeloid key readouts (classic)",
       colour = "Genotype", shape = "Lineage")

pdf(file.path(OUT, "fig_combo_keyreadouts.pdf"), width = 7.2, height = 4.0, useDingbats = FALSE)
print(p_combo)
dev.off()
message("wrote fig_combo_keyreadouts.pdf")

message("All enhanced panels written.")
