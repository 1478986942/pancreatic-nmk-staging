# Practice figures — Nature-style from real project data
# Data: methods_paper/figs_reproduction/plot_tables/
# Libraries: ggplot2, patchwork, ggpubr (all base R)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

OUT <- "E:/pancreatic_cancer/submit/painting_libraries/practice_figs"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

COL_CTRL  <- "#4d4d4d"
COL_B     <- "#1f77b4"
COL_E     <- "#d62728"
THEME_NATURE <- function() {
  theme_minimal(base_size=13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour="#ededed", linewidth=0.15),
          axis.line = element_line(colour="black", linewidth=0.22),
          axis.text  = element_text(colour="black"),
          plot.margin = margin(8,8,8,8))
}

# =====================================================================
# Exercise 1: SERPIN detection proportion — bar + points (Fig 2 style)
# =====================================================================
df_serpin <- read.delim("methods_paper/figs_reproduction/plot_tables/serpin_detection.tsv")
df_serpin$genotype_f <- factor(df_serpin$genotype,
                                levels=c("control_KO","Serpinb2_KO","Serpine1_KO"))

mean_df <- aggregate(mal_serpine1_pct ~ genotype_f, df_serpin, mean)
names(mean_df)[2] <- "mean"
pts_df  <- data.frame(genotype_f=df_serpin$genotype_f,
                      value   =df_serpin$mal_serpine1_pct)
colors1 <- c(control_KO=COL_CTRL, Serpinb2_KO=COL_B, Serpine1_KO=COL_E)

p1 <- ggplot(mean_df, aes(genotype_f, mean, fill=genotype_f)) +
  geom_col(width=0.5, colour="black", linewidth=0.26) +
  scale_fill_manual(values=colors1, guide="none") +
  geom_point(data=pts_df, aes(genotype_f, value), shape=21, fill="white",
             colour="black", stroke=0.26, size=1.1,
             position=position_jitter(width=0.08, seed=7)) +
  annotate("text", x=0.5, y=max(mean_df$mean)*1.18, label="a",
           hjust=0, vjust=1, size=3.2, fontface="bold") +
  labs(y="SERPINE1+ malignant cells (%)", x=NULL) +
  THEME_NATURE() +
  scale_y_continuous(expand=c(0.02,0))
print(p1)
ggsave(file.path(OUT,"ex1_serpin_detection.png"), p1, width=4.2, height=3.5, dpi=300)


# =====================================================================
# Exercise 2: CAF subset composition — grouped bar (Fig 2f style)
# =====================================================================
df_caf <- read.delim("methods_paper/figs_reproduction/plot_tables/caf_subsets.tsv")
cat("CAF columns:", paste(names(df_caf), collapse=", "), "\n")
# Columns after read.delim: LRCC15..MyCAF, Matrix.CAF, MyCAF, apCAF, iCAF
caf_types_orig <- c("LRCC15..MyCAF","Matrix.CAF","MyCAF","apCAF","iCAF")
caf_labels     <- c("LRCC15+ MyCAF","Matrix CAF","MyCAF","apCAF","iCAF")
df_caf <- df_caf[df_caf$chemistry == "Flex", ]
long_df <- do.call(rbind, mapply(function(orig, label) {
  data.frame(caf_type=factor(label, levels=caf_labels),
             genotype_f=factor(df_caf$genotype,
                               levels=c("control_KO","Serpinb2_KO","Serpine1_KO")),
             pct=round(as.numeric(df_caf[[orig]]), 1))
}, orig=caf_types_orig, label=caf_labels, SIMPLIFY=FALSE))

p2 <- ggplot(long_df, aes(caf_type, pct, fill=genotype_f)) +
  geom_bar(stat="identity", position=position_dodge(0.82),
           width=0.72, colour="black", linewidth=0.22) +
  scale_fill_manual(values=colors1) +
  annotate("text", x=0.35, y=70, label="b",
           hjust=0, vjust=1, size=3.2, fontface="bold") +
  labs(x=NULL, y="CAF subset proportion (%)", fill=NULL) +
  THEME_NATURE() +
  theme(axis.text.x = element_text(size=9, angle=22, hjust=1),
        legend.position = "bottom",
        legend.key.size = unit(0.22, "cm"))
print(p2)
ggsave(file.path(OUT,"ex2_caf_subsets.png"), p2, width=5.5, height=3.8, dpi=300)


# =====================================================================
# Exercise 3: UMAP — myeloid cells sampled to 8k (Fig 3g style)
# =====================================================================
cat("Loading UMAP coords ...\n")
df_umap <- read.delim("methods_paper/figs_reproduction/plot_tables/myeloid_umap.tsv",
                       stringsAsFactors=FALSE)
set.seed(42)
n_sample <- min(8000, nrow(df_umap))
sub_umap <- df_umap[sample(nrow(df_umap), n_sample), ]
umap_colors <- c("#e74c3c","#3498db","#2ecc71","#f39c12",
                 "#9b59b6","#1abc9c","#e67e22","#7f8c8d",
                 "#ecf0f1","#95a5a6")

p3 <- ggplot(sub_umap, aes(UMAP1, UMAP2, colour=subset)) +
  geom_point(size=0.1, alpha=0.75) +
  scale_colour_manual(values=umap_colors) +
  annotate("text", x=min(sub_umap$UMAP1)+0.5, y=max(sub_umap$UMAP2)-0.5,
           label="c", size=3.2, fontface="bold", colour="black") +
  theme_void(base_size=13) +
  theme(legend.position = "bottom",
        legend.text    = element_text(size=9),
        legend.title   = element_blank(),
        legend.key.size= unit(0.2, "cm"),
        plot.margin    = margin(6,6,6,6))
print(p3)
ggsave(file.path(OUT,"ex3_myeloid_umap.png"), p3, width=5, height=4.2, dpi=300)


# =====================================================================
# Exercise 4: Combine three panels into one Figure-like layout
# =====================================================================
combined <- (p1 | p2) / p3
combined <- combined +
  plot_layout(heights=c(1, 1.3)) +
  plot_annotation(tag_levels="A", title="Practice: Nature-style figure assembly",
                  subtitle="Data from GSE289138 serpin_repro/processed")
print(combined)
ggsave(file.path(OUT,"ex4_combined.png"), combined, width=9, height=8, dpi=300)


# =====================================================================
# Exercise 5: Boxplot with ggpubr significance brackets
# Requires: install.packages("ggpubr")
# =====================================================================
if (requireNamespace("ggpubr", quietly=TRUE)) {
suppressPackageStartupMessages({ library(ggpubr) })

p5 <- ggboxplot(df_serpin, x="genotype_f", y="mal_serpine1_pct",
                fill="genotype_f", color="black", palette=colors1,
                add="jitter", jitter.width=0.15, point.shape=21, point.fill="white",
                font.label=list(face="bold", size=10)) +
  stat_compare_means(method="wilcox.test", label="p.signif",
                     comparisons=list(c("control_KO","Serpinb2_KO"),
                                      c("control_KO","Serpine1_KO"))) +
  scale_fill_manual(values=colors1, guide="none") +
  labs(x=NULL, y="SERPINE1+ malignant %") +
  THEME_NATURE()
print(p5)
ggsave(file.path(OUT,"ex5_boxplot_ggpubr.png"), p5, width=4.5, height=4, dpi=300)
cat("Exercise 5: saved ex5_boxplot_ggpubr.png\n")
} else {
  cat("SKIP ex5: ggpubr not installed. Run: install.packages('ggpubr')\n")
}


# =====================================================================
# Exercise 6: Simple dot/violin from serpin_detection (classic Nature style)
# =====================================================================
p6 <- ggplot(df_serpin, aes(genotype_f, mal_serpine1_pct, fill=genotype_f)) +
  geom_violin(alpha=0.4, colour="grey50", trim=FALSE, linewidth=0.2) +
  geom_boxplot(width=0.2, outlier.shape=NA, fill=NA) +
  geom_point(position=position_jitter(width=0.12, seed=3),
             shape=21, fill="white", colour="black", stroke=0.25, size=1.3) +
  scale_fill_manual(values=colors1, guide="none") +
  annotate("text", x=0.5, y=max(df_serpin$mal_serpine1_pct)*1.15,
           label="d", size=3.2, fontface="bold", hjust=0, vjust=1) +
  labs(y="SERPINE1+ malignant cells (%)", x=NULL) +
  THEME_NATURE() +
  ylim(0, max(df_serpin$mal_serpine1_pct)*1.2)
print(p6)
ggsave(file.path(OUT,"ex6_violin_dotplot.png"), p6, width=4.2, height=4, dpi=300)


cat("\nDone. All figures saved to:", OUT, "\n")
list.files(OUT)
