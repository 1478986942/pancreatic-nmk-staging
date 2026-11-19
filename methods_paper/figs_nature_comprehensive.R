# =====================================================================
# Comprehensive enhanced figures — all available data, varied chart types
# Run: Rscript methods_paper/figs_nature_comprehensive.R
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(gridExtra); library(grid)
})

ROOT   <- "E:/pancreatic_cancer/submit"
PT     <- file.path(ROOT, "methods_paper", "figs_reproduction", "plot_tables")
SRP    <- file.path(ROOT, "serpin_repro", "processed")
OUT    <- file.path(ROOT, "methods_paper", "figs_nature")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

BS <- 6.6; AX <- 5.2; TI <- 7.8; ROW_H <- 1.62
`%||%` <- function(a,b) if(is.null(a)) b else a
rd <- function(f) read.delim(file.path(PT,f), check.names=FALSE, stringsAsFactors=FALSE)
rd_srp <- function(f) read.delim(file.path(SRP,f), sep="\t", stringsAsFactors=FALSE)

GENO_LAB  <- c(control_KO="Control", Serpinb2_KO="Serpinb2 KO", Serpine1_KO="Serpine1 KO")
GENO_COL  <- c(control_KO="#4d4d4d", Serpinb2_KO="#1f77b4", Serpine1_KO="#d62728")
CHEM_LAB  <- c(classic="Classic (3')", Flex="Flex")

p_theme <- function(title="", xlab="", ylab="", ...) {
  theme_minimal(base_size=BS) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_line(colour="#ededed",linewidth=0.15),
          axis.text=element_text(size=AX-0.2), axis.title=element_text(size=AX+0.1),
          plot.title=element_text(size=TI,face="bold",hjust=0),
          plot.title.position="plot", plot.margin=margin(3,4,3,4), ...)
}

# ════════════════════════════════════════════════════════════════════
# Panel 1: Cluster composition stacked bar (all 6 libraries, 24 clusters)
# ════════════════════════════════════════════════════════════════════
message("Panel 1: Cluster composition ...")
comp_raw <- rd_srp("cluster_sample_composition.tsv")
rownames(comp_raw) <- comp_raw[,1]; comp_raw <- comp_raw[,-1]
comp <- as.matrix(comp_raw); comp[] <- as.numeric(comp)
col_totals <- colSums(comp)
comp_prop <- sweep(comp, 2, col_totals, "/")
comp_long <- as.data.frame(as.table(comp_prop))
names(comp_long) <- c("Cluster","Sample","Fraction")
comp_long$Cluster <- factor(comp_long$Cluster, levels=rownames(comp))
comp_long$Sample  <- factor(comp_long$Sample, levels=colnames(comp))
sample_lab <- sapply(strsplit(as.character(comp_long$Sample),"_"),
                     function(x) paste0(GENO_LAB[x[2]]," (",CHEM_LAB[x[3]],")"))
comp_long$Sample_lab <- factor(sample_lab, levels=unique(sample_lab))
pal_cluster <- setNames(scales::hue_pal()(min(24,nrow(comp))), rownames(comp))

p1 <- ggplot(comp_long[comp_long$Fraction > 0.003,], aes(x=Sample_lab, y=Fraction, fill=Cluster)) +
  geom_col(width=0.85, colour="white", linewidth=0.12) +
  scale_fill_manual(values=pal_cluster, guide=guide_legend(nrow=4, keywidth=0.4, keyheight=0.15, title=NULL)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.02))) +
  coord_flip() +
  labs(x=NULL, y="Proportion of total cells", title="a  Leiden cluster composition across 6 libraries") +
  p_theme()

pdf(file.path(OUT, "enhanced_01_cluster_composition.pdf"), width=9.0, height=3.4, useDingbats=FALSE)
print(p1)
dev.off()

# ════════════════════════════════════════════════════════════════════
# Panel 2: T-cell composition line plot (classic only, n=3)
# ════════════════════════════════════════════════════════════════════
message("Panel 2: T-cell composition ...")
tc <- rd("tcell_subsets.tsv")[rd("tcell_subsets.tsv")$chemistry=="classic",]
subset_cols <- c("CD8 Teff","CD8 Tex term","CD8 Tex prog","CD8 Tex int","CD8 Tn-like","CD4 Teff Gzmk+","CD4 Tn-like","CD4 Treg")
tc_long <- do.call(rbind, lapply(subset_cols, function(s) {
  data.frame(genotype=tc$genotype, subset=s, value=tc[[s]])
}))
tc_long$genotype <- factor(tc_long$genotype, levels=names(GENO_LAB))
tc_long$subset   <- factor(tc_long$subset, levels=subset_cols)
tc_colors <- setNames(c("#1f77b4","#d62728","#ff7f0e","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f"), subset_cols)

p2 <- ggplot(tc_long, aes(x=subset, y=value, group=genotype, colour=genotype)) +
  geom_line(linewidth=0.7) +
  geom_point(size=2.8, stroke=0.4) +
  scale_colour_manual(values=GENO_COL) +
  ylim(c(0, max(tc_long$value)*1.15)) +
  labs(x=NULL, y="T-cell fraction (%)", title="b  T-cell subset profile by genotype (classic)") +
  p_theme(legend.position="right", legend.key.size=unit(0.22,"cm")) +
  guides(colour=guide_legend(ncol=1))

pdf(file.path(OUT, "enhanced_02_tcell_profile.pdf"), width=7.2, height=4.0, useDingbats=FALSE)
print(p2)
dev.off()

# ════════════════════════════════════════════════════════════════════
# Panel 3: Myeloid heatmap (6 samples x 9 subsets)
# ════════════════════════════════════════════════════════════════════
message("Panel 3: Myeloid heatmap ...")
my <- rd("myeloid_subsets.tsv")
my_cols <- c("Arg1+","Cxcl9+","Spp1/Marco+","Proliferating","C1q+ homeostatic","DC-like","Folr2+","Il1b+","Monocyte-like")
tmp <- data.frame(genotype=my$genotype, chemistry=my$chemistry)
for(col in my_cols) tmp <- cbind(tmp, as.numeric(my[[col]]))
rownames(tmp) <- paste0(tmp$genotype,"/",tmp$chemistry)
my_mat <- as.matrix(round(tmp[,-c(1,2)],1))

my_long <- as.data.frame(as.table(my_mat))
names(my_long) <- c("Gen.Chem","Subset","Value")
my_long$Subset <- factor(my_long$Subset, levels=my_cols)
my_long$Gen.Chem <- factor(my_long$Gen.Chem, levels=rownames(my_mat))

p3 <- ggplot(my_long, aes(Subset, Gen.Chem, fill=Value)) +
  geom_tile(colour="white", linewidth=0.2) +
  geom_text(aes(label=sprintf("%.1f", Value)), size=2.6, colour="grey25") +
  scale_fill_gradientn(colours=c("#f5f5f5","#b3d9ff","#ffb347","#d62728"),
                       limits=c(0,max(my_long$Value)), na.value="#fafafa", name="Fraction (%)") +
  labs(x=NULL, y=NULL, title="c  Myeloid subset fractions — 6 libraries, classic & Flex") +
  p_theme(axis.text.x=element_text(size=AX-0.4,angle=25,hjust=1),
          axis.text.y=element_text(size=AX-0.3))

pdf(file.path(OUT, "enhanced_03_myeloid_heatmap.pdf"), width=7.2, height=ROW_H*2.2, useDingbats=FALSE)
print(p3)
dev.off()

# ════════════════════════════════════════════════════════════════════
# Panel 4: Effect-size ratio (genotype vs chemistry)
# ════════════════════════════════════════════════════════════════════
message("Panel 4: Effect ratio ...")
ef <- rd("effect_sizes.tsv")
READOUT_LAB <- c(pct_TNK="T/NK %", pct_Macrophage="Macrophage %", pct_CAF="CAF %",
                 tnk_exhaust="T/NK exhaustion", tnk_effector="T/NK effector",
                 tnk_naive="T/NK naive", mac_spp1_marco_pct="SPP1+MARCO+ macro%",
                 mac_cxcl9_pct="CXCL9+ macro%", mal_serpine1_pct="SERPINE1+ mal epi",
                 mal_serpinb2_pct="SERPINB2+ mal epi")
ef$readout_lab <- READOUT_LAB[ef$readout]
ef$ratio <- ef$genotype_eff / ef$chemistry_eff
ef$dir   <- ifelse(ef$ratio < 1, "Chemistry dominates", "Genotype dominates")
ef$dir   <- factor(ef$dir, levels=c("Chemistry dominates","Genotype dominates"))

p4 <- ggplot(ef, aes(x=readout_lab, y=ratio, fill=dir)) +
  geom_col(width=0.65, colour="black", linewidth=0.18) +
  geom_hline(yintercept=1, linetype="dashed", colour="#888", linewidth=0.25) +
  geom_text(aes(label=sprintf("%.2f",ratio)), vjust=-0.35, size=2.5) +
  scale_fill_manual(values=c("Chemistry dominates"="#9ecae1","Genotype dominates"="#f2a3a3")) +
  coord_flip() +
  labs(x=NULL, y="Genotype effect / Chemistry effect", title="d  Effect-size ratio: how much does chemistry matter?") +
  p_theme() +
  annotate("text", x=1, y=max(ef$ratio)*1.08, label="ratio < 1 = chemistry effect larger than genotype",
           size=3, colour="grey50", fontface="italic")

pdf(file.path(OUT, "enhanced_04_effect_ratio.pdf"), width=7.2, height=4.8, useDingbats=FALSE)
print(p4)
dev.off()

# ════════════════════════════════════════════════════════════════════
# Panel 5: SERPIN mosaic detection — two side-by-side bar charts
# ════════════════════════════════════════════════════════════════════
message("Panel 5: SERPIN mosaic ...")
sd <- rd("serpin_detection.tsv")
sd <- sd[order(sd$chemistry, match(sd$genotype,names(GENO_LAB))),]
sd$label <- paste0(GENO_LAB[sd$genotype]," [",sd$chemistry,"]")
sd$label <- factor(sd$label, levels=sd$label)

p5a <- ggplot(sd, aes(x=label, y=mal_serpine1_pct, fill=genotype)) +
  geom_col(position=position_dodge(0.75), width=0.65, colour="black", linewidth=0.22) +
  geom_text(aes(label=sprintf("%.1f%%", mal_serpine1_pct), group=genotype),
            position=position_dodge(0.75), vjust=-0.3, size=2.8) +
  scale_fill_manual(values=GENO_COL, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18)), limits=c(0,42)) +
  annotate("rect", xmin=2.5, xmax=5.5, ymin=0, ymax=Inf, fill="red", alpha=0.05) +
  annotate("text", x=4, y=40, label="SERPINE1 still detected in Serpine1-KO (mosaic)",
           size=3.2, colour="#b30000", fontface="bold") +
  labs(x=NULL, y="SERPINE1+ malignant cells (%)", title="e  SERPINE1 persistence in Serpine1-KO") +
  p_theme(axis.text.x=element_text(size=6,angle=30,hjust=1))

p5b <- ggplot(sd, aes(x=label, y=mal_serpinb2_pct, fill=genotype)) +
  geom_col(position=position_dodge(0.75), width=0.65, colour="black", linewidth=0.22) +
  geom_text(aes(label=sprintf("%.2f%%", mal_serpinb2_pct), group=genotype),
            position=position_dodge(0.75), vjust=-0.3, size=2.5) +
  scale_fill_manual(values=GENO_COL, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,2)), limits=c(0,max(sd$mal_serpinb2_pct)*1.5)) +
  labs(x=NULL, y="SERPINB2+ malignant cells (%)", title="f  SERPINB2 absent from all libraries (confirmed KO)") +
  p_theme(axis.text.x=element_text(size=6,angle=30,hjust=1))

pdf(file.path(OUT, "enhanced_05_serpin_mosaic.pdf"), width=7.2, height=4.6, useDingbats=FALSE)
grid.arrange(p5a, p5b, nrow=2)
dev.off()

# ════════════════════════════════════════════════════════════════════
# Panel 6: Lineage comparison dot plot (T cell vs Myeloid)
# ════════════════════════════════════════════════════════════════════
message("Panel 6: Lineage comparison ...")
key_markers <- list(T_cell=c("CD8 Teff","CD8 Tex term","CD8 Tex prog","CD8 Tex int"),
                    Myeloid=c("Cxcl9+","Spp1/Marco+","Proliferating","Arg1+"))
combo_rows <- list()
for(lineage in names(key_markers)) {
  src <- if(lineage=="T_cell") tc[tc$chemistry=="classic",] else my[my$chemistry=="classic",]
  for(mk in key_markers[[lineage]]) {
    combo_rows[[length(combo_rows)+1]] <- data.frame(
      lineage=lineage, marker=mk, value=src[[mk]],
      genotype=src$genotype, stringsAsFactors=FALSE)
  }
}
comb_df <- do.call(rbind, combo_rows)
comb_df$genotype <- factor(comb_df$genotype, levels=names(GENO_LAB))
comb_df$marker <- factor(comb_df$marker, levels=c(unlist(key_markers)))

p6 <- ggplot(comb_df, aes(x=value, y=marker, colour=genotype, shape=lineage)) +
  geom_point(size=3.2, stroke=0.4) +
  facet_wrap(~lineage, scales="free_y", ncol=1) +
  scale_colour_manual(values=GENO_COL) +
  scale_shape_manual(values=c(T_cell=16, Myeloid=17)) +
  scale_x_continuous(expand=expansion(mult=c(0.05,0.2)), limits=c(-2,62)) +
  labs(x="Fraction (%)", y=NULL, title="g  Key T-cell vs Myeloid readouts (classic only)") +
  p_theme(legend.position="right")

pdf(file.path(OUT, "enhanced_06_lineage_comparison.pdf"), width=7.2, height=5.0, useDingbats=FALSE)
print(p6)
dev.off()

# ════════════════════════════════════════════════════════════════════
# Panel 7: UMAP overview — cell type / genotype / chemistry
# ════════════════════════════════════════════════════════════════════
message("Panel 7: UMAP ...")
um <- rd("umap_coords.tsv")
set.seed(1)
if(nrow(um)>10000) um <- um[sample(nrow(um),10000),]
um$genotype_f <- factor(um$genotype, levels=names(GENO_LAB))
um$condition_f <- factor(ifelse(um$condition=="classic","classic (3')","Flex"),
                         levels=c("classic (3')","Flex"))
anno_lv <- sort(unique(um$anno))
anno_pal <- setNames(scales::hue_pal()(length(anno_lv)), anno_lv)

umap_base <- function(col, pal, title) {
  ggplot(um, aes(UMAP1, UMAP2, colour=.data[[col]])) +
    geom_point(size=0.12, alpha=0.8) +
    scale_colour_manual(values=pal) +
    labs(title=title, x=NULL, y=NULL) +
    theme_void(base_size=BS) +
    theme(legend.position="bottom", legend.text=element_text(size=5),
          legend.key.size=unit(0.12,"cm"),
          plot.title=element_text(size=TI*0.9,face="bold",hjust=0))
}

p7a <- umap_base("anno", anno_pal, "Cell type")
p7b <- umap_base("genotype_f", GENO_COL, "Genotype")
p7c <- umap_base("condition_f", c("classic (3')"="#1f77b4","Flex"="#d62728"), "Chemistry")

pdf(file.path(OUT, "enhanced_07_umap_overview.pdf"), width=9.5, height=3.0, useDingbats=FALSE)
grid.arrange(p7a, p7b, p7c, ncol=3)
dev.off()

message("All comprehensive enhanced panels written to: ", OUT)
