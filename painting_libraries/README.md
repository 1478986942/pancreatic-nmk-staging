# Open-source bioinformatics plotting libraries — reference index
# Cloned for learning Nature-style figure production from Falcomatà et al. 2026

## R Libraries (used in figs_nature.R)

| Repo | Primary Use | Key Functions to Study |
|------|------------|----------------------|
| **ComplexHeatmap** | Multi-panel heatmaps, clustered matrices with side annotations | `Heatmap()`, `HeatmapList()`, `annotate_heatmap()` |
| **circlize** | Circular plots (chord diagrams, circular heatmaps like ED Fig 1e) | `circos.initialize()`, `circos.trackplots()`, `circos.network()` |
| **SCP** | Single-cell pipeline visualization: UMAP contours, stacked bar, trend heatmap | `umapPlot()`, `stackedBarPlot()`, `trendLinePlot()` |
| **scRNAtoolVis** | One-click publication figures: volcano, corner axes UMAP, bubble plot | `jjVolcano()`, `cornerAxes()`, `bubblePlot()` |
| **ClusterGVis** | Expression trends + left-side heatmap + right-side GO/KEGG annotation | `drawClusters()` |
| **enrichplot** | Pathway enrichment: cnetplot, gseavez, dotplot, emaplot | `cnetplot()`, `gseaplot2()`, `dotplot()` |
| **survminer** | Survival curves with risk tables and log-rank p-values | `ggsurvplot()`, `add_risktable()` |
| **ggpubr** | Quick significance brackets on box/violin/bar plots | `stat_compare_means()`, `ggboxplot()` |
| **patchwork** | Multi-panel composition with strict axis alignment | `+` operator, `plot_layout()`, `plot_annotation()` |
| **aplot** | Layer-based panel assembly (bar + density strip, etc.) | `insert_panels()`, `insert_xaxes_graphs()` |
| **ggsci** | Journal-specific color palettes (Nature npg, Science aaas) | `scale_colour_npg()`, `scale_fill_jama()` |

## Python Libraries (for spatial / single-cell ecosystems)

| Repo | Primary Use | Key Functions to Study |
|------|------------|----------------------|
| **marsilea** | Declarative composite heatmap (R-like ComplexHeatmap in Python) | `Marsilea().add_row()` with dendrograms, bar plots as side annotations |
| **PyComplexHeatmap** | Python port of ComplexHeatmap for scATAC/scRNA large annotation matrices | `ComplexHeatmap()`, `split` parameter, `rowAnnotation` |
| **spatialdata-plot** | Spatial transcriptomics overlays (Xenium/Visium + multiplex IF) | `plot_points()`, `plot_images()`, coordinate-level overlay |
| **plotnine** | ggplot2 syntax in Python — learn layer grammar without switching languages | `ggplot()`, `geom_point()`, `facet_wrap()` equivalent |

## Biological Illustration Assets (free BioRender alternatives)

| Repo | What it has | How to use |
|------|------------|-----------|
| **bioicons** | 2000+ SVG biological icons (cells, mice, antibodies, DNA, slides) | Download SVG → drag into Inkscape/draw.io for workflow schematics |
| **SERVIER/bioart** | CC BY 4.0 medical illustrations (tissues, organs, instruments) | Cite Servier Medical Art; same source used in the Nature paper |

---

## Learning Path (4 weeks)

### Week 1: R basics → Nature bar plots
1. Read `ComplexHeatmap/R/Heatmap.R` (30 min skim for API surface)
2. Re-implement `nature_bar()` from scratch using only ggplot2
3. Add significance brackets with `ggpubr::stat_compare_means()`
4. Output: reproduce Fig 2d exactly

### Week 2: Multi-panel layout + survival
1. Learn `patchwork` by combining 4 sub-plots from Fig 2 into one PDF
2. Reproduce KM curve from Survival object using `survminer::ggsurvplot()`
3. Try `aplot` to stack a bar chart under a heat map (Fig 1h style)

### Week 3: Spatial & single-cell viz
1. Pick any public scRNA-seq dataset (GSE120575 is ICB + CD45+ immune cells only, good exercise)
2. Draw UMAP with `scRNAtoolVis::cornerAxes()` — practice reading the vignette
3. Draw pathway enrichment with `enrichplot::cnetplot()` — match Fig 3c style

### Week 4: Composite heatmap mastery
1. Clone ComplexHeatmap examples from GitHub: `ComplexHeatmap/extdata/` or Bioconductor gallery
2. Recreate a simplified version of ED Fig 1e (circular heatmap) using `circlize`
3. Build your own multi-annotation heatmap combining gene expression + clinical metadata

---

## Minimal Working Example — start here

```r
# === 1. Install (one-time) ===
# BiocManager::install(c("ComplexHeatmap", "circlize", "enrichplot", "survminer", "patchwork"))
# install.packages(c("ggpubr", "aplot", "ggsci", "SCP"))
# pip install marsilea PyComplexHeatmap plotnine spatialdata-plot

# === 2. Nature-style bar + points (from Fig 2d) ===
library(ggplot2)
library(patchwork)

ctrl <- c(0.57, 0.60, 0.48, 0.62, 0.52, 0.37, 0.24)
b2   <- c(0.11, 0.29, 0.28, 0.35, 0.25, 0.08, 0.15)
e1   <- c(0.09, 0.10, 0.09, 0.13, 0.08, 0.13, 0.10)

colors <- c(Control="#4d4d4d", Serpinb2KO="#1f77b4", Serpine1KO="#d62728")

mn  <- data.frame(group=factor(names(colors), levels=names(colors)), mean=sapply(list(ctrl,b2,e1), mean))
pts <- do.call(rbind, lapply(names(colors), function(g) data.frame(group=g, val=c(ctrl,b2,e1)[
  if(g=="Control") 1:7 else if(g=="Serpinb2KO") 8:14 else 15:21])))
pts$group <- factor(pts$group, levels=names(colors))

p <- ggplot(mn, aes(group, mean, fill=group)) +
  geom_col(width=0.55, colour="black", linewidth=0.28) +
  scale_fill_manual(values=colors, guide="none") +
  geom_point(data=pts, aes(group, val), shape=21, fill="white", colour="black", stroke=0.28, size=1.2,
             position=position_jitter(width=0.085, seed=1)) +
  annotate("text", x=0.5, y=max(unlist(list(ctrl,b2,e1))*1.12), label="d", hjust=0, vjust=1, size=2.9, fontface="bold") +
  labs(y="Tumour mass (g)", title=NULL) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="#ededed"),
        axis.line=element_line(linewidth=0.24),
        axis.text.x=element_text(angle=32, hjust=1, size=11))

ggsave("exercise_bar.png", p, width=5, height=4, dpi=300)
```
