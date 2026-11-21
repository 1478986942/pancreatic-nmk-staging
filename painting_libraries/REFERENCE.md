# Bioinformatics Plotting Libraries — Quick Reference

All cloned shallow (depth=1). Each is a standalone reference library for learning Nature-style figure production.

## R Packages (install via BiocManager / install.packages)

```r
BiocManager::install(c("ComplexHeatmap", "circlize", "enrichplot",
                        "SCP", "ClusterGVis", "survminer"))
install.packages(c("ggpubr", "patchwork", "aplot", "ggsci"))
```

### Which one to open first?

| 想学什么 | 打开哪个仓库的 vignette |
|---------|----------------------|
| 热图 + 聚类树 + 行/列注释 | `ComplexHeatmap/vignettes/` |
| 圆形/弦图（对应 ED Fig 1e）| `circlize/demo/` |
| UMAP + 轮廓线 + 气泡图 | `SCP/vignettes/SCP_vignette.Rmd` |
| 火山图 + 左下角坐标轴 | `scRNAtoolVis/vignettes/scRNAtoolVis.Rmd` |
| 趋势线 + 左右热图+GO注释 | `ClusterGVis/R/drawClusters.R` |
| 通路网络图（对应 Fig 3c 风格）| `enrichplot/R/cnetplot.R` |
| Kaplan-Meier + 风险表 | `survminer/R/ggsurvplot.R` |
| 显著性括号一键添加 | `ggpubr/R/stat_compare_means.R` |
| 多图拼版严格对齐 | `patchwork/R/layout.R` |

## Python Packages (pip install)

```bash
pip install marsilea PyComplexHeatmap plotnine
# spatialdata-plot requires scverse ecosystem:
# pip install spatialdata spatialdata-plot
```

## Biological Illustration Assets

| 资源 | 地址 | 用途 |
|------|------|------|
| **bioicons** | `painting_libraries/bioicons/` | 2000+ SVG图标，拖进Inkscape直接用 |
| **Servier Medical Art** | https://smart.servier.com | CC BY 4.0医学插画，Nature论文引用来源 |

---

## Nature Figure 配色速查

```r
C_CTRL  <- "#4d4d4d"   # 灰色：对照组
C_B     <- "#1f77b4"   # 蓝色：第一处理组
C_E     <- "#d62728"   # 红色：第二处理组
C_MIX   <- "#f2a541"   # 橙色：混合/联合
C_NET   <- "#9467bd"   # 紫色：第三组
C_OLIVE <- "#8c564b"   # 棕色：第四组
C_TEAL  <- "#17becf"   # 青色：第五组
```

取自 ColorBrewer Set1，Nature 编辑默认接受。
