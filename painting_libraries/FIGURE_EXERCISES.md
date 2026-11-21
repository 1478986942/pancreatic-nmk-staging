# 可画的图清单 — 对应库 + 可用数据 + 练习建议

> 以 GSE202051、GSE300595、GSE289138 的已有产出为练习数据源，配合 17 个绘图库。

---

## 一、R 包可画图一览（按图型分组）

### 1. 柱状图 + 散点叠加（最常用）
**库**：基础 ggplot2 + `ggpubr::stat_compare_means()`  
**对应原文**：Fig 2d/e/h/i/k, Fig 3e/f/j/l, Fig 4a/d/e/f/j, Fig 5e/f/g/h  
**练习数据**：`methods_paper/figs_reproduction/plot_tables/myeloid_subsets.tsv`（SERPIN+巨噬细胞比例，每个患者一个值）  
```r
library(ggpubr)
p <- ggplot(df, aes(group, value, fill=group)) +
  geom_col(width=0.55, colour="black", linewidth=0.28) +
  geom_point(position=jitter(width=0.08), shape=21,
             fill="white", colour="black", stroke=0.28, size=1.2) +
  stat_compare_means(method="wilcox", label="p.signif")
```

### 2. Kaplan-Meier 生存曲线
**库**：`survminer::ggsurvplot()`  
**对应原文**：Fig 2l/n, ED 5e/l  
**练习数据**：用 GSE202051 的患者级 SERPIN 状态（高/低）构造生存变量（无真实随访，可造模拟数据练习）  
```r
library(survminer)
fit <- survfit(Surv(days, event) ~ serpin_group, data=df)
ggsurvplot(fit, data=df, risk.table=TRUE, pval=TRUE, conf.int=FALSE)
```

### 3. UMAP / tSNE 降维可视化
**库**：SCP（`FeatureDimPlot`）、scRNAtoolVis（`cornerAxes`）、基础 ggplot2  
**对应原文**：Fig 3g, ED 2e, ED 5f/g, ED 7c/f, ED 8p  
**练习数据**：`methods_paper/figs_reproduction/plot_tables/tcell_umap.tsv`（16万行，CD8 T细胞UMAP坐标）或 `myeloid_umap.tsv`（82万行）  
```r
df <- read.delim("figs_reproduction/plot_tables/myeloid_umap.tsv")
ggplot(df, aes(UMAP1, UMAP2, colour=cluster)) +
  geom_point(size=0.1, alpha=0.6) + theme_void()
```

### 4. 热图（矩阵 + 行/列注释）
**库**：`ComplexHeatmap`（首选）、enrichplot::heatplot（通路→基因热图）  
**对应原文**：Fig 3c（通路活性热图）、ED 1e（圆形热图）  
**练习数据**：可以造一个简单的表达矩阵，或直接读 `marco_corrected_test.tsv`  
```r
library(ComplexHeatmap)
mat <- matrix(rnorm(100), 10, 10)
ha_row <- rowAnnotation(cluster=cut(row.names(mat), 3), col=list(cluster=c(rep(1:3))))
Heatmap(mat, name="expr", row_annotation=ha_row,
        cluster_rows=TRUE, cluster_cols=TRUE,
        col=colorRampPalette(c("blue","white","red"))(100))
```

### 5. 气泡图（富集结果）
**库**：`enrichplot::dotplot()`  
**对应原文**：Fig 5c风格（通路富集气泡图）  
**练习数据**：用 `ClusterGVis` 或 enrichplot mock 对象  
```r
dotplot(enrich_result, showCategory=15, color=pAdjustMethod, size=GeneRatio)
```

### 6. 火山图（差异基因）
**库**：scRNAtoolVis（`jjVolcano`）、enrichplot::volplot  
**对应原文**：Fig 3d volcano  
```r
jjVolcano(de_results, x="log2FC", y="-log10(padj)", threshold=c(1, 0.05))
```

### 7. 通路网络图（Gene-Category Network）
**库**：`enrichplot::cnetplot()`  
**对应原文**：Fig 3c 风格的通路-基因关联  
```r
cnetplot(enrich_result, categorySize="pvalue", foldChange=fc_vector)
```

### 8. 山脊图 / GSEA 曲线
**库**：`enrichplot::ridgeplot()`  
**对应原文**：GSEA 富集曲线  
```r
ridgeplot(gsea_result, showCategory=20, fill="p.adjust")
```

### 9. 弦图 / 圆形聚类热图
**库**：`circlize`  
**对应原文**：ED Fig 1e 圆形热图  
```r
library(circlize)
circos.initialize With factions(...)
circos.trackHist(...) # 直方图
circos.trackLines(...) # 弦
circos.clear()
```

### 10. 堆叠柱状图 / 占比图
**库**：SCP（`stackedBarPlot`）、基础 ggplot2 (`position=position_fill`)  
**对应原文**：Fig 2f（T细胞亚群相对丰度）  
```r
ggplot(long_df, aes(x=subset, y=value, fill=genotype)) +
  geom_bar(position="fill", stat="identity")
```

### 11. 趋势折线 + 左右注释条
**库**：ClusterGVis（`drawClusters()`）  
**对应原文**：拟时序表达趋势  
```r
drawClusters(cluster_df, lineages=lineage_info, features=top_genes)
```

### 12. 多图拼版（严格对齐）
**库**：`patchwork` + `aplot`  
**对应原文**：任意多面板 Figure 组装  
```r
library(patchwork)
p1 + p2 / p3  # 上面两个，下面一个大panel
p1 | (p2 / p3)  # 左边一个，右边上下两个
```

---

## 二、Python 包可画图一览

| 图型 | 库 | 对应原文 | 练习数据 |
|------|----|---------|---------|
| 复合热图（行列注释+聚类树）| `marsilea` | Fig 3c风格 | `patient_level_serpin_marco.tsv` |
| ComplexHeatmap风格大矩阵 | `PyComplexHeatmap` | ED 1e | 同上 |
| 空间转录组叠加（Xenium/Visium）| `spatialdata-plot` | Fig 3b Xenium+t-CyCIF | 需下载 Zenodo 数据 |
| ggplot2语法（Python）| `plotnine` | 所有ggplot2图 | 任何TSV |

**marsilea 示例**：
```python
from marsilea import Marsilea, RowSideColors, ColumnSideColors
m = Marsilea()
m.add_matrix(expr_matrix)
m.add_row_cluster()
m.add_column_side_col(ColumnSideColors(meta["condition"]))
m.draw(output="heatmap.pdf")
```

---

## 三、已有数据 ↔ 推荐练习顺序

```
数据文件                                    可练的图                          难度
─────────────────────────────────────────────────────────────────────
tcell_umap.tsv (16万行)                    UMAP降维可视化                     ★★
myeloid_umap.tsv (82万行)                  UMAP降维（采样后）                 ★★★
caf_subsets.tsv                            柱状+散点（CAF比例）               ★
endo_subsets.tsv                           柱状+散点（内皮比例）              ★
serpin_detection.tsv                       比例条形图                         ★
effect_sizes.tsv (700B)                    森林图 / 效应量条形图              ★★
umap_coords.tsv (全细胞)                   UMAP + corner axes                ★★
marco_corrected_test.tsv                   热图                             ★★★
patient_level_serpin_marco.tsv             患者层聚类图                     ★★★★
GSE202051_qc_filtered.h5ad                 完整scRNA流程可视化               ★★★★★
```

---

## 四、三套完整练习方案（选一套开始）

### 方案A：从最简单的开始（30分钟）
目标：画出第1张 Nature 风格的柱状+散点图

```r
# 用 caf_subsets.tsv
df <- read.delim("methods_paper/figs_reproduction/plot_tables/caf_subsets.tsv")
ggplot(df, aes(group, value, fill=group)) +
  geom_col(width=0.55, colour="black", linewidth=0.28) +
  scale_fill_manual(values=c(Control="#4d4d4d", KO="#1f77b4"), guide="none") +
  geom_point(position=jitter(width=0.08, seed=42), shape=21,
             fill="white", colour="black", stroke=0.2, size=1) +
  labs(y="CAF proportion (%)", x=NULL) +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.y=element_line(colour="#ededed"),
        axis.text.x=element_text(angle=32, hjust=1))
ggsave("exercise_caf.png", width=4.5, height=3.5, dpi=300)
```

### 方案B：UMAP可视化（1小时）
目标：画出 myeloid 细胞的 UMAP，模拟 Fig 3g

```r
df <- read.delim("methods_paper/figs_reproduction/plot_tables/myeloid_umap.tsv")
# 随机采样1万点避免PDF爆炸
set.seed(1)
sub <- df[sample(nrow(df), min(10000, nrow(df))), ]
ggplot(sub, aes(UMAP1, UMAP2, colour=cluster)) +
  geom_point(size=0.12, alpha=0.7) +
  theme_void() +
  theme(legend.position="bottom", legend.text=element_text(size=10)) +
  scale_colour_brewer(palette="Set1")
ggsave("exercise_umap_myeloid.png", width=5, height=4.5, dpi=300)
```

### 方案C：综合练习——仿照 Fig 2 组合（2小时）
目标：把 CAF占比柱状图 + UMAP + 热图拼成一张 multi-panel figure

```r
library(patchwork)
# p1: CAF比例柱状图
# p2: UMAP图
# p3: 简单热图（用 serpin_detection.tsv 做）
# 组合：p1 | (p2 / p3)
# 加标题：plot_annotation(tag_levels="a")
```

---

## 五、下一步

想先从哪套方案开始？我可以：
1. **直接帮你跑方案A**（生成一张 CAF 柱状图 PNG/PDF）
2. **带你看 UMAP 的代码细节**（采样策略、theme设置）
3. **先安装这些 R 包**（需要 BiocManager）然后动手练习
