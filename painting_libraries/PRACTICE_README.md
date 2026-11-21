# GSE202051 练习数据说明

**数据**: Human PDAC snRNA-seq, Shen et al. Nature Medicine 2021
**位置**: `gse202051_hum/processed/GSE202051_qc_filtered.h5ad`
**规模**: 80,000 核 / 22,164 基因 / 44 患者 / 18 naive / 25 treated

---

## 可用变量一览

### 细胞类型注释
| 列名 | 含义 | 取值数 | 用法 |
|------|------|--------|------|
| `broad_celltypes` | 大类: Epithelial/Fibroblast/Immune/Endothelial/Endocrine/Schwann/unknown（注意是小写 `unknown`；**恶性细胞标为 `Epithelial`，数据中不存在 `Malignant` 标签**） | 7 | Fig B/Umap基础色 |
| `new_celltypes` | 细分类: 45亚型 (Malignant/Epid-Ductal/etc.) | 45 | 进阶聚类分析 |
| `Level 3 Annotation` | 三级注释 | 39 | 与 new_celltypes 对应 |
| `leiden` | Seurat 聚类 | 31 | 参考聚类编号 |
| `leiden_scvi` | SCVI 聚类 | 17 | 去批次后的聚类 |

### 治疗/反应信息
| 列名 | 含义 | 取值 | 用法 |
|------|------|------|------|
| `status` | 采样时状态 | treated / naive | Fig A 堆叠柱X轴 |
| `response` | 对新辅助治疗的反应 | Moderate/Poor/Minimal response, Untreated | Fig C/D 关键分组 |
| `pid` | 患者ID | 43个 | 合并patient-level统计 |
| `sampleid` | 样本ID | 44个 | 单样本QC用 |

### 分子评分
| 列名 | 含义 | 范围 | 用法 |
|------|------|------|------|
| `Moffitt_basal` | 基底样评分 | -0.1~0.4 | Fig E dotplot分组依据 |
| `Moffitt_classical` | 经典样评分 | -0.2~0.3 | 同上 |
| `Bailey_squamous` | 鳞状评分 | -0.1~0.2 | 同上 |
| `Bailey_progenitor` | 祖细胞评分 | -0.04~0.26 | 同上 |
| `Collison_QM` | 增殖QM评分 | -0.14~0.1 | 同上 |
| `COLLISON_CLASSICAL` | 经典分型 | -0.18~0.75 | 同上 |
| `MALIGNANT CELLS` | 恶性程度评分 | -0.03~0.24 | 肿瘤细胞鉴定 |
| `cnv_score` | CNV 负荷 | ~0.01~0.02 | CNV阳性判断 |

### CAF 亚型评分（多个坐标系）
| 列名 | 含义 |
|------|------|
| `Tuveson_iCAF` / `Tuveson_mCAF` | Tuveson 标准 CAF 分类 |
| `Neuzillet_CAFa/b/c/d` | Neuzillet 四型 CAF |
| `Davidson_CAF1/2/3` | Davidson 三型 CAF |
| `Pan_Stellate` / `Quiescent_Stellate` / `Activated_Stellate` | 胰腺星状细胞谱系 |
| `iCAF` / `myCAF` / `apCAF` | 简化CAF分类 |

---

## UMAP 坐标

| obsm 键 | 算法 | 推荐用途 |
|---------|------|---------|
| `X_umap` | Seurat UMAP | 默认，彩色面板用此 |
| `X_scVI` | SCVI 去批次 UMAP | 跨批次比较用 |
| `X_pca` | PCA | 快速预览 |
| `X_harmony` | Harmony 校正 | 批次校正后探索 |
| `X_cnv_umap` | CNV-aware | 只看非恶性细胞时 |

---

## 练习脚本说明

`practice_gse202051.py` 生成 9 张图（A–H + A4 合成图）：

```
gse202051_A_celltype_composition.png        双堆叠柱：左=每类型 naive/treated 占比，右=每状态细胞类型构成
gse202051_B_umap_broadcelltype.png          UMAP 按 broad_celltypes 着色
gse202051_C_umap_immune_response.png        免疫细胞 UMAP × response
gse202051_D_serpin_fraction_by_response.png SERPIN 信号 violin + 患者级统计
gse202051_E_marker_dotplot.png              marker dotplot：行=细胞类型，列=top基因，点=%表达，色=z-score
gse202051_F_marker_heatmap.png              34 marker × 7 谱系 热图（seaborn）
gse202051_G_umap_newcelltypes.png           45 亚型 UMAP：谱系定色相 + 谱系内明度分档 + 数字标注
gse202051_H_heatmap_complexheatmap.png      F 的 R/ComplexHeatmap 重绘版
gse202051_Figure1_A4.png                    A–F 按英寸排版的 A4 整图
heatmap_zscore_matrix.tsv                   → 供 R 读取的 z-score 矩阵
heatmap_row_annotation.tsv                  → 基因 → 谱系
heatmap_col_annotation.tsv                  → 谱系 → 核数
```

`practice_heatmap_complexheatmap.R` 只读取上述三个 TSV 重绘，**不重算任何统计**，
因此 Python 与 R 两版热图可比。

### 面板 D 的统计口径（重要）

- 推断单位是**患者**，不是细胞核。每患者先取 malignant（`broad_celltypes=='Epithelial'`）
  核的 SERPIN 信号均值，再在患者层面比较。
- 每组患者数：Moderate 7 / Minimal 10 / Poor 6（仅纳入 ≥15 个恶性核的患者）。
- **Kruskal–Wallis p = 0.339**；两两双侧 Mann–Whitney：
  Moderate vs Poor p = 0.445，Moderate vs Minimal p = 0.813，Minimal vs Poor p = 0.118。
- **结论：未检出差异。** 不得换检验、改阈值或事后分层去凑显著。
  名义 p 值本身也是探索性的，且未做多重比较校正。

### 练习中踩到的坑（已修）

1. A 面板归一化方向：按行归一的百分比不能拿去画按列堆叠的图（曾出现 353%）。
2. `min(p, 1e-300)` 会把所有 p 压到同一下限，log 轴整片消失；应为 `max(p, 1e-300)`。
3. 45 类 UMAP 用 `tab20` 相邻色会把 Ductal 与 Malignant 画成同一浅青；
   改为谱系定色相 + 谱系内明度分档。
4. `imshow` 拼图保持等比会在每行下方留大片空白；改为裁剪白边后按真实宽高比
   以英寸定位 axes。
5. `matplotlib` 没有 `'grey20'`（那是 ggplot 写法），应为 `'#333333'`。
