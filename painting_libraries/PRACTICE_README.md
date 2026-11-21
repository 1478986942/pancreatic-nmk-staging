# GSE202051 练习数据说明

**数据**: Human PDAC snRNA-seq, Shen et al. Nature Medicine 2021
**位置**: `gse202051_hum/processed/GSE202051_qc_filtered.h5ad`
**规模**: 80,000 核 / 22,164 基因 / 44 患者 / 18 naive / 25 treated

---

## 可用变量一览

### 细胞类型注释
| 列名 | 含义 | 取值数 | 用法 |
|------|------|--------|------|
| `broad_celltypes` | 大类: Malignant/Immune/Fibroblast/Endothelial/Endocrine/Unknown | 7 | Fig B/Umap基础色 |
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

`practice_gse202051.py` 生成 6 张练习图：

```
gse202051_A_celltype_composition.png      宽柱+高柱组合，细胞类型占比
gse202051_B_umap_broadcelltype.png        UMAP 按 broad_celltypes 着色
gse202051_C_umap_immune_response.png      免疫细胞 UMAP × response
gse202051_D_serpin_fraction_by_response.png SERPIN信号 violin + jitter
gse202051_E_marker_dotplot.png            Top marker log2FC vs p_adj 散点
gse202051_F_marker_heatmap.png            Top 8 markers × cell types heatmap
```
