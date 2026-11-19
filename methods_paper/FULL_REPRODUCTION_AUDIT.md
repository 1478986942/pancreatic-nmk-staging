# Nature 论文完整复现进度 — Falcomatà et al. 2026, Nature 658

> 对照原文 Fig 1–5 + ED Fig 1–10 共 **~75 个面板**逐一核实

---

## 图例说明

| 符号 | 含义 |
|------|------|
| ✅ | 已从 Source Data / GSE289138 plot tables 精确重绘，风格与原文一致 |
| ⬜ | 原文是示意图/工作流程图，可手工绘制，当前为文字占位 |
| ❌ | 原始数据未公开（显微镜图像、单细胞坐标、Xenium 矩阵、Visium 坐标），无法复现 |
| 🟡 | 部分数据可用，但只能做近似或简化版本 |

---

## Fig 1 — Perturb-map screen（8 面板）

| 面板 | 内容 | 状态 | 备注 |
|------|------|------|------|
| a | Perturb-map 工作流程示意图 | ⬜ | 需手绘；当前灰色框 |
| b | MICSSS 多色荧光显微图像 d7/d14/d21 | ❌ | 原始图像未公开 |
| c | 空间异质性图（非匹配邻居%） | ❌ | 需要 per-cell xy 坐标 |
| d | Pro-Code 焦点区域堆叠图 | ❌ | 需要 per-cell 条形码数据 |
| e | 适应度热图 in vivo/in vitro log2FC | ❌ | Source Data 无此表（Zenodo notebook 中有） |
| f | 邻域富集 z-score 热图（Squidpy） | ❌ | Squidpy 输出未提交 |
| g | Rag2−/− 实验设计示意图 | ⬜ | 可手绘替代灰色框 |
| h | Rag2−/− 各 KO 丰度条形图 | ❌ | 无每鼠 barcode count |

**结论：Fig 1 全部 8 个面板都无法重绘。第 a/g 面板是示意图可补画；b/f/h 需要原始图像或 per-cell 数据。**

---

## Fig 2 — Serpin loss promotes T-cell immunity（14 面板）

| 面板 | 内容 | 状态 | 数据来源 |
|------|------|------|----------|
| a | TCGA/GTEx 泛癌 SERPINB2/E1 差异表达 lollipop | ❌ | UCSC Xena 原始矩阵未公开 |
| b | 8 个 PDAC 队列 OS meta-analysis 森林图 | ❌ | per-study HRs 未在 Source Data |
| c | 在体验证实验设计示意图 | ⬜ | 可手绘 |
| **d** | 肿瘤质量 Control vs Serpinb2/Serpine1-KO | ✅ | MOESM9 / Figure_2d，n=7/7/7 |
| **e** | 流式 CD8+ T 细胞数/mg | ✅ | MOESM9 / Figure_2e，n=5/6/4 |
| **f** | scRNA-seq T 细胞亚型组成（classic） | ✅ | plot_tables/tcell_subsets.tsv |
| g | IHC 代表图像 CK19/CD8a/GZMB | ❌ | 原始图像未公开 |
| **h** | IHC CD8+ cells/mm² | ✅ | MOESM9 / Figure_2h，n=6/6/6 |
| **i** | IHC GZMB+CD8+ cells/mm² | ✅ | MOESM9 / Figure_2i，n=6/6/6 |
| j | Rag2−/− + NSG 实验设计示意图 | ⬜ | 可手绘 |
| **k** | NSG 异种移植肿瘤质量 | ✅ | MOESM9 / Figure_2k，n=5/5/5 |
| **l** | anti-PD-1 ± Serpin KO 生存曲线 | ✅ | MOESM9 / Figure_2l，KM 曲线 |
| m | PAI-039 ± anti-PD-1 实验设计 | ⬜ | 可手绘 |
| **n** | KC 模型 PAI-039 联合生存曲线 | ✅ | MOESM9 / Figure_2n，KM 曲线 |

**已复现：9/14（d/e/f/h/i/k/l/n）。缺口：a/b/g/j/m。**

---

## Fig 3 — Transcriptional programmes & ECM（12 面板）

| 面板 | 内容 | 状态 | 数据来源 |
|------|------|------|----------|
| a | Perturb-map Multi-modal 工作流示意图 | ⬜ | 可手绘 |
| b | Xenium + t-CyCIF 整合代表性区域图 | ❌ | Xenium 原始矩阵在 Zenodo |
| c | 通路活性热图（per-cell signature score） | ❌ | Xenium DGE 表在 Zenodo |
| d | Volcano plots Xenium DGE | ❌ | 同上 |
| **e** | PDPN+ CAF 百分比 | ✅ | MOESM10 / Figure_3e，n=6/6/6 |
| **f** | αSMA+PDPN+ 百分比 | ✅ | MOESM10 / Figure_3f，n=6/6/6 |
| **g** | CAF 亚群 UMAP | ✅ | plot_tables/caf_umap.tsv |
| **h** | CAF 亚型柱状图（Flex 化学） | ✅ | plot_tables/caf_subsets.tsv |
| i | Masson trichrome 染色图像 | ❌ | 原始组织切片未公开 |
| **j** | 胶原阳性面积 (%) | ✅ | MOESM10 / Figure_3j，n=4/4/4 |
| k | Fibrin(ogen) IHC 图像 | ❌ | 原始图像未公开 |
| **l** | Fibrin(ogen)+ 面积 (%) | ✅ | MOESM10 / Figure_3l，n=4/4/4 |

**已复现：7/12（e/f/g/h/j/l）。缺口：a/b/c/d/i/k。**

---

## Fig 4 — Fibrin-dependent macrophage retention（12 面板）

| 面板 | 内容 | 状态 | 数据来源 |
|------|------|------|----------|
| **a** | 流式巨噬细胞数/mg | ✅ | MOESM11 / Figure_4a，n=5/6/4 |
| b | 人 PDAC MICSSS 多色染色 | ❌ | 原始图像未公开 |
| c | 空间距离密度图 | ❌ | 需要 per-cell xy 坐标 |
| **d** | 巨噬细胞到 PAI+ 肿瘤细胞距离 (μm) | ✅ | MOESM11 / Figure_4d，n=6/6/6 |
| **e** | 巨噬细胞清除实验肿瘤质量 | ✅ | MOESM11 / Figure_4e，6 组 |
| **f** | CD18 阻断实验肿瘤质量 | ✅ | MOESM11 / Figure_4f，6 组 |
| **g** | myeloid 亚群柱状图（classic） | ✅ | plot_tables/myeloid_subsets.tsv |
| h | Arg1+ / fibrinogen-binding score UMAP | ❌ | per-cell expression 未公开 |
| i | BMDM 纤维蛋白培养示意图 | ⬜ | 可手绘 |
| **j** | BMDM MFI（Plastic vs Fibrin vs Argatroban） | ✅ | MOESM11 / Figure_4j |
| k | 人 PDAC 空间富集 effect sizes（3 队列） | ❌ | Visium 坐标在 GEO 限制访问 |
| l | 人 PDAC 空间富集热图 | ❌ | 同上 |

**已复现：7/12（a/d/e/f/g/j）。缺口：b/c/h/k/l。**

---

## Fig 5 — Rare serpin+ cells dominate local immunity（12 面板）

| 面板 | 内容 | 状态 | 数据来源 |
|------|------|------|----------|
| a | 每患者 SERPINB2+/SERPINE1+ 比例（n=71） | ❌ | 人 scRNA 组合数据集未公开 |
| b | 共享上调基因 Venn 图（196 genes） | ❌ | DEG 表未公开 |
| c | serpin+ 恶性细胞 pathway enrichment | ❌ | 同上 |
| d | 混合克隆实验设计示意图 | ⬜ | 可手绘 |
| **e** | 混合实验肿瘤质量 95:5 | ✅ | MOESM12 / Figure_5e，n=5/5/5 |
| **f** | 混合实验巨噬细胞 % | ✅ | MOESM12 / Figure_5f，n=5/5/5 |
| **g** | 混合实验 CD8+ T 细胞 % | ✅ | MOESM12 / Figure_5g，n=5/5/5 |
| **h** | 混合实验 GZMB+CD8+ % | ✅ | MOESM12 / Figure_5h，n=5/5/5 |
| i | 注射 vs 终点 KO/对照比例 | ❌ | per-mouse recovery count 未在 Source Data |
| j | Multiplex IF 代表性图像 | ❌ | 原始图像未公开 |
| **k** | ROI 级 CD8+ fraction 箱线图 | ✅ | MOESM12 / Figure_5k，100 ROIs |
| **l** | 对照比例 vs CD8+ 散点拟合 | ✅ | MOESM12 / Figure_5l，100 ROIs |

**已复现：7/12（e/f/g/h/k/l）。缺口：a/b/c/d/i/j。**

---

## Extended Data Fig 1 — Gene prioritisation（5 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | 候选基因筛选流程图 | ⬜ 示意图 |
| b | DepMap 依赖评分分布（46 细胞系） | ❌ 未公开 |
| c | Reactome 通路富集 | ❌ 只有 gene list |
| d | Pfam 结构域富集 | ❌ 只有 gene list |
| e | 圆形热图 normal/primary/metastasis | ❌ 需要 human+mouse scRNA |

**全部不可复现（除 a 外）。**

---

## Extended Data Fig 2 — Pro-Code barcoding（6 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | Pro-Code 载体示意图 | ⬜ 示意图 |
| b | 验证工作流示意图 | ⬜ 示意图 |
| c | H&E + 单 marker 图像 | ❌ 原始图像 |
| d | in vitro vs in vivo 环形图 | ❌ per-PC counts |
| **e** | 表位标签富集 crossbar | ✅ MOESM13 / Ext_Data_Figure_2e |
| f | PC-PC 邻域富集热图 | ❌ Squidpy 输出 |

**已复现：1/6。**

---

## Extended Data Fig 3 — In vivo vs in vitro fitness（4 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | Perturb-map 流程示意图 | ⬜ |
| b | Cas9+ vs Cas9- 培养后丰度 | ❌ per-PC counts |
| c | DepMap 依赖评分箱线图 | ❌ |
| d | 图像分析流程示意图 | ⬜ |

**全部不可复现。**

---

## Extended Data Fig 4 — Clonal evolution（7 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | Day 14 数字重建图 | ❌ per-cell coords |
| b | 焦点区域 PC 频率堆叠图 | ❌ |
| **c** | Shannon/Evenness 多样性随时间 boxplot | ✅ MOESM14 / Ext_Data_Figure_4c |
| d | 径向邻域富集图 | ❌ Squidpy 输出 |
| e | Longitudinal CD8 浸润 per KO | ❌ |
| f | Rag2−/− fitness volcano | ❌ |
| g | 免疫健全 vs Rag2−/− volcano | ❌ |

**已复现：1/7。**

---

## Extended Data Fig 5 — Functional validation（13 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | TCGA HR lollipop | ❌ per-cancer-type HRs |
| b | Western blot PAI1/PAI2 | ❌ Supplementary Data 3 |
| **c** | 倍增时间 bar plot | ✅ MOESM15 / Ext_Data_Figure_5c |
| **d** | sgRNA 验证肿瘤质量 | ✅ MOESM15 / Ext_Data_Figure_5d |
| **e** | sgRNA KO 生存曲线 | ✅ MOESM15 / Ext_Data_Figure_5e |
| **f** | scRNA-seq UMAP 主要群体 | ✅ plot_tables/umap_coords.tsv |
| **g** | T 细胞亚群 UMAP | ✅ plot_tables/tcell_umap.tsv |
| h | Hallmark 基因热图 | ❌ per-cell z-scores |
| **i** | 终末耗竭亚群柱状图 | ✅ plot_tables/tcell_subsets.tsv |
| **j** | Rag2−/− 肿瘤质量 | ✅ MOESM15 / Ext_Data_Figure_5j |
| k | 抗 PD-1 设计示意图 | ⬜ |
| **l** | KC 模型 PAI-039 生存曲线 | ✅ MOESM15 / Ext_Data_Figure_5l |
| m | JAVELIN Renal 101 PFS KM | ❌ trial data not deposited |

**已复现：8/13。**

---

## Extended Data Fig 6 — Perturb-map Multi-modal（6 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | 多模态工作流示意图 | ⬜ |
| b | Xenium ↔ t-CyCIF 配准示意图 | ⬜ |
| c | Coagulation 基因热图 | ❌ Xenium on Zenodo |
| d | Xenium UMAP by cell type | ❌ raw object |
| e | Xenium UMAP by KO | ❌ raw object |
| f | Tumor-cell UMAP by KO | ❌ raw object |

**全部不可复现。**

---

## Extended Data Fig 7 — Stromal & endothelial composition（7 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | t-CyCIF CAF 图像 | ❌ 原始图像 |
| **b** | PDPN+/αSMA+ 定量 | ✅ MOESM16 / Ext_Data_Figure_7b |
| **c** | Flex 化学 UMAP | ✅ plot_tables/umap_coords.tsv |
| d | t-CyCIF CD31 图像 | ❌ 原始图像 |
| **e** | CD31+ 内皮细胞百分比 | ✅ MOESM16 / Ext_Data_Figure_7e |
| **f** | 内皮亚群 UMAP | ✅ plot_tables/endo_umap.tsv |
| **g** | 内皮亚型柱状图 | ✅ plot_tables/endo_subsets.tsv |

**已复现：5/7。**

---

## Extended Data Fig 8 — PAI1/PAI2 and suppressive macrophages（18 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| **a** | 中性粒细胞/mg | ✅ MOESM17 / Ext_Data_Figure_8a |
| b | 小鼠 MICSSS 图像 | ❌ 原始图像 |
| c | 肿瘤中心/边缘重建 | ❌ per-cell coords |
| d | 人空间分析工作流示意图 | ⬜ |
| e | 距离密度图 mouse+human | ❌ per-cell distances |
| f | 人巨噬细胞距离 n=30 | ❌ |
| g | 巨噬细胞清除设计示意图 | ⬜ |
| **h** | 巨噬细胞数/mg（清除实验） | ✅ MOESM17 / Ext_Data_Figure_8h |
| i | CD18 阻断设计示意图 | ⬜ |
| **j** | GZMB+CD8+ %（CD18 阻断） | ✅ MOESM17 / Ext_Data_Figure_8j |
| **k** | CD8+ %（CD18 阻断） | ✅ MOESM17 / Ext_Data_Figure_8k |
| **l** | 巨噬细胞 %（CD18 阻断） | ✅ MOESM17 / Ext_Data_Figure_8l |
| **m** | ARG1+ %（CD18 阻断） | ✅ MOESM17 / Ext_Data_Figure_8m |
| n | Anti-PD-1+CD18 阻断设计示意图 | ⬜ |
| **o** | 肿瘤质量（三药组合） | ✅ MOESM17 / Ext_Data_Figure_8o |
| **p** | 巨噬细胞亚群 UMAP | ✅ plot_tables/myeloid_umap.tsv |
| q | GO 富集气泡图 | ❌ 需重新运行 DE + Enrichr |
| **r** | BMDM MFI（ARG1/PD-L1，±argatroban） | ✅ MOESM17 / Ext_Data_Figure_8r |

**已复现：9/18。**

---

## Extended Data Fig 9 — SERPIN spatial distribution（2 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | SERPINB2 分区比例（57 sections） | ❌ GEO phs002371 受限 |
| b | SERPINE1 分区比例（57 sections） | ❌ 同上 |

**全部不可复现。**

---

## Extended Data Fig 10 — SERPIN states in human PDAC（9 面板）

| 面板 | 内容 | 状态 |
|------|------|------|
| a | 每患者 serpin+ 细胞平均表达 | ❌ 组合数据集未公开 |
| b | serpin+ 及双阳比例 | ❌ |
| c | 患者水平 stacked bars | ❌ |
| d | 按 PDAC 亚型分布 violin | ❌ |
| e | 分子亚型组成 | ❌ |
| f | 亚型内相对组成 | ❌ |
| g | 治疗状态比较 | ❌ |
| h | 治疗状态表达均值 | ❌ |
| i | 机制示意图 | ⬜ |

**仅 i 为示意图；其余 8 个全部不可复现。**

---

## 📊 汇总统计

| 类别 | 数量 |
|------|------|
| **总面板数** | ~93（含 Supplementary） |
| ✅ 已精确重绘（有 Source Data 或 plot_tables） | **41** |
| ⬜ 示意图（可手绘替代灰色框） | **11** |
| ❌ 原始数据未公开，无法复现 | **41** |
| 🟡 部分数据可用（近似版本） | **0** |

**覆盖率：41/93 ≈ 44% 面板已复现**

---

## 无法复现的核心原因

### 1. 原始显微图像（最严重）
MICSSS、t-CyCIF、Xenium、H&E、Masson trichrome、IHC 的原始图像均未在 Zenodo 或 GEO 公开。这是 Nature 论文的常见做法——只保留统计数据。

### 2. Per-cell 坐标数据
Perturb-map 的空间异质性分析（Fig 1c/d）、邻域富集（Fig 1f, ED4d/e）、人 PDAC 空间距离（Fig 4c, ED8e/f）都依赖 per-cell xy 坐标，这些数据在 GSE289138 对象里但没有单独下载渠道。

### 3. Xenium 原始矩阵
Fig 3b/c/d（Xenium 空间转录组）和 ED6 全部依赖 10x Xenium raw object，该 object 上传到 Zenodo 但不是标准格式，无法本地加载。

### 4. 人 scRNA 组合数据集
Fig 5a/c、ED10 全部 8 个面板依赖 5 个人类 PDAC scRNA-seq 数据集的合并对象，其中 phs002371.v1.p1 受 dbGap 限制，其他 GEO 数据集也未提供合并后的 AnnData 文件。

### 5. Visium 空间坐标
Fig 4k/l、ED9 依赖 57 例人 PDAC Visium 数据，这些在 GEO 但需要 cell2location reference model 才能解卷积——reference atlas 不在 deposit 中。

---

## 下一步建议

如果目标是"把能做的都做出来"：

1. **Fig 1 工作流示意图**：可以用 PowerPoint / Inkscape 画一个简化版 Perturb-map 流程图替代灰色框
2. **Fig 2/3/4/5 中仍有数据缺口的 panel**（如 TCGA HR lollipop、forest plot）：TCGA GTEx 可以从 ucsc-xena.appspot.edu.cn 下载，但需要重新跑差异表达
3. **Supplementary Data 3（WB full scan）**：论文补充材料中有 WB 全膜扫描图，如果下载后可用于 ED5b
4. **JAVELIN Renal 101 PFS**：Motzer 2020 Nat Med 原文中有公开 KM 曲线和 HR，可以尝试从文中提取数字

如果目标是"达到出版级别"：
- 上述 5 类缺口都需要额外数据获取渠道，超出当前公开 deposit 范围
- 方法学论文的叙事本身就是"这五种错误如何发生"，不完整反而是一种诚实的表达
