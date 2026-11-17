# 项目约束：SERPIN–MARCO ICB 耐药轴 / 方法学纠错

本文件是本项目的常驻指令，对所有会话生效，优先于任何默认行为。
所有相对路径均相对于本文件所在目录。

**当前阶段：审计与数据排查。不产出任何新的生物学结论。**

---

## 1. 强制响应范式

在输出代码或执行任何 shell 命令之前，回复开头必须先列出四项：

1. **当前任务** —— 一句话概括。
2. **最终目标** —— 明确产出物。
3. **操作位置** —— 精确到脚本；优先修改已有文件，而不是新建。
4. **结果产物与保存位置** —— 给出路径。

---

## 2. 冻结令（最高优先级）

**GSE289138 的下游机制分析已全部完成（genotype 重跑，2026-11-14）。**

- `GSE289138_processed.h5ad` 已按 genotype 重跑，`treatment` 列已删除，
  `genotype` 列已加入。`step3_annotate_cnv.py` 和 `step4_myeloid.py` 已在
  干净的 processed object 上重跑。
- `audit_mouse_confound.py` 已生成 per-library readouts 和
  genotype-vs-chemistry effect sizes（见 `CORRECTIONS.md` C6 和
  `methods_paper/MANUSCRIPT_A_RESULTS.md` §6.6）。
- 冻结原因见 `CORRECTIONS.md` C4：该队列**没有 treatment 分组**。
  六个库是 control KO / Serpinb2 KO / Serpine1 KO × 两种化学
  （3' classic 与 Flex），每个基因型只有 1 个肿瘤，且两个化学是同一肿瘤的
  技术重复。`treatment` 列是预处理脚本里硬编码伪造的。
- 任何下游机制分析（CellChat、NicheNet、RNA velocity、CellRank、轨迹推断、
  空间邻域富集、配体-受体）仍然禁止，因为 n = 1 肿瘤/基因型，cell-level
  检验无效。

`serpin_repro/processed/QUARANTINE.md` 是该冻结令的物理标记，不要删除。

---

## 3. 拒绝"强行自洽"

- 发现组间无差异（例如 PDAC 队列 SERPIN 比例 p = 1.0），就报告无差异。
- 发现元数据共线（例如 treatment 与 condition 完全重合），就报告共线，
  并把它当作数据缺陷，而不是用模型去"校正"掉。
- **严禁**为了得到显著结果而：更换随机种子、调整 p 值阈值、逐次更换统计
  方法、事后增加分层、或改用单侧检验。
- 阴性结果与"该问题在此数据中不可回答"都是有效产出，必须写清楚。

---

## 4. 目录职责

| 目录 / 文件 | 允许做什么 |
|---|---|
| `CORRECTIONS.md` | 最高逻辑法典。新发现的数据谬误按 C5、C6… 追加，格式与 C4 一致（错误表现 → 真相 → 为何无声 → 后果 → 可自动捕捉的不变式 → 文件）。 |
| `methods_paper/` | 路径 A（方法学论文）的脚本与草稿。只放审计、量化、纠错类代码。 |
| `cohort_audit/` | 泛癌种 ICB 队列的**元数据爬取、临床维度核查、靶基因覆盖盲测**。此阶段**不做**生物学差异分析。 |
| `serpin_repro/processed/` | 预处理数据唯一存放地。`sample_metadata_CORRECTED.tsv` 是唯一元数据真理。 |
| `gse202051_hum/`、`gse300595/` | 已完成的分析；改动前先读 `CORRECTIONS.md`。 |

**禁止**在本项目任何位置生成含 `treatment`、`anti-PD-1`、`Vehicle` 等
伪造字段的衍生表。

---

## 5. 元数据断言（强制）

处理任何元数据、临床表、样本表时，必须写硬断言，不能只靠肉眼：

```python
# 1) 设计列不得与技术列共线
assert not (meta["treatment"] == meta["chemistry"]).all(), \
    "Fatal: treatment is collinear with chemistry"

# 2) 不允许凭空出现的分组列
assert "treatment" not in meta.columns, \
    "Fatal: fake treatment column detected"

# 3) 比例的分母必须是身份，不是巧合
assert mask.dtype == bool and mask.sum() < 0.5 * n_total, \
    "Fatal: mask covers most/all cells"

# 4) 请求的基因必须在面板里
assert set(genes) <= set(panel), f"missing from panel: {set(genes) - set(panel)}"
```

元数据**必须**从数据源自身的 sample 表解析（GEO series matrix 的
`!Sample_title` / `!Sample_characteristics_ch1`），不得硬编码。

---

## 6. 已知正确事实（防止再次写错）

- **GSE289138**：control KO / Serpinb2 KO / Serpine1 KO × {classic, Flex}，
  n = 1 肿瘤/基因型，**无治疗臂**。唯一可用结论是描述性存在性。
  **Genotype 重跑已完成**（2026-11-14）：`GSE289138_processed.h5ad` 带有
  `genotype` 列、无 `treatment` 列；`audit_mouse_confound.py` 已生成
  per-library readouts 和 genotype-vs-chemistry effect sizes
  （见 `CORRECTIONS.md` C6）。方向性核查（按论文所用化学分层，子集级
  结果见 `methods_paper/reproduction_report.tex`，已编译 24 页 PDF；
   报告中 13 张图全部由 R 4.6.1 + ggplot2 重绘（含 6 张 scRNA 图，
   底层表由 `audit_mouse_confound.py` 导出到
   `methods_paper/figs_reproduction/plot_tables/`），见
   `methods_paper/figs_nature.R` 与 `methods_paper/figs_nature/`）：
  T/NK exhaustion 在 classic 中两个 KO 均降低（Tex 总量 1.54×/1.75×；
  Flex 中反向）；CD8 Teff 在 NK 排除后升高 1.45×/1.46×（一致；全
  T/NK signature 因 NK 稀释而持平）；SPP1+MARCO+ 巨噬细胞 classic 中
  2.18× 降低（一致；Flex 不复现）；CAF 比例 Serpine1-KO 略降（一致）；
  CAF 亚型"基本保持"不复现；Serpine1 检出在 Serpine1-KO 库中持续存在
  （mosaic），敲除不可确认。总计：8 条可测 claim 中 7 条方向一致、
  2 条失败、1 条不可表征。
- **GSE202051**：人 PDAC snRNA，43 例（18 naive / 25 treated），224,988 核。
  SERPIN⁺ 比例 41.0% → 55.1%，p = 0.065 —— 只能说"未检出下降"。
- **GSE199102**：GeoMx DSP，22 例，235 个物理 ROI，570 段；
  面板**无 SERPINE1**。正确单位是"分区内 ROI + 患者层聚类"。
- **GSE300595**：snMultiome，12 例（5 R / 7 NR），35,115 核；
  SERPIN 比例 R vs NR p = 1.000；0 个 DAR 过 padj < 0.25。
- **队列 accession（已核查）**：IMvigor210 / Mariathasan_2018（n=348，ICB-Portal 开放）、
  Braun 2020 ccRCC（n=311 RNA-seq，ICB-Portal 有 `expression.tsv`）、GSE91061（109）、
  GSE78220（28）、PRJEB23709（73）。
- **单细胞 ICB 队列**：GSE120575（48 样本，17R/31NR，但**仅 16,291 个 CD45⁺ 免疫细胞，
  无肿瘤细胞**）、GSE123813（BCC/SCC，免疫/T 细胞分选，**无 response 列**）、
  GSE115978（7,186 细胞，含 Mal + Macrophage，**无疗效标注**）。
  **没有任何公共单细胞队列同时具备肿瘤 + 髓系 + 疗效。**
- **空间 ICB**：GSE300793（57 例 III 期黑色素瘤，辅助 anti-PD1 前取样），
  终点是 **2 年复发：33 RFS / 24 relapsed**；GEO 无 response 列，逐样本映射只在
  补充 PDF。GSE300445/446 是其 Visium / snRNA 子集。
- **ICB 下载路线**：用 **ICB-Portal**（`ngdc.cncb.ac.cn/icb/file/datasets/`，无需注册，
  29 数据集，每数据集 clinical / expression / landscape TSV + Rdata）。
  **ICBatlas 不可用**（页面与 API 均返回 405，Aliyun WAF）。
- **IMvigor210 基线预实验（已完成）**：68 R / 230 NR（50 NE 排除），
  `SERPINE1` p=0.377、`SPP1` p=0.582、`MARCO` p=0.201，BH q≥0.566；
  三者均为 **NR 更高**（Cliff's δ = −0.04…−0.10），方向一致但**统计为零**。
  不得据此宣称"bulk 表达可预测 ICB 耐药"，也不得换队列逐基因捞显著。
- **PDAC ICB 的硬缺口**：唯一带 RECIST 的队列 GSE248014 **只测肿瘤细胞**
  （5R/40NR）；唯一含髓系的 GSE197613 是分选样本且无 RECIST。
  因此"PDAC 中的 SERPIN⁺–髓系龛"在现有公共数据中不可检验。
- 详细核查见 `methods_paper/ICB_COHORT_AUDIT.md`。

---

## 7. 不得引用的旧文件

| 不得引用 | 替代 |
|---|---|
| `gse202051_hum/cross_platform_consensus.tsv` | `cross_platform_consensus_CORRECTED.tsv` |
| `gse202051_hum/processed_dsp/serpin_proximity_comparison.tsv` | `methods_paper/compartment_patient_clustered.tsv` |
| `gse202051_hum/processed/step3_summary.tsv` 的 MARCO 列 | `gse202051_hum/processed/marco_corrected_test.tsv` |
| `GSE289138_annotated.h5ad` 的 `treatment` 列 | `serpin_repro/processed/sample_metadata_CORRECTED.tsv` |

---

## 8. 统计纪律

- 推断单位必须是**患者**，不是细胞、不是段、不是 ROI。
- 空间数据先报 ICC 与设计效应，再选检验。
- 报告给定样本量下**可达的最小 p 值**（如 5 vs 7 时双侧 Mann-Whitney 下限为 0.002525）。
- 多重比较必须说明；名义 p 值一律标注为探索性。
- 跨平台一致性必须先核对**测的是不是同一个基因**。

---

## 9. 编码环境

- Python：`E:\pancreatic_cancer\submit\serpin_repro\venv\Scripts\python.exe`
- Windows / PowerShell：`python -c` 带嵌套引号会失败，**写 .py 文件再执行**。
- 不要在 `serpin_repro/venv/` 下做递归搜索或列目录（7 万+ 文件）。
- `D:\programs\research\AGENTS.md` 是指向本文件的符号链接，
  `D:\programs\research\.opencode` 是指向本目录 `.opencode` 的目录联接。
  不要在两处分别维护内容。

---

## 10. 当前主线的边界

允许：
- 方法学论文（路径 A）的补充分析与图表。
- `cohort_audit/` 中的队列可得性、样本量、疗效标签、靶基因覆盖核查。

不允许（需用户明确批准）：
- 任何新的生物学差异分析。
- 任何基于 GSE289138 的下游机制分析（CellChat、NicheNet、RNA velocity、
  CellRank、轨迹推断、空间邻域富集、配体-受体）。genotype 重跑已完成
  （2026-11-14），描述性 readouts 和化学效应量化已写入
  `CORRECTIONS.md` C6；cell-level 检验仍然禁止（n = 1 肿瘤/基因型）。
- 任何"泛癌种 SERPIN–髓系轴与 ICB 耐药"的实质性建模。

---

## 11. 操作边界：固定目录大纲，只改不建

**默认禁止新建文件。** 一切工作通过修改既有文件完成。

### 固定大纲（白名单）

```
submit/
├── AGENTS.md
├── CORRECTIONS.md
├── LIMITATIONS.md
├── RESULTS_RECONCILED.md
├── cross_platform_consensus_CORRECTED.tsv
├── methods_paper/
├── cohort_audit/
├── gse202051_hum/
├── gse300595/
└── serpin_repro/
```

不在上表中的目录一律不进入；上表中的具体文件以现有内容为准，不新增同级文件。

### 复用规则

| 要做的事 | 改哪个既有文件 |
|---|---|
| 新的队列核查 | `methods_paper/audit_icb_cohorts.py`（追加 accession，不改接口） |
| 新的核查结论 | `methods_paper/ICB_COHORT_AUDIT.md`（追加行） |
| 新发现的数据缺陷 | `CORRECTIONS.md`（追加 C5、C6…，格式同 C4） |
| 方法学论文补充分析 | `methods_paper/` 下既有脚本 |
| 队列工作区说明 | `cohort_audit/README.md` |

### 唯一例外

用户明确说"可以新建"时才允许创建新文件；否则一律复用上表。

### 为什么

C1–C4 的共同背景是"同一问题存在多套实现，结果互相矛盾"。
每多一个脚本副本，就多一个不一致的来源。只改不建是对这条教训的直接执行。
