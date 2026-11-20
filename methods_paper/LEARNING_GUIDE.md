# 论文学习方法论 — 以 Falcomatà et al. Nature 2026 为例

## 一、论文结构拆解（漏斗式叙事）

### 1.1 标题与摘要的关系

**标题**：`A serpin–myeloid axis in pancreatic cancer heterogeneity and immune evasion`
- 关键词：serpin（机制分子）、myeloid（效应细胞）、heterogeneity（现象）、immune evasion（表型）
- 四个词组成完整的因果链：**基因 → 细胞 → 现象 → 功能**

**摘要三段式结构**：
```
第1段：背景 + 核心问题
  "PDAC是免疫抑制的微环境，但肿瘤存在异质性——问题是：免疫抵抗是全局属性还是空间受限的niche？"

第2段：方法 + 主要发现
  "用Perturb-map做空间功能基因组学筛选→锁定serpin家族→证明SERPINE1/SERPINB2是主调控因子"

第3段：机制 + 意义
  "通过fibrin ECN基质 niche 保留并编程巨噬细胞→排斥CD8+ T细胞→联合抗PD-1增效"
```

**学习点**：摘要必须回答三个问题——**为什么做？做了什么？意味着什么？** 每个问题一句话。

---

### 1.2 Results 小节设计原则

原文13个Results小节，每个小节对应一个明确的**科学问题**：

| 小节 | 回答的问题 | 核心图表 |
|------|-----------|---------|
| 1. Tumour-derived TME regulators | **谁**在调控TME？ | ED Fig 1a-e（筛选流程） |
| 2. Extracellular factors shape clonal fate | **哪个基因**影响克隆命运？ | Fig 1a-h（Perturb-map筛选） |
| 3. Niche modelling precedes clone selection | 是niche先变还是克隆先变？ | Fig 1f, ED Fig 4（时间顺序） |
| 4. Serpinb2 and Serpine1 mediate immune evasion | 这些基因**做什么**？ | Fig 2a-i（表型验证） |
| 5. ...promote CD8 T cell exclusion | 具体怎么排斥T细胞？ | Fig 2d-i（流式/IHC） |
| 6. ...drive anti-PD-1 resistance | 能不能联合ICB？ | Fig 2j-l,n（生存曲线） |
| 7. ...promote suppressive programmes | 转录层面发生了什么？ | Fig 3c-d（通路热图） |
| 8. ...regulate ECM remodelling | 结构层面发生了什么？ | Fig 3e-f,i-l（ECM定量） |
| 9. TAMs mediate serpin immunosuppression | 哪个细胞执行抑制？ | Fig 4a-f（巨噬细胞） |
| 10. Serpins shape TAM states via fibrin niches | 怎么编程巨噬细胞？ | Fig 4g-j（极化实验） |
| 11. Serpin niches contain suppressive TAMs | 人肿瘤里有这个现象吗？ | Fig 4k,l（空间验证） |
| 12. SERPINB2/E1 define human PDAC states | 人的肿瘤细胞长什么样？ | Fig 5a-c（scRNA亚群） |
| 13. Serpin-expressing cells dominate local immunity | 少量阳性细胞够吗？ | Fig 5d-l（混合实验） |

**学习点**：每个小节 = 一个**独立的科学问题**，下一节的标题由上一节的结论自然引出。没有一节是"凑数"的。

---

### 1.3 Discussion 六段法

```
段落1：机制总结
  "我们发现了什么新东西"——基因驱动的局部免疫邻域重塑先于克隆选择

段落2：概念升华（ECM作为指导性组织者）
  "这改变了什么认知"——ECM不是物理屏障，是指令性支架

段落3：与已有研究的关系
  "这和别人的工作是什么关系"—— Stromal PAI1 vs tumor-derived PAI1/PAI2

段落4：更广的背景（凝血与免疫驱动）
  "这在更大范围里意味着什么"——PDAC凝血激活高发，PAI1可能是免疫驱动的选择压力

段落5：核心原则提炼（稀有细胞主导）
  "最关键的 Principle 是什么"——少量serpin+细胞足以驱动TME重塑

段落6：局限性与未来方向
  "还有什么不知道"——筛选限于curated因子、未直接测量纤溶活性等
```

**学习点**：Discussion 不是结果复述，而是**概念升华 + 文献对话 + 原则提炼 + 坦诚局限**。

---

## 二、Figure 设计规范

### 2.1 Nature 风格的核心特征

| 特征 | 具体表现 | 实现方法 |
|------|---------|---------|
| **灰底白网格** | `panel.grid.major = element_line(colour="#ededed", linewidth=0.18)` | theme_minimal + 自定义 |
| **黑色坐标轴** | `axis.line = element_line(linewidth=0.24)` | 显式绘制轴线 |
| **字体统一** | 全部用 Helvetica/Arial 风格的无衬线体 | R 默认已满足 |
| **panel letter** | 左上角 bold 小写字母 (a, b, c...) | `annotate("text", x=0.025, y=0.965, label="a")` |
| **p值标注** | 图内或图注中标注，不省略 | `geom_text(aes(label=sprintf("p=%.3g", p)))` |
| **n值透明** | 每个panel说明样本量 | legend 或 caption |
| **颜色规范** | Ctrl=灰色, KO1=蓝色, KO2=红色, Mix=橙色 | 固定调色板 `#4d4d4d #1f77b4 #d62728 #f2a541` |

### 2.2 每类图的适用场景

```r
# ── 柱状图+散点叠加（单个变量多组比较）──
# 适用：Fig 2d/e/h/i/k, Fig 3e/f/j/l, Fig 4a/d/e/f/j, Fig 5e/f/g/h
nature_bar <- function(groups, ylab, letter, pairs, claimed=NULL) {
  ggplot(mn, aes(group, mean, fill=group)) +
    geom_col(width=0.55, colour="black", linewidth=0.28) +      # 柱子：均值
    geom_point(position=jitter, shape=21, fill="white",          # 点：个体值
               colour="black", stroke=0.28, size=1.2) +
    stat_summary(fun=mean, geom="errorbar",                     # 误差线
                 width=0.2, linewidth=0.3)
}

# ── Kaplan-Meier 生存曲线 ──
# 适用：Fig 2l/n, ED 5e/l
km_panel <- function(surv_data, n_per_group, letter, cols, ltypes) {
  ggplot(km_df, aes(x=days, y=survival, colour=group, linetype=group)) +
    geom_step(linewidth=0.5) +                                    # 阶梯曲线
    geom_segment(aes(xend=days, yend=0.02),                      # risk table tick marks
                 colour="grey40", linewidth=0.3)
}

# ── UMAP/tSNE 降维可视化 ──
# 适用：Fig 3g, ED 2e, ED 5f/g, ED 7c/f, ED 8p
umap_panel <- function(df, col, letter, pal, pt=0.12) {
  ggplot(df, aes(UMAP1, UMAP2, colour=eval(col))) +
    geom_point(size=pt, alpha=0.8) +                              # 点要小且半透明
    scale_colour_manual(values=pal) +
    theme_void()                                                  # 去掉网格
}

# ── 分组柱状图（两个维度：亚型 × 基因型）──
# 适用：Fig 2f, Fig 4g, ED 5i
subset_bars <- function(df, subsets, letter, ylab) {
  ggplot(long_df, aes(subset, value, fill=genotype)) +
    geom_col(position=position_dodge(0.78), width=0.72,           # dodge 不要1.0
             colour="black", linewidth=0.24)
}

# ── 热图（矩阵数据）──
# 适用：Fig 1f(原版), Fig 3c(原版), Fig 4l(原版)
heatmap_panel <- function(mat, row_anno, col_anno) {
  pheatmap::pheatmap(mat,
    cluster_rows=TRUE, cluster_cols=TRUE,                        # 默认聚类
    annotation_row=row_anno, annotation_col=col_anno,            # 行/列注释
    color=colorRampPalette(c("blue","white","red"))(100),        # Nature蓝白红
    fontsize_row=6, fontsize_col=7)                               # 小字体
}

# ── Volcano plot ──
# 适用：Fig 3d(原版), ED 4f/g(原版)
volcano_plot <- function(de_results, sig_thresh=0.05, fct_thresh=1) {
  ggplot(de_results, aes(log2FC, -log10(pval))) +
    geom_point(alpha=0.4, size=0.8) +                             # 所有基因淡色
    geom_point(data=sig_genes, colour="red", size=1.2) +          # 显著基因红色
    geom_hline(yintercept=-log10(sig_thresh), linetype="dashed")+
    geom_vline(xintercept=c(-fct_thresh, fct_thresh), linetype="dashed")
}
```

### 2.3 Figure panel 排列规则

```
Nature 标准布局（7.2英寸画布）：

┌────────────────────────────────────────────┐
│ a   │  b   │  c   │                         │  ← 一行最多3-4个panel
├─────┼──────┼──────┤                         │
│  d   │  e   │  f   │                         │
├─────┴──────┴──────┤                         │
│         g           │  ← 占满整行的大图       │
└─────────────────────┘                        │
                                             │
Figure Legend:                                │
a Workflow schematic. b-d Quantification...   │  ← 从左到右，从上到下读
```

**关键规则**：
1. 字母顺序 = 阅读顺序（先左后右，先上后下）
2. 大panel放下面，小组件放上面
3. 相关panel放在一起（如 Fig 2d+i 都是IHC定量，紧挨着）
4.  schematics/图例放第一行第一个位置

---

## 三、代码组织模式

### 3.1 本项目的代码结构（单层扁平式）

```
submit/
├── AGENTS.md                    # 项目宪法（最高优先级指令）
├── CORRECTIONS.md               # 错误清单 C1-C6
├── LIMITATIONS.md               # 限制声明
├── RESULTS_RECONCILED.md        # 最终结论
│
├── methods_paper/               # 方法学论文主体
│   ├── MANUSCRIPT_A_RESULTS.md  # 论文草稿
│   ├── figs_nature.R            # 绘图脚本（单一入口）
│   ├── audit_icb_cohorts.py     # ICB队列审计
│   ├── audit_mouse_confound.py  # 小鼠混杂审计
│   └── figs_nature/*.pdf        # 输出目录
│
├── gse202051_hum/              # GSE202051 分析
│   ├── step1_inspect.py         # 步骤1: 检查
│   ├── step2_analysis.py        # 步骤2: 分析
│   ├── step3_stats.py           # 步骤3: 统计
│   └── processed/               # 中间产物
│
├── gse300595/                  # GSE300595 分析
│   ├── step1_qc.py              # QC
│   ├── step2_rna.py             # RNA分析
│   ├── step5_atac.py            # ATAC分析
│   └── step6_bridge.py          # DSP桥接
│
└── serpin_repro/               # GSE289138 重建
    ├── step2_cluster_markers.py
    ├── step3_annotate_cnv.py
    ├── step4_myeloid.py
    └── processed/
        ├── sample_metadata_CORRECTED.tsv   # 唯一元数据真理
        └── QUARANTINE.md                   # 冻结标记
```

**学习点**：每个目录 = 一个独立分析流水线，step1→step2→step3 串行执行，输出进 `processed/`。

### 3.2 Python 脚本最佳实践（本项目遵循的）

```python
# ✅ 好的做法
BASE = r"E:\pancreatic_cancer\submit\gse202051_hum"   # 用常量定义路径
OUT = os.path.join(BASE, "processed")
os.makedirs(OUT, exist_ok=True)

adata = sc.read_h5ad(os.path.join(DATA, "GSE202051_total.h5ad"))

# 硬断言（关键！）
assert "treatment" not in adata.obs.columns, "Fatal: fake treatment column detected"
assert mask.dtype == bool and mask.sum() < 0.5 * n_total, \
    "Fatal: mask covers most/all cells"

# 每一步写日志
print(f"[STEP 3] Computing patient-level metrics...", flush=True)

# ❌ 避免的做法
df = pd.read_csv("/home/user/data/GSE202051/sample_metadata.tsv")  # 硬编码绝对路径
result = some_function(df)  # 没有断言验证
# 没有日志输出
```

### 3.3 R 绘图脚本模式

```r
# figs_nature.R 的结构
# ── Part 1: Setup ──────────────────────────────────
suppressPackageStartupMessages({ library(ggplot2); library(gridExtra) })
ROOT   <- "E:/pancreatic_cancer/submit"
OUT    <- file.path(ROOT, "methods_paper", "figs_nature")
COL_CTRL <- "#4d4d4d"; COL_B <- "#1f77b4"; COL_E <- "#d62728"

# ── Part 2: Helper functions ───────────────────────
nature_bar <- function(...) { ... }
km_panel <- function(...) { ... }
ph <- function(letter, reason) { ... }   # placeholder for missing data

# ── Part 3: Data loading ───────────────────────────
g2d <- sd_bar("fig2", "Figure_2d")
g2e <- sd_bar("fig2", "Figure_2e")
# ...

# ── Part 4: Build each figure ──────────────────────
fig1 <- list(ph("a", "..."), ph("b", "..."), ...)
compose(fig1, "fig1.pdf", ncol=4)

fig2 <- list(ph("a","..."), nature_bar(g2d,...), ...)
compose(fig2, "fig2.pdf", height=5*1.24+0.35)

# ── Part 5: Run all ────────────────────────────────
message("All PDFs written to: ", OUT)
```

**学习点**：
1. **Helpers 集中在顶部**——所有通用函数定义在一起
2. **数据来源清晰**——`sd_bar()` / `rd()` 函数封装读取逻辑
3. **Figure 组装为列表**——`list(panel_a, panel_b, ...)` 然后 `compose()`
4. **Placeholder 函数**——`ph()` 统一处理无法复现的面板

---

## 四、写作叙述技巧

### 4.1 如何用一段话引出整个Results小节

**原文 Fig 1 开头**：
> "Owing to its central role in PDAC progression, we sought to identify genes involved in orchestrating the PDAC TME."

**公式**：`[临床/生物学重要性] + [知识空白] + [本研究目标]`

**仿写练习**（假设你的研究）：
> "Pancreatic cancer exhibits profound spatial heterogeneity in immune infiltration, yet the cellular programs driving T-cell exclusion remain poorly defined."

### 4.2 如何在图注中报告统计

**原文 Fig 2d 图注**：
> "Tumour mass of control, Serpinb2-KO, and Serpine1-KO orthotopic tumours (n = 7 mice per group). Two-sided unpaired Welch's t-tests."

**模板**：
```
[测量内容] of [条件] ([n = X per group]). 
[统计检验类型]; [显著性标记]. 
Data are mean ± s.d.
```

### 4.3 Discussion 的过渡句

| 目的 | 模板句式 |
|------|---------|
| 连接结果与前文 | "Consistent with our Perturb-map findings," |
| 提出新解释 | "These findings suggest that..." |
| 对比前人工作 | "Whereas previous studies reported X, we found Y," |
| 指出局限 | "Several limitations should be considered." |
| 展望未来 | "Together, our findings establish a principle linking..." |

---

## 五、自测清单（学完这篇你应该能做到的）

### 图表能力
- [ ] 能用 ggplot2 画出 Nature 风格的柱状图+散点（nature_bar）
- [ ] 能画 Kaplan-Meier 生存曲线（km_panel）
- [ ] 能画分组柱状图（position_dodge）
- [ ] 能画 UMAP 降维图（theme_void + 小点）
- [ ] 能创建灰底白网格的 theme
- [ ] 能在图上加字母标号（annotate text）
- [ ] 能在图注里报告正确的统计方法

### 代码能力
- [ ] 能写一个批量读取多个 Source Data Excel 的函数
- [ ] 能从 GEO 下载并解析 series matrix 文件
- [ ] 能写 assert 断言检查元数据的合理性
- [ ] 能用 RMarkdown 或 Quarto 生成可复现的报告
- [ ] 能把一整套图组合成 multi-panel PDF

### 写作能力
- [ ] 能写出三段落摘要（背景→方法→意义）
- [ ] 能给每个Results小节起一个明确回答问题的标题
- [ ] 能写出包含统计信息的图注
- [ ] 能在Discussion里坦诚地列出局限性

---

## 六、下一步学习建议

### 短期（本周）
1. 把这个项目的 `figs_nature.R` 当教材逐行读，理解每个 helper function
2. 用自己的数据画一张相似的柱状图+散点图
3. 写一篇模拟摘要（模仿原文三段式）

### 中期（本月）
1. 找一个公开 scRNA-seq 数据集，完整走一遍 QC→聚类→marker→作图的流程
2. 学习使用 `compartment_patient_clustered.py` 中的患者聚类推断方法
3. 模仿本文的 CORRECTIONS.md 格式，记录你自己分析中发现的错误

### 长期
1. 尝试投稿 Briefings in Bioinformatics 的方法学论文
2. 建立自己的 "错误清单" 文档习惯（这是本项目最大的方法论贡献）
3. 学习如何写 reproducible research 的 README
