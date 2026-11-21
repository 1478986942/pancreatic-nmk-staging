# Nature 论文学习方法 — Falcomatà et al. 2026 解剖报告

> 目标：从这篇论文提取可迁移的写作/画图/代码组织范式，用于你自己的生信论文

---

## 一、整体结构：漏斗 + 闭环

### 1.1 标题公式

```
[分子名]–[细胞类型] axis in [疾病] heterogeneity and [表型]
```

**解析**：
- `serpin–myeloid axis` = 分子轴（两个关键词）
- `pancreatic cancer heterogeneity` = 疾病背景 + 核心现象
- `immune evasion` = 最终生物学意义

**你的仿写模板**：
```
[X基因]–[Y细胞] axis in [Z癌症] progression and [治疗响应/耐药]
```

---

### 1.2 Abstract 三段式结构（严格对应 Results 章节）

| Abstract段落 | 内容功能 | 字数比例 | 包含要素 |
|-------------|---------|---------|---------|
| **第1段** | 背景 + 知识空白 + 本文要解决什么问题 | 25% | 疾病重要性、当前局限、核心科学问题 |
| **第2段** | 方法 + 主要发现（按 Results 顺序） | 55% | 技术平台→筛选结果→因果验证→机制→人体验证 |
| **第3段** | 概念升华 + 临床意义 | 20% | 新原则/新认知 + 治疗启示 |

**原文逐句拆解**：

```
Paragraph 1 (背景):
  Sentence 1: PDAC is characterized by X, yet Y.        ← 已知事实+矛盾
  Sentence 2: This raises the question of A or B?       ← 提出核心问题
  → 两句定调，不废话

Paragraph 2 (方法+发现):
  Sentence 1: Here, using X platform, we determined Y.  ← 方法+目标
  Sentence 2: This analysis revealed Z.                 ← 关键发现1（筛选）
  Sentence 3: We identify A and B as C.                 ← 关键发现2（锁定因子）
  Sentence 4: These factors promote D via E mechanism.  ← 机制链条
  Sentence 5: Loss of A/B or pharmacological inhibition improves F and synergizes with G.  ← 治疗意义
  → 五句话走完完整证据链，每句一个层次

Paragraph 3 (升华):
  Sentence 1: These findings identify X as local spatial organizers...    ← 新概念
  Sentence 2: Linking tumour-intrinsic heterogeneity to...                ← 连接两个领域
  → 一句话拔高，建立跨领域桥梁
```

**你的仿写模板**：
```
[ disease ] is characterized by [ feature ], yet [ paradox ].
This raises the question of [ scientific question ]?

Here, using [ technology ], we determined [ goal ].
This analysis revealed [ key finding 1 ].
We identify [ gene/protein ] as [ role ].
These [ molecules ] promote [ phenotype ] via [ mechanism ].
[ Intervention ] improves [ outcome ] and synergizes with [ existing therapy ].

These findings identify [ concept ] as [ new principle ],
linking [ domain A ] to [ domain B ].
```

---

### 1.3 Results 小节的排列逻辑

原文13个小节，**每个小节回答一个问题，且下一节的标题由上一节的结论自然引出**：

```
Fig 1: Who are the regulators?        → 筛选出 serpin 家族
    ↓ "These serpins shape the TME"
Fig 2: What do they do to tumours?   → KO 后肿瘤生长受抑，T细胞浸润增加
    ↓ "This effect requires adaptive immunity"
Fig 3: How do they reprogramme cells? → 转录程序改变（ECM/纤溶通路下调）
    ↓ "The ECM change involves fibrin"
Fig 4: Which cell mediates this?      → 巨噬细胞是效应执行者
    ↓ "This happens in human tumours too"
Fig 5: Is it conserved in humans?     → 人源数据验证 + 稀有细胞主导原则
```

**核心规则**：
1. 每个小节标题 = 一个完整的问题+答案陈述句
2. 相邻小节之间用过渡句连接（在段落开头或结尾）
3. 不设独立"Methods"小节在Results里——方法分散在各小节的方法描述中
4. 一张主图（Fig）通常对应1-2个小节，Extended Data补充细节

---

## 二、Figure 设计规范（可直接复用）

### 2.1 配色方案（Nature 标准）

```r
# 永远不要自己发明颜色，用这套
C_CTRL  <- "#4d4d4d"   # 灰色：对照组
C_B     <- "#1f77b4"   # 蓝色：第一个KO/处理组
C_E     <- "#d62728"   # 红色：第二个KO/处理组
C_MIX   <- "#f2a541"   # 橙色：混合/联合处理

# 衍生色（变浅用于相同基因的不同条件）
lighten <- function(col, amt=0.45) {
  v <- col2rgb(col)/255
  v <- v + (1-v) * amt
  rgb(v[1], v[2], v[3])
}
```

**为什么这四种颜色**：
- 灰/蓝/红是 ColorBrewer Set1 的前三个，被 ggplot2 默认采用，读者最熟悉
- 橙色列为第四组，与前三者区分度足够
- 所有图严格遵守同一套配色，读者无需重新学习图例

---

### 2.2 每类图的固定模板

#### 模板A：柱状图 + 个体散点（最常用）

```r
nature_bar <- function(groups, ylab, letter, pairs=NULL, claimed=NULL) {
  ggplot(mn, aes(group, mean, fill=group)) +
    geom_col(width=0.55, colour="black", linewidth=0.28) +     # 柱子：黑框55%宽
    scale_fill_manual(values=colors, guide="none") +
    geom_point(position=jitter(width=0.085),                   # 点：白色填充黑色边框
               shape=21, fill="white", colour="black", 
               stroke=0.28, size=1.2) +
    annotate("text", x=0.5, y=max_y*1.12,                     # 左上角字母标号
             label=letter, hjust=0, vjust=1, 
             size=2.9, fontface="bold") +
    labs(y=ylab, title=letter) +
    theme_minimal(base_size=14) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.y=element_line(colour="#ededed"),   # 浅灰横线
          axis.line=element_line(linewidth=0.24),              # 黑色坐标轴
          axis.text.x=element_text(angle=32, hjust=1))         # x标签旋转32°
}
```

**关键参数记忆**：
- `width=0.55`：柱子不要太胖，留白给散点
- `stroke=0.28`：散点边框粗细，太粗会盖住柱子
- `angle=32`：x轴标签不要90度（难看），32度是最佳角度
- `colour="#ededed"`：网格线颜色要淡，不能抢眼

---

#### 模板B：Kaplan-Meier 生存曲线

```r
km_panel <- function(surv_data, n_per_group, letter, cols, ltypes) {
  ggplot(km_df, aes(x=days, y=survival, colour=group, linetype=group)) +
    geom_step(linewidth=0.5) +                                  # 阶梯线 0.5pt
    scale_colour_manual(values=cols) +
    scale_linetype_manual(values=ltypes) +                     # solid/dashed 区分
    scale_y_continuous(breaks=seq(0,100,25), limits=c(0,101)) + # y轴25%一格
    annotate("text", x=max_x*0.7, y=0.3,                       # 图内标注 log-rank p
             label=paste0("p =", p_value), size=3, 
             colour="grey30", fontface="italic") +
    labs(x="Days", y="Survival (%)", title=letter) +
    theme(legend.position="bottom",                             # 图例放底部
          legend.key.size=unit(0.26,"cm"),
          plot.margin=margin(6,3,2,3))                          # 上方留更多空间给图例
}
```

**技巧**：
- 用 `linetype`（实线/虚线）同时编码两组信息：治疗 vs 对照 + 基因型
- 图例放底部而不是右侧，因为生存曲线需要横向空间展示时间轴
- p值直接标在图内，不在图注——方便读者快速捕捉显著性

---

#### 模板C：UMAP/tSNE 降维图

```r
umap_panel <- function(df, col, letter, pal, pt=0.12) {
  ggplot(df, aes(UMAP1, UMAP2, colour=.data[[col]])) +
    geom_point(size=pt, alpha=0.8) +                            # 点要小(0.12)且半透明
    scale_colour_manual(values=pal) +
    theme_void(base_size=14) +                                  # 完全去掉背景
    theme(legend.position="bottom",                            # 图例放底部
          legend.text=element_text(size=leg_cex),
          plot.title=element_text(size=TI, face="bold", hjust=0),
          plot.margin=margin(1,2,1,1))                         # 边距最小化
}
```

**要点**：
- `size=0.12` + `alpha=0.8`：点太多时会重叠，半透明能看到密度
- `theme_void()`：零网格零坐标轴，只保留点和颜色
- 如果点超过1万个，先 `sample(n=10000)` 再画——否则PDF会爆

---

#### 模板D：热图（通路活性/表达矩阵）

```r
heatmap_panel <- function(mat, row_anno, col_anno, letter) {
  pheatmap::pheatmap(mat,
    cluster_rows=TRUE, cluster_cols=TRUE,                    # 默认聚类
    annotation_row=row_anno, annotation_col=col_anno,        # 行/列注释条
    color=colorRampPalette(c("#becbe6","#ffffff","#f07d7d"))(100),  # Nature蓝白红
    fontsize_row=5, fontsize_col=6,                           # 小字体
    cutree_rows=2, cutree_cols=3,                             # 可选：强制分群数
    main=letter)
}
```

**配色选择**：
- 负→正：`blue → white → red`（经典，适用于log2FC）
- 低→高：`lightgray → blue → orange → red`（适用于百分比/丰度）
- **永远不要用 RdBu（红-蓝 diverging）**，打印后灰色阶分不清

---

### 2.3 Panel 布局规则

```
单张 Figure 的布局黄金法则：

┌──────────┬──────────┬──────────┐
│    a     │    b     │    c     │  ← 第一行：3个小组件
├──────────┼──────────┴──────────┤
│    d     │         e            │  ← 第二行：左小右大
├──────────┴──────────────────────┤
│              f                  │  ← 第三行：大图占满
└─────────────────────────────────┘

原则：
1. 从左到右、从上到下阅读（符合文字习惯）
2. 相关的panel紧挨着（如 d 和 e 讲同一个实验的不同方面）
3. 大图放下面（读者的视线最终落在最重要的结果上）
4. schematic 永远放第一个位置（a 面板）
```

---

### 2.4 图注（Figure Legend）写法

**原文 Fig 2 图注结构**：
```
Fig. 2 | Loss of Serpinb2 and Serpine1 promotes T cell-mediated anti-tumour immunity.

a, Differential expression... [1-2句说明数据和方法]

b, Meta-analysis of... [n=1,176 patients from 8 PDAC cohorts].

c, Experimental design for in vivo validation... [1句说清楚实验设置]

d, Tumour mass... (n = 7 mice per group). [n值必须标注]

e, Intratumoural CD8+ T cells... [统计检验方法: Two-sided unpaired Welch's t-tests]

f, Relative abundance of T cell subsets... [亚群定义说明]

g, Representative pseudocoloured immunohistochemistry... Scale bars, 125 µm. [比例尺]

h,i, Quantification of CD8+ (h) and GZMB+CD8+ (i) cells... [n值 + 统计方法]

j, Experimental design for analysis in Rag2−/− and NSG mice.

k, Tumour mass... (n = 5 mice per group).

l, Kaplan–Meier survival... [log-rank test, n值]

m, Experimental design for evaluation of PAI-039...

n, Kaplan–Meier survival analysis... [四组n值列出]

d,e,h,i,k, Two-sided unpaired Welch's t-tests. 
l,n, Two-sided log-rank (Mantel–Cox) tests.
Dots represent individual mice; bars indicate mean ± s.d.
```

**图注写作模板**：
```
Fig. X | [一句话总结这张图的核心结论]

[a] [数据描述] ([统计检验], n = X). [关键数字].

[b] [方法简述]. [n值]. 

[c] [实验设计]. [比例尺].

[d] [结果]. [统计学显著性标记, 如 *p < 0.05].

[e,f] [子图关系说明]. 

[g] [全局说明: 误差线含义, 数据表示方式].
```

---

## 三、Methods 写作技巧

### 3.1 方法部分的三层结构

原文 Methods 分五个大块，每块内部都是**目的→材料→操作→分析**的顺序：

```
Selection and prioritization of extracellular tumour-derived genes
  → 目的：找到调控TME的分泌因子
  → 材料：human/mouse scRNA-seq, DepMap, Compartments database
  → 操作：DE分析 → 过滤细胞外蛋白 → 取DepMap低频依赖基因
  → 分析：Reactome/Pfam富集 → 文献验证

Lentiviral vector construction and production
  → 目的：构建PC/CRISPR文库
  → 材料：Brie CRISPR library oligos, PC/CRISPR lentiviral vector (Addgene)
  → 操作：退火→连接→转化→筛选→测序验证
  → 分析：Sanger sequencing confirmation

Cell culture
  → 目的：培养KPC/KC细胞系
  → 材料：DMEM + 10% FBS, 100 U/ml penicillin/streptomycin
  → 操作：37°C, 5% CO2, 传代≤30次
  → 分析：Mycoplasma testing negative

Mouse experiments
  → 目的：原位移植模型
  → 材料：C57BL/6J, Rag2−/−, NSG mice (Jackson Laboratory)
  → 操作：50,000 cells orthotopically implanted into pancreas
  → 分析：ARRIVE guidelines, IACUC approved

Flow cytometry / CyTOF / MICSSS / CyCIF / Xenium
  → 每个技术单独一段
  → 标准结构：样本制备 → 抗体/探针 → 仪器 → 分析软件 → 定量方法
```

### 3.2 Methods 写作原则

**原则1：让他人能复现，但不啰嗦**
```
❌ 不好："Cells were cultured in Dulbecco's Modified Eagle Medium 
         (Thermo Fisher Scientific, cat. no. 11965092) supplemented 
         with 10% fetal bovine serum (Sigma-Aldrich, F2442) and 
         100 U/ml penicillin/streptomycin (Thermo Fisher Scientific, 
         15140122) at 37°C in a humidified incubator with 5% CO2."

✅ 好："Cells were maintained in DMEM supplemented with 10% FBS and 
       1% penicillin/streptomycin at 37°C, 5% CO2."
```

**原则2：商业试剂给出厂家+货号**
```
✅ "anti-mouse PD-1 (CD279) antibody (clone 29F.1A12; BioXcell, BE0273)"
```

**原则3：统计方法放在对应图注，_methods_里只放总体原则**
```
Methods: "Statistical analyses were performed using GraphPad Prism 10 
          or R v4.2.2. Details are provided in figure legends."

Figure legend: "Two-sided unpaired Welch's t-tests were used..."
```

---

## 四、Discussion 六段模板

```
Paragraph 1: Summary of key findings
  "Our study identifies X as a dominant regulator of Y in Z context."
  → 一句话概括全篇核心发现，不用重复细节

Paragraph 2: Mechanistic insight
  "Rather than acting through X, our data show that Y operates via Z."
  → 解释"怎么做到的"，提出新的机制理解

Paragraph 3: Comparison with prior work
  "Although previous studies reported X, our findings reveal Y..."
  → 对话前人工作：支持/拓展/修正

Paragraph 4: Broader implications
  "These results have implications for X field, suggesting Y..."
  → 拉到更广阔的领域谈意义

Paragraph 5: Core principle提炼
  "Notably, our work establishes a principle that X..."
  → 这句话是全文最高光时刻，要让读者记住

Paragraph 6: Limitations and future directions
  "Several limitations should be considered... Future studies should..."
  → 主动承认不足比审稿人挑出来强一百倍
```

---

## 五、代码组织模式（供你仿效）

### 5.1 单层扁平目录结构

```
project/
├── AGENTS.md              # 项目宪法（置顶，优先级最高）
├── CORRECTIONS.md         # 错误记录（每次修bug都追加）
├── LIMITATIONS.md         # 限制声明（永远诚实）
├── RESULTS_RECONCILED.md  # 校正后的结论
│
├── cohort_A/              # 每个队列一个目录
│   ├── step1_qc.py
│   ├── step2_analysis.py
│   └── processed/
│       ├── final_table.tsv
│       └── fig_main.pdf
│
├── cohort_B/
│   └── ...
│
└── methods_paper/         # 论文相关
    ├── manuscript.md
    ├── figures.R           # 单一绘图入口
    └── figs_output/
```

### 5.2 脚本命名规范

```
step1_inspect.py      # 第一步：探索性检查
step2_analysis.py     # 第二步：主分析
step3_stats.py        # 第三步：统计分析
step4_visualize.py    # 第四步：作图

diag_*.py             # 诊断脚本（非主线，只跑一次）
fix_*.py              # 修复脚本（纠正之前发现的错误）
audit_*.py            # 审计脚本（交叉验证）
```

### 5.3 关键函数封装模式

```python
# 1. 路径管理：所有路径集中在顶部常量
BASE = r"E:\path\to\data"
OUT  = os.path.join(BASE, "processed")
FIG  = os.path.join(BASE, "figures")

# 2. 数据加载：统一的读取函数
def load_h5ad(path):
    adata = sc.read_h5ad(path)
    assert "patient_id" in adata.obs.columns, \
        f"Fatal: missing patient_id in {path}"
    return adata

# 3. 断言先行：每个关键步骤前检查
def compute_fraction(mask, denominator):
    assert mask.dtype == bool, "mask must be boolean"
    assert mask.sum() > 0, "mask has no True values"
    assert denominator > mask.sum(), \
        f"denominator ({denominator}) < mask sum ({mask.sum()})"
    return mask.sum() / denominator

# 4. 日志输出：每步有明确标记
print(f"[STEP 3] Computing patient-level metrics...", flush=True)
# ... analysis ...
print(f"  Done. Output: {out_path}", flush=True)
```

---

## 六、自检清单（发paper前必过）

### 图表层面
- [ ] 每张图有一个明确的"一句话结论"写在标题或图注第一句
- [ ] 所有数值型axis都有单位（%, cells/mg, days...）
- [ ] 每组n值明确写出（"n = 5 mice per group"）
- [ ] 统计检验方法在图注或Methods中说明
- [ ] p值标注格式统一（*p<0.05, **p<0.01, ***p<0.001 或在图注写明实际p值）
- [ ] 比例尺标在图像上（显微照片必须有）
- [ ] 颜色对色盲友好（避免红绿直接对比，用红蓝或添加纹理）

### 写作层面
- [ ] Abstract 三段各自回答了：Why? What? So what?
- [ ] 每个 Results 小节标题是一个完整陈述句
- [ ] Discussion 没有简单重复 Results
- [ ] 局限性单独成段，主动交代而非被动暴露
- [ ] 所有缩写首次出现时全称+缩写并列
- [ ] 图注中的 "Data are representative of..." 只在真正代表时才用

### 代码层面
- [ ] 每个 .py 脚本顶部有清晰的功能说明 docstring
- [ ] 硬编码的路径改用常量
- [ ] 关键中间结果保存为 .tsv/.csv（不只是 .h5ad）
- [ ] 至少在一个地方写了 assert 检查
- [ ] print() 语句有 `[STEP X]` 标记，便于调试

---

## 七、立即可以动手练的三个练习

### 练习1：重写 Fig 2d 的图注
用原文 Fig 2d 的数据（Control: 7只鼠，Serpinb2-KO: 7只，Serpine1-KO: 7只；均值±SD；Welch's t-test p=0.001和p=0.0003），写一段符合 Nature 规范的图注。

### 练习2：设计你研究的 Figure 1
假设你在做一个"某基因在肺癌免疫治疗响应中的作用"研究，Figure 1 应该包含哪几个panel？每个panel回答什么问题？画出草图（文字描述即可）。

### 练习3：写一段 Discussion 第一段
模仿原文 Discussion 第一段的语气和结构，把你最近一个分析的结论写成一段话（3-5句）。

---

需要我把某个部分展开详细讲解，或者带你逐个完成这三个练习吗？