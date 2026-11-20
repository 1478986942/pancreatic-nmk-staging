# Learning Exercises — Practice with Real Code

## Exercise 1: Recreate a Nature-style Bar Plot

用以下模拟数据，画一张类似 Fig 2d 的柱状图：

```
Control:   [0.57, 0.60, 0.48, 0.62, 0.52, 0.37, 0.24]
Serpinb2:  [0.11, 0.29, 0.28, 0.35, 0.25, 0.08, 0.15]
Serpine1:  [0.09, 0.10, 0.09, 0.13, 0.08, 0.13, 0.10]
```

要求：
- 灰色柱（Control）、蓝色柱（Serpinb2-KO）、红色柱（Serpine1-KO）
- 白色圆点表示每个小鼠的值（带 jitter）
- 黑色边框
- x轴标签旋转32度
- y轴有标题 "Tumour mass (g)"
- 左上角标注 "d"
- 底部显示 p 值比较

## Exercise 2: Write a Hard Assert

下面这段代码有 bug，找出所有缺少断言的地方并加上：

```python
import pandas as pd

# Load metadata
meta = pd.read_csv("sample_metadata.tsv", sep="\t")

# Filter for treated patients
treated = meta[meta["status"] == "treated"]

# Compute SERPIN fraction
serpin_frac = treated["serpin_positive"] / treated["total_malignant"]

# Save result
out = serpin_frac.mean()
print(f"Mean SERPIN fraction in treated: {out}")
```

提示：至少加 3 个 assert，分别检查：
1. 处理列是否存在
2. 分母不为零
3. 过滤后的样本量合理

## Exercise 3: Design a Results Section Title

原文 Results 小节标题的规则：
- 陈述句，回答一个明确的问题
- 不用冒号，不用问句
- 动词通常是：identify / reveal / show / demonstrate / promote / drive

练习：把你的研究结果改成符合这个风格的标题。

你的原始想法：
- "ICB治疗前后PDAC肿瘤微环境的单细胞变化"
- " SERPIN在不同癌症类型中的表达差异分析"
- "基于空间转录组的肿瘤免疫微环境研究"

改成 Nature 风格：
1. _________________________________________________
2. _________________________________________________
3. _________________________________________________

## Exercise 4: Write an Abstract Paragraph

模仿原文第一段，写一段关于"胰腺癌免疫治疗耐药"的背景介绍（2-3句话）。

原文模板：
> "Pancreatic ductal carcinoma (PDAC) is characterized by a highly immunosuppressive, extracellular matrix-rich microenvironment, yet tumours display marked heterogeneity. This raises the question of whether immune resistance is a global tumour property or is organized within spatially restricted niches."

你的版本：
> "[疾病] is characterized by [关键特征], yet [矛盾/空白]. This raises the question of [核心科学问题]."

## Exercise 5: Create a Placeholder Function

原文有一个 `ph()` 函数用来生成灰色占位框。自己写一个类似的，但要支持：
- 面板字母（"a", "b"...）
- 原因文本（可以换行）
- 可选的颜色主题（灰色/黄色/红色表示不同严重程度）

---

## 答案参考（Exercise 1 代码框架）

```r
library(ggplot2)

# Data
ctrl <- c(0.57, 0.60, 0.48, 0.62, 0.52, 0.37, 0.24)
b2   <- c(0.11, 0.29, 0.28, 0.35, 0.25, 0.08, 0.15)
e1   <- c(0.09, 0.10, 0.09, 0.13, 0.08, 0.13, 0.10)

groups <- list(Control=ctrl, `Serpinb2 KO`=b2, `Serpine1 KO`=e1)
names(groups)

# Build mean dataframe
mn <- data.frame(
  group = factor(names(groups), levels=names(groups)),
  mean  = sapply(groups, mean)
)

# Build points dataframe
pts <- do.call(rbind, lapply(names(groups), function(g) {
  data.frame(group=g, value=groups[[g]])
}))
pts$group <- factor(pts$group, levels=names(groups))

# Plot
ggplot(mn, aes(group, mean, fill=group)) +
  geom_col(width=0.55, colour="black", linewidth=0.28) +
  scale_fill_manual(values=c("#4d4d4d","#1f77b4","#d62728"), guide="none") +
  geom_point(data=pts, aes(group, value), inherit.aes=FALSE,
             shape=21, fill="white", colour="black", stroke=0.28, size=1.2,
             position=position_jitter(width=0.085, height=0, seed=1)) +
  annotate("text", x=0.5, y=max(unlist(groups))*1.15, 
           label="d", hjust=0, vjust=1, size=3.5, fontface="bold") +
  labs(x=NULL, y="Tumour mass (g)") +
  theme_minimal(base_size=14) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="#ededed", linewidth=0.15),
        axis.text.x=element_text(size=11, angle=32, hjust=1, vjust=1),
        axis.text.y=element_text(size=11),
        axis.title.y=element_text(size=11),
        plot.margin=margin(10,10,10,10))
```
