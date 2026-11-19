# 开源项目状态报告

**日期**: 2026-11-19  
**提交**: `b32f462 Establish reproducible SERPIN-MARCO audit package`

---

## 📦 仓库位置

### Gitee（主仓库）
```
https://gitee.com/msdyy/projects_information-biology
```
- 分支: `master`
- 状态: ✅ 已推送同步

### GitHub（镜像）
```
https://github.com/1478986942/pancreatic-nmk-staging
```
- 分支: `master`
- 状态: ⚠️ 已覆盖（原内容为 N-M-K Framework）

---

## ✅ 已完成的工作

| 项目 | 详情 |
|------|------|
| **代码仓库建立** | 首提 commit，175 文件，35.4 MB |
| **Git 规范化** | 全局用户配置、LF 换行、提交规范 |
| **数据过滤** | .gitignore 排除 h5ad/tsv/venv 等大文件 |
| **README 文档** | 项目定位、科学结论、目录结构 |
| **审查清单** | CORRECTIONS.md（C1-C6）、LIMITATIONS.md |
| **论文草稿** | MANUSCRIPT_A_RESULTS.md（完整方法学论文） |
| **复核报告** | RESULTS_RECONCILED.md、reproduction_report.pdf |
| **队列审计** | ICB_COHORT_AUDIT.md（公开数据集可得性分析） |
| **图表重绘** | methods_paper/figs_nature/*.pdf（13 张图） |
| **代码检验** | 41 个 .py 文件 py_compile 通过 |

---

## ❌ 未完成事项

| 任务 | 优先级 | 说明 |
|------|--------|------|
| **选定目标期刊** | 🔴 高 | Briefings in Bioinformatics / Computers in Biology and Medicine |
| **投稿格式调整** | 🔴 高 | elsarticle.cls、插图规格、参考文献格式 |
| **Supplementary 整理** | 🟡 中 | 打包所有 .tsv/.pdf 为独立 supplement |
| **添加 License** | 🟢 低 | 建议 MIT 或 CC-BY 4.0 |
| **注册 DOI** | 🟢 低 | Zenodo 或 figshare |
| **预印本发布** | 🟡 中 | bioRxiv / BMC Bioinformatics |

---

## 🎯 投稿策略建议

### 核心卖点
- **"Five silent failures"** — 五种错误都产生合理数字且无异常，属于系统性风险
- **具体不变量设计** — 每个错误都有对应的硬断言可提前捕获
- **负结果即结果** — IMvigor210 阴性对照证明 SERPIN-髓系轴在 bulk RNA-seq 不可测

### 最佳投递期刊
1. **Briefings in Bioinformatics** (IF ≈ 8.5)
   - 最匹配：方法学 + 数据批判性分析
   - 开放审核，审稿周期约 6-8 周
   
2. **Genome Biology** (IF ≈ 10)
   - 更大影响力但需要更强叙事包装

3. **Nature Methods** (IF ≈ 40)
   - 不适合当前体量，除非扩展为多实验室案例集

### 文章定位
这不是"生物学发现"论文，而是**方法论教训论文**：
> "当公共组学数据被重新分析时，五种不同的静默失败如何改变结论"

---

## 📂 当前仓库结构

```
submit/
├── AGENTS.md                  # 项目约束与冻结令
├── CORRECTIONS.md             # C1-C6 错误清单（含修复方案）
├── LIMITATIONS.md             # 17 条限制声明
├── RESULTS_RECONCILED.md      # 校正后结论汇总表
├── README.md                  # 项目总览
├── methods_paper/
│   ├── MANUSCRIPT_A_RESULTS.md    # 完整论文草稿
│   ├── ICB_COHORT_AUDIT.md        # ICB 队列审计
│   ├── reproduction_report.pdf    # GSE289138 复现报告
│   └── figs_nature/               # 重绘的 13 张 Nature 图
├── gse202051_hum/
│   └── processed/                 # GSE202051 修正后表
├── gse300595/
│   └── processed/                 # GSE300595 患者级摘要
└── serpin_repro/
    ├── step3_annotate_cnv.py      # CNV 注释脚本
    ├── step4_myeloid.py           # 髓系分析脚本
    └── processed/                 # 修正后元数据
```

---

## 🚀 立即可执行的操作

1. **复制仓库** → `git clone https://gitee.com/msdyy/projects_information-biology.git`
2. **阅读论文草稿** → `methods_paper/MANUSCRIPT_A_RESULTS.md`
3. **查看错误清单** → `CORRECTIONS.md`
4. **核对限制** → `LIMITATIONS.md`

如需进一步整理投稿材料，请告知我具体方向。
