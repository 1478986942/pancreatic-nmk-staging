# cohort_audit — 泛癌种 ICB 队列核查工作区

**范围：只做元数据。不做生物学差异分析。**

本目录回答三个问题，除此之外不写任何分析代码：

1. 该队列是否有 **R/NR（或 CR/PR/SD/PD）** 随访标签？
2. 核心基因 **SERPINE1 / SPP1 / MARCO** 在该队列的测序深度/面板中是否可测？
3. 该队列是否同时包含**肿瘤细胞与髓系细胞**（决定能否算共现结构）？

## 已完成

`../methods_paper/ICB_COHORT_AUDIT.md` —— 首轮核查结论：

- 泛癌种 bulk 约 800 例带疗效标注，可直接使用。
- PDAC ICB 存在硬缺口：唯一带 RECIST 的队列只测肿瘤细胞。
- 空间 ICB 仅黑色素瘤，且 n 小。

## 已清（4/4）

| 项 | 结论 | 证据 |
|---|---|---|
| GSE120575 是否捕获肿瘤细胞 | **否**。16,291 个 CD45⁺ 免疫细胞；GEO 无细胞注释（suppl 里的 `patient_ID` 文件是空 SRA 模板） | GEO summary + suppl 文件检查 |
| GSE123813 series matrix 404 | **不是坏链接，是本来就不提供 matrix**。suppl 只有免疫/T 细胞元数据，`sort = CD45+ CD3+`，无 response 列 | `GSE123813_*_metadata.txt.gz` |
| ICBatlas / ICB-Portal 下载入口 | **ICBatlas 不可用**（自身模板与 API 均返回 405，Aliyun WAF）；**ICB-Portal 完全开放** | 浏览器实测 |
| GSE300793 response | **有终点，但不在 GEO**：2 年复发状态，33 RFS / 24 relapsed；逐样本映射在 39.8 MB 补充 PDF | JITC 2026, PMID 42315253 |

## 门户

| 门户 | URL | 注册 | 格式 | 结论 |
|---|---|---|---|---|
| ICB-Portal | `https://ngdc.cncb.ac.cn/icb/resources` | 不需要 | 每数据集 `clinical` / `expression` / `landscape` 的 `.tsv` + `.Rdata`；整包 `ICB_data.zip` | **首选**。29 数据集，含 IMvigor210（`Mariathasan_2018`）、Braun_2020、Riaz、Hugo、Gide、Motzer |
| ICBatlas | `https://guolab.wchscu.cn/ICBatlas/` | — | — | **不可用**。页面模板 405，API 被 Aliyun WAF 拦截 |
| IMvigor210 原 R 包 | GitHub 维护版 | 不需要 | R 包 | 可用，但 ICB-Portal 的 TSV 更直接 |
| Braun 2020 补充表 4 | Nature Medicine 2020 | 不需要 | 补充表 | 可用，但 ICB-Portal 已提供 `expression.tsv` |

已实测可下载：

- `Mariathasan_2018_clinical.tsv`：31 列，含 `Best Confirmed Overall Response`、
  `binaryResponse`、`Resp_NoResp`、`os`、`censOS`、`Treatment`、`Pre_On`。
- `Braun_2020_clinical.tsv`：125 列，含 `ORR`、`Benefit`、`PFS`、`OS`、
  `Tumor_Shrinkage`；311 例。

## 命名约定

- 脚本：`audit_<cohort>.py`，只读元数据，不下载表达矩阵。
- 输出：`<cohort>_metadata.tsv`，必须含 `patient_id`、`response`、
  `assay`、`n_samples` 四列，缺列即报错。
- 任何脚本开头必须有：
  ```python
  assert "treatment" not in meta.columns or meta.attrs.get("source") == "geo"
  ```
