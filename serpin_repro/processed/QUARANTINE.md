# QUARANTINE — GSE289138

**状态：genotype 重跑已完成，冻结解除（2026-11-14）。**

`GSE289138_processed.h5ad` 已按 genotype 重跑，`treatment` 列已删除，
`genotype` 列已加入。`step3_annotate_cnv.py` 和 `step4_myeloid.py` 已在
干净的 processed object 上重跑。`audit_mouse_confound.py` 已生成
per-library readouts 和 genotype-vs-chemistry effect sizes（见
`CORRECTIONS.md` C6 和 `methods_paper/MANUSCRIPT_A_RESULTS.md` §6.6）。

注意：`GSE289138_annotated.h5ad` 仍带有旧的 `treatment` 列（来自缓存），
但它已重跑过一次，其 `anno`/`malignant`/`leiden_epi` 列与新的 processed
object 一致。审计脚本通过 barcode join 从 annotated 转移这些列，因此
不影响结果。

## 唯一可用的元数据

`sample_metadata_CORRECTED.tsv`（同目录）：

| filename | tumor_id | chemistry | genotype |
|---|---|---|---|
| GSM8784900_CF111 | CF111 | classic | control_KO |
| GSM8784901_CF112 | CF112 | classic | Serpinb2_KO |
| GSM8784902_CF113 | CF113 | classic | Serpine1_KO |
| GSM9575309_CF111_Flex | CF111 | Flex | control_KO |
| GSM9575310_CF112_Flex | CF112 | Flex | Serpinb2_KO |
| GSM9575311_CF113_Flex | CF113 | Flex | Serpine1_KO |

n = 1 肿瘤/基因型；两个化学是同一肿瘤的技术重复。

## 解除条件（已全部满足，2026-11-14）

1. ✅ 按 genotype 重跑预处理，删除 `treatment` 列，加入 `genotype` 列。
2. ✅ 重跑后确认 `assert "treatment" not in adata.obs`（在
   `pdac_sc_preprocessing.py` 的 concat 步骤中已实现）。
3. ✅ 用户明确解除冻结（2026-11-14 指令"解冻，按 genotype 重跑"）。

在此之前，GSE289138 只能用于描述性陈述（SERPINE1⁺ 恶性细胞与
SPP1/MARCO 巨噬细胞在 KPC 肿瘤中存在），不能用于任何组间推断。

现在 GSE289138 可用于：
- 描述性 genotype 比较（n = 1 肿瘤/基因型，2 个化学技术重复）
- 化学效应量化（classic vs Flex 的 readout 差异）
- 方向性检查（与原文 claim 的对比）
- 不得用于：任何 cell-level p-value、任何治疗效应声明
