# ICB cohort availability audit — metadata only

Date: 2026-11-14
Purpose: decide whether the pan-cancer "SERPIN⁺ tumour / myeloid co-occurrence"
hypothesis can be started at all, before writing any analysis code.

**No expression data were downloaded and no analysis was run.** Everything below
comes from GEO series-matrix headers (via `!Sample_characteristics_ch1`),
E-utilities summaries, and published data-availability statements.

---

## Tier definitions

The hypothesis needs three things simultaneously: tumour cells, myeloid cells,
and an ICB response label. Cohorts are graded by which of these they actually
provide.

| tier | what it supports | requirement |
|---|---|---|
| **A** | cell-level co-occurrence within a sample | whole-tumour scRNA/snRNA containing tumour **and** myeloid, plus response |
| **B** | sample-level co-expression / co-infiltration | bulk transcriptome with response, n large |
| **C** | neither | sorted populations, tumour-only, or no response label |

Physical proximity (the "niche" in the spatial sense) needs spatial data and is
listed separately.

---

## Tier B — bulk, response-labelled, adequately powered

| cohort | cancer | n | therapy | response field | access |
|---|---|---|---|---|---|
| IMvigor210 | urothelial | 348 | atezolizumab | CR/PR/SD/PD | open — ICB-Portal `Mariathasan_2018_{clinical,expression}.tsv`; also R package `IMvigor210CoreBiologies` |
| Braun 2020 | ccRCC | 311 RNA-seq | nivolumab (CM-009/010/025) | ORR / Benefit / PFS / OS | open — ICB-Portal `Braun_2020_{clinical,expression}.tsv`; WES is EGA-restricted but not needed |
| GSE91061 (Riaz) | melanoma | 109 (51 pre, 58 on) | nivolumab | PD 48 / SD 34 / PR+CR 23 | open |
| GSE78220 (Hugo) | melanoma | 28 | pembrolizumab | PD 13 / PR 10 / CR 5 | open |
| PRJEB23709 (Gide) | melanoma | 73 | anti-PD-1 ± anti-CTLA-4 | response | open (ENA) |
| ICB-Portal (benchmark) | 9 cancers | >1,400, 29 datasets | mixed | harmonised | **open, no registration** — `ngdc.cncb.ac.cn/icb/file/datasets/<name>_{clinical,expression,landscape}.tsv` + `.Rdata`; bulk zips under `/icb/file/` |
| ICBatlas | 9 cancers | 1,515, 25 studies | mixed | harmonised | **not usable** — own page templates return HTTP 405 and the JSON API is blocked by an Aliyun WAF challenge, verified in a real browser |

**Usable n without any data-access negotiation: ~800 samples across melanoma,
RCC and urothelial, all with a response label.**

## Tier A — single-cell with tumour + myeloid

| cohort | cancer | n | compartments present | response | access |
|---|---|---|---|---|---|
| GSE115978 (Jerby-Arnon) | melanoma | 7,186 cells, 32 tumours | Mal 2,018, Macrophage 420, T, CAF, Endo | **none** (naive vs post only) | open |
| GSE120575 (Sade-Feldman) | melanoma | 48 samples, 16,291 cells | **CD45⁺ immune cells only — no tumour cells** | 17 R / 31 NR | open; GEO suppl is TPM + an empty SRA template, no cell-type annotation |
| GSE123813 (Yost) | BCC/SCC | 86 | **CD45⁺/CD3⁺ immune (T-cell) sorted — no tumour cells** | not in GEO | no series matrix (404 by design); suppl metadata has `cell.id, patient, treatment, sort, cluster, UMAP1, UMAP2`, no response column |
| GSE159115 | RCC | 14 | tumour + benign | **none** | open |

Only GSE115978 demonstrably contains both malignant cells and macrophages — and
it has no response label. Every response-labelled single-cell cohort found so
far is CD45⁺ immune-sorted (GSE120575, GSE123813). **No dataset pairs
whole-tumour single-cell composition with an ICB response endpoint at usable n.**

## PDAC ICB — all Tier C for the niche question

| cohort | n | therapy | response | what is actually sequenced | why it fails |
|---|---|---|---|---|---|
| GSE248014 | 45 samples / 25 patients | entinostat + nivolumab | **RECIST: 5 R / 40 NR** | `cell type: tumor cells` for all 45 | **tumour cells only** — no myeloid, no T |
| GSE197613 | 91 | neoadjuvant anti-PD-1 | arm A/B/Control, OS >2y/<2y | FACS-sorted CD11b / CD19 / CD4 / CD8 | sorted populations, no RECIST |
| GSE313101 | 37 | FOLFIRINOX ± nivolumab | arm only | bulk, pre/post | no response label |
| GSE279781 | 30 | CD137 agonist + anti-PD-1 | arm A/B/C | TILs only | no tumour, no myeloid |
| GSE277116 | 28 | neoadjuvant immunotherapy | path response 8 R / 20 NA | tumour / lymph node | 8 labelled responses |

This is the decisive result for the PDAC arm: **the one PDAC ICB dataset with
RECIST labels contains only tumour cells, and the one with myeloid cells is
sorted and has no RECIST.** The niche hypothesis cannot be tested in PDAC ICB
data that exists publicly today.

## Spatial ICB

| cohort | cancer | n | ICB | note |
|---|---|---|---|---|
| GSE300793 (+ GSE300445/446) | melanoma | 57 + 4 + 2 | adjuvant anti-PD-1 | **endpoint = 2-year relapse status: 33 relapse-free / 24 relapsed** (stage III, pre-treatment). GEO has no response field; the per-sample mapping is only in the 39.8 MB supplement PDF (JITC 2026) |
| GSE272362 | PDAC | 30 | none | spatial, no ICB |
| GSE335452 | PDAC | 15 | unclear | scRNA + spatial, CD8 exhaustion |
| GSE277116 | PDAC | 28 | neoadjuvant | TLS focus, 8 responses |

**No PDAC spatial dataset with ICB response exists publicly.** Physical-proximity
claims are therefore possible only in melanoma, at small n.

---

## Verdict

1. **The hypothesis is startable, but only at the sample level.** With ~800
   response-labelled bulk samples across three cancers, we can test whether a
   SERPIN⁺-tumour / myeloid co-infiltration axis associates with ICB resistance,
   and whether PDAC sits at the saturated end of that axis.
2. **Cell-level co-occurrence is not available at all.** GSE115978 has tumour
   and myeloid cells but no response; GSE120575 and GSE123813 have response but
   are immune-sorted and contain no tumour cells.
3. **Physical proximity is not testable outside melanoma.** No PDAC spatial ICB
   data exists.
4. **PDAC ICB cannot carry the niche claim.** The RECIST-labelled cohort is
   tumour-only.

## Consequence for the claim and the target journal

The defensible claim is:

> A SERPIN⁺-tumour / myeloid co-infiltration axis is associated with ICB
> resistance across melanoma, RCC and urothelial carcinoma; PDAC sits at the
> constitutive extreme of this axis, which is why the axis has no dynamic range
> there and does not stratify response.

Not: "SERPIN⁺ cells organise a physical exclusion niche conserved across
cancers." The spatial part would need melanoma-only evidence.

Novelty risk: SERPINE1/PAI-1 association with ICB resistance is already
reported in melanoma and RCC. A single-gene association will not clear a 10+
reviewer. The two defensible novelties are:

- the **co-infiltration structure** (SERPIN + SPP1/MARCO together) rather than
  SERPINE1 alone;
- the **PDAC exception** — we already have the null (GSE300595, p = 1.0), and it
  becomes a prediction of the model rather than a failure.

Realistic ceiling if the pan-cancer association is clean and the PDAC exception
holds: **8–12** (JITC / CCR / Advanced Science). If the association is weak:
5–7. Reaching 15+ would need spatial validation in a second cancer type or
functional data, neither of which is available.

## Open items — resolved

All four gating items are closed. None reverses the verdict; two of them
strengthen it.

### 1. GSE120575 has no tumour cells — Tier A loses its only response-labelled candidate

The GEO summary states the study profiled **16,291 immune cells** from 48
melanoma samples. The supplementary files are a TPM matrix plus
`GSE120575_patient_ID_single_cells.txt.gz`, which is an **unfilled SRA metadata
template** (35 columns, every row a `#` comment). There is no cell-type
annotation in GEO. The cohort is CD45⁺ immune-only; it cannot supply the tumour
half of a tumour–myeloid co-occurrence measure.

Consequence: **no public single-cell cohort pairs whole-tumour composition with
an ICB response label.** Tier A is empty for this hypothesis.

### 2. GSE123813 has no series matrix by design, and is immune-sorted

The 404 is not a broken accession: the series publishes supplementary files
only. `GSE123813_bcc_all_metadata.txt.gz` (1.38 MB) has columns
`cell.id, patient, treatment, sort, cluster, UMAP1, UMAP2`, with
`sort = CD45+ CD3+`. The `scc_` and `bcc_tcell_` metadata files are T-cell
subsets. **No response column exists in any GEO file**; response lives in the
Yost et al. 2019 supplement.

Consequence: same as item 1 — immune-only, no co-occurrence, and the response
label needs an external table.

### 3. ICBatlas is unreachable; ICB-Portal is the working bulk route

ICBatlas (`guolab.wchscu.cn/ICBatlas/`) returns HTTP 405 for its own Angular
templates (`/static/ICBatlas/pages/*.html`) and for
`/ICBatlas/api/document_list`, including when called from a real browser; the
API sits behind an Aliyun WAF `acw_sc__v2` challenge. It is not usable
programmatically or interactively.

ICB-Portal (`ngdc.cncb.ac.cn/icb/resources`) works with no registration and no
login. It exposes:

| resource | URL pattern |
|---|---|
| bulk bundle | `https://ngdc.cncb.ac.cn/icb/file/ICB_data.zip` (also `marker.zip`, `Results.zip`, `figures.zip`, `src.zip`) |
| per dataset | `https://ngdc.cncb.ac.cn/icb/file/datasets/<name>_clinical.tsv` |
| | `..._expression.tsv`, `..._landscape.tsv`, `..._.Rdata` |

29 datasets across three pages: Braun_2020; Gide_PRE/COMBINE/MONO/Paired_R/
Paired_NR; Hugo_PRE_2016; Jung_2019; Kim_2018; Lee_PRE_2020; Lee_Paired_R/NR;
Liu_2019; Liu_EXPOSURE/NAIVE; **Mariathasan_2018 (IMvigor210)**; MGH_Paired_R/NR;
MGH_PRE_2021; Miao_2018; Motzer_2020; Nathanson_PRE_2017; Riaz_PRE/EXPOSURE/
NAIVE/Paired_R/Paired_NR; Snyder_2017; VanAllen_2015.

Verified downloads:

- `Braun_2020_clinical.tsv` — 125 columns, includes `ORR`, `Benefit`, `PFS`,
  `OS`, `Tumor_Shrinkage`; 311 patients.
- `Mariathasan_2018_clinical.tsv` — 31 columns, includes
  `Best Confirmed Overall Response`, `binaryResponse`, `Resp_NoResp`, `os`,
  `censOS`, `Treatment`, `Pre_On`; IMvigor210.

Consequence: **the Tier B ammunition is one download, not a data-access
negotiation.** This is the route to use.

### 4. GSE300793 has a response endpoint, but only in the paper supplement

The series matrix carries only `tissue: Cutaneous melanoma`. The associated
paper (JITC 2026, PMID 42315253) defines the endpoint explicitly: 57 primary
cutaneous melanomas, stage III, sampled before adjuvant anti-PD-1,
**stratified by 2-year relapse status — 33 relapse-free, 24 relapsed**. GEO has
no response column. The only supplementary file is
`GSE300793_all_expression_data_withlast.rds` (40 MB), a compressed S4/list
object with no readable metadata strings; `pyreadr` cannot decode it and this
machine has no R. The per-sample relapse mapping is inside the 39.8 MB
`jitc-14-6-s001.pdf`, which PMC serves only behind a proof-of-work download
wall.

Consequence: the cohort is usable in principle, but **its response labels
require manual extraction from the PDF**. It is not a turnkey spatial
validation set.

### What this does to the verdict

Unchanged in direction, sharper in form:

- Tier B (~800 labelled bulk samples) is fully open and is the only realistic
  route.
- Tier A is empty: every response-labelled single-cell cohort is immune-sorted.
- The melanoma spatial cohort has labels, but not in machine-readable form.
- PDAC remains untestable for the niche question.

The next step is therefore still the bounded IMvigor210 pilot — now with
`Mariathasan_2018_clinical.tsv` and `Mariathasan_2018_expression.tsv` from
ICB-Portal, no R package required.

---

## IMvigor210 baseline pilot — pre-specified, single cohort

Script: `methods_paper/audit_icb_cohorts.py --audit-imvigor`
Data: ICB-Portal `Mariathasan_2018_{clinical,expression}.tsv` +
`Mariathasan_2018.Rdata` (gene names only exist in the Rdata — see
`CORRECTIONS.md` C5).

**Design, fixed before the result was seen**

- one cohort: IMvigor210 / Mariathasan 2018, urothelial, atezolizumab;
- baseline only: `Pre_On == "Pre"` (all 348 samples; asserted, not assumed);
- response: RECIST CR/PR (R) vs SD/PD (NR); 50 NE excluded → **68 R / 230 NR**;
- primary family: `SERPINE1`, `SPP1`, `MARCO`, two-sided Mann-Whitney U,
  Benjamini-Hochberg across the three;
- effect size: Cliff's delta (+1 = higher in R, −1 = higher in NR);
- secondary: mean z-score of the three genes, nominal only, exploratory;
- descriptive: pairwise Spearman rho among the three genes.

**Primary results**

| gene | median R | median NR | mean R | mean NR | U | p | q (BH) | Cliff's δ |
|---|---|---|---|---|---|---|---|---|
| SERPINE1 | 202.53 | 238.45 | 389.07 | 372.96 | 7268 | 0.377 | 0.566 | −0.071 |
| SPP1 | 233.97 | 269.58 | 671.64 | 807.69 | 7476 | 0.582 | 0.582 | −0.044 |
| MARCO | 1.73 | 2.87 | 7.75 | 9.82 | 7021 | 0.201 | 0.566 | −0.102 |

**Secondary (exploratory).** Mean z-score of the three genes: median R
−0.319 vs NR −0.255, U = 7035, p = 0.209, δ = −0.100.

**Descriptive co-variation** (labelled n = 298): SERPINE1–SPP1 ρ = +0.402,
SERPINE1–MARCO ρ = +0.434, SPP1–MARCO ρ = +0.532.

**Verdict: directionally consistent, statistically null.**

All three genes and the composite are higher in non-responders — the predicted
direction — but none reaches nominal significance, and BH q is ≥ 0.57. The
effect sizes are small (|δ| = 0.04–0.10). With 68 vs 230 the cohort can resolve
moderate effects (|δ| ≈ 0.2) but not small ones; the observed effects are
therefore compatible with a weak or absent bulk association, and this dataset
cannot distinguish those two.

Two details that must not be smoothed over:

- **SERPINE1's mean is higher in responders** (389.07 vs 372.96) while its
  median and Cliff's delta are higher in non-responders. The mean is tail-driven;
  the rank statistic was the pre-specified one. Do not report the mean as a
  direction.
- The three genes **do co-vary moderately** (ρ = 0.40–0.53), so the
  "co-expression module" premise is not contradicted — but co-variation is not
  association with response.

**What this does to the project**

- It does **not** support "SERPIN⁺/MARCO⁺ bulk expression predicts ICB
  resistance" as a stand-alone claim.
- It does **not** refute it either. Bulk TPM is the wrong instrument for a
  spatial co-occurrence hypothesis, and this is one cancer type.
- For the methods paper it is a useful negative: an axis that is visually
  obvious in single-cell and spatial data does not collapse into a bulk
  single-gene difference. That is precisely the "bulk cannot substitute for
  spatial co-occurrence" point.

**Do not** escalate this into a multi-cohort gene-by-gene fishing exercise. If
the axis is tested further, it needs a pre-specified composite in an
independent cohort with the same endpoint.
