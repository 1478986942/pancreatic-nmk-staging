# Reconciled results

Date: 2026-11-14

**Read `CORRECTIONS.md` and `LIMITATIONS.md` first.** Three defects (DSP
pseudo-replication, DSP SERPINB2-only axis, GSE202051 MARCO denominator) changed
two headline claims. This file states only what the corrected analyses support.

---

## 1. The question

Does neoadjuvant chemotherapy preserve or eliminate the SERPIN⁺ tumour /
MARCO⁺ macrophage niche in human PDAC, and is there a myeloid chromatin basis
for the niche that could explain preservation?

## 2. Bottom line

1. **The SERPIN⁺ malignant programme persists after chemotherapy**, and if
   anything is more common in treated patients — but the difference is **not
   significant** (41.0 % → 55.1 % of malignant nuclei, p = 0.065, n = 18 vs 25).
2. **Preservation of the MARCO⁺ macrophage niche is NOT supported.** With the
   corrected denominator, the MARCO⁺ fraction of macrophages is 8.2 % in naive
   vs 4.7 % in treated patients (p = 0.76) — a non-significant *decrease*, not
   an increase. The earlier "preserved" reading came from dividing by all nuclei
   while macrophage density was rising.
3. **The spatial T-cell claim is neither confirmed as published nor refuted.** The
   original p = 0.0063 is invalid because it pools three tissue compartments with
   different baselines. The correctly stratified, patient-clustered analysis gives
   a **weaker but directionally consistent** result: SERPIN-high *epithelial*
   regions have +2.15 T-cell/NK signal (permutation p = 0.034, signed-rank
   p = 0.048), carried mainly by untreated patients. There is no evidence of
   T-cell exclusion, but the effect is borderline and compartment-specific.
4. **No epigenomic specificity can be demonstrated.** Zero DARs survive FDR, and
   the responder-enriched and non-responder-enriched accessibility programmes
   are statistically indistinguishable in the spatial data (ρ = 0.88–0.90).
5. **The one comparison that reaches p < 0.05 is myeloid abundance, not the
   SERPIN axis**: responders have more myeloid nuclei than non-responders
   (13.9 % vs 3.9 %, p = 0.048, δ = +0.714, n = 5 vs 7). This is exploratory and
   on the edge of what n = 12 can resolve.

Nothing in this project survives multiple-testing correction. Treat all of the
above as directional, hypothesis-generating observations.

---

## 3. GSE202051 — snRNA, 224,988 nuclei, 43 patients

Definitions mirror the source study: `SERPIN_score` = mean of `SERPINE1` +
`SERPINB2`; median split **within malignant epithelium**.

### SERPIN⁺ malignant programme

| group | n patients | mean SERPIN⁺ fraction of malignant nuclei |
|---|---|---|
| naive | 18 | 41.0 % |
| treated | 25 | 55.1 % |

Mann-Whitney p = 0.065; Cliff's δ = −0.338 (computed naive vs treated, so
negative = naive lower). Direction is *increase* after chemotherapy; magnitude
is not significant.

### MARCO⁺ macrophage compartment (corrected)

Denominator is the dataset's own `Level 3 == "Macrophage"` label (10,140
nuclei), not all nuclei.

| metric | naive | treated | p |
|---|---|---|---|
| MARCO⁺ / macrophage | 8.20 % | 4.71 % | 0.758 |
| MARCO⁺ / myeloid | 7.15 % | 3.93 % | 0.912 |
| macrophage fraction of nuclei | 3.54 % | 6.30 % | 0.143 |
| (old, incorrect) MARCO⁺ per 100 nuclei | 0.291 | 0.472 | 0.990 |

MARCO is detected in 764 / 224,988 nuclei (0.34 %), so this compartment is at
the edge of measurability. The most that can be said: **macrophage density rises
modestly and the MARCO⁺ fraction falls modestly after chemotherapy; neither is
significant.**

### Interpretation

The SERPIN⁺ tumour-cell state is not eliminated by chemotherapy, which is the
one part of the original premise that holds — weakly. The MARCO⁺ macrophage
half of the premise does not hold once measured correctly.

---

## 4. GSE199102 — GeoMx DSP, 194 physical ROIs (570 segments), 22 patients

**Axis: `SERPINB2` only.** `SERPINE1` is absent from the panel.

### The T-cell claim: three units, three answers

The original analysis pooled all 570 segments across Immune / Epithelial / CAF and
compared SERPIN-high vs SERPIN-low. That pools three tissues with different
baselines and is not a valid contrast. Aggregating everything to 194 ROIs is not
valid either, because it averages those same three baselines together.

| unit of analysis | T-cell p |
|---|---|
| pool 570 segments (original) | 0.0063 |
| average to 194 physical ROIs | 0.685 |
| **compartment-stratified, patient-clustered** | **0.034 (permutation)** |

Correct analysis — one observation per ROI within a compartment, inference
clustered at patient level (epithelial, within-compartment SERPINB2 split):

| stratum | SERPIN-high | SERPIN-low | diff | p permutation | p signed-rank |
|---|---|---|---|---|---|
| all | 13.13 | 10.98 | +2.15 | **0.034** | **0.048** |
| naive | 11.76 | 9.56 | +2.21 | 0.074 | 0.067 |
| treated | 15.07 | 14.27 | +0.80 | 0.210 | 0.297 |

CAF is not reliable (pooled +2.52 but median per-patient difference −0.07,
signed-rank p = 0.395 — outlier-driven). Immune compartment is null
(p = 0.88). ICC of the T-cell score across segments is 0.000 (design effect
1.00), and segment-level type I error is 6.8 %, so clustering was never the main
defect — compartment pooling was.

### Myeloid signal

| analysis | SERPIN-high | SERPIN-low | p (patient-clustered) |
|---|---|---|---|
| myeloid (CD163/CD68), epithelial | 23.70 | 23.46 | 0.813 |
| myeloid, CAF | 27.52 | 29.03 | 0.535 |
| myeloid, Immune | 74.18 | 78.82 | 0.545 |
| myeloid, Immune, treated | 62.63 | 75.18 | 0.231 |

No compartment shows a significant decrease in SERPIN-high regions once patients
are the unit of inference. The borderline "immune exclusion" reading from the
pooled analysis does not survive.

### Interpretation

There is no sign of T-cell exclusion from SERPIN-high regions. In epithelium the
T-cell signal is modestly *higher* in SERPIN-high ROIs (+2.15), but this is
borderline after patient clustering and absent in treated patients. This is a
weaker statement than "T-cell rich niches", and it rests on a 2-gene myeloid and
3-gene T-cell score.

**Multiplicity.** The compartment table holds 18 contrasts (3 compartments x
2 scores x 3 strata). The epithelial p = 0.034 is nominal; it does not survive
any multiplicity correction. Treat it as exploratory and directionally
consistent, not as a significant finding.

---

## 5. GSE300595 — snMultiome, 12 patients (5 R / 7 NR), 35,115 nuclei

All samples are post-treatment. No pre-treatment biopsies exist.

### Data quality

- 35,115 nuclei pass (paper: 30,840). The gap is explained by four fragment-based
  filters that cannot be computed from a count matrix alone.
- Myeloid compartment = RNA cluster 7, 2,364 nuclei, 16 subclusters.
- ATAC: 21,769 peaks detected in ≥5 % of myeloid nuclei; 18,671 tested.

### Responders vs non-responders

| measure | responders | non-responders | p | δ |
|---|---|---|---|---|
| myeloid fraction of nuclei | 13.9 % | 3.9 % | **0.048** | +0.714 |
| LAM fraction of myeloid | 42.4 % | 31.1 % | 0.202 | +0.486 |
| LAM module score | +0.040 | −0.133 | 0.106 | +0.600 |
| SERPIN⁺ fraction of malignant | 48.5 % | 42.9 % | 1.000 | −0.029 |

The myeloid-abundance difference is the only p < 0.05 result, and it is on the
edge of the resolution of n = 5 vs 7. Note also that **SERPIN fraction does not
distinguish responders from non-responders at all** in this cohort (p = 1.0), so
the cohort does not reproduce the assumption that high SERPIN goes with poor
response.

### Chromatin accessibility

- **0 peaks** reach padj < 0.25 (and 0 at padj < 0.1) among 18,671 tested.
- Programmes were therefore built by effect-size ranking with unadjusted
  p < 0.05: 370 non-responder-enriched peaks → 366 genes (226 present on the DSP
  panel); 662 responder-enriched peaks → 620 genes (389 on panel).
- The most extreme peaks sit exactly at the permutation floor for 5 vs 7
  (p = 0.002525), i.e. they cannot become more significant regardless of effect
  size.
- 349 / 370 NR peaks lie within 100 kb of a protein-coding TSS, so gene
  assignment is broadly plausible despite the nearest-TSS approximation.

### Why the epigenomic question cannot be answered here

The intended test was: project the non-responder myeloid accessibility programme
onto spatial ROIs and ask whether it is enriched in the SERPIN-high niche. It
was run, and it fails on specificity:

| test | result |
|---|---|
| NR programme in SERPIN-high vs low (segment) | p < 0.0001, δ = +0.208 |
| responder programme, same comparison | p = 0.0046, δ = +0.137 |
| NR vs responder programme correlation | **ρ = 0.883 (segment), 0.895 (ROI)** |
| NR programme, ROI level, naive stratum | p = 0.144 |
| NR programme, ROI level, treated stratum | p = 0.441 |
| partial correlation controlling epithelial+CAF+myeloid | ρ 0.187 → 0.189 (unchanged) |
| patient NR ATAC score vs SERPIN fraction (n = 12) | ρ = −0.014, p = 0.97 |
| patient NR ATAC score vs myeloid fraction (n = 12) | **ρ = 0.811, p = 0.0014** |

Both directional programmes rise together in SERPIN-high regions, so the score
is tracking a shared axis — most plausibly myeloid/transcriptional content —
rather than a response-specific myeloid state. Consistent with this, the
patient-level accessibility score tracks myeloid abundance, not SERPIN status.

**Conclusion: the data do not provide independent spatial support for a
non-responder-specific myeloid chromatin programme in the SERPIN-high niche.**

---

## 6. What the three cohorts jointly support

| claim | status |
|---|---|
| SERPIN⁺ malignant programme persists after chemotherapy | weak support (p = 0.065, direction increase) |
| MARCO⁺ macrophage niche is preserved | **not supported** (corrected: p = 0.76, direction decrease) |
| SERPIN-high niches are T-cell rich | weak/exploratory (epithelial +2.15, patient-clustered p = 0.034–0.048; naive-driven) |
| SERPIN-high niches are immune-excluded | **not supported** (no myeloid decrease survives patient clustering) |
| Non-responder myeloid chromatin programme is enriched in SERPIN-high niche | **not supported** (programmes indistinguishable, ρ = 0.88) |
| Responders have more myeloid nuclei | exploratory support (p = 0.048, n = 5 vs 7) |
| Any of the above is FDR-significant | **no** |

## 7. What would be needed to answer the original question

1. **Paired pre/post biopsies** from the same patients — the single biggest gap.
   Without them no chemotherapy effect is identifiable, only cross-sectional
   contrasts.
2. **A platform that measures `SERPINE1`, `MARCO`, `SPP1` and `MRC1` together**,
   or at least two platforms sharing the same SERPIN gene. Currently no two
   cohorts measure the same SERPIN axis.
3. **Fragments for the ATAC data**, enabling peak calling, TSS enrichment,
   nucleosome signal and chromVAR/motif analysis — i.e. an actual mechanistic
   chromatin readout rather than nearest-gene annotation.
4. **Larger response groups** (the 5 vs 7 split cannot distinguish a real
   moderate effect from noise, and pins the best attainable p at 0.0025).
5. **Spatial data with histology labels and ≥1 myeloid gene per programme
   state**, so that a myeloid-specific (not tumour-content) signal can be
   isolated.

Until then, the defensible position is narrow: chemotherapy does not clearly
erase the SERPIN⁺ tumour state, the MARCO⁺ macrophage state is not demonstrably
preserved, and no chromatin-level explanation for niche persistence can be
substantiated with these data.
