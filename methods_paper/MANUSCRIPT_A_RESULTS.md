# Path A — Methods paper (full closed-loop draft)

**Working title**

> *Five silent failures in reanalysis of public tumour–immune data, one
> ecological fallacy, and one reconstruction gap: unit of analysis, panel
> coverage, denominator, fabricated metadata, identifier loss, and the
> undeposited upstream pipeline*

**Target venues.** *Briefings in Bioinformatics* (primary, IF ~8.5);
*Computers in Biology and Medicine* (fallback, IF ~7.0); *Journal of
Translational Medicine* (short perspective).

**Thesis (one sentence).** Every error that changed our conclusions produced a
plausible number and raised no exception, so the lesson is not "be careful" but
"assert specific invariants against the study design" — and even after all five
are fixed, a spatially defined co-occurrence hypothesis does not survive
transfer to bulk expression, which is the sixth lesson. The seventh lesson
comes from the GSE289138 reconstruction: the source repository for the
Falcomatà 2026 *Nature* paper contains figure-level notebooks only; the
upstream processing that produced its intermediate h5ad files is not
deposited, so any reproduction is a reconstruction, not an exact replay, and
the chemistry effect must be quantified alongside the genotype effect before
any directional claim is made.

---

## Abstract

Public single-cell, spatial and bulk datasets are routinely reanalysed to test
hypotheses they were not designed for, and such reanalyses rarely report how the
original analytical choices shaped the result. Over a four-cohort audit
(snRNA-seq n = 43 patients / 224,988 nuclei; GeoMx DSP n = 22 patients / 194
physical ROIs / 570 segments; snMultiome n = 12 patients / 35,115 nuclei; mouse
snRNA-seq n = 6 libraries), we found five defects that each changed a conclusion
while producing no software error. (i) **Unit of analysis**: one spatial
contrast gave p = 0.0063, p = 0.685 and p = 0.034 depending on whether segments
were pooled, averaged to ROIs, or stratified by compartment with patient-clustered
inference; the first two are wrong for different reasons. (ii) **Panel
coverage**: `SERPINE1` is absent from the 9,009-gene GeoMx panel, so a two-gene
SERPIN score silently became a one-gene `SERPINB2` score, and three platforms
appeared to agree on an axis they did not measure. (iii) **Denominator
definition**: a cell-type mask built from continuous score columns with
`.any(axis=1)` matched all 224,988 nuclei, so a per-macrophage fraction was in
fact per nucleus, reversing the direction of the apparent treatment effect.
(iv) **Fabricated metadata**: a hardcoded `treatment` column in a mouse cohort
invented an anti-PD-1 contrast that does not exist; the real design is three
genotypes × two library chemistries, and the "treatment effect" was probe-capture
chemistry. (v) **Identifier loss**: the ICB-Portal expression matrix has no
gene-identifier column, so the standard `index_col=0` read yields a 348 × 347
sample matrix with every gene silently missing. We then ran one pre-specified,
single-cohort negative control: in IMvigor210 (68 responders / 230
non-responders, baseline only), `SERPINE1`, `SPP1` and `MARCO` were all higher in
non-responders — the predicted direction — but none was significant (p = 0.20 to
0.58; BH q ≥ 0.57), and the three-gene composite was likewise null (p = 0.209).
The genes do co-vary (Spearman ρ = 0.40 to 0.53), so the co-expression premise
survives; what does not survive is the leap from spatial co-occurrence to bulk
abundance. We give one invariant per failure and a five-item checklist that would
have caught all five automatically.

---

## Key points

1. Five defects that changed conclusions produced valid numbers and no exception;
   each was caught only by an external check against the study design.
2. The most severe is not statistical but documentary: a fabricated `treatment`
   column that is perfectly collinear with sequencing chemistry.
3. The same spatial contrast has three defensible-looking answers under three
   units; only the compartment-stratified, patient-clustered unit matches both
   the biology and the sampling design.
4. A spatially defined co-occurrence axis does not collapse into a bulk
   single-gene difference: a pre-specified IMvigor210 test was directionally
   consistent and statistically null.
5. Reporting minimum attainable p, gene coverage, denominator identity, metadata
   provenance and identifier axes is cheaper than any of the errors.

---

## 1. Introduction

Reanalysis of public omics data is now a primary mode of discovery. A typical
study downloads a deposited matrix, applies a modern pipeline, and reports
differential expression, cell–cell communication, trajectory or spatial
neighbourhood enrichment. What is almost never reported is the *sensitivity of
the conclusion to analytical choices that the deposit does not constrain*: which
level of the sampling hierarchy is treated as the unit, whether a requested gene
was actually measured, what the denominator of a fraction is, where the design
columns came from, and whether the matrix's row axis is even genes.

These are not exotic. In our own audit they occurred in four consecutive
projects, in four different assay types, and in every case the erroneous pipeline
ran to completion and emitted publishable-looking output. The failure mode is
therefore not "a bug" but "a valid computation of the wrong quantity", which is
strictly harder to detect: type checks, shape checks and range checks all pass.

This paper has two parts. The first documents five such failures with the exact
numbers before and after correction, and states for each one a single assertion
that would have caught it at load time. The second reports a pre-specified
negative control that tests whether the biological axis motivating the whole
project — a SERPIN-high, MARCO-positive myeloid niche associated with immune
evasion — is detectable in bulk expression at all. It is not. We report that
null result as a result, because it delimits what the preceding single-cell and
spatial observations can legitimately be used to claim.

The contribution is not a new method. It is a compact, reproducible failure
catalogue with invariants, plus a demonstration that correcting the five
analytical defects does not rescue a claim that was never identifiable in the
data to begin with.

---

## 2. Data and scope

| cohort | assay | n | structure | role in this paper |
|---|---|---|---|---|
| GSE202051 | snRNA-seq | 43 patients (18 naive / 25 treated), 224,988 nuclei | patient-level | C3 denominator failure |
| GSE199102 | GeoMx DSP, 9,009-gene panel | 22 patients, 235 physical ROIs, 570 segments | segment ⊂ ROI ⊂ patient | C1 unit failure, C2 panel failure |
| GSE300595 | 10x Multiome (ATAC + RNA) | 12 patients (5 R / 7 NR), 35,115 nuclei | patient-level | resolution floor, cross-platform check |
| GSE289138 | mouse snRNA-seq | 6 libraries, 3 genotypes × 2 chemistries, n = 1 tumour/genotype | library = genotype × chemistry | C4 fabricated metadata |
| IMvigor210 (Mariathasan 2018) | bulk RNA-seq, urothelial, atezolizumab | 348 baseline (68 R / 230 NR / 50 NE) | patient-level | pre-specified negative control |

No patients overlap across the human cohorts. All cross-cohort comparisons are
consistency checks, not validations. GSE289138 is a mouse orthotopic model and is
used only to demonstrate the metadata failure; it supports no treatment
comparison.

The audit was not a systematic survey. The five failures are those that occurred
in one laboratory's four projects over a defined period. We make no prevalence
claim; a prevalence estimate would require a screened corpus with a pre-registered
coding scheme.

---

## 3. Failure 1 — the unit of analysis (C1)

### 3.1 What happened

A single biological question — "is T-cell signal higher in SERPIN-high spatial
regions?" — produced three defensible-looking answers from the same 570
observations:

| analysis unit | T-cell p |
|---|---|
| pool all 570 segments | **0.0063** |
| average to 194 physical ROIs | **0.685** |
| compartment-stratified, patient-clustered | **0.034** |

All three are computable in a few lines. Only the third is defensible, and the
reason is structural rather than statistical.

### 3.2 Root cause

The design is **crossed, not nested**, and this is invisible unless checked:

- 194 physical ROIs, each split into at most one segment per compartment
  (Immune / Epithelial / CAF). Verified: **0 duplicated `(ROI, compartment)`
  pairs**.
- Within a compartment, each observation is therefore a *distinct ROI*.
- Across compartments, segments of one ROI are *different tissues*, not
  replicates.

Two consequences follow:

1. **Pooling compartments** mixes tissues with different baselines (T-cell score
   Immune ≈ 25, Epithelial ≈ 12). This is a confounding error, not primarily a
   replication error.
2. **Averaging to ROIs** injects the between-compartment variance into every
   observation and destroys compartment-specific effects.

Clustering, the error we *expected*, was minor. One-way ICC of the T-cell score
across segments nested in ROIs is **0.000** (design effect 1.00); for SERPINB2 it
is 0.168 (design effect 1.33). A null simulation (2,000 permutations with SERPIN
status assigned at ROI level) gives segment-level type I error rates of:

| score | type I, segment level | type I, ROI level |
|---|---|---|
| T-cell | 0.068 | 0.054 |
| Myeloid | 0.019 | 0.044 |
| SERPINB2 | 0.095 | 0.054 |
| NR accessibility programme | 0.103 | 0.049 |
| R accessibility programme | 0.069 | 0.054 |

Segment-level testing is therefore mildly liberal (6.8–10.3 %), **not** enough to
explain a move from p = 0.0063 to p = 0.685. Pseudo-replication is over-diagnosed;
design crossing is under-diagnosed.

### 3.3 Corrected analysis

One observation per ROI within a compartment; inference clustered at patient
level (22 patients, 7–15 ROIs each). Two independent patient-aware tests:
signed-rank on per-patient paired differences, and within-patient label
permutation (5,000 draws).

Epithelial compartment, within-compartment SERPINB2 split:

| stratum | mean high | mean low | diff | p pooled | p signed-rank | p permutation |
|---|---|---|---|---|---|---|
| all | 13.13 | 10.98 | +2.15 | 0.00046 | **0.048** | **0.034** |
| naive | 11.76 | 9.56 | +2.21 | 0.00113 | 0.067 | 0.074 |
| treated | 15.07 | 14.27 | +0.80 | 0.370 | 0.297 | 0.210 |

CAF illustrates a **second** trap in the same table: pooled diff +2.52
(p = 0.00025) but median per-patient difference −0.07 and signed-rank p = 0.395 —
the pooled effect is outlier-driven. Reporting only the pooled test would have
been wrong in the opposite direction. Immune compartment: null (permutation
p = 0.88). Myeloid: no compartment shows a significant *decrease* in SERPIN-high
regions after patient clustering (Immune all: diff −4.65, signed-rank p = 0.62).

All 18 contrasts in the compartment table (3 compartments × 2 scores × 3 strata)
are uncorrected. The epithelial p = 0.034 is nominal and does not survive
multiplicity correction; it is reported as exploratory.

### 3.4 Sidebar — the resolution floor

The same project contains a smaller but cleaner arithmetic trap. With 5
responders vs 7 non-responders, the minimum attainable two-sided Mann-Whitney p
is

```
2 / C(12,5) = 2 / 792 = 0.002525
```

Several differentially accessible peaks sat at exactly this value, i.e. maximal
rank separation. Among 18,671 peaks tested, **none** reached padj < 0.25, and no
gene reached padj < 0.1 in pseudobulk RNA. Reporting a p at the floor as
"significant" would be arithmetically impossible to strengthen — the statistic
has no room left. Any test statistic at the floor should be flagged as such
rather than interpreted as strong evidence.

### 3.5 Invariant that would have caught it

> For every observation, assert the full nesting path
> (`segment → ROI → patient`) and refuse to test a contrast that crosses a level
> at which the predictor varies in the opposite direction to the outcome
> baseline. Report ICC and design effect before choosing the test. Report the
> minimum attainable p for the group sizes.

---

## 4. Failure 2 — silent gene dropout from a targeted panel (C2)

### 4.1 What happened

A scoring function selected genes by filtering a requested list against the
matrix:

```python
pick = lambda names: [g for g in names if g in matrix.columns]
SERPIN_G = pick(["SERPINE1", "SERPINB2"])
```

On the 9,009-gene GeoMx panel, `SERPINE1` is **absent** and `SERPINB2` is
present, so this silently became a **one-gene `SERPINB2` score** — while the
manuscript text described a two-gene SERPIN axis. The panel does contain
`SERPINA1`, `SERPINA3`, `SERPINA5`, `SERPINB1`, `SERPINB2`, `SERPINB5`,
`SERPINB6`, `SERPINB8`, `SERPINF1`, `SERPING1`, `SERPINH1` — but not `SERPINE1`.

### 4.2 The resulting false consistency

| platform | dataset | genes actually used |
|---|---|---|
| snRNA | GSE202051 | `SERPINE1` + `SERPINB2` (SERPINB2 ~undetectable → effectively `SERPINE1`) |
| snMultiome | GSE300595 | `SERPINE1` only (`SERPINB2` in 11 / 35,115 nuclei, 0.03 %) |
| GeoMx DSP | GSE199102 | **`SERPINB2` only** |

Three cohorts, three different molecules, described as one axis. In GSE300595 the
two genes are near mutually exclusive in detectability, so "SERPIN-high" means
opposite genes on different platforms. The apparent cross-platform agreement was
an artefact of a filtering convention.

Handling after correction: DSP results are labelled by their true axis, a
sensitivity analysis using the mean of all 11 panel SERPIN genes is reported
alongside, and cross-platform statements refer to "the SERPIN-high compartment"
without implying that the same gene was measured.

### 4.3 Invariant that would have caught it

> Never filter a requested gene list silently. Assert
> `requested ⊆ panel` at load time, print a coverage table, and fail loudly on
> any missing gene that contributes to a reported score.

A one-line `assert set(genes) <= set(panel), missing` would have prevented this.
Targeted panels are usually treated as "the transcriptome"; coverage is an
analysis parameter, not metadata.

---

## 5. Failure 3 — denominator from boolean coercion (C3)

### 5.1 What happened

A myeloid mask was built from per-cell-type *score* columns:

```python
myeloid_cols = [c for c in obs.columns
                if any(x in c.lower() for x in ["monocyte","macrophage","m0","m1","m2"])
                and "score" not in c.lower()]
myeloid_mask = obs[myeloid_cols].any(axis=1)
```

The matched columns (`Macrophage`, `Monocytes_1/2`, `M0`, `M1`, `M2`) are
**continuous floats**, not booleans. `.any(axis=1)` treats any non-zero value as
True, so the mask covered **all 224,988 nuclei**. The diagnostic signature was
present in the output and not noticed:

```
n_myeloid_total = 224,988   # == total cells
```

A "MARCO⁺ fraction of macrophages" was therefore MARCO⁺ **per nucleus**.

### 5.2 Effect: the direction reverses

Corrected denominators come from the dataset's own annotation (`Level 1 ==
"Myeloid"` → 11,938 nuclei; `Level 3 == "Macrophage"` → 10,140 nuclei), with
MARCO positivity from the raw count column.

| metric | naive (n = 18) | treated (n = 25) | Cliff's δ | p |
|---|---|---|---|---|
| **MARCO⁺ / macrophage** (correct) | **8.20 %** | **4.71 %** | +0.058 | **0.758** |
| MARCO⁺ / myeloid (correct) | 7.15 % | 3.93 % | +0.022 | 0.912 |
| MARCO⁺ per 100 nuclei (old metric) | 0.291 | 0.472 | −0.004 | 0.990 |
| macrophage fraction of nuclei | 3.54 % | 6.30 % | — | 0.143 |

Per nucleus, MARCO detection *rises* with treatment (0.29 → 0.47 per 100 nuclei).
But macrophage density also rises (3.5 % → 6.3 %), and once normalised to
macrophages the MARCO⁺ fraction *falls* (8.2 % → 4.7 %). Neither is significant.
MARCO is detected in only 764 / 224,988 nuclei (0.34 %), so power is negligible.

The original claim — "chemotherapy does not eliminate the MARCO⁺ macrophage
niche" — is **not supported**. The corrected statement is that macrophage density
increases modestly and the MARCO⁺ fraction of macrophages decreases modestly
after chemotherapy, neither significantly, with detection too sparse for a firm
conclusion in either direction. The apparent "preservation" was an artefact:
rising density divided by an all-nuclei denominator.

### 5.3 Invariant that would have caught it

> Assert that a denominator is an identity, not a coincidence:
> `mask.dtype == bool` and `mask.sum() < 0.5 * n_total`, and evaluate the mask
> against a categorical annotation rather than a score column.

A unit test such as `assert n_masked < 0.5 * n_total` fails immediately here.
Boolean coercion of `float` score matrices is a recurring API footgun in
single-cell tooling, and it fails in the direction that makes a fraction look
stable.

---

## 6. Failure 4 — fabricated treatment metadata (C4)

### 6.1 What happened

The mouse cohort was read with a hardcoded metadata table in
`serpin_repro/pdac_sc_preprocessing.py`:

```python
"GSM8784900_CF111": {"sample": "CF111", "condition": "classic", "treatment": "Vehicle"},
...
"GSM9575309_CF111_Flex": {"sample": "CF111", "condition": "Flex", "treatment": "anti-PD-1"},
```

That table was the sole source of the `treatment` column. `step5_axis.py`
compared the two values and its docstring described the output as "effect of
anti-PD-1 (Flex) vs vehicle (classic)".

### 6.2 What the deposit actually contains

The GEO series matrix header (`GSE289138_series_matrix.txt.gz`) reads:

| GSM | `!Sample_title` |
|---|---|
| GSM8784900 | Pancreatic tumor and TME cells, **control KO** |
| GSM8784901 | Pancreatic tumor and TME cells, **Serpinb2 KO** |
| GSM8784902 | Pancreatic tumor and TME cells, **Serpine1 KO** |
| GSM9575309 | Pancreatic tumor and TME cells, **control KO - Flex** |
| GSM9575310 | Pancreatic tumor and TME cells, **Serpinb2 KO - Flex** |
| GSM9575311 | Pancreatic tumor and TME cells, **Serpine1 KO - Flex** |

`!Sample_characteristics_ch1` gives only `cell line: KPC FC1245 orthotopically
injected`. **There is no treatment field anywhere in the deposit.**

The real design is **3 genotypes × 2 library chemistries (3' classic vs Flex)**,
with **one tumour per genotype**; the two chemistries are technical replicates of
the same tumour, not different arms.

### 6.3 Why it was silent

The fabricated column produced entirely plausible numbers and every test ran
without error. The "effect" it reported was in fact the chemistry difference:

| readout | classic | Flex |
|---|---|---|
| median genes per cell | 1,458 | 3,185 |
| `Gapdh` detection | 91.7 % | **0.0 %** |
| `Col1a1` detection | 5.6 % | 88.1 % |
| `Krt19` detection | 14.4 % | 82.7 % |
| features in the library | 32,285 | 19,059 (13,192 classic-only) |
| `Serpine1`⁺ epithelial nuclei | 12.8–26.2 % | 45.9–51.2 % |

The last row is exactly the kind of number that would have been reported as
"anti-PD-1 expands the SERPINE1⁺ compartment". It is a probe-capture effect: Flex
is a targeted probe chemistry and detects many genes at far higher rates than 3'
classic sequencing. `treatment` was perfectly collinear with `condition` (0 cells
off-diagonal), so no amount of modelling could separate them.

### 6.4 Consequence

No anti-PD-1 claim can be made from GSE289138. The "mouse reproduction of the
Falcomata treatment effect" does not exist. The cohort can support only a
**descriptive genotype comparison** (control KO vs Serpine1 KO vs Serpinb2 KO)
with n = 1 tumour per genotype and two technical chemistry replicates — and even
that cannot confirm knockdown, because Perturb-map tumours are mosaic: `Serpine1`
detection is in fact *highest* in the Serpine1 KO library (0.46 vs 0.21 mean
counts per epithelial cell, classic), which is what a mosaic tumour with a
minority of edited clones looks like.

The existing `GSE289138_annotated.h5ad` still carries the invalid `treatment`
column and must not be used until the pipeline is re-run by genotype.
`serpin_repro/processed/QUARANTINE.md` is the physical freeze marker.

The pipeline was re-run by genotype (`pdac_sc_preprocessing.py` with
`enforce_metadata()`, then `step3_annotate_cnv.py` and `step4_myeloid.py`).
The resulting `GSE289138_processed.h5ad` carries `genotype` and no `treatment`.
`serpin_repro/processed/QUARANTINE.md` is updated to reflect the genotype re-run.
`methods_paper/audit_mouse_confound.py` then computed per-library readouts and
descriptive genotype-vs-chemistry effect sizes; the key numbers are reported in
`CORRECTIONS.md` C6 and in the manuscript's §6.6.

### 6.6 Reproduction audit: per-library readouts and genotype-vs-chemistry effect sizes

After the genotype re-run, `audit_mouse_confound.py` computed the following
per-library readouts (descriptive only; no cell-level p-values):

| readout | classic control_KO | classic Serpinb2_KO | classic Serpine1_KO | Flex control_KO | Flex Serpinb2_KO | Flex Serpine1_KO |
|---|---|---|---|---|---|---|
| T/NK % | 20.25 | 15.00 | 16.84 | 6.22 | 11.22 | 7.29 |
| Macrophage % | 19.13 | 29.81 | 29.23 | 15.57 | 19.99 | 17.84 |
| CAF % | 5.87 | 4.34 | 5.12 | 11.26 | 17.39 | 10.96 |
| Epithelial % | 25.79 | 15.09 | 12.57 | 58.17 | 37.44 | 52.94 |
| Neutrophil % | 14.55 | 13.84 | 18.82 | 2.48 | 3.40 | 5.15 |
| T/NK exhaustion | 1.869 | 1.100 | 1.472 | 0.557 | 1.047 | 0.891 |
| T/NK effector | 1.745 | 1.550 | 1.769 | 0.459 | 1.067 | 0.849 |
| T/NK naive | 0.607 | 0.600 | 0.667 | 0.272 | 0.626 | 0.355 |
| SPP1+MARCO+ macro % | 2.30 | 0.98 | 0.44 | 2.35 | 2.65 | 2.56 |
| Cxcl9+ macro % | 39.06 | 33.11 | 46.11 | 30.31 | 40.71 | 40.00 |
| Serpine1+ mal epi % | 10.32 | 13.86 | 23.67 | 35.15 | 23.35 | 28.44 |
| Serpinb2+ mal epi % | 0.00 | 0.00 | 0.00 | 0.08 | 0.06 | 0.17 |

**Genotype vs chemistry effect sizes (descriptive).**

| readout | genotype_eff | chemistry_eff | ratio |
|---|---|---|---|
| T/NK % | 5.125 | 9.120 | 0.56 |
| Macrophage % | 7.550 | 8.257 | 0.91 |
| CAF % | 3.980 | 8.093 | 0.49 |
| T/NK exhaustion | 0.629 | 0.649 | 0.97 |
| T/NK effector | 0.413 | 0.896 | 0.46 |
| T/NK naive | 0.211 | 0.224 | 0.94 |
| SPP1+MARCO+ macro % | 1.080 | 1.280 | 0.84 |
| Cxcl9+ macro % | 11.700 | 7.487 | 1.56 |
| Serpine1+ mal epi % | 12.575 | 13.030 | 0.97 |
| Serpinb2+ mal epi % | 0.055 | 0.103 | 0.53 |

**Directional checks against the paper's claims.**

- T/NK exhaustion (classic; paper: lower in Serpinb2-KO and Serpine1-KO):
  control_KO 1.869, Serpinb2_KO 1.100, Serpine1_KO 1.472. Both KOs are
  lower than control — direction consistent with the paper (signature
  reductions of 1.70× and 1.27×; the paper reports 3.5× and 1.7× for
  Tex-term subset fractions in Fig. 2f — different metric, same sign).
  In Flex the sign reverses (0.557 vs 1.047 / 0.891): chemistry flips
  this readout, which is why the pooled object must never be tested.
- T/NK effector (classic; paper: higher in Serpinb2-KO and Serpine1-KO):
  the all-T/NK signature is flat (1.745 / 1.550 / 1.769) because NK cells
  dilute it (NK/other = 25.7/47.0/38.9% of classic T/NK). At the paper's
  own granularity, the CD8 Teff fraction with NK excluded RISES
  (33.79 → 48.86/49.38%, 1.45×/1.46×) — supported at subset level.

**Subset-level reproduction (added 2026-11-15; full tables in
`methods_paper/reproduction_report.tex`, compiled 24 pp. including the
source-data audit and the R-redrawn figures). ** T cells
(classic only, n = 3,017, fractions renormalised over CD4⁺CD8⁺ T cells as
in Fig 2f): CD8 Tex total 41.55 → 26.97/23.70% (1.54×/1.75× reductions;
paper reports 3.5×/1.7× for the terminal subset); CD8 Teff 33.79 →
48.86/49.38%. Macrophages (classic): SPP1⁺/MARCO⁺ 12.12 → 5.57% (2.18×
reduction in Serpine1-KO), Cxcl9⁺ 3.02 → 6.05/7.82%, proliferating
10.89 → 5.80/6.05% — three Fig 4g claims reproduce; Arg1⁺ moves opposite
(13.60 → 21.67/23.22%). None of the macrophage claims reproduce in Flex.
CAFs (Flex-only): the "subtype composition largely preserved" claim is NOT
reproduced (matrix CAF 62.06 → 13.33/24.78%; LRCC15⁺ 0 → 1.64/22.56%).
Updated verdict: **7 of 8 testable claims directionally consistent in the
intended chemistry; 2 fail; 1 not characterisable (SERPINB2⁺ malignant
state, 0.00–0.17% detection).**
- SPP1+MARCO+ macrophages (classic; paper: lower in Serpine1-KO):
  control_KO 2.30, Serpinb2_KO 0.98, Serpine1_KO 0.44. Both KOs lower;
  direction consistent with the paper.
- CAF fraction (Flex; paper: lower in Serpine1-KO):
  control_KO 11.26, Serpinb2_KO 17.39, Serpine1_KO 10.96. Serpine1-KO is
  *lower* than control; Serpinb2-KO is higher.
- SERPINE1+ malignant epithelial: classic 10.32–23.67 %; Flex 23.35–35.15 %.
  In classic the Serpine1-KO library is highest (23.67 %); in Flex the
  control library is highest (35.15 %), Serpine1-KO second (28.44 %).
  The correct statement is that Serpine1 detection *persists* in the
  Serpine1-KO libraries — expected for a mosaic Perturb-map tumour with a
  minority of wild-type clones, and the reason knockdown cannot be
  confirmed from this deposit.

**Interpretation.**

1. With one tumour per genotype, no cell-level test is valid. The unit of
   inference is the tumour (n = 1 per genotype), so differences are
   descriptive only.
2. For most readouts, the chemistry effect is comparable to or larger than
   the genotype effect (ratios 0.46–0.97). A genotype claim from a pooled
   analysis is not identifiable for those readouts.
3. Directional checks in the intended chemistry (classic for T cells,
   Flex for stroma): T/NK exhaustion consistent (Tex total 1.54×/1.75×
   reductions); T/NK effector consistent at CD8-subset level (Teff
   1.45×/1.46× increases; the flat all-T/NK signature was NK-diluted);
   SPP1⁺/MARCO⁺ macrophages consistent in classic (2.18× reduction);
   CAF fraction consistent (Serpine1-KO marginally lower); CAF subtype
   preservation NOT reproduced.
4. The source repository contains figure-level notebooks only; the upstream
   processing that produced its intermediate h5ad is not deposited. This is
   a reconstruction, not an exact replay.
5. **Nature source-data audit (added 2026-11-15).** The non-deposited in
   vivo / flow / IHC / survival panels were re-plotted and p-values
   recomputed from the Nature source-data spreadsheets (MOESM9--17,
   archived in `serpin_repro/data/Nature_SourceData/`). Recomputed Welch
   p-values track the claimed numbers closely; the 2l log-rank contrasts
   recompute to values *smaller* than claimed (0.000470/0.000181 vs
    0.0002/0.0005). Fig. 2n claimed set {0.0003, 0.0024, 0.0143}
    re-maps as: 0.0003 = Vehicle vs combo (0.0002903); 0.0024 = PAI-039
    vs combo (0.002402); 0.0143 = Vehicle vs PAI-039 (0.01429) -- the
    legend attributes 0.0003 to Vehicle vs PAI-039, but that contrast
    actually gives 0.0143. The 3j/3l p-value set is confirmed but
    permuted in attribution. Fig. 5e claimed 0.0578 = Serpine1-KO vs
    mix (0.05779), not a Control contrast. Fig. 1 (Perturb-map screen)
    has no source-data file -- a documented reporting gap.
    Full details: `methods_paper/reproduction_report.tex`
    Section `sec:sourcadata`; figure re-plots in
    `methods_paper/figs_reproduction/figG--figJ_sd_*.png`.

### 6.5 Invariant that would have caught it

> Never hardcode per-sample metadata. Parse it from the deposit's own sample
> table, then assert that (i) every value in each design column appears in the
> parsed source, and (ii) no design column is perfectly collinear with a
> technical column (chemistry, batch, run date).

Here, either check fails immediately: `treatment` has no counterpart in the
series matrix, and `treatment` is 100 % determined by `condition`.

---

## 7. Failure 5 — identifier loss in a public bulk matrix (C5)

### 7.1 What happened

The standard bulk-matrix convention is
`pd.read_csv(path, sep="\t", index_col=0)`, i.e. first column = gene symbol.
Applied to ICB-Portal's `Mariathasan_2018_expression.tsv`, this yields a
**348 × 347 matrix of samples**.

The TSV is **18,021 rows × 348 columns** and the first column is a **sample ID**
(`SAMf2ce197162ce`), not a gene. There is no gene-identifier column at all:

```
line 1: n_fields=348  first3=['SAMf2ce197162ce', 'SAM698d8d76b934', ...]
line 2: n_fields=348  first3=['0.740221204305847', '3.63390276783659', ...]
```

Gene names exist only in `Mariathasan_2018.Rdata`, as `TPM` `dim_0`; sample IDs
are `TPM` `dim_1`. The Rdata has 18,021 rows but **18,020 unique symbols** — one
symbol is duplicated and must be averaged before testing.

### 7.2 Why it was silent

`read_csv` succeeds, the shape is plausible, and `index_col=0` produces a
valid-looking numeric frame. A gene lookup then fails and a naive script prints
"gene not found, skipped" for every target; a script that indexes rows
positionally would analyse arbitrary genes with no error at all. The failure is
not an exception — it is a silent redefinition of the row axis.

### 7.3 Handling

Gene names were recovered from the Rdata, and the TSV was proved to be the same
matrix before its names were used:

```python
assert list(tsv.columns) == rdata_sample_ids
assert tsv.shape == rdata_matrix.shape
assert np.allclose(tsv.values, rdata_matrix)
```

Only then were the Rdata names attached, and the duplicated symbol was averaged.
The same check must be applied to every other `*_expression.tsv` in ICB-Portal
until proven otherwise.

### 7.4 Invariant that would have caught it

> Assert that the expression matrix has a gene-identifier axis. If it does not,
> parse the authoritative container for names and assert
> `np.allclose(tsv.values, authoritative_matrix)` and
> `list(tsv.columns) == authoritative_sample_ids` before any gene is touched.
> Also assert the response coding against observed values (`Pre`, not `pre`;
> `CR/PR`/`SD/PD`, not `R`/`NR`) instead of a guessed regex.

---

## 8. The ecological fallacy — a pre-specified negative control (IMvigor210)

### 8.1 Why this test, and why only one

The five failures above are analytical. They do not by themselves test whether
the motivating biology is real. The motivating claim is a *spatial
co-occurrence*: SERPIN-high tumour regions are enriched for a MARCO-positive
myeloid compartment that excludes T cells. Bulk RNA-seq cannot measure
co-occurrence; it measures a homogenate. The question this section asks is
therefore narrow and pre-specified:

> If that axis is real, is its bulk expression at least directionally associated
> with ICB non-response in an independent cohort?

One cohort, one contrast, three pre-specified genes, one pre-specified statistic.
No cohort shopping, no gene substitution, no threshold adjustment.

### 8.2 Design, fixed before the result was seen

- **Cohort**: IMvigor210 / Mariathasan 2018, urothelial carcinoma, atezolizumab.
- **Baseline only**: `Pre_On == "Pre"` — all 348 samples, asserted rather than
  assumed.
- **Response**: RECIST CR/PR (R) vs SD/PD (NR); 50 NE excluded → **68 R / 230
  NR**.
- **Primary family**: `SERPINE1`, `SPP1`, `MARCO`, two-sided Mann-Whitney U,
  Benjamini-Hochberg across the three.
- **Effect size**: Cliff's δ (+1 = higher in responders, −1 = higher in
  non-responders).
- **Secondary (exploratory)**: mean z-score of the three genes, nominal only.
- **Descriptive**: pairwise Spearman ρ among the three genes.

### 8.3 Primary results

| gene | median R | median NR | mean R | mean NR | U | p | q (BH) | Cliff's δ |
|---|---|---|---|---|---|---|---|---|
| SERPINE1 | 202.53 | 238.45 | 389.07 | 372.96 | 7268 | 0.377 | 0.566 | −0.071 |
| SPP1 | 233.97 | 269.58 | 671.64 | 807.69 | 7476 | 0.582 | 0.582 | −0.044 |
| MARCO | 1.73 | 2.87 | 7.75 | 9.82 | 7021 | 0.201 | 0.566 | −0.102 |

**Secondary (exploratory).** Mean z-score of the three genes: median R −0.319 vs
NR −0.255, U = 7035, **p = 0.209**, δ = −0.100.

**Descriptive co-variation** (labelled n = 298): SERPINE1–SPP1 ρ = +0.402;
SERPINE1–MARCO ρ = +0.434; SPP1–MARCO ρ = +0.532.

### 8.4 Verdict: directionally consistent, statistically null

All three genes and the composite are higher in non-responders — the predicted
direction — but none reaches nominal significance, and BH q is ≥ 0.57. The effect
sizes are small (|δ| = 0.04–0.10). With 68 vs 230 the cohort can resolve moderate
effects (|δ| ≈ 0.2) but not small ones; the observed effects are therefore
compatible with a weak or absent bulk association, and this dataset cannot
distinguish those two.

Two details that must not be smoothed over:

- **SERPINE1's mean is higher in responders** (389.07 vs 372.96) while its median
  and Cliff's δ are higher in non-responders. The mean is tail-driven; the rank
  statistic was the pre-specified one. The mean must not be reported as a
  direction.
- The three genes **do co-vary moderately** (ρ = 0.40–0.53), so the
  "co-expression module" premise is not contradicted — but co-variation is not
  association with response.

### 8.5 What this does and does not establish

- It does **not** support "SERPIN⁺/MARCO⁺ bulk expression predicts ICB
  resistance" as a stand-alone claim.
- It does **not** refute the spatial hypothesis. Bulk TPM is the wrong instrument
  for a co-occurrence claim, and this is one cancer type with one drug.
- It does establish a boundary: a pattern that is visually obvious in
  single-cell and spatial data does not automatically appear as a bulk
  single-gene difference. Any future test needs a pre-specified composite in an
  independent cohort with the same endpoint, reported whatever the result.

This is the sixth lesson of the paper, and the only one that is biological rather
than computational: **the ecological fallacy is not avoided by using better
statistics on the wrong measurement scale.**

---

## 9. Synthesis — why all six were silent, and what to assert

| failure | silent because | cheap invariant |
|---|---|---|
| 1. unit of analysis | every unit yields a valid number | print ICC / design effect; assert the nesting path; report minimum attainable p |
| 2. panel dropout | list comprehension filters without error | `assert requested ⊆ panel`; coverage table |
| 3. denominator | `.any()` on floats returns a mask | `assert mask.dtype == bool`; `mask.sum() << n`; use an annotation |
| 4. fabricated metadata | hardcoded dict looks like data | parse the deposit's sample table; assert no design × technical collinearity |
| 5. identifier loss | `read_csv` succeeds with a plausible shape | assert a gene-identifier axis; prove the matrix against its authoritative container |
| 6. reconstruction gap | source repo deposits figure-level notebooks only; upstream h5ad not shared | re-run the pipeline by genotype; report per-library readouts and genotype-vs-chemistry effect sizes descriptively; never claim a treatment effect |

The common structure: each failure preserves *types* and *ranges*, so nothing
raises, and each is detectable only by checking the analysis against the **study
design** rather than against the data.

### Checklist for reanalyses (proposed)

1. Draw the sampling hierarchy explicitly; test at the level the hypothesis
   refers to; report ICC and design effect before choosing the test.
2. Print a gene-coverage table for every score and fail on silent dropout.
3. State every denominator as a cell count and assert it against a categorical
   annotation, not a score column.
4. Compute the minimum attainable p for the group sizes and flag statistics at
   the floor.
5. Parse all design metadata from the deposit itself; assert that no design
   column is collinear with a technical column, and that every expression matrix
   has a verified identifier axis.
6. Report the same contrast under ≥ 2 defensible units; if conclusions differ,
   the unit is doing the work.
7. For any claim that crosses measurement scales (spatial → bulk), state the
   scale at which it is falsifiable and pre-specify one test.

Items 1–6 are cheap: each is a few lines at load time. Item 7 is the expensive
one, and it is the one that decides whether a paper has a result or only a
pipeline.

---

## 10. Discussion

### 10.1 Relation to existing methodological literature

The statistical content of Failures 1 and 3 is not new. Cluster-robust and
multilevel inference, and the distinction between crossed and nested designs, are
standard. What is new is the observation that in spatial-omics practice the
design is *neither* cleanly nested nor cleanly crossed, and that the field's
default unit (the segment, or the spot) is chosen by the assay's geometry rather
than by the hypothesis. Likewise, denominator errors from boolean coercion are
ordinary programming mistakes; the point is that single-cell APIs make them easy
and silent because score matrices and annotation matrices share a shape and a
dtype.

Failure 4 belongs to a different category: data provenance. It is not a
statistical error and cannot be fixed by a better model, because the contrast it
encodes does not exist. The generalisable warning is that hardcoded per-sample
dictionaries are indistinguishable from parsed metadata once they are in a
dataframe, and that collinearity with a technical variable is the fingerprint of
a fabricated design column.

Failure 5 is an interoperability defect specific to public portals that publish
"the same" matrix in two containers (TSV and Rdata) with different amounts of
metadata. The TSV is lossy; the Rdata is authoritative. Assuming the TSV is
self-describing is the trap.

Failure 6 (the reconstruction audit of GSE289138) adds a fourth category:
**code availability as a design constraint.** The source repository for the
Falcomatà 2026 *Nature* paper contains figure-level notebooks only. The upstream
processing that produced its intermediate h5ad files is not deposited. The
reproduction is therefore a reconstruction, not an exact replay. This is not a
statistical error; it is a reporting gap that affects what claims the deposit
can support. The practical consequence is that any genotype comparison from
GSE289138 is descriptive only (n = 1 tumour per genotype), and the chemistry
effect must be quantified alongside the genotype effect before any directional
claim is made.

### 10.2 Why this matters for the reanalysis literature

A large fraction of published reanalyses report a p-value without reporting the
unit, the panel coverage, the denominator, the provenance of design columns, or
the identifier axis. Each omission is individually defensible and collectively
fatal: the same dataset can yield opposite conclusions under choices that the
paper does not disclose. Our three answers to one spatial question (p = 0.0063,
0.685, 0.034) are a compact demonstration that "reproducible" is not the same as
"identified".

### 10.3 The negative control as a discipline

The IMvigor210 test was designed to be able to fail, and it did. Reporting it
serves two purposes. First, it prevents the five corrected analyses from being
retrospectively assembled into a claim they do not support. Second, it makes the
scale argument concrete: even with 298 labelled patients, the bulk association is
|δ| ≤ 0.10. A spatial co-occurrence hypothesis is not a bulk expression
hypothesis, and no amount of cohort-swapping changes that.

### 10.4 What would change the conclusion

The spatial claim would be strengthened by direct, patient-level measurement of
co-occurrence in tissue with both tumour and myeloid compartments and a response
annotation. No public single-cell cohort currently provides all three: every
response-labelled single-cell ICB cohort we audited is CD45⁺ immune-sorted
(GSE120575: 16,291 CD45⁺ cells, no tumour cells; GSE123813: sorted, no response
column; GSE115978: no response annotation). The PDAC-specific gap is the same:
the only cohort with RECIST (GSE248014) is tumour-cell-only (5 R / 40 NR), and
the only cohort with myeloid cells (GSE197613) is sorted and has no RECIST. This
is a data-availability fact, not a statistical one, and it is why the project
stops here.

---

## 11. Limitations

1. **Case series, not a survey.** Five failures from one laboratory's projects;
   no prevalence claim is possible without a screened corpus.
2. **One disease area.** PDAC and urothelial carcinoma. The failure mechanisms
   are assay-generic, but the specific numbers are not.
3. **The negative control is one cohort, one drug, one endpoint.** It cannot
   distinguish "no bulk association" from "small bulk association below
   |δ| ≈ 0.2".
4. **SERPINE1 mean and median disagree in direction.** The rank statistic was
   pre-specified; the mean is reported only to document the tail.
5. **No FDR-significant result survives anywhere in the project.** GSE202051
   SERPIN fraction p = 0.065; corrected MARCO fraction p = 0.758; GSE300595
   pseudobulk DE 0 genes at padj < 0.1; DAR 0 peaks at padj < 0.25; LAM fraction
   p = 0.202; IMvigor210 q ≥ 0.57. All conclusions are directional and
   hypothesis-generating.
6. **GSE289138 supports no treatment comparison at all** and only a descriptive,
   n = 1 per genotype comparison; knockdown cannot be confirmed in mosaic
   tumours. The source repository contains figure-level notebooks only; the
   upstream processing is not deposited, so any reproduction is a
   reconstruction, not an exact replay.
7. **The DSP bridge is non-specific**: NR and R accessibility programmes are
   indistinguishable in DSP (ρ = 0.883 segment-level, 0.895 ROI-level), so DSP
   cannot test the direction of the myeloid programme.
8. **The patient-level ATAC score is partly a cellularity proxy** (ρ = 0.811 with
   myeloid fraction), and the R-vs-NR peak set was selected on that contrast.
9. **Nearest-TSS gene assignment is approximate**; a ≤ 100 kb version is reported
   alongside and gives the same conclusions.
10. **Cross-platform statements concern compartments, not genes**, after C2.

Full itemised constraints are in `LIMITATIONS.md`.

---

## 12. Target venues and framing

| venue | fit | framing |
|---|---|---|
| *Briefings in Bioinformatics* (IF ~8.5) | primary | review/methodological audit; five failures + invariants + negative control |
| *Computers in Biology and Medicine* (IF ~7.0) | fallback | computational reproducibility and data-integrity angle |
| *Journal of Translational Medicine* | short perspective | "what public ICB reanalyses can and cannot show" |

The manuscript is deliberately not framed as a biological discovery. Its claim is
methodological and negative, and the negative control is part of the evidence,
not an appendix.

---

## 13. Data and code availability

- **Public data**: GSE202051, GSE199102, GSE300595, GSE289138; IMvigor210 via
  ICB-Portal (`ngdc.cncb.ac.cn/icb/file/datasets/`), files
  `Mariathasan_2018_{clinical,expression}.tsv` and `Mariathasan_2018.Rdata`.
- **Code**: `methods_paper/audit_icb_cohorts.py` (`--inspect`, `--suppl`,
  `--rds`, `--rds-strings`, `--audit-imvigor`);
  `methods_paper/check_design.py`; `methods_paper/design_effect.py`;
  `methods_paper/compartment_patient_clustered.py`;
  `gse300595/step6_bridge.py`; `gse300595/step6c_dsp_recheck.py`;
  `gse202051_hum/fix_marco_denominator.py`;
  `methods_paper/audit_mouse_confound.py` (incl. subset-level functions and
  `--source-data` flag for MOESM9--17 re-audit);
  `serpin_repro/data/Nature_SourceData/` (13 MOESM xlsx files);
  `methods_paper/figs_reproduction/figG--figJ_sd_*.png` (matplotlib SD re-plots);
  `methods_paper/figs_nature.R` (R 4.6.1 + ggplot2 Nature-style re-draw of
  Fig 2d–n, 3e–l, 4a/d/e/f/j, 5e–l, ED 5c/d/j/l, ED 7b/e, ED 8a–r, and the
  six scRNA-seq figures figA–F from `plot_tables/*.tsv`);
  `methods_paper/figs_nature/*.pdf` (thirteen vector figures; all figures in
  the report are R-generated);
  `methods_paper/figs_reproduction/source_data_audit.tsv` (recomputed p table);
  `methods_paper/reproduction_report.tex` (paper-style reproduction-audit
  report of the GSE289138 scRNA-seq component; compiled PDF 24 pp.);
  `methods_paper/verify_mouse_genotype.py`.
- **Correction log**: `CORRECTIONS.md` (C1–C6, plus superseded-file table).
- **Constraints**: `LIMITATIONS.md`; cohort audit: `methods_paper/ICB_COHORT_AUDIT.md`.
- **Frozen artefact**: `serpin_repro/processed/QUARANTINE.md` marks the mouse
  cohort; `GSE289138_annotated.h5ad` carries the old `treatment` column (from
  cache) but its `anno`/`malignant`/`leiden_epi` columns are consistent with
  the genotype re-run. The clean object is `GSE289138_processed.h5ad`.

---

## 14. Figures

1. **Design schematic** (segment ⊂ ROI ⊂ patient) annotated with the three
   p-values for the same contrast.
2. **ICC / design effect per score**, showing T-cell ICC = 0.000 and SERPINB2
   ICC = 0.168.
3. **Null-simulation type I error**, segment vs ROI, per score.
4. **Denominator reversal for MARCO**: per-nucleus vs per-macrophage, naive vs
   treated.
5. **Panel coverage across the three platforms**, highlighting `SERPINE1` absent
   from GeoMx.
6. **IMvigor210 negative control**: per-gene distributions (R vs NR), BH q, and
   the composite z-score.

---

## Superseded content

This file replaces the earlier three-failure skeleton. The earlier version's
"Failure 4 (secondary) — rank tests at the resolution floor" is retained here as
§3.4 (sidebar) rather than as a standalone failure, because it is an arithmetic
constraint on interpretation, not an independent defect. C4 (fabricated
metadata), C5 (identifier loss) and the IMvigor210 negative control are new in
this version.
