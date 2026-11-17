# Corrected analyses — four defects found and fixed

Date: 2026-11-14
Scope: GSE202051 (snRNA), GSE199102 (GeoMx DSP), GSE300595 (snMultiome),
GSE289138 (mouse snRNA)

Four errors were found: three while building the GSE300595 → DSP bridge, and a
fourth (C4) while auditing what the mouse cohort could support. All four change
previously reported or previously planned conclusions. Two of them affected
headline claims. This file supersedes any earlier summary that conflicts with it.

---

## C1. DSP: the unit of analysis was wrong — and so was the first correction

This entry was revised twice. The final version supersedes both earlier readings.

### What was done before

GSE199102 was analysed as 570 independent ROIs, pooling all compartments and
comparing SERPIN-high vs SERPIN-low segments. Result:
T-cell/NK 15.05 vs 12.78, delta +2.28, **p = 0.0063** (naive).

### What is structurally true

The 570 GeoMx observations are **194 physical ROIs** split into geometric
compartments:

| segments per physical ROI | ROIs |
|---|---|
| 3 | 193 |
| 2 | 37 |
| 1 | 5 |

Critically, verified directly: **each physical ROI contributes at most ONE
segment per compartment** (0 duplicated `(ROI, compartment)` pairs). So

- *within* a compartment, every observation is a distinct physical ROI;
- *across* compartments, segments from the same ROI are different tissues
  (Immune / Epithelial / CAF), not replicates.

### My first correction was also wrong

I first aggregated every ROI to a single mean and re-tested. That gave
p = 0.685 and I concluded the T-cell claim was an artifact. **That conclusion
was itself an overcorrection.** Two diagnostics show why:

1. **Clustering was never the main problem.** ICC of the T-cell score across
   segments nested in ROIs is **0.000** (design effect 1.00). A null simulation
   (2,000 label permutations at ROI level) gives a segment-level type I error of
   **6.8 %** at alpha = 0.05 — mildly liberal, not enough to turn p = 0.0063 into
   noise.
2. **Averaging compartments is a category error.** T-cell baseline differs by
   compartment (Immune ~25, Epithelial ~12). Averaging three compartments per ROI
   injects large between-compartment variance into every ROI value and swamps a
   modest within-compartment effect.

### Correct analysis: compartment-stratified, patient-clustered

One observation per ROI inside a compartment; inference clustered at patient
level (22 patients, 7–15 ROIs each). Two independent patient-aware tests:
per-patient signed-rank on paired differences, and within-patient label
permutation (5,000 draws).

Epithelial compartment, within-compartment SERPINB2 split:

| stratum | mean high | mean low | diff | p pooled | p signed-rank | p permutation |
|---|---|---|---|---|---|---|
| all | 13.13 | 10.98 | +2.15 | 0.00046 | **0.048** | **0.034** |
| naive | 11.76 | 9.56 | +2.21 | 0.00113 | 0.067 | 0.074 |
| treated | 15.07 | 14.27 | +0.80 | 0.370 | 0.297 | 0.210 |

CAF compartment is **not** reliable: pooled diff +2.52 (p = 0.00025) but the
median per-patient difference is −0.07 and signed-rank p = 0.395. The pooled
effect is driven by a few outlier patients. Immune compartment: null
(permutation p = 0.88).

Myeloid: no compartment shows a significant *decrease* in SERPIN-high regions
after patient clustering (Immune all: diff −4.65, signed-rank p = 0.62).

### Final conclusion for C1

- The original **p = 0.0063 is not interpretable**: it pools three compartments
  with different baselines.
- The original **direction survives, weakened**: SERPIN-high *epithelial* regions
  have modestly higher T-cell/NK signal (**+2.15**, patient-clustered
  p = 0.034–0.048), carried mainly by untreated patients (naive-only p = 0.074;
  treated null).
- There is **no evidence of T-cell exclusion**, but the effect is borderline,
  compartment-specific and exploratory — not the strong "T-cell rich" claim
  originally made.
- Myeloid signal is not significantly lower anywhere after clustering; the
  borderline "immune exclusion" reading does not hold up.
- **Multiplicity:** the compartment table contains 18 contrasts
  (3 compartments x 2 scores x 3 strata). The epithelial p = 0.034 is **nominal**
  and does not survive any multiplicity correction. It is reported as an
  exploratory, directionally consistent signal — not a finding.

### The methods lesson

One dataset, three defensible units, three answers:

| unit of analysis | T-cell p |
|---|---|
| pool 570 segments | 0.0063 |
| average to 194 ROIs | 0.685 |
| compartment-stratified, patient-clustered | 0.034 (permutation) |

None of the naive choices is right. The unit must match the biological question
(compartment-specific) and the sampling design (ROI nested in patient). This is
the centrepiece of `methods_paper/`.

**Files.** `methods_paper/check_design.py`, `methods_paper/design_effect.py`,
`methods_paper/compartment_patient_clustered.py`,
`methods_paper/compartment_patient_clustered.tsv`,
`methods_paper/type_I_simulation.tsv`, `gse300595/step6c_dsp_recheck.py`.

---

## C2. The DSP "SERPIN score" was SERPINB2 only — SERPINE1 is not on the panel

**What was done before.** The DSP script scored
`SERPIN_G = pick(["SERPINE1","SERPINB2"])`. Because `pick()` silently drops
absent genes, this resolved to a **one-gene score**.

**What is actually true.** `SERPINE1` is **absent** from the GeoMx WTA panel
(9,009 genes). Verified directly. The panel contains `SERPINA1`, `SERPINA3`,
`SERPINA5`, `SERPINB1`, `SERPINB2`, `SERPINB5`, `SERPINB6`, `SERPINB8`,
`SERPINF1`, `SERPING1`, `SERPINH1` — but no `SERPINE1`.

**Consequence: the three platforms do not measure one SERPIN axis.**

| platform | dataset | SERPIN genes actually used |
|---|---|---|
| snRNA | GSE202051 | `SERPINE1` + `SERPINB2` (SERPINB2 ~undetectable → effectively `SERPINE1`) |
| snMultiome | GSE300595 | `SERPINE1` only (`SERPINB2` in 11 / 35,115 nuclei) |
| GeoMx DSP | GSE199102 | **`SERPINB2` only** |

**Handling.** DSP results are now labelled by their true axis. A sensitivity
analysis using the mean of all 11 panel SERPIN genes is reported alongside.
Cross-platform statements refer to "the SERPIN-high compartment", and never
imply that the same gene was measured.

**Files.** `gse300595/diag_dsp_panel.py`, `gse300595/step6_bridge.py`.

---

## C3. GSE202051 "MARCO+ macrophage fraction" used the wrong denominator

**What was done before.** `MARCO_frac` was intended as MARCO⁺ / macrophages and
was used to conclude the MARCO⁺ macrophage niche is preserved.

**What was actually computed.** The myeloid mask was built as

```python
myeloid_cols = [c for c in obs.columns
                if any(x in c.lower() for x in ["monocyte","macrophage","m0","m1","m2"])
                and "score" not in c.lower()]
myeloid_mask = obs[myeloid_cols].any(axis=1)
```

The matched columns (`Macrophage`, `Monocytes_1/2`, `M0`, `M1`, `M2`) are
**continuous float scores**, not booleans. `.any(axis=1)` treats any non-zero
float as True, so the mask covered **all 224,988 nuclei**
(`n_myeloid_total = 224988 = total cells`). `MARCO_frac` was therefore
MARCO⁺ **per nucleus**, not per macrophage.

**Corrected version.** Denominator taken from the dataset's own annotation
(`Level 1 == "Myeloid"` → 11,938 nuclei; `Level 3 == "Macrophage"` → 10,140
nuclei), with MARCO positivity from the raw count column.

| metric | naive (n=18) | treated (n=25) | Cliff's δ | p |
|---|---|---|---|---|
| **MARCO⁺ / macrophage** (correct) | **8.20 %** | **4.71 %** | +0.058 | **0.758** |
| MARCO⁺ / myeloid (correct) | 7.15 % | 3.93 % | +0.022 | 0.912 |
| MARCO⁺ per 100 nuclei (old metric) | 0.291 | 0.472 | −0.004 | 0.990 |
| macrophage fraction of nuclei | 3.54 % | 6.30 % | — | 0.143 |

**The direction reverses.** Per nucleus, MARCO detection *rises* with treatment
(0.29 → 0.47 per 100 nuclei). But macrophage density also rises (3.5 % → 6.3 %),
and once normalised to macrophages the MARCO⁺ fraction *falls* (8.2 % → 4.7 %).
Neither is significant. MARCO is detected in only 764 / 224,988 nuclei (0.34 %),
so power is negligible.

**Conclusion.** "Chemotherapy does not eliminate the MARCO⁺ macrophage niche" is
**not supported** by this data. The corrected statement is: *macrophage density
increases modestly and the MARCO⁺ fraction of macrophages decreases modestly
after chemotherapy; neither change is statistically significant, and MARCO
detection is too sparse to support a firm conclusion in either direction.* The
apparent "preservation" was an artifact (rising density divided by an
all-nuclei denominator).

**Files.** `gse202051_hum/fix_marco_denominator.py`,
`gse202051_hum/processed/patient_marco_corrected.tsv`,
`gse202051_hum/processed/marco_corrected_test.tsv`.

---

## C4. GSE289138: a fabricated `treatment` column, and a treatment effect that does not exist

Found 2026-11-14, after C1–C3. This is the most severe of the four because it
invents a contrast rather than mis-analysing a real one.

**What was done before.** The mouse cohort was read with a hardcoded metadata
table in `serpin_repro/pdac_sc_preprocessing.py`:

```python
"GSM8784900_CF111": {"sample": "CF111", "condition": "classic", "treatment": "Vehicle"},
...
"GSM9575309_CF111_Flex": {"sample": "CF111", "condition": "Flex", "treatment": "anti-PD-1"},
```

That table was the sole source of the `treatment` column. `step5_axis.py`
compared the two values and its docstring described the output as
"effect of anti-PD-1 (Flex) vs vehicle (classic)".

**What the deposit actually contains.** The GEO series matrix header
(`GSE289138_series_matrix.txt.gz`) reads:

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

So the real design is **3 genotypes x 2 library chemistries (3' classic vs
Flex)**, with **one tumor per genotype**; the two chemistries are technical
replicates of the same tumor, not different arms.

**Why it was silent.** The fabricated column produced entirely plausible
numbers and every test ran without error. The "effect" it reported was in fact
the chemistry difference:

| readout | classic | Flex |
|---|---|---|
| median genes per cell | 1,458 | 3,185 |
| `Gapdh` detection | 91.7 % | **0.0 %** |
| `Col1a1` detection | 5.6 % | 88.1 % |
| `Krt19` detection | 14.4 % | 82.7 % |
| features in the library | 32,285 | 19,059 (13,192 classic-only) |
| Serpine1⁺ epithelial nuclei | 12.8–26.2 % | 45.9–51.2 % |

The last row is exactly the kind of number that would have been reported as
"anti-PD-1 expands the SERPINE1⁺ compartment". It is a probe-capture effect:
Flex is a targeted probe chemistry and detects many genes at far higher rates
than 3' classic sequencing. `treatment` was perfectly collinear with
`condition` (0 cells off-diagonal), so no amount of modelling could separate
them.

**Consequence.** No anti-PD-1 claim can be made from GSE289138. The "mouse
reproduction of the Falcomata treatment effect" does not exist. The cohort can
support only a **descriptive genotype comparison** (control KO vs Serpine1 KO
vs Serpinb2 KO) with n = 1 tumor per genotype and two technical chemistry
replicates — and even that cannot confirm knockdown, because Perturb-map
tumors are mosaic: `Serpine1` detection is in fact *highest* in the Serpine1
KO library (0.46 vs 0.21 mean counts per epithelial cell, classic), which is
what a mosaic tumor with a minority of edited clones looks like.

**Invariant that would have caught it.**

> Never hardcode per-sample metadata. Parse it from the deposit's own sample
> table, then assert that (i) every value in each design column appears in the
> parsed source, and (ii) no design column is perfectly collinear with a
> technical column (chemistry, batch, run date).

Here, either check fails immediately: `treatment` has no counterpart in the
series matrix, and `treatment` is 100 % determined by `condition`.

**Files.** `serpin_repro/pdac_sc_preprocessing.py` (fixed: `treatment` removed,
`genotype` added, with an assert that a tumor maps to one genotype),
`serpin_repro/processed/sample_metadata_CORRECTED.tsv`,
`methods_paper/audit_mouse_confound.py`, `methods_paper/verify_mouse_genotype.py`.
The existing `GSE289138_annotated.h5ad` still carries the invalid `treatment`
column; it must not be used until the pipeline is re-run.

---

## C5. ICB-Portal expression TSV has no gene-identifier column

Found 2026-11-14 while setting up the IMvigor210 pilot. Lower severity than
C1–C4 (it breaks loudly if you check), but it is the same class of defect: a
file that loads without error and is not what its name implies.

**What was done before.** The standard bulk-matrix convention is
`pd.read_csv(path, sep="\t", index_col=0)`, i.e. first column = gene symbol.
Applied to `Mariathasan_2018_expression.tsv` this yields a 348 x 347 matrix.

**What the file actually contains.** The TSV is **18,021 rows x 348 columns**
and the first column is a **sample ID** (`SAMf2ce197162ce`), not a gene. There
is no gene-identifier column at all:

```
line 1: n_fields=348  first3=['SAMf2ce197162ce', 'SAM698d8d76b934', ...]
line 2: n_fields=348  first3=['0.740221204305847', '3.63390276783659', ...]
```

Gene names exist only in `Mariathasan_2018.Rdata`, as `TPM` `dim_0`; sample IDs
are `TPM` `dim_1`. The Rdata also shows 18,021 rows but **18,020 unique
symbols** — one symbol is duplicated and must be averaged before testing.

**Why it was silent.** `read_csv` succeeds, the shape is plausible, and
`index_col=0` produces a valid-looking numeric frame. A gene lookup then fails
and a naive script prints "gene not found, skipped" for every target; a script
that indexes rows positionally would analyse arbitrary genes with no error.

**Consequence.** The expression TSV alone cannot support any gene-level
analysis. Row order is only meaningful when it is proved equal to the Rdata
matrix. The same applies to every other `*_expression.tsv` in ICB-Portal until
checked.

**Invariant that would have caught it.**

> Assert that the expression matrix has a gene-identifier axis. If it does not,
> parse the Rdata for names and assert
> `np.allclose(tsv.values, rdata_matrix)` and
> `list(tsv.columns) == rdata_sample_ids` before any gene is touched. Also
> assert the response coding against observed values (`Pre`, not `pre`;
> `CR/PR`/`SD/PD`, not `R`/`NR`) instead of a guessed regex.

**Files.** `methods_paper/audit_icb_cohorts.py` (`audit_imvigor210`),
`methods_paper/ICB_COHORT_AUDIT.md`.

---

## C6. Reproduction audit: GSE289138 genotype re-run

Found 2026-11-14, after C4, by explicitly re-running the pipeline on the
genotype-corrected design.

**What was done.** `pdac_sc_preprocessing.py` was re-run with
`enforce_metadata()` and the `treatment` column removed. The output
`GSE289138_processed.h5ad` now carries `genotype` and no `treatment`.
`step3_annotate_cnv.py` and `step4_myeloid.py` were re-run on the clean
object. `audit_mouse_confound.py` then computed per-library readouts and
genotype-vs-chemistry effect sizes.

**Per-library readouts (from audit_mouse_confound.py).**

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
  control_KO 1.869, Serpinb2_KO 1.100, Serpine1_KO 1.472.
  Both KOs are lower than control; direction consistent with the paper
  (signature reductions of 1.70x and 1.27x; the paper reports 3.5x and
  1.7x for Tex-term subset fractions in Fig. 2f -- different metric, same
  sign). In Flex the sign reverses (control 0.557, Serpinb2 1.047,
  Serpine1 0.891): chemistry flips this readout, which is why the pooled
  object must never be tested.
- T/NK effector (classic; paper: higher in Serpinb2-KO and Serpine1-KO):
  signature over ALL classic T/NK cells is flat (control 1.745,
  Serpinb2 1.550, Serpine1 1.769) -- but that signature is diluted by NK
  cells (NK/other = 25.7/47.0/38.9% of classic T/NK). At the paper's own
  granularity, the CD8 Teff fraction with NK excluded RISES
  (33.79 -> 48.86/49.38%, i.e. 1.45x/1.46x increases). Supported at
  subset level; the flat signature is a denominator artefact, not a
  refutation.

**Subset-level reproduction (per-cell classification / subclustering; added
2026-11-15).** T cells restricted to classic (the Fig. 2f platform;
n = 3,017), fractions renormalised over CD4+CD8 T cells as in Fig. 2f:

| subset | control | Serpinb2-KO | Serpine1-KO |
|---|---|---|---|
| CD8 Teff | 33.79 | 48.86 | 49.38 |
| CD8 Tex total (prog+int+term) | 41.55 | 26.97 | 23.70 |
| CD4 Treg | 5.02 | 1.60 | 4.48 |

Tex total falls 1.54x / 1.75x in the two KOs (paper: 3.5x / 1.7x for the
terminal subset); CD8 Teff rises 1.45x / 1.46x (paper: concomitant
increase). The prog/int/term split is rank-based and unstable at n = 1
(all "term" cells landed in the control library), so only the combined Tex
fraction is interpreted.

Macrophages (12,390-cell myeloid subset, leiden res 2, module-labelled):
classic reproduces three Fig. 4g claims -- Spp1/Marco+ 12.12 -> 5.57%
(2.18x reduction in Serpine1-KO), Cxcl9+ 3.02 -> 6.05/7.82% (increase),
proliferating 10.89 -> 5.80/6.05% (reduction); Arg1+ moves opposite
(13.60 -> 21.67/23.22%). Flex reproduces none of these (Spp1/Marco+
28.78 -> 20.46/27.92%). The Fig. 4g platform is not stated in the source
Methods.

CAFs (Flex-only subclustering, n = 4,570): CAF fraction reproduces
marginally (control 11.26 vs Serpine1-KO 10.96; Serb2 17.39, compatible
with the paper's Serpinb2 null), but the "CAF subtype composition largely
preserved" claim is NOT reproduced: matrix CAF 62.06 -> 13.33/24.78%,
LRCC15+ MyCAF 0 -> 1.64/22.56%. n = 1 per genotype prevents separating
genotype from tumour identity.

**Updated directional verdict: 7 of 8 testable claims are directionally
consistent in the intended chemistry; 2 fail (CAF subtype preservation,
Arg1+ mirror); 1 is not characterisable (SERPINB2+ malignant state,
0.00-0.17% detection).**
- SPP1+MARCO+ macrophages (classic; paper: lower in Serpine1-KO):
  control_KO 2.30, Serpinb2_KO 0.98, Serpine1_KO 0.44.
  Both KOs lower; direction consistent with the paper.
- CAF fraction (Flex; paper: lower in Serpine1-KO):
  control_KO 11.26, Serpinb2_KO 17.39, Serpine1_KO 10.96.
  Serpine1-KO is *lower* than control; Serpinb2-KO is higher.
- SERPINE1+ malignant epithelial: classic 10.32–23.67 %; Flex 23.35–35.15 %.
  In classic the Serpine1-KO library is highest (23.67 % vs 10.32/13.86);
  in Flex the control library is highest (35.15 %), Serpine1-KO second
  (28.44 %). The correct statement is that Serpine1 detection *persists*
  in the Serpine1-KO libraries — expected for a mosaic Perturb-map tumour
  with a minority of wild-type clones, and the reason knockdown cannot be
  confirmed from this deposit.

**Interpretation.**

1. With one tumour per genotype, no cell-level test is valid. The unit of
   inference is the tumour (n = 1 per genotype), so differences are
   descriptive only.
2. For most readouts, the chemistry effect is comparable to or larger than
   the genotype effect (ratios 0.46–0.97). For Cxcl9+ macrophages the
   chemistry effect is smaller (ratio 1.56), but that readout is not a
   core paper claim.
3. Directional checks, run in the chemistry the paper used for each
   analysis (classic for T cells, Flex for stroma):
   - T/NK exhaustion: both KOs lower in classic — consistent; sign
     reverses in Flex, so the pooled object must never be tested.
   - T/NK effector: flat as an all-T/NK signature (NK-diluted) but rising
     at the CD8-subset level with NK excluded (33.79 -> 48.86/49.38%) —
     consistent.
   - SPP1+MARCO+ macrophages: both KOs lower in classic — consistent
     (2.18x reduction in Serpine1-KO); not reproduced in Flex.
   - CAF fraction: Serpine1-KO marginally lower (consistent); Serpinb2-KO
     higher, matching the paper's own null for Serpinb2 CAFs.
   - CAF subtype preservation: NOT reproduced (matrix CAF 62.06 ->
     13.33/24.78%).
4. The source repository contains figure-level notebooks only; the upstream
   processing that produced its intermediate h5ad is not deposited. This is
   a reconstruction, not an exact replay.

  5. **Nature source-data audit (added 2026-11-15).** MOESM9--17 downloaded
     and archived to `serpin_repro/data/Nature_SourceData/` (13 xlsx). The
     non-deposited in vivo / flow / IHC / survival panels were re-plotted and
     p-values recomputed (Welch $t$; no-censoring log-rank where applicable).
     Findings: recomputed values track the claimed numbers closely
     (Fig.~2d 0.001023 vs 0.001; NSG mass 0.004317/0.02754 matching
     0.0043/0.0275 to four digits; Fig.~2l log-rank 0.000470/0.000181 vs
      0.0002/0.0005 -- recomputed \emph{smaller}). Fig.~2n claimed set
      {0.0003, 0.0024, 0.0143} re-maps as: 0.0003 = Vehicle vs combo
      (0.0002903); 0.0024 = PAI-039 vs combo (0.002402); 0.0143 = Vehicle
      vs PAI-039 (0.01429) -- the legend attributes 0.0003 to Vehicle vs
      PAI-039, but that contrast actually gives 0.0143. The 3j/3l p-value
      set {0.0145, 0.0282, 0.0482, 0.0564} is confirmed but permuted in
      panel attribution. Fig.~5e claimed 0.0578 = Serpine1-KO vs mix
      (0.05779), not a Control contrast. Fig.~1 (Perturb-map screen)
      has no source-data file -- a documented reporting gap.
      See `methods_paper/reproduction_report.tex` Section `sec:sourcadata`.

**Files.** `serpin_repro/pdac_sc_preprocessing.py` (genotype re-run),
`serpin_repro/step3_annotate_cnv.py`, `serpin_repro/step4_myeloid.py`,
`methods_paper/audit_mouse_confound.py` (incl. subset-level functions
and `--source-data` flag),
`serpin_repro/data/Nature_SourceData/` (13 MOESM xlsx files),
`methods_paper/figs_reproduction/figG--figJ_sd_*.png` (matplotlib SD re-plots),
`methods_paper/figs_nature.R` (R 4.6.1 + ggplot2; Nature-style re-draw of
every non-deposited panel -- `fig2`, `fig3`, `fig4`, `fig5`, `fig5kl`,
`ed5_7`, `ed8` -- plus the six scRNA-seq figures `figA`--`figF` read from
`figs_reproduction/plot_tables/*.tsv`),
`methods_paper/figs_nature/*.pdf` (the thirteen redrawn vector figures;
every figure in the report is now R-generated),
`methods_paper/figs_reproduction/source_data_audit.tsv` (recomputed p table),
`methods_paper/reproduction_report.tex` (paper-style reproduction-audit
report, compiled PDF 24 pp.),
`serpin_repro/processed/GSE289138_processed.h5ad`,
`serpin_repro/processed/GSE289138_annotated.h5ad` (re-written, still has
`treatment` from the old cache; the processed object is clean).

---

## Superseded files

| superseded | replaced by |
|---|---|
| `gse202051_hum/cross_platform_consensus.tsv` | `cross_platform_consensus_CORRECTED.tsv` (this directory) |
| `gse202051_hum/processed_dsp/serpin_proximity_comparison.tsv` | `methods_paper/compartment_patient_clustered.tsv` |
| `gse202051_hum/processed/step3_summary.tsv` MARCO columns | `gse202051_hum/processed/marco_corrected_test.tsv` |
| `serpin_repro/processed/GSE289138_annotated.h5ad` `treatment` column | `serpin_repro/processed/sample_metadata_CORRECTED.tsv` (genotype-based) |
| `methods_paper/figs_reproduction/figA--figJ_*.png` (matplotlib, report figures) | `methods_paper/figs_nature/*.pdf` (R/ggplot2; the PNGs are still emitted by `audit_mouse_confound.py` as a cross-check but are no longer cited by the report) |

The old `cross_platform_consensus.tsv` states "SERPIN-high ROIs have MORE
T-cells, not fewer". That file's p = 0.0063 is not interpretable (it pools
compartments). The direction is partially supported by the corrected analysis —
SERPIN-high **epithelial** regions show +2.15 T-cell signal at
patient-clustered p = 0.034–0.048 — but the original file must not be cited:
its statistic is wrong even though its sign happens to be right.

## What is NOT changed

- **GSE202051 SERPIN result** is unaffected: 41.0 % → 55.1 % of malignant
  epithelial nuclei are SERPIN-high (naive → treated), `p = 0.065`,
  Cliff's δ = −0.338 (δ computed as naive vs treated). This was always a
  patient-level test on 43 patients (18 vs 25), not pseudo-replicated.
- **GSE300595 absence of FDR-significant DARs** is unaffected.
- **Zero patient overlap** across the three cohorts is unaffected; all
  cross-cohort statements remain consistency-only, never paired.
