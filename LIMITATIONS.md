# Limitations

Consolidated limitations for the three-cohort analysis. Every item below is a
real constraint on what the data can support; none is a formality.

---

## 1. Causal inference is impossible by design

- **GSE300595**: all 12 patients are **post-treatment** resection specimens.
  There are **no matched pre-treatment biopsies** in the series. Every
  comparison is therefore *responder vs non-responder at one time point*, not
  *before vs after chemotherapy*. Any statement of the form "chemotherapy
  changes X" is unavailable from this dataset; the only supportable form is
  "X differs between responders and non-responders after treatment".
- **GSE199102**: naive vs treated comparison is between patients, not within
  patients.
- **GSE202051**: naive vs treated is also between patients (18 naive, 25
  treated), with no within-patient pairing.

## 2. Zero patient overlap across cohorts

The three cohorts share **no patients** (verified). Consequences:

- No paired snRNA ↔ DSP analysis is possible; this was formally abandoned.
- No cross-platform sample matching, no shared patient-level covariates.
- All cross-cohort results are **consistency checks**, which can fail to
  replicate for platform or cohort reasons unrelated to biology.
- Nothing here can be meta-analysed as if the cohorts were independent samples
  of one population.

## 3. The three platforms do not measure the same SERPIN axis

| platform | SERPIN genes available and used |
|---|---|
| GSE202051 snRNA | `SERPINE1` + `SERPINB2`; SERPINB2 is ~undetectable, so effectively `SERPINE1` |
| GSE300595 snMultiome | `SERPINE1` only; `SERPINB2` detected in 11 / 35,115 nuclei (0.03 %) |
| GSE199102 GeoMx DSP | **`SERPINB2` only**; `SERPINE1` is not on the 9,009-gene panel |

"Increased/decreased SERPIN signal" therefore means different genes on
different platforms. This is not fixable post hoc and must be stated whenever
platforms are compared.

## 4. GSE300595: fragment-level QC is impossible

GEO distributes **only** `filtered_feature_bc_matrix.h5` for this series. There
is no `fragments.tsv.gz`, no BAM, and no `RAW.tar`. Therefore these four
metrics from the source paper **cannot be computed**:

- `nucleosome_signal`
- `TSS.enrichment`
- `blacklist_fraction`
- `pct_reads_in_peaks`

Only matrix-derivable metrics were applied: `nCount_RNA`, `nFeature_RNA`,
`percent.mt`, `percent.ribo`, `nCount_ATAC`, `nFeature_ATAC`, `percent.mt.atac`.

**Consequence.** 35,115 nuclei were retained vs 30,840 in the paper. The gap is
expected — the uncomputable filters remove cells on the basis of fragment
geometry, which is invisible in a count matrix. Per-sample pass rates ranged
41–86 %. The two datasets are therefore **not** the same cell set, and cell
counts should not be compared to the paper directly.

## 5. GSE300595: no peak calling, no motif analysis, no chromVAR

- With no fragments, de novo peak calling is impossible. The peak space is the
  **union of the 12 cellranger-arc peak sets** (293,174 peaks, median width
  975 bp, 297.5 Mb), filtered to peaks detected in ≥5 % of myeloid nuclei
  (21,769 peaks).
- `snapatac2` could not be installed (the `macs3` wheel fails to build on
  Windows), but this is moot — it requires fragments that do not exist here.
- **No genome FASTA was obtained**, so motif enrichment, TF footprinting and
  chromVAR are all unavailable. Chromatin results therefore cannot be linked to
  transcription-factor activity mechanistically. DARs are annotated only to the
  nearest protein-coding gene.
- LSI on binarized TF-IDF peaks explains little variance (component 0 = 0.019),
  which is normal for sparse binary peak data but limits cluster resolution.

## 6. GSE300595: statistics are severely underpowered

- Only **5 responders vs 7 non-responders**.
- The minimum attainable two-sided Mann-Whitney p with 5 vs 7 is **0.002525**.
  Several DAR peaks sit exactly at this floor, i.e. they are the most extreme
  rank separation possible and cannot be more significant regardless of effect
  size.
- **No peak survives FDR** (0 at padj < 0.25 and 0 at padj < 0.1 among 18,671
  tested peaks). The same is true of the pseudobulk RNA DE (0 genes at
  padj < 0.1).
- All GSE300595 results are **exploratory / hypothesis-generating**. The
  apparent DAR "programmes" are effect-size rankings with an unadjusted p < 0.05
  filter, not discoveries.

## 7. GSE300595: myeloid abundance is extremely sample-skewed

Myeloid nuclei per sample range from 21 (pc149) to 563 (pc131), and myeloid
fraction from 0.7 % (pc123) to 24.2 % (pc102). The responders happen to include
the two highest-myeloid samples. This is not balanced by design and cannot be
balanced post hoc; patient-level medians are reported precisely to avoid
letting one sample dominate, but the imbalance remains a confound for any
abundance-dependent measure.

## 8. GSE300595: myeloid pseudobulk is contaminated by ambient epithelial RNA

Responder-elevated genes in the myeloid pseudobulk include `MUC6`, `HNF1B`,
`PROM1` and `SERPINE1` — epithelial products. This is ambient contamination
typical of snRNA-seq from solid tumour, plus doublets. It means:

- the myeloid "responder" signature is partly epithelial signal;
- `SERPINE1` in a myeloid pseudobulk cannot be read as macrophage SERPINE1;
- DE directions for epithelial genes in the myeloid compartment are not
  interpretable as myeloid biology.

## 9. GSE300595: SERPINB2 and MARCO are effectively unmeasurable

- `SERPINB2`: 11 / 35,115 nuclei (0.03 %). The `SERPIN_score` is therefore
  driven entirely by `SERPINE1`, despite the gene list matching GSE202051.
- `MARCO`: detected in 0.57 % of nuclei. Any MARCO-based statement in
  GSE300595 lacks power and should be treated as qualitative.
- `M1_like` module scores are ≈ 0 across all myeloid cells, so the "M1" state is
  nominal and weakly separated; state labels rest mainly on LAM / M2 / DC axes.

## 10. GSE300595: myeloid states were assigned by cross-cluster z-score

Absolute module scores were degenerate (the `M2_like` module was high in every
subcluster, because `CD163`/`MRC1` are broadly expressed by macrophages here).
States were therefore assigned by z-scoring each module **across subclusters**
and taking the most distinguishing module, with a threshold of z > 0.5; clusters
below threshold were labelled `Macrophage_other`. This yields relative, not
absolute, states: a cluster labelled `LAM` means "most LAM-like among these
clusters", not "LAM by an external standard".

## 11. GeoMx DSP: severe panel gaps for the myeloid programme

The panel is a whole-transcriptome assay (9,009 genes) but is missing most genes
that define the LAM/MARCO niche of interest:

| gene | on panel |
|---|---|
| `SERPINE1` | **no** |
| `MARCO` | **no** |
| `SPP1` | **no** |
| `MRC1` | **no** |
| `GPNMB` | **no** |
| `TREM2` | **no** |
| `LPL` | **no** |
| `PPARG` | **no** |
| `CD8A` | **no** |
| `GNLY` | **no** |
| `PECAM1`, `CD31` | **no** |
| `EPCAM`, `KRT19` | **no** |
| `CD163`, `CD68` | yes |
| `CD3D`, `CD3E`, `NKG7` | yes |
| `CD74` | **no** |

Consequences: the myeloid score is `CD163`+`CD68` only; the T-cell score is
`CD3D`+`CD3E`+`NKG7` (no CD8A, no GNLY, so CD8 vs NK cannot be separated); the
LAM programme can only be approximated by `APOE`, `APOC1`, `CLIP4`, `CD163`,
`C1QA`, `C1QB`, `LYZ`, `CSF1R`, `CTSD`; epithelial content must be estimated
from `KRT8/KRT18/MUC1/MSLN/CLDN4/LGALS4/TSPAN8/KRT17`.

## 12. GeoMx DSP: the correct unit of analysis is neither segments nor ROIs

194 physical ROIs were split into ≤1 segment per compartment (Immune /
Epithelial / CAF). Two distinct mistakes are possible here, and both were made
at different points:

- **Pooling compartments** (570 segments at once) mixes three tissues with
  different baselines. This is the original error and it is a confounding error,
  not mainly a clustering error: ICC of the T-cell score across segments is
  **0.000** (design effect 1.00) and segment-level type I error is 6.8 %, so
  treating segments as independent is only mildly liberal.
- **Averaging to one value per ROI** injects the between-compartment variance
  into every observation and destroys compartment-specific effects (T-cell
  p went from 0.0063 to 0.685).

The defensible unit is **one observation per ROI within a compartment, with
inference clustered at patient level** (22 patients, 7–15 ROIs each). That gives
T-cell p = 0.034 in epithelium. See `CORRECTIONS.md` C1 and
`methods_paper/compartment_patient_clustered.tsv`.

## 13. GeoMx DSP: histology labels are missing

All 570 matched segments have **no histology annotation** (589 / 658 rows in
`segment_props` are NaN; only 69 are labelled). The available `Segment` field
(`Immune` / `Epithelial` / `CAF`) is a *geometric compartment*, chosen by the
assay, not a pathological classification. Tumour vs stroma vs normal duct
cannot be distinguished.

## 14. The DSP "bridge" is non-specific

The decisive limitation of the GSE300595 → DSP bridge: the non-responder and
responder accessibility programmes are **not distinguishable** in DSP.

| comparison | rho |
|---|---|
| NR programme vs R programme, segment level | 0.883 |
| NR programme vs R programme, physical-ROI level | 0.895 |

Both programmes score higher in SERPIN-high regions, with similar or larger
effects for the **responder** programme. Partial correlation controlling for
epithelial/CAF/myeloid content does not change the NR association
(rho 0.187 → 0.189), because the shared axis is intrinsic to both gene sets,
not a compositional confound. At physical-ROI level the NR association also
loses significance within treatment strata (naive p = 0.14, treated p = 0.44).

**Therefore the DSP data cannot test the direction of the myeloid chromatin
programme.** Any claim that the NR accessibility programme is specifically
enriched in the SERPIN-high niche is unsupported.

## 15. Patient-level ATAC score measures myeloid abundance, not SERPIN status

Across the 12 patients, the NR accessibility score does not correlate with
SERPIN fraction (`rho = -0.014, p = 0.97`) but correlates strongly with myeloid
fraction (`rho = 0.811, p = 0.0014`). The score is therefore largely a
cellularity/depth proxy, and the responder vs non-responder separation in
`NR_minus_R` (p = 0.0025) is partly circular — those peaks were selected for
that contrast.

## 16. Gene assignment by nearest TSS is approximate

DARs were annotated to the nearest protein-coding gene TSS, ignoring strand,
enhancer-gene looping and orientation. Most peaks are acceptable (349 / 370 NR
peaks lie within 100 kb of a TSS), but some assignments are implausible
(e.g. a chr17 peak assigned to `KCNJ18` at 829 kb). Distant assignments add
noise to every gene-level programme; a ≤100 kb version (`NR_near100kb`) is
reported alongside and gives the same conclusions.

## 17. No FDR-significant result survives anywhere in this project

Summarised for clarity:

| analysis | test | significant after FDR? |
|---|---|---|
| GSE202051 SERPIN fraction | patient-level MWU, p = 0.065 | no |
| GSE202051 MARCO fraction (corrected) | patient-level MWU, p = 0.758 | no |
| GSE300595 myeloid pseudobulk DE | 0 genes padj < 0.1 | no |
| GSE300595 myeloid DAR | 0 peaks padj < 0.25 | no |
| GSE300595 LAM fraction R vs NR | p = 0.202 | no |
| IMvigor210 SERPINE1 baseline | MWU, p = 0.377, q = 0.566 | no |
| IMvigor210 SPP1 baseline | MWU, p = 0.582, q = 0.582 | no |
| IMvigor210 MARCO baseline | MWU, p = 0.201, q = 0.566 | no |
| IMvigor210 three-gene z-score | MWU, p = 0.209 (exploratory) | no |

The project's conclusions are therefore **directional and hypothesis-generating**.
They should be presented as such, and none should be described as a demonstrated
effect.

## 18. Public bulk matrices may not be self-describing

ICB-Portal's `Mariathasan_2018_expression.tsv` has **no gene-identifier column**:
18,021 rows × 348 columns, first column a sample ID. The conventional
`pd.read_csv(..., index_col=0)` therefore yields a 348 × 347 *sample* matrix, and
every gene lookup silently fails or, worse, a positional lookup analyses the
wrong genes. Gene names exist only in `Mariathasan_2018.Rdata` (`TPM` `dim_0`),
which also shows 18,021 rows but 18,020 unique symbols (one duplicate, averaged).
Any ICB-Portal expression TSV must be proved equal to its authoritative container
(`np.allclose`, matching sample order) before names are attached. See
`CORRECTIONS.md` C5.

## 19. Bulk expression cannot stand in for spatial co-occurrence

The IMvigor210 pilot (68 R / 230 NR, baseline only) returned a clean negative:
`SERPINE1` p = 0.377, `SPP1` p = 0.582, `MARCO` p = 0.201, BH q ≥ 0.566; the
three-gene z-score p = 0.209. All effects point in the predicted direction
(higher in non-responders) but |Cliff's δ| is only 0.04–0.10, and with 68 vs 230
the cohort can resolve moderate effects (|δ| ≈ 0.2) but not small ones. The three
genes do co-vary (ρ = 0.40–0.53), so the co-expression premise survives; the leap
from co-expression to response association does not.

Consequences that constrain the manuscript:

- No claim that bulk SERPIN/SPP1/MARCO expression predicts ICB resistance.
- No claim that the spatial hypothesis is refuted; bulk TPM is the wrong
  instrument for a co-occurrence claim, and this is one cancer type and one drug.
- No escalation to other cohorts with gene-by-gene testing. A future test must
  use a pre-specified composite in an independent cohort with the same endpoint,
  and must report the result whatever it is.
- The SERPINE1 **mean** is higher in responders (389.07 vs 372.96) while its
  median and Cliff's δ are higher in non-responders; the mean is tail-driven and
  must not be reported as a direction.

## 20. Software / environment constraints

- `snapatac2` unavailable (macs3 build failure on Windows); irrelevant here but
  means no ArchR/Signac-equivalent workflow was possible.
- `harmonypy` pinned to 0.0.9 (pure Python); 0.0.10 wheel build fails.
- Analyses run on Python 3.10 with scanpy 1.11.5; results are reproducible from
  the scripts in `gse300595/` and `gse202051_hum/`.
