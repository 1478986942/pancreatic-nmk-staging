# Pancreatic cancer SERPIN–MARCO audit

This repository contains the reproducible **audit and correction** package for the
SERPIN–MARCO / immune-checkpoint-blockade analysis project. It is deliberately
not a repository of new biological claims: the current phase is metadata,
metadata-provenance, unit-of-analysis, panel-coverage, denominator and
confounding audits.

## Current scientific status

- **GSE289138** contains three genotypes (`control KO`, `Serpinb2 KO`,
  `Serpine1 KO`) crossed with two library chemistries (`classic`, `Flex`), with
  **one tumour per genotype and no treatment arm**. Cell-level hypothesis tests
  and downstream mechanism analyses are therefore prohibited; only descriptive
  genotype/chemistry readouts are admissible.
- **GSE202051** weakly supports persistence of the SERPIN-positive malignant
  programme after chemotherapy (41.0% to 55.1%, nominal `p = 0.065`), but does
  not support preservation of the MARCO-positive macrophage fraction after
  correcting the denominator (`p = 0.758`).
- **GSE199102** requires compartment-stratified, patient-clustered inference.
  The pooled spatial result is invalid; the corrected epithelial T-cell signal
  is borderline and does not survive multiplicity correction.
- **GSE300595** yields zero FDR-significant DARs and no SERPIN response
  association (`p = 1.000`). The ATAC-to-spatial bridge tracks myeloid abundance
  and does not identify a response-specific programme.
- The pre-specified **IMvigor210** baseline bulk test is statistically null for
  `SERPINE1`, `SPP1` and `MARCO`; no bulk-predictive claim is supported.
- No public cohort currently permits cell-level testing of a PDAC
  SERPIN-positive / myeloid niche with an ICB-response endpoint.

## Read first

1. [`AGENTS.md`](AGENTS.md) — standing project constraints and analysis freeze.
2. [`CORRECTIONS.md`](CORRECTIONS.md) — C1–C6 error catalogue and invariants.
3. [`LIMITATIONS.md`](LIMITATIONS.md) — design and identifiability limits.
4. [`RESULTS_RECONCILED.md`](RESULTS_RECONCILED.md) — only the corrected,
   currently supported conclusions.
5. [`methods_paper/MANUSCRIPT_A_RESULTS.md`](methods_paper/MANUSCRIPT_A_RESULTS.md)
   — methods-paper draft.
6. [`methods_paper/ICB_COHORT_AUDIT.md`](methods_paper/ICB_COHORT_AUDIT.md)
   — cohort availability and response-label audit.
7. [`methods_paper/reproduction_report.pdf`](methods_paper/reproduction_report.pdf)
   — GSE289138 reconstruction and source-data audit report.

## Repository layout

| Path | Role |
|---|---|
| `methods_paper/` | audit scripts, result tables, manuscript draft, figures and reproduction report |
| `cohort_audit/` | cohort availability workspace documentation |
| `gse202051_hum/` | GSE202051 and GSE199102 analysis scripts plus corrected summary tables |
| `gse300595/` | GSE300595 analysis scripts plus patient-level audit summaries |
| `serpin_repro/` | GSE289138 preprocessing/audit scripts, corrected metadata and quarantine marker |

## Data and reproducibility policy

Raw matrices, processed H5AD objects, virtual environments, nested upstream
repositories and unrelated training artefacts are intentionally excluded from
Git. The included summary tables are the patient-level or descriptive outputs
needed to audit the reported numbers. Public source accessions are:

- GSE202051
- GSE199102
- GSE300595
- GSE289138
- IMvigor210 / Mariathasan 2018 (ICB-Portal)

The local analysis environment used for the mouse workflow is documented in
`AGENTS.md`. Do not commit access tokens, passwords, raw patient data or
credentials.

## Provenance note

The GSE289138 analysis is a reconstruction rather than an exact replay: the
upstream source repository deposits figure-level notebooks but not the
processing pipeline that generated its intermediate H5AD files. The chemistry
effect is therefore quantified alongside the genotype effect before any
directional interpretation.
