# Reproduced Figure Status — Nature paper (Falcomata et al. 2026)

Run with: `Rscript methods_paper/figs_nature.R`  
Output directory: `methods_paper/figs_nature/`

---

## Data availability by panel

| Panel | Source | Status |
|-------|--------|--------|
| **Fig 1** all panels | Perturb-map / MICSSS / DepMap | ⬜ All placeholders — no source data deposited (workflow schematics only) |
| **Fig 2d** tumour mass, IgG | MOESM9 / Figure_2d | ✅ Bar + dot plot from deposited data |
| **Fig 2e** CD8+/mg flow | MOESM9 / Figure_2e | ✅ Bar + dot plot from deposited data |
| **Fig 2f** T-cell subsets | GSE289138 tcell_subsets.tsv | ✅ Grouped bar chart (classic chemistry) |
| **Fig 2h** CD8+/mm2 IHC | MOESM9 / Figure_2h | ✅ Bar + dot plot |
| **Fig 2i** GZMB+CD8+/mm2 | MOESM9 / Figure_2i | ✅ Bar + dot plot |
| **Fig 2k** tumour mass NSG | MOESM9 / Figure_2k | ✅ Bar + dot plot |
| **Fig 2l** survival IgG/aPD1 | MOESM9 / Figure_2l | ✅ KM curves, log-rank claims annotated |
| **Fig 2n** survival vehicle/PAI-039 | MOESM9 / Figure_2n | ✅ KM curves, p-values annotated |
| **Fig 2a,b,c,g,j,m** | External / drawing | ⬜ Placeholders |
| **Fig 3e,f,j,l** ECM & fibrin quantification | MOESM10 / Figure_3{e,f,j,l} | ✅ Bar + dot plots |
| **Fig 3g** CAF UMAP | plot_tables/caf_umap.tsv | ✅ UMAP |
| **Fig 3h** CAF subsets (Flex) | plot_tables/caf_subsets.tsv | ✅ Grouped bars |
| **Fig 3a,b,d,i,k** | Xenium images & heatmaps | ⬜ Xenium data on Zenodo only, not in Source Data |
| **Fig 4a** macrophage count | MOESM11 / Figure_4a | ✅ Bar + dot plot |
| **Fig 4d** spatial distance to PAI cells | MOESM11 / Figure_4d | ✅ Bar + dot plot |
| **Fig 4e** macrophage depletion tumour mass | MOESM11 / Figure_4e | ✅ Bar plot with pal_for colours |
| **Fig 4f** CD18-blockade tumour mass | MOESM11 / Figure_4f | ✅ Bar plot |
| **Fig 4g** myeloid subset fractions | plot_tables/myeloid_subsets.tsv | ✅ Grouped bars (classic chemistry) |
| **Fig 4j** BMDM MFI plastic vs fibrin | MOESM11 / Figure_4j | ✅ Grouped bar |
| **Fig 4b,c,h,i,k,l** | Raw images, per-cell scores, human spatial | ⬜ Not available in public deposits |
| **Fig 5e–h** mixing experiment | MOESM12 / Figure_5{e,f,g,h} | ✅ Bar + dot plots |
| **Fig 5k** ROI CD8 fraction boxplot | MOESM12 / Figure_5k | ✅ Boxplot + jittered points |
| **Fig 5l** control fraction vs CD8 scatter | MOESM12 / Figure_5l | ✅ Scatter + lm fit |
| **Fig 5a–d,i,j** | Human compendium / drawings | ⬜ Not deposited |
| **ED2e** Pro-Code tag enrichment | MOESM13 / Ext_Data_Figure_2e | ✅ Crossbar + jitter |
| **ED4c** clonal diversity over time | MOESM14 / Ext_Data_Figure_4c | ✅ Boxplots, d7/d14/d21 |
| **ED5c** doubling time | MOESM15 / Ext_Data_Figure_5c | ✅ Bar plot |
| **ED5d** tumour mass sgRNA validation | MOESM15 / Ext_Data_Figure_5d | ✅ Bar plot |
| **ED5e** survival sgRNA KOs | MOESM15 / Ext_Data_Figure_5e | ✅ KM curve |
| **ED5f** scRNA-seq UMAP major pops | plot_tables/umap_coords.tsv | ✅ UMAP (subset 12k pts) |
| **ED5g** T-cell subset UMAP | plot_tables/tcell_umap.tsv | ✅ UMAP |
| **ED5i** terminal exhausted subsets | plot_tables/tcell_subsets.tsv | ✅ Grouped bars |
| **ED5j** Rag2-/- tumour mass | MOESM15 / Ext_Data_Figure_5j | ✅ Bar plot |
| **ED5l** KC survival PAI-039 combo | MOESM15 / Ext_Data_Figure_5l | ✅ KM curve |
| **ED7b** PDPN+/aSMA+ (%) | MOESM16 / Ext_Data_Figure_7b | ✅ Bar plot |
| **ED7e** CD31+ endothelial (%) | MOESM16 / Ext_Data_Figure_7e | ✅ Bar plot |
| **ED7c** Flex UMAP | plot_tables/umap_coords.tsv | ✅ UMAP |
| **ED7f** endothelial subset UMAP | plot_tables/endo_umap.tsv | ✅ UMAP |
| **ED7g** endothelial subsets | plot_tables/endo_subsets.tsv | ✅ Grouped bars |
| **ED8a** neutrophils/mg | MOESM17 / Ext_Data_Figure_8a | ✅ Bar plot |
| **ED8h** macrophages/mg | MOESM17 / Ext_Data_Figure_8h | ✅ Bar plot |
| **ED8j–m,o,r** various flow/IHC quantifications | MOESM17 | ✅ Bar plots |
| **ED8p** macrophage subset UMAP | plot_tables/myeloid_umap.tsv | ✅ UMAP |
| **ED1, ED3, ED6, ED9, ED10** | External / schematic | ⬜ All placeholders |

---

## Legend for symbols

- ✅ = redrawn from deposited Source Data spreadsheet or GSE289138 plot tables
- ⬜ = placeholder grey box with reason stated (no underlying data accessible)
- Note: Some Xenium/spatial-deconvolution outputs are on Zenodo but NOT in Source Data spreadsheets

---

## Supplementary panels (not part of Nature manuscript)

| File | Content |
|------|---------|
| `figS1_umap.pdf` | GSE289138 UMAP: cell type / genotype / chemistry (3 panels side-by-side) |
| `figS2_serpin.pdf` | SERPINE1/SERPINB2 detection % in malignant epithelium by genotype × chemistry |
| `figS3_effects.pdf` | Genotype-vs-chemistry effect-size ratio per readout |

---

## How to re-run

```bash
cd E:/pancreatic_cancer/submit
"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" methods_paper/figs_nature.R
```

All output goes to `methods_paper/figs_nature/*.pdf`.
