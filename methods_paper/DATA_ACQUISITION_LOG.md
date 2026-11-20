# Data Acquisition Log — What's Missing and How to Get It

**Date**: 2026-11-20  
**Goal**: Fill remaining ~34 unreproducible panels with publicly accessible data or substitutes

---

## Gap Analysis (by panel)

### Fig 1 — All 8 panels missing
| Panel | Needed | Source Attempted | Status |
|-------|--------|------------------|--------|
| a | Perturb-map workflow schematic | Drawing needed | ⬜ No data, needs hand-drawing |
| b | MICSSS multiplex images d7/d14/d21 | Nature source repo Zenodo | ❌ Not in deposited files |
| c | Spatial heterogeneity maps | Per-cell xy coordinates | ❌ In Zenodo notebook only |
| d | Pro-Code focal region composition | Per-cell barcode counts | ❌ Same |
| e | Fitness heatmap log2FC | Per-gene in vivo/in vitro values | ❌ Not in MOESM9–17 |
| f | Neighbourhood enrichment z-scores | Squidpy output | ❌ Computed on author servers |
| g | Rag2−/− design schematic | Drawing | ⬜ Needs hand-drawing |
| h | KO abundance in Rag2−/− | Per-mouse barcode counts | ❌ Not deposited |

### Fig 2 — Panels a, b, g, j, m missing
| Panel | Needed | Source | Status |
|-------|--------|--------|--------|
| **a** | TCGA/GTEx pan-cancer lollipop | UCSC Xena + GTEX | ✅ Supplement A created |
| **b** | Forest plot HR meta-analysis | Published cohort HRs | ✅ Supplement B created |
| g | IHC representative images | Author microscope images | ❌ Never deposited |
| j | Rag2−/− + NSG design schematic | Drawing | ⬜ Needs hand-drawing |
| m | PAI-039 ± anti-PD-1 design | Drawing | ⬜ Needs hand-drawing |

### Fig 3 — Panels a, b, c, d, i, k missing
| Panel | Needed | Source | Status |
|-------|--------|--------|--------|
| a | Multi-modal workflow schematic | Drawing | ⬜ Needs hand-drawing |
| b | Xenium + t-CyCIF region image | 10x Genomics raw object | ❌ On Zenodo as .tar.gz |
| c | Pathway activity heatmap | Xenium DGE per cell | ❌ Zenodo only |
| d | Volcano plots Xenium DGE | Same | ❌ Zenodo only |
| i | Masson trichrome H&E images | Histology slides | ❌ Not deposited |
| k | Fibrin(ogen) IHC images | Same | ❌ Not deposited |

### Fig 4 — Panels b, c, h, k, l missing
| Panel | Needed | Source | Status |
|-------|--------|--------|--------|
| b | Human PDAC MICSSS image | Author IF images | ❌ Not deposited |
| c | Spatial proximity density plots | Per-cell distances | ❌ Requires Xenium coords |
| h | Arg1 / fibrinogen-binding UMAP | Per-cell gene expression | ❌ Zenodo raw object |
| k | Human spatial enrichment effect sizes | Cell2location deconvolution | ❌ Needs reference atlas |
| l | Human spatial enrichment heatmap | Same | ❌ Same |

### Fig 5 — Panels a, b, c, d, i, j missing
| Panel | Needed | Source | Status |
|-------|--------|--------|--------|
| a | SERPIN fraction per patient (n=71) | 5 human scRNA compendium | ❌ Includes dbGap data |
| b | Venn diagram shared DE genes | DEG tables | ❌ Not deposited |
| c | Pathway enrichment serpin+ cells | GSEA results | ❌ Not deposited |
| d | Clonal mixing schematic | Drawing | ⬜ Needs hand-drawing |
| i | Injected vs recovered fraction | Per-mouse recovery counts | ❌ Not in Source Data |
| j | Multiplex IF representative images | Raw imaging | ❌ Not deposited |

### ED Figs gaps summary
| Figure | Missing panels | Primary blocker |
|--------|---------------|-----------------|
| ED1 | a-e (all except a) | DepMap raw scores not deposited; circular heatmap needs scRNA input |
| ED2 | a-d, f | Per-PC counts not deposited; image not deposited |
| ED3 | a,b,c | Perturb-map raw counts not deposited |
| ED4 | a,b,d,e,f,g | Per-cell coords + Squidpy outputs not deposited |
| ED5 | a,b,h,m | WB blots in Supplementary Data 3 (separate download); TCGA done in supplement A; JAVELIN done in supplement C |
| ED6 | a-f | All Xenium-related — Zenodo raw object required |
| ED7 | a,d | IF images not deposited |
| ED8 | b,c,e,f,q | Per-cell coords + raw images |
| ED9 | a,b | Visium spatial data restricted |
| ED10 | a-h | Human scRNA compendium with dbGap restrictions |

---

## Action Plan — What Can Be Done Now

### Phase 1: Immediate (no new downloads needed)
1. **Hand-draw schematics** for all ⬜ panels (Fig 1a/g, Fig 2c/j/m, Fig 3a, Fig 4i, Fig 5d, ED series)
   - Tools: PowerPoint, Inkscape, or BioRender free tier
   - Simple flow diagrams using the methods section text

2. **Add Supplement F**: Download Zenodo Xenium object to extract Fig 3b/c/d and ED6 d-f
   - Zenodo deposit: https://zenodo.org/records/18882157
   - File format: likely AnnData .h5ad or similar
   - Need R/Python environment to read it

### Phase 2: External data retrieval
3. **Supplement G**: TCGA survival data for ED5m (JAVELIN Renal 101) — already done
4. **Supplement H**: DepMap CRISPR dependency scores
   - URL: https://depmap.org/portal/download/
   - File: `CRISPR(Gene Dependency)public25Q3.csv.zip`
   - Can recreate ED1b

5. **Supplement I**: Pan-cancer survival from TCGA/GTEx directly
   - Already attempted via UCSC Xena — failed due to SSL
   - Alternative: use `TCGAbiolinks` R package or download from GDC portal

### Phase 3: Hard blockers (require author cooperation)
6. **MICSSS/Xenium raw images** — Contact Brian Brown lab for deposition
7. **Per-cell coordinates** — Available in Zenodo but need to parse Xenium output format
8. **Human scRNA compendium** — dbGap controlled access; requires IRB application

---

## Quick Win: Zenodo Xenium Object

The Nature paper deposits Xenium data on Zenodo:
```
https://zenodo.org/records/18882157
```

This likely contains:
- Fig 3b: Representative Xenium+t-CyCIF region
- Fig 3c: Pathway activity heatmap  
- Fig 3d: Volcano plots
- ED6 d-f: UMAP visualizations

Action: Download and inspect the Zenodo archive contents.

---

## Quick Win: DepMap Direct Download

```r
# R code to fetch DepMap dependency scores
download.file("https:// depmap.org/portal/download/data/public/25Q3/CRISPR/GeneDependency/CRISPR(Gene Dependency)public25Q3.csv.gz",
              "depmap_dependency.csv.gz", mode="wb")
dep <- read.csv(gzfile("depmap_dependency.csv.gz"))
# Filter for PDAC cell lines (46 lines from paper)
```

This would recreate ED1b directly from source.

---

## Next Steps Decision Points

1. **Hand-draw schematics** — Quick win, ~2 hours work
2. **Download Zenodo Xenium object** — Medium effort, unlocks 6 panels
3. **Fetch DepMap scores** — Easy, fills ED1b
4. **Contact authors for images** — Hardest path, may never happen

Recommend starting with #1 and #3 (quick wins), then tackle #2.
