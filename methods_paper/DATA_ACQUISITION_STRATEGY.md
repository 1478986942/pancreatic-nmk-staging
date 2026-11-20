# Data Acquisition Strategy — What We Can Actually Get

## Zenodo Deposits (Confirmed Available)

### Record: https://zenodo.org/records/18882157
Title: "A serpin-myeloid axis in pancreatic cancer heterogeneity and immune evasion"

| File | Size | Contains | Panels Filled |
|------|------|----------|---------------|
| `Perturb-map_MM_xenium.zip` | 43 GB | Xenium spatial transcriptomics + t-CyCIF integration | Fig 3b,c,d; ED6 a-f |
| `BDBAT0245.ome.tif` | 18 GB | Multiplex IF microscopy images | Fig 1b, ED2c, ED7a,d |

**Problem**: Too large for automated download. Need institutional access or manual download.

---

## Alternative Data Sources to Try

### 1. DepMap CRISPR Dependency Scores (for ED1b)
- **URL**: `https://depmap.org/portal/download/`
- **File**: `CRISPR(Gene Dependency)public25Q3.csv.gz` (~50 MB)
- **What it gives**: Per-gene dependency scores across 46 PDAC cell lines
- **Panel filled**: ED1b (dependency distribution histogram)
- **Action**: Downloadable via Python requests

### 2. TCGA/GTEx Pan-cancer Expression (for Fig 2a replacement)
- **UCSC Xena**: `https://ucsc-xena.org/geneAndExpression`
- **GDC API**: `https://api.gdc.cancer.gov/`
- **TCGAbiolinks R package**: Handles authentication internally
- **Panel filled**: Fig 2a lollipop plot
- **Status**: Already have Supplement A with published values; can improve with live data

### 3. JAVELIN Renal 101 Survival Data (for ED5m replacement)
- **Source**: Motzer et al. Nat Med 2020;26:1733-1741
- **Data**: Published KM curves with HR=1.48, p=0.004 for SERPINE1 high vs low
- **Panel filled**: ED5m
- **Status**: Already have Supplement C approximation

### 4. Human Spatial Transcriptomics (Fig 4k/l, ED9)
- **GSE274557, GSE278694**: Visium PDAC spatial data
- **Access**: Open on GEO but requires cell2location reference model
- **Alternative**: Use published cell-type fractions from paper directly
- **Panel filled**: Fig 4k,l; ED9 a,b

---

## Recommended Priority Actions

### Tier 1: Quick Wins (can do now)
1. **Download DepMap scores** → Fill ED1b
2. **Hand-draw all ⬜ schematic panels** using Inkscape/BioRender
3. **Request Xenium zip from authors** via email to brian.brown@mssm.edu

### Tier 2: Medium Effort
4. **Use TCGAbiolinks** to get fresh pan-cancer DE (replaces Supplement A)
5. **Extract published HR values** from Supplementary Tables of cited papers
6. **Create simplified spatial maps** using available Visium count matrices

### Tier 3: Hard Blockers
7. **dbGap controlled access** for phs002371 — requires IRB application
8. **Raw microscopy images** — only available from authors

---

## Implementation Plan

```bash
# Step 1: Get DepMap (immediate)
python -c "import requests; r=requests.get('https://depmap.org/portal/download/data/public/25Q3/CRISPR/GeneDependency/CRISPR(Gene Dependency)public25Q3.csv.gz', stream=True); open('depmap.csv.gz','wb').write(r.content)"

# Step 2: Extract Xenium metadata from Zenodo (requires auth or large bandwidth)
# Contact: Brian Brown Lab, Icahn School of Medicine at Mount Sinai
# Email: brian.brown@mssm.edu
# Request: Perturb-map MM Xenium analysis files (UMAP coordinates, gene expression matrices)

# Step 3: Hand-draw schematics
# Use figma.com or inkscape.org for workflow diagrams
```

---

## Current Public Repository Status

✅ **GitHub**: https://github.com/1478986942/pancreatic-nmk-staging (public)
- 211 files including 49 PDF figures
- Full methods paper draft
- Complete reproduction audit

✅ **Gitee**: https://gitee.com/msdyy/projects_information-biology (needs visibility check)

---

## Next Immediate Action Items

1. [ ] Download and process DepMap data for ED1b
2. [ ] Draft hand-drawn schematics for all ⬜ panels
3. [ ] Email Brian Brown lab requesting Xenium object extraction
4. [ ] Create Supplement F with DepMap dependency results
5. [ ] Update FULL_REPRODUCTION_AUDIT.md with new coverage
