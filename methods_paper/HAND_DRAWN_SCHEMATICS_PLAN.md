# Hand-Drawn Schematics — Action Plan

## Panels Requiring Hand-Drawn Diagrams

These panels are schematic/workflow illustrations in the original paper. We can recreate them from the methods section text.

### Fig 1a — Perturb-map Workflow
**Original**: KPC cells → PC/CRISPR library transduction → orthotopic implantation → imaging → analysis pipeline  
**Recreate from**: Methods "Lentiviral vector construction and production" + "Perturb-map workflow" sections

### Fig 1g — Rag2−/− Experimental Design  
**Original**: KPC_PC/CRISPR cells → orthotopic transplant into Rag2−/− → tumour harvest day 14  
**Recreate from**: Methods "Mouse experiments" section

### Fig 2c — In Vivo Validation Design
**Original**: Control/Serpinb2-KO/Serpine1-KO KPC cells → orthotopic implantation → flow/IHC/survival analysis  
**Recreate from**: Methods "Cell culture" + "Mouse experiments" sections

### Fig 2j — Rag2−/− and NSG Design
**Original**: Serpin KO cells → transplantation into immunodeficient hosts → tumour growth comparison  
**Recreate from**: Same as above

### Fig 2m — PAI-039 ± anti-PD-1 Design
**Original**: KC/KPC tumours → treatment with Vehicle/PAI-039/anti-PD-1/combo → survival  
**Recreate from**: Methods "For ICB" and "For PAI1 inhibition studies" sections

### Fig 3a — Perturb-map Multi-modal Workflow
**Original**: Perturb-map → Xenium ST → t-CyCIF → integration → spatial DE  
**Recreate from**: Methods "Xenium spatial transcriptomics" + "Post-Xenium cyclic immunofluorescence"

### Fig 4i — BMDM Fibrin Culture Schematic
**Original**: BMDM on plastic vs fibrin-coated plates → conditioned medium → flow cytometry  
**Recreate from**: Methods "Fibrin co-culture assays"

### ED1a — Gene Prioritisation Workflow
**Original**: scRNA-seq → extracellular annotation → DepMap filtering → pathway enrichment → final gene list  
**Recreate from**: Methods "Selection and prioritization of extracellular tumour-derived genes"

### ED2a,b — Pro-Code Vector & Validation Workflow
**Original**: Pro-Code vector structure → library cloning → sorting → in vivo validation  
**Recreate from**: Methods "Lentiviral vector construction"

### ED3a,d — Perturb-map Workflow & Image Analysis
**Original**: Full Perturb-map pipeline → image processing → deconvolution → spatial analysis  
**Recreate from**: Methods sections on MICSSS + CyCIF

### ED5k — Anti-PD-1 + Serpin KO Design
**Original**: Serpin KO tumours → anti-PD-1 treatment → survival comparison  
**Recreate from**: Fig 2l legend + Methods

### ED6a,b — Multi-modal Workflow
**Original**: Perturb-map → Xenium → t-CyCIF registration → integration  
**Recreate from**: Methods "Xenium-CyCIF integration"

### ED7a,d — t-CyCIF Images
**Original**: PDPN/αSMA/PDGFRα staining → quantification; CD31 staining  
**Recreate from**: Methods "CyCIF" section

### ED8g,i,n — Treatment Design Schematics
**Original**: Macrophage depletion protocol; CD18 blockade protocol; combo treatment schedule  
**Recreate from**: Methods "For macrophage depletion" and "For CD18-blocking"

### ED10i — Proposed Model Schematic
**Original**: SERPIN → fibrin → αMβ2 integrin → macrophage retention → immunosuppression  
**Recreate from**: Discussion section paragraph 1-2

---

## Recommended Tools

| Tool | Cost | Best For |
|------|------|----------|
| **BioRender** | Free tier (60 figures) | Biological pathways, cell diagrams |
| **Inkscape** | Free | Technical schematics, workflow diagrams |
| **PowerPoint** | Included | Simple box-and-arrow flows |
| **draw.io** | Free | Process diagrams |

---

## Priority Order

1. **Fig 1a** (main workflow) — most cited figure element
2. **Fig 2c/j/m** (experimental designs) — critical for Results narrative
3. **ED1a** (prioritisation workflow) — establishes screening logic
4. **All others** — supplementary context

---

## Template Structure for Each Schematic

```
[Box: Cell type/condition] → [Arrow: procedure] → [Box: readout/measurement]
     ↓                                            ↓
  "KPC cells"                              "Flow cytometry"
  "Serpinb2-KO"                            "CD8+ T cell quantification"
  
Scale bar / Legend: n=x mice per group
Statistical test: Two-sided Welch's t-test
```

---

## Next Steps

1. Create BioRender account (free) or install Inkscape
2. Draft Fig 1a using Methods section text as reference
3. Batch-create all remaining schematics in same style
4. Export as PDF at Nature specifications (300 dpi minimum)
5. Replace grey placeholder boxes in figs_nature.R with actual figures

Estimated time: 4-6 hours total for all ~15 schematics.
