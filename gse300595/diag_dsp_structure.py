import os, gzip
import numpy as np
import pandas as pd

DSP = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"
seg = pd.read_csv(gzip.open(os.path.join(DSP, "segment_props.txt.gz"), "rt"), sep="\t")
cnt = pd.read_csv(gzip.open(os.path.join(DSP, "q3norm_counts.txt.gz"), "rt"),
                  sep="\t", index_col=0)
panel = set(cnt.index)

print("segments:", len(seg))
print("unique (Slide_name, ROI_number):",
      seg.groupby(["Slide_name", "ROI_number"]).ngroups)
print("unique (Patient, ROI_number):",
      seg.groupby(["Patient", "ROI_number"]).ngroups)
print("\nsegments per physical ROI (top):")
print(seg.groupby(["Slide_name", "ROI_number"]).size().value_counts().to_dict())
print("\nunique Patient:", seg["Patient"].nunique())
print("\nSegment x Patient (head):")
print(seg.groupby(["Patient", "Segment"]).size().unstack(fill_value=0).head(8).to_string())

markers = {
 "epithelial": ["KRT8","KRT18","KRT7","CEACAM6","S100P","MUC1","TFF1","MSLN","CLDN4","LGALS4","TSPAN8","KRT17","KRT5"],
 "CAF":        ["COL1A1","COL1A2","COL3A1","DCN","LUM","FAP","FBLN1","PDGFRB","ACTA2","TAGLN","POSTN","COL5A1"],
 "myeloid":    ["CD68","CD163","CSF1R","LYZ","C1QA","C1QB","CD14","ITGAM","TYROBP","FCER1G","AIF1","CD4"],
 "Tcell":      ["CD3D","CD3E","CD2","TRBC1","TRBC2","IL7R","NKG7","KLRD1","PRF1","CD8B","CD247"],
 "endo":       ["VWF","CDH5","CLDN5","ERG","PTPRB","KDR","ENG"],
}
for k, gs in markers.items():
    have = [g for g in gs if g in panel]
    miss = [g for g in gs if g not in panel]
    print(f"\n{k}: have={have}\n   missing={miss}")
