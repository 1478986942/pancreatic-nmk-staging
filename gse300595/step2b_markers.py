"""Step 2b - quantify canonical marker panels per Leiden cluster."""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import PROC

PANELS = {
    "Malignant_epi": ["EPCAM", "KRT19", "KRT18", "KRT8", "MUC1", "CEACAM6",
                      "S100P", "TFF1", "TSPAN8", "MSLN", "CLDN4", "LGALS4"],
    "Ductal_normal": ["CFTR", "SLC4A4", "PKHD1", "ANXA4", "KRT7", "BICC1"],
    "Acinar":        ["PRSS1", "CPA1", "CELA3A", "PNLIP", "CTRB1", "REG1A", "REG1B", "GP2"],
    "Endocrine":     ["INS", "GCG", "SST", "PPY", "CHGA", "CHGB", "SLC30A8"],
    "Tuft":          ["POU2F3", "TRPM5", "PLCG2", "AVIL", "IL25"],
    "Fibroblast":    ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1", "LAMA2"],
    "Myofibro_SMC":  ["ACTA2", "TAGLN", "MYH11", "DES", "CNN1"],
    "Pericyte":      ["PDGFRB", "RGS5", "CSPG4", "NOTCH3", "KCNJ8"],
    "Endothelial":   ["PECAM1", "VWF", "CDH5", "CLDN5", "ERG", "PTPRB"],
    "Schwann":       ["PLP1", "S100B", "SOX10", "MPZ", "NRXN1", "CDH19"],
    "T_NK":          ["CD3D", "CD3E", "CD2", "TRBC1", "TRBC2", "IL7R",
                      "NKG7", "GNLY", "KLRD1", "PRF1"],
    "B":             ["MS4A1", "CD79A", "CD79B", "BANK1", "BLK", "IGKC"],
    "Plasma":        ["MZB1", "JCHAIN", "SDC1", "XBP1", "IGHG1"],
    "Mast":          ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
    "Myeloid":       ["CD68", "CD14", "LYZ", "AIF1", "C1QA", "CSF1R",
                      "ITGAM", "TYROBP", "FCER1G", "MRC1", "CD163", "SPI1"],
    "pDC":           ["LILRA4", "CLEC4C", "IRF7", "GZMB", "TCF4"],
    "cDC":           ["CD1C", "CLEC9A", "LAMP3", "FLT3", "ID2", "BATF3"],
    "Proliferating": ["MKI67", "TOP2A", "CCNB1", "RRM2", "PCNA"],
    "Neural_like":   ["NRXN1", "CACNA1A", "KCNMB2", "RIMS2", "SOX5", "DPP6"],
}

adata = sc.read_h5ad(os.path.join(PROC, "rna_clustered.h5ad"))
raw = adata.raw.to_adata()          # all genes, log-normalised
clusters = adata.obs["leiden"].astype(str)

rows = {}
for name, genes in PANELS.items():
    present = [g for g in genes if g in raw.var_names]
    sub = raw[:, present]
    df = pd.DataFrame(sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X,
                      columns=present, index=raw.obs_names)
    rows[name] = df.groupby(clusters).mean().mean(axis=1)

score_df = pd.DataFrame(rows).T
score_df.index.name = "panel"
score_df.to_csv(os.path.join(PROC, "cluster_panel_scores.tsv"), sep="\t")

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 100)
print("Cluster mean log-normalised expression, averaged over panel genes")
print((score_df * 100).round(0).astype(int).to_string())

# also print a few individual diagnostic genes
diag = ["EPCAM", "KRT19", "CEACAM6", "S100P", "CFTR", "PRSS1", "INS",
        "COL1A1", "PECAM1", "PTPRC", "CD68", "CD3E", "MS4A1", "CPA3", "PLP1"]
diag = [g for g in diag if g in raw.var_names]
sub = raw[:, diag]
df = pd.DataFrame(sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X,
                  columns=diag, index=raw.obs_names)
print("\nSelected diagnostic genes (x100):")
print((df.groupby(clusters).mean() * 100).round(0).astype(int).to_string())
