"""Step 2c - diagnostic expression of lineage + SERPIN/myeloid program genes."""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import PROC

adata = sc.read_h5ad(os.path.join(PROC, "rna_clustered.h5ad"))
raw = adata.raw.to_adata()
clusters = adata.obs["leiden"].astype(str)

GENES = [
    # islet / neural
    "GCG", "SST", "PPY", "CHGA", "CHGB", "NEUROD1", "NKX2-2", "ISL1",
    "SOX10", "PLP1", "S100B", "MPZ", "NRXN1", "KCNMB2", "CACNA1A",
    # epithelial / malignant
    "KRT19", "CEACAM6", "S100P", "MUC1", "TFF1", "MSLN", "CFTR", "MUC6", "MUC5B",
    # acinar
    "PRSS1", "CPA1", "PNLIP", "CTRB1", "REG1A",
    # SERPIN axis + myeloid program
    "SERPINE1", "SPP1", "MARCO", "APOE", "APOC1", "GPNMB", "TREM2",
    "CD163", "MRC1", "C1QA", "LPL", "PPARG", "CLIP4", "RGS1", "SGK1",
    "ITGAV", "PADI2", "ADAM19", "SLC11A1", "CD68", "AIF1", "C1QB",
    "CD14", "FCGR3A", "NOS2", "IL1B", "CXCL9", "CXCL10", "CD80",
    "FLT3", "AFF3", "CD1C", "CLEC9A", "LAMP3",
    # immune
    "PTPRC", "CD3E", "MS4A1", "MZB1", "CPA3", "TPSAB1",
]
GENES = [g for g in GENES if g in raw.var_names]
print(f"{len(GENES)}/{len(GENES)} genes present: {GENES}\n")

sub = raw[:, GENES]
X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
df = pd.DataFrame(X, columns=GENES, index=raw.obs_names)
means = df.groupby(clusters).mean() * 100
frac = df.groupby(clusters).apply(lambda d: (d > 0).mean()) * 100

pd.set_option("display.width", 400)
pd.set_option("display.max_columns", 200)
pd.set_option("display.max_rows", 200)

print("=== mean log-normalised expression x100 ===")
print(means.round(0).astype(int).to_string())
print("\n=== percent expressing ===")
print(frac.round(0).astype(int).to_string())

means.to_csv(os.path.join(PROC, "cluster_diagnostic_means.tsv"), sep="\t")
frac.to_csv(os.path.join(PROC, "cluster_diagnostic_frac.tsv"), sep="\t")
