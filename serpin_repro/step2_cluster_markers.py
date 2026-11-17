"""
Step 2: Cluster marker export + manual annotation basis (mirrors repo's
2_human_scRNAseq_atlas_2/3 notebooks: annotate leiden clusters by markers).

Writes processed/cluster_top_markers.tsv and prints per-cluster top genes so the
annotation map can be defined explicitly (repo style: dict leiden -> cell type).
"""

import os

import numpy as np
import pandas as pd
import scanpy as sc

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "processed")
FIG_DIR = os.path.join(BASE, "figures")

adata = sc.read_h5ad(os.path.join(OUT_DIR, "GSE289138_processed.h5ad"))
print(adata)
print(f"cells={adata.n_obs} genes={adata.n_vars}")

res = adata.uns["rank_genes_groups"]
groups = list(res["names"].dtype.names)

rows = []
for g in groups:
    names = list(res["names"][g])
    scores = list(res["scores"][g])
    lfc = list(res["logfoldchanges"][g])
    padj = list(res["pvals_adj"][g])
    for i in range(len(names)):
        rows.append(
            {
                "leiden": g,
                "rank": i + 1,
                "gene": names[i],
                "score": scores[i],
                "log2FC": lfc[i],
                "padj": padj[i],
            }
        )

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "cluster_top_markers.tsv"), sep="\t", index=False)

# top 15 per cluster for quick manual inspection
print("\n===== TOP 15 MARKERS PER LEIDEN CLUSTER =====")
for g in groups:
    sub = df[(df.leiden == g) & (df["rank"] <= 15)].sort_values("rank")
    n = int((adata.obs["leiden"] == g).sum())
    print(f"\ncluster {g} (n={n}): " + ", ".join(sub["gene"].tolist()))

# also dump canonical marker expression per cluster for annotation sanity
CANON = {
    "Epithelial": ["Epcam", "Krt8", "Krt18", "Krt19", "Krt5", "Krt6a", "S100p"],
    "T/NK": ["Cd3e", "Cd3d", "Trbc2", "Cd8a", "Nkg7", "Gzma", "Gzmb", "Prf1"],
    "Macrophage": ["Adgre1", "C1qa", "C1qb", "C1qc", "Cd68", "Lyz2", "Marco", "Spp1", "Mrc1", "Ctss"],
    "Monocyte": ["Ly6c2", "Ccr2", "Plac8", "Chil3", "Ms4a4c", "Ifitm6"],
    "B cell": ["Cd79a", "Cd79b", "Ms4a1", "Ebf1", "Ighm"],
    "Plasma cell": ["Mzb1", "Jchain", "Sdc1", "Xbp1", "Derl3"],
    "CAF/Stromal": ["Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Dcn", "Lum", "Acta2", "Postn"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Esam", "Egfl7", "Cldn5"],
    "Mast cell": ["Kit", "Cpa3", "Mcpt4", "Fcer1a", "Ms4a2"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg", "Csf3r", "Mmp8"],
    "pDC": ["Siglech", "Bst2", "Ccr9", "Irf8"],
    "DC": ["Flt3", "Cd74", "H2-Ab1", "Cd83", "Fscn1"],
    "RBC": ["Hba-a1", "Hba-a2", "Hbb-bs", "Alas2"],
}
norm = adata.layers["scVI_normalized"]
var_idx = {g: i for i, g in enumerate(adata.var_names)}
print("\n===== CANONICAL MARKER MEAN EXPRESSION PER CLUSTER =====")
recs = []
for g in groups:
    mask = (adata.obs["leiden"] == g).to_numpy()
    rec = {"leiden": g, "n": int(mask.sum())}
    for ct, genes in CANON.items():
        present = [x for x in genes if x in var_idx]
        if present:
            vals = np.asarray(norm[mask][:, [var_idx[x] for x in present]]).mean()
            rec[ct] = round(float(vals), 2)
    recs.append(rec)
canon_df = pd.DataFrame(recs).sort_values("n", ascending=False)
canon_df.to_csv(os.path.join(OUT_DIR, "cluster_canonical_marker_scores.tsv"), sep="\t", index=False)
print(canon_df.to_string(index=False))

print("\n===== SAMPLE COMPOSITION PER CLUSTER =====")
comp = pd.crosstab(adata.obs["leiden"], adata.obs["filename"])
print(comp.to_string())
comp.to_csv(os.path.join(OUT_DIR, "cluster_sample_composition.tsv"), sep="\t")
