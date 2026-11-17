"""
Step 3: Manual annotation (repo-style leiden -> cell type map) + epithelial
inferCNV (mirrors 2_human_scRNAseq_atlas_2_infercnv.ipynb / _3_annotations).

Outputs:
  processed/GSE289138_annotated.h5ad
  processed/cnv_scores.tsv
  figures/umap_annotated.png, umap_cnv.png
"""

import os

import matplotlib

matplotlib.use("Agg")
import infercnvpy as cnv
import numpy as np
import pandas as pd
import scanpy as sc

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "processed")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Manual annotation derived from per-cluster marker analysis (step2 output)
ANNO = {
    "0": "Macrophage",
    "1": "Epithelial",
    "2": "T/NK",
    "3": "Epithelial",
    "4": "Neutrophil",
    "5": "Epithelial",
    "6": "Macrophage",
    "7": "CAF",
    "8": "Epithelial",
    "9": "Epithelial",
    "10": "Macrophage",
    "11": "Epithelial",
    "12": "T/NK",
    "13": "Endothelial",
    "14": "CAF",
    "15": "B cell",
    "16": "Epithelial",
    "17": "CAF",
    "18": "Monocyte",
    "19": "DC",
    "20": "Macrophage",
    "21": "Exocrine",
    "22": "Neutrophil",
    "23": "RBC",
}

def main():
    adata = sc.read_h5ad(os.path.join(OUT_DIR, "GSE289138_processed.h5ad"))
    adata.obs["anno"] = adata.obs["leiden"].astype(str).map(ANNO).fillna("Unknown")
    print(adata.obs["anno"].value_counts().to_string())

    sc.pl.umap(adata, color=["anno", "leiden", "condition"], save="_umap_annotated.png", ncols=3)

    # ---------------------------------------------------------------- inferCNV
    # reference = non-epithelial stromal/immune cells; target = epithelial
    epi_mask = (adata.obs["anno"] == "Epithelial").to_numpy()
    ref_mask = adata.obs["anno"].isin(["CAF", "Endothelial"]).to_numpy()
    print(f"epithelial={epi_mask.sum()} reference={ref_mask.sum()}")

    epi = adata[epi_mask].copy()
    # inferCNV works on log-normalized counts; use scVI-normalized layer for stability
    norm = adata.layers["scVI_normalized"]
    epi.X = norm[epi_mask].copy()

    ref_adata = adata[ref_mask].copy()
    ref_adata.X = norm[ref_mask].copy()

    # --- attach mouse genomic coordinates (Ensembl BioMart) ---
    pos_path = os.path.join(BASE, "data", "mouse_gene_pos.tsv")
    pos = pd.read_csv(pos_path, sep="\t", header=None, names=["ensembl", "symbol", "chrom", "start", "end"])
    pos = pos.dropna(subset=["chrom", "start", "end"])
    # normalise chromosome naming to Ensembl style used by infercnvpy (chr1..chr19, chrX)
    def norm_chrom(c):
        c = str(c)
        if c.startswith("chr"):
            return c
        if c == "MT":
            return "chrM"
        if c.isdigit():
            return f"chr{c}"
        return f"chr{c}"

    pos["chrom"] = pos["chrom"].map(norm_chrom)
    pos = pos.drop_duplicates(subset=["symbol"])
    for obj in (epi, ref_adata):
        obj.var["chromosome"] = obj.var_names.map(pos.set_index("symbol")["chrom"])
        obj.var["start"] = obj.var_names.map(pos.set_index("symbol")["start"])
        obj.var["end"] = obj.var_names.map(pos.set_index("symbol")["end"])
        for c in ["chromosome", "start", "end"]:
            obj.var[c] = obj.var[c].astype(object)
    epi = epi[:, epi.var["chromosome"].notna() & epi.var["start"].notna()].copy()
    ref_adata = ref_adata[:, ref_adata.var["chromosome"].notna() & ref_adata.var["start"].notna()].copy()
    epi.var["start"] = epi.var["start"].astype(int)
    epi.var["end"] = epi.var["end"].astype(int)
    ref_adata.var["start"] = ref_adata.var["start"].astype(int)
    ref_adata.var["end"] = ref_adata.var["end"].astype(int)
    print(f"epithelial with coords: {epi.n_vars} genes; reference: {ref_adata.n_vars} genes")

    cnv.tl.infercnv(
        epi,
        reference=ref_adata.X,
        window_size=250,
        step=50,
        exclude_chromosomes=["chrX", "chrY"],
        n_jobs=1,
        layer=None,
        key_added="cnv",
        inplace=True,
    )
    # infercnvpy stores the CNV matrix in obsm["X_cnv"]; derive a per-cell score
    pca_out = sc.pp.pca(epi.obsm["X_cnv"], n_comps=30, svd_solver="arpack", zero_center=False)
    epi.obsm["cnv_pca"] = pca_out if isinstance(pca_out, np.ndarray) else epi.obsm["X_pca"]
    sc.pp.neighbors(epi, use_rep="cnv_pca", n_neighbors=15, key_added="cnv_neighbors")
    cnv.tl.leiden(epi, resolution=1.0, key_added="cnv_leiden")
    cnv.tl.cnv_score(epi, groupby="cnv_leiden", use_rep="cnv", key_added="cnv_score")

    # malignant = high CNV score epithelial
    cnv_vals = epi.obs["cnv_score"].to_numpy()
    thr = float(np.percentile(cnv_vals, 60))
    epi.obs["malignant"] = np.where(cnv_vals > thr, "malignant", "non-malignant")
    print(f"CNV threshold(60th pct)={thr:.4f}")
    print(epi.obs["malignant"].value_counts().to_string())

    epi.obs["leiden_epi"] = epi.obs["leiden"].astype(str)
    sc.tl.leiden(epi, resolution=0.3, key_added="leiden_epi")
    sc.tl.rank_genes_groups(epi, "leiden_epi", method="wilcoxon")
    print("\n--- epithelial subclusters ---")
    print(epi.obs["leiden_epi"].value_counts().sort_index().to_string())

    sc.pl.umap(
        epi,
        color=["cnv_score", "malignant", "leiden_epi", "Serpine1", "Serpinb2", "Spp1"],
        save="_umap_cnv_epithelial.png",
        ncols=3,
    )

    # annotate back into main object
    adata.obs["cnv_score"] = np.nan
    adata.obs.loc[epi.obs_names, "cnv_score"] = epi.obs["cnv_score"]
    adata.obs["malignant"] = "NA"
    adata.obs.loc[epi.obs_names, "malignant"] = epi.obs["malignant"]
    adata.obs["leiden_epi"] = "NA"
    adata.obs.loc[epi.obs_names, "leiden_epi"] = epi.obs["leiden_epi"]

    # epithelial cluster composition
    epi_comp = pd.crosstab(epi.obs["leiden_epi"], epi.obs["malignant"])
    print("\n--- CNV by epithelial subcluster ---")
    print(epi_comp.to_string())

    # Serpin expression per epithelial subcluster (core paper readout)
    var_idx = {g: i for i, g in enumerate(epi.var_names)}
    rows = []
    for cl in sorted(epi.obs["leiden_epi"].unique()):
        m = (epi.obs["leiden_epi"] == cl).to_numpy()
        rec = {
            "leiden_epi": cl,
            "n": int(m.sum()),
            "malignant_frac": round(float((epi.obs["malignant"][m] == "malignant").mean()), 3),
        }
        for gene in ["Serpine1", "Serpinb2", "Spp1", "Marco", "Krt19", "Onecut2"]:
            if gene in var_idx:
                rec[gene] = round(float(np.asarray(epi.X[m, var_idx[gene]]).mean()), 3)
        rows.append(rec)
    serpin_epi = pd.DataFrame(rows)
    serpin_epi.to_csv(os.path.join(OUT_DIR, "epithelial_serpin_by_cluster.tsv"), sep="\t", index=False)
    print("\n--- Serpin expression per epithelial subcluster ---")
    print(serpin_epi.to_string(index=False))

    # per-cell Serpine1/B2 status among malignant epithelial cells (count-based)
    epi_counts = adata.layers["counts"][epi_mask]
    epi_counts = epi_counts[:, [adata.var_names.get_loc(g) for g in epi.var_names]]
    mal = (epi.obs["malignant"] == "malignant").to_numpy()
    for gene in ["Serpine1", "Serpinb2"]:
        if gene in var_idx:
            vals = np.asarray(epi_counts[:, var_idx[gene]].todense()).ravel()
        else:
            continue
        for label, mask in [("malignant", mal), ("non-malignant", ~mal)]:
            pos = (vals >= 2) & mask
            print(f"{gene}+ (>=2 UMI) {label} epithelial: {int(pos.sum())} / {int(mask.sum())} "
                  f"({pos.sum() / max(mask.sum(), 1) * 100:.1f}%)")

    adata.write(os.path.join(OUT_DIR, "GSE289138_annotated.h5ad"))
    epi.write(os.path.join(OUT_DIR, "epithelial_cnv.h5ad"))
    print("\nWROTE annotated h5ad + epithelial cnv h5ad")


if __name__ == "__main__":
    main()
