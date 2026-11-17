"""
Step 4: Myeloid subset re-clustering (repo: adata_myeloid, leiden resolution=2)
to resolve the SPP1+MARCO+ macrophage niche, plus quantitative Serpin-axis readout.

Also fixes the Serpine positivity statistic (must use count-based detection, not
value > 0 on a continuous normalised matrix).
"""

import os

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "processed")
FIG_DIR = os.path.join(BASE, "figures")


def main():
    adata = sc.read_h5ad(os.path.join(OUT_DIR, "GSE289138_annotated.h5ad"))
    print(adata)

    # ------------------------------------------------ raw count detection
    counts = adata.layers["counts"]
    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    for gene in ["Serpine1", "Serpinb2", "Spp1", "Marco"]:
        if gene in var_idx:
            col = np.asarray(counts[:, var_idx[gene]].todense()).ravel()
            adata.obs[f"{gene}_detected"] = col > 0
            adata.obs[f"{gene}_count"] = col

    print("\n===== detection rates by annotation =====")
    rows = []
    for ct in adata.obs["anno"].unique():
        m = (adata.obs["anno"] == ct).to_numpy()
        rec = {"anno": ct, "n": int(m.sum())}
        for gene in ["Serpine1", "Serpinb2", "Spp1", "Marco"]:
            rec[f"{gene}_pct"] = round(float(adata.obs[f"{gene}_detected"][m].mean() * 100), 1)
        rows.append(rec)
    det_df = pd.DataFrame(rows).sort_values("n", ascending=False)
    det_df.to_csv(os.path.join(OUT_DIR, "serpin_detection_by_celltype.tsv"), sep="\t", index=False)
    print(det_df.to_string(index=False))

    # ------------------------------------------------ malignant epithelial Serpin
    mal = (adata.obs["malignant"] == "malignant").to_numpy()
    epi = (adata.obs["anno"] == "Epithelial").to_numpy()
    print("\n===== Serpin detection in malignant vs non-malignant epithelial =====")
    for label, mask in [("malignant", epi & mal), ("non-malignant epi", epi & ~mal)]:
        n = int(mask.sum())
        s1 = float(adata.obs["Serpine1_detected"][mask].mean() * 100)
        s2 = float(adata.obs["Serpinb2_detected"][mask].mean() * 100)
        print(f"  {label}: n={n}  Serpine1+={s1:.1f}%  Serpinb2+={s2:.1f}%")

    # ------------------------------------------------ myeloid subset
    my_mask = adata.obs["anno"].isin(["Macrophage", "Monocyte", "DC"]).to_numpy()
    my = adata[my_mask].copy()
    print(f"\n===== myeloid subset: {my.n_obs} cells =====")
    print(my.obs["anno"].value_counts().to_string())

    my.X = my.layers["scVI_normalized"]
    sc.pp.highly_variable_genes(
        my, flavor="seurat_v3", n_top_genes=3000, layer="counts", subset=True, batch_key="filename"
    )
    sc.pp.normalize_total(my, target_sum=1e4)
    sc.pp.log1p(my)
    my.raw = my

    sc.pp.neighbors(my, use_rep="X_scVI")
    sc.tl.umap(my)
    sc.tl.leiden(my, resolution=2, key_added="leiden_myeloid")
    sc.tl.rank_genes_groups(my, "leiden_myeloid", method="wilcoxon")
    print("\n--- myeloid subclusters ---")
    print(my.obs["leiden_myeloid"].value_counts().sort_index().to_string())

    # per-subcluster marker + serpin readout
    norm = my.layers["scVI_normalized"]
    mvar = {g: i for i, g in enumerate(my.var_names)}
    rows = []
    for cl in sorted(my.obs["leiden_myeloid"].unique()):
        m = (my.obs["leiden_myeloid"] == cl).to_numpy()
        rec = {"leiden_myeloid": cl, "n": int(m.sum())}
        for gene in ["Spp1", "Marco", "Serpine1", "Serpinb2", "C1qa", "C1qb", "C1qc",
                     "Lyz2", "Cd68", "Mrc1", "Arg1", "Ccr2", "Ly6c2", "Cd83", "Cd86", "H2-Oa"]:
            if gene in mvar:
                rec[gene] = round(float(np.asarray(norm[m, mvar[gene]]).mean()), 2)
        rows.append(rec)
    my_df = pd.DataFrame(rows)
    my_df.to_csv(os.path.join(OUT_DIR, "myeloid_subcluster_serpin.tsv"), sep="\t", index=False)
    print("\n--- myeloid subcluster marker expression ---")
    print(my_df.to_string(index=False))

    # top markers per myeloid subcluster for manual naming
    res = my.uns["rank_genes_groups"]
    print("\n===== TOP 12 MARKERS PER MYELOID SUBCLUSTER =====")
    for g in res["names"].dtype.names:
        print(f"\ncluster {g} (n={int((my.obs['leiden_myeloid']==g).sum())}): "
              + ", ".join(list(res["names"][g])[:12]))

    # SPP1+MARCO+ double positive macrophages (count-based detection, not normalised values)
    mcounts = my.layers["counts"]
    if "Spp1" in mvar and "Marco" in mvar:
        sp = np.asarray(mcounts[:, mvar["Spp1"]].todense()).ravel()
        ma = np.asarray(mcounts[:, mvar["Marco"]].todense()).ravel()
        # require robust expression: >=2 UMI for both genes
        both = (sp >= 2) & (ma >= 2)
        print(f"\nSPP1+MARCO+ myeloid cells (>=2 UMI both): {int(both.sum())} / {my.n_obs}")
        my.obs["SPP1_MARCO_pp"] = both
        sp_df = (
            my.obs.groupby("leiden_myeloid", observed=True)["SPP1_MARCO_pp"]
            .agg(["sum", "size"])
            .rename(columns={"sum": "spp1_marco_pp", "size": "n"})
        )
        sp_df["frac"] = (sp_df["spp1_marco_pp"] / sp_df["n"] * 100).round(1)
        print(sp_df.to_string())
        sp_df.to_csv(os.path.join(OUT_DIR, "myeloid_spp1_marco_pp.tsv"), sep="\t")

    sc.pl.umap(my, color=["anno", "leiden_myeloid", "Spp1", "Marco", "Serpine1"], save="_umap_myeloid.png", ncols=3)

    my.write(os.path.join(OUT_DIR, "myeloid.h5ad"))
    adata.write(os.path.join(OUT_DIR, "GSE289138_annotated.h5ad"))
    print("\nWROTE myeloid.h5ad + updated annotated h5ad")


if __name__ == "__main__":
    main()
