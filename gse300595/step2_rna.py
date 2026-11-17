"""GSE300595 step 2 - RNA preprocessing, Harmony integration, clustering.

Reads processed/rna_qc.h5ad (raw counts, QC-passed nuclei from step1),
normalises, selects HVGs, runs PCA + Harmony (batch = sample), builds the
neighbour graph and Leiden clusters, and writes:

    processed/rna_clustered.h5ad
    processed/leiden_markers.tsv      top markers per Leiden cluster
    processed/cluster_composition.tsv cluster x sample / response counts

Annotation is done separately (step2b) after inspecting the markers, so that
cell-type labels are driven by the data rather than by hard-coded guesses.
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import BASE, PROC, FIG

sc.settings.verbosity = 1


def main():
    t0 = time.time()
    print(f"scanpy {sc.__version__}", flush=True)
    adata = sc.read_h5ad(os.path.join(PROC, "rna_qc.h5ad"))
    print(f"Loaded {adata.shape[0]} nuclei x {adata.shape[1]} genes", flush=True)

    adata.layers["counts"] = adata.X.copy()

    sc.pp.filter_genes(adata, min_cells=10)
    print(f"After min_cells=10 gene filter: {adata.shape[1]} genes", flush=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample",
                                flavor="seurat")
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"HVG selected: {n_hvg}", flush=True)

    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    print(f"PCA done ({time.time()-t0:.0f}s)", flush=True)

    # ---- Harmony batch correction over the 12 samples ----------------------
    try:
        sc.external.pp.harmony_integrate(adata, "sample", max_iter_harmony=20)
        rep = "X_pca_harmony"
        print(f"Harmony integration done ({time.time()-t0:.0f}s)", flush=True)
    except Exception as e:                                    # pragma: no cover
        print(f"[WARN] Harmony failed ({e}); falling back to uncorrected PCA", flush=True)
        adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"]
        rep = "X_pca_harmony"

    sc.pp.neighbors(adata, n_neighbors=20, use_rep=rep, random_state=0)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden", random_state=0,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=0)

    n_cl = adata.obs["leiden"].nunique()
    print(f"Leiden clusters: {n_cl}", flush=True)
    print(adata.obs["leiden"].value_counts().sort_index().to_string(), flush=True)

    # ---- cluster markers ----------------------------------------------------
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(os.path.join(PROC, "leiden_markers.tsv"), sep="\t", index=False)

    # compact top-20 table for quick inspection
    top = (markers.sort_values(["group", "scores"], ascending=[True, False])
                  .groupby("group").head(20))
    top.to_csv(os.path.join(PROC, "leiden_top20.tsv"), sep="\t", index=False)

    # ---- composition --------------------------------------------------------
    comp = (adata.obs.groupby(["leiden", "response"], observed=True)
                  .size().unstack(fill_value=0))
    comp["total"] = comp.sum(axis=1)
    comp.to_csv(os.path.join(PROC, "cluster_composition.tsv"), sep="\t")
    print("\nCluster x response:\n" + comp.to_string(), flush=True)

    by_sample = (adata.obs.groupby(["leiden", "sample"], observed=True)
                       .size().unstack(fill_value=0))
    by_sample.to_csv(os.path.join(PROC, "cluster_by_sample.tsv"), sep="\t")

    adata.write_h5ad(os.path.join(PROC, "rna_clustered.h5ad"), compression="gzip")
    print(f"\n[saved] rna_clustered.h5ad  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
