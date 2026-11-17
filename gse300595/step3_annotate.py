"""GSE300595 step 3 - cell-type annotation + malignancy (CNV) calling.

1. Assign a coarse cell type to each Leiden cluster using the canonical marker
   panels quantified in step2b/2c.
2. Run a compact inferCNV-style copy-number score on the epithelial clusters to
   separate malignant epithelium from normal duct / acinar / islet cells.
   Method: cluster-mean expression, genes sorted by genomic position, 101-gene
   moving average, centred on the immune-cell reference, burden = mean |dev|.
3. Reproduce the GSE202051 SERPIN definition exactly:
       SERPIN_score = sc.tl.score_genes(["SERPINE1", "SERPINB2"])
       SERPIN_high  = score > median within malignant epithelium
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.ndimage import uniform_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import BASE, PROC

GTF = os.path.join(BASE, "ref", "gencode.v44.basic.annotation.gtf.gz")

# coarse cluster labels from step2b/2c marker quantification
CLUSTER_LABEL = {
    "0": "Ductal", "1": "Ductal", "2": "Ductal",
    "3": "Endocrine_alpha_pp", "4": "Endocrine_alpha",
    "5": "Endocrine_beta", "12": "Endocrine_delta",
    "6": "Tuft",
    "7": "Myeloid",
    "8": "Epithelial_malignant", "13": "Epithelial_malignant_prolif",
    "9": "Ductal_mucinous", "25": "Ductal_mucinous",
    "10": "Endothelial", "11": "Fibroblast", "18": "Schwann",
    "19": "Pericyte_SMC",
    "14": "Plasma", "15": "T_NK", "17": "B", "20": "Mast",
    "16": "Acinar_ductal", "24": "Acinar_ductal", "26": "Acinar_ductal",
    "21": "Acinar", "22": "Acinar", "23": "Acinar",
}
BROAD = {
    "Ductal": "Epithelial", "Ductal_mucinous": "Epithelial",
    "Epithelial_malignant": "Epithelial", "Epithelial_malignant_prolif": "Epithelial",
    "Acinar": "Acinar", "Acinar_ductal": "Acinar",
    "Endocrine_alpha_pp": "Endocrine", "Endocrine_alpha": "Endocrine",
    "Endocrine_beta": "Endocrine", "Endocrine_delta": "Endocrine",
    "Tuft": "Epithelial",
    "Myeloid": "Immune", "Plasma": "Immune", "T_NK": "Immune",
    "B": "Immune", "Mast": "Immune",
    "Endothelial": "Stromal", "Fibroblast": "Stromal",
    "Schwann": "Stromal", "Pericyte_SMC": "Stromal",
}


def load_gene_positions():
    pos = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[2] != "gene":
                continue
            attrs = parts[8]
            name = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    name = kv.split(" ", 1)[1].strip('"')
                    break
            if name and name not in pos:
                pos[name] = (parts[0], int(parts[3]))
    return pos


def main():
    adata = sc.read_h5ad(os.path.join(PROC, "rna_clustered.h5ad"))
    print(f"Loaded {adata.shape}", flush=True)

    adata.obs["leiden"] = adata.obs["leiden"].astype(str)
    adata.obs["celltype"] = adata.obs["leiden"].map(CLUSTER_LABEL)
    adata.obs["broadtype"] = adata.obs["celltype"].map(BROAD)
    missing = adata.obs["celltype"].isna()
    if missing.any():
        print(f"[WARN] {missing.sum()} cells in unlabelled clusters: "
              f"{sorted(adata.obs.loc[missing,'leiden'].unique())}", flush=True)

    # ------------------------------------------------------------------ CNV
    print("\n=== CNV inference (inferCNV-style) ===", flush=True)
    raw = adata.raw.to_adata()
    pos = load_gene_positions()
    genes = []
    for g in raw.var_names:
        key = g.split("-")[0] if g not in pos else g
        if key in pos:
            genes.append((g, pos[key][0], pos[key][1]))
    gdf = pd.DataFrame(genes, columns=["gene", "chrom", "start"])
    gdf = gdf[gdf["chrom"].isin([f"chr{i}" for i in range(1, 23)] + ["chrX"])]
    gdf["chrom"] = pd.Categorical(gdf["chrom"],
                                  categories=[f"chr{i}" for i in range(1, 23)] + ["chrX"],
                                  ordered=True)
    gdf = gdf.sort_values(["chrom", "start"]).reset_index(drop=True)
    print(f"Genes with coordinates: {len(gdf)}", flush=True)

    X = raw[:, gdf["gene"].tolist()].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    clusters = adata.obs["leiden"].astype(str).values
    uniq = sorted(pd.unique(clusters), key=int)
    M = np.vstack([X[clusters == c].mean(axis=0) for c in uniq])
    Mdf = pd.DataFrame(M, index=uniq, columns=gdf["gene"].tolist())

    # smooth along the genome within each chromosome
    smooth = np.empty_like(M)
    for chrom, idx in gdf.groupby("chrom", observed=True).indices.items():
        sub = Mdf.iloc[:, idx].values
        smooth[:, idx] = uniform_filter1d(sub, size=101, axis=1, mode="nearest")

    immune_cl = [c for c in uniq if c in ("7", "14", "15", "17", "20")]
    ref = smooth[[uniq.index(c) for c in immune_cl]].mean(axis=0)
    cnv = smooth - ref

    burden = np.abs(cnv).mean(axis=1)
    # chromosome-level mean deviation, for interpretation
    chr_dev = {}
    for chrom, idx in gdf.groupby("chrom", observed=True).indices.items():
        chr_dev[chrom] = cnv[:, idx].mean(axis=1)
    chr_df = pd.DataFrame(chr_dev, index=uniq)

    cnv_df = pd.DataFrame({
        "leiden": uniq,
        "celltype": [CLUSTER_LABEL.get(c, "?") for c in uniq],
        "n_cells": [int((clusters == c).sum()) for c in uniq],
        "cnv_burden": burden,
    })
    cnv_df.to_csv(os.path.join(PROC, "cnv_burden_by_cluster.tsv"), sep="\t", index=False)
    chr_df.to_csv(os.path.join(PROC, "cnv_chrom_deviation.tsv"), sep="\t")

    pd.set_option("display.width", 200)
    print(cnv_df.sort_values("cnv_burden", ascending=False).to_string(index=False),
          flush=True)

    # ---- malignancy calling -------------------------------------------------
    # PRIMARY (mirrors GSE202051 exactly): among epithelial cells, a malignancy
    # score (malignant-epithelial marker programme) above the epithelial median.
    epi_types = ["Ductal", "Ductal_mucinous", "Epithelial_malignant",
                 "Epithelial_malignant_prolif", "Acinar_ductal"]
    adata.obs["is_epi"] = adata.obs["celltype"].isin(epi_types).values
    epi_mask = adata.obs["is_epi"].values

    mal_genes = [g for g in ["KRT19", "CEACAM6", "S100P", "MSLN", "TFF1",
                             "MUC1", "EPCAM", "KRT8", "KRT18"]
                 if g in adata.raw.var_names]
    sc.tl.score_genes(adata, mal_genes, score_name="malignancy_score",
                      use_raw=True, random_state=0)
    mal_thr = float(np.median(adata.obs.loc[epi_mask, "malignancy_score"]))
    adata.obs["is_malignant"] = False
    adata.obs.loc[epi_mask, "is_malignant"] = (
        adata.obs.loc[epi_mask, "malignancy_score"] > mal_thr)
    print(f"\nMalignancy score genes: {mal_genes}", flush=True)
    print(f"Epithelial cells: {int(epi_mask.sum())}, threshold={mal_thr:.4f}, "
          f"malignant={int(adata.obs['is_malignant'].sum())}", flush=True)

    # SECONDARY (sensitivity): cluster-level CNV burden above non-epithelial
    # background (mean + 1 SD).  Reported alongside, not used for primary calls.
    bg = cnv_df[~cnv_df["celltype"].isin(epi_types + ["Tuft"])]["cnv_burden"]
    bg_mean, bg_sd = float(bg.mean()), float(bg.std())
    cnv_thr = bg_mean + bg_sd
    epi = cnv_df[cnv_df["celltype"].isin(epi_types)].copy()
    epi["cnv_call"] = np.where(epi["cnv_burden"] > cnv_thr, "malignant", "normal")
    print(f"\nNon-epithelial CNV background: mean={bg_mean:.4f} sd={bg_sd:.4f} "
          f"-> threshold={cnv_thr:.4f}", flush=True)
    print(epi.sort_values("cnv_burden", ascending=False).to_string(index=False),
          flush=True)
    cnv_mal = set(epi.loc[epi["cnv_call"] == "malignant", "leiden"])
    adata.obs["is_malignant_cnv"] = adata.obs["leiden"].isin(cnv_mal)
    print(f"CNV-malignant clusters (sensitivity): {sorted(cnv_mal, key=int)}", flush=True)

    mal_mask = adata.obs["is_malignant"].values

    # ------------------------------------------------- SERPIN (GSE202051 def.)
    serpin_genes = [g for g in ["SERPINE1", "SERPINB2"] if g in adata.raw.var_names]
    print(f"\nSERPIN genes available in full matrix: {serpin_genes}", flush=True)
    # SERPINB2 is essentially undetectable in these nuclear profiles (~0.03 % of
    # nuclei); score_genes therefore effectively reduces to SERPINE1.  Kept in
    # the gene list for exact parity with the GSE202051 definition.
    sc.tl.score_genes(adata, serpin_genes, score_name="SERPIN_score",
                      use_raw=True, random_state=0)

    adata.obs["SERPIN_status"] = "non_malignant"
    thr_s = float(np.median(adata.obs.loc[mal_mask, "SERPIN_score"]))
    adata.obs.loc[mal_mask, "SERPIN_status"] = np.where(
        adata.obs.loc[mal_mask, "SERPIN_score"] > thr_s, "SERPIN_high", "SERPIN_low")
    print(f"SERPIN threshold (median of malignant): {thr_s:.4f}", flush=True)
    print(adata.obs["SERPIN_status"].value_counts().to_string(), flush=True)

    # ------------------------------- per-patient malignant/SERPIN fractions
    rows = []
    for (s, resp), grp in adata.obs.groupby(["sample", "response"], observed=True):
        mal = grp[grp["is_malignant"].astype(bool)]
        mal_cnv = grp[grp["is_malignant_cnv"].astype(bool)]
        rows.append(dict(
            sample=s, response=resp, n_cells=len(grp),
            n_malignant=len(mal),
            n_SERPIN_high=int((mal["SERPIN_status"] == "SERPIN_high").sum()),
            SERPIN_frac=(mal["SERPIN_status"] == "SERPIN_high").mean() if len(mal) else np.nan,
            mean_SERPIN_score=mal["SERPIN_score"].mean() if len(mal) else np.nan,
            n_malignant_cnv=len(mal_cnv),
            SERPIN_frac_cnv=(mal_cnv["SERPIN_score"] >
                             np.median(adata.obs.loc[adata.obs["is_malignant_cnv"],
                                                     "SERPIN_score"])).mean()
                            if len(mal_cnv) else np.nan,
            n_myeloid=int((grp["celltype"] == "Myeloid").sum()),
            myeloid_frac=(grp["celltype"] == "Myeloid").mean(),
        ))
    pat = pd.DataFrame(rows).sort_values(["response", "sample"])
    pat.to_csv(os.path.join(PROC, "patient_serpin_myeloid.tsv"), sep="\t", index=False)
    print("\nPer-sample malignant / SERPIN / myeloid:", flush=True)
    print(pat.to_string(index=False), flush=True)

    adata.obs["celltype"] = adata.obs["celltype"].astype("category")
    adata.obs["broadtype"] = adata.obs["broadtype"].astype("category")
    adata.write_h5ad(os.path.join(PROC, "rna_annotated.h5ad"), compression="gzip")
    print("\n[saved] rna_annotated.h5ad", flush=True)


if __name__ == "__main__":
    main()
