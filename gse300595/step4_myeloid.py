"""GSE300595 step 4 - myeloid subclustering, states, and response contrast.

Subsets the myeloid compartment (Leiden cluster 7), re-clusters it from raw
counts, annotates macrophage states (LAM / M1-like / M2-like / DC-like /
proliferating) from module scores, and contrasts responders (n=5) vs
non-responders (n=7).

Statistics are patient-level (12 patients): cells within a patient are not
independent.  Per-patient state fractions and module means are compared with
Mann-Whitney U; per-gene testing uses patient pseudobulk.  With n=5 vs n=7 the
result is descriptive / hypothesis-generating, not confirmatory.
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import BASE, PROC

MODULES = {
    "LAM":     ["APOE", "APOC1", "GPNMB", "TREM2", "SPP1", "LPL", "PPARG",
                "CLIP4", "CTSD", "FABP5"],
    "M1_like": ["NOS2", "IL1B", "CXCL9", "CXCL10", "CXCL11", "CD80", "SOCS3",
                "IL6", "TNF"],
    "M2_like": ["CD163", "MRC1", "MSR1", "RGS1", "SGK1", "ITGAV", "ADAM19",
                "PADI2", "MAF", "F13A1"],
    "DC_like": ["FLT3", "AFF3", "CD1C", "CLEC9A", "LAMP3", "ID2", "BATF3",
                "IRF8"],
    "Prolif":  ["MKI67", "TOP2A", "CCNB1", "RRM2", "PCNA"],
}
STATE_MODULES = ["LAM", "M1_like", "M2_like", "DC_like", "Prolif"]


def mwu(a, b, label):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        print(f"  {label:28s}: insufficient n", flush=True)
        return np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Cliff's delta
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    delta = (gt - lt) / (len(a) * len(b))
    print(f"  {label:28s}: R={np.median(a):+.3f} (n={len(a)})  "
          f"NR={np.median(b):+.3f} (n={len(b)})  p={p:.4f}  delta={delta:+.3f}",
          flush=True)
    return p


def main():
    ann = sc.read_h5ad(os.path.join(PROC, "rna_annotated.h5ad"))
    my_ids = ann.obs_names[ann.obs["celltype"] == "Myeloid"]
    print(f"Myeloid nuclei in cluster 7: {len(my_ids)}", flush=True)

    rna = sc.read_h5ad(os.path.join(PROC, "rna_qc.h5ad"))
    adata = rna[my_ids].copy()
    del rna
    adata.obs["response"] = ann.obs.loc[adata.obs_names, "response"].values
    adata.obs["sample"] = ann.obs.loc[adata.obs_names, "sample"].values
    adata.obs["patient"] = ann.obs.loc[adata.obs_names, "patient"].values

    adata.layers["counts"] = adata.X.copy()
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # keep full-gene raw counts for pseudobulk (before HVG subsetting)
    counts_all = adata.layers["counts"].copy()
    genes_all = adata.var_names.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=1500, batch_key="sample")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=40, svd_solver="arpack")
    sc.external.pp.harmony_integrate(adata, "sample", max_iter_harmony=30)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony", random_state=0)
    sc.tl.leiden(adata, resolution=0.6, key_added="myeloid_cluster",
                 random_state=0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=0)
    print(f"Myeloid subclusters: {adata.obs['myeloid_cluster'].nunique()}", flush=True)

    for name, genes in MODULES.items():
        present = [g for g in genes if g in adata.raw.var_names]
        sc.tl.score_genes(adata, present, score_name=f"score_{name}",
                          use_raw=True, random_state=0)

    # ---- annotate states from module scores ------------------------------
    # Module scores are not comparable in absolute terms (M2 markers such as
    # CD163/MRC1 are broadly expressed by all macrophages here), so each module
    # is z-scored ACROSS subclusters and the cluster adopts its most
    # distinguishing module.
    sc_cols = [f"score_{m}" for m in STATE_MODULES]
    cl_scores = adata.obs.groupby("myeloid_cluster", observed=True)[sc_cols].mean()
    print("\n=== module scores by subcluster ===", flush=True)
    print(cl_scores.round(3).to_string(), flush=True)

    z = (cl_scores - cl_scores.mean()) / cl_scores.std(ddof=0)
    print("\n=== module scores z-scored across subclusters ===", flush=True)
    print(z.round(2).to_string(), flush=True)

    state = {}
    for cl, row in z.iterrows():
        best = row.idxmax().replace("score_", "")
        state[cl] = best if row.max() > 0.5 else "Macrophage_other"
    adata.obs["myeloid_state"] = adata.obs["myeloid_cluster"].map(state)

    print("\n=== cluster -> state mapping ===", flush=True)
    n_map = adata.obs.groupby(["myeloid_cluster", "myeloid_state"],
                              observed=True).size()
    print(n_map.to_string(), flush=True)

    # ---- subtype markers --------------------------------------------------
    marker_genes = ["APOE", "APOC1", "GPNMB", "TREM2", "SPP1", "LPL", "PPARG",
                    "CLIP4", "CD163", "MRC1", "RGS1", "SGK1", "ITGAV", "PADI2",
                    "ADAM19", "NOS2", "IL1B", "CXCL9", "CXCL10", "CD80",
                    "FLT3", "AFF3", "CD1C", "CLEC9A", "LAMP3", "HLA-DRA",
                    "CD74", "C1QA", "MKI67", "SLC11A1", "MARCO", "LYVE1",
                    "CX3CR1", "FCGR3A", "CD14"]
    marker_genes = [g for g in marker_genes if g in adata.raw.var_names]
    raw = adata.raw.to_adata()
    sub = raw[:, marker_genes]
    Xm = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    mdf = pd.DataFrame(Xm, columns=marker_genes, index=raw.obs_names)
    mdf["state"] = adata.obs["myeloid_state"].values
    means = mdf.groupby("state", observed=True)[marker_genes].mean() * 100
    pd.set_option("display.width", 400)
    pd.set_option("display.max_columns", 200)
    print("\n=== mean expression x100 by myeloid STATE ===", flush=True)
    print(means.round(0).astype(int).to_string(), flush=True)
    means.to_csv(os.path.join(PROC, "myeloid_state_marker_means.tsv"), sep="\t")

    # ---- per-patient state composition -----------------------------------
    pat_rows = []
    for s, grp in adata.obs.groupby("sample", observed=True):
        d = dict(sample=s, response=grp["response"].iloc[0], n_myeloid=len(grp))
        vc = grp["myeloid_state"].value_counts()
        for st in sorted(adata.obs["myeloid_state"].unique()):
            d[f"frac_{st}"] = vc.get(st, 0) / len(grp)
        for c in sc_cols:
            d[f"{c}_mean"] = grp[c].mean()
        pat_rows.append(d)
    pat = pd.DataFrame(pat_rows).sort_values(["response", "sample"])
    pat.to_csv(os.path.join(PROC, "myeloid_patient_states.tsv"), sep="\t", index=False)
    print("\n=== per-patient myeloid state fractions ===", flush=True)
    print(pat.to_string(index=False), flush=True)

    # ---- responder vs non-responder tests (patient level) ----------------
    print("\n=== RESPONDER vs NON-RESPONDER (patient-level Mann-Whitney) ===",
          flush=True)
    R = pat[pat["response"] == "responder"]
    NR = pat[pat["response"] == "nonresponder"]
    for st in sorted(adata.obs["myeloid_state"].unique()):
        mwu(R[f"frac_{st}"], NR[f"frac_{st}"], f"fraction {st}")
    print("", flush=True)
    for c in sc_cols:
        mwu(R[f"{c}_mean"], NR[f"{c}_mean"], f"module {c}")

    # ---- pseudobulk DE ----------------------------------------------------
    print("\n=== pseudobulk differential expression (R vs NR, myeloid) ===",
          flush=True)
    counts = pd.DataFrame(counts_all.toarray() if hasattr(counts_all, "toarray")
                          else np.asarray(counts_all),
                          index=adata.obs_names, columns=genes_all)
    samples = adata.obs["sample"].values
    uniq = sorted(pd.unique(samples))
    pb = np.vstack([counts.loc[samples == s].sum(axis=0).values for s in uniq])
    pb_df = pd.DataFrame(pb, index=uniq, columns=genes_all)
    pb_df.to_csv(os.path.join(PROC, "myeloid_pseudobulk_counts.tsv"), sep="\t")

    cpm = pb_df.div(pb_df.sum(axis=1), axis=0) * 1e6
    logcpm = np.log1p(cpm)
    samp_resp = pd.Series({s: adata.obs.loc[adata.obs["sample"] == s,
                                             "response"].iloc[0] for s in uniq})
    r_s = list(samp_resp[samp_resp == "responder"].index)
    nr_s = list(samp_resp[samp_resp == "nonresponder"].index)
    print(f"  responders n={len(r_s)}, non-responders n={len(nr_s)}", flush=True)

    det = (cpm > 1).sum(axis=0) >= 6          # detected in >=6/12 patients
    expr = logcpm.loc[r_s + nr_s, det]
    res = []
    for g in expr.columns:
        v1 = expr.loc[r_s, g].values
        v0 = expr.loc[nr_s, g].values
        try:
            u, p = stats.mannwhitneyu(v1, v0, alternative="two-sided")
        except ValueError:
            continue
        m1 = cpm.loc[r_s, g].mean() + 1.0     # pseudocount of 1 CPM
        m0 = cpm.loc[nr_s, g].mean() + 1.0
        res.append(dict(gene=g, cpm_R=m1 - 1, cpm_NR=m0 - 1,
                        log2FC=np.log2(m1 / m0), p=p))
    de = pd.DataFrame(res)
    de["padj"] = np.nan
    try:
        from statsmodels.stats.multitest import multipletests
        de["padj"] = multipletests(de["p"], method="fdr_bh")[1]
    except Exception:
        pass
    de = de.sort_values("p")
    de.to_csv(os.path.join(PROC, "myeloid_pseudobulk_DE.tsv"), sep="\t", index=False)

    cols = ["gene", "cpm_R", "cpm_NR", "log2FC", "p", "padj"]
    print("\n  Top 20 higher in RESPONDERS:", flush=True)
    print(de.sort_values("log2FC", ascending=False).head(20)[cols]
          .to_string(index=False), flush=True)
    print("\n  Top 20 higher in NON-RESPONDERS:", flush=True)
    print(de.sort_values("log2FC", ascending=True).head(20)[cols]
          .to_string(index=False), flush=True)
    print(f"\n  Genes with padj<0.1: {int((de['padj'] < 0.1).sum())} "
          f"| padj<0.25: {int((de['padj'] < 0.25).sum())}", flush=True)

    adata.write_h5ad(os.path.join(PROC, "myeloid_annotated.h5ad"), compression="gzip")
    print("\n[saved] myeloid_annotated.h5ad", flush=True)


if __name__ == "__main__":
    main()
