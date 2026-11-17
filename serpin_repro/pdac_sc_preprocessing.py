"""
GSE289138 (KPC mouse PDAC, Falcomata et al. Nature 2026) ScRNAseq preprocessing.

Pipeline mirrors Falcomata_PDAC_2026/ScRNAseq/1_human_scRNAseq_example_pre-process.ipynb:
  1. per-sample 10x mtx -> h5ad
  2. SOLO doublet detection via scVI (per sample)
  3. MAD-based outlier QC (total_counts / n_genes / top20 pct) + mt outlier
  4. concat all samples, HVG 6000 (seurat_v3, on counts layer)
  5. scVI integration with batch_key + continuous covariates
  6. neighbors -> UMAP -> leiden -> rank_genes_groups
  7. marker-based annotation for mouse (Serpine1/Serpinb2 axis, SPP1+MARCO+ macrophages)

Environment: venv (Python 3.10) + scanpy 1.11.5 + scvi-tools 1.3.3 + torch 2.11.0+cu128
Run: python pdac_sc_preprocessing.py
"""

import os
import re
import sys
import time
import warnings

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

warnings.filterwarnings("ignore")

BASE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE, "data", "GSE289138")
OUT_DIR = os.path.join(BASE, "processed")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

try:
    import torch

    DEVICE = "gpu" if torch.cuda.is_available() else "cpu"
except Exception:
    DEVICE = "cpu"

sc.settings.verbosity = 2
sc.settings.figdir = FIG_DIR
sc.set_figure_params(figsize=(5, 4), dpi=100, facecolor="white", frameon=False)
scvi.settings.seed = 0

# CellTypist expects human symbols; mouse data keeps mouse symbols.
RUNTIME = {"device": DEVICE}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------
# Metadata for GSE289138.
#
# CORRECTED 2026-11-14.  The previous version of this table was wrong: it
# labelled the classic libraries "Vehicle" and the Flex libraries "anti-PD-1",
# and wrote both into an obs column called `treatment`.  GSE289138 contains no
# treatment arms at all.  Authoritative source: the GEO series matrix header
# (GSE289138_series_matrix.txt.gz):
#
#   !Sample_title  "Pancreatic tumor and TME cells, control KO"
#                  "Pancreatic tumor and TME cells, Serpinb2 KO"
#                  "Pancreatic tumor and TME cells, Serpine1 KO"
#                  "Pancreatic tumor and TME cells, control KO - Flex"
#                  "Pancreatic tumor and TME cells, Serpinb2 KO - Flex"
#                  "Pancreatic tumor and TME cells, Serpine1 KO - Flex"
#   !Sample_characteristics_ch1  "cell line: KPC FC1245 orthotopically injected"
#
# So the design is 3 genotypes x 2 library chemistries (3' classic vs Flex),
# one tumor per genotype, and the two chemistries are technical replicates of
# the same tumor.  There is NO anti-PD-1 contrast.  Any Vehicle-vs-anti-PD-1
# result produced from this cohort is invalid; see CORRECTIONS.md C4.
# --------------------------------------------------------------------------
SAMPLE_META = {
    "GSM8784900_CF111": {"sample": "CF111", "condition": "classic", "genotype": "control_KO"},
    "GSM8784901_CF112": {"sample": "CF112", "condition": "classic", "genotype": "Serpinb2_KO"},
    "GSM8784902_CF113": {"sample": "CF113", "condition": "classic", "genotype": "Serpine1_KO"},
    "GSM9575309_CF111_Flex": {"sample": "CF111", "condition": "Flex", "genotype": "control_KO"},
    "GSM9575310_CF112_Flex": {"sample": "CF112", "condition": "Flex", "genotype": "Serpinb2_KO"},
    "GSM9575311_CF113_Flex": {"sample": "CF113", "condition": "Flex", "genotype": "Serpine1_KO"},
}

# A tumor ID must map to exactly one genotype; if this ever fails, the table
# is wrong and no downstream analysis should run.
_GENO_BY_SAMPLE = {}
for _lib, _m in SAMPLE_META.items():
    _prev = _GENO_BY_SAMPLE.setdefault(_m["sample"], _m["genotype"])
    assert _prev == _m["genotype"], (
        f"sample {_m['sample']} assigned two genotypes: {_prev} vs {_m['genotype']}"
    )
assert "treatment" not in next(iter(SAMPLE_META.values())), "GSE289138 has no treatment arms"


def enforce_metadata(adata, name):
    """Re-assert the corrected design on a cached per-sample object.

    The per-sample cache was written by an earlier run whose `read_samples`
    injected a fabricated `treatment` column.  The cache is reused for speed,
    so the invalid column has to be removed explicitly instead of inherited,
    and the genotype/chemistry labels have to be re-derived from SAMPLE_META.
    """
    meta = SAMPLE_META[name]
    if "treatment" in adata.obs.columns:
        adata.obs.drop(columns=["treatment"], inplace=True)
    assert "treatment" not in adata.obs.columns, f"treatment survived in {name}"
    assert (adata.obs["sample"].astype(str) == meta["sample"]).all(), \
        f"{name}: sample label mismatch"
    assert (adata.obs["condition"].astype(str) == meta["condition"]).all(), \
        f"{name}: chemistry label mismatch"
    adata.obs["genotype"] = meta["genotype"]
    return adata


# --------------------------------------------------------------------------
# 1. Read each sample
# --------------------------------------------------------------------------
def read_samples():
    adatas = {}
    for name in sorted(SAMPLE_META):
        path = os.path.join(DATA_DIR, name)
        if not os.path.isdir(path):
            log(f"SKIP missing: {name}")
            continue
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
        adata.var_names_make_unique()
        meta = SAMPLE_META[name]
        adata.obs["filename"] = name
        adata.obs["sample"] = meta["sample"]
        adata.obs["condition"] = meta["condition"]   # library chemistry, NOT treatment
        adata.obs["genotype"] = meta["genotype"]
        adata.obs_names = [f"{name}_{bc}" for bc in adata.obs_names]
        adatas[name] = adata
        log(f"read {name}: {adata.n_obs} cells x {adata.n_vars} genes")
    return adatas


# --------------------------------------------------------------------------
# 2. SOLO doublet detection (per sample), mirroring notebook
# --------------------------------------------------------------------------
def run_solo(adata, name):
    log(f"SOLO doublet detection on {name} ({adata.n_obs} cells)")
    a = adata.copy()
    sc.pp.filter_genes(a, min_cells=10)
    sc.pp.filter_genes(a, max_counts=100000)
    sc.pp.filter_cells(a, min_genes=10)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, subset=True, flavor="seurat_v3")
    scvi.model.SCVI.setup_anndata(a)
    vae = scvi.model.SCVI(a)
    vae.train(accelerator=DEVICE, devices=1, max_epochs=100)
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train(accelerator=DEVICE, devices=1)
    soft = solo.predict(soft=True)
    hard = solo.predict(soft=False)
    hard = np.asarray(hard).ravel()
    prob_d = (
        np.asarray(soft["doublet"]).ravel()
        if "doublet" in soft.columns
        else np.zeros(len(hard))
    )
    prob_s = (
        np.asarray(soft["singlet"]).ravel()
        if "singlet" in soft.columns
        else np.zeros(len(hard))
    )
    pred = np.where(hard == "doublet", "doublet", "singlet")
    # robust criterion: predicted doublet AND (probability margin > 0 OR high doublet prob)
    margin = prob_d - prob_s
    keep_doublet = (pred == "doublet") & ((margin > 0) | (prob_d > 0.5))
    doublet_idx = np.asarray(soft.index)[keep_doublet]
    log(
        f"  {name}: predicted_doublet={int((pred == 'doublet').sum())} "
        f"confirmed_doublet={len(doublet_idx)}"
    )
    adata.obs["doublet"] = adata.obs.index.isin(doublet_idx)
    adata.obs["doublet_prob"] = prob_d
    return adata


# --------------------------------------------------------------------------
# 3. QC with MAD outliers (mouse gene prefixes)
# --------------------------------------------------------------------------
def is_outlier(adata, metric, nmads):
    m = adata.obs[metric]
    return (m < np.median(m) - nmads * 1.4826 * np.median(np.abs(m - np.median(m)))) | (
        np.median(m) + nmads * 1.4826 * np.median(np.abs(m - np.median(m))) < m
    )


def qc_sample(adata, name):
    sc.pp.filter_cells(adata, min_genes=200)
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.match(r"^Rp[sl]\d")
    adata.var["hb"] = adata.var_names.str.match(r"^Hb[ab]-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True
    )
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 20)
    n0 = adata.n_obs
    adata = adata[(~adata.obs["outlier"]) & (~adata.obs["mt_outlier"])].copy()
    log(f"  {name}: QC {n0} -> {adata.n_obs} cells")
    return adata


# --------------------------------------------------------------------------
# 4. Annotation helpers (mouse symbols)
# --------------------------------------------------------------------------
MARKERS = {
    "Epithelial": ["Epcam", "Krt8", "Krt18", "Krt19", "Krt7", "Krt5", "Krt6a", "S100p", "S100a6", "Krt20"],
    "T/NK": ["Cd3e", "Cd3d", "Cd3g", "Trbc1", "Trbc2", "Cd8a", "Cd4", "Nkg7", "Klrd1", "Gzma", "Gzmb", "Prf1", "Ifng"],
    "Macrophage": ["Adgre1", "C1qa", "C1qb", "C1qc", "Cd68", "Lyz2", "Marco", "Spp1", "Lgals3", "Cd163", "Mrc1", "Ctss"],
    "Monocyte": ["Ly6c2", "Ccr2", "Plac8", "Chil3", "F13a1", "Ms4a4c", "Ifitm6"],
    "B cell": ["Cd79a", "Cd79b", "Ms4a1", "Ebf1", "Cd19", "Ighm", "Ighd"],
    "Plasma cell": ["Mzb1", "Jchain", "Sdc1", "Xbp1", "Derl3", "Igkc"],
    "CAF/Stromal": ["Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Dcn", "Lum", "Acta2", "Tagln", "Postn", "Ddr2"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Esam", "Egfl7", "Emcn", "Cldn5"],
    "Mast cell": ["Kit", "Cpa3", "Mcpt4", "Fcer1a", "Ms4a2", "Gata2"],
    "Neutrophil": ["S100a8", "S100a9", "Ly6g", "Retnlg", "Csf3r", "Mmp8", "Mmp9"],
    "pDC": ["Siglech", "Bst2", "Ccr9", "Irf8", "Irf7", "Cox6a2"],
    "RBC": ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2"],
    "DC": ["Flt3", "Cd74", "H2-Ab1", "Cd83", "Ccr7", "Fscn1"],
}


def score_and_annotate(adata):
    for ct, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.tl.score_genes(adata, present, score_name=f"score_{ct}", random_state=0)
    score_cols = [f"score_{c}" for c in MARKERS if f"score_{c}" in adata.obs.columns]
    score_mat = adata.obs[score_cols].to_numpy()
    score_mat = np.nan_to_num(score_mat, nan=-np.inf)
    best = np.argmax(score_mat, axis=1)
    adata.obs["anno_auto"] = [score_cols[i].replace("score_", "") for i in best]
    adata.obs["anno_auto_score"] = score_mat[np.arange(len(best)), best]
    return adata


def main():
    log(f"Device: {DEVICE}")
    adatas = read_samples()
    if not adatas:
        log("No samples found, aborting")
        sys.exit(1)

    # --- per-sample SOLO + QC (cached) ---
    cache_dir = os.path.join(OUT_DIR, "per_sample_cache")
    os.makedirs(cache_dir, exist_ok=True)
    processed = []
    for name, adata in adatas.items():
        cache_path = os.path.join(cache_dir, f"{name}_qc.h5ad")
        if os.path.exists(cache_path):
            log(f"=== {name}: loading cached QC result ===")
            processed.append(enforce_metadata(sc.read_h5ad(cache_path), name))
            continue
        log(f"=== processing {name} ===")
        adata = run_solo(adata, name)
        adata = qc_sample(adata, name)
        adata.write(cache_path)
        processed.append(adata)
        del adata
        import gc

        gc.collect()

    adata = ad.concat(processed, join="outer", label="batch", index_unique="-")
    adata.var_names_make_unique()
    del processed
    import gc

    gc.collect()
    assert "treatment" not in adata.obs.columns, \
        "Fatal: fabricated treatment column survived concat"
    assert set(adata.obs["genotype"].astype(str)) == {
        "control_KO", "Serpinb2_KO", "Serpine1_KO"}, \
        "Fatal: unexpected genotype labels after concat"
    log(f"CONCAT: {adata.n_obs} cells x {adata.n_vars} genes")
    log("Doublet/QC summary:")
    for f in ["filename", "condition", "genotype", "sample"]:
        print(adata.obs[f].value_counts().to_string())
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", color="condition", save="_qc_mt.png")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="condition", save="_qc_genes.png")

    # --- integration ---
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=6000, layer="counts", subset=True, batch_key="filename"
    )
    log(f"HVG selected: {adata.n_vars}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="filename",
        continuous_covariate_keys=["pct_counts_mt", "total_counts", "pct_counts_ribo"],
    )
    model_dir = os.path.join(OUT_DIR, "scvi_model")
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
    if os.path.exists(os.path.join(model_dir, "model.pt")):
        log("Loading cached scVI model from " + model_dir)
        model = scvi.model.SCVI.load(model_dir, adata=adata)
    else:
        log("Training scVI on " + DEVICE)
        model.train(accelerator=DEVICE, devices=1, max_epochs=300, early_stopping=True)
        model.save(model_dir, overwrite=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scVI_normalized"] = model.get_normalized_expression(library_size=1e4)
    log("scVI latent extracted")

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1, key_added="leiden")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    log("Leiden clusters: " + str(adata.obs["leiden"].value_counts().sort_index().to_dict()))

    adata = score_and_annotate(adata)

    sc.pl.umap(adata, color=["condition", "filename", "leiden", "anno_auto"], save="_scvi_overview.png", ncols=2)
    sc.pl.umap(adata, color=["Serpine1", "Serpinb2", "Spp1", "Marco"], save="_serpin_axis.png", ncols=4)
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby="leiden", save="_marker_heatmap.png")

    out_h5ad = os.path.join(OUT_DIR, "GSE289138_processed.h5ad")
    adata.write(out_h5ad)
    log(f"WROTE {out_h5ad}")

    # per-cluster Serpin summary (the paper's core readout)
    summary = (
        adata.obs.groupby(["leiden", "anno_auto"], observed=True)
        .agg(n=("leiden", "size"))
        .reset_index()
    )
    norm = adata.layers["scVI_normalized"]
    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    for gene in ["Serpine1", "Serpinb2", "Spp1", "Marco"]:
        if gene in var_idx:
            col = np.asarray(norm[:, var_idx[gene]]).ravel()
            tmp = adata.obs[["leiden", "anno_auto"]].copy()
            tmp["expr"] = col
            g = tmp.groupby(["leiden", "anno_auto"], observed=True)["expr"].agg(["mean", "size"])
            g.columns = [f"{gene}_{c}" for c in g.columns]
            summary = summary.merge(g.reset_index(), on=["leiden", "anno_auto"], how="left")
    summary = summary.sort_values("n", ascending=False)
    summary.to_csv(os.path.join(OUT_DIR, "cluster_serpin_summary.csv"), index=False)
    log("Serpin summary:")
    print(summary.to_string(index=False))

    log("=== DONE ===")


if __name__ == "__main__":
    main()
