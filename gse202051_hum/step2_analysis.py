"""
GSE202051 snRNA-seq step 2 — QC filtering, scVI integration, SERPIN/MARCO axis quantification.

Builds on step1. Outputs:
  - processed/GSE202051_qc_filtered.h5ad          (filtered + normalized)
  - processed/scvi_model/                         (scVI model)
  - processed/patient_level_serpin_marco.tsv      (per-patient SERPIN+/MARCO+ fractions)
  - processed/comparison_summary.tsv              (naive vs treated summary)
"""

import os, sys, time, warnings
warnings.filterwarnings("ignore")

import torch
torch.set_float32_matmul_precision('high')
import pandas as pd
import scanpy as sc
import scvi
import anndata as ad

BASE    = os.path.dirname(os.path.abspath(__file__))
DATA    = os.path.join(BASE, "data")
OUT     = os.path.join(BASE, "processed")
FIG     = os.path.join(BASE, "figures")
for d in [OUT, FIG]:
    os.makedirs(d, exist_ok=True)

H5AD   = os.path.join(DATA, "GSE202051_total.h5ad")
PARENT = r"E:\pancreatic_cancer\submit\gse202051_hum"

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
print("[STEP 2] Loading data...", flush=True)
t0 = time.time()
adata = sc.read_h5ad(H5AD)
print(f"  Loaded: {adata.n_obs:,} x {adata.n_vars:,}  ({time.time()-t0:.1f}s)", flush=True)

# ---------------------------------------------------------------------------
# 2a. Patient-level status mapping  (pid → unified treatment)
# ---------------------------------------------------------------------------
# Map original treatment_status → unified: naive / treated
TRT_MAP = {
    "Untreated":  "naive",
    "CRT":        "treated",
    "CRTl":       "treated",
    "CRTx":       "treated",
    "CRTn":       "treated",
    "RT":         "treated",
    "GART":       "treated",
    "CRTln":      "treated",
}
adata.obs["status"] = adata.obs["treatment_status"].map(TRT_MAP).fillna("unknown")
print(f"\n  status distribution:")
print(f"    naive : {(adata.obs['status']=='naive').sum():,} cells  "
      f"[{(adata.obs['status']=='naive').sum()/len(adata)*100:.1f}%]")
print(f"    treated: {(adata.obs['status']=='treated').sum():,} cells  "
      f"[{(adata.obs['status']=='treated').sum()/len(adata)*100:.1f}%]")

# Per-patient status
pid_status = adata.obs.groupby("pid")["status"].first()
print(f"  Patients: naive={int((pid_status=='naive').sum())}, treated={int((pid_status=='treated').sum())}")

# ---------------------------------------------------------------------------
# 2b. QC (snRNA-seq thresholds)
# ---------------------------------------------------------------------------
# QC metrics are pre-computed in the original h5ad; reuse them directly
qc_cols_present = [c for c in ["n_genes_by_counts","total_counts","pct_counts_mt"] if c in adata.obs.columns]
if qc_cols_present:
    print(f"  Using pre-computed QC cols: {qc_cols_present}", flush=True)
    for c in qc_cols_present:
        s = adata.obs[c]
        print(f"    {c}: median={s.median():.0f}  Q1={s.quantile(.25):.0f}  Q3={s.quantile(.75):.0f}")
else:
    # Compute manually
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    print("  QC computed de novo", flush=True)

# Apply snRNA-seq-friendly thresholds
adata = adata[(adata.obs["n_genes_by_counts"] >= 200) &
              (adata.obs["n_genes_by_counts"] <= 6000)].copy()
print(f"  After QC: {adata.n_obs:,} nuclei")

# ---------------------------------------------------------------------------
# 2c. Normalize, HVG
# ---------------------------------------------------------------------------
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="seurat_v3",
    layer="counts", batch_key="pid",
)
print(f"  HVGs: {adata.n_vars} selected", flush=True)

# ---------------------------------------------------------------------------
# 2d. scVI integration  (batch=pid, covariate=status to avoid confounding)
# ---------------------------------------------------------------------------
# Subsample to ~80K cells for faster scVI (preserves patient/treatment balance)
import numpy as np
np.random.seed(0)
n_target = 80000
if adata.n_obs > n_target:
    idx = np.random.choice(adata.n_obs, n_target, replace=False)
    adata = adata[idx].copy()
    print(f"  Subsampled to {adata.n_obs:,} cells for scVI", flush=True)

print(f"\n  Setting up scVI ...", flush=True)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="pid",
    categorical_covariate_keys=["status"],   # protect treatment signal from batch correction
)
model = scvi.model.SCVI(
    adata, n_layers=2, n_latent=30,
    gene_likelihood="nb",
)
print("  Training scVI...", flush=True)
train_t0 = time.time()
model.train(max_epochs=100, early_stopping=True, early_stopping_patience=15,
            batch_size=8192)
print(f"  scVI trained in {time.time()-train_t0:.0f}s", flush=True)

adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_scvi")
print(f"  UMAP + Leiden done. Clusters: {adata.obs['leiden_scvi'].nunique()}", flush=True)

# Save scVI model + filtered adata
model.save(os.path.join(OUT, "scvi_model"), overwrite=True)
adata.write(os.path.join(OUT, "GSE202051_qc_filtered.h5ad"))
print(f"  Saved to {OUT}", flush=True)

# ---------------------------------------------------------------------------
# 2e. Define SERPIN-high malignant epithelial cells
# ---------------------------------------------------------------------------
# Malignant epithelial = high Moffitt_basal + Bailey_progenitor score + EPCAM/KRT high
epi_genes = [g for g in ["EPCAM","KRT8","KRT18","KRT19"] if g in adata.var_names]
malignant_markers = ["Moffitt_basal", "Bailey_progenitor", "Collison_QM"]
existing_scores = [g for g in malignant_markers if g in adata.obs.columns]
if existing_scores:
    adata.obs["malignancy_score"] = adata.obs[existing_scores].mean(axis=1)
else:
    # Compute de novo
    mal_genes = []
    for glist in [["MAL2", "KRT8", "EPCAM"], ["VIM", "ACTA2"]]:
        mal_genes += [g for g in glist if g in adata.var_names]
    if mal_genes:
        sc.tl.score_genes(adata, mal_genes, score_name="malignancy_score")
    else:
        adata.obs["malignancy_score"] = 0.0

# Epithelial mask (use broad_celltypes annotation from original paper)
if "broad_celltypes" in adata.obs.columns:
    epi_mask = adata.obs["broad_celltypes"] == "Epithelial"
else:
    epi_mask = pd.Series(False, index=adata.obs_names)

mal_mask = epi_mask & (adata.obs["malignancy_score"] > adata.obs["malignancy_score"].median())
print(f"\n  Epithelial: {epi_mask.sum():,}  |  Malignant-ish: {mal_mask.sum():,}")

# SERPINE1 / SERPINB2 scoring (on ALL cells; will filter to malignant epi later)
serpin_genes = [g for g in ["SERPINE1", "SERPINB2"] if g in adata.var_names]
if serpin_genes:
    sc.tl.score_genes(adata, serpin_genes, score_name="SERPIN_score")
    thr = float(adata[mal_mask].obs["SERPIN_score"].median())
    adata.obs.loc[mal_mask, "SERPIN_status"] = np.where(
        adata.obs.loc[mal_mask, "SERPIN_score"] > thr, "SERPIN_high", "SERPIN_low"
    )
else:
    adata.obs["SERPIN_status"] = "non_malignant"
adata.obs["SERPIN_status"] = adata.obs["SERPIN_status"].fillna("non_malignant")

# Print SERPIN stats
n_serpin_high = (adata.obs["SERPIN_status"] == "SERPIN_high").sum()
print(f"  SERPIN_high malignant cells: {n_serpin_high:,}  "
      f"[{n_serpin_high/mal_mask.sum()*100:.1f}% of malignant epithelial]")

# ---------------------------------------------------------------------------
# 2f. Myeloid/MARCO definition
# ---------------------------------------------------------------------------
myeloid_genes = [g for g in ["CD68","CD163","MARCO","C1QA","CSF1R"] if g in adata.var_names]
immune_mask = adata.obs["broad_celltypes"] == "Immune" if "broad_celltypes" in adata.obs.columns \
    else pd.Series(False, index=adata.obs_names)
if myeloid_genes and immune_mask.any():
    # Score on ALL cells (view bug fix)
    sc.tl.score_genes(adata, myeloid_genes, score_name="Myeloid_score")
    mye_thr = float(adata[immune_mask].obs["Myeloid_score"].median())
    myeloid_mask = immune_mask & (adata.obs["Myeloid_score"] > mye_thr)
else:
    # Use pre-computed scores if available
    mye_cols = [c for c in adata.obs.columns if c in
                ["Monocytes_1","Monocytes_2","Macrophage","M0","M1","M2"]]
    if mye_cols:
        myeloid_mask = adata.obs[mye_cols].any(axis=1)
    else:
        myeloid_mask = pd.Series(False, index=adata.obs_names)

print(f"  Myeloid cells: {myeloid_mask.sum():,}")

# MARCO expression (log1p space after normalization)
marco_gene = "MARCO" if "MARCO" in adata.var_names else None
if marco_gene:
    marco_expr = adata[:, marco_gene].X.toarray().ravel()
    adata.obs["MARCO_expr"] = marco_expr
    # Threshold: top quartile of myeloid cells
    mye_expr = marco_expr[myeloid_mask]
    mar_thr = np.percentile(mye_expr, 75) if len(mye_expr) > 10 else mye_expr.median()
    adata.obs["MARCO_pos"] = (marco_expr > mar_thr) & myeloid_mask
    n_marco_pos = int(adata.obs["MARCO_pos"].sum())
    print(f"  MARCO+ myeloid cells: {n_marco_pos:,}  [{n_marco_pos/myeloid_mask.sum()*100:.1f}% of myeloid]")
else:
    adata.obs["MARCO_pos"] = False

# ---------------------------------------------------------------------------
# 2g. Patient-level SERPIN/MARCO fractions → comparison
# ---------------------------------------------------------------------------
print("\n  Computing patient-level fractions...", flush=True)
patient_rows = []
for pid, grp in adata.obs.groupby("pid"):
    st = grp["status"].iloc[0]
    n_ep = grp["is_epithelial"].sum() if "is_epithelial" in grp.columns else mal_mask[grp.index].sum()
    n_mal = grp["SERPIN_status"].isin(["SERPIN_high","SERPIN_low"]).sum()
    serpin_hi = (grp["SERPIN_status"] == "SERPIN_high").sum()
    serpin_frac = serpin_hi / n_mal if n_mal > 0 else np.nan
    n_imm = grp["is_myeloid"].sum() if "is_myeloid" in grp.columns else myeloid_mask[grp.index].sum()
    marco_pos = grp["MARCO_pos"].sum()
    marco_frac = marco_pos / n_imm if n_imm > 0 else np.nan
    patient_rows.append({
        "pid": pid,
        "status": st,
        "n_cells": len(grp),
        "n_malignant_epi": n_mal,
        "n_SERPIN_high": serpin_hi,
        "SERPIN_frac": round(serpin_frac, 4),
        "n_myeloid": n_imm,
        "n_MARCO_pos": marco_pos,
        "MARCO_frac": round(marco_frac, 4),
    })

pat_df = pd.DataFrame(patient_rows)
pat_df.to_csv(os.path.join(OUT, "patient_level_serpin_marco.tsv"), sep="\t", index=False)
print(f"  {len(pat_df)} patients loaded")
print(pat_df.groupby("status")[["SERPIN_frac","MARCO_frac","n_SERPIN_high","n_MARCO_pos"]].agg(["mean","median","count"]).to_string())

# Wilcoxon on patient-level fractions (paired where possible)
naive_pids = pat_df[pat_df["status"]=="naive"]["pid"].tolist()
treated_pids = pat_df[pat_df["status"]=="treated"]["pid"].tolist()
from scipy import stats

def paired_test(naive_vals, treated_vals):
    common = set(naive_vals.index) & set(treated_vals.index)
    if len(common) < 3: return None, None, "insufficient_pairs"
    v1 = naive_vals.loc[list(common)]
    v2 = treated_vals.loc[list(common)]
    stat, p = stats.wilcoxon(v1.values, v2.values)
    return stat, p, f"n_pair={len(common)}"

print("\n=== STEP 2 SUMMARY ===", flush=True)
print(f"Total cells (QC-filtered): {len(adata):,}", flush=True)
print(f"Malignant epithelial: {mal_mask.sum():,}", flush=True)
print(f"SERPIN-high: {(adata.obs['SERPIN_status']=='SERPIN_high').sum():,}", flush=True)
print(f"Myeloid: {myeloid_mask.sum():,}", flush=True)
print(f"MARCO+: {int(adata.obs['MARCO_pos'].sum()):,}", flush=True)
print(f"\nSaved: patient_level_serpin_marco.tsv, GSE202051_qc_filtered.h5ad, scvi_model/",
      flush=True)
print("=== STEP 2 COMPLETE ===", flush=True)
