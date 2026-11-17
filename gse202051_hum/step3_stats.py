"""
GSE202051 step 3 — apply trained scVI to full dataset + formal stats.

Loads the scVI model trained on 80K subsample, applies it to full 220K cells,
then computes patient-level SERPIN/MARCO fractions with statistical tests.
"""

import os, sys, time, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy import stats

BASE   = r"E:\pancreatic_cancer\submit\gse202051_hum"
DATA   = os.path.join(BASE, "data")
OUT    = os.path.join(BASE, "processed")
FIG    = os.path.join(BASE, "figures")

FULL_H5AD = os.path.join(DATA, "GSE202051_total.h5ad")
SCVI_MODEL_DIR = os.path.join(OUT, "scvi_model")

# Treatment mapping
TRT_MAP = {
    "Untreated": "naive",
    "CRT": "treated", "CRTl": "treated", "CRTx": "treated",
    "CRTn": "treated", "RT": "treated", "GART": "treated", "CRTln": "treated",
}
MARKER_GENES = {
    "SERPINE1": "SERPINE1",
    "SERPINB2": "SERPINB2",
    "MARCO": "MARCO",
    "SPP1": "SPP1",
    "CD68": "CD68",
    "CD163": "CD163",
}

print("[STEP 3] Loading full dataset...", flush=True)
t0 = time.time()
adata_full = sc.read_h5ad(FULL_H5AD)
print(f"  Full dataset: {adata_full.n_obs:,} x {adata_full.n_vars:,}  ({time.time()-t0:.1f}s)", flush=True)

# Add unified status column
if "status" not in adata_full.obs.columns or not {"naive","treated"}.intersection(adata_full.obs["status"].unique()):
    adata_full.obs["status"] = adata_full.obs["treatment_status"].map(TRT_MAP).fillna("unknown")
    pid_status = adata_full.obs.groupby("pid")["status"].first()
    print(f"  Status re-mapped: naive={int((pid_status=='naive').sum())} pts, treated={int((pid_status=='treated').sum())} pts", flush=True)

# --- Load HVG mask from step2 ---
# Step 2 used n_top_genes=2000 seurat_v3; load the filtered h5ad to get HVG indices
h5ad_filtered = os.path.join(OUT, "GSE202051_qc_filtered.h5ad")
if os.path.exists(h5ad_filtered):
    adata_sub = sc.read_h5ad(h5ad_filtered)
    # Get HVG names
    hvgs = set(adata_sub.var_names).intersection(adata_full.var_names)
    print(f"\n  Found HVGs from step2: {len(hvgs)} genes", flush=True)
else:
    hvgs = None

# --- Apply scVI model to full dataset (no training, just transform) ---
print("\n[STEP 3] Applying scVI latent to full dataset ...", flush=True)
if os.path.exists(SCVI_MODEL_DIR):
    # Re-setup anndata with same params
    scvi.model.SCVI.setup_anndata(
        adata_full, layer="counts", batch_key="pid",
        categorical_covariate_keys=["status"] if "status" in adata_full.obs.columns else None,
    )
    model = scvi.model.SCVI.load(SCVI_MODEL_DIR, adata=adata_full)
    adata_full.obsm["X_scVI"] = model.get_latent_representation()
    print(f"  scVI latent computed: {adata_full.obsm['X_scVI'].shape}", flush=True)

    # Neighbors + UMAP
    sc.pp.neighbors(adata_full, use_rep="X_scVI", n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_full)
    sc.tl.leiden(adata_full, resolution=0.5, key_added="leiden_scvi")
    print(f"  UMAP + Leiden done. Clusters: {adata_full.obs['leiden_scvi'].nunique()}", flush=True)
else:
    print("  WARNING: scVI model not found, skipping dimensionality reduction", flush=True)

# --- Define masks ---
epi_mask = adata_full.obs["broad_celltypes"] == "Epithelial" if "broad_celltypes" in adata_full.obs.columns \
    else pd.Series(False, index=adata_full.obs_names)
mal_thr = adata_full.obs["Moffitt_basal"].median() if "Moffitt_basal" in adata_full.obs.columns else 0
mal_mask = epi_mask & (adata_full.obs["Moffitt_basal"] > mal_thr) if "Moffitt_basal" in adata_full.obs.columns \
    else epi_mask
myeloid_cols = [c for c in adata_full.obs.columns if any(x in c.lower()
    for x in ["monocyte","macrophage","m0","m1","m2"]) and "score" not in c.lower()]
if not myeloid_cols:
    # Fallback: use broad_celltypes or immune-related score columns
    myeloid_cols = [c for c in ["Monocytes_1","Monocytes_2","Macrophage","M0","M1","M2"]
                    if c in adata_full.obs.columns]
myeloid_mask = adata_full.obs[myeloid_cols].any(axis=1) if myeloid_cols else pd.Series(False, index=adata_full.obs_names)

print(f"\n  Epithelial: {epi_mask.sum():,}  Malignant-ish: {mal_mask.sum():,}  Myeloid: {myeloid_mask.sum():,}", flush=True)

# --- SERPIN scoring ---
serpin_genes = [g for g in ["SERPINE1","SERPINB2"] if g in adata_full.var_names]
marco_gene = "MARCO" if "MARCO" in adata_full.var_names else None
spp1_gene  = "SPP1"  if "SPP1"  in adata_full.var_names else None

if serpin_genes:
    sc.tl.score_genes(adata_full, serpin_genes, score_name="SERPIN_score")
    thr = float(adata_full[mal_mask].obs["SERPIN_score"].median())
    adata_full.obs.loc[mal_mask, "SERPIN_status"] = np.where(
        adata_full.obs.loc[mal_mask, "SERPIN_score"] > thr, "SERPIN_high", "SERPIN_low"
    )
    adata_full.obs["SERPIN_status"] = adata_full.obs["SERPIN_status"].fillna("non_malignant")
    print(f"\n  SERPIN threshold (malignant epi median): {thr:.4f}")
else:
    adata_full.obs["SERPIN_status"] = "non_malignant"
    print("  WARNING: SERPIN genes not found!", flush=True)

# MARCO expression in myeloid cells
if marco_gene:
    marco_expr = np.asarray(adata_full[:, marco_gene].X.toarray()).ravel()
    adata_full.obs["MARCO_expr"] = marco_expr
    mye_expr = marco_expr[myeloid_mask]
    mar_thr = float(np.percentile(mye_expr, 75)) if len(mye_expr) > 10 else float(mye_expr.median())
    adata_full.obs["MARCO_pos"] = (marco_expr > mar_thr) & myeloid_mask
    print(f"  MARCO threshold (75th pct of myeloid): {mar_thr:.4f}")
else:
    adata_full.obs["MARCO_pos"] = False
    print("  WARNING: MARCO gene not found!", flush=True)

# SPP1 co-expression
spp1_pos = None
if spp1_gene:
    spp1_expr = np.asarray(adata_full[:, spp1_gene].X.toarray()).ravel()
    adata_full.obs["SPP1_expr"] = spp1_expr
    spp1_thr = float(np.percentile(mye_expr if myeloid_mask.any() else spp1_expr, 75))
    adata_full.obs["SPP1_pos"] = (spp1_expr > spp1_thr)
    spp1_pos = adata_full.obs["SPP1_pos"] & adata_full.obs["MARCO_pos"]

# --- Patient-level aggregation ---
print("\n=== Patient-Level SERPIN–MARCO Axis ===", flush=True)
patient_rows = []
for pid, grp in adata_full.obs.groupby("pid"):
    st = grp["status"].iloc[0]
    n_cells = len(grp)
    n_mal   = int(mal_mask[grp.index].sum())
    n_serpin_hi = int((grp["SERPIN_status"] == "SERPIN_high").sum())
    serpin_frac = n_serpin_hi / n_mal if n_mal > 0 else np.nan
    n_myel  = int(myeloid_mask[grp.index].sum())
    n_marco = int(grp["MARCO_pos"].sum())
    marco_frac = n_marco / n_myel if n_myel > 0 else np.nan
    n_spp1_marco = int(spp1_pos.loc[grp.index].sum()) if spp1_pos is not None else 0
    patient_rows.append({
        "pid": pid, "status": st, "n_cells": n_cells,
        "n_malignant_epi": n_mal, "n_SERPIN_high": n_serpin_hi,
        "SERPIN_frac": round(serpin_frac, 4),
        "n_myeloid": n_myel, "n_MARCO_pos": n_marco,
        "MARCO_frac": round(marco_frac, 4),
        "n_SPP1_MARCO_pp": n_spp1_marco,
    })

pat_df = pd.DataFrame(patient_rows)
pat_df.to_csv(os.path.join(OUT, "patient_level_serpin_marco_full.tsv"), sep="\t", index=False)

print(f"\n  Total patients: {len(pat_df)}")
print(f"  Naive: {(pat_df['status']=='naive').sum()}  Treated: {(pat_df['status']=='treated').sum()}")
print(f"\n  --- SERPIN_frac (malignant epi) ---")
for st in ["naive","treated"]:
    sub = pat_df[pat_df["status"]==st]
    valid = sub[sub["SERPIN_frac"].notna()]
    print(f"    {st} (n={len(valid)}): mean={valid['SERPIN_frac'].mean():.4f}  "
          f"median={valid['SERPIN_frac'].median():.4f}  "
          f"range=[{valid['SERPIN_frac'].min():.4f}, {valid['SERPIN_frac'].max():.4f}]")

print(f"\n  --- MARCO_frac (myeloid) ---")
for st in ["naive","treated"]:
    sub = pat_df[pat_df["status"]==st]
    valid = sub[sub["MARCO_frac"].notna()]
    print(f"    {st} (n={len(valid)}): mean={valid['MARCO_frac'].mean():.4f}  "
          f"median={valid['MARCO_frac'].median():.4f}  "
          f"range=[{valid['MARCO_frac'].min():.4f}, {valid['MARCO_frac'].max():.4f}]")

# Wilcoxon rank-sum (unpaired, n=18 vs n=25)
serpin_naive = pat_df[pat_df["status"]=="naive"]["SERPIN_frac"].dropna()
serpin_treated = pat_df[pat_df["status"]=="treated"]["SERPIN_frac"].dropna()
marco_naive  = pat_df[pat_df["status"]=="naive"]["MARCO_frac"].dropna()
marco_treated = pat_df[pat_df["status"]=="treated"]["MARCO_frac"].dropna()

stat_s, p_s = stats.mannwhitneyu(serpin_naive, serpin_treated, alternative="two-sided")
stat_m, p_m = stats.mannwhitneyu(marco_naive,  marco_treated, alternative="two-sided")

print(f"\n  Mann-Whitney U tests:")
print(f"    SERPIN_frac:     U={stat_s:.0f}  p={p_s:.4f}  (naive n={len(serpin_naive)}, treated n={len(serpin_treated)})")
print(f"    MARCO_frac:      U={stat_m:.0f}  p={p_m:.4f}  (naive n={len(marco_naive)},  treated n={len(marco_treated)})")

# Effect size (Cliff's delta)
def cliff_delta(a, b):
    n = len(a) * len(b)
    if n == 0: return np.nan
    greater = sum(x > y for x in a for y in b)
    less    = sum(x < y for x in a for y in b)
    return (greater - less) / n

es_serpin = cliff_delta(serpin_naive, serpin_treated)
es_marco  = cliff_delta(marco_naive,  marco_treated)
print(f"    Cliff's delta SERPIN_frac: {es_serpin:.3f}  (near 0 = no effect)")
print(f"    Cliff's delta MARCO_frac:  {es_marco:.3f}")

# Save final summary
summary = {
    "dataset": "GSE202051_human_snRNA",
    "total_cells": int(adata_full.n_obs),
    "total_patients": int(len(pat_df)),
    "naive_patients": int((pat_df["status"]=="naive").sum()),
    "treated_patients": int((pat_df["status"]=="treated").sum()),
    "n_malignant_epi_cells": int(mal_mask.sum()),
    "n_serpin_high_total": int((pat_df["n_SERPIN_high"]).sum()),
    "n_myeloid_total": int(myeloid_mask.sum()),
    "n_marco_pos_total": int(adata_full.obs["MARCO_pos"].sum()),
    "SERPIN_frac_naive_mean": round(float(serpin_naive.mean()), 4),
    "SERPIN_frac_treated_mean": round(float(serpin_treated.mean()), 4),
    "SERPIN_U": float(stat_s), "SERPIN_p": float(p_s), "SERPIN_cliff_delta": round(es_serpin, 4),
    "MARCO_frac_naive_mean": round(float(marco_naive.mean()), 4),
    "MARCO_frac_treated_mean": round(float(marco_treated.mean()), 4),
    "MARCO_U": float(stat_m), "MARCO_p": float(p_m), "MARCO_cliff_delta": round(es_marco, 4),
}
pd.DataFrame([summary]).to_csv(os.path.join(OUT, "step3_summary.tsv"), sep="\t", index=False)
print(f"\n[OK] step3_summary.tsv written")
print(f"[OK] patient_level_serpin_marco_full.tsv written")
print("\n=== STEP 3 COMPLETE ===", flush=True)
