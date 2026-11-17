"""
GSE202051 snRNA-seq step 1 — data inspection & metadata reconstruction.

Paper: Hwang et al., Nat Genet 2022
GEO: GSE202051, GSE199102
Repository: github.com/karthikj89/humanpdac

This script loads the total h5ad, inspects key columns, and writes
step1_stats.tsv + sample_metadata_clean.tsv for downstream steps.
"""

import os, sys, time, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc

BASE    = os.path.dirname(os.path.abspath(__file__))
DATA    = os.path.join(BASE, "data")
OUT     = os.path.join(BASE, "processed")
FIG     = os.path.join(BASE, "figures")
for d in [DATA, OUT, FIG]:
    os.makedirs(d, exist_ok=True)

H5AD = os.path.join(DATA, "GSE202051_total.h5ad")
if not os.path.exists(H5AD):
    gz = os.path.join(DATA, "GSE202051_total.h5ad.gz")
    if not os.path.exists(gz):
        print("[ERROR] No h5ad or h5ad.gz found in data/", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] Decompressing {gz} ...", flush=True)
    t0 = time.time()
    import gzip
    with gzip.open(gz, "rb") as gz_f, open(H5AD, "wb") as out_f:
        for i, chunk in enumerate(iter(lambda: gz_f.read(16 << 20), b"")):
            out_f.write(chunk)
            if i % 30 == 0:
                print(f"  {i} chunks ({round(i*16/1024, 1)} GB)", flush=True)
    print(f"[OK] Decompressed in {time.time()-t0:.0f}s -> {os.path.getsize(H5AD)/1e9:.2f} GB", flush=True)

# ---------------------------------------------------------------------------
# Load (backed mode first for quick inspect, then full load)
# ---------------------------------------------------------------------------
print("\n[STEP 1] Loading AnnData (backed)...", flush=True)
t0 = time.time()
adata = sc.read_h5ad(H5AD, backed="r")
print(f"  Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes  ({time.time()-t0:.1f}s)", flush=True)
print(f"  Layers: {list(adata.layers.keys())}", flush=True)

# ---------------------------------------------------------------------------
# Key metadata inspection
# ---------------------------------------------------------------------------
print("\n=== Treatment distribution ===", flush=True)
for col in ["treatment_status", "new_treatment"]:
    if col in adata.obs.columns:
        vc = adata.obs[col].value_counts()
        print(f"\n  {col}:")
        for k, v in vc.items():
            print(f"    {str(k):25s}  n={v:>7,d}  ({v/len(adata)*100:.1f}%)")

# Patient-level treatment (each pid has one treatment)
if "pid" in adata.obs.columns and "treatment_status" in adata.obs.columns:
    pid_trt = adata.obs.groupby("pid")["treatment_status"].first()
    print(f"\n  Patients per treatment:")
    for k, v in pid_trt.value_counts().items():
        print(f"    {str(k):25s}  {v} patients")

print("\n=== Broad cell type distribution ===", flush=True)
if "broad_celltypes" in adata.obs.columns:
    vc = adata.obs["broad_celltypes"].value_counts()
    for k, v in vc.items():
        print(f"  {str(k):25s}  n={v:>7,d}  ({v/len(adata)*100:.1f}%)")

print("\n=== leiden clusters ===", flush=True)
if "leiden" in adata.obs.columns:
    vc = adata.obs["leiden"].value_counts().sort_index()
    for k, v in vc.items():
        print(f"  cluster {k:>3s}  n={v:>6,d}")

print("\n=== QC metrics (quick scan) ===", flush=True)
sc.pp.calculate_qc_metrics(adata, inplace=True)
qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
for c in qc_cols:
    if c in adata.obs.columns:
        s = adata.obs[c]
        print(f"  {c:25s}  median={s.median():.0f}  "
              f"Q1={s.quantile(0.25):.0f}  Q3={s.quantile(0.75):.0f}  "
              f"min={s.min():.0f}  max={s.max():.0f}")

# ---------------------------------------------------------------------------
# SERPIN / MARCO gene presence check
# ---------------------------------------------------------------------------
print("\n=== Key gene presence ===", flush=True)
key_genes = ["SERPINE1", "SERPINB2", "MARCO", "SPP1", "CD68", "CD163",
             "EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "PECAM1"]
gene_info = []
for g in key_genes:
    present = g in adata.var_names
    expr_vals = None
    pct = None
    if present:
        idx = list(adata.var_names).index(g)
        counts = np.asarray(adata[:, idx].X.todense()).ravel() if hasattr(adata.X, 'todense') else np.asarray(adata.X[:, idx].toarray()).ravel()
        expr_vals = float(counts.mean())
        pct = float((counts > 0).mean() * 100)
    gene_info.append({
        "gene": g,
        "present": present,
        "mean_expr": round(expr_vals, 4) if expr_vals is not None else None,
        "pct_cells": round(pct, 2) if pct is not None else None,
    })
pd.DataFrame(gene_info).to_csv(os.path.join(OUT, "step1_gene_presence.tsv"), sep="\t", index=False)
print(pd.DataFrame(gene_info).to_string(index=False))

# ---------------------------------------------------------------------------
# Save summary stats
# ---------------------------------------------------------------------------
stats = {
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_patients": int(adata.obs["pid"].nunique()) if "pid" in adata.obs.columns else None,
    "has_counts_layer": "counts" in adata.layers,
    "has_raw": hasattr(adata, "raw"),
    "n_leiden_clusters": int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs.columns else None,
    "time_s": round(time.time() - t0, 1),
}
pd.DataFrame([stats]).to_csv(os.path.join(OUT, "step1_stats.tsv"), sep="\t", index=False)

print(f"\n[OK] step1_stats.tsv written")
print(f"[OK] step1_gene_presence.tsv written")
print(f"\n=== STEP 1 COMPLETE ===")
