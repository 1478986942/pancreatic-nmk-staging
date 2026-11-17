"""Design check for the mouse cohort GSE289138: what can n=3 support?"""
import pandas as pd
import anndata as ad

p = r"E:\pancreatic_cancer\submit\serpin_repro\processed\GSE289138_annotated.h5ad"
a = ad.read_h5ad(p, backed="r")
obs = a.obs[["filename", "sample", "condition", "treatment"]].copy()
print("=== sample x condition x treatment (cells) ===")
print(obs.groupby(["sample", "condition", "treatment"], observed=True).size().to_string())

print("\n=== cells per filename ===")
print(obs.groupby("filename", observed=True).size().to_string())

print("\n=== annotated cell types per sample ===")
obs2 = a.obs[["sample", "anno", "malignant"]].copy()
print(obs2.groupby(["sample", "anno"], observed=True).size().unstack(fill_value=0).to_string())

print("\n=== malignant epithelial cells per sample ===")
epi = obs2[obs2["anno"] == "Epithelial"]
print(epi.groupby(["sample", "malignant"], observed=True).size().unstack(fill_value=0).to_string())

print("\n=== Serpine1 / Serpinb2 / Spp1 / Marco detection rate by anno ===")
cols = ["Serpine1_detected", "Serpinb2_detected", "Spp1_detected", "Marco_detected"]
det = a.obs[["anno"] + cols].copy()
for c in cols:
    det[c] = det[c].astype(str).eq("True")
print(det.groupby("anno", observed=True)[cols].mean().round(4).to_string())
a.file.close()
