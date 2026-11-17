"""Decisive test of the metadata error in GSE289138.

GEO sample titles (esummary):
  GSM8784900 CF111 (classic) = control KO      GSM9575309 CF111_Flex = control KO
  GSM8784901 CF112 (classic) = Serpinb2 KO     GSM9575310 CF112_Flex = Serpinb2 KO
  GSM8784902 CF113 (classic) = Serpine1 KO     GSM9575311 CF113_Flex = Serpine1 KO

Our preprocessing hardcoded the SAME six libraries as
  classic = Vehicle, Flex = anti-PD-1

If the GEO titles are right, Serpine1 must be near-absent in CF113 and
Serpinb2 near-absent in CF112, in BOTH chemistries.
"""
import numpy as np
import pandas as pd
import anndata as ad

p = r"E:\pancreatic_cancer\submit\serpin_repro\processed\GSE289138_annotated.h5ad"
a = ad.read_h5ad(p, backed="r")

def det(g):
    x = a[:, g].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    return pd.Series(x > 0, index=a.obs_names)

obs = a.obs
print("=== per-library detection rate (%) of the two perturbed genes ===")
tab = pd.DataFrame({
    "Serpine1": det("Serpine1"),
    "Serpinb2": det("Serpinb2"),
    "Spp1": det("Spp1"),
    "Marco": det("Marco"),
})
tab["lib"] = obs["filename"].values
tab["chem"] = obs["condition"].values
tab["sample"] = obs["sample"].values
res = tab.groupby(["sample", "chem"], observed=True)[["Serpine1", "Serpinb2", "Spp1", "Marco"]].mean() * 100
print(res.round(2).to_string())

print("\n=== epithelial cells only, mean counts of perturbed gene ===")
epi = obs["anno"].astype(str).eq("Epithelial").values
for g in ["Serpine1", "Serpinb2"]:
    x = a[:, g].X
    x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    s = pd.Series(x[epi], index=obs.loc[epi, "filename"].values)
    print(f"  {g}: per-library mean count among epithelial cells")
    print("   ", s.groupby(level=0).mean().round(3).to_dict())

print("\n=== interpretation ===")
print("If Serpine1 ~0 in CF113 and Serpinb2 ~0 in CF112, the libraries are GENOTYPE KOs,")
print("not Vehicle/anti-PD-1, and the hardcoded 'treatment' column is fabricated.")
a.file.close()
