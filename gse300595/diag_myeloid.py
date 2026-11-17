import os, sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, r"E:\pancreatic_cancer\submit\gse300595")
from io_utils import PROC

a = sc.read_h5ad(os.path.join(PROC, "rna_qc.h5ad"))
print("rna_qc genes:", a.shape[1])
for g in ["SERPINE1", "SERPINB2", "SERPINB1", "SERPINH1", "MARCO", "SPP1"]:
    if g in a.var_names:
        v = a[:, g].X
        v = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.asarray(v).ravel()
        print(f"{g:10s} present  detected_in={int((v>0).sum()):6d}  "
              f"({100*(v>0).mean():5.2f}%)  mean={v.mean():.3f}")
    else:
        print(f"{g:10s} ABSENT from count matrix")

b = sc.read_h5ad(os.path.join(PROC, "rna_annotated.h5ad"))
raw = b.raw.to_adata()
print("\nannotated raw genes:", raw.shape[1])
print("SERPINB2 in raw:", "SERPINB2" in raw.var_names)

# per-sample myeloid marker detection
mark = ["CD68", "AIF1", "C1QA", "CSF1R", "CD14", "LYZ", "ITGAM", "PTPRC"]
mark = [g for g in mark if g in raw.var_names]
sub = raw[:, mark]
X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
df = pd.DataFrame(X, columns=mark, index=raw.obs_names)
df["sample"] = b.obs["sample"].values
df["leiden"] = b.obs["leiden"].values
df["myeloid_score"] = df[mark].mean(axis=1)
g = df.groupby("sample")
print("\nPer-sample cells with myeloid_score>0.1 (log-normalised mean of 8 markers):")
res = pd.DataFrame({
    "n_cells": g.size(),
    "n_myeloid_like": g.apply(lambda d: int((d["myeloid_score"] > 0.1).sum()), include_groups=False),
})
res["pct"] = (100 * res["n_myeloid_like"] / res["n_cells"]).round(1)
print(res.to_string())

print("\nWhere do myeloid-like cells sit by leiden (top clusters per sample)?")
my = df[df["myeloid_score"] > 0.1]
print(my.groupby(["sample", "leiden"]).size().unstack(fill_value=0).to_string())
