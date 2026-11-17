"""Recompute the MARCO axis in GSE202051 with the CORRECT denominator.

Why this exists
---------------
step3_stats.py built its myeloid mask as

    myeloid_cols = [c for c in obs.columns
                    if any(x in c.lower() for x in
                           ["monocyte","macrophage","m0","m1","m2"])
                    and "score" not in c.lower()]
    myeloid_mask = obs[myeloid_cols].any(axis=1)

Those matched columns (Macrophage, Monocytes_1/2, M0/M1/M2) are continuous
float scores, not booleans.  `.any(axis=1)` treats any non-zero float as True,
so myeloid_mask covered ALL 224,988 nuclei and

    MARCO_frac = n_MARCO_pos / n_myeloid   ==   MARCO+ / ALL nuclei

The label "MARCO+ macrophage fraction" was therefore never computed.  Here the
denominator is the dataset's own annotation (Level 1 == 'Myeloid' and
Level 3 == 'Macrophage'), and MARCO positivity is taken from the raw counts
column streamed out of the CSR matrix.
"""
import os
import sys
import numpy as np
import pandas as pd
import h5py
import anndata as ad
from scipy import stats

ADATA = r"E:\pancreatic_cancer\submit\gse202051_hum\data\GSE202051_total.h5ad"
OUT = r"E:\pancreatic_cancer\submit\gse202051_hum\processed"


def decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def stream_gene_column(path, gene):
    """Return (n_cells,) dense vector of raw counts for one gene, block-wise."""
    f = h5py.File(path, "r")
    var_names = decode(f["var"]["_index"][:])
    if gene not in var_names:
        raise KeyError(f"{gene} not in var_names")
    j = int(np.where(var_names == gene)[0][0])
    print(f"  {gene} is column {j} of {len(var_names)}", flush=True)

    indptr = f["X"]["indptr"][:]
    n = len(indptr) - 1
    out = np.zeros(n, dtype=np.float32)

    data_ds = f["X"]["data"]
    idx_ds = f["X"]["indices"]

    BLOCK = 20000
    for start in range(0, n, BLOCK):
        stop = min(start + BLOCK, n)
        a, b = int(indptr[start]), int(indptr[stop])
        if b <= a:
            continue
        idx = idx_ds[a:b]
        dat = data_ds[a:b]
        # row id for every nonzero in this block
        rows = np.repeat(np.arange(start, stop), np.diff(indptr[start:stop + 1]))
        m = idx == j
        if m.any():
            out[rows[m]] = dat[m]
        print(f"    rows {start}-{stop}  ({b-a} nnz)", end="\r", flush=True)
    f.close()
    print(f"\n  {gene}: detected in {int((out > 0).sum())} / {n} nuclei", flush=True)
    return out


def main():
    print("[1] reading obs annotation", flush=True)
    a = ad.read_h5ad(ADATA, backed="r")
    obs = a.obs[["pid", "treatment_status", "Level 1 Annotation",
                 "Level 3 Annotation"]].copy()
    a.file.close()
    print(f"  {len(obs)} nuclei, {obs['pid'].nunique()} patients", flush=True)
    print("  treatment_status:", obs["treatment_status"].value_counts().to_dict(),
          flush=True)

    print("\n[2] streaming MARCO counts", flush=True)
    marco = stream_gene_column(ADATA, "MARCO")

    obs["MARCO_expr"] = marco
    obs["MARCO_pos"] = marco > 0

    lvl1 = obs["Level 1 Annotation"].astype(str)
    lvl3 = obs["Level 3 Annotation"].astype(str)
    obs["is_myeloid"] = lvl1 == "Myeloid"
    obs["is_macrophage"] = lvl3 == "Macrophage"
    print(f"\n  Myeloid (Level 1): {int(obs['is_myeloid'].sum())}", flush=True)
    print(f"  Macrophage (Level 3): {int(obs['is_macrophage'].sum())}", flush=True)

    status_map = {"Untreated": "naive", "Treated": "treated"}
    obs["status"] = obs["treatment_status"].map(status_map)
    if obs["status"].isna().any():
        print("  treatment_status values:", obs["treatment_status"].unique(), flush=True)
        obs["status"] = obs["treatment_status"].astype(str).str.lower().map(
            lambda s: "naive" if "untreat" in s or "naive" in s else "treated")

    print("\n[3] per-patient aggregation", flush=True)
    rows = []
    for pid, grp in obs.groupby("pid"):
        mac = grp[grp["is_macrophage"]]
        mye = grp[grp["is_myeloid"]]
        rows.append(dict(
            pid=pid, status=grp["status"].iloc[0],
            n_cells=len(grp),
            n_macrophage=len(mac),
            n_myeloid=len(mye),
            marco_pos_mac=int(mac["MARCO_pos"].sum()),
            marco_pos_mye=int(mye["MARCO_pos"].sum()),
            MARCO_frac_of_macrophage=(mac["MARCO_pos"].mean()
                                      if len(mac) else np.nan),
            MARCO_frac_of_myeloid=(mye["MARCO_pos"].mean()
                                   if len(mye) else np.nan),
            MARCO_per_100_nuclei=100.0 * grp["MARCO_pos"].mean(),
        ))
    pat = pd.DataFrame(rows).sort_values(["status", "pid"])
    pat.to_csv(os.path.join(OUT, "patient_marco_corrected.tsv"), sep="\t",
               index=False)

    pd.set_option("display.width", 220)
    print(pat.round(4).to_string(index=False), flush=True)

    print("\n[4] naive vs treated (patient level, Mann-Whitney)", flush=True)
    res = []
    for metric in ["MARCO_frac_of_macrophage", "MARCO_frac_of_myeloid",
                   "MARCO_per_100_nuclei"]:
        nv = pat.loc[pat["status"] == "naive", metric].dropna().values
        tv = pat.loc[pat["status"] == "treated", metric].dropna().values
        if len(nv) < 2 or len(tv) < 2:
            continue
        u, p = stats.mannwhitneyu(nv, tv, alternative="two-sided")
        gt = sum((x > y) for x in nv for y in tv)
        lt = sum((x < y) for x in nv for y in tv)
        delta = (gt - lt) / (len(nv) * len(tv))
        res.append(dict(metric=metric, n_naive=len(nv), n_treated=len(tv),
                        naive_mean=round(float(np.mean(nv)), 4),
                        treated_mean=round(float(np.mean(tv)), 4),
                        naive_median=round(float(np.median(nv)), 4),
                        treated_median=round(float(np.median(tv)), 4),
                        cliffs_delta=round(float(delta), 3),
                        p=round(float(p), 4)))
        print(f"  {metric:28s} naive={np.mean(nv):.4f} treated={np.mean(tv):.4f} "
              f"delta={delta:+.3f} p={p:.4f}", flush=True)

    df = pd.DataFrame(res)
    df.to_csv(os.path.join(OUT, "marco_corrected_test.tsv"), sep="\t", index=False)

    print("\n[5] macrophage abundance also re-checked", flush=True)
    pat["macrophage_frac_of_nuclei"] = pat["n_macrophage"] / pat["n_cells"]
    nv = pat.loc[pat["status"] == "naive", "macrophage_frac_of_nuclei"]
    tv = pat.loc[pat["status"] == "treated", "macrophage_frac_of_nuclei"]
    u, p = stats.mannwhitneyu(nv.dropna(), tv.dropna(), alternative="two-sided")
    print(f"  macrophage fraction of nuclei: naive={nv.mean():.4f} "
          f"treated={tv.mean():.4f} p={p:.4f}", flush=True)

    print("\n[saved] patient_marco_corrected.tsv, marco_corrected_test.tsv",
          flush=True)


if __name__ == "__main__":
    main()
