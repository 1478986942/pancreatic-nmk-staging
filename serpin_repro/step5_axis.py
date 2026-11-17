"""
Step 5: Serpin-myeloid axis quantification.

CORRECTED 2026-11-14 (see CORRECTIONS.md C4).  This script previously described
its classic-vs-Flex comparison as "effect of anti-PD-1 (Flex) vs vehicle
(classic)".  GSE289138 has no treatment arms: the six libraries are
control KO / Serpinb2 KO / Serpine1 KO, each sequenced with two chemistries
(3' classic and Flex).  `condition` is therefore a LIBRARY CHEMISTRY, and any
comparison across it is a platform comparison, not a treatment comparison.

Do not report a treatment effect from this script.
"""

import os

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE, "processed")
FIG_DIR = os.path.join(BASE, "figures")


def main():
    adata = sc.read_h5ad(os.path.join(OUT_DIR, "GSE289138_annotated.h5ad"))
    my = sc.read_h5ad(os.path.join(OUT_DIR, "myeloid.h5ad"))

    counts = adata.layers["counts"]
    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    for gene in ["Serpine1", "Serpinb2", "Spp1", "Marco"]:
        col = np.asarray(counts[:, var_idx[gene]].todense()).ravel()
        adata.obs[f"{gene}_detected"] = col >= 2

    # ------------------------------------------------ 1. malignant epithelial Serpin
    print("=" * 70)
    print("1. SERPINE1/SERPINB2 IN MALIGNANT EPITHELIAL CELLS")
    print("=" * 70)
    epi = adata.obs["anno"] == "Epithelial"
    mal = adata.obs["malignant"] == "malignant"
    for label, mask in [("malignant", epi & mal), ("non-malignant epi", epi & ~mal)]:
        n = int(mask.sum())
        s1 = adata.obs["Serpine1_detected"][mask].mean() * 100
        s2 = adata.obs["Serpinb2_detected"][mask].mean() * 100
        print(f"  {label:20s} n={n:6d}  Serpine1+={s1:5.1f}%  Serpinb2+={s2:5.1f}%")

    # Serpine1-high malignant epithelial subpopulation
    print("\n  --- SERPINE1-high malignant epithelial (>=2 UMI) per condition ---")
    rows = []
    for cond in ["classic", "Flex"]:
        m = epi & mal & (adata.obs["condition"] == cond)
        n = int(m.sum())
        rows.append({
            "condition": cond,
            "n_malignant_epi": n,
            "Serpine1_pct": round(adata.obs["Serpine1_detected"][m].mean() * 100, 1),
            "Serpinb2_pct": round(adata.obs["Serpinb2_detected"][m].mean() * 100, 1),
        })
    cond_df = pd.DataFrame(rows)
    print(cond_df.to_string(index=False))

    # ------------------------------------------------ 2. SERPINB2+ macrophage niche
    print("\n" + "=" * 70)
    print("2. SERPINB2+ / MARCO+ MACROPHAGE NICHE (myeloid cluster 20)")
    print("=" * 70)
    mycounts = my.layers["counts"]
    mvar = {g: i for i, g in enumerate(my.var_names)}
    my.obs["Serpinb2_detected"] = np.asarray(mycounts[:, mvar["Serpinb2"]].todense()).ravel() >= 2
    my.obs["Marco_detected"] = np.asarray(mycounts[:, mvar["Marco"]].todense()).ravel() >= 2
    my.obs["Spp1_detected"] = np.asarray(mycounts[:, mvar["Spp1"]].todense()).ravel() >= 2
    my.obs["SPP1_MARCO_pp"] = my.obs["Spp1_detected"] & my.obs["Marco_detected"]

    niche = my.obs["leiden_myeloid"] == "20"
    print(f"  myeloid cluster 20: n={int(niche.sum())}")
    print(f"    Serpinb2+ : {my.obs['Serpinb2_detected'][niche].mean() * 100:.1f}%")
    print(f"    Marco+   : {my.obs['Marco_detected'][niche].mean() * 100:.1f}%")
    print(f"    SPP1+MARCO+ : {my.obs['SPP1_MARCO_pp'][niche].mean() * 100:.1f}%")
    others = ~niche
    print(f"  other myeloid clusters: Serpinb2+ {my.obs['Serpinb2_detected'][others].mean() * 100:.1f}%, "
          f"Marco+ {my.obs['Marco_detected'][others].mean() * 100:.1f}%")

    # condition comparison within niche
    print("\n  --- niche occupancy by condition ---")
    rows = []
    for cond in ["classic", "Flex"]:
        m = my.obs["condition"] == cond
        total = int(m.sum())
        rows.append({
            "condition": cond,
            "n_myeloid": total,
            "n_niche": int((niche & m).sum()),
            "niche_pct": round((niche & m).sum() / total * 100, 2),
            "spp1_marco_pp_n": int((my.obs["SPP1_MARCO_pp"] & m).sum()),
            "serpinb2_pos_n": int((my.obs["Serpinb2_detected"] & m).sum()),
        })
    niche_df = pd.DataFrame(rows)
    print(niche_df.to_string(index=False))

    # statistical test on niche fraction
    tab = pd.crosstab(my.obs["condition"][niche.values], my.obs["condition"])
    cont = np.array([
        [int((niche & (my.obs["condition"] == c)).sum()),
         int(((~niche) & (my.obs["condition"] == c)).sum())]
        for c in ["classic", "Flex"]
    ])
    chi2, p, dof, expected = stats.chi2_contingency(cont)
    print(f"\n  chi2 test niche fraction classic vs Flex: chi2={chi2:.2f} p={p:.4f}")

    # per-sample niche fraction
    print("\n  --- per sample ---")
    ps = (
        my.obs.assign(niche=niche.values)
        .groupby(["filename", "condition"], observed=True)["niche"]
        .agg(["sum", "size"])
    )
    ps["niche_pct"] = (ps["sum"] / ps["size"] * 100).round(2)
    print(ps.to_string())
    ps.to_csv(os.path.join(OUT_DIR, "niche_fraction_by_sample.tsv"), sep="\t")

    # ------------------------------------------------ 3. epithelial Serpin state vs niche
    print("\n" + "=" * 70)
    print("3. TUMOUR-INTRINSIC SERPIN STATE vs IMMUNE EXCLUSION")
    print("=" * 70)
    # T/NK abundance as immune-exclusion readout
    tnk = (adata.obs["anno"] == "T/NK").to_numpy()
    cond_tnk = pd.crosstab(adata.obs["condition"], tnk)
    print("  T/NK cell fraction by condition:")
    for c in ["classic", "Flex"]:
        n = int((adata.obs["condition"] == c).sum())
        print(f"    {c}: {tnk[(adata.obs['condition'] == c).to_numpy()].sum()} / {n} "
              f"({tnk[(adata.obs['condition'] == c).to_numpy()].mean() * 100:.1f}%)")

    # Serpine1-high malignant epithelial fraction per sample
    print("\n  Serpine1-high malignant fraction per sample:")
    s1hi = epi & mal & adata.obs["Serpine1_detected"]
    ss = (
        adata.obs.assign(s1hi=s1hi)
        .groupby(["filename", "condition"], observed=True)["s1hi"]
        .agg(["sum", "size"])
    )
    ss["pct"] = (ss["sum"] / ss["size"] * 100).round(2)
    print(ss.to_string())
    ss.to_csv(os.path.join(OUT_DIR, "serpine1_malignant_fraction_by_sample.tsv"), sep="\t")

    # ------------------------------------------------ 4. figures
    sc.settings.figdir = FIG_DIR
    os.makedirs(FIG_DIR, exist_ok=True)
    sc.pl.umap(
        adata,
        color=["anno", "malignant", "Serpine1_detected", "Serpinb2_detected", "condition"],
        save="_umap_final_overview.png",
        ncols=3,
    )
    sc.pl.umap(
        my,
        color=["leiden_myeloid", "Serpinb2_detected", "Marco_detected", "SPP1_MARCO_pp", "condition"],
        save="_umap_myeloid_final.png",
        ncols=3,
    )
    print("\nfigures written to", FIG_DIR)


if __name__ == "__main__":
    main()
