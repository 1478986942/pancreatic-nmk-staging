"""Reproduction audit of GSE289138 (Falcomata et al., Nature 2026).

GSE289138 is the source study's own mouse scRNA-seq deposit.  Its design is
3 genotypes x 2 library chemistries (3' classic vs Flex), **one tumour per
genotype**; the two chemistries are technical replicates of the same tumour.
There is no treatment arm.  See CORRECTIONS.md C4.

This script does two things:

1. Re-asserts the corrected design on the genotype re-run
   (`GSE289138_processed.h5ad` must carry `genotype`, `condition` = chemistry,
   and must not carry `treatment`).
2. Reproduces, descriptively, the readouts the source paper derives from this
   deposit, and quantifies how much of each readout is genotype versus
   chemistry.

It deliberately reports **no cell-level p-values**.  With one tumour per
genotype, a cell-level test is pseudo-replicated; the unit of inference would
have to be the tumour, and there is only one per genotype.

Run:
    python audit_mouse_confound.py
"""

import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "serpin_repro", "processed")
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "figs_reproduction")
PROCESSED = os.path.join(PROC, "GSE289138_processed.h5ad")
ANNOTATED = os.path.join(PROC, "GSE289138_annotated.h5ad")
META = os.path.join(PROC, "sample_metadata_CORRECTED.tsv")

GENOTYPES = ["control_KO", "Serpinb2_KO", "Serpine1_KO"]
CHEMISTRIES = ["classic", "Flex"]

# Signature genes.  These are readouts, not findings: they approximate the
# subsets the paper names, and are reported as mean log-normalised expression
# per library, never as a cell-level test.
TNK_EXHAUST = ["Pdcd1", "Havcr2", "Tox", "Lag3", "Tigit", "Ctla4", "Entpd1"]
TNK_EFF = ["Gzmb", "Prf1", "Ifng", "Nkg7", "Ccl5", "Gzmk"]
TNK_NAIVE = ["Sell", "Il7r", "Tcf7", "Slamf6", "Lef1"]
MAC_SUPPRESS = ["Spp1", "Marco", "Arg1", "C1qa", "C1qb", "C1qc"]
MAC_CXCL9 = ["Cxcl9", "Cxcl10", "Cxcl11"]


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------
# 1. Metadata assertions
# --------------------------------------------------------------------------
def load_metadata():
    m = pd.read_csv(META, sep="\t")
    assert m["filename"].is_unique, "Fatal: duplicated library in metadata"
    assert set(m["genotype"]) == set(GENOTYPES), \
        f"Fatal: unexpected genotypes {set(m['genotype'])}"
    assert set(m["chemistry"]) == set(CHEMISTRIES), \
        f"Fatal: unexpected chemistries {set(m['chemistry'])}"
    # One tumour per genotype: this is the fact that forbids cell-level tests.
    per_geno = m.groupby("genotype")["tumor_id"].nunique()
    assert (per_geno == 1).all(), f"Fatal: >1 tumour per genotype: {per_geno.to_dict()}"
    return m.set_index("filename")


def load_processed(meta):
    a = sc.read_h5ad(PROCESSED)
    assert "treatment" not in a.obs.columns, \
        "Fatal: fabricated treatment column survived the genotype re-run"
    for col in ["genotype", "condition", "sample"]:
        assert col in a.obs.columns, f"Fatal: {col} missing from processed object"

    # Re-derive every design label from the authoritative table, not from the
    # object, and fail if they disagree.
    src = {"sample": "tumor_id", "condition": "chemistry", "genotype": "genotype"}
    for col, mcol in src.items():
        want = meta[mcol].reindex(a.obs["filename"].astype(str)).to_numpy()
        got = a.obs[col].astype(str).to_numpy()
        assert (got == want).all(), f"Fatal: {col} disagrees with metadata table"

    assert set(a.obs["genotype"].astype(str)) == set(GENOTYPES)
    return a


def attach_annotation(a):
    """Transfer the leiden->anno map and CNV calls from the annotated object.

    The cell set is identical (both objects derive from the same per-sample
    cache), so the transfer is a join on barcode, not a re-clustering.
    """
    if not os.path.exists(ANNOTATED):
        raise SystemExit(
            "annotated object missing; run serpin_repro/step3_annotate_cnv.py first")
    b = ad.read_h5ad(ANNOTATED, backed="r")
    missing = set(a.obs_names) - set(b.obs_names)
    assert not missing, f"Fatal: {len(missing)} processed cells absent from annotated"
    for col in ["anno", "malignant", "leiden_epi"]:
        if col in b.obs.columns:
            a.obs[col] = b.obs[col].reindex(a.obs_names).astype(str).to_numpy()
    b.file.close()
    assert "anno" in a.obs.columns, "Fatal: annotation transfer failed"
    return a


# --------------------------------------------------------------------------
# 2. Readouts
# --------------------------------------------------------------------------
def _dense(x):
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def signature(a, genes, mask):
    present = [g for g in genes if g in a.var_names]
    if not present:
        return np.full(int(mask.sum()), np.nan)
    idx = [a.var_names.get_loc(g) for g in present]
    sub = _dense(a.X[mask][:, idx]).astype(float)
    return sub.mean(axis=1)


def detection(a, gene, mask):
    if gene not in a.var_names:
        return np.full(int(mask.sum()), np.nan)
    col = _dense(a.layers["counts"][mask][:, a.var_names.get_loc(gene)]).ravel()
    return (col >= 2).astype(float)


def per_library(a, meta):
    rows = []
    anno = a.obs["anno"].astype(str).to_numpy()
    mal = a.obs["malignant"].astype(str).to_numpy()
    is_tnk = anno == "T/NK"
    is_mac = anno == "Macrophage"
    is_caf = anno == "CAF"
    is_epi = anno == "Epithelial"
    is_mal = is_epi & (mal == "malignant")

    for lib, g in a.obs.groupby("filename", observed=True):
        idx = g.index
        mask_lib = a.obs["filename"].astype(str).to_numpy() == str(lib)
        n = int(mask_lib.sum())
        rec = {
            "filename": str(lib),
            "tumour": meta.loc[str(lib), "tumor_id"],
            "genotype": meta.loc[str(lib), "genotype"],
            "chemistry": meta.loc[str(lib), "chemistry"],
            "n_cells": n,
        }

        # composition (% of all cells)
        for name, mask in [("TNK", is_tnk), ("Macrophage", is_mac),
                           ("CAF", is_caf), ("Epithelial", is_epi),
                           ("Neutrophil", anno == "Neutrophil")]:
            rec[f"pct_{name}"] = round(float((mask & mask_lib).sum()) / n * 100, 2)

        # T/NK state signatures (mean log-normalised expression among T/NK)
        tnk_lib = is_tnk & mask_lib
        for name, genes in [("exhaust", TNK_EXHAUST), ("effector", TNK_EFF),
                            ("naive", TNK_NAIVE)]:
            v = signature(a, genes, tnk_lib)
            rec[f"tnk_{name}"] = round(float(np.nanmean(v)), 3) if len(v) else np.nan

        # Macrophage state (% of macrophages, count-based)
        mac_lib = is_mac & mask_lib
        rec["mac_spp1_marco_pct"] = round(
            float(np.nanmean(detection(a, "Spp1", mac_lib)
                             * detection(a, "Marco", mac_lib))) * 100, 2
        ) if mac_lib.sum() else np.nan
        rec["mac_cxcl9_pct"] = round(
            float(np.nanmean(detection(a, "Cxcl9", mac_lib))) * 100, 2
        ) if mac_lib.sum() else np.nan

        # Malignant epithelial SERPIN detection (% of malignant epithelial)
        mal_lib = is_mal & mask_lib
        for gene, key in [("Serpine1", "mal_serpine1_pct"),
                          ("Serpinb2", "mal_serpinb2_pct")]:
            rec[key] = round(float(np.nanmean(detection(a, gene, mal_lib))) * 100, 2) \
                if mal_lib.sum() else np.nan

        rows.append(rec)

    tab = pd.DataFrame(rows).sort_values(["chemistry", "genotype"]).reset_index(drop=True)
    return tab


def effects(tab, col):
    """Descriptive genotype and chemistry effect sizes, in the column's units.

    genotype effect = mean over chemistries of (max - min across genotypes)
    chemistry effect = mean over genotypes of |classic - Flex| for one tumour
    """
    geno = []
    for chem in CHEMISTRIES:
        v = tab.loc[tab["chemistry"] == chem, col].dropna()
        geno.append(float(v.max() - v.min()) if len(v) > 1 else np.nan)
    chem = []
    for gt in GENOTYPES:
        v = tab.loc[tab["genotype"] == gt].set_index("chemistry")[col]
        if set(CHEMISTRIES) <= set(v.index):
            chem.append(abs(float(v["classic"]) - float(v["Flex"])))
    return float(np.nanmean(geno)), float(np.nanmean(chem)) if chem else np.nan


READOUTS = [
    "pct_TNK", "pct_Macrophage", "pct_CAF",
    "tnk_exhaust", "tnk_effector", "tnk_naive",
    "mac_spp1_marco_pct", "mac_cxcl9_pct",
    "mal_serpine1_pct", "mal_serpinb2_pct",
]


# --------------------------------------------------------------------------
# 3. Subset-level readouts -- analogues of Fig. 2f (T-cell subsets),
#    Fig. 4g (macrophage subsets), Fig. 3g/h + ED Fig. 7c (CAF subsets).
#    Same discipline as above: per-library descriptive fractions, and no
#    cell-level p-values anywhere.
# --------------------------------------------------------------------------

def _module_means_by_cluster(sub, key, modules):
    """Mean expression of each gene module, per level of `key`."""
    cl_vec = sub.obs[key].astype(str).to_numpy()
    order = sorted(set(cl_vec), key=lambda x: int(x) if x.isdigit() else 10**6)
    out = {}
    for cl in order:
        m = cl_vec == cl
        row = {}
        for name, genes in modules.items():
            vals = []
            for g in genes:
                if g in sub.var_names:
                    x = sub.X[m][:, sub.var_names.get_loc(g)]
                    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
                    vals.append(float(np.asarray(x).ravel().mean()))
            row[name] = float(np.mean(vals)) if vals else np.nan
        out[cl] = row
    return out


def _argmax_labels(mm, modules):
    """Label every cluster by its most z-prominent module."""
    names = list(modules.keys())
    mat = np.array([[mm[cl][n] for n in names] for cl in mm], dtype=float)
    mu = np.nanmean(mat, axis=0)
    sd = np.nanstd(mat, axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    z = (mat - mu) / sd
    return {cl: names[int(np.nanargmax(z[i]))] for i, cl in enumerate(mm)}


def _frac_table(sub, label_col, meta):
    lib = sub.obs["filename"].astype(str).to_numpy()
    lab = sub.obs[label_col].astype(str).to_numpy()
    rows = []
    for lib_id in meta.index:
        m = lib == lib_id
        n = int(m.sum())
        if n == 0:
            continue
        row = {"filename": lib_id, "genotype": meta.loc[lib_id, "genotype"],
               "chemistry": meta.loc[lib_id, "chemistry"], "n": n}
        for s in sorted(set(lab)):
            row[s] = round(float((lab[m] == s).sum()) / n * 100, 2)
        rows.append(row)
    return pd.DataFrame(rows)


def _cell_module_scores(sub, modules):
    """Per-cell mean log-normalised expression of each module."""
    out = {}
    for name, genes in modules.items():
        gs = [g for g in genes if g in sub.var_names]
        if not gs:
            out[name] = np.full(sub.n_obs, np.nan)
            continue
        idx = [sub.var_names.get_loc(g) for g in gs]
        X = sub.X[:, idx]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        out[name] = X.astype(float).mean(axis=1)
    return pd.DataFrame(out, index=sub.obs_names)


def tnk_subsets(a, meta):
    """T-cell subsets of Fig. 2f, by per-cell module-score classification.

    Restricted to the classic (fresh 5' GEX) libraries: that is the platform
    the source paper used for the T-cell analysis, and pooling the
    chemistries puts chemistry-specific clusters into the graph.

    Classification rules (stated up front, deterministic):
      CD8 (Cd8a >= Cd4):
        Tex      if exhaust-module > effector-module and > naive-module;
                 Tex cells split into prog / int / term by tertiles of the
                 TCF7 module (Tcf7, Sell, Il7r): high = prog, low = term.
        Tn-like  else if naive-module > effector-module.
        Teff     otherwise.
      CD4 (Cd4 > Cd8a):
        Treg     if Foxp3 > 0.
        Tn-like  else if naive-module > effector-module.
        Teff Gzmk+ otherwise (the paper's only CD4 effector subset).
      Neither marker detected: NK/other.
    """
    sub = a[(a.obs["anno"].astype(str) == "T/NK")
            & (a.obs["condition"].astype(str) == "classic")].copy()

    def gene_vec(g):
        if g not in sub.var_names:
            return np.zeros(sub.n_obs)
        x = sub.X[:, sub.var_names.get_loc(g)]
        x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
        return np.asarray(x).ravel().astype(float)

    cd8a, cd4, nkg7 = gene_vec("Cd8a"), gene_vec("Cd4"), gene_vec("Nkg7")
    mods = _cell_module_scores(sub, {
        "naive": ["Sell", "Tcf7", "Ccr7", "Il7r"],
        "eff": ["Gzmk", "Gzmb", "Prf1", "Nkg7"],
        "exh": ["Pdcd1", "Havcr2", "Tox", "Lag3", "Entpd1"],
        "tcf7m": ["Tcf7", "Sell", "Il7r"],
    })

    lab = np.empty(sub.n_obs, dtype=object)
    # lineage gates require the marker to be detected at all; cells with
    # neither Cd8a nor Cd4 are NK/other (the paper's Fig. 2f panel excludes
    # NK cells, so this bucket is reported but not part of the 8 subsets)
    is_cd8 = (cd8a > cd4) & (cd8a > 0)
    is_cd4 = (cd4 > cd8a) & (cd4 > 0)
    foxp3 = gene_vec("Foxp3")

    cd8_idx = np.where(is_cd8)[0]
    exh_dom = (mods.loc[:, "exh"].to_numpy() > mods.loc[:, "eff"].to_numpy()) & \
              (mods.loc[:, "exh"].to_numpy() > mods.loc[:, "naive"].to_numpy())
    tex_idx = cd8_idx[exh_dom[cd8_idx]]
    lab[:] = "NK/other"
    lab[is_cd8 & ~exh_dom
        & (mods.loc[:, "naive"].to_numpy() <= mods.loc[:, "eff"].to_numpy())] = "CD8 Teff"
    lab[is_cd8 & ~exh_dom
        & (mods.loc[:, "naive"].to_numpy() > mods.loc[:, "eff"].to_numpy())] = "CD8 Tn-like"
    lab[is_cd4 & (foxp3 > 0)] = "CD4 Treg"
    rest_cd4 = is_cd4 & (foxp3 <= 0) & (lab != "CD4 Treg")
    lab[rest_cd4 & (mods.loc[:, "naive"].to_numpy() >= mods.loc[:, "eff"].to_numpy())] = "CD4 Tn-like"
    lab[rest_cd4 & (mods.loc[:, "naive"].to_numpy() < mods.loc[:, "eff"].to_numpy())] = "CD4 Teff Gzmk+"
    lab[~is_cd8 & ~is_cd4] = "NK/other"

    # Tex tertiles by the TCF7 module: rank-based positional split (robust
    # to the zero-inflated ties that break value-percentile tertiles)
    if len(tex_idx) >= 3:
        t = mods.loc[:, "tcf7m"].to_numpy()[tex_idx]
        order = np.argsort(t, kind="stable")  # ascending TCF7 module
        n = len(tex_idx)
        n_term = n // 3
        n_prog = n // 3
        for k, pos in enumerate(order):
            if k < n_term:
                lab[tex_idx[pos]] = "CD8 Tex term"
            elif k >= n - n_prog:
                lab[tex_idx[pos]] = "CD8 Tex prog"
            else:
                lab[tex_idx[pos]] = "CD8 Tex int"
    elif len(tex_idx) in (1, 2):
        for i in tex_idx:
            lab[i] = "CD8 Tex term"

    sub.obs["subset"] = pd.Categorical(lab)
    log("T/NK per-cell classification counts:")
    log(sub.obs["subset"].value_counts().to_string())
    _dump_subset_umap(sub, "tcell_umap.tsv")
    return _frac_table(sub, "subset", meta), {}


def macro_subsets(meta):
    """Myeloid subclustering (step4 leiden); approximate Fig. 4g states."""
    import scanpy as sc
    my = sc.read_h5ad(os.path.join(PROC, "myeloid.h5ad"))
    modules = {
        "Spp1/Marco+": ["Spp1", "Marco"],
        "Cxcl9+": ["Cxcl9", "Cxcl10"],
        "Proliferating": ["Mki67", "Top2a", "Stmn1"],
        "Arg1+": ["Arg1"],
        "Il1b+": ["Il1b"],
        "Folr2+": ["Folr2", "Trem2"],
        "Monocyte-like": ["Ly6c2", "Ccr2", "Plac8"],
        "C1q+ homeostatic": ["C1qa", "C1qb", "C1qc"],
        "DC-like": ["H2-Ab1", "H2-Aa", "Flt3"],
    }
    mm = _module_means_by_cluster(my, "leiden_myeloid", modules)
    labels = _argmax_labels(mm, modules)
    my.obs["subset"] = pd.Categorical(
        [labels[c] for c in my.obs["leiden_myeloid"].astype(str)])
    log("Myeloid cluster -> subset map: " + str(labels))
    _dump_subset_umap(my, "myeloid_umap.tsv")
    return _frac_table(my, "subset", meta), labels


def caf_subsets(a, meta):
    """CAF subclustering; approximate Fig. 3g/h + ED Fig. 7c (Flex stroma).

    Restricted to the Flex libraries: the source paper subclustered CAFs from
    the Flex protocol specifically (ED Fig. 7c), and pooled-chemistry CAF
    clustering produced chemistry-pure clusters with unstable labels.
    """
    import scanpy as sc
    sub = a[(a.obs["anno"].astype(str) == "CAF")
            & (a.obs["condition"].astype(str) == "Flex")].copy()
    sc.tl.pca(sub, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=15)
    sc.tl.leiden(sub, resolution=1.0, key_added="leiden_caf")
    modules = {
        "iCAF": ["Il6", "Cxcl1", "Cxcl2", "Pdgfra"],
        "MyCAF": ["Acta2", "Tagln", "Myh11"],
        "apCAF": ["H2-Ab1", "H2-Aa", "Cd74"],
        "LRCC15+ MyCAF": ["Lrcc15", "Htra1"],
        "Matrix CAF": ["Col1a1", "Col1a2", "Dcn", "Lum"],
    }
    mm = _module_means_by_cluster(sub, "leiden_caf", modules)
    labels = _argmax_labels(mm, modules)
    sub.obs["subset"] = pd.Categorical(
        [labels[c] for c in sub.obs["leiden_caf"].astype(str)])
    log("CAF cluster -> subset map: " + str(labels))
    _dump_subset_umap(sub, "caf_umap.tsv")
    return _frac_table(sub, "subset", meta), labels


def _dump_subset_umap(sub, fname, label_col="subset"):
    """Write per-cell UMAP coordinates + compartment label for one subset.

    The R figure script (methods_paper/figs_nature.R) needs these to redraw the
    compartment UMAP panels (Fig. 3g, ED Fig. 5g, 7c/f, 8p) from exactly the
    cells that produced the proportion tables. Descriptive, no statistics.
    """
    out = os.path.join(FIG_DIR, "plot_tables")
    os.makedirs(out, exist_ok=True)
    if "X_umap" not in sub.obsm:
        log("  no X_umap on %s -- %s not written" % (label_col, fname))
        return
    xy = np.asarray(sub.obsm["X_umap"], dtype=float)
    cols = [c for c in (label_col, "genotype", "condition")
            if c in sub.obs.columns]
    df = sub.obs[cols].copy()
    df["UMAP1"] = xy[:, 0]
    df["UMAP2"] = xy[:, 1]
    df.to_csv(os.path.join(out, fname), sep="\t", index=False)
    log("  wrote %s (%d cells)" % (fname, len(df)))


def endo_subsets(a, meta):
    """Endothelial subclustering; approximate ED Fig. 7f,g (Flex stroma).

    Same descriptive scheme as caf_subsets: Leiden on the Flex endothelial
    compartment, each cluster named by the argmax of three canonical vascular
    modules. Returns (None, {}) when the compartment is too small.
    """
    import scanpy as sc
    sub = a[(a.obs["anno"].astype(str) == "Endothelial")
            & (a.obs["condition"].astype(str) == "Flex")].copy()
    if sub.n_obs < 60:
        log("  endothelial Flex compartment too small (%d) -- skipped"
            % sub.n_obs)
        return None, {}
    sc.tl.pca(sub, n_comps=min(30, sub.n_obs - 1), svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=15)
    sc.tl.leiden(sub, resolution=1.0, key_added="leiden_endo")
    modules = {
        "Capillary": ["Gpihbp1", "Car4", "Cd36", "Fabp4"],
        "Venous": ["Ackr1", "Nr2f2", "Vwf", "Sele"],
        "Lymphatic": ["Prox1", "Lyve1", "Pdpn", "Ccl21a"],
    }
    mm = _module_means_by_cluster(sub, "leiden_endo", modules)
    labels = _argmax_labels(mm, modules)
    sub.obs["subset"] = pd.Categorical(
        [labels[c] for c in sub.obs["leiden_endo"].astype(str)])
    log("Endothelial cluster -> subset map: " + str(labels))
    _dump_subset_umap(sub, "endo_umap.tsv")
    frac = _frac_table(sub, "subset", meta)
    frac.to_csv(os.path.join(FIG_DIR, "plot_tables", "endo_subsets.tsv"),
                sep="\t", index=False)
    return frac, labels


def export_plot_tables(a, meta, tab, tnk_renorm, mac_frac, caf_frac, eff):
    """Export the tables behind the six scRNA figures so they can be redrawn
    in R (methods_paper/figs_nature.R). No p-values, per-library only."""
    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, "plot_tables")
    os.makedirs(out, exist_ok=True)

    tnk_renorm.to_csv(os.path.join(out, "tcell_subsets.tsv"),
                      sep="\t", index=False)
    mac_frac.to_csv(os.path.join(out, "myeloid_subsets.tsv"),
                    sep="\t", index=False)
    caf_frac.to_csv(os.path.join(out, "caf_subsets.tsv"),
                    sep="\t", index=False)
    tab[["filename", "genotype", "chemistry",
         "mal_serpine1_pct", "mal_serpinb2_pct"]].to_csv(
        os.path.join(out, "serpin_detection.tsv"), sep="\t", index=False)
    eff.to_csv(os.path.join(out, "effect_sizes.tsv"), sep="\t", index=False)

    # UMAP: downsample so the file stays small enough for ggplot2
    import pandas as pd
    xy = np.asarray(a.obsm["X_umap"], dtype=float)
    obs = a.obs[["anno", "genotype", "condition"]].copy()
    obs["UMAP1"] = xy[:, 0]
    obs["UMAP2"] = xy[:, 1]
    rng = np.random.default_rng(0)
    keep = rng.choice(obs.shape[0], size=min(25000, obs.shape[0]),
                      replace=False)
    obs.iloc[np.sort(keep)].to_csv(
        os.path.join(out, "umap_coords.tsv"), sep="\t", index=False)
    log("plot tables written to " + out)


def make_figures(a, meta, tab, tnk_renorm, mac_frac, caf_frac, eff):
    """Six descriptive figures for reproduction_report.tex. No p-values."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    os.makedirs(FIG_DIR, exist_ok=True)
    geno_order = ["control_KO", "Serpinb2_KO", "Serpine1_KO"]
    geno_colors = {"control_KO": "#4d4d4d", "Serpinb2_KO": "#1f77b4",
                   "Serpine1_KO": "#d62728"}
    w = 0.26

    # --- Fig A: UMAP overview (ED Fig. 5f analogue) ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, col, title in zip(
            axes, ["anno", "genotype", "condition"],
            ["cell type (annotation)", "genotype", "library chemistry"]):
        sc.pl.umap(a, color=col, ax=ax, show=False, frameon=False,
                   title=title)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figA_umap_overview.png"), dpi=200)
    plt.close(fig)
    log("wrote figA_umap_overview.png")

    # --- Fig B: T-cell subsets (Fig. 2f analogue), classic, renormalised ---
    t = tnk_renorm[tnk_renorm["chemistry"] == "classic"].set_index("genotype")
    subsets = [c for c in tnk_renorm.columns
               if c not in ("filename", "genotype", "chemistry", "n")]
    mat = t[subsets].T.loc[subsets, geno_order]
    x = np.arange(len(subsets))
    fig, ax = plt.subplots(figsize=(9.6, 4.4))
    for i, g in enumerate(geno_order):
        ax.bar(x + (i - 1) * w, mat[g].to_numpy(), w, label=g,
               color=geno_colors[g])
    ax.set_xticks(x)
    ax.set_xticklabels(subsets, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("% of CD4+CD8+ T cells (classic)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figB_tcell_subsets.png"), dpi=200)
    plt.close(fig)
    log("wrote figB_tcell_subsets.png")

    # --- Fig C: myeloid subsets (Fig. 4g analogue), classic vs Flex ---
    key = ["Spp1/Marco+", "Cxcl9+", "Proliferating", "Arg1+",
           "C1q+ homeostatic", "Folr2+"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for ax, chem in zip(axes, ["classic", "Flex"]):
        m = mac_frac[mac_frac["chemistry"] == chem].set_index("genotype")
        mm = m[key].T.loc[key, geno_order]
        x = np.arange(len(key))
        for i, g in enumerate(geno_order):
            ax.bar(x + (i - 1) * w, mm[g].to_numpy(), w, label=g,
                   color=geno_colors[g])
        ax.set_xticks(x)
        ax.set_xticklabels(key, rotation=30, ha="right", fontsize=8)
        ax.set_title(chem)
        ax.set_ylabel("% of myeloid cells")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figC_macro_subsets.png"), dpi=200)
    plt.close(fig)
    log("wrote figC_macro_subsets.png")

    # --- Fig D: CAF subsets (ED Fig. 7c analogue), Flex, stacked ---
    c = caf_frac[caf_frac["chemistry"] == "Flex"].set_index("genotype")
    csub = [s for s in caf_frac.columns
            if s not in ("filename", "genotype", "chemistry", "n")]
    cm = c[csub].T.loc[csub, geno_order]
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    bottom = np.zeros(len(geno_order))
    cmap = plt.get_cmap("tab20")
    for i, s in enumerate(csub):
        ax.bar(geno_order, cm.loc[s].to_numpy(), bottom=bottom, label=s,
               color=cmap(i % 20))
        bottom = bottom + cm.loc[s].to_numpy()
    ax.set_ylabel("% of Flex CAFs")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figD_caf_subsets.png"), dpi=200)
    plt.close(fig)
    log("wrote figD_caf_subsets.png")

    # --- Fig E: serpin detection in malignant epithelium (mosaic) ---
    libs = tab.sort_values(["chemistry", "genotype"]).reset_index(drop=True)
    x = np.arange(len(libs))
    wb = 0.38
    fig, ax = plt.subplots(figsize=(9, 4.4))
    ax.bar(x - wb / 2, libs["mal_serpine1_pct"], wb, color="black",
           label="Serpine1+")
    ax.bar(x + wb / 2, libs["mal_serpinb2_pct"], wb, color="black",
           alpha=0.3, hatch="//", label="Serpinb2+")
    ymax = float(libs["mal_serpine1_pct"].max()) * 1.25
    ax.set_ylim(0, ymax)
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("GSM", "").replace("_", "\n")
                        for s in libs["filename"]], fontsize=7)
    ax.set_ylabel("% malignant epithelial (>=2 UMI)")
    n_classic = int((libs["chemistry"] == "classic").sum())
    ax.axvline(n_classic - 0.5, color="k", ls=":", lw=1)
    ax.text((n_classic - 1) / 2, ymax * 0.94, "classic", ha="center")
    ax.text(n_classic + (len(libs) - n_classic - 1) / 2, ymax * 0.94,
            "Flex", ha="center")
    handles = [plt.Rectangle((0, 0), 1, 1, color=geno_colors[g])
               for g in geno_order]
    ax.legend(handles, geno_order, frameon=False, fontsize=8,
              loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figE_serpin_mosaic.png"), dpi=200)
    plt.close(fig)
    log("wrote figE_serpin_mosaic.png")

    # --- Fig F: genotype vs chemistry effect ratio ---
    e = eff.sort_values("ratio_geno_over_chem")
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.barh(e["readout"], e["ratio_geno_over_chem"].to_numpy(),
            color="#5a5a5a")
    ax.axvline(1.0, color="#d62728", ls="--", lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("genotype effect / chemistry effect  (< 1: chemistry "
                  "dominates)", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figF_effect_ratio.png"), dpi=200)
    plt.close(fig)
    log("wrote figF_effect_ratio.png")


# --------------------------------------------------------------------------
# 4. Source-data audit of the non-deposited panels (added 2026-11-15).
#
# Nature publishes per-figure Source Data spreadsheets (MOESM9-17) with the
# per-mouse values behind Figs 2-5 and ED Figs 4/5/7/8.  Fig 1 and ED
# Figs 1/3/6/9/10 have NO source data file (the Perturb-map screen itself).
# Here we re-plot the deposited per-mouse values and recompute the tests the
# legends claim (Welch's t; log-rank where the sheet encodes all events).
# --------------------------------------------------------------------------

SD_DIR = os.path.join(ROOT, "serpin_repro", "data", "Nature_SourceData")
SD_FILES = {
    "fig2": "41586_2026_11002_MOESM9_ESM.xlsx",
    "fig3": "41586_2026_11002_MOESM10_ESM.xlsx",
    "fig4": "41586_2026_11002_MOESM11_ESM.xlsx",
    "fig5": "41586_2026_11002_MOESM12_ESM.xlsx",
    "ed4": "41586_2026_11002_MOESM14_ESM.xlsx",
    "ed5": "41586_2026_11002_MOESM15_ESM.xlsx",
    "ed7": "41586_2026_11002_MOESM16_ESM.xlsx",
    "ed8": "41586_2026_11002_MOESM17_ESM.xlsx",
}


def _sd_groups(fn_key, sheet):
    """Parse a Source Data sheet -> {group: per-mouse values}. Header row 1,
    data from row 2; survival sheets handled separately."""
    import pandas as pd
    xl = pd.ExcelFile(os.path.join(SD_DIR, SD_FILES[fn_key]))
    df = xl.parse(sheet, header=None)
    header = [None if pd.isna(v) else str(v).replace("\xa0", " ").strip()
              for v in df.iloc[1]]
    groups = {}
    for j, name in enumerate(header):
        if name is None or name in ("", "Survival"):
            continue
        vals = pd.to_numeric(df.iloc[2:, j], errors="coerce").dropna()
        if len(vals):
            groups[name] = vals.astype(float).to_numpy()
    return groups


def _sd_survival(fn_key, sheet):
    """Parse a survival sheet: col 1 = day, group columns flag events."""
    import pandas as pd
    xl = pd.ExcelFile(os.path.join(SD_DIR, SD_FILES[fn_key]))
    df = xl.parse(sheet, header=None)
    header = [None if pd.isna(v) else str(v).replace("\xa0", " ").strip()
              for v in df.iloc[1]]
    out = {}
    for j, name in enumerate(header):
        if name is None or name in ("", "Survival"):
            continue
        rows = df.iloc[2:, [1, j]].apply(pd.to_numeric, errors="coerce")
        events = rows[(rows[j] == 1) & rows[1].notna()][1].astype(float).tolist()
        out[name] = np.array(events)
    return out


def _welch_p(a, b):
    from scipy import stats
    return float(stats.ttest_ind(a, b, equal_var=False)[1])


def _logrank_nocensor(t1, t2):
    """Mantel-Cox log-rank with no censoring (Source Data sheets encode one
    row per event and the row counts equal the legend n per group)."""
    from scipy.stats import chi2
    t1, t2 = np.asarray(t1), np.asarray(t2)
    times = np.sort(np.unique(np.concatenate([t1, t2])))
    O1 = E1 = V = 0.0
    for t in times:
        at1 = int((t1 >= t).sum())
        at2 = int((t2 >= t).sum())
        n_tot = at1 + at2
        if n_tot == 0:
            continue
        d1 = int((t1 == t).sum())
        d_tot = d1 + int((t2 == t).sum())
        E1 += d_tot * at1 / n_tot
        O1 += d1
        if n_tot > 1:
            V += d_tot * (at1 / n_tot) * (at2 / n_tot) * (n_tot - d_tot) / (n_tot - 1)
    z = (O1 - E1) / np.sqrt(V) if V > 0 else 0.0
    return float(chi2.sf(z ** 2, 1))


def _sd_panel(ax, groups, order, ylabel, title, pairs, claimed,
              colors=None, logy=False):
    """Grouped bar + per-mouse points; annotate Welch's p (recomputed)."""
    import matplotlib.pyplot as plt
    if colors is None:
        colors = {"Control": "#4d4d4d", "Serpinb2 KO": "#1f77b4",
                  "Serpine1 KO": "#d62728"}
    xs = np.arange(len(order))
    for i, g in enumerate(order):
        v = groups[g]
        ax.bar(i, float(np.mean(v)), 0.62, color=colors.get(g, "#999999"),
               alpha=0.85)
        ax.scatter(np.full(len(v), i) + np.random.uniform(-0.12, 0.12, len(v)),
                   v, s=14, color="white", edgecolor="black", zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(order, rotation=20, ha="right", fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=8)
    if logy:
        ax.set_yscale("log")
    for k, (a, b) in enumerate(pairs):
        if a not in groups or b not in groups:
            continue
        p = _welch_p(groups[a], groups[b])
        y = max(groups[a].max(), groups[b].max())
        ax.text(len(order) - 0.4, y * (1.35 - 0.18 * k),
                "%s vs %s: p=%.4g%s" % (b, a, p,
                                        claimed.get((a, b), "")),
                fontsize=6.5)
    ax.tick_params(labelsize=7)


def source_data_audit():
    """Re-plot + re-test the panels whose Source Data Nature deposits."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    os.makedirs(FIG_DIR, exist_ok=True)
    log("=" * 84)
    log("Source-data audit (Nature MOESM9-17; per-mouse values)")
    log("=" * 84)
    results = []

    # ---- Fig 2d: tumour mass -------------------------------------------
    g = _sd_groups("fig2", "Figure_2d")
    order = ["Control", "Serpinb2 KO", "Serpine1 KO"]
    pB = _welch_p(g["Control"], g["Serpinb2 KO"])
    pE = _welch_p(g["Control"], g["Serpine1 KO"])
    results.append(("Fig 2d tumour mass", "C vs B", "Welch", "0.001", "%.4g" % pB))
    results.append(("Fig 2d tumour mass", "C vs E", "Welch", "0.0003", "%.4g" % pE))
    results.append(("Fig 2d burden ratio", "B/C, E/C", "means",
                    ">50% cut", "%.0f%%, %.0f%%" % (
                        100 * np.mean(g["Serpinb2 KO"]) / np.mean(g["Control"]),
                        100 * np.mean(g["Serpine1 KO"]) / np.mean(g["Control"]))))

    # ---- Fig 2e: CD8 flow ----------------------------------------------
    g = _sd_groups("fig2", "Figure_2e")
    pB = _welch_p(g["Control"], g["Serpinb2 KO"])
    pE = _welch_p(g["Control"], g["Serpine1 KO"])
    results.append(("Fig 2e CD8 flow", "C vs B", "Welch", "0.0436", "%.4g" % pB))
    results.append(("Fig 2e CD8 flow", "C vs E", "Welch", "-", "%.4g" % pE))
    results.append(("Fig 2e fold", "B/C, E/C", "means", ">2x",
                    "%.2fx, %.2fx" % (
                        np.mean(g["Serpinb2 KO"]) / np.mean(g["Control"]),
                        np.mean(g["Serpine1 KO"]) / np.mean(g["Control"]))))

    # ---- Fig 2h / 2i: IHC counts ---------------------------------------
    for sheet, label, claimed in [("Figure_2h", "CD8 IHC", "0.0064"),
                                  ("Figure_2i", "GZMB IHC", "0.0018/0.041/0.0194")]:
        g = _sd_groups("fig2", sheet)
        names = list(g.keys())
        for other in names[1:]:
            results.append((f"Fig {label}", f"C vs {other}", "Welch", claimed,
                            "%.4g" % _welch_p(g["Control"], g[other])))

    # ---- Fig 2k: NSG mass; 4a: macrophage flow; 3e/f/j/l; 5e-h ----------
    for fn_key, sheet, label, claimed in [
            ("fig2", "Figure_2k", "NSG mass", "0.0043/0.0275"),
            ("fig4", "Figure_4a", "macrophage flow", "0.0219"),
            ("fig3", "Figure_3e", "PDPN+ CAF", "0.2218 (B) / 0.002 (E)"),
            ("fig3", "Figure_3f", "aSMA+PDPN CAF", "0.7267 (B) / 0.0359 (E)"),
            ("fig3", "Figure_3j", "trichrome", "0.0564 (B) / 0.0282 (E)"),
            ("fig3", "Figure_3l", "fibrin area", "0.0482 (B) / 0.0145 (E)"),
            ("fig5", "Figure_5e", "mix mass", "0.0578"),
            ("fig5", "Figure_5f", "mix macrophages", "0.0005"),
            ("fig5", "Figure_5g", "mix CD8", "0.9146"),
            ("fig5", "Figure_5h", "mix GZMB", "0.0064")]:
        g = _sd_groups(fn_key, sheet)
        names = list(g.keys())
        base = names[0]
        for other in names[1:]:
            results.append((f"{label}", f"{base} vs {other}", "Welch", claimed,
                            "%.4g" % _welch_p(g[base], g[other])))

    # ---- Survival: 2l (6 groups), 2n (4 groups) ------------------------
    s = _sd_survival("fig2", "Figure_2l")
    results.append(("Fig 2l median (d)", "Control aPD1", "median", "24",
                    "%.0f" % np.median(s["Control aPD1"])))
    for gname in ["Serpinb2 KO aPD1", "Serpine1 KO aPD1"]:
        results.append(("Fig 2l median (d)", gname, "median", "47",
                        "%.0f" % np.median(s[gname])))
        p = _logrank_nocensor(s["Control aPD1"], s[gname])
        results.append(("Fig 2l log-rank", "aPD1: C vs %s" % gname.split()[0],
                        "log-rank", "0.0002/0.0005", "%.4g" % p))
    s = _sd_survival("fig2", "Figure_2n")
    base = "Vehicle control"
    for gname in ["PAI-039", "aPD1", "PAI-039+aPD1"]:
        p = _logrank_nocensor(s[base], s[gname])
        results.append(("Fig 2n log-rank", f"{base} vs {gname}", "log-rank",
                        "0.0003/0.0024/0.0143", "%.4g" % p))
    # remaining pairwise contrasts (claimed 0.0024 may be mono-therapy vs combo)
    results.append(("Fig 2n log-rank", "aPD1 vs PAI-039+aPD1", "log-rank",
                    "0.0003/0.0024/0.0143",
                    "%.4g" % _logrank_nocensor(s["aPD1"], s["PAI-039+aPD1"])))
    results.append(("Fig 2n log-rank", "PAI-039 vs PAI-039+aPD1", "log-rank",
                    "0.0003/0.0024/0.0143",
                    "%.4g" % _logrank_nocensor(s["PAI-039"],
                                               s["PAI-039+aPD1"])))

    # ---- ED 5d / ED 7b: secondary panels -------------------------------
    g = _sd_groups("ed5", "Ext_Data_Figure_5d")
    names = list(g.keys())
    for other in names[1:]:
        results.append(("ED 5d D14 weight", f"{names[0]} vs {other}", "Welch",
                        "-", "%.4g" % _welch_p(g[names[0]], g[other])))

    res = pd.DataFrame(results, columns=["panel", "contrast", "test",
                                         "claimed p", "recomputed"])
    log(res.to_string(index=False))
    res.to_csv(os.path.join(FIG_DIR, "source_data_audit.tsv"), sep="\t",
               index=False)
    log("wrote source_data_audit.tsv")

    # ------------------------------- figures ----------------------------
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
    g = _sd_groups("fig2", "Figure_2d")
    _sd_panel(axes[0], g, list(g.keys()), "tumour mass (g)", "Fig 2d",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .001)",
               ("Control", "Serpine1 KO"): " (claim 3e-4)"})
    g = _sd_groups("fig2", "Figure_2e")
    _sd_panel(axes[1], g, list(g.keys()), "CD8+ per mg", "Fig 2e",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .044)"})
    g = _sd_groups("fig2", "Figure_2h")
    _sd_panel(axes[2], g, list(g.keys()), "CD8+ per mm2", "Fig 2h",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .0064)"})
    g = _sd_groups("fig2", "Figure_2i")
    _sd_panel(axes[3], g, list(g.keys()), "GZMB+CD8+ per mm2", "Fig 2i",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .0018)"})
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figG_sd_fig2_barpanels.png"), dpi=200)
    plt.close(fig)
    log("wrote figG_sd_fig2_barpanels.png")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))
    km_colors = {"Control": "#4d4d4d", "Serpinb2 KO": "#1f77b4",
                 "Serpine1 KO": "#d62728"}
    for ax, (groups, treat, panel) in zip(
            axes[:2],
            [(_sd_survival("fig2", "Figure_2l"), "IgG", "Fig 2l left"),
             (_sd_survival("fig2", "Figure_2l"), "aPD1", "Fig 2l right")]):
        for gname, t in groups.items():
            if not gname.endswith(treat):
                continue
            t = np.sort(t)
            ys = 100.0 * (1.0 - np.arange(1, len(t) + 1) / len(t))
            ax.step(np.concatenate([[0], t]), np.concatenate([[100], ys]),
                    where="post", color=km_colors[gname.rsplit(" " + treat)[0]],
                    label=gname)
        ax.set_xlabel("days"); ax.set_ylabel("survival (%)")
        ax.set_title(panel, fontsize=9); ax.legend(fontsize=6)
        ax.set_ylim(0, 105)
    ax = axes[2]
    for gname, t in _sd_survival("fig2", "Figure_2n").items():
        t = np.sort(t)
        ys = 100.0 * (1.0 - np.arange(1, len(t) + 1) / len(t))
        ax.step(np.concatenate([[0], t]), np.concatenate([[100], ys]),
                where="post", label=gname)
    ax.set_xlabel("days"); ax.set_ylabel("survival (%)")
    ax.set_title("Fig 2n (KC, PAI-039/aPD1)", fontsize=9)
    ax.legend(fontsize=6); ax.set_ylim(0, 105)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figH_sd_survival.png"), dpi=200)
    plt.close(fig)
    log("wrote figH_sd_survival.png")

    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
    g = _sd_groups("fig4", "Figure_4a")
    _sd_panel(axes[0], g, list(g.keys()), "macrophages per mg",
              "Fig 4a", [("Control", "Serpinb2 KO"),
                         ("Control", "Serpine1 KO")],
              {("Control", "Serpine1 KO"): " (claim .0219)"}, logy=True)
    g = _sd_groups("fig3", "Figure_3e")
    _sd_panel(axes[1], g, list(g.keys()), "PDPN+ (%)", "Fig 3e",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .2218)",
               ("Control", "Serpine1 KO"): " (claim .002)"})
    g = _sd_groups("fig3", "Figure_3j")
    _sd_panel(axes[2], g, list(g.keys()), "trichrome+ (%)", "Fig 3j",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .0564)",
               ("Control", "Serpine1 KO"): " (claim .0282)"})
    g = _sd_groups("fig3", "Figure_3l")
    _sd_panel(axes[3], g, list(g.keys()), "fibrin(ogen)+ (%)", "Fig 3l",
              [("Control", "Serpinb2 KO"), ("Control", "Serpine1 KO")],
              {("Control", "Serpinb2 KO"): " (claim .0482)",
               ("Control", "Serpine1 KO"): " (claim .0145)"})
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figI_sd_fig4_fig3.png"), dpi=200)
    plt.close(fig)
    log("wrote figI_sd_fig4_fig3.png")

    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
    for i, (sheet, lab) in enumerate(
            [("Figure_5e", "tumour mass (g)"),
             ("Figure_5f", "macrophages (%)"),
             ("Figure_5g", "CD8+ T cells (%)"),
             ("Figure_5h", "GZMB+CD8+ (%)")]):
        g = _sd_groups("fig5", sheet)
        names = list(g.keys())
        colors5 = {names[0]: "#4d4d4d", names[1]: "#d62728"}
        if len(names) > 2:
            colors5[names[2]] = "#f2a541"
        _sd_panel(axes[i], g, names, lab, sheet.replace("Figure_", "Fig "),
                  [(names[1], names[2])] if len(names) > 2 else [],
                  {}, colors=colors5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "figJ_sd_fig5_mixing.png"), dpi=200)
    plt.close(fig)
    log("wrote figJ_sd_fig5_mixing.png")

    log("source-data audit done")


def main():
    meta = load_metadata()
    a = load_processed(meta)
    a = attach_annotation(a)

    log("=" * 84)
    log("GSE289138 reproduction audit — genotype re-run")
    log("=" * 84)
    log(f"cells: {a.n_obs:,}   libraries: {a.obs['filename'].nunique()}   "
        f"tumours: {a.obs['sample'].nunique()}   genotypes: {a.obs['genotype'].nunique()}")
    log("design: 3 genotypes x 2 chemistries, ONE tumour per genotype")
    log("unit of inference: tumour (n = 1 per genotype) -> no cell-level p-values")
    log("")

    tab = per_library(a, meta)
    log("--- per-library readouts ---")
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        log(tab.to_string(index=False))
    log("")

    log("--- genotype vs chemistry effect size (descriptive) ---")
    log(f"{'readout':<24}{'genotype_eff':>14}{'chemistry_eff':>16}{'ratio':>9}")
    eff_rows = []
    for col in READOUTS:
        ge, ce = effects(tab, col)
        ratio = ge / ce if ce and np.isfinite(ce) and ce > 0 else np.nan
        eff_rows.append({"readout": col, "genotype_eff": ge,
                         "chemistry_eff": ce, "ratio_geno_over_chem": ratio})
        log(f"{col:<24}{ge:>14.3f}{ce:>16.3f}{ratio:>9.2f}")
    eff = pd.DataFrame(eff_rows)
    log("")

    # Direct checks of the source paper's directional claims.
    log("--- directional checks against the paper's claims ---")
    classic = tab[tab["chemistry"] == "classic"].set_index("genotype")
    flex = tab[tab["chemistry"] == "Flex"].set_index("genotype")

    def direction(frame, col, higher_in):
        vals = frame[col]
        return {gt: round(float(vals.get(gt, np.nan)), 3) for gt in GENOTYPES}

    log("T/NK exhaustion (classic; paper: lower in Serpinb2-KO and Serpine1-KO)")
    log("  " + str(direction(classic, "tnk_exhaust", None)))
    log("T/NK effector (classic; paper: higher in Serpinb2-KO and Serpine1-KO)")
    log("  " + str(direction(classic, "tnk_effector", None)))
    log("SPP1+MARCO+ macrophages (classic; paper: lower in Serpine1-KO)")
    log("  " + str(direction(classic, "mac_spp1_marco_pct", None)))
    log("CAF fraction (Flex; paper: lower in Serpine1-KO)")
    log("  " + str(direction(flex, "pct_CAF", None)))
    log("SERPINE1+ malignant epithelial (paper: rare, genotype-dependent)")
    log("  classic " + str(direction(classic, "mal_serpine1_pct", None)))
    log("  Flex    " + str(direction(flex, "mal_serpine1_pct", None)))
    log("")

    log("--- interpretation ---")
    log("1. Every readout is reported per library. With one tumour per genotype,")
    log("   differences cannot be separated from tumour identity.")
    log("2. Where the chemistry effect exceeds the genotype effect, a genotype")
    log("   claim from a pooled analysis is not identifiable.")
    log("3. The source repository contains figure-level notebooks only; the")
    log("   upstream processing that produced its intermediate h5ad is not")
    log("   deposited, so this is a reconstruction, not an exact replay.")

    # ---- subset-level reproduction (Fig. 2f / 4g / 3g-h analogues) ----
    log("")
    log("--- subset-level reproduction (Fig. 2f / 4g / 3g-h analogues) ---")
    tnk_frac, tnk_lab = tnk_subsets(a, meta)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        log(tnk_frac.to_string(index=False))
    # Fig. 2f's denominator excludes NK cells; renormalise over CD4+CD8 T cells
    tnk_renorm = tnk_frac.copy()
    sub_cols = [c for c in tnk_renorm.columns if c not in
                ("filename", "genotype", "chemistry", "n", "NK/other")]
    for c in sub_cols:
        tnk_renorm[c] = tnk_renorm[c] / (100.0 - tnk_renorm["NK/other"]) * 100.0
    tnk_renorm = tnk_renorm.drop(columns=["NK/other"]).round(2)
    log("T-cell subsets renormalised over CD4+CD8 T cells (NK/other excluded):")
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        log(tnk_renorm.to_string(index=False))
    mac_frac, mac_lab = macro_subsets(meta)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        log(mac_frac.to_string(index=False))
    caf_frac, caf_lab = caf_subsets(a, meta)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        log(caf_frac.to_string(index=False))
    endo_frac, endo_lab = endo_subsets(a, meta)
    if endo_frac is not None:
        with pd.option_context("display.width", 220, "display.max_columns", 30):
            log(endo_frac.to_string(index=False))

    log("--- subset directional checks (per chemistry) ---")
    t = tnk_frac[tnk_frac["chemistry"] == "classic"]
    for col in ["CD8 Teff", "CD8 Tex prog", "CD8 Tex int", "CD8 Tex term",
                "CD4 Treg"]:
        if col in t.columns:
            log("  T cells classic %s: %s" % (
                col, {r["genotype"]: r[col] for _, r in t.iterrows()}))
    m = mac_frac[mac_frac["chemistry"] == "Flex"]
    for col in ["Spp1/Marco+", "Cxcl9+", "Proliferating", "Arg1+"]:
        if col in m.columns:
            log("  macrophages Flex %s: %s" % (
                col, {r["genotype"]: r[col] for _, r in m.iterrows()}))
    c = caf_frac[caf_frac["chemistry"] == "Flex"]
    for col in ["iCAF", "MyCAF", "apCAF", "LRCC15+ MyCAF", "Matrix CAF"]:
        if col in c.columns:
            log("  CAF Flex %s: %s" % (
                col, {r["genotype"]: r[col] for _, r in c.iterrows()}))

    # ---- figures for reproduction_report.tex ----
    log("")
    log("--- figures for reproduction_report.tex ---")
    export_plot_tables(a, meta, tab, tnk_renorm, mac_frac, caf_frac, eff)
    make_figures(a, meta, tab, tnk_renorm, mac_frac, caf_frac, eff)
    log("figures written to " + FIG_DIR)


if __name__ == "__main__":
    if "--source-data" in sys.argv:
        source_data_audit()
    else:
        main()
