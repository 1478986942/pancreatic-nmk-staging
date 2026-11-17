"""GSE300595 step 5 - myeloid ATAC: peak matrix assembly, LSI, DAR, annotation.

Builds a common peak space from the 12 per-sample cellranger-arc peak sets by
merging overlapping intervals across samples, then assembles the peak count
matrix for the MYELOID nuclei only (identified from RNA in step 4).

No fragments.tsv.gz is distributed for GSE300595, so peak calling is not
possible here; the union of the per-sample peak sets supplied by cellranger-arc
is used instead.  Consequently no chromVAR/JASPAR motif enrichment is computed
(that would require the genome sequence and fragment-level insertions).  DARs
are instead annotated to their nearest gene, which is what the downstream
bridge to the SERPIN-high niche needs.

Differential accessibility: patient-level pseudobulk (12 patients), CPM +
log1p, Mann-Whitney per peak (n=5 responders vs 7 non-responders).
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from scipy import stats
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import BASE, PROC, SAMPLES, h5_path, read_multiome_h5, split_features

GTF = os.path.join(BASE, "ref", "gencode.v44.basic.annotation.gtf.gz")


def load_genes():
    """gene -> (chrom, start, end, strand)"""
    g = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != "gene":
                continue
            name = None
            for kv in p[8].split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    name = kv.split(" ", 1)[1].strip('"')
                    break
            if name and name not in g:
                g[name] = (p[0], int(p[3]), int(p[4]), p[6])
    return g


def main():
    # ---- myeloid barcodes per sample --------------------------------------
    my = sc.read_h5ad(os.path.join(PROC, "myeloid_annotated.h5ad"))
    my_obs = my.obs[["sample", "response", "myeloid_state"]].copy()
    bc_by_sample = {s: set(g.index.str.split("_", n=1).str[1])
                    for s, g in my_obs.groupby("sample", observed=True)}
    print("Myeloid cells per sample:", flush=True)
    print(my_obs.groupby("sample", observed=True).size().to_string(), flush=True)

    # ---- merged peak space ------------------------------------------------
    peaks = pd.read_csv(os.path.join(PROC, "merged_peaks.bed"), sep="\t")
    print(f"\nMerged peak space: {len(peaks)}", flush=True)

    peak_offsets = {}
    for chrom, grp in peaks.groupby("chrom", sort=False):
        idx = grp.index.values
        order = np.argsort(grp["start"].values)
        peak_offsets[chrom] = (grp["start"].values[order], idx[order])

    blocks = []
    obs_rows = []
    for s in SAMPLES:
        X, barcodes, ftype, fname, fid = read_multiome_h5(h5_path(s))
        _, Xp, _, pm = split_features(X, ftype)
        del X
        names = fname[pm]
        keep = np.isin(barcodes, list(bc_by_sample[s]))
        if keep.sum() == 0:
            print(f"  [{s}] no myeloid cells", flush=True)
            continue
        Xp = Xp[keep].tocsr()
        cells = barcodes[keep]
        print(f"  [{s}] {Xp.shape[0]} myeloid nuclei x {Xp.shape[1]} peaks", flush=True)

        # map this sample's peaks to the merged index
        chroms = np.array([n.split(":", 1)[0] for n in names])
        mapping = np.full(len(names), -1, dtype=np.int64)
        for chrom, (starts, gidx) in peak_offsets.items():
            sel = np.where(chroms == chrom)[0]
            if len(sel) == 0:
                continue
            ends = np.array([int(names[i].rsplit("-", 1)[1]) for i in sel])
            j = np.searchsorted(starts, np.array(
                [int(names[i].rsplit(":", 1)[1].split("-")[0]) for i in sel]),
                side="right") - 1
            valid = (j >= 0)
            jj = np.clip(j, 0, len(starts) - 1)
            valid &= (ends <= peaks.loc[gidx[jj], "end"].values)
            mapping[sel[valid]] = gidx[jj[valid]]

        ok = mapping >= 0
        M = sp.csr_matrix(
            (np.ones(ok.sum(), dtype=np.float32),
             (np.where(ok)[0], mapping[ok])),
            shape=(len(names), len(peaks)))
        Xm = (Xp @ M).tocsr()
        blocks.append(Xm)
        obs_rows.append(pd.DataFrame(
            {"cell_id": [f"{s}_{b}" for b in cells]}, index=[f"{s}_{b}" for b in cells]))

    X_all = sp.vstack(blocks, format="csr")
    obs = pd.concat(obs_rows)
    obs = obs.join(my_obs, how="left")
    print(f"\nAssembled myeloid ATAC: {X_all.shape}", flush=True)

    # ---- peak filtering ---------------------------------------------------
    n_cells = X_all.shape[0]
    detect = np.asarray((X_all > 0).sum(axis=0)).ravel()
    keep_peak = detect >= max(10, int(0.05 * n_cells))
    print(f"Peaks detected in >=5% of myeloid nuclei: {int(keep_peak.sum())}", flush=True)

    atac = ad.AnnData(X=X_all[:, keep_peak].tocsr(), obs=obs,
                      var=pd.DataFrame(index=peaks.loc[keep_peak, "peak_id"].values))
    atac.var["chrom"] = peaks.loc[keep_peak, "chrom"].values
    atac.var["start"] = peaks.loc[keep_peak, "start"].values
    atac.var["end"] = peaks.loc[keep_peak, "end"].values
    atac.layers["counts"] = atac.X.copy()

    # ---- TF-IDF + LSI -----------------------------------------------------
    Xbin = (atac.X > 0).astype(np.float32)
    tf = sp.diags(1.0 / np.maximum(np.asarray(Xbin.sum(axis=1)).ravel(), 1)) @ Xbin
    idf = np.log1p(n_cells / (1.0 + np.asarray(Xbin.sum(axis=0)).ravel()))
    tfidf = normalize(tf.multiply(idf), norm="l2", axis=1)
    svd = TruncatedSVD(n_components=40, random_state=0)
    lsi = svd.fit_transform(tfidf)
    atac.obsm["X_lsi"] = lsi[:, 1:]          # drop depth-dominated component 0
    print(f"LSI explained variance (first 10): "
          f"{np.round(svd.explained_variance_ratio_[:10], 3)}", flush=True)

    sc.pp.neighbors(atac, use_rep="X_lsi", n_neighbors=15, random_state=0)
    sc.tl.leiden(atac, resolution=0.5, key_added="atac_cluster",
                 random_state=0, flavor="igraph", n_iterations=2, directed=False)
    print(f"ATAC clusters: {atac.obs['atac_cluster'].nunique()}", flush=True)
    print(atac.obs.groupby(["atac_cluster", "myeloid_state"],
                           observed=True).size().to_string(), flush=True)

    # ---- DAR: patient pseudobulk -----------------------------------------
    samples = atac.obs["sample"].values
    uniq = sorted(pd.unique(samples))
    pb = np.vstack([np.asarray(atac.X[samples == s].sum(axis=0)).ravel()
                    for s in uniq])
    pb_df = pd.DataFrame(pb, index=uniq, columns=atac.var_names)
    pb_df.to_csv(os.path.join(PROC, "myeloid_atac_pseudobulk_counts.tsv"), sep="\t")

    cpm = pb_df.div(pb_df.sum(axis=1).replace(0, np.nan), axis=0) * 1e6
    logcpm = np.log1p(cpm.fillna(0))

    samp_resp = pd.Series({s: atac.obs.loc[atac.obs["sample"] == s,
                                           "response"].iloc[0] for s in uniq})
    r_s = list(samp_resp[samp_resp == "responder"].index)
    nr_s = list(samp_resp[samp_resp == "nonresponder"].index)

    det = (cpm > 5).sum(axis=0) >= 6
    expr = logcpm.loc[r_s + nr_s, det]
    print(f"\nPeaks tested for DAR: {expr.shape[1]} "
          f"(R n={len(r_s)}, NR n={len(nr_s)})", flush=True)

    res = []
    for pk in expr.columns:
        v1 = expr.loc[r_s, pk].values
        v0 = expr.loc[nr_s, pk].values
        if np.all(v1 == 0) and np.all(v0 == 0):
            continue
        try:
            u, p = stats.mannwhitneyu(v1, v0, alternative="two-sided")
        except ValueError:
            continue
        m1 = cpm.loc[r_s, pk].mean() + 1.0
        m0 = cpm.loc[nr_s, pk].mean() + 1.0
        res.append(dict(peak=pk, cpm_R=m1 - 1, cpm_NR=m0 - 1,
                        log2FC=np.log2(m1 / m0), p=p))
    dar = pd.DataFrame(res)
    from statsmodels.stats.multitest import multipletests
    dar["padj"] = multipletests(dar["p"], method="fdr_bh")[1]
    dar = dar.merge(peaks[["peak_id", "chrom", "start", "end"]],
                    left_on="peak", right_on="peak_id", how="left")
    dar = dar.sort_values("p")
    dar.to_csv(os.path.join(PROC, "myeloid_DAR_pseudobulk.tsv"), sep="\t", index=False)

    print(f"  DAR padj<0.25: {int((dar['padj'] < 0.25).sum())} | "
          f"padj<0.1: {int((dar['padj'] < 0.1).sum())}", flush=True)
    print("\n  Top 20 peaks MORE ACCESSIBLE in RESPONDERS (uncorrected p):",
          flush=True)
    print(dar.sort_values("log2FC", ascending=False).head(20)
          [["peak", "chrom", "start", "end", "cpm_R", "cpm_NR", "log2FC", "p", "padj"]]
          .to_string(index=False), flush=True)
    print("\n  Top 20 peaks MORE ACCESSIBLE in NON-RESPONDERS (uncorrected p):",
          flush=True)
    print(dar.sort_values("log2FC", ascending=True).head(20)
          [["peak", "chrom", "start", "end", "cpm_R", "cpm_NR", "log2FC", "p", "padj"]]
          .to_string(index=False), flush=True)

    # ---- nearest-gene annotation -----------------------------------------
    print("\n=== nearest-gene annotation of DARs ===", flush=True)
    genes = load_genes()
    gdf = pd.DataFrame([(g, v[0], v[1], v[2]) for g, v in genes.items()],
                       columns=["gene", "chrom", "start", "end"])
    gchrom = {c: sub for c, sub in gdf.groupby("chrom", sort=False)}

    def nearest(chrom, s, e):
        if chrom not in gchrom:
            return None, np.nan
        sub = gchrom[chrom]
        mid = (s + e) // 2
        d = np.abs(((sub["start"].values + sub["end"].values) // 2) - mid)
        i = int(np.argmin(d))
        return sub["gene"].values[i], float(d[i])

    ann = dar.copy()
    ng = [nearest(r.chrom, r.start, r.end) if isinstance(r.chrom, str) else (None, np.nan)
          for r in ann.itertuples()]
    ann["nearest_gene"] = [a[0] for a in ng]
    ann["dist_to_tss_bp"] = [a[1] for a in ng]
    ann.to_csv(os.path.join(PROC, "myeloid_DAR_annotated.tsv"), sep="\t", index=False)

    # gene-level programme: top DAR genes per direction
    for direction, sub in [("RESPONDER", ann.sort_values("log2FC", ascending=False).head(300)),
                           ("NONRESPONDER", ann.sort_values("log2FC", ascending=True).head(300))]:
        genes_dir = [g for g in sub["nearest_gene"].dropna().unique() if g]
        pd.Series(genes_dir).to_csv(
            os.path.join(PROC, f"DAR_genes_{direction}.txt"),
            index=False, header=False)
        print(f"  {direction}: {len(genes_dir)} unique nearest genes "
              f"(top: {genes_dir[:12]})", flush=True)

    atac.write_h5ad(os.path.join(PROC, "myeloid_atac.h5ad"), compression="gzip")
    print("\n[saved] myeloid_atac.h5ad", flush=True)


if __name__ == "__main__":
    main()
