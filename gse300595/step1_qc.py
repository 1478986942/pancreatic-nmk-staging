"""GSE300595 step 1 - QC of the 12 10x Multiome (cellranger-arc) samples.

Paper thresholds (Ojetunde et al., Clin Transl Med 2026):
    nCount_ATAC > 400, nCount_RNA > 400
    nFeature_ATAC > 200, nFeature_RNA > 200
    percent.mt (RNA) < 1 %, percent.ribo < 2 %
    percent.mt.atac < 5 %
    nucleosome_signal < 1.5, TSS.enrichment 2-10,
    blacklist_fraction < 2 %, pct_reads_in_peaks > 10 %

IMPORTANT: GEO distributes ONLY filtered_feature_bc_matrix.h5 for this series
(no fragments.tsv.gz / BAM).  The four fragment-based metrics --
nucleosome_signal, TSS.enrichment, blacklist_fraction, pct_reads_in_peaks --
therefore CANNOT be computed here and are not applied.  Every metric that can
be derived from the count matrix is applied exactly as in the paper.

Outputs (processed/):
    rna_qc.h5ad        cells x genes, raw counts, QC-passed nuclei
    cell_qc.tsv        per-cell QC metrics + sample/response annotation
    qc_summary.tsv     per-sample nuclei counts before/after filtering
    merged_peaks.bed   union peak set used later for ATAC (step 3)
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import (BASE, DATA, PROC, SAMPLES, h5_path, read_multiome_h5,
                      split_features, parse_peak)

META = pd.read_csv(os.path.join(BASE, "metadata_samples.tsv"), sep="\t")
META = META.set_index("sample")

# ---- paper QC thresholds (matrix-computable subset) ------------------------
MIN_COUNT_RNA = 400
MIN_FEAT_RNA = 200
MIN_COUNT_ATAC = 400
MIN_FEAT_ATAC = 200
MAX_PCT_MT = 1.0
MAX_PCT_MT_ATAC = 5.0
MAX_PCT_RIBO = 2.0

THRESHOLDS = dict(
    nCount_RNA=f"> {MIN_COUNT_RNA}", nFeature_RNA=f"> {MIN_FEAT_RNA}",
    percent_mt=f"< {MAX_PCT_MT}%", percent_ribo=f"< {MAX_PCT_RIBO}%",
    nCount_ATAC=f"> {MIN_COUNT_ATAC}", nFeature_ATAC=f"> {MIN_FEAT_ATAC}",
    percent_mt_atac=f"< {MAX_PCT_MT_ATAC}%",
    NOT_COMPUTABLE="nucleosome_signal, TSS.enrichment, blacklist_fraction, "
                   "pct_reads_in_peaks (require fragments.tsv.gz, not distributed)",
)


def main():
    t_start = time.time()
    adatas = []
    qc_rows = []
    all_peaks = []          # (chrom, start, end) over ALL cells, used for union

    for s in SAMPLES:
        t0 = time.time()
        path = h5_path(s)
        X, barcodes, ftype, fname, fid = read_multiome_h5(path)
        Xg, Xp, gm, pm = split_features(X, ftype)
        gene_names = fname[gm]
        peak_names = fname[pm]

        # ---- RNA QC metrics -------------------------------------------------
        n_count_rna = np.asarray(Xg.sum(axis=1)).ravel()
        n_feat_rna = np.asarray((Xg > 0).sum(axis=1)).ravel()

        is_mt = np.char.startswith(gene_names.astype(str), "MT-")
        is_ribo = (np.char.startswith(gene_names.astype(str), "RPS") |
                   np.char.startswith(gene_names.astype(str), "RPL"))
        mt_counts = np.asarray(Xg[:, is_mt].sum(axis=1)).ravel()
        ribo_counts = np.asarray(Xg[:, is_ribo].sum(axis=1)).ravel()
        pct_mt = np.where(n_count_rna > 0, 100.0 * mt_counts / n_count_rna, 0.0)
        pct_ribo = np.where(n_count_rna > 0, 100.0 * ribo_counts / n_count_rna, 0.0)

        # ---- ATAC QC metrics ------------------------------------------------
        n_count_atac = np.asarray(Xp.sum(axis=1)).ravel()
        n_feat_atac = np.asarray((Xp > 0).sum(axis=1)).ravel()
        chroms = np.array([p.split(":", 1)[0] for p in peak_names])
        mt_peak_mask = chroms == "chrM"
        mt_atac_counts = np.asarray(Xp[:, mt_peak_mask].sum(axis=1)).ravel()
        pct_mt_atac = np.where(n_count_atac > 0,
                               100.0 * mt_atac_counts / n_count_atac, 0.0)

        # ---- joint filter (one barcode set for both modalities) -------------
        keep = ((n_count_rna > MIN_COUNT_RNA) & (n_feat_rna > MIN_FEAT_RNA) &
                (pct_mt < MAX_PCT_MT) & (pct_ribo < MAX_PCT_RIBO) &
                (n_count_atac > MIN_COUNT_ATAC) & (n_feat_atac > MIN_FEAT_ATAC) &
                (pct_mt_atac < MAX_PCT_MT_ATAC))

        qc_rows.append(dict(
            sample=s, patient=META.loc[s, "patient"],
            response=META.loc[s, "response"], treatment=META.loc[s, "treatment"],
            n_barcodes=len(barcodes), n_pass=int(keep.sum()),
            pass_pct=round(100.0 * keep.sum() / len(barcodes), 2),
            median_nCount_RNA=int(np.median(n_count_rna[keep])) if keep.sum() else 0,
            median_nFeature_RNA=int(np.median(n_feat_rna[keep])) if keep.sum() else 0,
            median_nCount_ATAC=int(np.median(n_count_atac[keep])) if keep.sum() else 0,
            median_nFeature_ATAC=int(np.median(n_feat_atac[keep])) if keep.sum() else 0,
            median_pct_mt=round(float(np.median(pct_mt[keep])), 3) if keep.sum() else 0,
            median_pct_ribo=round(float(np.median(pct_ribo[keep])), 3) if keep.sum() else 0,
        ))

        # ---- keep QC-passed RNA ---------------------------------------------
        Xg_keep = Xg[keep].tocsr()
        obs = pd.DataFrame(dict(
            barcode=barcodes[keep],
            sample=s,
            patient=META.loc[s, "patient"],
            response=META.loc[s, "response"],
            treatment=META.loc[s, "treatment"],
            nCount_RNA=n_count_rna[keep],
            nFeature_RNA=n_feat_rna[keep],
            nCount_ATAC=n_count_atac[keep],
            nFeature_ATAC=n_feat_atac[keep],
            percent_mt=pct_mt[keep],
            percent_ribo=pct_ribo[keep],
            percent_mt_atac=pct_mt_atac[keep],
        ))
        obs["cell_id"] = obs["sample"] + "_" + obs["barcode"]

        a = ad.AnnData(X=Xg_keep, obs=obs.set_index("cell_id"),
                       var=pd.DataFrame(index=pd.Index(gene_names, name="gene")))
        a.var_names_make_unique()
        adatas.append(a)

        # union peak set uses all peaks (not only QC-passed cells)
        for p in peak_names:
            c, st, en = parse_peak(p)
            all_peaks.append((c, st, en))

        print(f"[{s}] barcodes={len(barcodes):6d} pass={int(keep.sum()):6d} "
              f"rna_nnz/cell={a.X.nnz//max(1,a.n_obs):5d} "
              f"({time.time()-t0:.1f}s)", flush=True)

    # ---- concatenate RNA ----------------------------------------------------
    rna = ad.concat(adatas, join="outer", merge="same", index_unique=None)
    rna.obs_names_make_unique()
    print(f"\nConcat RNA: {rna.shape[0]} nuclei x {rna.shape[1]} genes", flush=True)
    rna.write_h5ad(os.path.join(PROC, "rna_qc.h5ad"), compression="gzip")

    # ---- per-cell QC table --------------------------------------------------
    cell_qc = rna.obs.copy()
    cell_qc.to_csv(os.path.join(PROC, "cell_qc.tsv"), sep="\t")

    # ---- per-sample QC summary ---------------------------------------------
    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(os.path.join(PROC, "qc_summary.tsv"), sep="\t", index=False)
    print("\n" + qc_df.to_string(index=False), flush=True)
    print(f"\nTOTAL nuclei retained: {int(qc_df['n_pass'].sum())}", flush=True)

    # ---- merged union peak set ---------------------------------------------
    peaks = pd.DataFrame(all_peaks, columns=["chrom", "start", "end"])
    peaks = peaks.sort_values(["chrom", "start", "end"], kind="mergesort")

    merged = []
    for chrom, grp in peaks.groupby("chrom", sort=False):
        starts = grp["start"].values
        ends = grp["end"].values
        cur_s, cur_e = starts[0], ends[0]
        for st, en in zip(starts[1:], ends[1:]):
            if st <= cur_e:                       # overlapping -> merge
                if en > cur_e:
                    cur_e = en
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = st, en
        merged.append((chrom, cur_s, cur_e))

    merged_df = pd.DataFrame(merged, columns=["chrom", "start", "end"])
    merged_df["peak_id"] = (merged_df["chrom"] + ":" +
                            merged_df["start"].astype(str) + "-" +
                            merged_df["end"].astype(str))
    merged_df["width"] = merged_df["end"] - merged_df["start"]
    merged_df.to_csv(os.path.join(PROC, "merged_peaks.bed"), sep="\t", index=False)

    print(f"\nUnion peak set: {len(merged_df)} peaks "
          f"(from {len(peaks)} raw peaks across 12 samples)", flush=True)
    print(f"  width median={merged_df['width'].median():.0f} bp, "
          f"total={merged_df['width'].sum()/1e6:.1f} Mb", flush=True)

    # ---- thresholds record --------------------------------------------------
    with open(os.path.join(PROC, "qc_thresholds.txt"), "w") as fh:
        for k, v in THRESHOLDS.items():
            fh.write(f"{k}\t{v}\n")

    print(f"\nDone in {time.time()-t_start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
