"""Shared IO helpers for the GSE300595 10x Multiome (cellranger-arc) h5 files.

The GEO supplementary files are `filtered_feature_bc_matrix.h5` produced by
cellranger-arc 2.0.2 (GRCh38).  Each file stores ONE matrix whose rows are the
concatenated Gene Expression genes and ATAC peaks; `/matrix/features/feature_type`
distinguishes them ("Gene Expression" vs "Peaks").

Only count matrices are distributed (no fragments.tsv.gz), therefore all
fragment-based QC metrics (TSS enrichment, nucleosome signal, blacklist
fraction, fraction of reads in peaks) CANNOT be computed for this dataset.
"""
import os
import numpy as np
import h5py
from scipy.sparse import csc_matrix, csr_matrix

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data")
PROC = os.path.join(BASE, "processed")
FIG = os.path.join(BASE, "figures")
for _d in (DATA, PROC, FIG):
    os.makedirs(_d, exist_ok=True)

SAMPLES = ["pc99", "pc102", "pc123", "pc131", "pc149", "pc158",
           "pc159", "pc167", "pc169", "pc172", "pc173", "pc187"]


def h5_path(sample):
    return os.path.join(DATA, f"GSE300595_{sample}.filtered_feature_bc_matrix.h5")


def read_multiome_h5(path):
    """Read a cellranger-arc filtered_feature_bc_matrix.h5.

    Returns
    -------
    X : scipy sparse (n_cells x n_features), features = genes then peaks
    barcodes : np.ndarray[str]
    feature_type : np.ndarray[str]
    feature_name : np.ndarray[str]
    feature_id : np.ndarray[str]
    """
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        shape = tuple(int(x) for x in m["shape"][:])
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]

        def _decode(dset):
            arr = dset[:]
            return np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x)
                             for x in arr])

        barcodes = _decode(m["barcodes"])
        feat = m["features"]
        feature_type = _decode(feat["feature_type"])
        feature_name = _decode(feat["name"])
        feature_id = _decode(feat["id"])

    # cellranger stores features x cells; figure out the orientation from indptr
    if len(indptr) == shape[1] + 1:
        X_fc = csc_matrix((data, indices, indptr), shape=shape)   # features x cells
    elif len(indptr) == shape[0] + 1:
        X_fc = csr_matrix((data, indices, indptr), shape=shape)   # features x cells
    else:
        raise ValueError(f"Unrecognized matrix layout in {path}: "
                         f"shape={shape}, len(indptr)={len(indptr)}")

    return X_fc.T.tocsr(), barcodes, feature_type, feature_name, feature_id


def split_features(X, feature_type):
    """Split the combined feature matrix into (genes, peaks)."""
    gex_mask = feature_type == "Gene Expression"
    peak_mask = feature_type == "Peaks"
    return X[:, gex_mask], X[:, peak_mask], gex_mask, peak_mask


def parse_peak(name):
    """'chr1:10000-10500' -> ('chr1', 10000, 10500)."""
    chrom, rest = name.rsplit(":", 1)
    start, end = rest.split("-")
    return chrom, int(start), int(end)
