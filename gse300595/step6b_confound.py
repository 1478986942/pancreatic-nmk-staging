"""GSE300595 step 6b - confound control for the DSP bridge.

Two problems surfaced in step 6 that must be resolved before any claim is made:

1. PSEUDO-REPLICATION.  The 570 "ROIs" are 235 physical ROIs that GeoMx split
   into geometric segments (193 ROIs -> 3 segments, 37 -> 2, 5 -> 1).  Treating
   570 segments as independent inflates significance.  Every test is therefore
   repeated on ROI-level aggregates (one value per physical ROI).

2. NON-SPECIFICITY.  The NR-enriched and R-enriched myeloid accessibility
   programmes BOTH score higher in SERPIN-high segments.  If both move together,
   the module score is tracking a generic axis (tumour/stroma content) rather
   than a response-specific myeloid programme.  This script quantifies that with
   partial Spearman correlations and pairwise programme correlation.

Nothing here is confirmatory: n=22 patients, cross-cohort, no shared patients.
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import PROC

DSP_DIR = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"
TRT_MAP = {"Untreated": "naive", "CRT": "treated",
           "CRTL": "treated", "CRTLN": "treated"}

CELL_SCORES = {
    "epithelial": ["KRT8", "KRT18", "MUC1", "MSLN", "CLDN4", "LGALS4",
                   "TSPAN8", "KRT17"],
    "CAF": ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "FBLN1", "PDGFRB",
            "ACTA2", "POSTN", "COL5A1"],
    "myeloid": ["CD68", "CD163", "CSF1R", "LYZ", "C1QA", "C1QB", "CD14",
                "ITGAM", "TYROBP", "FCER1G", "AIF1"],
    "Tcell": ["CD3D", "CD3E", "CD2", "NKG7", "PRF1"],
    "endo": ["VWF", "CDH5", "CLDN5", "ENG"],
}


def partial_spearman(x, y, z):
    """Spearman of x,y controlling for z (rank-based partial correlation)."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    n = len(x)
    if n < 8:
        return np.nan, np.nan, n
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    # residualise ranks on z
    def resid(a):
        b = np.polyfit(rz, a, 1)
        return a - np.polyval(b, rz)
    ex, ey = resid(rx), resid(ry)
    r = np.corrcoef(ex, ey)[0, 1]
    # t-test with df = n - 3 (two controls -> one here)
    dfree = n - 3
    t = r * np.sqrt(dfree / max(1e-12, 1 - r * r))
    p = 2 * stats.t.sf(abs(t), dfree)
    return r, p, n


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def main():
    cnt = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "q3norm_counts.txt.gz"), "rt"),
                      sep="\t", index_col=0)
    seg = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "segment_props.txt.gz"), "rt"),
                      sep="\t")
    seg["ROI_id"] = seg["Sample_ID"].str.replace("-", ".", regex=False)
    common = [r for r in cnt.columns if r in set(seg["ROI_id"])]
    meta = seg.set_index("ROI_id").loc[common].copy()
    X = np.log2(cnt.loc[:, common].T + 1.0)
    Z = ((X - X.mean()) / X.std().replace(0, np.nan)).fillna(0.0)
    panel = set(Z.columns)

    meta["status"] = meta["Treatment"].map(TRT_MAP).fillna("other")
    meta["phys_roi"] = meta["Slide_name"] + "_" + meta["ROI_number"].astype(str)

    # ---- cell-type context scores ----------------------------------------
    for name, genes in CELL_SCORES.items():
        have = [g for g in genes if g in panel]
        meta[f"ctx_{name}"] = Z[have].mean(axis=1)
    print("context scores built:", [c for c in meta.columns if c.startswith("ctx_")],
          flush=True)
    print("\ncontext score mean by Segment:", flush=True)
    print(meta.groupby("Segment")[[c for c in meta.columns if c.startswith("ctx_")]]
          .mean().round(3).to_string(), flush=True)

    # ---- programmes -------------------------------------------------------
    dar = pd.read_csv(os.path.join(PROC, "myeloid_DAR_annotated.tsv"), sep="\t")
    nr_genes = [l.strip() for l in
                open(os.path.join(PROC, "DAR_pcgenes_NONRESPONDER.txt")) if l.strip()]
    r_genes = [l.strip() for l in
               open(os.path.join(PROC, "DAR_pcgenes_RESPONDER.txt")) if l.strip()]
    nr_have = [g for g in nr_genes if g in panel]
    r_have = [g for g in r_genes if g in panel]
    meta["NR_score"] = Z[nr_have].mean(axis=1)
    meta["R_score"] = Z[r_have].mean(axis=1)
    meta["SERPINB2"] = Z["SERPINB2"]
    fam = [g for g in ["SERPINA1", "SERPINA3", "SERPINA5", "SERPINB1", "SERPINB2",
                       "SERPINB5", "SERPINB6", "SERPINB8", "SERPINF1", "SERPING1",
                       "SERPINH1"] if g in panel]
    meta["SERPIN_family"] = Z[fam].mean(axis=1)

    # ---- (2) are the two programmes actually distinct? -------------------
    rho, p = stats.spearmanr(meta["NR_score"], meta["R_score"])
    print(f"\n=== NR vs R programme correlation (segment level) ===", flush=True)
    print(f"  Spearman rho = {rho:.3f}  p = {p:.3g}  n = {len(meta)}", flush=True)
    print("  -> if rho is high, the two 'programmes' are not distinguishable\n"
          "     and DSP module scores are dominated by a shared axis.", flush=True)

    # ---- (3) partial correlation controlling for cell-type context -------
    print("=== Spearman of programme vs SERPINB2 (raw and partial) ===", flush=True)
    rows = []
    for prog in ["NR_score", "R_score"]:
        for serpin in ["SERPINB2", "SERPIN_family"]:
            for ctrl_name, ctrls in [
                    ("none", []),
                    ("epithelial", ["ctx_epithelial"]),
                    ("CAF", ["ctx_CAF"]),
                    ("myeloid", ["ctx_myeloid"]),
                    ("epi+CAF", ["ctx_epithelial", "ctx_CAF"]),
                    ("epi+CAF+myeloid", ["ctx_epithelial", "ctx_CAF", "ctx_myeloid"])]:
                x = meta[prog].values
                y = meta[serpin].values
                if ctrls:
                    zc = meta[ctrls].mean(axis=1).values
                    r, pp, n = partial_spearman(x, y, zc)
                    rows.append(dict(programme=prog, serpin=serpin,
                                     controlled_for=ctrl_name,
                                     rho=round(float(r), 3) if np.isfinite(r) else np.nan,
                                     p=round(float(pp), 4) if np.isfinite(pp) else np.nan,
                                     n=n))
                else:
                    m = np.isfinite(x) & np.isfinite(y)
                    r, pp = stats.spearmanr(x[m], y[m])
                    rows.append(dict(programme=prog, serpin=serpin,
                                     controlled_for=ctrl_name,
                                     rho=round(float(r), 3),
                                     p=round(float(pp), 4), n=int(m.sum())))
    pdf = pd.DataFrame(rows)
    pdf.to_csv(os.path.join(PROC, "bridge_partial_correlation.tsv"),
               sep="\t", index=False)
    print(pdf.to_string(index=False), flush=True)

    # ---- (1) ROI-level aggregation to remove pseudo-replication ----------
    print("\n=== ROI-level aggregation (one value per physical ROI) ===",
          flush=True)
    agg = meta.groupby("phys_roi").agg(
        NR=("NR_score", "mean"), R=("R_score", "mean"),
        SERPINB2=("SERPINB2", "mean"), SERPIN_family=("SERPIN_family", "mean"),
        epi=("ctx_epithelial", "mean"), caf=("ctx_CAF", "mean"),
        mye=("ctx_myeloid", "mean"), status=("status", "first"),
        patient=("Patient", "first"), n_seg=("Segment", "size"))
    print(f"  physical ROIs: {len(agg)}  (from {len(meta)} segments)", flush=True)

    roi_rows = []
    for serpin in ["SERPINB2", "SERPIN_family"]:
        thr = float(agg[serpin].median())
        for prog in ["NR", "R"]:
            for stratum, sub in [("all", agg),
                                 ("naive", agg[agg["status"] == "naive"]),
                                 ("treated", agg[agg["status"] == "treated"])]:
                hi = sub.loc[sub[serpin] > thr, prog].values
                lo = sub.loc[sub[serpin] <= thr, prog].values
                if len(hi) < 5 or len(lo) < 5:
                    continue
                u, pp = stats.mannwhitneyu(hi, lo, alternative="two-sided")
                roi_rows.append(dict(serpin=serpin, programme=prog, stratum=stratum,
                                     n_high=len(hi), n_low=len(lo),
                                     med_high=round(float(np.median(hi)), 4),
                                     med_low=round(float(np.median(lo)), 4),
                                     delta=round(float(np.median(hi) - np.median(lo)), 4),
                                     cliffs_delta=round(float(cliffs_delta(hi, lo)), 3),
                                     p=round(float(pp), 4)))
    rdf = pd.DataFrame(roi_rows)
    rdf.to_csv(os.path.join(PROC, "bridge_roi_level.tsv"), sep="\t", index=False)
    print(rdf.to_string(index=False), flush=True)

    # ---- segment-level vs ROI-level p, side by side ----------------------
    print("\n=== effect of pseudo-replication on significance ===", flush=True)
    print("  segment level (n=570) vs ROI level (n~235), NR programme:", flush=True)
    for serpin in ["SERPINB2", "SERPIN_family"]:
        thr_s = float(meta[serpin].median())
        hi = meta.loc[meta[serpin] > thr_s, "NR_score"].values
        lo = meta.loc[meta[serpin] <= thr_s, "NR_score"].values
        _, p_seg = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        thr_r = float(agg[serpin].median())
        hi2 = agg.loc[agg[serpin] > thr_r, "NR"].values
        lo2 = agg.loc[agg[serpin] <= thr_r, "NR"].values
        _, p_roi = stats.mannwhitneyu(hi2, lo2, alternative="two-sided")
        print(f"    {serpin:14s}: p_segment={p_seg:.5f}  p_ROI={p_roi:.5f}", flush=True)

    # ---- does myeloid context explain the NR score? ---------------------
    print("\n=== Does the NR programme reflect myeloid content? ===", flush=True)
    for prog in ["NR", "R"]:
        r, pp = stats.spearmanr(agg[prog], agg["mye"])
        print(f"  {prog}_score vs myeloid context: rho={r:.3f} p={pp:.3g}", flush=True)
    r, pp = stats.spearmanr(agg["NR"], agg["R"])
    print(f"  NR vs R at ROI level: rho={r:.3f} p={pp:.3g}", flush=True)

    agg.to_csv(os.path.join(PROC, "bridge_roi_level_scores.tsv"), sep="\t")
    print("\n[saved] bridge_partial_correlation.tsv, bridge_roi_level.tsv, "
          "bridge_roi_level_scores.tsv", flush=True)


if __name__ == "__main__":
    main()
