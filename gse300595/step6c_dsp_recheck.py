"""GSE300595 step 6c - re-test the earlier GSE199102 claims at the correct unit.

The original GSE199102 analysis treated 570 GeoMx segments as independent ROIs.
They are not: 194 physical ROIs were split into 570 geometric segments
(Immune / Epithelial / CAF).  This script re-runs the two claims that the
project narrative leans on, at segment level and at physical-ROI level:

    C1  SERPIN-high regions have MORE T-cell/NK signal (naive)
    C2  SERPIN-high regions have LESS myeloid signal (naive)

and adds a cell-type-context partial correlation, because a segment's signal is
dominated by which compartment it samples.

Panel limits (unchanged): T-cell score = CD3D/CD3E/NKG7 (no CD8A, GNLY);
myeloid score = CD163/CD68 (no MARCO, SPP1, MRC1); SERPIN = SERPINB2 only
(SERPINE1 absent).
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

TCELL = ["CD3D", "CD3E", "NKG7"]
MYE = ["CD163", "CD68"]
EPI = ["KRT8", "KRT18", "MUC1", "MSLN", "CLDN4", "LGALS4", "TSPAN8", "KRT17"]
CAF = ["COL1A1", "COL1A2", "DCN", "LUM", "FAP", "FBLN1", "PDGFRB",
       "ACTA2", "POSTN", "COL5A1"]


def cliffs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def partial_spear(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    n = len(x)
    if n < 8:
        return np.nan, np.nan, n
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    def resid(a):
        return a - np.polyval(np.polyfit(rz, a, 1), rz)
    r = np.corrcoef(resid(rx), resid(ry))[0, 1]
    df = n - 3
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    return r, 2 * stats.t.sf(abs(t), df), n


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

    # scores as in the original script: mean of raw q3norm values (not z)
    raw = cnt.loc[:, common].T
    meta["TCELL"] = raw[[g for g in TCELL if g in panel]].mean(axis=1)
    meta["MYE"] = raw[[g for g in MYE if g in panel]].mean(axis=1)
    meta["SERPINB2"] = raw["SERPINB2"]
    meta["SERPINB2_z"] = Z["SERPINB2"]
    meta["epi"] = Z[[g for g in EPI if g in panel]].mean(axis=1)
    meta["caf"] = Z[[g for g in CAF if g in panel]].mean(axis=1)

    out = []

    # ---------------- C1/C2 at SEGMENT level (reproduce original) ----------
    print("=== SEGMENT level (n=570, as originally analysed) ===", flush=True)
    for score in ["TCELL", "MYE"]:
        for stratum in ["naive", "treated"]:
            sub = meta[meta["status"] == stratum]
            thr = float(meta["SERPINB2"].median())
            hi = sub.loc[sub["SERPINB2"] > thr, score].values
            lo = sub.loc[sub["SERPINB2"] <= thr, score].values
            u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
            print(f"  {score:6s} {stratum:8s}: hi={np.median(hi):7.3f} "
                  f"lo={np.median(lo):7.3f} delta={np.median(hi)-np.median(lo):+6.3f} "
                  f"p={p:.4f}", flush=True)
            out.append(dict(level="segment", score=score, stratum=stratum,
                            n_high=len(hi), n_low=len(lo),
                            med_high=round(float(np.median(hi)), 3),
                            med_low=round(float(np.median(lo)), 3),
                            delta=round(float(np.median(hi) - np.median(lo)), 3),
                            cliffs_delta=round(float(cliffs(hi, lo)), 3),
                            p=round(float(p), 4)))

    # ---------------- C1/C2 at PHYSICAL ROI level --------------------------
    print("\n=== PHYSICAL ROI level (n=194, no pseudo-replication) ===", flush=True)
    agg = meta.groupby("phys_roi").agg(
        TCELL=("TCELL", "mean"), MYE=("MYE", "mean"),
        SERPINB2=("SERPINB2", "mean"), SERPINB2_z=("SERPINB2_z", "mean"),
        epi=("epi", "mean"), caf=("caf", "mean"),
        status=("status", "first"), patient=("Patient", "first"))
    print(f"  physical ROIs: {len(agg)}", flush=True)
    thr_r = float(agg["SERPINB2"].median())
    for score in ["TCELL", "MYE"]:
        for stratum in ["naive", "treated"]:
            sub = agg[agg["status"] == stratum]
            hi = sub.loc[sub["SERPINB2"] > thr_r, score].values
            lo = sub.loc[sub["SERPINB2"] <= thr_r, score].values
            if len(hi) < 5 or len(lo) < 5:
                print(f"  {score:6s} {stratum:8s}: insufficient", flush=True)
                continue
            u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
            print(f"  {score:6s} {stratum:8s}: hi={np.median(hi):7.3f} "
                  f"lo={np.median(lo):7.3f} delta={np.median(hi)-np.median(lo):+6.3f} "
                  f"p={p:.4f}  (n_hi={len(hi)}, n_lo={len(lo)})", flush=True)
            out.append(dict(level="physical_ROI", score=score, stratum=stratum,
                            n_high=len(hi), n_low=len(lo),
                            med_high=round(float(np.median(hi)), 3),
                            med_low=round(float(np.median(lo)), 3),
                            delta=round(float(np.median(hi) - np.median(lo)), 3),
                            cliffs_delta=round(float(cliffs(hi, lo)), 3),
                            p=round(float(p), 4)))

    # ---------------- partial correlation at ROI level ---------------------
    print("\n=== ROI level: partial Spearman vs SERPINB2 ===", flush=True)
    for score in ["TCELL", "MYE"]:
        for ctrl, zc in [("none", None), ("epithelial", agg["epi"].values),
                         ("CAF", agg["caf"].values),
                         ("epi+CAF", (agg["epi"] + agg["caf"]).values)]:
            if zc is None:
                m = np.isfinite(agg[score]) & np.isfinite(agg["SERPINB2_z"])
                r, p = stats.spearmanr(agg.loc[m, score], agg.loc[m, "SERPINB2_z"])
                n = int(m.sum())
            else:
                r, p, n = partial_spear(agg[score].values,
                                        agg["SERPINB2_z"].values, zc)
            print(f"  {score:6s} controlled_for={ctrl:11s}: rho={r:+.3f} p={p:.4f} n={n}",
                  flush=True)
            out.append(dict(level="physical_ROI_partial", score=score,
                            stratum=f"ctrl={ctrl}", n_high=n, n_low=np.nan,
                            med_high=round(float(r), 3), med_low=np.nan,
                            delta=np.nan, cliffs_delta=np.nan,
                            p=round(float(p), 4)))

    # ---------------- within-compartment test ------------------------------
    print("\n=== Within Segment compartment (controls sampling bias) ===", flush=True)
    for segname in ["Immune", "Epithelial", "CAF"]:
        sub = meta[meta["Segment"] == segname]
        for score in ["TCELL", "MYE"]:
            hi = sub.loc[sub["SERPINB2"] > float(meta["SERPINB2"].median()),
                         score].values
            lo = sub.loc[sub["SERPINB2"] <= float(meta["SERPINB2"].median()),
                         score].values
            if len(hi) < 5 or len(lo) < 5:
                continue
            u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
            print(f"  {segname:11s} {score:6s}: hi={np.median(hi):7.3f} "
                  f"lo={np.median(lo):7.3f} delta={np.median(hi)-np.median(lo):+6.3f} "
                  f"p={p:.4f} (n_hi={len(hi)}, n_lo={len(lo)})", flush=True)
            out.append(dict(level=f"segment_within_{segname}", score=score,
                            stratum="all", n_high=len(hi), n_low=len(lo),
                            med_high=round(float(np.median(hi)), 3),
                            med_low=round(float(np.median(lo)), 3),
                            delta=round(float(np.median(hi) - np.median(lo)), 3),
                            cliffs_delta=round(float(cliffs(hi, lo)), 3),
                            p=round(float(p), 4)))

    df = pd.DataFrame(out)
    df.to_csv(os.path.join(PROC, "bridge_dsp_claims_recheck.tsv"), sep="\t",
              index=False)
    print("\n[saved] bridge_dsp_claims_recheck.tsv", flush=True)


if __name__ == "__main__":
    main()
