"""Quantify the pseudo-replication defect (CORRECTIONS.md C1) as a methods result.

Purpose
-------
The correction file reports one before/after pair (segment-level p=0.0063 vs
physical-ROI-level p=0.685).  A methods paper needs the mechanism quantified, not
just exhibited.  This script computes:

  1. Intraclass correlation (ICC) of each score across GeoMx segments nested in
     physical ROIs -> how much variance is between-ROI vs within-ROI.
  2. Design effect (1 + (m-1)*ICC) and effective sample size.
  3. A null simulation: assign SERPIN status at the PHYSICAL ROI level, carry it
     to segments, then test at segment level.  If segments were independent the
     type I error rate would be 5 %.  The observed rate is the inflation.

Null simulation is the key number: it separates "our p-value moved" from
"segment-level testing is systematically anti-conservative".
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import stats

DSP_DIR = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"
OUT = r"E:\pancreatic_cancer\submit\methods_paper"
os.makedirs(OUT, exist_ok=True)

TCELL = ["CD3D", "CD3E", "NKG7"]
MYE = ["CD163", "CD68"]
TRT_MAP = {"Untreated": "naive", "CRT": "treated",
           "CRTL": "treated", "CRTLN": "treated"}

B = 2000
RNG = np.random.default_rng(20261114)


def icc_oneway(values, groups):
    """One-way random-effects ICC(1) with unbalanced-cluster correction."""
    df = pd.DataFrame({"y": values, "g": groups}).dropna()
    y = df["y"].values
    g = df["g"].values
    k = df["g"].nunique()
    N = len(df)
    if k < 2 or N - k < 1:
        return np.nan, np.nan, np.nan, np.nan
    grand = y.mean()
    grp = df.groupby("g")["y"]
    n_g = grp.size().values
    m_g = grp.mean().values
    ss_between = float(np.sum(n_g * (m_g - grand) ** 2))
    ss_within = float(np.sum((y - grp.transform("mean").values) ** 2))
    ms_b = ss_between / (k - 1)
    ms_w = ss_within / (N - k)
    m0 = (N - np.sum(n_g ** 2) / N) / (k - 1)
    denom = ms_b + (m0 - 1) * ms_w
    icc = (ms_b - ms_w) / denom if denom > 0 else np.nan
    icc = max(0.0, min(1.0, icc))
    deff = 1.0 + (m0 - 1) * icc
    n_eff = N / deff
    return icc, deff, n_eff, m0


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

    raw = cnt.loc[:, common].T
    scores = {
        "T-cell (CD3D/CD3E/NKG7)": raw[[g for g in TCELL if g in panel]].mean(axis=1),
        "Myeloid (CD163/CD68)": raw[[g for g in MYE if g in panel]].mean(axis=1),
        "SERPINB2": raw["SERPINB2"],
    }
    # add the two accessibility programmes if present
    for tag in ["NONRESPONDER", "RESPONDER"]:
        p = os.path.join(r"E:\pancreatic_cancer\submit\gse300595\processed",
                         f"DAR_pcgenes_{tag}.txt")
        if os.path.exists(p):
            genes = [l.strip() for l in open(p) if l.strip()]
            have = [g for g in genes if g in panel]
            if len(have) >= 5:
                scores[f"{tag} accessibility programme"] = Z[have].mean(axis=1)

    print(f"segments={len(meta)}  physical ROIs={meta['phys_roi'].nunique()}", flush=True)
    cs = meta.groupby("phys_roi").size()
    print(f"cluster size: mean={cs.mean():.2f} median={cs.median():.0f} "
          f"min={cs.min()} max={cs.max()}", flush=True)
    print(f"design effect at ICC observed below applies to N={len(meta)}\n", flush=True)

    # ---- (1)+(2) ICC / design effect -------------------------------------
    rows = []
    print("=== ICC and design effect (segments nested in physical ROI) ===",
          flush=True)
    for name, v in scores.items():
        icc, deff, n_eff, m0 = icc_oneway(v.values, meta["phys_roi"].values)
        rows.append(dict(score=name, n_segments=len(meta),
                         n_physical_ROI=int(meta["phys_roi"].nunique()),
                         mean_cluster_size=round(float(m0), 3),
                         ICC=round(float(icc), 4),
                         design_effect=round(float(deff), 3),
                         effective_n=round(float(n_eff), 1),
                         n_inflation=round(len(meta) / n_eff, 2) if n_eff else np.nan))
        print(f"  {name:34s} ICC={icc:6.3f}  deff={deff:5.2f}  "
              f"n_eff={n_eff:6.1f}  (nominal n={len(meta)}, "
              f"inflation x{len(meta)/n_eff:.2f})", flush=True)
    icc_df = pd.DataFrame(rows)
    icc_df.to_csv(os.path.join(OUT, "design_effect_icc.tsv"), sep="\t", index=False)

    # ---- (3) null simulation: type I error under ROI-level null ----------
    print("\n=== Type I error simulation (B=%d) ===" % B, flush=True)
    rois = meta["phys_roi"].values
    uniq = np.unique(rois)
    roi_index = {r: i for i, r in enumerate(uniq)}
    seg_roi_idx = np.array([roi_index[r] for r in rois])

    sim_rows = []
    for name, v in scores.items():
        v = v.values.astype(float)
        ok = np.isfinite(v)
        if ok.sum() < 50:
            continue
        roi_mean = np.array([v[ok & (seg_roi_idx == i)].mean() for i in range(len(uniq))])
        roi_n = np.array([int((ok & (seg_roi_idx == i)).sum()) for i in range(len(uniq))])
        valid_roi = np.where(roi_n >= 1)[0]

        n_high = len(valid_roi) // 2
        rej_seg = 0
        rej_roi = 0
        done = 0
        while done < B:
            perm = RNG.permutation(valid_roi)
            high_rois = set(perm[:n_high].tolist())
            is_high_seg = np.array([i in high_rois for i in seg_roi_idx])
            hs = v[ok & is_high_seg]
            ls = v[ok & ~is_high_seg]
            if len(hs) < 3 or len(ls) < 3:
                continue
            _, p_seg = stats.mannwhitneyu(hs, ls, alternative="two-sided")
            rh = [i for i in valid_roi if i in high_rois]
            rl = [i for i in valid_roi if i not in high_rois]
            _, p_roi = stats.mannwhitneyu(roi_mean[rh], roi_mean[rl],
                                          alternative="two-sided")
            if p_seg < 0.05:
                rej_seg += 1
            if p_roi < 0.05:
                rej_roi += 1
            done += 1

        sim_rows.append(dict(score=name, B=done,
                             type_I_segment_level=round(rej_seg / done, 4),
                             type_I_ROI_level=round(rej_roi / done, 4),
                             inflation_factor=round((rej_seg / done) /
                                                    max(1e-9, rej_roi / done), 2)))
        print(f"  {name:34s} alpha=0.05 -> segment={rej_seg/done:.3f}  "
              f"ROI={rej_roi/done:.3f}", flush=True)

    sim_df = pd.DataFrame(sim_rows)
    sim_df.to_csv(os.path.join(OUT, "type_I_simulation.tsv"), sep="\t", index=False)

    # ---- observed before/after for reference -----------------------------
    print("\n=== Observed segment vs ROI tests (SERPINB2 median split) ===",
          flush=True)
    obs_rows = []
    thr_seg = float(raw["SERPINB2"].median())
    agg = meta.groupby("phys_roi").agg(
        status=("status", "first"))
    for name, v in scores.items():
        v = v.values.astype(float)
        roi_mean = pd.Series(v, index=meta.index).groupby(
            meta["phys_roi"].values).mean()
        hi = v[raw["SERPINB2"].values > thr_seg]
        lo = v[raw["SERPINB2"].values <= thr_seg]
        _, p_seg = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        thr_roi = float(raw["SERPINB2"].groupby(meta["phys_roi"].values).mean().median())
        roi_serp = raw["SERPINB2"].groupby(meta["phys_roi"].values).mean()
        h2 = roi_mean.values[roi_serp.values > thr_roi]
        l2 = roi_mean.values[roi_serp.values <= thr_roi]
        _, p_roi = stats.mannwhitneyu(h2, l2, alternative="two-sided")
        obs_rows.append(dict(score=name, p_segment=round(float(p_seg), 5),
                             p_ROI=round(float(p_roi), 5),
                             p_ratio=round(float(p_roi / p_seg), 2) if p_seg > 0 else np.nan))
        print(f"  {name:34s} p_segment={p_seg:.5f}  p_ROI={p_roi:.5f}", flush=True)
    pd.DataFrame(obs_rows).to_csv(os.path.join(OUT, "observed_p_before_after.tsv"),
                                  sep="\t", index=False)

    print(f"\n[saved] {OUT}", flush=True)


if __name__ == "__main__":
    main()
