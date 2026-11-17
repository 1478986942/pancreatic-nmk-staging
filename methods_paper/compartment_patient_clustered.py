"""Compartment-stratified, patient-clustered re-analysis of GSE199102.

Design (verified in check_design.py):
  - 235 physical ROIs, each split into <=1 segment per compartment
    (Immune / Epithelial / CAF); 0 duplicated (ROI, compartment) pairs.
  - Therefore WITHIN a compartment every observation is a distinct ROI:
    there is no ROI-level pseudo-replication inside a stratum.
  - What remains is PATIENT clustering (7-15 ROIs per patient, 22 patients).

Consequence: the earlier "segment-level" test was wrong for a different reason
than assumed.  It pooled three different tissue compartments into one contrast,
which is a confounded comparison, not a clustered one.  The correct test is
compartment-stratified with patient-level inference.

Two patient-aware tests are reported:
  A. Per-patient paired difference: for each patient having both SERPIN-high and
     SERPIN-low ROIs in a compartment, take mean(high) - mean(low); then a
     one-sample Wilcoxon signed-rank across patients.
  B. Within-patient label permutation: shuffle SERPIN status among a patient's
     ROIs, recompute the pooled difference, build a null over B iterations.
"""
import os
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
B = 5000
RNG = np.random.default_rng(7)


def main():
    cnt = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "q3norm_counts.txt.gz"), "rt"),
                      sep="\t", index_col=0)
    seg = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "segment_props.txt.gz"), "rt"),
                      sep="\t")
    seg["ROI_id"] = seg["Sample_ID"].str.replace("-", ".", regex=False)
    common = [r for r in cnt.columns if r in set(seg["ROI_id"])]
    d = seg.set_index("ROI_id").loc[common].copy()
    raw = cnt.loc[:, common].T

    d["status"] = d["Treatment"].map(TRT_MAP).fillna("other")
    d["TCELL"] = raw[[g for g in TCELL if g in raw.columns]].mean(axis=1)
    d["MYE"] = raw[[g for g in MYE if g in raw.columns]].mean(axis=1)
    d["SERPINB2"] = raw["SERPINB2"]

    rows = []
    for compartment in ["Epithelial", "CAF", "Immune"]:
        for split in ["global", "within_compartment"]:
            for stratum in ["all", "naive", "treated"]:
                sub = d[(d["Segment"] == compartment)].copy()
                if stratum != "all":
                    sub = sub[sub["status"] == stratum]
                if len(sub) < 20:
                    continue
                if split == "global":
                    thr = float(d["SERPINB2"].median())
                else:
                    thr = float(sub["SERPINB2"].median())
                sub["grp"] = np.where(sub["SERPINB2"] > thr, "high", "low")
                sub = sub.dropna(subset=["SERPINB2"])

                for score in ["TCELL", "MYE"]:
                    hi = sub.loc[sub["grp"] == "high", score].values
                    lo = sub.loc[sub["grp"] == "low", score].values
                    if len(hi) < 5 or len(lo) < 5:
                        continue
                    # pooled (ignores patient clustering) - for comparison only
                    _, p_pool = stats.mannwhitneyu(hi, lo, alternative="two-sided")
                    obs_diff = float(np.mean(hi) - np.mean(lo))

                    # A. per-patient paired differences
                    diffs = []
                    for pat, g in sub.groupby("Patient"):
                        gh = g.loc[g["grp"] == "high", score].dropna().values
                        gl = g.loc[g["grp"] == "low", score].dropna().values
                        if len(gh) >= 1 and len(gl) >= 1:
                            diffs.append(float(np.mean(gh) - np.mean(gl)))
                    diffs = np.asarray(diffs)
                    if len(diffs) >= 5:
                        try:
                            _, p_sign = stats.wilcoxon(diffs, alternative="two-sided")
                        except ValueError:
                            p_sign = np.nan
                        med_diff = float(np.median(diffs))
                    else:
                        p_sign, med_diff = np.nan, np.nan

                    # B. within-patient permutation of labels
                    groups = [g.index.values for _, g in sub.groupby("Patient")]
                    vals = sub[score].values
                    lab = (sub["grp"] == "high").values
                    idx_of = {ix: i for i, ix in enumerate(sub.index)}
                    null = np.empty(B)
                    for b in range(B):
                        perm = lab.copy()
                        for gidx in groups:
                            pos = np.array([idx_of[ix] for ix in gidx])
                            perm[pos] = RNG.permutation(perm[pos])
                        h = vals[perm]
                        l = vals[~perm]
                        null[b] = h.mean() - l.mean()
                    p_perm = float((np.abs(null) >= abs(obs_diff)).mean())

                    rows.append(dict(
                        compartment=compartment, serpin_split=split,
                        stratum=stratum, score=score,
                        n_ROI=len(sub), n_high=len(hi), n_low=len(lo),
                        n_patients_paired=len(diffs) if len(diffs) else 0,
                        mean_high=round(float(np.mean(hi)), 3),
                        mean_low=round(float(np.mean(lo)), 3),
                        obs_diff=round(obs_diff, 3),
                        median_patient_diff=round(med_diff, 3) if np.isfinite(med_diff) else np.nan,
                        p_pooled=round(float(p_pool), 5),
                        p_patient_signedrank=round(float(p_sign), 5) if np.isfinite(p_sign) else np.nan,
                        p_patient_permutation=round(p_perm, 5)))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "compartment_patient_clustered.tsv"),
               sep="\t", index=False)

    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 60)
    show = res[res["serpin_split"] == "within_compartment"]
    print("=== within-compartment SERPIN split, patient-aware inference ===",
          flush=True)
    print(show.to_string(index=False), flush=True)

    print("\n=== global SERPIN split, patient-aware inference ===", flush=True)
    print(res[res["serpin_split"] == "global"].to_string(index=False), flush=True)

    print("\n[note] p_pooled ignores patient clustering and is shown only so the "
          "inflation from ignoring it is visible.", flush=True)
    print(f"[saved] {os.path.join(OUT, 'compartment_patient_clustered.tsv')}",
          flush=True)


if __name__ == "__main__":
    main()
