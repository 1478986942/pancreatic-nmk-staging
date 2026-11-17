"""GSE300595 step 6 - bridge the myeloid ATAC programmes to the spatial SERPIN niche.

QUESTION
    Are the non-responder (NR) enriched myeloid chromatin-accessibility programmes
    enriched in the SERPIN-high microenvironment, as measured independently in the
    GSE199102 GeoMx DSP cohort?

DESIGN / HONESTY CONSTRAINTS
    * GSE300595 (12 patients, all post-treatment) and GSE199102 (570 ROIs) share
      ZERO patients.  This is therefore a cross-cohort, cross-platform CONSISTENCY
      analysis, not a paired one, and it cannot support causal claims.
    * The three platforms do NOT measure the same SERPIN genes:
          GSE202051 snRNA : SERPINE1 + SERPINB2   (SERPINB2 ~0 -> effectively SERPINE1)
          GSE300595 snMultiome : SERPINE1 only     (SERPINB2 in 11/35115 nuclei)
          GSE199102 DSP   : SERPINB2 only          (SERPINE1 is NOT on the panel)
      The DSP SERPIN axis is reported as-is (SERPINB2), with a SERPIN-family
      mean as a sensitivity analysis.  It is NOT silently equated with SERPINE1.
    * DSP myeloid coverage is partial: MARCO, SPP1, MRC1, GPNMB, TREM2, LPL and
      PPARG are absent, so the LAM programme can only be approximated.
    * No peak survives FDR in the ATAC contrast (n=5 vs 7).  Programmes are
      rank/effect-size based and explicitly exploratory; the DSP bridge is the
      independent readout, and every test is reported with its effect size.
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import BASE, PROC

DSP_DIR = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"
DSP_OUT = r"E:\pancreatic_cancer\submit\gse202051_hum\processed_dsp"

TRT_MAP = {"Untreated": "naive", "CRT": "treated",
           "CRTL": "treated", "CRTLN": "treated"}

N_NULL = 2000
RNG = np.random.default_rng(0)


def cliffs_delta(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    gt = sum((x > y) for x in a for y in b)
    lt = sum((x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def mwu_block(df, score_col, group_col, hi_label, lo_label, tag):
    hi = df.loc[df[group_col] == hi_label, score_col].dropna().values
    lo = df.loc[df[group_col] == lo_label, score_col].dropna().values
    if len(hi) < 3 or len(lo) < 3:
        return None
    u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
    return dict(comparison=tag, score=score_col,
                n_high=len(hi), n_low=len(lo),
                med_high=round(float(np.median(hi)), 4),
                med_low=round(float(np.median(lo)), 4),
                delta_med=round(float(np.median(hi) - np.median(lo)), 4),
                cliffs_delta=round(float(cliffs_delta(hi, lo)), 3),
                p=round(float(p), 4))


def main():
    # ------------------------------------------------------------------ DSP load
    cnt = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "q3norm_counts.txt.gz"), "rt"),
                      sep="\t", index_col=0)
    seg = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "segment_props.txt.gz"), "rt"),
                      sep="\t")
    print(f"DSP matrix: {cnt.shape[0]} genes x {cnt.shape[1]} ROIs", flush=True)

    seg["ROI_id"] = seg["Sample_ID"].str.replace("-", ".", regex=False)
    common = [r for r in cnt.columns if r in set(seg["ROI_id"])]
    print(f"Matched ROIs: {len(common)}", flush=True)

    meta = seg.set_index("ROI_id").loc[common]
    X = cnt.loc[:, common].T                                  # ROIs x genes
    X = np.log2(X + 1.0)
    # z-score each gene across ROIs
    mu = X.mean(axis=0)
    sd = X.std(axis=0).replace(0, np.nan)
    Z = (X - mu) / sd
    Z = Z.fillna(0.0)
    panel = list(Z.columns)

    meta["status"] = meta["Treatment"].map(TRT_MAP).fillna("other")
    print("Treatment x status:", flush=True)
    print(meta.groupby(["status", "Treatment"], observed=True).size().to_string(),
          flush=True)
    print("\nSegment values:", meta["Segment"].value_counts().to_dict(), flush=True)

    # ------------------------------------------------------- SERPIN definitions
    # PRIMARY = SERPINB2 only.  SERPINE1 is absent from the panel, so this is the
    # only gene that the original DSP "SERPIN" score actually used.
    serpinb2 = Z["SERPINB2"]
    thr_b2 = float(serpinb2.median())
    meta["SERPIN_status"] = np.where(serpinb2 > thr_b2, "SERPIN_high", "SERPIN_low")
    print(f"\nPRIMARY SERPIN axis = SERPINB2  (median split at {thr_b2:.3f})", flush=True)
    print(meta.groupby(["status", "SERPIN_status"], observed=True).size().to_string(),
          flush=True)

    # SENSITIVITY = mean z of all panel SERPIN genes
    fam = [g for g in ["SERPINA1", "SERPINA3", "SERPINA5", "SERPINB1", "SERPINB2",
                       "SERPINB5", "SERPINB6", "SERPINB8", "SERPINF1", "SERPING1",
                       "SERPINH1"] if g in panel]
    meta["SERPIN_family"] = Z[fam].mean(axis=1)
    thr_fam = float(meta["SERPIN_family"].median())
    meta["SERPIN_family_status"] = np.where(meta["SERPIN_family"] > thr_fam,
                                            "SERPIN_high", "SERPIN_low")
    print(f"\nSENSITIVITY SERPIN axis = mean z of {fam}", flush=True)

    # ------------------------------------------------------------- programmes
    def load_prog(tag):
        p = os.path.join(PROC, f"DAR_pcgenes_{tag}.txt")
        return [l.strip() for l in open(p) if l.strip()]

    dar = pd.read_csv(os.path.join(PROC, "myeloid_DAR_annotated.tsv"), sep="\t")
    near = dar[(dar["dist_to_tss_bp"] <= 100_000) &
               (dar["log2FC"] < -0.5) & (dar["p"] < 0.05)]["nearest_pc_gene"]
    nr_genes_all = load_prog("NONRESPONDER")
    r_genes_all = load_prog("RESPONDER")
    nr_genes_near = [g for g in near.dropna().unique() if g]

    programmes = {}
    for name, genes in [("NR_all", nr_genes_all),
                        ("NR_near100kb", nr_genes_near),
                        ("R_all", r_genes_all)]:
        present = [g for g in genes if g in panel]
        programmes[name] = present
        print(f"programme {name:14s}: {len(genes):4d} genes -> {len(present):4d} on panel",
              flush=True)

    # LAM proxy from genes that DO exist on the panel
    lam_proxy = [g for g in ["APOE", "APOC1", "CLIP4", "CD163", "C1QA", "C1QB",
                             "LYZ", "CSF1R", "CTSD", "FABP5"] if g in panel]
    print(f"LAM proxy genes on panel: {lam_proxy}", flush=True)

    for name, genes in list(programmes.items()) + [("LAM_proxy", lam_proxy)]:
        if len(genes) >= 5:
            meta[name + "_score"] = Z[genes].mean(axis=1)

    # ------------------------------------------------- random-gene-set null
    # How unusual is a programme mean?  Compare against random gene sets of the
    # same size drawn from panel genes with comparable detection.
    sd_genes = [g for g in panel if Z[g].std() > 0]
    null = {}
    for name in ["NR_all", "NR_near100kb", "R_all"]:
        k = len(programmes[name])
        if k < 5:
            continue
        draws = []
        for _ in range(N_NULL):
            pick = RNG.choice(sd_genes, size=k, replace=False)
            draws.append(float(Z[pick].mean(axis=1).median()))
        null[name] = np.asarray(draws)
        print(f"null built for {name} (k={k}, {N_NULL} draws)", flush=True)

    # ------------------------------------------------------ DSP association
    rows = []
    for base_status in ["SERPIN_status", "SERPIN_family_status"]:
        label = "SERPINB2" if base_status == "SERPIN_status" else "SERPIN_family"
        for score in ["NR_all_score", "NR_near100kb_score", "R_all_score",
                      "LAM_proxy_score"]:
            if score not in meta.columns:
                continue
            for stratum, sub in [("all", meta),
                                 ("naive", meta[meta["status"] == "naive"]),
                                 ("treated", meta[meta["status"] == "treated"])]:
                r = mwu_block(sub, score, base_status, "SERPIN_high", "SERPIN_low",
                              f"{label} | {stratum}")
                if r:
                    r["serpin_axis"] = label
                    r["stratum"] = stratum
                    # empirical p against the random-gene-set null
                    if score.endswith("_score") and score[:-6] in null:
                        obs = r["delta_med"]
                        draws = null[score[:-6]]
                        r["null_empirical_p"] = round(
                            float((np.abs(draws) >= abs(obs)).mean()), 4)
                    rows.append(r)

    res = pd.DataFrame(rows)
    cols = ["serpin_axis", "stratum", "score", "n_high", "n_low",
            "med_high", "med_low", "delta_med", "cliffs_delta", "p",
            "null_empirical_p"]
    res = res[[c for c in cols if c in res.columns]]
    res.to_csv(os.path.join(PROC, "bridge_dsp_serpin_association.tsv"),
               sep="\t", index=False)

    pd.set_option("display.width", 220)
    print("\n=== DSP: programme score by SERPIN status ===", flush=True)
    print(res.to_string(index=False), flush=True)

    # --------------------------------------- continuous correlation in DSP
    corr_rows = []
    for score in ["NR_all_score", "NR_near100kb_score", "R_all_score",
                  "LAM_proxy_score"]:
        if score not in meta.columns:
            continue
        for stratum, sub in [("all", meta),
                             ("naive", meta[meta["status"] == "naive"]),
                             ("treated", meta[meta["status"] == "treated"])]:
            for sgene in ["SERPINB2", "SERPIN_family"]:
                if sgene == "SERPIN_family":
                    v = sub["SERPIN_family"]
                else:
                    v = Z.loc[sub.index, "SERPINB2"]
                m = sub[score].notna() & v.notna()
                if m.sum() < 10:
                    continue
                rho, p = stats.spearmanr(sub.loc[m, score], v[m])
                corr_rows.append(dict(stratum=stratum, score=score, against=sgene,
                                      n=int(m.sum()), rho=round(float(rho), 3),
                                      p=round(float(p), 4)))
    cdf = pd.DataFrame(corr_rows)
    cdf.to_csv(os.path.join(PROC, "bridge_dsp_continuous_corr.tsv"),
               sep="\t", index=False)
    print("\n=== DSP: continuous Spearman (programme vs SERPIN) ===", flush=True)
    print(cdf.to_string(index=False), flush=True)

    # -------------------------------------------- stratified by Segment
    seg_rows = []
    for segname in meta["Segment"].dropna().unique():
        sub = meta[meta["Segment"] == segname]
        if len(sub) < 30:
            continue
        for score in ["NR_all_score", "NR_near100kb_score", "R_all_score"]:
            if score not in meta.columns:
                continue
            r = mwu_block(sub, score, "SERPIN_status", "SERPIN_high", "SERPIN_low",
                          f"segment={segname}")
            if r:
                r["segment"] = segname
                seg_rows.append(r)
    sdf = pd.DataFrame(seg_rows)
    if not sdf.empty:
        sdf = sdf[["segment", "score", "n_high", "n_low", "med_high", "med_low",
                   "delta_med", "cliffs_delta", "p"]]
        sdf.to_csv(os.path.join(PROC, "bridge_dsp_by_segment.tsv"),
                   sep="\t", index=False)
        print("\n=== DSP: stratified by geometric Segment ===", flush=True)
        print(sdf.to_string(index=False), flush=True)

    # ------------------------------------------- patient-level ATAC vs SERPIN
    print("\n=== Patient level: ATAC programme vs snRNA SERPIN fraction (n=12) ===",
          flush=True)
    pb = pd.read_csv(os.path.join(PROC, "myeloid_atac_pseudobulk_counts.tsv"),
                     sep="\t", index_col=0)
    print(f"  ATAC pseudobulk: {pb.shape}", flush=True)

    pat = pd.read_csv(os.path.join(PROC, "patient_serpin_myeloid.tsv"), sep="\t")
    pat = pat.set_index("sample")

    # CTRL peaks = the R-enriched set, used only to normalise library composition
    dar_sets = {
        "NR": dar[(dar["log2FC"] < -0.5) & (dar["p"] < 0.05)]["peak"].tolist(),
        "R": dar[(dar["log2FC"] > 0.5) & (dar["p"] < 0.05)]["peak"].tolist(),
    }
    atac_rows = []
    for s in pb.index:
        d = dict(sample=s)
        for tag, peaks in dar_sets.items():
            present = [p for p in peaks if p in pb.columns]
            if present:
                d[f"{tag}_cpm"] = float(np.log1p(pb.loc[s, present]).mean())
        atac_rows.append(d)
    atac = pd.DataFrame(atac_rows).set_index("sample")
    atac["NR_minus_R"] = atac["NR_cpm"] - atac["R_cpm"]
    atac = atac.join(pat[["response", "SERPIN_frac", "mean_SERPIN_score",
                          "myeloid_frac"]])
    atac.to_csv(os.path.join(PROC, "bridge_patient_atac_vs_serpin.tsv"),
                sep="\t")
    print(atac.round(4).to_string(), flush=True)

    pt_rows = []
    for xcol in ["NR_cpm", "R_cpm", "NR_minus_R"]:
        for ycol in ["SERPIN_frac", "mean_SERPIN_score", "myeloid_frac"]:
            m = atac[xcol].notna() & atac[ycol].notna()
            if m.sum() < 6:
                continue
            rho, p = stats.spearmanr(atac.loc[m, xcol], atac.loc[m, ycol])
            pt_rows.append(dict(atac=xcol, serpin=ycol, n=int(m.sum()),
                                rho=round(float(rho), 3), p=round(float(p), 4)))
    ptdf = pd.DataFrame(pt_rows)
    ptdf.to_csv(os.path.join(PROC, "bridge_patient_spearman.tsv"),
                sep="\t", index=False)
    print("\n  Spearman (patient level):", flush=True)
    print(ptdf.to_string(index=False), flush=True)

    # responder vs non-responder sanity check on the ATAC scores
    print("\n  R vs NR on ATAC programme scores (sanity):", flush=True)
    for xcol in ["NR_cpm", "R_cpm", "NR_minus_R"]:
        r = atac.loc[atac["response"] == "responder", xcol].dropna()
        n = atac.loc[atac["response"] == "nonresponder", xcol].dropna()
        if len(r) >= 2 and len(n) >= 2:
            u, p = stats.mannwhitneyu(r, n, alternative="two-sided")
            print(f"    {xcol:12s} R={r.median():.3f} NR={n.median():.3f} "
                  f"p={p:.4f}", flush=True)

    print("\n[saved] bridge_*.tsv in processed/", flush=True)


if __name__ == "__main__":
    main()
