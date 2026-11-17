"""
GSE199102 GeoMx DSP spatial analysis - SERPIN niche proximity in human PDAC.
Tests whether SERPIN-high ROIs show higher MYELOID and lower TCELL signals
vs SERPIN-low ROIs, stratified by treatment (naive vs CRT).
All 570 matched ROIs are used (histology labels unavailable for this subset).
"""
import os, sys, gzip, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"E:\pancreatic_cancer\submit\serpin_repro\venv\Lib\site-packages")
import numpy as np
import pandas as pd
from scipy import stats

BASE  = r"E:\pancreatic_cancer\submit\gse202051_hum"
DATA  = os.path.join(BASE, "dsp_data")
OUT   = os.path.join(BASE, "processed_dsp")
for d in [OUT]:
    os.makedirs(d, exist_ok=True)

print("[STEP DSP] Loading...", flush=True)
t0 = time.time()

cnt  = pd.read_csv(gzip.open(os.path.join(DATA,"q3norm_counts.txt.gz"),"rt"), sep="\t", index_col=0)
seg  = pd.read_csv(gzip.open(os.path.join(DATA,"segment_props.txt.gz"),"rt"), sep="\t")
print(f"  Matrix {cnt.shape}  Segments {seg.shape[0]}  ({time.time()-t0:.1f}s)", flush=True)

# Match ROI IDs: count matrix uses dots, seg uses dashes
seg_roi = seg["Sample_ID"].str.replace("-",".").astype(str)
common  = [r for r in cnt.columns if r in set(seg_roi)]
mat     = cnt.loc[:, common].T
meta    = seg[["Sample_ID","Patient","Treatment","Histology",
               "X_coordinate","Y_coordinate","AOI_area","Nuclei"]].copy()
meta["ROI_id"] = meta["Sample_ID"].str.replace("-",".")
meta = meta.merge(pd.DataFrame({"ROI_id": common}), on="ROI_id", how="inner")

merged = meta[["ROI_id","Patient","Treatment","Histology",
               "X_coordinate","Y_coordinate","AOI_area","Nuclei"]].reset_index(drop=True)
m_ord  = mat.loc[meta["ROI_id"].values].reset_index(drop=True)
for c in m_ord.columns:
    merged[c] = m_ord[c].values
print(f"  Merged: {len(merged)} ROIs  ({time.time()-t0:.1f}s)", flush=True)

TRT_MAP = {"Untreated":"naive","CRT":"treated","CRTL":"treated","CRTLN":"treated"}
merged["status"] = merged["Treatment"].map(TRT_MAP).fillna("other")
print(f"  naive={ (merged['status']=='naive').sum() }  treated={ (merged['status']=='treated').sum() }", flush=True)

def pick(names): return [g for g in names if g in merged.columns]
SERPIN_G  = pick(["SERPINE1","SERPINB2"])
MYEO_G    = pick(["MARCO","CD163","CD68","SPP1"])
TCELL_G   = pick(["CD3D","CD3E","CD8A","GNLY","NKG7"])
ENDO_G    = pick(["PECAM1","VWF","CD31"])
MYEO_BR_G = pick(["CD14","CSF1R","CTSS","LST1"])

print("\n=== Panel coverage ===", flush=True)
for name, genes in [("SERPIN",SERPIN_G),("MYEO",MYEO_G),("TCELL",TCELL_G),
                     ("ENDO",ENDO_G),("MYEO_BR",MYEO_BR_G)]:
    print(f"  {name}: {genes or ['NONE']}", flush=True)

for genes, prefix in [(SERPIN_G,"SERPIN"),(MYEO_G,"MYELOID"),(TCELL_G,"TCELL"),
                       (ENDO_G,"ENDO"),(MYEO_BR_G,"MYELOID_BR")]:
    if not genes:
        merged[f"{prefix}_score"] = np.nan; continue
    merged[f"{prefix}_score"] = merged[genes].mean(axis=1)
    s = merged[f"{prefix}_score"]
    print(f"  {prefix:10s}: n={s.notna().sum():4d} med={s.median():.2f} mean={s.mean():.2f}", flush=True)

roi_all = merged[merged["SERPIN_score"].notna()].copy()
thr     = float(roi_all["SERPIN_score"].median())
roi_all["SERPIN_status"] = np.where(roi_all["SERPIN_score"] > thr, "SERPIN_high", "SERPIN_low")
print(f"\n=== ROIs analyzed: n={len(roi_all)}  thr={thr:.2f}  HIGH={(roi_all['SERPIN_status']=='SERPIN_high').sum()} LOW={(roi_all['SERPIN_status']=='SERPIN_low').sum()} ===", flush=True)

results = []
for st in ["naive","treated"]:
    sub = roi_all[roi_all["status"]==st]
    n_hi = int((sub["SERPIN_status"]=="SERPIN_high").sum())
    n_lo = int((sub["SERPIN_status"]=="SERPIN_low").sum())
    if n_hi < 3 or n_lo < 3:
        print(f"  {st}: skip (hi={n_hi}, lo={n_lo})", flush=True); continue
    hi, lo = sub[sub["SERPIN_status"]=="SERPIN_high"], sub[sub["SERPIN_status"]=="SERPIN_low"]
    for sc, label in [("MYELOID_score","Myeloid signal"),("TCELL_score","T-cell/NK"),
                      ("ENDO_score","Endothelial"),("MYELOID_BR_score","Broad myeloid")]:
        vh = hi[sc].dropna().values; vl = lo[sc].dropna().values
        if len(vh)<3 or len(vl)<3: continue
        stat, p = stats.mannwhitneyu(vh,vl,alternative="two-sided")
        results.append({"group":st,"n_HIGH":len(vh),"n_LO":len(vl),"score":label,
                         "med_HIGH":round(float(np.median(vh)),3),"med_LO":round(float(np.median(vl)),3),
                         "delta":round(float(np.median(vh)-np.median(vl)),3),
                         "U":round(float(stat),1),"p":round(float(p),4)})

rdf = pd.DataFrame(results)
print("\n=== SERPIN-high vs SERPIN-low comparison ===", flush=True)
if not rdf.empty:
    print(rdf.to_string(index=False), flush=True)
else:
    print("  No comparisons", flush=True)
rdf.to_csv(os.path.join(OUT,"serpin_proximity_comparison.tsv"),sep="\t",index=False)

print("\n=== SERPIN-high niche: naive vs treated ===", flush=True)
sh = roi_all[roi_all["SERPIN_status"]=="SERPIN_high"]
for sc, label in [("MYELOID_score","Myeloid"),("TCELL_score","T-cell"),
                  ("SERPIN_score","SERPIN"),("ENDO_score","Endothelial")]:
    nv = sh[sh["status"]=="naive"][sc].dropna()
    tv = sh[sh["status"]=="treated"][sc].dropna()
    if len(nv)>=3 and len(tv)>=3:
        stat,p = stats.mannwhitneyu(nv,tv,alternative="two-sided")
        print(f"  {label:15s} naive(n={len(nv):3d}) med={nv.median():.2f}  "
              f"trtd(n={len(tv):3d}) med={tv.median():.2f}  U={stat:.0f} p={p:.4f}", flush=True)
    else:
        print(f"  {label:15s} insufficient (n_naive={len(nv)}, n_trtd={len(tv)})", flush=True)

pt = roi_all.pivot_table(values=["SERPIN_score","MYELOID_score","TCELL_score"],
                          index=["status","SERPIN_status"], aggfunc="mean")
print("\n=== Mean scores by status x SERPIN_status ===", flush=True)
print(pt.round(3).to_string(), flush=True)

pat = roi_all.groupby(["Patient","status"]).agg(
    n_ROIs=("ROI_id","count"),
    n_SERPIN_high=("SERPIN_status", lambda x:(x=="SERPIN_high").sum()),
    SERPIN_frac=("SERPIN_status", lambda x:(x=="SERPIN_high").mean()),
    mean_MYEO=("MYELOID_score","mean"),
    mean_TCELL=("TCELL_score","mean"),
    mean_SERPIN=("SERPIN_score","mean"),
).reset_index()
pat.to_csv(os.path.join(OUT,"dsp_patient_level.tsv"),sep="\t",index=False)
print(f"\nPatient-level DSP: {len(pat)} patients", flush=True)
print(pat.to_string(index=False), flush=True)

summ = {"dataset":"GSE199102_GeoMx_DSP","total_ROIs":len(merged),
        "analyzed_ROIs":len(roi_all),
        "naive_ROIs":int((roi_all["status"]=="naive").sum()),
        "treated_ROIs":int((roi_all["status"]=="treated").sum()),
        "SERPIN_threshold":round(thr,4),
        "SERPIN_high_ROIs":int((roi_all["SERPIN_status"]=="SERPIN_high").sum()),
        "SERPIN_low_ROIs":int((roi_all["SERPIN_status"]=="SERPIN_low").sum()),
        "SERPIN_genes":str(SERPIN_G),"MYEO_genes":str(MYEO_G),"TCELL_genes":str(TCELL_G)}
pd.DataFrame([summ]).to_csv(os.path.join(OUT,"step_dsp_summary.tsv"),sep="\t",index=False)
print(f"\n[OK] step_dsp_summary.tsv", flush=True)
print(f"[OK] serpin_proximity_comparison.tsv", flush=True)
print(f"[OK] dsp_patient_level.tsv", flush=True)
print("\n=== DSP STEP COMPLETE ===", flush=True)
