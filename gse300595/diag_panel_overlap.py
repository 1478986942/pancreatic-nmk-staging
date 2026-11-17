import os, gzip
import numpy as np
import pandas as pd

DSP = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"
PROC = r"E:\pancreatic_cancer\submit\gse300595\processed"

cnt = pd.read_csv(gzip.open(os.path.join(DSP, "q3norm_counts.txt.gz"), "rt"),
                  sep="\t", index_col=0)
print("DSP matrix:", cnt.shape)
panel = list(cnt.index)
print("first 30 panel genes:", panel[:30])

for tag in ["NONRESPONDER", "RESPONDER", "NONRESPONDER_top200", "RESPONDER_top200"]:
    p = os.path.join(PROC, f"DAR_pcgenes_{tag}.txt")
    if not os.path.exists(p):
        print(f"{tag}: MISSING"); continue
    genes = [l.strip() for l in open(p) if l.strip()]
    hit = [g for g in genes if g in panel]
    print(f"{tag:22s}: n={len(genes):4d}  overlap_with_DSP={len(hit):3d}  {hit[:40]}")

print("\nPanel genes containing SERP/MARCO/SPP1/CD:",
      [g for g in panel if any(k in g for k in ("SERP", "MARCO", "SPP1", "CD16", "CD68", "CD163"))])

seg = pd.read_csv(gzip.open(os.path.join(DSP, "segment_props.txt.gz"), "rt"), sep="\t")
print("\nsegment_props columns:", list(seg.columns))
print(seg.head(3).to_string())
print("\nTreatment values:", seg["Treatment"].value_counts().to_dict())
print("Histology values:", seg["Histology"].value_counts(dropna=False).head(10).to_dict())
