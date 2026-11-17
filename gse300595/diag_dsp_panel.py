import os, gzip
import pandas as pd

DSP = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"
cnt = pd.read_csv(gzip.open(os.path.join(DSP, "q3norm_counts.txt.gz"), "rt"),
                  sep="\t", index_col=0)
panel = set(cnt.index)

probe = ["SERPINE1", "SERPINB2", "SERPINA1", "SERPINA3", "SERPINB1",
         "SERPINF1", "SERPINH1", "MARCO", "SPP1", "CD163", "CD68", "MRC1",
         "CD3D", "CD3E", "CD8A", "GNLY", "NKG7", "PECAM1", "VWF", "CD31",
         "APOE", "APOC1", "GPNMB", "TREM2", "LPL", "PPARG", "CLIP4",
         "C1QA", "C1QB", "HLA-DRA", "CD74", "FLT3", "CD1C", "LYZ", "CSF1R",
         "PTPRC", "COL1A1", "EPCAM", "KRT19"]

print(f"{'gene':10s} {'in DSP panel':>14s}")
for g in probe:
    print(f"{g:10s} {'YES' if g in panel else '-- NO --':>14s}")

# also check for any SERPINE alias
print("\nAny name starting with SERPINE:", sorted([g for g in panel if g.startswith("SERPINE")]))
print("Any name containing MARCO:", sorted([g for g in panel if "MARCO" in g]))
