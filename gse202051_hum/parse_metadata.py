"""Parse GSE202051 series matrix → clean sample_metadata.tsv."""
import gzip, re, os, pandas as pd

DATA = r"E:\pancreatic_cancer\submit\gse202051_hum\data"
txt_path = os.path.join(DATA, "series_matrix.txt")
text = open(txt_path, encoding="utf-8", errors="replace").read()
lines = [l for l in text.split("\n") if l.startswith("!Sample_")]

def split_soft(line):
    return [p.strip('"').strip() for p in line[1:].split("\t")]

gsms      = split_soft(next(l for l in lines if l.startswith("!Sample_geo_accession")))
titles    = split_soft(next(l for l in lines if l.startswith("!Sample_title")))
ctypes    = split_soft(next((l for l in lines if "cell type:" in l), None)) or [""] * len(gsms)
sources   = split_soft(next((l for l in lines if "source_name" in l), None)) or [""] * len(gsms)

# Collect all treatment lines (there can be multiple per sample in GEO)
trt_map = {}
for tl in (l for l in lines if "treatment:" in l):
    vals = split_soft(tl)
    for gsm, v in zip(gsms, vals[1:]):
        trt_map[gsm] = re.sub(r"^treatment:\s*", "", v).strip('"')

rows = []
for i, gsm in enumerate(gsms):
    trt_raw = trt_map.get(gsm, "unknown")
    # Normalize to naive / treated
    if trt_raw == "Untreated":
        status = "naive"
    else:
        status = "treated"
    # Extract patient ID from title
    m = re.search(r"(PDAC_[UT]_\d+)", titles[i] if i < len(titles) else "")
    pid = m.group(1) if m else gsm
    rows.append({
        "GSM": gsm,
        "patient_id": pid,
        "title": (titles[i][:60] if i < len(titles) else ""),
        "cell_type": ctypes[i] if i < len(ctypes) else "",
        "source_name": sources[i] if i < len(sources) else "",
        "raw_treatment": trt_raw,
        "status": status,
    })

df = pd.DataFrame(rows)
print(f"Total samples: {len(df)}")
print(f"\nStatus distribution:")
print(df["status"].value_counts().to_string())
print(f"\nRaw treatment breakdown:")
print(df["raw_treatment"].value_counts().to_string())
print(f"\nSample types:")
print(df[df["cell_type"].str.contains("Normal|Organoid", case=False, na=False)].head())
# Filter to tumor samples only for analysis
tumor = df[df["cell_type"].str.contains("Tumor", na=False)].copy()
print(f"\nTumor samples: {len(tumor)}")
print(tumor["status"].value_counts().to_string())
tumor.to_csv(os.path.join(DATA, "sample_metadata.tsv"), sep="\t", index=False)
print(f"\nSaved → sample_metadata.tsv ({len(df)} rows, {len(tumor)} tumor-only)")
