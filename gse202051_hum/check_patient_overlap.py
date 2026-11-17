"""
Patient mapping between GSE202051 (snRNA, U/T IDs) and GSE199102 (DSP, MGH IDs).
Hwang et al. Nat Genet 2022 used both platforms on the same ~40 patients.
snRNA uses PDAC_U_#/PDAC_T_## convention; DSP uses MGH##### convention.

Strategy:
  1. Check if snRNA h5ad contains any MGH-like fields
  2. Use sample_metadata.tsv title strings to find MGH mentions
  3. If no direct mapping found, treat as two independent cohorts from same study
"""
import os, sys, re, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"E:\pancreatic_cancer\submit\serpin_repro\venv\Lib\site-packages")
import pandas as pd
import numpy as np

BASE = r"E:\pancreatic_cancer\submit\gse202051_hum"
DATA = os.path.join(BASE, "data")
OUT  = os.path.join(BASE, "processed")

# ── 1. Load snRNA patient data ────────────────────────────────────────────────
print("=== snRNA patients (GSE202051) ===", flush=True)
pat_snRNA = pd.read_csv(os.path.join(OUT, "patient_level_serpin_marco_full.tsv"), sep="\t")
snRNA_ids = set(pat_snRNA["pid"].astype(str))
print(f"  Total: {len(pat_snRNA)} patients", flush=True)
print(f"  naive (U):  {sorted(pat_snRNA[pat_snRNA['status']=='naive']['pid'].astype(str))}", flush=True)
print(f"  treated(T): {sorted(pat_snRNA[pat_snRNA['status']=='treated']['pid'].astype(str))}", flush=True)

# Check sample_metadata for any MGH references
meta = pd.read_csv(os.path.join(DATA, "sample_metadata.tsv"), sep="\t")
mg_in_title = meta[meta["title"].str.contains("MGH", case=False, na=False)]
print(f"\n  Titles containing 'MGH': {len(mg_in_title)} / {len(meta)}", flush=True)
if len(mg_in_title) > 0:
    print(mg_in_title[["GSM","patient_id","title"]].to_string(), flush=True)

# ── 2. Load DSP patient data ──────────────────────────────────────────────────
print("\n=== DSP patients (GSE199102) ===", flush=True)
pat_dsp = pd.read_csv(os.path.join(BASE, "processed_dsp", "dsp_patient_level.tsv"), sep="\t")
dsp_ids = set(pat_dsp["Patient"].astype(str))
print(f"  Total: {pat_dsp['Patient'].nunique()} patients", flush=True)
print(f"  naive:  {sorted(pat_dsp[pat_dsp['status']=='naive']['Patient'].astype(str).unique())}", flush=True)
print(f"  treated:{sorted(pat_dsp[pat_dsp['status']=='treated']['Patient'].astype(str).unique())}", flush=True)

# ── 3. Attempt to map via number patterns ────────────────────────────────────
print("\n=== Mapping attempt ===", flush=True)

# Hypothesis: MGH#### may encode year + sequence number
# e.g., MGH003 = early patient, MGH2xxx = later added patients
# Let's check if the numeric part correlates with treatment
def num_part(mid):
    m = re.search(r'(\d+)', mid)
    return int(m.group(1)) if m else 0

pat_dsp["MGH_num"] = pat_dsp["Patient"].apply(num_part)
print(f"\n  DSP patient MGH numbers vs treatment:", flush=True)
for st in ["naive","treated"]:
    sub = pat_dsp[pat_dsp["status"]==st]
    nums = sorted(sub["MGH_num"].unique())
    print(f"    {st}: {[f'MGH{n:03d}' for n in nums]}", flush=True)

# Pattern observation: MGH003-011 are all naive (small numbers)
# MGH2462-2997 are mixed (larger numbers), let's check
naive_mgh = sorted(pat_dsp[pat_dsp["status"]=="naive"]["MGH_num"].unique())
trtd_mgh  = sorted(pat_dsp[pat_dsp["status"]=="treated"]["MGH_num"].unique())
print(f"\n  Naive MGH nums:     {naive_mgh}", flush=True)
print(f"  Treated MGH nums:   {trtd_mgh}", flush=True)

# Check snRNA patient count vs DSP
# Paper says 43 total, 22 in DSP. So ~21 snRNA-only or non-overlapping.
# Key question: do any U/T IDs match MGH IDs?
overlap = snRNA_ids & dsp_ids
print(f"\n  Direct ID overlap: {len(overlap)}", flush=True)

# ── 4. Try mapping by ordering assumption ────────────────────────────────────
# If we assume first-named patients correspond:
# snRNA naive U1-U18 → likely maps to some subset of MGH naive (13 pts)
# snRNA treated T1-T25 → likely maps to some subset of MGH treated (7 pts)
# But this is a guess without explicit mapping table.

# Let's see if there's a known mapping in supplementary files
print("\n=== Searching for mapping table ===", flush=True)
for fname in os.listdir(os.path.join(DATA)):
    fpath = os.path.join(DATA, fname)
    try:
        with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()[:5000]
        # Look for lines with both MGH and PDAC/U/T patterns
        for line in content.split('\n'):
            if ('MGH' in line or 'PDAC_U' in line or 'PDAC_T' in line) and \
               ('U1' in line or 'T1' in line or 'patient' in line.lower()):
                print(f"  {fname}: {line[:200]}", flush=True)
    except:
        pass

print("\n=== CONCLUSION ===", flush=True)
if not overlap:
    print("No shared patient IDs between datasets.", flush=True)
    print("Treat as INDEPENDENT cohorts from the same Hwang et al. study.", flush=True)
    print("Cross-platform consistency检验 is appropriate (not paired analysis).", flush=True)
else:
    print(f"{len(overlap)} overlapping patients found!", flush=True)
