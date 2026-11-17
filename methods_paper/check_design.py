"""Check whether the GeoMx design is a crossed or nested structure.

Critical question raised by design_effect.py: ICC for the T-cell score is ~0, so
pooling segments is NOT strongly pseudo-replicated.  Why, then, does the p-value
collapse when aggregating to physical ROI?

Hypothesis: each physical ROI contributes at most ONE segment per compartment
(Immune / Epithelial / CAF).  If so:
  - POOLING compartments mixes three different tissues and is confounded,
  - WITHIN a compartment, each observation is a distinct physical ROI, i.e. not
    pseudo-replicated at ROI level at all.
The only remaining clustering is at PATIENT level.
"""
import os
import gzip
import numpy as np
import pandas as pd

DSP_DIR = r"E:\pancreatic_cancer\submit\gse202051_hum\dsp_data"

seg = pd.read_csv(gzip.open(os.path.join(DSP_DIR, "segment_props.txt.gz"), "rt"),
                  sep="\t")
seg["phys_roi"] = seg["Slide_name"] + "_" + seg["ROI_number"].astype(str)

dup = seg.groupby(["phys_roi", "Segment"]).size()
print("duplicated (physical ROI, compartment) pairs:",
      int((dup > 1).sum()), "of", len(dup))
print("max segments per (ROI, compartment):", int(dup.max()))

print("\nsegments per ROI:", seg.groupby("phys_roi").size().value_counts().to_dict())
print("patients:", seg["Patient"].nunique())
print("ROIs per patient:", seg.groupby("Patient")["phys_roi"].nunique().describe().round(2).to_dict())

# is every ROI inside exactly one patient?
per_roi_pat = seg.groupby("phys_roi")["Patient"].nunique()
print("ROIs spanning >1 patient:", int((per_roi_pat > 1).sum()))

print("\ncompartment counts per patient (all 22):")
print(seg.groupby(["Patient", "Segment"]).size().unstack(fill_value=0).to_string())
