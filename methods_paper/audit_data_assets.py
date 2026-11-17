"""Audit what mouse / spatial resources actually exist on disk.

Answers, factually, the question: does the proposed three-layer
'mouse causal -> human ICB -> spatial niche' chain have data behind it?
"""
import os
import glob
import json
import anndata as ad

BASE = r"E:\pancreatic_cancer\submit\serpin_repro"
PROC = os.path.join(BASE, "processed")
DATA = os.path.join(BASE, "data")

print("=== processed h5ad files ===")
for p in sorted(glob.glob(os.path.join(PROC, "*.h5ad"))):
    try:
        a = ad.read_h5ad(p, backed="r")
        print(f"\n{p}")
        print(f"  shape: {a.n_obs} cells x {a.n_vars} genes")
        print(f"  obs cols: {list(a.obs.columns)}")
        for c in a.obs.columns:
            nu = a.obs[c].nunique(dropna=True)
            if nu <= 30:
                print(f"    {c}: {sorted(map(str, a.obs[c].dropna().unique()))[:30]}")
        print(f"  obsm: {list(a.obsm.keys())}")
        print(f"  uns keys: {list(a.uns.keys())[:15]}")
        if "spatial" in a.obsm:
            print("  HAS SPATIAL COORDS")
        a.file.close()
    except Exception as e:
        print(f"\n{p}\n  ERROR {type(e).__name__}: {e}")

print("\n=== data/ tree (non-venv) ===")
for root, dirs, files in os.walk(DATA):
    dirs[:] = [d for d in dirs if d not in ("venv", ".git")]
    lvl = root.replace(DATA, "").count(os.sep)
    if lvl > 3:
        continue
    print("  " * lvl + os.path.basename(root) + "/")
    for f in files[:12]:
        fp = os.path.join(root, f)
        print("  " * (lvl + 1) + f + f"  ({os.path.getsize(fp)/1e6:.1f} MB)")

print("\n=== any spatial / visium / xenium artifacts anywhere under submit ===")
pat = ["*visium*", "*Visium*", "*xenium*", "*Xenium*", "*spatial*", "*Spatial*"]
hits = []
for p in pat:
    hits += glob.glob(os.path.join(r"E:\pancreatic_cancer\submit", "**", p), recursive=True)
seen = set()
for h in hits:
    if "venv" in h or "node_modules" in h or ".git" in h:
        continue
    if h in seen:
        continue
    seen.add(h)
    print("  ", h)
