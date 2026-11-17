import pandas as pd
import os

PROC = r"E:\pancreatic_cancer\submit\gse300595\processed"

df = pd.read_csv(os.path.join(PROC, "leiden_top20.tsv"), sep="\t")
counts = pd.read_csv(os.path.join(PROC, "cluster_composition.tsv"), sep="\t", index_col=0)

for g, sub in df.groupby("group", sort=True):
    genes = ", ".join(sub["names"].head(14).tolist())
    n = int(counts.loc[int(g), "total"]) if int(g) in counts.index else -1
    resp = ""
    if int(g) in counts.index:
        resp = f"R={int(counts.loc[int(g),'responder'])} NR={int(counts.loc[int(g),'nonresponder'])}"
    print(f"{str(g):>3s} n={n:5d} [{resp:>16s}]  {genes}")
