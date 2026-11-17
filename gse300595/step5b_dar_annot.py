"""Step 5b - annotate DARs to protein-coding genes and build direction-specific
accessibility programmes.

Peak-level tests do not survive FDR at n=5 vs 7, so the DAR table is used to
build rank-based PROGRAMMES (top peaks by effect size with unadjusted p<0.05).
These are explicitly exploratory; the bridge in step 6 is the independent test.
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import BASE, PROC

GTF = os.path.join(BASE, "ref", "gencode.v44.basic.annotation.gtf.gz")


def load_pc_genes():
    g = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != "gene":
                continue
            name = gtype = None
            for kv in p[8].split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    name = kv.split(" ", 1)[1].strip('"')
                elif kv.startswith("gene_type "):
                    gtype = kv.split(" ", 1)[1].strip('"')
            if name and gtype == "protein_coding" and name not in g:
                g[name] = (p[0], int(p[3]), int(p[4]))
    return g


def main():
    dar = pd.read_csv(os.path.join(PROC, "myeloid_DAR_pseudobulk.tsv"), sep="\t")
    print(f"DAR rows: {len(dar)}", flush=True)

    genes = load_pc_genes()
    print(f"Protein-coding genes with coordinates: {len(genes)}", flush=True)
    gdf = pd.DataFrame([(g, v[0], v[1], v[2]) for g, v in genes.items()],
                       columns=["gene", "chrom", "start", "end"])
    gchrom = {c: sub.reset_index(drop=True) for c, sub in gdf.groupby("chrom")}

    def nearest(chrom, s, e):
        if not isinstance(chrom, str) or chrom not in gchrom:
            return None, np.nan
        sub = gchrom[chrom]
        mid = (s + e) // 2
        tss = sub["start"].values            # gene start ~ TSS on + strand
        d = np.abs(tss - mid)
        i = int(np.argmin(d))
        return sub["gene"].values[i], float(d[i])

    ng = [nearest(r.chrom, r.start, r.end) for r in dar.itertuples()]
    dar["nearest_pc_gene"] = [a[0] for a in ng]
    dar["dist_to_tss_bp"] = [a[1] for a in ng]
    dar.to_csv(os.path.join(PROC, "myeloid_DAR_annotated.tsv"), sep="\t", index=False)

    cols = ["peak", "chrom", "start", "end", "cpm_R", "cpm_NR", "log2FC",
            "p", "nearest_pc_gene", "dist_to_tss_bp"]
    print("\n=== Top 25 peaks MORE ACCESSIBLE IN NON-RESPONDERS ===", flush=True)
    print(dar.sort_values("log2FC").head(25)[cols].to_string(index=False), flush=True)
    print("\n=== Top 25 peaks MORE ACCESSIBLE IN RESPONDERS ===", flush=True)
    print(dar.sort_values("log2FC", ascending=False).head(25)[cols]
          .to_string(index=False), flush=True)

    # ---- programmes -------------------------------------------------------
    # effect-size ranking with a permissive unadjusted filter, NOT FDR
    sel_nr = dar[(dar["log2FC"] < -0.5) & (dar["p"] < 0.05)]
    sel_r = dar[(dar["log2FC"] > 0.5) & (dar["p"] < 0.05)]
    print(f"\nNR-enriched peaks (log2FC<-0.5, p<0.05): {len(sel_nr)}", flush=True)
    print(f"R-enriched  peaks (log2FC> 0.5, p<0.05): {len(sel_r)}", flush=True)

    for tag, sub in [("NONRESPONDER", sel_nr), ("RESPONDER", sel_r)]:
        sub = sub.sort_values("p")
        genes_dir = [g for g in sub["nearest_pc_gene"].dropna().unique() if g]
        pd.Series(genes_dir).to_csv(
            os.path.join(PROC, f"DAR_pcgenes_{tag}.txt"), index=False, header=False)
        print(f"\n{tag}: {len(genes_dir)} unique protein-coding genes", flush=True)
        print(f"  {genes_dir[:40]}", flush=True)

    # also a wider, purely effect-size-based list (top 200 each way)
    for tag, sub in [("NONRESPONDER_top200", dar.sort_values("log2FC").head(200)),
                     ("RESPONDER_top200", dar.sort_values("log2FC", ascending=False).head(200))]:
        genes_dir = [g for g in sub["nearest_pc_gene"].dropna().unique() if g]
        pd.Series(genes_dir).to_csv(
            os.path.join(PROC, f"DAR_pcgenes_{tag}.txt"), index=False, header=False)
        print(f"{tag}: {len(genes_dir)} genes", flush=True)

    print("\n[saved] myeloid_DAR_annotated.tsv + DAR_pcgenes_*.txt", flush=True)


if __name__ == "__main__":
    main()
