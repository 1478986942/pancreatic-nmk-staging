"""Metadata-only audit of candidate ICB cohorts.

Two modes:

  python audit_icb_cohorts.py GSE248014 GSE197613 ...
      read each GEO series-matrix header (no expression) and summarise every
      !Sample_characteristics key with value counts.

  python audit_icb_cohorts.py --inspect <file> [<file> ...]
      read a downloaded supplementary metadata file and report its columns,
      the first rows, and any response / cell-type columns.

  python audit_icb_cohorts.py --suppl GSE120575 ...
      list the supplementary files published for each series.

  python audit_icb_cohorts.py --audit-imvigor [data_dir]
      run the single-cohort IMvigor210 baseline pilot (Mann-Whitney U on
      SERPINE1 / SPP1 / MARCO, R vs NR). Downloads from ICB-Portal if absent.

This answers: does the cohort have a response annotation, at what n, and does
it contain a cell-type label that could identify tumour cells.
"""
import gzip
import os
import sys
import urllib.request
from collections import Counter, defaultdict

CACHE = r"C:\Users\Administrator\AppData\Local\Temp\opencode"

ACC = ["GSE248014", "GSE197613", "GSE313101", "GSE279781"]


def list_suppl(acc):
    pre = acc[: max(3, len(acc) - 3)]
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}nnn/{acc}/suppl/")
    try:
        html = urllib.request.urlopen(url, timeout=60).read().decode(
            "utf-8", "replace")
    except Exception as e:
        print(f"  ERROR {type(e).__name__}: {e}")
        return
    names = []
    for chunk in html.split('href="')[1:]:
        name = chunk.split('"')[0]
        if name.startswith(("?", "/", "..", "http")):
            continue
        names.append(name)
    print(f"  {len(names)} supplementary files")
    for n in names:
        print(f"    {n}")


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def inspect(path):
    print("=" * 78)
    print(path)
    try:
        with _open(path) as f:
            head = [next(f) for _ in range(4)]
    except StopIteration:
        pass
    except Exception as e:
        print(f"  ERROR {type(e).__name__}: {e}")
        return
    sep = "\t" if head[0].count("\t") >= head[0].count(",") else ","
    cols = [c.strip().strip('"') for c in head[0].rstrip("\n").split(sep)]
    print(f"  delimiter={'TAB' if sep == chr(9) else 'COMMA'}  n_columns={len(cols)}")
    print(f"  columns: {cols[:40]}")
    for i, line in enumerate(head[1:4], 1):
        print(f"  row{i}: {line.rstrip()[:200]}")
    low = [c.lower() for c in cols]
    hits = [c for c, l in zip(cols, low)
            if any(k in l for k in ("respon", "benefit", "recist", "outcome",
                                    "celltype", "cell_type", "cell type",
                                    "annotation", "cluster", "malig", "tumor",
                                    "patient", "sample"))]
    print(f"  relevant columns: {hits}")
    print()


def inspect_rds(path):
    """Report the columns of a .rds object without loading it into R.

    Works for data.frame / matrix payloads. A Seurat object is a nested S4
    list and pyreadr will refuse it; that failure is itself the finding and is
    reported rather than hidden.
    """
    print("=" * 78)
    print(path)
    try:
        import pyreadr
    except ImportError:
        print("  ERROR pyreadr not installed")
        return
    try:
        res = pyreadr.read_r(path)
    except Exception as e:
        print(f"  ERROR {type(e).__name__}: {str(e)[:300]}")
        print("  -> not a flat data.frame; colData must be read with R")
        return
    for key, obj in res.items():
        print(f"  object={key!r} type={type(obj).__name__}")
        print(f"  shape={getattr(obj, 'shape', None)}")
        cols = [str(c) for c in getattr(obj, "columns", [])]
        print(f"  n_columns={len(cols)}")
        print(f"  columns={cols[:80]}")
        low = [c.lower() for c in cols]
        hits = [c for c, l in zip(cols, low)
                if any(k in l for k in ("respon", "benefit", "recist",
                                        "outcome", "patient", "sample",
                                        "treat", "time", "status"))]
        print(f"  relevant columns: {hits}")
        if hits:
            print(obj[hits].head(6).to_string()[:1500])
    print()


def inspect_rds_strings(path, keywords=None):
    """Scan an .rds binary for printable strings that look like metadata keys.

    An RDS file serialises symbol names as length-prefixed strings, so column
    names such as `response` or `RECIST` survive in the raw bytes even when the
    object is an S4/Seurat tree that pyreadr cannot decode. This reports those
    names; it does not read or analyse expression values.
    """
    import re
    if keywords is None:
        keywords = ("respon", "recist", "benefit", "relapse", "progress",
                    "outcome", "efficacy", "nonresp", "responder", "pfs",
                    "survival", "status", "patient", "sample", "treat")
    print("=" * 78)
    print(path)
    try:
        with open(path, "rb") as f:
            raw = f.read()
    except Exception as e:
        print(f"  ERROR {type(e).__name__}: {e}")
        return
    # saveRDS defaults to gzip; bzip2 and xz are also possible.
    magic = raw[:6]
    try:
        if magic[:2] == b"\x1f\x8b":
            import gzip as _gz
            raw = _gz.decompress(raw)
            print("  compression = gzip")
        elif magic[:3] == b"BZh":
            import bz2
            raw = bz2.decompress(raw)
            print("  compression = bzip2")
        elif magic[:6] == b"\xfd7zXZ\x00":
            import lzma
            raw = lzma.decompress(raw)
            print("  compression = xz")
        else:
            print("  compression = none")
    except Exception as e:
        print(f"  decompress ERROR {type(e).__name__}: {str(e)[:200]}")
        return
    runs = re.findall(rb"[\x20-\x7e]{4,80}", raw)
    seen = set()
    hits = []
    for run in runs:
        s = run.decode("ascii", "replace")
        low = s.lower()
        if any(k in low for k in keywords) and s not in seen:
            seen.add(s)
            hits.append(s)
    print(f"  file bytes = {len(raw)}")
    print(f"  printable runs = {len(runs)}   keyword hits = {len(hits)}")
    for s in hits[:80]:
        print(f"    {s}")
    print()


def _bh(pvals):
    """Benjamini-Hochberg q-values, no statsmodels dependency."""
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m, dtype=float)
    out[order] = np.minimum(q, 1.0)
    return out


IMVIGOR_BASE = "https://ngdc.cncb.ac.cn/icb/file/datasets/"
IMVIGOR_FILES = {
    "clinical": "Mariathasan_2018_clinical.tsv",
    "expression": "Mariathasan_2018_expression.tsv",
    "rdata": "Mariathasan_2018.Rdata",
}


def _fetch_imvigor(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    paths = {}
    for key, name in IMVIGOR_FILES.items():
        path = os.path.join(data_dir, name)
        if not os.path.exists(path) or os.path.getsize(path) < 1000:
            print(f"  downloading {name} ...")
            urllib.request.urlretrieve(IMVIGOR_BASE + name, path)
        paths[key] = path
    return paths


def audit_imvigor210(data_dir=CACHE):
    """IMvigor210 baseline pilot: three pre-specified genes, R vs NR.

    Design, fixed before looking at the result:
      - one cohort (Mariathasan_2018 / IMvigor210, urothelial, atezolizumab)
      - baseline only (`Pre_On == "Pre"`)
      - response = RECIST CR/PR vs SD/PD; NE excluded
      - primary family = SERPINE1, SPP1, MARCO; two-sided Mann-Whitney U;
        Benjamini-Hochberg across the three
      - secondary = mean z-score of the three genes, nominal only, exploratory
      - descriptive = pairwise Spearman rho among the three genes

    The ICB-Portal expression TSV has **no gene-identifier column**. Gene names
    live in the `.Rdata` (`TPM` dim_0). This function therefore reads the Rdata
    for names and cross-checks its values against the TSV.
    """
    import warnings

    import numpy as np
    import pandas as pd
    from scipy.stats import mannwhitneyu, spearmanr

    paths = _fetch_imvigor(data_dir)

    print("=" * 78)
    print("IMvigor210 / Mariathasan_2018 — baseline bulk pilot")
    print("=" * 78)

    clin = pd.read_csv(paths["clinical"], sep="\t")
    assert "Sample" in clin.columns, "Fatal: clinical table has no Sample column"
    assert clin["Sample"].is_unique, "Fatal: duplicated sample IDs in clinical table"
    assert "Pre_On" in clin.columns, "Fatal: no Pre_On column — cannot isolate baseline"
    assert "Treatment" in clin.columns, "Fatal: no Treatment column"

    pre_vals = set(clin["Pre_On"].dropna().astype(str).str.strip())
    assert pre_vals == {"Pre"}, f"Fatal: non-baseline samples present: {pre_vals}"
    treat_vals = set(clin["Treatment"].dropna().astype(str).str.strip())
    assert treat_vals == {"anti-PDL1"}, f"Fatal: mixed therapy: {treat_vals}"

    # Gene names come from the Rdata; the TSV alone cannot be mapped.
    try:
        import rdata
    except ImportError as exc:
        raise SystemExit(
            "rdata is required to recover gene names: pip install rdata") from exc
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        root = rdata.conversion.convert(
            rdata.parser.parse_file(paths["rdata"]))["Mariathasan_2018"]
    tpm = root["TPM"]
    genes = [str(g) for g in tpm.coords["dim_0"].values]
    samples = [str(s) for s in tpm.coords["dim_1"].values]
    mat = np.asarray(tpm.values, dtype=float)

    assert len(set(samples)) == len(samples), "Fatal: duplicated sample IDs in TPM"
    assert set(samples) == set(clin["Sample"]), \
        "Fatal: TPM and clinical sample sets differ"

    # Some TPM builds keep one row per transcript, so a symbol can repeat.
    # Count first; aggregate only after the TSV cross-check, so aggregation
    # cannot hide a mismatch.
    from collections import Counter
    gene_counts = Counter(genes)
    n_dup_symbols = sum(1 for c in gene_counts.values() if c > 1)
    print(f"  gene rows={len(genes)}  unique symbols={len(gene_counts)}"
          f"  duplicated symbols={n_dup_symbols}")

    # The TSV is the same matrix with row names dropped. Prove it before using
    # the Rdata names on it.
    tsv = pd.read_csv(paths["expression"], sep="\t")
    assert list(tsv.columns) == samples, \
        "Fatal: expression TSV column order differs from Rdata TPM"
    assert tsv.shape == mat.shape, \
        f"Fatal: TSV {tsv.shape} vs TPM {mat.shape}"
    assert np.allclose(tsv.values.astype(float), mat, rtol=1e-6, atol=1e-9,
                       equal_nan=True), \
        "Fatal: expression TSV values differ from Rdata TPM"

    if n_dup_symbols:
        frame = pd.DataFrame(mat, index=genes)
        agg = frame.groupby(level=0, sort=False).mean()
        genes = [str(g) for g in agg.index]
        mat = agg.to_numpy()
        print(f"  duplicate symbols averaged by mean -> {mat.shape[0]} genes")

    clin = clin.set_index("Sample").loc[samples]

    resp_col = "Resp_NoResp" if "Resp_NoResp" in clin.columns else "binaryResponse"
    assert resp_col in clin.columns, "Fatal: no response column"
    resp = clin[resp_col].astype("string").str.strip()
    codes = set(resp.dropna().unique())
    assert codes == {"Response", "No_Response"}, \
        f"Fatal: unexpected response coding: {codes}"

    is_r = resp.eq("Response").fillna(False).to_numpy(dtype=bool)
    is_nr = resp.eq("No_Response").fillna(False).to_numpy(dtype=bool)
    n_ne = int(resp.isna().sum())
    assert not (is_r & is_nr).any(), "Fatal: a sample is both R and NR"
    assert is_r.sum() + is_nr.sum() + n_ne == len(resp), \
        "Fatal: response labels do not partition the cohort"

    print(f"  samples total        : {len(resp)}")
    print(f"  baseline (Pre_On=Pre): {len(resp)}  [asserted]")
    print(f"  treatment            : {sorted(treat_vals)}")
    print(f"  Responder (CR/PR)    : {int(is_r.sum())}")
    print(f"  Non-responder (SD/PD): {int(is_nr.sum())}")
    print(f"  excluded (NE/NA)     : {n_ne}")
    print()

    targets = ["SERPINE1", "SPP1", "MARCO"]
    pos = {g: i for i, g in enumerate(genes)}
    missing = [g for g in targets if g not in pos]
    assert not missing, f"Fatal: genes absent from matrix: {missing}"
    dup_targets = [g for g in targets if gene_counts.get(g, 0) > 1]
    if dup_targets:
        print(f"  NOTE: target(s) had multiple transcript rows, averaged: {dup_targets}")

    rows = []
    for g in targets:
        vals = mat[pos[g]]
        r = vals[is_r]
        nr = vals[is_nr]
        assert np.isfinite(r).all() and np.isfinite(nr).all(), \
            f"Fatal: non-finite values for {g}"
        u, p = mannwhitneyu(r, nr, alternative="two-sided")
        # Cliff's delta: +1 = every responder above every non-responder,
        # -1 = every responder below. This matches the median direction.
        delta = 2.0 * u / (r.size * nr.size) - 1.0
        rows.append({
            "gene": g,
            "median_R": float(np.median(r)),
            "median_NR": float(np.median(nr)),
            "mean_R": float(r.mean()),
            "mean_NR": float(nr.mean()),
            "U": float(u),
            "p": float(p),
            "delta": float(delta),
        })

    q = _bh([row["p"] for row in rows])
    for row, qv in zip(rows, q):
        row["q"] = float(qv)

    print("-" * 78)
    print("  primary: three pre-specified genes, Mann-Whitney U (two-sided)")
    print("-" * 78)
    header = (f"  {'gene':<9}{'med_R':>9}{'med_NR':>10}"
              f"{'mean_R':>10}{'mean_NR':>10}{'U':>10}{'p':>11}{'q_BH':>10}{'delta':>8}")
    print(header)
    for row in rows:
        print(f"  {row['gene']:<9}{row['median_R']:>9.2f}{row['median_NR']:>10.2f}"
              f"{row['mean_R']:>10.2f}{row['mean_NR']:>10.2f}{row['U']:>10.0f}"
              f"{row['p']:>11.4g}{row['q']:>10.4g}{row['delta']:>8.3f}")
    print()
    print("  delta = Cliff's delta; +1 = higher in responders,")
    print("          -1 = higher in non-responders.")
    print()

    # Secondary, exploratory: the composite the hypothesis actually names.
    labelled = is_r | is_nr
    z = np.zeros((len(targets), int(labelled.sum())), dtype=float)
    for i, g in enumerate(targets):
        v = mat[pos[g]][labelled]
        sd = v.std(ddof=0)
        z[i] = (v - v.mean()) / sd if sd > 0 else 0.0
    comp = z.mean(axis=0)
    cr = comp[is_r[labelled]]
    cnr = comp[is_nr[labelled]]
    cu, cp = mannwhitneyu(cr, cnr, alternative="two-sided")
    cdelta = 2.0 * cu / (cr.size * cnr.size) - 1.0
    print("-" * 78)
    print("  secondary (exploratory): mean z-score of the three genes")
    print("-" * 78)
    print(f"  median_R={np.median(cr):.3f}  median_NR={np.median(cnr):.3f}"
          f"  U={cu:.0f}  p={cp:.4g}  delta={cdelta:.3f}")
    print()

    # Descriptive: do the three genes even co-vary?
    print("-" * 78)
    print("  descriptive: pairwise Spearman rho among the three genes (labelled n)")
    print("-" * 78)
    for i in range(len(targets)):
        for j in range(i + 1, len(targets)):
            rho, _ = spearmanr(mat[pos[targets[i]]][labelled],
                               mat[pos[targets[j]]][labelled])
            print(f"  {targets[i]:<9} vs {targets[j]:<9} rho={rho:+.3f}")
    print()

    return rows, {"n_R": int(is_r.sum()), "n_NR": int(is_nr.sum()),
                  "n_NE": n_ne, "composite_p": float(cp),
                  "composite_delta": float(cdelta)}


def fetch_header(acc):
    pre = acc[: max(3, len(acc) - 3)]
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{pre}nnn/{acc}"
           f"/matrix/{acc}_series_matrix.txt.gz")
    path = os.path.join(CACHE, f"{acc}_sm.txt.gz")
    if not os.path.exists(path) or os.path.getsize(path) < 200:
        urllib.request.urlretrieve(url, path)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        return f.read()


def parse(text):
    titles = []
    chars = defaultdict(list)
    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue
        parts = line.split("\t")
        key = parts[0]
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key == "!Sample_title":
            titles = vals
        elif key == "!Sample_characteristics_ch1":
            for v in vals:
                if ":" in v:
                    k, _, val = v.partition(":")
                    chars[k.strip()].append(val.strip())
                elif v:
                    chars["_raw"].append(v)
    return titles, chars


def main(argv):
    if argv and argv[0] == "--inspect":
        for p in argv[1:]:
            inspect(p)
        return
    if argv and argv[0] == "--suppl":
        for acc in argv[1:]:
            print("=" * 78)
            print(acc)
            list_suppl(acc)
            print()
        return
    if argv and argv[0] == "--rds":
        for p in argv[1:]:
            inspect_rds(p)
        return
    if argv and argv[0] == "--rds-strings":
        for p in argv[1:]:
            inspect_rds_strings(p)
        return
    if argv and argv[0] == "--audit-imvigor":
        data_dir = argv[1] if len(argv) > 1 else CACHE
        audit_imvigor210(data_dir)
        return

    for acc in (argv or ACC):
        print("=" * 78)
        print(acc)
        try:
            titles, chars = parse(fetch_header(acc))
        except Exception as e:
            print("  ERROR", type(e).__name__, e)
            continue
        print(f"  n_samples = {len(titles)}")
        if titles:
            print(f"  first titles: {titles[:3]}")
        for k, v in chars.items():
            if not v:
                continue
            c = Counter(v)
            if len(c) <= 12:
                print(f"  [{k}] " + ", ".join(f"{a}={b}" for a, b in c.most_common()))
            else:
                print(f"  [{k}] {len(c)} distinct; top: "
                      + ", ".join(f"{a}={b}" for a, b in c.most_common(6)))
        print()


if __name__ == "__main__":
    main(sys.argv[1:])
