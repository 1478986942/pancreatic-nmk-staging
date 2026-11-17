"""Search GEO for PDAC datasets with chemotherapy +/- spatial or single-cell.

Queries NCBI E-utilities (esearch/esummary) and reports candidate accessions.
"""

import json
import sys
import time
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = {
    "PDAC_scRNA_chemo": (
        'pancreatic[All Fields] AND ("single cell"[All Fields] OR scRNA[All Fields] '
        'OR "single-cell"[All Fields]) AND (chemotherapy[All Fields] OR gemcitabine[All Fields] '
        'OR "nab-paclitaxel"[All Fields] OR abraxane[All Fields] OR FOLFIRINOX[All Fields] '
        'OR "neoadjuvant"[All Fields])'
    ),
    "PDAC_Visium": (
        'pancreatic[All Fields] AND (visium[All Fields] OR "spatial transcriptomics"[All Fields] '
        'OR "spatially resolved"[All Fields])'
    ),
    "PDAC_neoadjuvant_any": (
        'pancreatic[All Fields] AND (neoadjuvant[All Fields] OR "pre-treatment"[All Fields] '
        'OR "pretreatment"[All Fields] OR "treatment naive"[All Fields] '
        'OR "untreated"[All Fields])'
    ),
    "pancreas_scRNA_general": (
        'pancreatic[All Fields] AND ("single cell"[All Fields] OR "single-cell"[All Fields])'
    ),
    "pancreas_spatial_general": (
        'pancreas[All Fields] AND (visium[All Fields] OR "spatial transcriptomics"[All Fields])'
    ),
}


def get(url):
    with urllib.request.urlopen(url, timeout=90) as r:
        return r.read().decode("utf-8", errors="replace")


def esearch(term, n=60):
    url = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "term": term, "retmax": n, "retmode": "json"}
    )
    data = json.loads(get(url))
    return data["esearchresult"].get("idlist", []), int(data["esearchresult"].get("count", 0))


def esummary(ids):
    if not ids:
        return []
    url = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "id": ",".join(ids), "retmode": "json"}
    )
    data = json.loads(get(url))
    out = []
    for uid in data["result"]["uids"]:
        r = data["result"][uid]
        out.append(
            {
                "gse": r.get("accession"),
                "title": r.get("title", "")[:200],
                "n_samples": r.get("n_samples"),
                "type": r.get("gdstype"),
                "taxon": r.get("taxon"),
                "pubdate": r.get("pubdate"),
                "summary": r.get("summary", "")[:400],
            }
        )
    return out


for label, term in QUERIES.items():
    ids, count = esearch(term)
    print(f"\n{'=' * 78}\n### {label}   hits={count}  retrieved={len(ids)}\n{'=' * 78}")
    for rec in esummary(ids):
        print(f"\n[{rec['gse']}] n={rec['n_samples']} type={rec['type']} {rec['pubdate']}")
        print(f"  title: {rec['title']}")
        print(f"  {rec['summary'][:300]}")
    time.sleep(1)
