"""Download the 12 GSE300595 10x Multiome filtered_feature_bc_matrix.h5 files.

Resumable: skips files already present with the expected size.
"""
import os
import sys
import time
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data")
os.makedirs(DATA, exist_ok=True)

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300595/suppl/"
SAMPLES = ["pc99", "pc102", "pc123", "pc131", "pc149", "pc158",
           "pc159", "pc167", "pc169", "pc172", "pc173", "pc187"]

# expected sizes in bytes (approx, from FTP listing); used only for sanity check
EXPECTED_MB = {"pc99": 32, "pc102": 12, "pc123": 30, "pc131": 39,
               "pc149": 11, "pc158": 37, "pc159": 12, "pc167": 31,
               "pc169": 26, "pc172": 7.7, "pc173": 32, "pc187": 15}


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}TB"


def download(name):
    fname = f"GSE300595_{name}.filtered_feature_bc_matrix.h5"
    dst = os.path.join(DATA, fname)
    if os.path.exists(dst) and os.path.getsize(dst) > 1_000_000:
        print(f"[skip] {fname} already present ({human(os.path.getsize(dst))})", flush=True)
        return
    url = FTP + fname
    print(f"[get ] {fname}", flush=True)
    t0 = time.time()
    tmp = dst + ".part"
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
        total = int(r.headers.get("Content-Length", 0))
        done = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            done += len(chunk)
            pct = 100.0 * done / total if total else 0
            print(f"       {human(done)} / {human(total)} ({pct:5.1f}%)", end="\r", flush=True)
    os.replace(tmp, dst)
    print(f"\n       done in {time.time()-t0:.1f}s -> {human(os.path.getsize(dst))}", flush=True)


if __name__ == "__main__":
    t0 = time.time()
    for s in SAMPLES:
        try:
            download(s)
        except Exception as e:
            print(f"[FAIL] {s}: {e}", flush=True)
    print(f"\nAll done in {time.time()-t0:.1f}s", flush=True)
    total = sum(os.path.getsize(os.path.join(DATA, f))
                for f in os.listdir(DATA) if f.endswith(".h5"))
    print(f"Total h5 bytes: {human(total)}", flush=True)
