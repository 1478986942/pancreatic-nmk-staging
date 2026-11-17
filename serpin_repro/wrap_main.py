import re

p = r"E:\pancreatic_cancer\submit\serpin_repro\step3_annotate_cnv.py"
s = open(p, encoding="utf-8").read()
marker = 'adata = sc.read_h5ad(os.path.join(OUT_DIR, "GSE289138_processed.h5ad"))'
if "def main():" in s:
    print("already wrapped")
else:
    i = s.index(marker)
    head = s[:i]
    body = s[i:]
    body = "\n".join("    " + ln if ln.strip() else ln for ln in body.split("\n"))
    new = head + "def main():\n" + body + "\n\nif __name__ == \"__main__\":\n    main()\n"
    open(p, "w", encoding="utf-8").write(new)
    print("wrapped main()")
