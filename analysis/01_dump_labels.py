# -*- coding: utf-8 -*-
"""Dump variable names + labels for the 8 working-layer files (metadata only)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, json
import pyreadstat

FILES = {
    "CHARLS": os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta"),
    "ELSA":   os.path.join(DL, "2. ELSA 英国/Working_data/elsa.dta"),
    "HRS":    os.path.join(DL, "3. HRS  美国/HRS_美国/Working_data/hrs.dta"),
    "KLoSA":  os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Working_data/klosa.dta"),
    "LASI":   os.path.join(DL, "5. LASI 印度/LASI_印度/Working_data/lasi.dta"),
    "MHAS":   os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Working_data/mhas.dta"),
    "SHARE":  os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Working_data/share.dta"),
    "CLHLS":  os.path.join(DL, "CLHLS/Working_data/clhls08_18.dta"),
}

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
os.makedirs(OUT, exist_ok=True)

allmeta = {}
lines = []
for db, path in FILES.items():
    if not os.path.exists(path):
        lines.append(f"## {db}  MISSING: {path}")
        continue
    df, meta = pyreadstat.read_dta(path, metadataonly=True)
    cols = list(meta.column_names)
    labs = list(meta.column_labels)
    allmeta[db] = {"n_rows": meta.number_rows, "n_cols": meta.number_columns,
                   "vars": [{"v": c, "l": l} for c, l in zip(cols, labs)]}
    lines.append(f"\n{'='*90}\n## {db}  rows={meta.number_rows}  cols={meta.number_columns}\n{'='*90}")
    for c, l in zip(cols, labs):
        lines.append(f"{c}\t{l}")

with open(os.path.join(OUT, "01_var_labels.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

with open(os.path.join(OUT, "01_var_labels.json"), "w", encoding="utf-8") as f:
    json.dump(allmeta, f, ensure_ascii=False, indent=1)

print("done", {k: (v["n_rows"], v["n_cols"]) for k, v in allmeta.items()})
