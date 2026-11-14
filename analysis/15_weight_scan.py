import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
# -*- coding: utf-8 -*-
import pyreadstat, os, re
FILES = {
 "CHARLS_h": os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Raw_data/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta"),
 "ELSA_h":   os.path.join(DL, "2. ELSA 英国/Raw_data/Harmonized ELSA/h_elsa_g3.dta"),
 "HRS_h":    os.path.join(DL, "3. HRS  美国/HRS_美国/Raw_data/Gateway Harmonized HRS/H_HRS_d.dta"),
 "KLoSA_h":  os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Raw_data/H_KLoSA/H_KLoSA_e2.dta"),
 "LASI_h":   os.path.join(DL, "5. LASI 印度/LASI_印度/Raw_data/Harmonized LASI (A.3)/H_LASI_a3.dta"),
 "MHAS_h":   os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Raw_data/Harmonized MHAS File/H_MHAS_c2.dta"),
 "SHARE_h":  os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Raw_data/Harmonized SHARE/H_SHARE_f2.dta"),
}
for k,p in FILES.items():
    if not os.path.exists(p):
        print(k,"MISSING",p); continue
    df,meta=pyreadstat.read_dta(p, metadataonly=True)
    hits=[(c,l) for c,l in zip(meta.column_names, meta.column_labels) if re.search(r'weight|wt', str(c)+" "+str(l), re.I)]
    print(f"\n## {k} rows={meta.number_rows} cols={meta.number_columns}  weight-like={len(hits)}")
    for c,l in hits[:14]:
        print(f"   {c:14s} {l}")
