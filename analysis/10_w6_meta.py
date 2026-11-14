import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
# -*- coding: utf-8 -*-
import pyreadstat, os, json
base=os.path.join(RAW, "CHARLS")
mods=["Health_Status_and_Functioning.dta","Health_Care_Insurance.dta","Sample_Infor.dta"]
out=[]
for y in ["2021","2022","2023"]:
    for m in mods:
        p=os.path.join(base,y,m)
        if not os.path.exists(p): 
            out.append(f"\n## {y}/{m} MISSING"); continue
        df,meta=pyreadstat.read_dta(p,metadataonly=True)
        out.append(f"\n{'='*80}\n## {y}/{m}  rows={meta.number_rows} cols={meta.number_columns}\n{'='*80}")
        for c,l in zip(meta.column_names, meta.column_labels):
            out.append(f"{c}\t{l}")
open(os.path.join(OUT,"10_w6_meta.txt"),"w",encoding="utf-8").write("\n".join(out))
print("ok")
