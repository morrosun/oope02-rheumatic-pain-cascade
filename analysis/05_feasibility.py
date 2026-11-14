# -*- coding: utf-8 -*-
"""Feasibility probe: pain / diagnosis / treatment cascade variables, by cohort."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, json
import pyreadstat
import pandas as pd

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
FILES = {
 "CHARLS": (os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta"),
            ["ID","wave","age","ragender","raeducl","arthre","rxarthr","rxarthr_c",
             "da042s1","da042s2","da042s3","da042s4","da042s5","da042s7","da042s8","da042s9",
             "da042s10","da042s11","da042s12","da042s13","da042s15",
             "pain","adlab_c","iadl","iwstat","radyear"]),
 "KLoSA":  (os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Working_data/klosa.dta"),
            ["pid","wave","agey","ragender","raeducl","arthre","rxarthr",
             "pain_1","pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13",
             "painfr","painhlv","paina","adltot6","iwstat","radyear"]),
 "LASI":   (os.path.join(DL, "5. LASI 印度/LASI_印度/Working_data/lasi.dta"),
            ["prim_key","r1agey","ragender","raeducl","r1arthre","r1rxarthr","r1jointp","r1painfr","r1paina",
             "r1adltot6","r1iadltot_l"]),
 "MHAS":   (os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Working_data/mhas.dta"),
            ["rahhidnp","wave","agey","ragender","raeducl","arthre","rxarthr","painfr","painlv","paina",
             "adltot6","iwstat","radyear"]),
 "HRS":    (os.path.join(DL, "3. HRS  美国/HRS_美国/Working_data/hrs.dta"),
            ["hhidpn","wave","ragey_b","ragender","raeducl","arthre","rxarthr","rxpain","painfr","painlv","paina",
             "adl6a","iadl5a","iwstat","radyear"]),
 "SHARE":  (os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Working_data/share.dta"),
            ["mergeid","wave","country","agey","ragender","raeducl","arthre","rxinflm","rxpain","rxosteo",
             "pain_s","painfr","painlv","adltot6"]),
 "ELSA":   (os.path.join(DL, "2. ELSA 英国/Working_data/elsa.dta"),
            ["idauniqc","wave","agey","ragender","raeducl","arthre","rxosteo","painfr","painlv","adltot6"]),
}

lines=[]
def w(s=""):
    lines.append(s); print(s, flush=True)

for db,(path,cols) in FILES.items():
    w("\n" + "="*84); w(f"## {db}  {os.path.basename(path)}")
    if not os.path.exists(path):
        w("MISSING"); continue
    df, meta = pyreadstat.read_dta(path, usecols=cols)
    cols = [c for c in cols if c in df.columns]
    w(f"rows={len(df)}  cols_loaded={df.shape[1]}")
    if "wave" in df: w("waves: " + str(df["wave"].value_counts().sort_index().to_dict()))
    for c in cols:
        if c in ("ID","pid","prim_key","rahhidnp","hhidpn","mergeid","wave","country"): continue
        s = df[c]
        nm = int(s.notna().sum())
        if nm==0:
            w(f"  {c:14s} ALL MISSING"); continue
        vc = s.value_counts(dropna=True).head(5)
        w(f"  {c:14s} n={nm:7d}  " + ", ".join(f"{k}:{v}" for k,v in vc.items()))
    # cascade cross-tabs
    try:
        if "arthre" in df and "painfr" in df:
            t = pd.crosstab(df["painfr"].fillna(-1), df["arthre"].fillna(-1))
            w("  painfr x arthre:"); w(t.to_string())
        if "arthre" in df and "pain_s" in df:
            t = pd.crosstab(df["pain_s"].fillna(-1), df["arthre"].fillna(-1)); w("  pain_s x arthre:"); w(t.to_string())
        if "r1arthre" in df and "r1jointp" in df:
            t = pd.crosstab(df["r1jointp"].fillna(-1), df["r1arthre"].fillna(-1)); w("  r1jointp x r1arthre:"); w(t.to_string())
        if "arthre" in df and "rxarthr" in df:
            sub = df[df["arthre"]==1]
            w(f"  among arthre==1 (n={len(sub)}): rxarthr dist = " + str(sub["rxarthr"].value_counts(dropna=False).head(6).to_dict()))
    except Exception as e:
        w(f"  crosstab ERR {e}")

with open(os.path.join(OUT,"05_feasibility.txt"),"w",encoding="utf-8") as f:
    f.write("\n".join(lines))
print("WROTE 05_feasibility.txt")
