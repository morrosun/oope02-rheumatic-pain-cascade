# -*- coding: utf-8 -*-
"""16: 从 harmonized 层提取人-波权重，统一成 (cohort, ID, wave, weight)"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pyreadstat, pandas as pd, numpy as np

CFG = [
 ("CHARLS", os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Raw_data/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta"),
  "ID", {w: f"r{w}wtresp" for w in range(1,6)}),
 ("ELSA",   os.path.join(DL, "2. ELSA 英国/Raw_data/Harmonized ELSA/h_elsa_g3.dta"),
  "idauniqc", {w: f"r{w}cwtresp" for w in range(1,10)}),
 ("KLoSA",  os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Raw_data/H_KLoSA/H_KLoSA_e2.dta"),
  "pid", {w: f"r{w}wtresp" for w in range(1,9)}),
 ("LASI",   os.path.join(DL, "5. LASI 印度/LASI_印度/Raw_data/Harmonized LASI (A.3)/H_LASI_a3.dta"),
  "prim_key", {1: "r1wtresp"}),
 ("MHAS",   os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Raw_data/Harmonized MHAS File/H_MHAS_c2.dta"),
  "rahhidnp", {w: f"r{w}wtresp" for w in range(1,6)}),
 ("SHARE",  os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Raw_data/Harmonized SHARE/H_SHARE_f2.dta"),
  "mergeid", {w: f"r{w}wtresp" for w in [1,2,4,5,6,7,8]}),
]

frames=[]
for cohort, path, idv, wmap in CFG:
    cols=[idv]+list(wmap.values())
    df,_=pyreadstat.read_dta(path, usecols=cols, metadataonly=False)
    present={w:v for w,v in wmap.items() if v in df.columns}
    missing=[v for v in wmap.values() if v not in df.columns]
    if missing: print(f"{cohort}: 缺列 {missing}，已跳过")
    long = df[[idv]+list(present.values())].melt(id_vars=[idv], value_vars=list(present.values()),
                                                 var_name="wvar", value_name="weight")
    inv={v:k for k,v in present.items()}
    long["wave"]=long["wvar"].map(inv)
    long=long.rename(columns={idv:"ID"})[["ID","wave","weight"]]
    long["cohort"]=cohort
    n=long["weight"].notna().sum()
    print(f"{cohort:7s} waves={sorted(present)} rows={len(long):7d}  非缺失权重={n:7d}  median={long['weight'].median():.1f}", flush=True)
    frames.append(long)

# HRS 用工作层的 wtresp
hrs,_=pyreadstat.read_dta(os.path.join(DL, "3. HRS  美国/HRS_美国/Working_data/hrs.dta"),
                          usecols=["hhidpn","wave","wtresp"])
hrs=hrs.rename(columns={"hhidpn":"ID","wtresp":"weight"}); hrs["cohort"]="HRS"
print(f"HRS     rows={len(hrs):7d}  非缺失权重={hrs['weight'].notna().sum():7d}  median={hrs['weight'].median():.1f}")
frames.append(hrs)

w=pd.concat(frames, ignore_index=True)
w["ID"]=w["ID"].astype(str)
w.to_csv(os.path.join(OUT,"16_weights.csv"), index=False, encoding="utf-8-sig")
print("WROTE 16_weights.csv", w.shape)
