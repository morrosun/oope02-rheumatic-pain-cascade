import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
# -*- coding: utf-8 -*-
import pyreadstat, pandas as pd, numpy as np, os
p=os.path.join(RAW, "CHARLS/2021/Health_Status_and_Functioning.dta")
cols=["ID","da011_15_","da012_15_","da015_15_","da016_15_","da016_1_15_","da016_2_15_",
      "da021_15_","da069","da071","da072","da074","da066","da067","da068",
      "db001","db003","db052","da078","da079","proxy"]
df,_=pyreadstat.read_dta(p,usecols=cols)
L=[]
L.append(f"rows={len(df)}")
for c in cols[1:]:
    s=df[c]; nm=int(s.notna().sum())
    vc=s.value_counts(dropna=True).head(8)
    L.append(f"{c:14s} n={nm:6d}  " + ", ".join(f"{k}:{v}" for k,v in vc.items()))
# 交叉：新诊断 × 治疗
sub=df[df["da011_15_"]==1]
L.append(f"\n新诊断关节炎 n={len(sub)}")
L.append("  da021_15_ (治疗方式) 分布: " + str(sub["da021_15_"].value_counts(dropna=False).head(8).to_dict()))
L.append("  有疼痛(da069==1) 比例: " + str(round((sub["da069"]==1).mean()*100,1)) + "%")
pain=df[df["da069"]==1]
L.append(f"\n有疼痛 n={len(pain)}")
L.append("  其中新诊断关节炎 %: " + str(round((pain["da011_15_"]==1).mean()*100,1)))
L.append("  其中在用止痛药(非处方) %: " + str(round((pain["da071"]==1).mean()*100,1)))
L.append("  其中在用阿片类 %: " + str(round((pain["da072"]==1).mean()*100,1)))
open(os.path.join(OUT,"11_w6_probe.txt"),"w",encoding="utf-8").write("\n".join(L))
print("\n".join(L))
