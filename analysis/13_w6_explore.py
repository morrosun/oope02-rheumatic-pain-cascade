# -*- coding: utf-8 -*-
"""Step 0-a: 探查 W6 三年包结构、权重、ID 格式与 da011×da012 跳转逻辑"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pyreadstat, pandas as pd, numpy as np

BASE = os.environ.get("AGING_DB_RAW", os.path.join(ROOT, "data_raw"))
L=[]
def w(s=""):
    L.append(str(s)); print(s, flush=True)

# 1) Sample_Infor
si,_ = pyreadstat.read_dta(os.path.join(BASE,"2023","Sample_Infor.dta"))
w("## Sample_Infor")
w(f"rows={len(si)}  uniqueID={si['ID'].nunique()}")
w("res_wave: " + str(si['res_wave'].value_counts().sort_index().to_dict()))
w("crosssection: " + str(si['crosssection'].value_counts().to_dict()))
w("died: " + str(si['died'].value_counts().to_dict()))
w("ID dtype=" + str(si['ID'].dtype) + "  head=" + str(si['ID'].head(3).tolist()))

# 2) Weights
wt,mt = pyreadstat.read_dta(os.path.join(BASE,"2023","Weights.dta"))
w("\n## Weights")
w(f"rows={len(wt)}  cols={mt.number_columns}")
for c,l in zip(mt.column_names, mt.column_labels):
    s=wt[c]
    desc=""
    if pd.api.types.is_numeric_dtype(s) and s.notna().any():
        desc = "min=%.3f med=%.3f max=%.3f" % (s.min(), s.median(), s.max())
    w(f"  {c:22s} {str(l)[:44]:46s} n={s.notna().sum():6d} dtype={s.dtype} {desc} " + str(s.dropna().head(3).tolist()))
w("  weight ID unique=" + str(wt['ID'].nunique()) if 'ID' in wt.columns else "  no ID")

# 3) 三年包 da011×da012 交叉
frames=[]
for y in ["2021","2022","2023"]:
    p=os.path.join(BASE,y,"Health_Status_and_Functioning.dta")
    cols=["ID","da011_15_","da012_15_","da015_15_","da016_1_15_","da016_2_15_",
          "da021_15_","da069","da071","da072","da074","proxy"]
    df,_=pyreadstat.read_dta(p,usecols=cols)
    df["res_wave"]=int(y)
    frames.append(df)
    w(f"\n## {y}: rows={len(df)} uniqueID={df['ID'].nunique()}")
    w("da011_15_ x da012_15_:")
    w(pd.crosstab(df["da011_15_"].fillna(-1), df["da012_15_"].fillna(-1)).to_string())
    w("da021_15_ by da011_15_:")
    w(pd.crosstab(df["da011_15_"].fillna(-1), df["da021_15_"].fillna(-1)).to_string())

allw=pd.concat(frames, ignore_index=True)
w("\n## 合并后（未去重）")
w(f"rows={len(allw)}  uniqueID={allw['ID'].nunique()}")
dup = allw[allw.duplicated('ID', keep=False)].sort_values('ID')
w(f"重复 ID 数={dup['ID'].nunique()}  记录数={len(dup)}")
w("重复者的 res_wave 组合（前 10 组）:")
w(dup.groupby('ID')['res_wave'].apply(lambda s: tuple(sorted(s))).value_counts().head(10).to_string())

# 4) ID 与成品层 charls.dta 的格式
ch,_ = pyreadstat.read_dta(os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta"),
                           usecols=["ID","wave","arthre","rxarthr"])
w("\n## charls.dta")
w(f"ID dtype={ch['ID'].dtype}  head={ch['ID'].head(3).tolist()}")
ch5 = ch[ch["wave"]==5][["ID","arthre"]].copy()
w(f"w5 rows={len(ch5)}  arthre: {ch5['arthre'].value_counts(dropna=False).to_dict()}")
ids_w6=set(allw['ID'].astype(str)); ids_ch=set(ch5['ID'].astype(str))
w(f"W6 并集 ID={len(ids_w6)}  charls w5 ID={len(ids_ch)}  匹配={len(ids_w6 & ids_ch)}  匹配率={len(ids_w6 & ids_ch)/len(ids_w6)*100:.1f}%")

open(os.path.join(OUT,"13_w6_explore.txt"),"w",encoding="utf-8").write("\n".join(L))
