# -*- coding: utf-8 -*-
"""
27: 把 CHARLS 第六轮（2021–2023）并入第 1 步的加权级联（作为"任何疼痛"层的对照行）
输出：27_charls_w6_row.csv / 27_compare.md
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, warnings
import numpy as np, pandas as pd, pyreadstat
warnings.filterwarnings("ignore")
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
BASE = os.environ.get("AGING_DB_RAW", os.path.join(ROOT, "data_raw"))
CH=os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta")
def num(s): return pd.to_numeric(s,errors="coerce")

STD_W={1:5370,2:4550,3:3720,4:2970,5:2210,6:1520,7:910,8:630}
EDGES=[50,55,60,65,70,75,80,85,200]
def ageg(s): return pd.cut(s,bins=EDGES,right=False,labels=list(range(1,9))).astype(float)

# ---------- 派生集（14 已生成） ----------
d=pd.read_csv(os.path.join(OUT,"14_w6_derived.csv"),dtype={"ID":str})
print("W6 派生集",d.shape, flush=True)

# ---------- 权重 ----------
wt,_=pyreadstat.read_dta(os.path.join(BASE,"2023","Weights.dta"),usecols=["ID","INDV_weight_post2"])
wt["ID"]=wt["ID"].astype(str); wt=wt.rename(columns={"INDV_weight_post2":"weight"})
d=d.merge(wt,on="ID",how="left")

# ---------- 年龄：w5 年龄 + 年差；缺失用出生年 ----------
w5,_=pyreadstat.read_dta(CH,usecols=["ID","wave","age"])
w5=w5[w5["wave"]==5][["ID","age"]].rename(columns={"age":"age5"}); w5["ID"]=w5["ID"].astype(str)
d=d.merge(w5,on="ID",how="left")
d["age_est"]=num(d["age5"])+(d["res_wave"]-2020)

by=[]
for y in [2021,2022,2023]:
    p=os.path.join(BASE,str(y),"Demographics.dta")
    if not os.path.exists(p): continue
    x,_=pyreadstat.read_dta(p,usecols=["ID","ba004_1"])
    x["ID"]=x["ID"].astype(str); x["byear"]=num(x["ba004_1"]); x["res_wave"]=y
    by.append(x[["ID","byear"]])
by=pd.concat(by).drop_duplicates("ID")
d=d.merge(by,on="ID",how="left")
d["age_est"]=d["age_est"].where(d["age_est"].notna(), d["res_wave"]-d["byear"])
print("年龄可得：",d["age_est"].notna().sum(),"/",len(d), flush=True)

d=d[(d["age_est"]>=50)&(d["age_est"]<=105)].copy()
d["age"]=d["age_est"]; d["ageg"]=ageg(d["age"])
d["dx"]=num(d["ever_w6"]); d["tx"]=num(d["treat"]); d["adl"]=num(d["adl_any"])
d["pain_any"]=num(d["pain_any"]); d["pain_sev"]=num(d["pain_sev"])
print("分析样本",len(d),"有权重",int(d["weight"].notna().sum()), flush=True)

def wrate(df,y):
    m=df[[y,"weight"]].dropna()
    if len(m)==0 or m["weight"].sum()==0: return np.nan
    return np.average(m[y],weights=m["weight"])
def std(df,y):
    dd=df[[y,"weight","ageg"]].dropna(); num_=0; den=0
    for g,wg in STD_W.items():
        s=dd[dd["ageg"]==g]
        if len(s)==0 or s["weight"].sum()==0: continue
        num_+=wg*np.average(s[y],weights=s["weight"]); den+=wg
    return num_/den if den else np.nan

rows=[]
def add(label,painvar):
    s=d[d[painvar].notna()].copy()
    p=s[s[painvar]==1]
    dxp=p[p["dx"].notna()]
    txp=dxp[dxp["dx"]==1]
    txp=txp[txp["tx"].notna()]
    rows.append(dict(定义=label,n=len(s),
        疼痛率=round(wrate(s,painvar)*100,1),疼痛率标化=round(std(s,painvar)*100,1),
        诊断率=round(wrate(dxp,"dx")*100,1),诊断率标化=round(std(dxp,"dx")*100,1),
        治疗率=round(wrate(txp,"tx")*100,1) if len(txp) else np.nan,
        治疗率标化=round(std(txp,"tx")*100,1) if len(txp) else np.nan,
        n_dx=len(dxp), n_tx=len(txp)))
add("体痛 ≥ 一点", "pain_any")
add("体痛 ≥ 相当多", "pain_sev")

res=pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"27_charls_w6_row.csv"),index=False,encoding="utf-8-sig")
print(res.to_string(index=False))

# 四状态 ADL
d["state"]=np.select(
  [d["pain_any"]==0,(d["pain_any"]==1)&(d["dx"]==0),(d["pain_any"]==1)&(d["dx"]==1)&(d["tx"]==0),(d["pain_any"]==1)&(d["dx"]==1)&(d["tx"]==1)],
  ["无痛","痛·未诊断","痛·诊断未治疗","痛·诊断在治疗"],default=None)
st=d.dropna(subset=["state"]).groupby("state").apply(
    lambda g: pd.Series({"n":len(g),"ADL受限%":round(wrate(g.assign(adl=g["adl"]),"adl")*100,1)}), include_groups=False).reset_index()
st.to_csv(os.path.join(OUT,"27_charls_w6_state.csv"),index=False,encoding="utf-8-sig")
print("\n四状态："); print(st.to_string(index=False))
