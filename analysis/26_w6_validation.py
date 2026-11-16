# -*- coding: utf-8 -*-
"""
26: 用 CHARLS 第六轮 da011_15_（直接问"自上次访问以来是否新诊断"）独立验证 A4
设计：基线 w5(2020) 有疼痛且未诊断 → W6(2021-2023) 是否新诊断 → W6 ADL（调整基线 ADL）
输出：26_w6_validation.txt / .csv
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, warnings
import numpy as np, pandas as pd, pyreadstat, statsmodels.api as sm
warnings.filterwarnings("ignore")
BASE = os.environ.get("AGING_DB_RAW", os.path.join(ROOT, "data_raw"))
CH=os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta")
def num(s): return pd.to_numeric(s,errors="coerce")

L=[]
def w(s=""):
    L.append(str(s)); print(s,flush=True)

# ---------- 1) 基线 w5 ----------
SITES=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10",
       "da042s11","da042s12","da042s13","da042s14","da042s15"]
b,_=pyreadstat.read_dta(CH,usecols=["ID","wave","age","ragender","raeducl","hatotfa","hrural","arthre","adlab_c",
                                    "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"]+SITES)
b=b[b["wave"]==5].copy()
b["ID"]=b["ID"].astype(str)
b["pain5"]=(b[SITES]==1).sum(axis=1).astype(float)
b.loc[b[SITES].isna().all(axis=1),"pain5"]=np.nan
b["nsite5"]=b["pain5"]
b["pain5"]=np.where(b["pain5"].notna(),(b["pain5"]>0).astype(float),np.nan)
b["dx5"]=num(b["arthre"]).where(lambda s:s.isin([0,1]))
b["adl5"]=np.where(num(b["adlab_c"]).notna(),(num(b["adlab_c"])>=1).astype(float),np.nan)
b["age5"]=num(b["age"]); b["male"]=num(b["ragender"])
ed=num(b["raeducl"]); b["edu_low"]=(ed==1).astype(float); b.loc[ed.isna(),"edu_low"]=np.nan
b["edu_mid"]=(ed==2).astype(float); b.loc[ed.isna(),"edu_mid"]=np.nan
b["rural"]=num(b["hrural"])
cc=[x for x in ["hibpe","diabe","hearte","stroke","cancre","lunge"] if x in b.columns]
b["comorb"]=(b[cc]==1).sum(axis=1).astype(float); b.loc[b[cc].isna().sum(axis=1)>3,"comorb"]=np.nan
b["smoke"]=num(b["smoken"])
b["wq"]=np.nan
if num(b["hatotfa"]).notna().sum()>200:
    b["wq"]=pd.qcut(num(b["hatotfa"]).rank(method="first"),5,labels=False)+1
w(f"w5 基线 {len(b)} 人；有疼痛 {int((b['pain5']==1).sum())}；未诊断 {int((b['dx5']==0).sum())}")

# ---------- 2) W6 ----------
cols=["ID","da011_15_","da012_15_","da021_15_","da069","proxy",
      "db001","db003","db005","db007","db009","db011","db012","db014","db016","db018","db020","db022","db052"]
fr=[]
for y in [2021,2022,2023]:
    d,_=pyreadstat.read_dta(os.path.join(BASE,str(y),"Health_Status_and_Functioning.dta"),usecols=cols)
    d["res_wave"]=y; fr.append(d)
w6=pd.concat(fr,ignore_index=True)
w6["ID"]=w6["ID"].astype(str)
w6["_nm"]=w6[cols[1:]].notna().sum(axis=1)
w6=w6.sort_values(["ID","res_wave","_nm"],ascending=[True,False,False]).drop_duplicates("ID",keep="first")
# Follow-up ADL uses the same six items as baseline (adlab_c)
ADL6=["db001","db003","db005","db007","db009","db011"]
IADL6=["db012","db014","db016","db018","db020","db022"]
w6["adl6"]=(w6[ADL6]>=2).any(axis=1).astype(float)
w6.loc[w6[ADL6].isna().all(axis=1),"adl6"]=np.nan
# 原 13 项版（ADL6+IADL6+行走）保留作敏感性
adl_cols=ADL6+IADL6+["db052"]
w6["adl13"]=(w6[adl_cols]>=2).any(axis=1).astype(float)
w6.loc[w6[adl_cols].isna().all(axis=1),"adl13"]=np.nan
w6["newdx"]=num(w6["da011_15_"]).map({1:1,2:0})
w6["pain6"]=np.where(num(w6["da069"]).between(2,5),1,np.where(num(w6["da069"])==1,0,np.nan))
w(f"W6 去重 {len(w6)} 人；新诊断 {int((w6['newdx']==1).sum())}；da011 有值 {int(w6['newdx'].notna().sum())}")

# ---------- 3) 合并 ----------
m=b.merge(w6[["ID","res_wave","newdx","pain6","adl6","adl13"]],on="ID",how="inner")
w(f"\n与 W6 匹配 {len(m)} 人")
w("基线 dx5 × W6 da011 有值：")
w(pd.crosstab(m["dx5"].fillna(-1), m["newdx"].notna()).to_string())

s=m[(m["dx5"]==0)&(m["pain5"]==1)&m["newdx"].notna()].copy()
w(f"\n分析样本（w5 有疼痛且未诊断，W6 有 da011）：{len(s)}")
rate=s["newdx"].mean()
w(f"两年内新诊断率：{rate*100:.1f}%（{int(s['newdx'].sum())}/{len(s)}）")

# ---------- 4) 模型 ----------
s["age10"]=(s["age5"]-50)/10
sd=s["nsite5"].std(); mu=s["nsite5"].mean(); s["sev_z"]=(s["nsite5"]-mu)/sd
base=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]

def fit(d,y,terms,tag):
    dd=d.dropna(subset=[y]+terms)
    if len(dd)<100: return
    X=sm.add_constant(dd[terms],has_constant="add")
    try:
        r=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd["ID"]})
    except Exception as e:
        w(f"  拟合失败 {tag}: {e}"); return
    w(f"\n[{tag}] n={len(dd)}  pseudoR2={r.prsquared:.3f}")
    for t in terms:
        if t not in r.params: continue
        b0,se=r.params[t],r.bse[t]
        w(f"   {t:9s} OR={np.exp(b0):7.3f} ({np.exp(b0-1.96*se):.3f}-{np.exp(b0+1.96*se):.3f})  p={r.pvalues[t]:.4g}")
    return r

# 验证一：谁会新诊断（应复现 男性<1、年龄>1、教育/财富≈1）
fit(s,"newdx",[t for t in base if s[t].notna().sum()>100],"M1 新诊断的决定因素")

# 验证二：新诊断 → W6 ADL（调整基线 ADL 与疼痛部位数）
q=s.dropna(subset=["adl6","adl5"]).copy()
fit(q,"adl6",["newdx","adl5"]+[t for t in base+["sev_z"] if q[t].notna().sum()>100],"M2 新诊断 → W6 ADL（6 项，与基线同构念）")
# 敏感性：13 项（ADL6+IADL6+行走）
q2=s.dropna(subset=["adl13","adl5"]).copy()
fit(q2,"adl13",["newdx","adl5"]+[t for t in base+["sev_z"] if q2[t].notna().sum()>100],"M2b 敏感性：13 项结局")

# 与 A4 的 CHARLS 估计对比
a4=pd.read_csv(os.path.join(OUT,"24_incident.csv"))
c=a4[(a4.cohort=="CHARLS")&(a4.term=="y_dx")]
if len(c):
    c=c.iloc[0]
    w(f"\nA4（0→1 推断）CHARLS：OR={c['OR']:.3f} ({c['lo']:.3f}-{c['hi']:.3f}), n={int(c['n'])}")
w(f"本次（直接问新诊断）CHARLS：见上方 M2 的 newdx")

with open(os.path.join(OUT,"26_w6_validation.txt"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("\nWROTE 26_w6_validation.txt")
