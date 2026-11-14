# -*- coding: utf-8 -*-
"""29: 复盘审计 —— 权重匹配、参照组污染、样本流失、口径一致性"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
L=[]; w=lambda s: (L.append(str(s)), print(s, flush=True))

W=pd.read_csv(os.path.join(OUT,"16_weights.csv"),dtype={"ID":str})

# ---------- 1. HRS wtresp 是否按波正确 ----------
w("\n" + "="*70)
w("【1】HRS working_data 的 wtresp 是否按波正确（每波非缺失率与中位数）")
w("="*70)
h,_=pyreadstat.read_dta(os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),
                        usecols=["hhidpn","wave","wtresp","ragey_b"])
g=h.groupby("wave").agg(rows=("wtresp","size"), nonmiss=("wtresp","count"),
                        med=("wtresp","median"), age_med=("ragey_b","median"))
g["pct"]=(g["nonmiss"]/g["rows"]*100).round(1)
w(g.to_string())

# ---------- 2. 各库权重匹配率（以 17 的构建口径） ----------
w("\n" + "="*70)
w("【2】17_cascade 各库权重 merge 匹配率")
w("="*70)
CFG=[("CHARLS",os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),
      "ID",["ID","wave","age","arthre"],lambda d:d[d.wave.isin([1,2,3,4])]),
     ("ELSA",os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),"idauniqc",
      ["idauniqc","wave","agey","arthre"],lambda d:d),
     ("HRS",os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),"hhidpn",
      ["hhidpn","wave","ragey_b","arthre"],lambda d:d),
     ("KLoSA",os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),"pid",
      ["pid","wave","agey","arthre"],lambda d:d),
     ("LASI",os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),"prim_key",
      ["prim_key","r1agey","r1arthre"],lambda d:d),
     ("MHAS",os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),"rahhidnp",
      ["rahhidnp","wave","agey","arthre"],lambda d:d),
     ("SHARE",os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),"mergeid",
      ["mergeid","wave","agey","arthre"],lambda d:d)]
for c,p,idc,cols,filt in CFG:
    df,_=pyreadstat.read_dta(p,usecols=cols)
    if idc!="ID": df=df.rename(columns={idc:"ID"})
    d=filt(df.copy()); d["ID"]=d["ID"].astype(str)
    if "wave" not in d.columns: d["wave"]=1
    ww=W[W.cohort==c][["ID","wave","weight"]].copy(); ww["wave"]=ww["wave"].astype(float)
    d["wave"]=d["wave"].astype(float)
    m=d.merge(ww,on=["ID","wave"],how="left")
    hit=m["weight"].notna().sum()
    w(f"{c:7s} rows={len(d):7d}  merge匹配={hit:7d} ({hit/len(d)*100:5.1f}%)  唯一ID={d.ID.nunique():7d}  "
      f"重复(ID,wave)={int(d.duplicated(['ID','wave']).sum()):5d}")
    if hit/len(d) < 0.5:
        byw=d.assign(hit=m["weight"].notna().values).groupby("wave")["hit"].agg(["size","sum"])
        byw["pct"]=(byw["sum"]/byw["size"]*100).round(1)
        w("        逐波匹配率: " + byw["pct"].to_dict().__str__())

# ---------- 3. 参照组污染：无痛但已诊断 ----------
w("\n" + "="*70)
w("【3】A3 四状态参照组污染：pain==0 但 arthre==1 的比例")
w("="*70)
# 用 17 已生成的四状态无法看出，重建简化版
CH=pyreadstat.read_dta(os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),
                       usecols=["ID","wave","age","arthre","da042s12","da042s2","da042s3","da042s4","da042s5",
                                "da042s8","da042s9","da042s10","da042s11","da042s13","da042s14","da042s15"])[0]
CH=CH[(CH.age>=50)&(CH.age<=105)]
pcols=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"]
CH["pain"]=(CH[pcols]==1).any(axis=1).astype(float); CH.loc[CH[pcols].isna().all(axis=1),"pain"]=np.nan
nop=CH[CH.pain==0]
w(f"CHARLS 无痛者 n={len(nop)}；其中 arthre==1 = {int((nop.arthre==1).sum())} ({(nop.arthre==1).mean()*100:.1f}%)")

# ---------- 4. 各库 adl / 结局变量取值分布 ----------
w("\n" + "="*70)
w("【4】ADL 结局变量取值分布（确认 >=1 的二分是否等价）")
w("="*70)
ADL={"CHARLS":(os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),"adlab_c","ID"),
     "ELSA":(os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),"adltot6","idauniqc"),
     "HRS":(os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),"adl6a","hhidpn"),
     "LASI":(os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),"r1adltot6","prim_key"),
     "MHAS":(os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),"adltot6","rahhidnp")}
for c,(p,col,idc) in ADL.items():
    s=pyreadstat.read_dta(p,usecols=[col])[0][col]
    w(f"{c:7s} {col:10s} n={s.notna().sum():6d} 分布={s.value_counts().sort_index().head(8).to_dict()}  mean={s.mean():.2f}")

# ---------- 5. 模型样本流失 ----------
w("\n" + "="*70)
w("【5】A2 模型 dropna 后的样本流失（诊断模型：有疼痛者）")
w("="*70)
try:
    det=pd.read_csv(os.path.join(OUT,"19_determinants.csv"))
    d0=det[det.stage=="D_dx_given_pain"][["panel","cohort","n"]].drop_duplicates()
    w(d0.to_string(index=False))
except Exception as e:
    w(f"读取失败 {e}")

open(os.path.join(OUT,"29_audit.txt"),"w",encoding="utf-8").write("\n".join(L))
print("\nWROTE 29_audit.txt")
