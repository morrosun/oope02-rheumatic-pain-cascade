# -*- coding: utf-8 -*-
"""
37: 四状态描述表的修正版 —— 把「无痛但已诊断」从「无痛」里拆出来
输出：37_four_state_v2.csv / 37_four_state_v2.md
说明：加权（harmonized 人-波权重），未做年龄标化（与原 17_four_state 口径一致，
      便于逐项对照；标化版见 31 的回归结果）。
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
def num(s): return pd.to_numeric(s,errors="coerce")
W=pd.read_csv(os.path.join(OUT,"16_weights.csv"),dtype={"ID":str})[["cohort","ID","wave","weight"]]

SPEC={
"CHARLS":dict(p=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",
  cols=["ID","wave","age","arthre","rxarthr_c","adlab_c",
        "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
  mode="list",dx="arthre",tx="rxarthr_c",age="age",adl="adlab_c",waves=[1,2,3,4]),
"ELSA":dict(p=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",
  cols=["idauniqc","wave","agey","arthre","painfr","adltot6"],
  mode="single",pain="painfr",dx="arthre",tx=None,age="agey",adl="adltot6",waves=None),
"HRS":dict(p=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",
  cols=["hhidpn","wave","ragey_b","arthre","rxarthr","painfr","adl6a"],
  mode="single",pain="painfr",dx="arthre",tx="rxarthr",age="ragey_b",adl="adl6a",waves=None),
"KLoSA":dict(p=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",
  cols=["pid","wave","agey","arthre","rxarthr",
        "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13",
        "dressb","bathb","eatb","toiltb","bedb_k","urinb"],
  mode="list",dx="arthre",tx="rxarthr",age="agey",adl="ADL6",waves=None),
"LASI":dict(p=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",
  cols=["prim_key","r1agey","r1arthre","r1rxarthr","r1jointp","r1adltot6"],
  mode="joint",pain="r1jointp",dx="r1arthre",tx="r1rxarthr",age="r1agey",adl="r1adltot6",waves=None),
"MHAS":dict(p=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",
  cols=["rahhidnp","wave","agey","arthre","rxarthr","painfr","adltot6"],
  mode="single",pain="painfr",dx="arthre",tx="rxarthr",age="agey",adl="adltot6",waves=None),
"SHARE":dict(p=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",
  cols=["mergeid","wave","agey","arthre","rxpain","painfr","pain_s","painlv",
        "walkra","dressa","batha","eata","beda","toilta"],
  mode="single",pain="painfr",dx="arthre",tx="rxpain",age="agey",adl="ADL6",waves=[5,6,7,8]),
}
def load(c,s):
    df,_=pyreadstat.read_dta(s["p"],usecols=s["cols"])
    if s["idc"]!="ID": df=df.rename(columns={s["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if s["waves"] is not None: d=d[d["wave"].isin(s["waves"])]
    d["ID"]=d["ID"].astype(str)
    d["age"]=num(d[s["age"]]); d=d[(d.age>=50)&(d.age<=105)].copy()
    if s["mode"]=="list":
        pc=[x for x in s["cols"] if x.startswith("da042s") or x.startswith("pain_")]
        pv=(d[pc]==1).any(axis=1).astype(float); pv[d[pc].isna().all(axis=1)]=np.nan
        d["pain"]=pv
    else:
        pv=num(d[s["pain"]]); d["pain"]=pv.where(pv.isin([0,1]))
    d["dx"]=num(d[s["dx"]]).where(lambda x:x.isin([0,1]))
    d["tx"]=num(d[s["tx"]]).where(lambda x:x.isin([0,1])) if s["tx"] else np.nan
    if s["adl"]=="ADL6":
        ac=["dressb","bathb","eatb","toiltb","bedb_k","urinb"] if c=="KLoSA" else ["walkra","dressa","batha","eata","beda","toilta"]
        d["adl"]=(d[ac]==1).any(axis=1).astype(float); d.loc[d[ac].isna().all(axis=1),"adl"]=np.nan
    else:
        a=num(d[s["adl"]]); d["adl"]=np.where(a.notna(),(a>=1).astype(float),np.nan)
    d["wave"]=d["wave"].astype(float)
    d=d.merge(W[W.cohort==c][["ID","wave","weight"]],on=["ID","wave"],how="left")
    return d

def state(d,has_tx):
    if has_tx:
        s=np.select([(d.pain==0)&(d.dx==0),(d.pain==1)&(d.dx==0),(d.pain==1)&(d.dx==1)&(d.tx==0),(d.pain==1)&(d.dx==1)&(d.tx==1),(d.pain==0)&(d.dx==1)],
                    ["0_no_pain_undx","1_pain_undx","2_pain_dx_untx","3_pain_dx_tx","4_no_pain_dx"],default=None)
    else:
        s=np.select([(d.pain==0)&(d.dx==0),(d.pain==1)&(d.dx==0),(d.pain==1)&(d.dx==1),(d.pain==0)&(d.dx==1)],
                    ["0_no_pain_undx","1_pain_undx","2_pain_dx","4_no_pain_dx"],default=None)
    return s

def wrate(g):
    m=g[["adl","weight"]].dropna()
    if len(m)==0 or m["weight"].sum()==0: return np.nan
    return np.average(m["adl"],weights=m["weight"])*100

rows=[]
for c,s in SPEC.items():
    d=load(c,s); d=d[d["weight"].notna()&d["pain"].notna()].copy()
    has_tx=bool(s["tx"]) and c!="SHARE"
    d["state"]=state(d,has_tx)
    dd=d.dropna(subset=["state"])
    for st,g in dd.groupby("state"):
        rows.append(dict(cohort=c,state=st,n=len(g),adl_pct=round(wrate(g),1)))
    # 旧口径：无痛（含已诊断）
    old=d[d["pain"]==0]
    rows.append(dict(cohort=c,state="OLD_no_pain_combined",n=len(old),adl_pct=round(wrate(old),1)))
res=pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"37_four_state_v2.csv"),index=False,encoding="utf-8-sig")
print(res.pivot_table(index="cohort",columns="state",values="adl_pct",aggfunc="first").to_string())
print()
print(res.pivot_table(index="cohort",columns="state",values="n",aggfunc="first").to_string())
