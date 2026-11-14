# -*- coding: utf-8 -*-
"""
41: 疼痛缺失是否与关节炎状态相关（选择偏倚的方向与量级）
对每个库-面板，比较「纳入」与「疼痛缺失」两组的关节炎率，
并给出疼痛率与诊断覆盖率的上下界。
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
L=[]; w=lambda s:(L.append(str(s)),print(s,flush=True))
W=pd.read_csv(os.path.join(OUT,"16_weights.csv"),dtype={"ID":str})[["cohort","ID","wave","weight"]]
def num(s): return pd.to_numeric(s,errors="coerce")

SPEC={
"CHARLS":dict(p=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",age="age",
  cols=["ID","wave","age","arthre"],sites=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
  wB=[1,2,3,4],modeB="list"),
"ELSA":dict(p=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",age="agey",
  cols=["idauniqc","wave","agey","arthre","painfr"],sites=[],wA=None,modeB="single",painA="painfr"),
"HRS":dict(p=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",age="ragey_b",
  cols=["hhidpn","wave","ragey_b","arthre","painfr"],sites=[],wA=None,modeB="single",painA="painfr"),
"KLoSA":dict(p=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",age="agey",
  cols=["pid","wave","agey","arthre"],sites=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
  wB=None,modeB="list"),
"LASI":dict(p=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",age="r1agey",
  cols=["prim_key","r1agey","r1arthre","r1jointp"],sites=[],wB=None,modeB="joint",painB="r1jointp"),
"MHAS":dict(p=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",age="agey",
  cols=["rahhidnp","wave","agey","arthre","painfr"],sites=[],wA=None,modeB="single",painA="painfr"),
"SHARE":dict(p=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",age="agey",
  cols=["mergeid","wave","agey","arthre","painfr","pain_s"],sites=[],wA=[5,6,7,8],wB=[1,2,4],modeB="single",
  painA="painfr",painB="pain_s"),
}
rows=[]
for c,s in SPEC.items():
    cols=list(dict.fromkeys(s["cols"]+s["sites"]))
    df,_=pyreadstat.read_dta(s["p"],usecols=cols)
    if s["idc"]!="ID": df=df.rename(columns={s["idc"]:"ID"})
    df["ID"]=df["ID"].astype(str)
    if "wave" not in df.columns: df["wave"]=1
    df["age"]=num(df[s["age"]])
    base=df[(df.age>=50)&(df.age<=105)].copy()
    dxcol="r1arthre" if c=="LASI" else "arthre"
    for tier in ["A","B"]:
        if tier=="A" and c in ["CHARLS","KLoSA","LASI"]: continue
        if tier=="B" and c in ["ELSA","HRS","MHAS"]: continue
        d=base.copy()
        wf=s.get("wA") if tier=="A" else s.get("wB")
        if wf is not None: d=d[d.wave.isin(wf)]
        if tier=="B":
            if s["modeB"]=="list":
                sub=d[s["sites"]]; pv=(sub==1).any(axis=1).astype(float); pv[sub.isna().all(axis=1)]=np.nan
                d["pain"]=pv
            else:
                pv=num(d[s["painB"]]); d["pain"]=pv.where(pv.isin([0,1]))
        else:
            pv=num(d[s["painA"]]); d["pain"]=pv.where(pv.isin([0,1]))
        d["dx"]=num(d[dxcol]).where(lambda x:x.isin([0,1]))
        d["wave"]=d.wave.astype(float)
        d=d.merge(W[W.cohort==c][["ID","wave","weight"]],on=["ID","wave"],how="left")
        e=d[d.weight.notna()]
        inc=e[e.pain.notna()]; exc=e[e.pain.isna()]
        if len(inc)==0: continue
        k=int((inc.pain==1).sum()); n=len(inc)
        n_exc=len(exc)
        lo=(k+0)/(n+n_exc)*100; hi=(k+n_exc)/(n+n_exc)*100
        rows.append(dict(cohort=c,tier=tier,n_inc=n,n_exc=n_exc,
                         dx_inc=round(inc["dx"].mean()*100,1) if inc["dx"].notna().any() else np.nan,
                         dx_exc=round(exc["dx"].mean()*100,1) if (n_exc and exc["dx"].notna().any()) else np.nan,
                         painrate=round(k/n*100,1),
                         bound=f"{lo:.1f}–{hi:.1f}",
                         gap=round((inc["dx"].mean()-exc["dx"].mean())*100,1) if (n_exc and exc["dx"].notna().any()) else np.nan))
r=pd.DataFrame(rows).sort_values("gap",ascending=False)
r.to_csv(os.path.join(OUT,"41_pain_missing_bias.csv"),index=False,encoding="utf-8-sig")
w(r.to_string(index=False))
w("\ndx_inc/dx_exc = 纳入组与疼痛缺失组的关节炎率（%）；gap = 两者之差")
w("bound = 把缺失者全按无痛/全按有痛补回时，疼痛率的取值区间")
w("\n→ gap 越大，说明疼痛缺失与关节炎状态越不独立，选择偏倚的方向越明确。")
open(os.path.join(OUT,"41_pain_missing_bias.txt"),"w",encoding="utf-8").write("\n".join(L))
print("\nWROTE 41_pain_missing_bias.txt")
