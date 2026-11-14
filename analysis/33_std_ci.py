# -*- coding: utf-8 -*-
"""
33: P9 —— 年龄标准化率的 95% CI（整群 bootstrap）
设计：按「群」重抽（CHARLS=communityID；其余=个体首次观测的家户 ID），
      B=1000；每次重算各年龄组的加权率再按 WHO 2000–2025 标准人口加权。
      CI 取 bootstrap 百分位（2.5 / 97.5）。
输出：33_std_ci.csv / 33_std_ci.md
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
B=1000
RNG=np.random.default_rng(20261114)

STD_W={1:5370,2:4550,3:3720,4:2970,5:2210,6:1520,7:910,8:630}
EDGES=[50,55,60,65,70,75,80,85,200]
def agegrp(s): return pd.cut(s,bins=EDGES,right=False,labels=list(range(1,9))).astype(float)
def anypos(df,cols):
    cols=[c for c in cols if c in df.columns]
    if not cols: return pd.Series(np.nan,index=df.index)
    sub=df[cols]; out=(sub==1).any(axis=1).astype(float)
    out[sub.isna().all(axis=1)]=np.nan; return out

W=pd.read_csv(os.path.join(OUT,"16_weights.csv"),dtype={"ID":str})[["cohort","ID","wave","weight"]]

SPEC={
 "CHARLS":dict(p=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",clu="communityID",
   base=["ID","communityID","wave","age","arthre","rxarthr_c"],
   sites=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
   pain_mode="list",dx="arthre",tx="rxarthr_c",age="age",waves=[1,2,3,4]),
 "ELSA":dict(p=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",clu="hhidc",
   base=["idauniqc","hhidc","wave","agey","arthre","painfr"],
   sites=[],pain_mode="single",pain="painfr",dx="arthre",tx=None,age="agey",waves=None),
 "HRS":dict(p=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",clu="hhhid",
   base=["hhidpn","hhhid","wave","ragey_b","arthre","rxarthr","painfr"],
   sites=[],pain_mode="single",pain="painfr",dx="arthre",tx="rxarthr",age="ragey_b",waves=None),
 "KLoSA":dict(p=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",clu="hhidc",
   base=["pid","hhidc","wave","agey","arthre","rxarthr"],
   sites=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
   pain_mode="list",dx="arthre",tx="rxarthr",age="agey",waves=None),
 "LASI":dict(p=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",clu="hhid",
   base=["prim_key","hhid","r1agey","r1arthre","r1rxarthr","r1jointp"],
   sites=[],pain_mode="joint",pain="r1jointp",dx="r1arthre",tx="r1rxarthr",age="r1agey",waves=None),
 "MHAS":dict(p=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",clu="hhidc",
   base=["rahhidnp","hhidc","wave","agey","arthre","rxarthr","painfr"],
   sites=[],pain_mode="single",pain="painfr",dx="arthre",tx="rxarthr",age="agey",waves=None),
 "SHARE":dict(p=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",clu="hhidc",
   base=["mergeid","hhidc","wave","agey","arthre","rxpain","painfr","pain_s"],
   sites=[],pain_mode="single",pain="painfr",dx="arthre",tx="rxpain",age="agey",waves=None),
}

def load(c,spec,tier):
    cols=list(dict.fromkeys(spec["base"]+spec["sites"]))
    df,_=pyreadstat.read_dta(spec["p"],usecols=cols)
    if spec["idc"]!="ID": df=df.rename(columns={spec["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if spec["waves"] is not None: d=d[d["wave"].isin(spec["waves"])]
    if tier=="B":
        if c=="SHARE": d=d[d["wave"].isin([1,2,4])]
    else:
        if c=="SHARE": d=d[d["wave"].isin([5,6,7,8])]
    d["ID"]=d["ID"].astype(str); d["CLU"]=d[spec["clu"]].astype(str)
    d["age"]=pd.to_numeric(d[spec["age"]],errors="coerce")
    d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    d["ageg"]=agegrp(d["age"])
    if tier=="B":
        if c=="LASI": pv=pd.to_numeric(d["r1jointp"],errors="coerce")
        elif spec["pain_mode"]=="list":
            sub=d[spec["sites"]]; pv=(sub==1).any(axis=1).astype(float); pv[sub.isna().all(axis=1)]=np.nan
        else: pv=pd.to_numeric(d["pain_s"],errors="coerce")
    else:
        if spec["pain_mode"]=="list":
            sub=d[spec["sites"]]; pv=(sub==1).any(axis=1).astype(float); pv[sub.isna().all(axis=1)]=np.nan
        else: pv=pd.to_numeric(d[spec["pain"]],errors="coerce")
    d["pain"]=np.where(pv.notna(),np.where(pv>0,1.0,0.0),np.nan)
    d["dx"]=pd.to_numeric(d[spec["dx"]],errors="coerce").where(lambda s:s.isin([0,1]))
    d["tx"]=pd.to_numeric(d[spec["tx"]],errors="coerce").where(lambda s:s.isin([0,1])) if spec["tx"] else np.nan
    d["wave"]=d["wave"].astype(float)
    d=d.merge(W[W.cohort==c][["ID","wave","weight"]],on=["ID","wave"],how="left")
    # 群固定：个体首次观测的家户
    first=(d.sort_values(["ID","wave"]).groupby("ID")["CLU"].first().rename("CLU_FIX"))
    d=d.merge(first,on="ID",how="left")
    return d[["cohort" if False else "ID","CLU_FIX","wave","ageg","pain","dx","tx","weight"]]

def boot_std(d,y,denom,cluster_col="CLU_FIX",B=B):
    """d 已含 y/denom 掩码；返回 (点估计, lo, hi, n, n_clusters)"""
    m=d[denom].dropna(subset=[y,"ageg","weight",cluster_col]).copy()
    if len(m)==0: return (np.nan,)*3+(0,0)
    # 预聚合到 (cluster, ageg)
    m["_yw"]=m[y]*m["weight"]
    g=m.groupby([cluster_col,"ageg"]).agg(w=("weight","sum"),yw=("_yw","sum"),n=(y,"size")).reset_index()
    clus=g[cluster_col].unique(); G=sorted(m["ageg"].dropna().unique())
    ci={c:i for i,c in enumerate(clus)}; gi={a:i for i,a in enumerate(G)}
    Wm=np.zeros((len(clus),len(G))); YWm=np.zeros_like(Wm)
    Wm[g[cluster_col].map(ci).values, g["ageg"].map(gi).values]=g["w"].values
    YWm[g[cluster_col].map(ci).values, g["ageg"].map(gi).values]=g["yw"].values
    wg=np.array([STD_W[a] for a in G],float)
    def std_from(Wm_,YWm_):
        num=YWm_.sum(0); den=Wm_.sum(0)
        ok=den>0
        if not ok.any(): return np.nan
        r=np.divide(num,den,out=np.zeros_like(num),where=ok)
        return np.sum(wg[ok]*r[ok])/np.sum(wg[ok])
    point=std_from(Wm,YWm)
    C=len(clus); reps=np.empty(B)
    for b in range(B):
        idx=RNG.integers(0,C,C)
        reps[b]=std_from(Wm[idx,:],YWm[idx,:])
    lo,hi=np.nanpercentile(reps,[2.5,97.5])
    return point,lo,hi,len(m),C

rows=[]
for tier in ["A","B"]:
    for c in (["HRS","MHAS","ELSA","SHARE"] if tier=="A" else ["CHARLS","KLoSA","LASI","SHARE"]):
        d=load(c,SPEC[c],tier)
        # 1. 疼痛率
        d["_one"]=1.0
        p,lo,hi,n,nc=boot_std(d,"pain",d["_one"].eq(1))
        rows.append(dict(tier=tier,cohort=c,metric="pain_rate",est=p*100,lo=lo*100,hi=hi*100,n=n,n_clusters=nc))
        # 2. 诊断覆盖（疼痛者中）
        q=d[(d["pain"]==1)].copy(); q["_one"]=1.0
        p,lo,hi,n,nc=boot_std(q,"dx",q["_one"].eq(1))
        rows.append(dict(tier=tier,cohort=c,metric="dx_given_pain",est=p*100,lo=lo*100,hi=hi*100,n=n,n_clusters=nc))
        # 3. 治疗覆盖（诊断且有疼痛）
        if d["tx"].notna().any():
            q=d[(d["pain"]==1)&(d["dx"]==1)].copy(); q["_one"]=1.0
            p,lo,hi,n,nc=boot_std(q,"tx",q["_one"].eq(1))
            rows.append(dict(tier=tier,cohort=c,metric="tx_given_dx_pain",est=p*100,lo=lo*100,hi=hi*100,n=n,n_clusters=nc))

res=pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"33_std_ci.csv"),index=False,encoding="utf-8-sig")
print(res.to_string(index=False))
