# -*- coding: utf-8 -*-
"""
31: A3 修正版 —— 修复两处结构性问题
  (1) 参照组洁净：参照 = 无痛【且未诊断】；"无痛已诊断"单列为 s4，不再污染参照
  (2) 构念对齐：主模型用三状态（无痛/痛未诊断/痛已诊断），7 库全部可比；
      四状态（拆治疗）仅在 5 个有 tx 变量的库内做，不与三状态库合并
输出：31_outcome_v2.csv / 31_pooled_v2.csv / A3修正版.md
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, warnings
import numpy as np, pandas as pd, pyreadstat, statsmodels.api as sm
warnings.filterwarnings("ignore")
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
def num(s): return pd.to_numeric(s,errors="coerce")

def dl(est,se):
    est=np.array(est,float); se=np.array(se,float)
    ok=~np.isnan(est)&~np.isnan(se)&(se>0); est,se=est[ok],se[ok]; k=len(est)
    if k==0: return None
    w=1/se**2; fe=np.sum(w*est)/np.sum(w); Q=np.sum(w*(est-fe)**2); df=k-1
    C=np.sum(w)-np.sum(w**2)/np.sum(w); tau2=max(0.0,(Q-df)/C) if C>0 else 0.0
    wr=1/(se**2+tau2); mu=np.sum(wr*est)/np.sum(wr); sem=np.sqrt(1/np.sum(wr))
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    return dict(k=k,OR=np.exp(mu),lo=np.exp(mu-1.96*sem),hi=np.exp(mu+1.96*sem),I2=I2,Q=Q)

def fit(d,y,terms):
    dd=d.dropna(subset=[y]+terms)
    if len(dd)<200 or dd[y].nunique()<2: return None,len(dd)
    X=sm.add_constant(dd[terms],has_constant="add")
    try:
        m=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd["ID"]})
    except Exception:
        return None,len(dd)
    return {t:(m.params[t],m.bse[t]) for t in terms if t in m.params},len(dd)

COD={
"CHARLS":dict(path=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",
  cols=["ID","wave","age","ragender","raeducl","hatotfa","hrural","arthre","rxarthr_c","adlab_c",
        "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="list",dx="arthre",tx="rxarthr_c",age="age",adl="adlab_c",sev="nsite",rural="hrural",wealth="hatotfa",waves=[1,2,3,4]),
"ELSA":dict(path=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",
  cols=["idauniqc","wave","agey","ragender","raeducl","hatotb","arthre","painfr","painlv","adltot6",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",tx=None,age="agey",adl="adltot6",sev="painlv",rural=None,wealth="hatotb"),
"HRS":dict(path=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",
  cols=["hhidpn","wave","ragey_b","ragender","raeducl","atotn","rural","arthre","rxarthr","painfr","painlv","adl6a",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",tx="rxarthr",age="ragey_b",adl="adl6a",sev="painlv",rural="rural",wealth="atotn"),
"KLoSA":dict(path=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",
  cols=["pid","wave","agey","ragender","raeducl","atotf","rural","arthre","rxarthr",
        "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13","painhlv",
        "dressb","bathb","eatb","toiltb","bedb_k","urinb",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="list",dx="arthre",tx="rxarthr",age="agey",adl="ADL6",sev="painhlv",rural="rural",wealth="atotf"),
"LASI":dict(path=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",
  cols=["prim_key","r1agey","ragender","raeducl","hh1atotb","hh1rural","r1arthre","r1rxarthr","r1jointp","r1painfrq",
        "r1adltot6","r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge","r1smoken"],
  pain=None,joint="r1jointp",dx="r1arthre",tx="r1rxarthr",age="r1agey",adl="r1adltot6",
  sev="r1painfrq",rural="hh1rural",wealth="hh1atotb"),
"MHAS":dict(path=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",
  cols=["rahhidnp","wave","agey","ragender","raeducl","atotb","rural","arthre","rxarthr","painfr","painlv","adltot6",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",tx="rxarthr",age="agey",adl="adltot6",sev="painlv",rural="rural",wealth="atotb"),
"SHARE":dict(path=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",
  cols=["mergeid","wave","agey","ragender","raeducl","hhatotn","rural","arthre","rxpain","painfr","pain_s","painlv",
        "walkra","dressa","batha","eata","beda","toilta",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",tx="rxpain",age="agey",adl="ADL6",sev="painlv",rural="rural",wealth="hhatotn"),
}

def load(c,cfg,joint=None,waves=None):
    df,_=pyreadstat.read_dta(cfg["path"],usecols=cfg["cols"])
    if cfg["idc"]!="ID": df=df.rename(columns={cfg["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if waves is not None: d=d[d["wave"].isin(waves)]
    d["ID"]=d["ID"].astype(str)
    d["age"]=num(d[cfg["age"]]); d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    if joint:
        pv=num(d[joint]); d["pain"]=pv.where(pv.isin([0,1])); d["sev_raw"]=num(d[cfg["sev"]])
    elif cfg["pain"]=="list":
        pc=[x for x in cfg["cols"] if x.startswith("da042s") or x.startswith("pain_")]
        cnt=(d[pc]==1).sum(axis=1).astype(float); cnt[d[pc].isna().all(axis=1)]=np.nan
        d["pain"]=np.where(cnt.notna(),(cnt>0).astype(float),np.nan); d["sev_raw"]=cnt
    elif isinstance(cfg["pain"],str):
        pv=num(d[cfg["pain"]]); d["pain"]=pv.where(pv.isin([0,1]))
        d["sev_raw"]=num(d[cfg["sev"]]) if cfg["sev"] in d.columns else np.nan
    d["dx"]=num(d[cfg["dx"]]).where(lambda s:s.isin([0,1]))
    d["tx"]=num(d[cfg["tx"]]).where(lambda s:s.isin([0,1])) if cfg["tx"] else np.nan
    if cfg["adl"]=="ADL6":
        cols=["dressb","bathb","eatb","toiltb","bedb_k","urinb"] if c=="KLoSA" else ["walkra","dressa","batha","eata","beda","toilta"]
        d["adl"]=(d[cols]==1).any(axis=1).astype(float); d.loc[d[cols].isna().all(axis=1),"adl"]=np.nan
    else:
        a=num(d[cfg["adl"]]); d["adl"]=np.where(a.notna(),(a>=1).astype(float),np.nan)
    d["age10"]=(d["age"]-50)/10; d["male"]=num(d["ragender"])
    ed=num(d["raeducl"]); d["edu_low"]=np.where(ed.notna(),(ed==1).astype(float),np.nan)
    d["edu_mid"]=np.where(ed.notna(),(ed==2).astype(float),np.nan)
    d["rural"]=num(d[cfg["rural"]]) if cfg["rural"] else np.nan
    cc=[x for x in ["hibpe","diabe","hearte","stroke","cancre","lunge"] if x in d.columns]
    if c=="LASI": cc=["r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge"]
    d["comorb"]=(d[cc]==1).sum(axis=1).astype(float); d.loc[d[cc].isna().sum(axis=1)>len(cc)//2,"comorb"]=np.nan
    smk="r1smoken" if c=="LASI" else "smoken"
    d["smoke"]=num(d[smk]) if smk in d.columns else np.nan
    # 财富五分位：与 A2 对齐，按 cohort×wave
    d["wq"]=np.nan
    if cfg["wealth"] in d.columns:
        for _,g in d.groupby("wave"):
            s=num(g[cfg["wealth"]])
            if s.notna().sum()<200: continue
            d.loc[g.index,"wq"]=pd.qcut(s.rank(method="first"),5,labels=False)+1
    s=pd.Series(d["sev_raw"],dtype=float)
    mu=s[d["pain"]==1].mean(); sd=s[d["pain"]==1].std()
    d["sev_z"]=(s-mu)/sd if sd and sd>0 else np.nan
    d["cohort"]=c; d["has_tx"]=bool(cfg["tx"]) and (c!="SHARE")
    return d

D={}
for c,cfg in COD.items():
    D[c]=load(c,cfg,joint=cfg.get("joint"),waves=cfg.get("waves"))
SHB=load("SHARE",COD["SHARE"],joint="pain_s",waves=[1,2,4])

BASE=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
rows=[]
PANELS={"A":["HRS","MHAS","ELSA","SHARE"],"B":["CHARLS","KLoSA","LASI","SHARE"]}

def mk(d,mode):
    """mode: 'three' = 无痛未诊断/痛未诊断/痛已诊断/无痛已诊断
             'four'  = 无痛未诊断/痛未诊断/痛已诊断未治疗/痛已诊断在治疗/无痛已诊断"""
    if mode=="three":
        s=np.select([(d.pain==0)&(d.dx==0),(d.pain==1)&(d.dx==0),(d.pain==1)&(d.dx==1),(d.pain==0)&(d.dx==1)],
                    [0,1,2,3],default=np.nan)
    else:
        s=np.select([(d.pain==0)&(d.dx==0),(d.pain==1)&(d.dx==0),
                     (d.pain==1)&(d.dx==1)&(d.tx==0),(d.pain==1)&(d.dx==1)&(d.tx==1),(d.pain==0)&(d.dx==1)],
                    [0,1,2,3,4],default=np.nan)
    return s

for panel,cohorts in PANELS.items():
    for c in cohorts:
        d=SHB if (c=="SHARE" and panel=="B") else D[c]
        d=d.copy()
        terms=[t for t in BASE if d[t].notna().sum()>100]
        for mode in ["three","four"]:
            if mode=="four" and not bool(d["has_tx"].iloc[0]): continue
            d["state"]=mk(d,mode)
            dd=d.dropna(subset=["adl","state"]).copy()
            ks=[1,2,3] if mode=="three" else [1,2,3,4]
            for k in ks: dd[f"s{k}"]=(dd.state==k).astype(float)
            st=[f"s{k}" for k in ks if dd[f"s{k}"].sum()>50]
            tt=[t for t in terms if dd[t].notna().sum()>100]
            r,n=fit(dd,"adl",st+tt)
            if r:
                for t in st:
                    b,se=r[t]; rows.append(dict(panel=panel,cohort=c,mode=mode,term=t,n=n,
                        OR=round(np.exp(b),3),lo=round(np.exp(b-1.96*se),3),hi=round(np.exp(b+1.96*se),3),b=b,se=se))
        # M2：疼痛者中 已诊断 vs 未诊断
        sub=d[(d.pain==1)&d.dx.notna()].copy()
        for tag,extra in [("M2_no_sev",[]),("M2_adj_sev",["sev_z"] if (d.sev_z.notna().sum()>100) else [])]:
            tt=[t for t in terms+extra+["dx"] if sub[t].notna().sum()>100]
            if extra==[] and tag=="M2_adj_sev": continue
            if d.sev_z.notna().sum()<=100 and tag=="M2_adj_sev": continue
            r,n=fit(sub,"adl",tt)
            if r and "dx" in r:
                b,se=r["dx"]; rows.append(dict(panel=panel,cohort=c,mode=tag,term="dx_vs_undx",n=n,
                    OR=round(np.exp(b),3),lo=round(np.exp(b-1.96*se),3),hi=round(np.exp(b+1.96*se),3),b=b,se=se))

res=pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"31_outcome_v2.csv"),index=False,encoding="utf-8-sig")

pool=[]
for (panel,mode,term),g in res.groupby(["panel","mode","term"]):
    r=dl(g["b"].values,g["se"].values)
    if r is None: continue
    pool.append(dict(panel=panel,mode=mode,term=term,k=r["k"],OR=round(r["OR"],3),
                     lo=round(r["lo"],3),hi=round(r["hi"],3),I2=round(r["I2"],1),
                     cohorts=",".join(sorted(g["cohort"]))))
pooled=pd.DataFrame(pool)
pooled.to_csv(os.path.join(OUT,"31_pooled_v2.csv"),index=False,encoding="utf-8-sig")
print("=== per-cohort ===")
print(res[["panel","cohort","mode","term","n","OR","lo","hi"]].to_string(index=False))
print("\n=== pooled ===")
print(pooled.to_string(index=False))
