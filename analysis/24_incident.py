# -*- coding: utf-8 -*-
"""
24: A4 —— 纵向级联：基线"有疼痛但未诊断"者，2 年内获得诊断的比例；
      以及获得诊断是否与 ADL 改善相关（滞后结局 + 基线 ADL 调整）
输出：24_incident.csv / 24_incident_pooled.csv / 纵向级联.md
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
    return dict(k=k,OR=np.exp(mu),lo=np.exp(mu-1.96*sem),hi=np.exp(mu+1.96*sem),I2=I2)

def fit(d,y,terms):
    dd=d.dropna(subset=[y]+terms)
    if len(dd)<150 or dd[y].nunique()<2: return None,len(dd)
    X=sm.add_constant(dd[terms],has_constant="add")
    try:
        m=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd["ID"]})
    except Exception:
        return None,len(dd)
    return {t:(m.params[t],m.bse[t]) for t in terms if t in m.params}, len(dd)

CFG={
"CHARLS":dict(p=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",
  cols=["ID","wave","iwy","age","ragender","raeducl","hatotfa","hrural","arthre","adlab_c",
        "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
  dx="arthre",adl="adlab_c",age="age",year="iwy",rural="hrural",wealth="hatotfa",waves=[1,2,3,4,5]),
"ELSA":dict(p=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",
  cols=["idauniqc","wave","iwindy","agey","ragender","raeducl","hatotb","arthre","painfr","painlv","adltot6",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",adl="adltot6",age="agey",year="iwindy",rural=None,wealth="hatotb"),
"HRS":dict(p=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",
  cols=["hhidpn","wave","iwendy","ragey_b","ragender","raeducl","atotn","rural","arthre","painfr","painlv","adl6a",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",adl="adl6a",age="ragey_b",year="iwendy",rural="rural",wealth="atotn"),
"KLoSA":dict(p=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",
  cols=["pid","wave","iwy","agey","ragender","raeducl","atotf","rural","arthre",
        "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13","painhlv",
        "dressb","bathb","eatb","toiltb","bedb_k","urinb","hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
  dx="arthre",adl="ADL6",age="agey",year="iwy",rural="rural",wealth="atotf"),
"MHAS":dict(p=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",
  cols=["rahhidnp","wave","iwy","agey","ragender","raeducl","atotb","rural","arthre","painfr","painlv","adltot6",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",adl="adltot6",age="agey",year="iwy",rural="rural",wealth="atotb"),
"SHARE":dict(p=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",
  cols=["mergeid","wave","iwy","agey","ragender","raeducl","hhatotn","rural","arthre","painfr","painlv",
        "walkra","dressa","batha","eata","beda","toilta","hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",adl="ADL6",age="agey",year="iwy",rural="rural",wealth="hhatotn",waves=[5,6,7,8]),
}

def load(c,cfg):
    df,_=pyreadstat.read_dta(cfg["p"],usecols=cfg["cols"])
    if cfg["idc"]!="ID": df=df.rename(columns={cfg["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if cfg.get("waves"): d=d[d["wave"].isin(cfg["waves"])]
    d["ID"]=d["ID"].astype(str)
    d["age"]=num(d[cfg["age"]]); d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    if isinstance(cfg["pain"],str):
        pv=num(d[cfg["pain"]]); d["pain"]=pv.where(pv.isin([0,1]))
        d["sev_raw"]=num(d[cfg["pain"].replace("painfr","painlv")]) if "painlv" in d.columns else np.nan
    else:
        pv=(d[cfg["pain"]]==1).sum(axis=1).astype(float); pv[d[cfg["pain"]].isna().all(axis=1)]=np.nan
        d["pain"]=np.where(pv.notna(),(pv>0).astype(float),np.nan); d["sev_raw"]=pv
    d["dx"]=num(d[cfg["dx"]]).where(lambda s:s.isin([0,1]))
    if cfg["adl"]=="ADL6":
        cols=["dressb","bathb","eatb","toiltb","bedb_k","urinb"] if c=="KLoSA" else ["walkra","dressa","batha","eata","beda","toilta"]
        d["adl"]=(d[cols]==1).any(axis=1).astype(float); d.loc[d[cols].isna().all(axis=1),"adl"]=np.nan
    else:
        a=num(d[cfg["adl"]]); d["adl"]=np.where(a.notna(),(a>=1).astype(float),np.nan)
    d["year"]=num(d[cfg["year"]])
    d["male"]=num(d["ragender"]); ed=num(d["raeducl"])
    d["edu_low"]=(ed==1).astype(float); d.loc[ed.isna(),"edu_low"]=np.nan
    d["edu_mid"]=(ed==2).astype(float); d.loc[ed.isna(),"edu_mid"]=np.nan
    d["rural"]=num(d[cfg["rural"]]) if cfg["rural"] else np.nan
    cc=[x for x in ["hibpe","diabe","hearte","stroke","cancre","lunge"] if x in d.columns]
    d["comorb"]=(d[cc]==1).sum(axis=1).astype(float); d.loc[d[cc].isna().sum(axis=1)>3,"comorb"]=np.nan
    d["smoke"]=num(d["smoken"]) if "smoken" in d.columns else np.nan
    # 【P10 修正】财富五分位按 cohort×wave 分层，与 A2（19）口径一致
    d["wq"]=np.nan
    if cfg["wealth"] in d.columns:
        for _,g in d.groupby("wave"):
            s=num(g[cfg["wealth"]])
            if s.notna().sum()<200: continue
            d.loc[g.index,"wq"]=pd.qcut(s.rank(method="first"),5,labels=False)+1
    # 【P10 修正】严重度 z 在疼痛者内标准化，与 A3（31）一致
    s=pd.Series(d["sev_raw"],dtype=float)
    mu=s[d["pain"]==1].mean(); sd=s[d["pain"]==1].std()
    d["sev_z"]=(s-mu)/sd if sd and sd>0 else np.nan
    d["age10"]=(d["age"]-50)/10
    return d[["ID","wave","year","age","age10","male","edu_low","edu_mid","wq","rural","comorb","smoke","sev_z","pain","dx","adl"]]

rows=[]; desc=[]
for c,cfg in CFG.items():
    d=load(c,cfg).sort_values(["ID","wave"])
    d["nxt"]=d.groupby("ID")["wave"].shift(-1)
    nxt=d.groupby("ID")[["pain","dx","adl","year"]].shift(-1)
    nxt.columns=["pain_n","dx_n","adl_n","year_n"]
    d=pd.concat([d,nxt],axis=1)
    d["gap"]=d["year_n"]-d["year"]
    p=d[(d["pain"]==1)&(d["dx"]==0)&(d["gap"]>0)&(d["gap"]<=3)&d["dx_n"].notna()].copy()
    p["y_dx"]=p["dx_n"]
    if len(p)<150: 
        print(c,"skip, n=",len(p)); continue
    rate=p["y_dx"].mean()
    desc.append(dict(cohort=c,n=len(p),rate=round(rate*100,1),
                     gap_med=p["gap"].median(),adl_n=int(p["adl_n"].notna().sum())))
    terms0=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
    ind=["y_dx"]
    terms=[t for t in terms0+ind if p[t].notna().sum()>100]
    r,n=fit(p,"y_dx",terms)
    if r and "y_dx" in r:
        pass
    # 诊断获得的决定因素（不含 y_dx 本身）
    r1,n1=fit(p,"y_dx",[t for t in terms0 if p[t].notna().sum()>100])
    # 获得诊断 → 后续 ADL（基线 ADL 调整）
    q=p.dropna(subset=["adl_n","adl"]).copy()
    r2,n2=fit(q,"adl_n",["y_dx","adl"]+[t for t in terms0+["sev_z"] if q[t].notna().sum()>100])
    for tag,rr,nn in [("acquire_rate",None,len(p)),("M_acquire",r1,n1),("M_adl",r2,n2)]:
        if tag=="acquire_rate": continue
        if rr is None: continue
        for t,(b,se) in rr.items():
            rows.append(dict(cohort=c,model=tag,term=t,n=nn,OR=round(np.exp(b),3),
                             lo=round(np.exp(b-1.96*se),3),hi=round(np.exp(b+1.96*se),3),b=b,se=se))
    print(f"{c}: base={len(p)} acquire={rate*100:.1f}% adl_model_n={n2}",flush=True)

pd.DataFrame(desc).to_csv(os.path.join(OUT,"24_incident_desc.csv"),index=False,encoding="utf-8-sig")
res=pd.DataFrame(rows); res.to_csv(os.path.join(OUT,"24_incident.csv"),index=False,encoding="utf-8-sig")
pool=[]
for (model,term),g in res.groupby(["model","term"]):
    r=dl(g["b"].values,g["se"].values)
    if r is None: continue
    pool.append(dict(model=model,term=term,k=r["k"],OR=round(r["OR"],3),lo=round(r["lo"],3),
                     hi=round(r["hi"],3),I2=round(r["I2"],1),cohorts=",".join(sorted(g["cohort"]))))
pooled=pd.DataFrame(pool); pooled.to_csv(os.path.join(OUT,"24_incident_pooled.csv"),index=False,encoding="utf-8-sig")
print("\n=== desc ==="); print(pd.DataFrame(desc).to_string(index=False))
print("\n=== pooled ==="); print(pooled.to_string(index=False))
