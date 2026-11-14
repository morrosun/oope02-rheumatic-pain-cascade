# -*- coding: utf-8 -*-
"""
22: A3 —— 四状态与 ADL 结局；以及"已诊断 vs 未诊断"（调整疼痛严重度）
输出：22_outcome_gradient.csv / 22_pooled.csv / 结局梯度.md
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
    ok=~np.isnan(est)&~np.isnan(se)&(se>0)
    est,se=est[ok],se[ok]; k=len(est)
    if k==0: return None
    w=1/se**2; fe=np.sum(w*est)/np.sum(w); Q=np.sum(w*(est-fe)**2); df=k-1
    C=np.sum(w)-np.sum(w**2)/np.sum(w); tau2=max(0.0,(Q-df)/C) if C>0 else 0.0
    wr=1/(se**2+tau2); mu=np.sum(wr*est)/np.sum(wr); sem=np.sqrt(1/np.sum(wr))
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    return dict(k=k,OR=np.exp(mu),lo=np.exp(mu-1.96*sem),hi=np.exp(mu+1.96*sem),I2=I2)

def fit(d,y,terms):
    dd=d.dropna(subset=[y]+terms)
    if len(dd)<200 or dd[y].nunique()<2: return None,len(dd)
    X=sm.add_constant(dd[terms],has_constant="add")
    try:
        m=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd["ID"]})
    except Exception:
        return None,len(dd)
    o={}
    for t in terms:
        if t=="const": continue
        b=m.params.get(t,np.nan); se=m.bse.get(t,np.nan)
        o[t]=(b,se)
    return o,len(dd)

COD={"CHARLS":dict(path=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",
      cols=["ID","wave","age","ragender","raeducl","hatotfa","hrural","arthre","rxarthr_c","adlab_c","iadl",
            "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
      pain=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
      dx="arthre",tx="rxarthr_c",age="age",adl="adlab_c",iadl="iadl",sev="nsite",rural="hrural",wealth="hatotfa",waves=[1,2,3,4]),
    "ELSA":dict(path=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",
      cols=["idauniqc","wave","agey","ragender","raeducl","hatotb","arthre","painfr","painlv","adltot6","iadltot2_e",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
      pain="painfr",dx="arthre",tx=None,age="agey",adl="adltot6",iadl="iadltot2_e",sev="painlv",rural=None,wealth="hatotb"),
    "HRS":dict(path=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",
      cols=["hhidpn","wave","ragey_b","ragender","raeducl","atotn","rural","arthre","rxarthr","painfr","painlv","adl6a","iadl5a",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
      pain="painfr",dx="arthre",tx="rxarthr",age="ragey_b",adl="adl6a",iadl="iadl5a",sev="painlv",rural="rural",wealth="atotn"),
    "KLoSA":dict(path=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",
      cols=["pid","wave","agey","ragender","raeducl","atotf","rural","arthre","rxarthr",
            "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13","painhlv",
            "dressb","bathb","eatb","toiltb","bedb_k","urinb","mealsb","shopb","medsb","moneyb","phoneb","transb",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
      pain=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
      dx="arthre",tx="rxarthr",age="agey",adl="ADL6",iadl="IADL6",sev="painhlv",rural="rural",wealth="atotf"),
    "LASI":dict(path=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",
      cols=["prim_key","r1agey","ragender","raeducl","hh1atotb","hh1rural","r1arthre","r1rxarthr","r1jointp","r1painfrq",
            "r1adltot6","r1iadltot_l","r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge","r1smoken"],
      pain=None,joint="r1jointp",dx="r1arthre",tx="r1rxarthr",age="r1agey",adl="r1adltot6",iadl="r1iadltot_l",
      sev="r1painfrq",rural="hh1rural",wealth="hh1atotb"),
    "MHAS":dict(path=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",
      cols=["rahhidnp","wave","agey","ragender","raeducl","atotb","rural","arthre","rxarthr","painfr","painlv","adltot6","iadlfour",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
      pain="painfr",dx="arthre",tx="rxarthr",age="agey",adl="adltot6",iadl="iadlfour",sev="painlv",rural="rural",wealth="atotb"),
    "SHARE":dict(path=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",
      cols=["mergeid","wave","agey","ragender","raeducl","hhatotn","rural","arthre","rxpain","painfr","pain_s","painlv",
            "walkra","dressa","batha","eata","beda","toilta","phonea","medsa","moneya","shopa","mealsa","mapa","housewka","leavhsa","laundrya",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
      pain="painfr",dx="arthre",tx="rxpain",age="agey",adl="ADL6",iadl="IADL9",sev="painlv",rural="rural",wealth="hhatotn"),
}

def load(c,cfg,joint=None,waves=None):
    df,_=pyreadstat.read_dta(cfg["path"],usecols=cfg["cols"])
    if cfg["idc"]!="ID": df=df.rename(columns={cfg["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if waves is not None: d=d[d["wave"].isin(waves)]
    d["ID"]=d["ID"].astype(str)
    d["age"]=num(d[cfg["age"]]); d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    # pain
    if joint:
        pv=num(d[joint]); d["pain"]=pv.where(pv.isin([0,1]))
    elif isinstance(cfg["pain"],str):
        pv=num(d[cfg["pain"]]); d["pain"]=pv.where(pv.isin([0,1]))
    else:
        pv=(d[cfg["pain"]]==1).sum(axis=1).astype(float); pv[d[cfg["pain"]].isna().all(axis=1)]=np.nan
        d["pain"]=np.where(pv.notna(),(pv>0).astype(float),np.nan)
    if (joint or isinstance(cfg["pain"],str)):
        d["sev_raw"]=num(d[cfg["sev"]]) if cfg["sev"] in d.columns else np.nan
    else:
        d["sev_raw"]=pv
    d["dx"]=num(d[cfg["dx"]]).where(lambda s:s.isin([0,1]))
    d["tx"]=num(d[cfg["tx"]]).where(lambda s:s.isin([0,1])) if cfg["tx"] else np.nan
    # ADL / IADL
    if cfg["adl"]=="ADL6":
        cols=["dressb","bathb","eatb","toiltb","bedb_k","urinb"] if c=="KLoSA" else ["walkra","dressa","batha","eata","beda","toilta"]
        d["adl"]=(d[cols]==1).any(axis=1).astype(float); d.loc[d[cols].isna().all(axis=1),"adl"]=np.nan
    else:
        a=num(d[cfg["adl"]]); d["adl"]=np.where(a.notna(),(a>=1).astype(float),np.nan)
    if cfg["iadl"]=="IADL6":
        cols=["mealsb","shopb","medsb","moneyb","phoneb","transb"]
        d["iadl"]=(d[cols]==1).any(axis=1).astype(float); d.loc[d[cols].isna().all(axis=1),"iadl"]=np.nan
    elif cfg["iadl"]=="IADL9":
        cols=["phonea","medsa","moneya","shopa","mealsa","mapa","housewka","leavhsa","laundrya"]
        d["iadl"]=(d[cols]==1).any(axis=1).astype(float); d.loc[d[cols].isna().all(axis=1),"iadl"]=np.nan
    else:
        a=num(d[cfg["iadl"]]); d["iadl"]=np.where(a.notna(),(a>=1).astype(float),np.nan)
    # 协变量
    d["age10"]=(d["age"]-50)/10
    d["male"]=num(d["ragender"])
    ed=num(d["raeducl"]); d["edu_low"]=(ed==1).astype(float); d.loc[ed.isna(),"edu_low"]=np.nan
    d["edu_mid"]=(ed==2).astype(float); d.loc[ed.isna(),"edu_mid"]=np.nan
    d["rural"]=num(d[cfg["rural"]]) if cfg["rural"] else np.nan
    cc=[x for x in ["hibpe","diabe","hearte","stroke","cancre","lunge"] if x in d.columns]
    if c=="LASI": cc=["r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge"]
    d["comorb"]=(d[cc]==1).sum(axis=1).astype(float); d.loc[d[cc].isna().sum(axis=1)>3,"comorb"]=np.nan
    smk="r1smoken" if c=="LASI" else "smoken"
    d["smoke"]=num(d[smk]) if smk in d.columns else np.nan
    d["wq"]=np.nan
    if cfg["wealth"] in d.columns:
        w=num(d[cfg["wealth"]]).rank(method="first"); d["wq"]=pd.qcut(w,5,labels=False)+1
    # 严重度 z（仅在疼痛者内）
    s=pd.Series(d["sev_raw"],dtype=float)
    mu=s[d["pain"]==1].mean(); sd=s[d["pain"]==1].std()
    d["sev_z"]=(s-mu)/sd if sd and sd>0 else np.nan
    d["cohort"]=c
    d["has_tx"]=bool(cfg["tx"]) and (c!="SHARE")
    return d

D={}
for c,cfg in COD.items():
    D[c]=load(c,cfg,joint=cfg.get("joint"),waves=cfg.get("waves"))
    print(c,"n",len(D[c]),"pain",int((D[c]['pain']==1).sum()),"adl",int(D[c]['adl'].notna().sum()),flush=True)
SHB=load("SHARE",COD["SHARE"],joint="pain_s",waves=[1,2,4])

def state(d):
    if bool(d["has_tx"].iloc[0]):
        s=np.select([d["pain"]==0,(d["pain"]==1)&(d["dx"]==0),(d["pain"]==1)&(d["dx"]==1)&(d["tx"]==0),(d["pain"]==1)&(d["dx"]==1)&(d["tx"]==1)],
                    [0,1,2,3],default=np.nan)
    else:
        s=np.select([d["pain"]==0,(d["pain"]==1)&(d["dx"]==0),(d["pain"]==1)&(d["dx"]==1)],[0,1,2],default=np.nan)
    return s

BASE=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
rows=[]
PANELS={"A":["HRS","MHAS","ELSA","SHARE"],"B":["CHARLS","KLoSA","LASI","SHARE"]}
for panel,cohorts in PANELS.items():
    for c in cohorts:
        d = SHB if (c=="SHARE" and panel=="B") else D[c]
        d=d.copy(); d["state"]=state(d)
        terms=[t for t in BASE if d[t].notna().sum()>100]
        # M1：一次拟合，state 哑变量，参照 = 无痛
        dd=d.dropna(subset=["adl","state"]).copy()
        for k in [1,2,3]:
            dd[f"s{k}"]=(dd["state"]==k).astype(float)
        st_terms=[f"s{k}" for k in [1,2,3] if dd[f"s{k}"].sum()>50]
        t2=[t for t in terms if dd[t].notna().sum()>100]
        r,n=fit(dd,"adl",st_terms+t2)
        if r:
            for t in st_terms:
                b,se=r[t]; rows.append(dict(panel=panel,cohort=c,model="M1_vs_nopain",term=t,n=n,
                                            OR=round(np.exp(b),3),lo=round(np.exp(b-1.96*se),3),hi=round(np.exp(b+1.96*se),3),
                                            b=b,se=se))
        # M2：疼痛者中 已诊断 vs 未诊断（调整严重度）
        sub=d[(d["pain"]==1)&d["dx"].notna()].copy(); sub["y"]=sub["dx"]
        t3=[t for t in terms+["sev_z"] if sub[t].notna().sum()>100]
        # M2 正确写法：adl ~ dx + covars (+sev)
        for tag,tt in [("M2_no_sev",[x for x in t3 if x!="sev_z"]),("M2_adj_sev",t3)]:
            r4,n4=fit(sub,"adl",["dx"]+tt)
            if r4 and "dx" in r4:
                b,se=r4["dx"]; rows.append(dict(panel=panel,cohort=c,model=tag,term="dx_vs_undx",n=n4,
                                                OR=round(np.exp(b),3),lo=round(np.exp(b-1.96*se),3),hi=round(np.exp(b+1.96*se),3),b=b,se=se))

res=pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"22_outcome_gradient.csv"),index=False,encoding="utf-8-sig")

pool=[]
for (panel,model,term),g in res.groupby(["panel","model","term"]):
    r=dl(g["b"].values,g["se"].values)
    if r is None: continue
    pool.append(dict(panel=panel,model=model,term=term,k=r["k"],OR=round(r["OR"],3),
                     lo=round(r["lo"],3),hi=round(r["hi"],3),I2=round(r["I2"],1),
                     cohorts=",".join(sorted(g["cohort"]))))
pooled=pd.DataFrame(pool)
pooled.to_csv(os.path.join(OUT,"22_pooled.csv"),index=False,encoding="utf-8-sig")
print("\n=== per-cohort ===")
print(res[["panel","cohort","model","term","n","OR","lo","hi"]].to_string(index=False))
print("\n=== pooled ===")
print(pooled.to_string(index=False))
