# -*- coding: utf-8 -*-
"""21: 敏感性 —— 把疼痛严重度加进诊断模型，看 edu_low 是否衰减"""
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

def num(s): return pd.to_numeric(s, errors="coerce")

def fit_or(d, y, terms):
    dd=d.dropna(subset=[y]+terms)
    if len(dd)<200 or dd[y].nunique()<2: return np.nan, len(dd)
    X=sm.add_constant(dd[terms],has_constant="add")
    try:
        m=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd["ID"]})
        return float(np.exp(m.params["edu_low"])), len(dd)
    except Exception:
        return np.nan, len(dd)

CFG=[
 # cohort, panel, path, idcol, age, dx, edu, wealth, rural, comorb, smoke, painbuilder, sev
 ("HRS","A",os.path.join(DL, "3. HRS  美国/HRS_美国/Working_data/hrs.dta"),"hhidpn","ragey_b","arthre","raeducl","atotn","rural",
  ["hibpe","diabe","hearte","stroke","cancre","lunge"],"smoken","painfr","painlv"),
 ("MHAS","A",os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Working_data/mhas.dta"),"rahhidnp","agey","arthre","raeducl","atotb","rural",
  ["hibpe","diabe","hearte","stroke","cancre","lunge"],"smoken","painfr","painlv"),
 ("ELSA","A",os.path.join(DL, "2. ELSA 英国/Working_data/elsa.dta"),"idauniqc","agey","arthre","raeducl","hatotb",None,
  ["hibpe","diabe","hearte","stroke","cancre","lunge"],"smoken","painfr","painlv"),
 ("SHARE","A",os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Working_data/share.dta"),"mergeid","agey","arthre","raeducl","hhatotn","rural",
  ["hibpe","diabe","hearte","stroke","cancre","lunge"],"smoken","painfr","painlv"),
 ("CHARLS","B",os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta"),"ID","age","arthre","raeducl","hatotfa","hrural",
  ["hibpe","diabe","hearte","stroke","cancre","lunge"],"smoken",
  ["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],"nsite"),
 ("KLoSA","B",os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Working_data/klosa.dta"),"pid","agey","arthre","raeducl","atotf","rural",
  ["hibpe","diabe","hearte","stroke","cancre","lunge"],"smoken",
  ["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],"painhlv"),
 ("LASI","B",os.path.join(DL, "5. LASI 印度/LASI_印度/Working_data/lasi.dta"),"prim_key","r1agey","r1arthre","raeducl","hh1atotb","hh1rural",
  ["r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge"],"r1smoken","r1jointp","r1painfrq"),
]

rows=[]
for (c,panel,path,idc,age,dx,edu,wth,rur,cm,smk,pain,sev) in CFG:
    cols=list(dict.fromkeys([idc,"ragender",age,dx,edu,wth,rur,smk]+[x for x in cm if x]+
                            ([pain] if isinstance(pain,str) else pain)+([sev] if sev!="nsite" else [])))
    cols=[x for x in cols if x]
    df,_=pyreadstat.read_dta(path,usecols=cols)
    if idc!="ID": df=df.rename(columns={idc:"ID"})
    d=df.copy(); d["ID"]=d["ID"].astype(str)
    d["age"]=num(d[age]); d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    if isinstance(pain,str):
        pv=num(d[pain])
        d["pain"]=pv.where(pv.isin([0,1]))
    else:
        pv=(d[pain]==1).sum(axis=1).astype(float)
        pv[d[pain].isna().all(axis=1)]=np.nan
        d["pain"]=np.where(pv.notna(),(pv>0).astype(float),np.nan)
    d["dx"]=num(d[dx]).where(lambda s: s.isin([0,1]))
    d["age10"]=(d["age"]-50)/10
    d["male"]=num(d["ragender"]) if "ragender" in d.columns else np.nan
    ed=num(d[edu]); d["edu_low"]=(ed==1).astype(float); d.loc[ed.isna(),"edu_low"]=np.nan
    d["edu_mid"]=(ed==2).astype(float); d.loc[ed.isna(),"edu_mid"]=np.nan
    d["rural"]=num(d[rur]) if rur else np.nan
    cc=[x for x in cm if x in d.columns]
    d["comorb"]=(d[cc]==1).sum(axis=1).astype(float); d.loc[d[cc].isna().sum(axis=1)>3,"comorb"]=np.nan
    d["smoke"]=num(d[smk]) if smk in d.columns else np.nan
    # 财富五分位
    d["wq"]=np.nan
    if wth in d.columns:
        s=num(d[wth]).rank(method="first")
        d["wq"]=pd.qcut(s,5,labels=False)+1
    # 严重度 z
    if sev=="nsite":
        sv=pv
    else:
        sv=num(d[sev]) if sev in d.columns else np.nan
    sd=np.nanstd(sv); mu=np.nanmean(sv)
    d["sev_z"]=(sv-mu)/sd if sd and sd>0 else np.nan
    sub=d[d["pain"]==1]
    base=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
    base=[t for t in base if sub[t].notna().sum()>100]
    o1,n1=fit_or(sub,"dx",base)
    o2,n2=fit_or(sub,"dx",base+["sev_z"])
    rows.append(dict(panel=panel,cohort=c,n=n1,OR_edu_base=round(o1,3) if pd.notna(o1) else np.nan,
                     OR_edu_sev=round(o2,3) if pd.notna(o2) else np.nan,
                     sev_var=sev, sev_sd=round(sd,2) if sd==sd else np.nan))
    print(rows[-1], flush=True)

res=pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"21_severity_sensitivity.csv"),index=False,encoding="utf-8-sig")
print(res.to_string(index=False))
