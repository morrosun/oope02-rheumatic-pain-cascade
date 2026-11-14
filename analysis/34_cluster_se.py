# -*- coding: utf-8 -*-
"""
34: P8 —— 聚类层级敏感性
对 A2 的诊断模型，分别用「个体」与「群（社区/家户）」聚类，
比较标准误比值（设计效应平方根）与 OR 的 95% CI 宽度。
输出：34_cluster_se.csv / 34_cluster_se.md
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

def anypos(df,cols):
    cols=[c for c in cols if c in df.columns]; sub=df[cols]
    out=(sub==1).any(axis=1).astype(float); out[sub.isna().all(axis=1)]=np.nan; return out

def comorb(d,names):
    cols=[c for c in names if c in d.columns]
    cnt=(d[cols]==1).sum(axis=1).astype(float)
    cnt[d[cols].isna().sum(axis=1)>len(cols)//2]=np.nan; return cnt

def add_wq(d,col):
    d["wq"]=np.nan
    if col not in d.columns: return d
    for _,g in d.groupby(["cohort","wave"]):
        s=g[col]
        if s.notna().sum()<200: continue
        d.loc[g.index,"wq"]=pd.qcut(s.rank(method="first"),5,labels=False)+1
    return d

def fit(d,y,terms,cluster):
    dd=d.dropna(subset=[y]+terms+[cluster]).copy()
    if len(dd)<200 or dd[y].nunique()<2: return None,len(dd),0
    X=sm.add_constant(dd[terms],has_constant="add")
    try:
        m=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd[cluster]})
    except Exception:
        return None,len(dd),0
    return {t:(m.params[t],m.bse[t]) for t in terms if t in m.params},len(dd),dd[cluster].nunique()

C6=["hibpe","diabe","hearte","stroke","cancre","lunge"]
CFG={
"CHARLS":dict(path=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",clu="communityID",
  cols=["ID","communityID","wave","age","ragender","raeducl","hatotfa","hrural","arthre",
        "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
  dx="arthre",age="age",edu="raeducl",wealth="hatotfa",rural="hrural",smoke="smoken",waves=[1,2,3,4]),
"ELSA":dict(path=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",clu="hhidc",
  cols=["idauniqc","hhidc","wave","agey","ragender","raeducl","hatotb","arthre","painfr",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="agey",edu="raeducl",wealth="hatotb",rural=None,smoke="smoken",waves=None),
"HRS":dict(path=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",clu="hhhid",
  cols=["hhidpn","hhhid","wave","ragey_b","ragender","raeducl","atotn","rural","arthre","painfr",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="ragey_b",edu="raeducl",wealth="atotn",rural="rural",smoke="smoken",waves=None),
"KLoSA":dict(path=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",clu="hhidc",
  cols=["pid","hhidc","wave","agey","ragender","raeducl","atotf","rural","arthre",
        "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
  dx="arthre",age="agey",edu="raeducl",wealth="atotf",rural="rural",smoke="smoken",waves=None),
"LASI":dict(path=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",clu="hhid",
  cols=["prim_key","hhid","r1agey","ragender","raeducl","hh1atotb","hh1rural","r1arthre","r1jointp",
        "r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge","r1smoken"],
  pain=None,joint="r1jointp",dx="r1arthre",age="r1agey",edu="raeducl",wealth="hh1atotb",rural="hh1rural",smoke="r1smoken",waves=None),
"MHAS":dict(path=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",clu="hhidc",
  cols=["rahhidnp","hhidc","wave","agey","ragender","raeducl","atotb","rural","arthre","painfr",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="agey",edu="raeducl",wealth="atotb",rural="rural",smoke="smoken",waves=None),
"SHARE":dict(path=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",clu="hhidc",
  cols=["mergeid","hhidc","wave","agey","ragender","raeducl","hhatotn","rural","arthre","painfr","pain_s",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="agey",edu="raeducl",wealth="hhatotn",rural="rural",smoke="smoken",waves=None),
}

def build(c,cfg,joint=None,waves=None):
    df,_=pyreadstat.read_dta(cfg["path"],usecols=cfg["cols"])
    if cfg["idc"]!="ID": df=df.rename(columns={cfg["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if waves is not None: d=d[d["wave"].isin(waves)]
    d["ID"]=d["ID"].astype(str); d["CLU"]=d[cfg["clu"]].astype(str)
    d["age"]=num(d[cfg["age"]]); d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    if joint: pv=num(d[joint])
    elif isinstance(cfg["pain"],str): pv=num(d[cfg["pain"]])
    else: pv=anypos(d,cfg["pain"])
    d["pain"]=pv.where(pv.isin([0,1]))
    d["dx"]=num(d[cfg["dx"]]).where(lambda s:s.isin([0,1]))
    d["age10"]=(d["age"]-50)/10; d["male"]=num(d["ragender"])
    ed=num(d[cfg["edu"]]); d["edu_low"]=np.where(ed.notna(),(ed==1).astype(float),np.nan)
    d["edu_mid"]=np.where(ed.notna(),(ed==2).astype(float),np.nan)
    d["rural"]=num(d[cfg["rural"]]) if cfg["rural"] else np.nan
    cc=["r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge"] if c=="LASI" else C6
    d["comorb"]=comorb(d,cc); d["smoke"]=num(d[cfg["smoke"]])
    d["cohort"]=c
    d=add_wq(d,cfg["wealth"])
    # 群固定为个体首次观测的家户
    first=d.sort_values(["ID","wave"]).groupby("ID")["CLU"].first().rename("CLU_FIX")
    d=d.merge(first,on="ID",how="left")
    return d

TERMS=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
res=[]
for c,cfg in CFG.items():
    d=build(c,cfg,joint=cfg.get("joint"),waves=cfg.get("waves"))
    sub=d[d["pain"]==1].copy()
    terms=[t for t in TERMS if sub[t].notna().sum()>100]
    r1,n1,k1=fit(sub,"dx",terms,"ID")
    r2,n2,k2=fit(sub,"dx",terms,"CLU_FIX")
    if r1 is None or r2 is None: 
        print(c,"fail"); continue
    for t in ["male","age10","edu_low"]:
        if t not in r1 or t not in r2: continue
        b1,se1=r1[t]; b2,se2=r2[t]
        res.append(dict(cohort=c,term=t,n=n1,n_clusters=k2,
                        OR_id=np.exp(b1),lo_id=np.exp(b1-1.96*se1),hi_id=np.exp(b1+1.96*se1),se_id=se1,
                        OR_clu=np.exp(b2),lo_clu=np.exp(b2-1.96*se2),hi_clu=np.exp(b2+1.96*se2),se_clu=se2,
                        se_ratio=se2/se1))
out=pd.DataFrame(res)
out.to_csv(os.path.join(OUT,"34_cluster_se.csv"),index=False,encoding="utf-8-sig")
print(out[["cohort","term","n","n_clusters","OR_id","lo_id","hi_id","OR_clu","lo_clu","hi_clu","se_ratio"]].round(3).to_string(index=False))
print("\nSE 比值（群/个体）中位数 = %.3f，最大 = %.3f" % (out.se_ratio.median(), out.se_ratio.max()))
