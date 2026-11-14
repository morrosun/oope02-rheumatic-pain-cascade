# -*- coding: utf-8 -*-
"""
38: 财富五分位的线性假设检验 —— "财富没有梯度"是否可能被连续编码掩盖
对每个库的诊断模型，把 wq 从连续改为 4 个哑变量（Q1 为参照），
做联合 Wald 检验；若 p>0.05 且各档 OR 跨 1，则"无梯度"成立。
输出：38_wq_linearity.csv / .md
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
    o=(sub==1).any(axis=1).astype(float); o[sub.isna().all(axis=1)]=np.nan; return o
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

C6=["hibpe","diabe","hearte","stroke","cancre","lunge"]
CFG={
"CHARLS":dict(p=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",
  cols=["ID","wave","age","ragender","raeducl","hatotfa","hrural","arthre",
        "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
  dx="arthre",age="age",edu="raeducl",wealth="hatotfa",rural="hrural",smoke="smoken",waves=[1,2,3,4]),
"ELSA":dict(p=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),idc="idauniqc",
  cols=["idauniqc","wave","agey","ragender","raeducl","hatotb","arthre","painfr","hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="agey",edu="raeducl",wealth="hatotb",rural=None,smoke="smoken",waves=None),
"HRS":dict(p=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",
  cols=["hhidpn","wave","ragey_b","ragender","raeducl","atotn","rural","arthre","painfr","hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="ragey_b",edu="raeducl",wealth="atotn",rural="rural",smoke="smoken",waves=None),
"KLoSA":dict(p=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),idc="pid",
  cols=["pid","wave","agey","ragender","raeducl","atotf","rural","arthre",
        "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
  dx="arthre",age="agey",edu="raeducl",wealth="atotf",rural="rural",smoke="smoken",waves=None),
"LASI":dict(p=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),idc="prim_key",
  cols=["prim_key","r1agey","ragender","raeducl","hh1atotb","hh1rural","r1arthre","r1jointp",
        "r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge","r1smoken"],
  pain=None,joint="r1jointp",dx="r1arthre",age="r1agey",edu="raeducl",wealth="hh1atotb",rural="hh1rural",smoke="r1smoken",waves=None),
"MHAS":dict(p=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",
  cols=["rahhidnp","wave","agey","ragender","raeducl","atotb","rural","arthre","painfr","hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="agey",edu="raeducl",wealth="atotb",rural="rural",smoke="smoken",waves=None),
"SHARE":dict(p=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),idc="mergeid",
  cols=["mergeid","wave","agey","ragender","raeducl","hhatotn","rural","arthre","painfr","pain_s",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr",dx="arthre",age="agey",edu="raeducl",wealth="hhatotn",rural="rural",smoke="smoken",waves=None),
}
def build(c,cfg,joint=None,waves=None):
    df,_=pyreadstat.read_dta(cfg["p"],usecols=cfg["cols"])
    if cfg["idc"]!="ID": df=df.rename(columns={cfg["idc"]:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if waves is not None: d=d[d["wave"].isin(waves)]
    d["ID"]=d["ID"].astype(str); d["age"]=num(d[cfg["age"]])
    d=d[(d.age>=50)&(d.age<=105)].copy()
    if joint: pv=num(d[joint])
    elif isinstance(cfg["pain"],str): pv=num(d[cfg["pain"]])
    else: pv=anypos(d,cfg["pain"])
    d["pain"]=pv.where(pv.isin([0,1]))
    d["dx"]=num(d[cfg["dx"]]).where(lambda s:s.isin([0,1]))
    d["age10"]=(d.age-50)/10; d["male"]=num(d["ragender"])
    ed=num(d[cfg["edu"]]); d["edu_low"]=np.where(ed.notna(),(ed==1).astype(float),np.nan)
    d["edu_mid"]=np.where(ed.notna(),(ed==2).astype(float),np.nan)
    d["rural"]=num(d[cfg["rural"]]) if cfg["rural"] else np.nan
    cc=["r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge"] if c=="LASI" else C6
    d["comorb"]=comorb(d,cc); d["smoke"]=num(d[cfg["smoke"]]); d["cohort"]=c
    return add_wq(d,cfg["wealth"])

BASE=["age10","male","edu_low","edu_mid","rural","comorb","smoke"]
rows=[]
for c,cfg in CFG.items():
    d=build(c,cfg,joint=cfg.get("joint"),waves=cfg.get("waves"))
    sub=d[(d.pain==1)].copy()
    for q in [2,3,4,5]:
        sub[f"wq{q}"]=(sub["wq"]==q).astype(float)
        sub.loc[sub["wq"].isna(),f"wq{q}"]=np.nan
    terms=[t for t in BASE if sub[t].notna().sum()>100]
    wqterms=[f"wq{q}" for q in [2,3,4,5]]
    # 连续模型
    d1=sub.dropna(subset=["dx","wq"]+terms)
    X1=sm.add_constant(d1[terms+["wq"]],has_constant="add")
    m1=sm.Logit(d1["dx"].astype(float),X1).fit(disp=0,cov_type="cluster",cov_kwds={"groups":d1["ID"]})
    b1,se1=m1.params["wq"],m1.bse["wq"]
    # 分类模型
    d2=sub.dropna(subset=["dx"]+terms+wqterms)
    X2=sm.add_constant(d2[terms+wqterms],has_constant="add")
    m2=sm.Logit(d2["dx"].astype(float),X2).fit(disp=0,cov_type="cluster",cov_kwds={"groups":d2["ID"]})
    R=np.zeros((4,len(m2.params))); cols=list(m2.params.index)
    for i,t in enumerate(wqterms): R[i,cols.index(t)]=1
    test=m2.wald_test(R,scalar=True)
    pval=float(test.pvalue)
    ors={t:round(float(np.exp(m2.params[t])),2) for t in wqterms}
    rows.append(dict(cohort=c,n=int(len(d2)),OR_linear=round(float(np.exp(b1)),3),
                     lo=round(float(np.exp(b1-1.96*se1)),3),hi=round(float(np.exp(b1+1.96*se1)),3),
                     Q2=ors["wq2"],Q3=ors["wq3"],Q4=ors["wq4"],Q5=ors["wq5"],
                     joint_p=round(pval,4)))
out=pd.DataFrame(rows)
out.to_csv(os.path.join(OUT,"38_wq_linearity.csv"),index=False,encoding="utf-8-sig")
print(out.to_string(index=False))
print("\n联合检验 p<0.05 的库数：", int((out.joint_p<0.05).sum()), "/", len(out))
