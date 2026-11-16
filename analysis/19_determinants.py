# -*- coding: utf-8 -*-
"""
19: Determinants of loss along the cascade (diagnosis among people with pain; treatment among those diagnosed and still in pain)
每库 logistic（个体聚类稳健），按面板随机效应合并 log OR。
输出：19_determinants.csv（逐库）、19_pooled.csv（合并）、级联流失决定因素.md
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, warnings
import numpy as np
import pandas as pd
import pyreadstat
import statsmodels.api as sm

warnings.filterwarnings("ignore")

# ---------------- 工具 ----------------
def anypos(df, cols):
    cols=[c for c in cols if c in df.columns]
    sub=df[cols]
    out=(sub==1).any(axis=1).astype(float)
    out[sub.isna().all(axis=1)]=np.nan
    return out

def comorb(d, names):
    cols=[c for c in names if c in d.columns]
    cnt=(d[cols]==1).sum(axis=1).astype(float)
    cnt[d[cols].isna().sum(axis=1) > len(cols)//2] = np.nan
    return cnt

def add_wq(d, col):
    d["wq"]=np.nan
    if col not in d.columns: return d
    for _,g in d.groupby(["cohort","wave"]):
        s=g[col]
        if s.notna().sum()<200: continue
        r=s.rank(method="first")
        d.loc[g.index,"wq"]=pd.qcut(r,5,labels=False)+1
    return d

def fit(d, y, terms):
    dd=d.dropna(subset=[y]+terms).copy()
    if len(dd)<200 or dd[y].nunique()<2: return None, len(dd)
    X=sm.add_constant(dd[terms], has_constant="add")
    try:
        m=sm.Logit(dd[y].astype(float), X).fit(disp=0, cov_type="cluster",
                                               cov_kwds={"groups": dd["ID"]})
    except Exception:
        return None, len(dd)
    out={}
    for t in terms:
        b=m.params.get(t, np.nan); se=m.bse.get(t, np.nan)
        out[t]=(np.exp(b), np.exp(b-1.96*se), np.exp(b+1.96*se), b, se)
    return out, len(dd)

def dl(est, se):
    """DerSimonian-Laird 随机效应合并 log OR"""
    est=np.array(est,float); se=np.array(se,float)
    ok=~np.isnan(est)&~np.isnan(se)&(se>0)
    est,se=est[ok],se[ok]; k=len(est)
    if k==0: return None
    w=1/se**2; fe=np.sum(w*est)/np.sum(w)
    Q=np.sum(w*(est-fe)**2); df=k-1
    C=np.sum(w)-np.sum(w**2)/np.sum(w)
    tau2=max(0.0,(Q-df)/C) if C>0 else 0.0
    wr=1/(se**2+tau2); mu=np.sum(wr*est)/np.sum(wr); sem=np.sqrt(1/np.sum(wr))
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    return dict(k=k, mu=mu, se=sem, OR=np.exp(mu), lo=np.exp(mu-1.96*sem), hi=np.exp(mu+1.96*sem),
                Q=Q, I2=I2, tau2=tau2)

# ---------------- 构建各库分析集 ----------------
def build(cohort, path, cols, pain, dx, tx, age, edu, wealth, rural, comorb_names, smoke, id_col, waves=None, joint=None):
    df,_=pyreadstat.read_dta(path,usecols=cols)
    if id_col!="ID": df=df.rename(columns={id_col:"ID"})
    d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    if waves is not None: d=d[d["wave"].isin(waves)]
    d["ID"]=d["ID"].astype(str)
    d["age"]=pd.to_numeric(d[age],errors="coerce")
    d=d[(d["age"]>=50)&(d["age"]<=105)].copy()
    if joint:
        pv=pd.to_numeric(d[joint],errors="coerce")
    elif isinstance(pain,str):
        pv=pd.to_numeric(d[pain],errors="coerce")
    else:
        pv=anypos(d,pain)
    d["pain"]=pv.where(pv.isin([0,1]))
    d["dx"]=pd.to_numeric(d[dx],errors="coerce").where(lambda s: s.isin([0,1])) if dx in d.columns else np.nan
    d["tx"]=pd.to_numeric(d[tx],errors="coerce").where(lambda s: s.isin([0,1])) if (tx and tx in d.columns) else np.nan
    d["age10"]=(d["age"]-50)/10
    d["male"]=pd.to_numeric(d["ragender"],errors="coerce") if "ragender" in d.columns else np.nan
    ed=pd.to_numeric(d[edu],errors="coerce") if edu in d.columns else pd.Series(np.nan,index=d.index)
    d["edu_low"]=(ed==1).astype(float); d.loc[ed.isna(),"edu_low"]=np.nan
    d["edu_mid"]=(ed==2).astype(float); d.loc[ed.isna(),"edu_mid"]=np.nan
    d["rural"]=pd.to_numeric(d[rural],errors="coerce") if (rural and rural in d.columns) else np.nan
    d["comorb"]=comorb(d, comorb_names)
    d["smoke"]=pd.to_numeric(d[smoke],errors="coerce") if (smoke and smoke in d.columns) else np.nan
    d["cohort"]=cohort
    d=add_wq(d, wealth)
    return d

C6=["hibpe","diabe","hearte","stroke","cancre","lunge"]
CFG={
"CHARLS": dict(path=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),
  cols=["ID","wave","age","ragender","raeducl","hatotfa","hrural","arthre","rxarthr_c",
        "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
  dx="arthre", tx="rxarthr_c", age="age", edu="raeducl", wealth="hatotfa", rural="hrural", smoke="smoken",
  id_col="ID", waves=[1,2,3,4]),
"ELSA": dict(path=os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),
  cols=["idauniqc","wave","agey","ragender","raeducl","hatotb","arthre","painfr",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr", dx="arthre", tx=None, age="agey", edu="raeducl", wealth="hatotb", rural=None, smoke="smoken",
  id_col="idauniqc"),
"HRS": dict(path=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),
  cols=["hhidpn","wave","ragey_b","ragender","raeducl","atotn","rural","arthre","rxarthr","painfr",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr", dx="arthre", tx="rxarthr", age="ragey_b", edu="raeducl", wealth="atotn", rural="rural", smoke="smoken",
  id_col="hhidpn"),
"KLoSA": dict(path=os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),
  cols=["pid","wave","agey","ragender","raeducl","atotf","rural","arthre","rxarthr",
        "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
  dx="arthre", tx="rxarthr", age="agey", edu="raeducl", wealth="atotf", rural="rural", smoke="smoken",
  id_col="pid"),
"LASI": dict(path=os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),
  cols=["prim_key","r1agey","ragender","raeducl","hh1atotb","hh1rural","r1arthre","r1rxarthr","r1jointp",
        "r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge","r1smoken"],
  pain=None, joint="r1jointp", dx="r1arthre", tx="r1rxarthr", age="r1agey", edu="raeducl",
  wealth="hh1atotb", rural="hh1rural", smoke="r1smoken", id_col="prim_key"),
"MHAS": dict(path=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),
  cols=["rahhidnp","wave","agey","ragender","raeducl","atotb","rural","arthre","rxarthr","painfr",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr", dx="arthre", tx="rxarthr", age="agey", edu="raeducl", wealth="atotb", rural="rural", smoke="smoken",
  id_col="rahhidnp"),
"SHARE": dict(path=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),
  cols=["mergeid","wave","agey","ragender","raeducl","hhatotn","rural","arthre","rxpain","painfr","pain_s",
        "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],
  pain="painfr", dx="arthre", tx="rxpain", age="agey", edu="raeducl", wealth="hhatotn", rural="rural", smoke="smoken",
  id_col="mergeid"),
}

D={}
for c,cfg in CFG.items():
    kw={k:v for k,v in cfg.items() if k not in ("path","cols","id_col")}
    if c=="LASI":
        d=build(c,cfg["path"],cfg["cols"],None,cfg["dx"],cfg["tx"],cfg["age"],cfg["edu"],
                cfg["wealth"],cfg["rural"],["r1hibpe","r1diabe","r1hearte","r1stroke","r1cancre","r1lunge"],
                cfg["smoke"],cfg["id_col"],joint=cfg["joint"])
    else:
        d=build(c,cfg["path"],cfg["cols"],cfg["pain"],cfg["dx"],cfg["tx"],cfg["age"],cfg["edu"],
                cfg["wealth"],cfg["rural"],C6,cfg["smoke"],cfg["id_col"],waves=cfg.get("waves"))
    D[c]=d
    print(f"{c}: n={len(d)} pain={int((d['pain']==1).sum())} dx={int((d['dx']==1).sum())} tx_nonmiss={int(d['tx'].notna().sum())}", flush=True)

# SHARE 关节痛面板（w1/2/4）
s_cfg=CFG["SHARE"]
shareB=build("SHARE",s_cfg["path"],s_cfg["cols"],None,s_cfg["dx"],s_cfg["tx"],s_cfg["age"],s_cfg["edu"],
             s_cfg["wealth"],s_cfg["rural"],C6,s_cfg["smoke"],s_cfg["id_col"],waves=[1,2,4],joint="pain_s")

TERMS=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
rows=[]
def run(cohort, d, panel, label):
    for stage, sub in [("D_dx_given_pain", d[d["pain"]==1]),
                       ("T_tx_given_dx_pain", d[(d["pain"]==1)&(d["dx"]==1)])]:
        if stage.startswith("T") and (d["tx"].notna().sum()<200 or cohort=="SHARE"):
            continue
        y = "dx" if stage.startswith("D") else "tx"
        terms=[t for t in TERMS if t in sub.columns and sub[t].notna().sum()>100]
        res, n = fit(sub, y, terms)
        if res is None:
            print(f"  [{cohort} {stage}] 拟合失败 n={n}"); continue
        for t,(orr,lo,hi,b,se) in res.items():
            rows.append(dict(panel=panel, cohort=cohort, stage=stage, term=t, n=n,
                             OR=round(orr,3), lo=round(lo,3), hi=round(hi,3), b=b, se=se))

PANEL_A=["HRS","MHAS","ELSA","SHARE"]
PANEL_B=["CHARLS","KLoSA","LASI","SHARE"]
for c in PANEL_A: run(c, D[c], "A_frequent", c)
for c in PANEL_B:
    run(c, shareB if c=="SHARE" else D[c], "B_joint", c)

det=pd.DataFrame(rows)
det.to_csv(os.path.join(OUT,"19_determinants.csv"),index=False,encoding="utf-8-sig")

# 合并
pool=[]
for (panel,stage,term),g in det.groupby(["panel","stage","term"]):
    r=dl(g["b"].values, g["se"].values)
    if r is None: continue
    pool.append(dict(panel=panel, stage=stage, term=term, k=r["k"], OR=round(r["OR"],3),
                     lo=round(r["lo"],3), hi=round(r["hi"],3), I2=round(r["I2"],1),
                     tau2=round(r["tau2"],3), cohorts=",".join(sorted(g["cohort"]))))
pooled=pd.DataFrame(pool)
pooled.to_csv(os.path.join(OUT,"19_pooled.csv"),index=False,encoding="utf-8-sig")
print("\n=== pooled ===")
print(pooled.to_string(index=False))
