# -*- coding: utf-8 -*-
"""30: 复盘审计之二 —— A3 参照组重定义、KLoSA 非单调、tx 缺失"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat, statsmodels.api as sm
import warnings; warnings.filterwarnings("ignore")
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
L=[]; w=lambda s:(L.append(str(s)),print(s,flush=True))
def num(s): return pd.to_numeric(s,errors="coerce")

# ---------- A. 各库 tx 在诊断者中的缺失率 ----------
w("="*72); w("【A】诊断者中治疗变量的缺失率（决定治疗层分母是否可信）"); w("="*72)
SPEC=[("CHARLS",os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),"ID",
       "arthre","rxarthr_c",["ID","wave","age","arthre","rxarthr_c"],None),
      ("HRS",os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),"hhidpn",
       "arthre","rxarthr",["hhidpn","wave","ragey_b","arthre","rxarthr"],None),
      ("KLoSA",os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),"pid",
       "arthre","rxarthr",["pid","wave","agey","arthre","rxarthr"],None),
      ("LASI",os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),"prim_key",
       "r1arthre","r1rxarthr",["prim_key","r1agey","r1arthre","r1rxarthr"],None),
      ("MHAS",os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),"rahhidnp",
       "arthre","rxarthr",["rahhidnp","wave","agey","arthre","rxarthr"],None),
      ("SHARE",os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),"mergeid",
       "arthre","rxpain",["mergeid","wave","agey","arthre","rxpain"],None)]
for c,p,idc,dxc,txc,cols,_ in SPEC:
    df,_=pyreadstat.read_dta(p,usecols=cols)
    if idc!="ID": df=df.rename(columns={idc:"ID"})
    age="r1agey" if c=="LASI" else ("ragey_b" if c=="HRS" else "agey")
    if c=="CHARLS": age="age"
    d=df[(df[age]>=50)&(df[age]<=105)]
    dx=num(d[dxc]).where(lambda s:s.isin([0,1])); tx=num(d[txc]).where(lambda s:s.isin([0,1]))
    n_dx=int((dx==1).sum()); n_tx=int(((dx==1)&tx.notna()).sum())
    byw=pd.DataFrame({"dx":dx,"tx":tx,"wave":d["wave"] if "wave" in d.columns else 1}).groupby("wave").apply(
        lambda g: pd.Series({"dx_n":int((g.dx==1).sum()),"tx_n":int(((g.dx==1)&g.tx.notna()).sum())}))
    byw["miss%"]=((byw.dx_n-byw.tx_n)/byw.dx_n.replace(0,np.nan)*100).round(1)
    w(f"{c:7s} 诊断者={n_dx:6d}  有tx={n_tx:6d}  缺失={100-n_tx/max(1,n_dx)*100:5.1f}%   逐波缺失%: {byw['miss%'].to_dict()}")

# ---------- B. A3 参照组重定义 ----------
w("\n"+"="*72); w("【B】A3 M1：把「无痛已诊断」从参照组剔出后的 OR 变化"); w("="*72)
CFG={"HRS":dict(p=os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),idc="hhidpn",
      cols=["hhidpn","wave","ragey_b","ragender","raeducl","atotn","rural","arthre","rxarthr","painfr","adl6a",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],pain="painfr",dx="arthre",tx="rxarthr",
      age="ragey_b",adl="adl6a",rural="rural",wealth="atotn"),
     "MHAS":dict(p=os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),idc="rahhidnp",
      cols=["rahhidnp","wave","agey","ragender","raeducl","atotb","rural","arthre","rxarthr","painfr","adltot6",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],pain="painfr",dx="arthre",tx="rxarthr",
      age="agey",adl="adltot6",rural="rural",wealth="atotb"),
     "CHARLS":dict(p=os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),idc="ID",
      cols=["ID","wave","age","ragender","raeducl","hatotfa","hrural","arthre","rxarthr_c","adlab_c",
            "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15",
            "hibpe","diabe","hearte","stroke","cancre","lunge","smoken"],pain="list",dx="arthre",tx="rxarthr_c",
      age="age",adl="adlab_c",rural="hrural",wealth="hatotfa")}
def fit(d,y,terms):
    dd=d.dropna(subset=[y]+terms)
    if len(dd)<200 or dd[y].nunique()<2: return None,len(dd)
    X=sm.add_constant(dd[terms],has_constant="add")
    try: m=sm.Logit(dd[y].astype(float),X).fit(disp=0,cov_type="cluster",cov_kwds={"groups":dd["ID"]})
    except Exception: return None,len(dd)
    return {t:(m.params[t],m.bse[t]) for t in terms if t in m.params},len(dd)
BASE=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]
for c,cfg in CFG.items():
    df,_=pyreadstat.read_dta(cfg["p"],usecols=cfg["cols"])
    d=df.rename(columns={cfg["idc"]:"ID"}) if cfg["idc"]!="ID" else df.copy()
    if "wave" not in d.columns: d["wave"]=1
    d["ID"]=d["ID"].astype(str); d["age"]=num(d[cfg["age"]]); d=d[(d.age>=50)&(d.age<=105)].copy()
    if cfg["pain"]=="list":
        pc=[x for x in cfg["cols"] if x.startswith("da042s")]
        d["pain"]=np.where((d[pc]==1).any(axis=1),1.0,np.nan); d.loc[d[pc].notna().any(axis=1),"pain"]=np.where((d.loc[d[pc].notna().any(axis=1),pc]==1).any(axis=1),1.0,0.0)
    else:
        d["pain"]=num(d[cfg["pain"]]).where(lambda s:s.isin([0,1]))
    d["dx"]=num(d[cfg["dx"]]).where(lambda s:s.isin([0,1]))
    d["tx"]=num(d[cfg["tx"]]).where(lambda s:s.isin([0,1]))
    a=num(d[cfg["adl"]]); d["adl"]=np.where(a.notna(),(a>=1).astype(float),np.nan)
    d["age10"]=(d.age-50)/10; d["male"]=num(d["ragender"])
    ed=num(d["raeducl"]); d["edu_low"]=np.where(ed.notna(),(ed==1).astype(float),np.nan); d["edu_mid"]=np.where(ed.notna(),(ed==2).astype(float),np.nan)
    d["rural"]=num(d[cfg["rural"]]); cc=[x for x in ["hibpe","diabe","hearte","stroke","cancre","lunge"] if x in d.columns]
    d["comorb"]=(d[cc]==1).sum(axis=1).astype(float); d.loc[d[cc].isna().sum(axis=1)>3,"comorb"]=np.nan
    d["smoke"]=num(d["smoken"]); d["wq"]=pd.qcut(num(d[cfg["wealth"]]).rank(method="first"),5,labels=False)+1
    d["has_tx"]=True
    terms=[t for t in BASE if d[t].notna().sum()>100]
    # 现方案：参照 = pain==0（含已诊断）
    d["state"]=np.select([d.pain==0,(d.pain==1)&(d.dx==0),(d.pain==1)&(d.dx==1)&(d.tx==0),(d.pain==1)&(d.dx==1)&(d.tx==1)],[0,1,2,3],default=np.nan)
    dd=d.dropna(subset=["adl","state"]).copy()
    for k in [1,2,3]: dd[f"s{k}"]=(dd.state==k).astype(float)
    st=[f"s{k}" for k in [1,2,3] if dd[f"s{k}"].sum()>50]
    r,n=fit(dd,"adl",st+terms)
    # 新方案：参照 = pain==0 & dx==0；无痛已诊断单独为 s4
    d["state2"]=np.select([(d.pain==0)&(d.dx==0),(d.pain==1)&(d.dx==0),(d.pain==1)&(d.dx==1)&(d.tx==0),(d.pain==1)&(d.dx==1)&(d.tx==1),(d.pain==0)&(d.dx==1)],[0,1,2,3,4],default=np.nan)
    d2=d.dropna(subset=["adl","state2"]).copy()
    for k in [1,2,3,4]: d2[f"t{k}"]=(d2.state2==k).astype(float)
    st2=[f"t{k}" for k in [1,2,3,4] if d2[f"t{k}"].sum()>50]
    r2,n2=fit(d2,"adl",st2+terms)
    w(f"\n{c}  现参照组 n={n}（其中无痛已诊断 {int(((d.pain==0)&(d.dx==1)).sum())}）  新参照组 n={n2}")
    if r:
        w("  现方案 OR: " + "  ".join(f"{k}={np.exp(v[0]):.2f}" for k,v in r.items() if k.startswith("s")))
    if r2:
        w("  新方案 OR: " + "  ".join(f"{k}={np.exp(v[0]):.2f}" for k,v in r2.items() if k.startswith("t")))

# ---------- C. KLoSA 非单调核对 ----------
w("\n"+"="*72); w("【C】KLoSA 四状态 ADL 原始率（核对 s2<1 是否真实）"); w("="*72)
kl,_=pyreadstat.read_dta(os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),
    usecols=["pid","wave","agey","arthre","rxarthr","pain_2","pain_3","pain_4","pain_5","pain_8","pain_9",
             "pain_10","pain_11","pain_12","pain_13","dressb","bathb","eatb","toiltb","bedb_k","urinb"])
pc=["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"]
kl["pain"]=np.where((kl[pc]==1).any(axis=1),1.0,0.0); kl.loc[kl[pc].isna().all(axis=1),"pain"]=np.nan
kl["dx"]=num(kl.arthre).where(lambda s:s.isin([0,1])); kl["tx"]=num(kl.rxarthr).where(lambda s:s.isin([0,1]))
ac=["dressb","bathb","eatb","toiltb","bedb_k","urinb"]
kl["adl"]=np.where((kl[ac]==1).any(axis=1),1.0,0.0); kl.loc[kl[ac].isna().all(axis=1),"adl"]=np.nan
kl["st"]=np.select([kl.pain==0,(kl.pain==1)&(kl.dx==0),(kl.pain==1)&(kl.dx==1)&(kl.tx==0),(kl.pain==1)&(kl.dx==1)&(kl.tx==1)],[0,1,2,3],default=np.nan)
g=kl.dropna(subset=["adl","st"]).groupby("st").agg(n=("adl","size"),adl=("adl","mean"),age=("agey","median"))
g["adl"]=(g["adl"]*100).round(1); w(g.to_string())
w("  → 若 s2 的 ADL 率低于 s1，说明「诊断未治疗」组确实比「未诊断」组更健康，非单调是真实的。")

open(os.path.join(OUT,"30_audit2.txt"),"w",encoding="utf-8").write("\n".join(L))
print("\nWROTE 30_audit2.txt")
