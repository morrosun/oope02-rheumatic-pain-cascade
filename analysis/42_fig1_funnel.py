# -*- coding: utf-8 -*-
"""
42: 图 1 数据 —— 同分母级联 + 「每 100 名疼痛者中最终接受治疗的人数」

口径（重要，需进图注）
--------------------
tx 只在「有治疗条目的波次」里可观测（HRS w5–11、MHAS w3–5、CHARLS w1–4、
KLoSA w3–8、LASI w1、SHARE 各面板对应波次；ELSA 无此条目）。
因此对有三段数据的库，整行限制在这些波次内，使三段同分母；
ELSA 只报两段。该限制对诊断覆盖率的影响 ≤1.1 个百分点（见 43_diag_txwaves.txt）。

四个量都在同一次整群 bootstrap 中计算：
  pain_rate        = 疼痛 / 全人群
  dx_given_pain    = 诊断 / 疼痛
  tx_given_dx_pain = 治疗 / (诊断且疼痛，且 tx 可观测)
  yield_per100     = 治疗 / 疼痛            ← 直接估计，不是两比例相乘
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
B   = 1000
RNG = np.random.default_rng(20261114)

STD_W = {1:5370,2:4550,3:3720,4:2970,5:2210,6:1520,7:910,8:630}
EDGES = [50,55,60,65,70,75,80,85,200]
def agegrp(s): return pd.cut(s,bins=EDGES,right=False,labels=list(range(1,9))).astype(float)

W = pd.read_csv(os.path.join(OUT,"16_weights.csv"),dtype={"ID":str})[["cohort","ID","wave","weight"]]

SPEC = {
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

def load(c, spec, tier):
    cols = list(dict.fromkeys(spec["base"] + spec["sites"]))
    df,_ = pyreadstat.read_dta(spec["p"], usecols=cols)
    if spec["idc"] != "ID": df = df.rename(columns={spec["idc"]:"ID"})
    d = df.copy()
    if "wave" not in d.columns: d["wave"] = 1
    if spec["waves"] is not None: d = d[d["wave"].isin(spec["waves"])]
    if c == "SHARE":
        d = d[d["wave"].isin([1,2,4] if tier=="B" else [5,6,7,8])]
    d["ID"] = d["ID"].astype(str); d["CLU"] = d[spec["clu"]].astype(str)
    d["age"] = pd.to_numeric(d[spec["age"]], errors="coerce")
    d = d[(d["age"]>=50)&(d["age"]<=105)].copy()
    d["ageg"] = agegrp(d["age"])
    if tier=="B" and c=="LASI":
        pv = pd.to_numeric(d["r1jointp"], errors="coerce")
    elif spec["pain_mode"]=="list":
        sub = d[spec["sites"]]
        pv = (sub==1).any(axis=1).astype(float); pv[sub.isna().all(axis=1)] = np.nan
    elif tier=="B":
        pv = pd.to_numeric(d["pain_s"], errors="coerce")
    else:
        pv = pd.to_numeric(d[spec["pain"]], errors="coerce")
    d["pain"] = np.where(pv.notna(), np.where(pv>0,1.0,0.0), np.nan)
    d["dx"] = pd.to_numeric(d[spec["dx"]], errors="coerce").where(lambda s:s.isin([0,1]))
    d["tx"] = (pd.to_numeric(d[spec["tx"]], errors="coerce").where(lambda s:s.isin([0,1]))
               if spec["tx"] else np.nan)
    d["wave"] = d["wave"].astype(float)
    d = d.merge(W[W.cohort==c][["ID","wave","weight"]], on=["ID","wave"], how="left")
    first = d.sort_values(["ID","wave"]).groupby("ID")["CLU"].first().rename("CLU_FIX")
    d = d.merge(first, on="ID", how="left")
    return d[["ID","CLU_FIX","wave","ageg","pain","dx","tx","weight"]]

def boot_multi(d, restrict_tx=True):
    d = d.copy()
    waves_tx = []
    if restrict_tx and d["tx"].notna().any():
        byw = d.groupby("wave")["tx"].apply(lambda s: s.notna().sum())
        waves_tx = sorted(byw[byw>0].index.tolist())
        d = d[d["wave"].isin(waves_tx)]
    m = d.dropna(subset=["pain","dx","ageg","weight","CLU_FIX"]).copy()
    if len(m)==0: return {}
    pain = (m["pain"]==1)
    dxp  = pain & (m["dx"]==1)
    txo  = dxp & m["tx"].notna()
    m["P"] = pain.astype(float)
    m["D"] = dxp.astype(float)
    m["T"] = (dxp & (m["tx"]==1)).astype(float)
    m["O"] = txo.astype(float)          # tx 可观测的诊断者
    for c in "PDTO": m["_"+c] = m[c]*m["weight"]
    g = m.groupby(["CLU_FIX","ageg"]).agg(
        w=("weight","sum"), P=("_P","sum"), D=("_D","sum"),
        T=("_T","sum"), O=("_O","sum")).reset_index()
    clus = g["CLU_FIX"].unique(); G = sorted(m["ageg"].dropna().unique())
    ci = {c:i for i,c in enumerate(clus)}; gi = {a:i for i,a in enumerate(G)}
    M = {k:np.zeros((len(clus),len(G))) for k in "WPDTO"}
    for k,col in zip("WPDTO",["w","P","D","T","O"]):
        M[k][g["CLU_FIX"].map(ci).values, g["ageg"].map(gi).values] = g[col].values
    wg = np.array([STD_W[a] for a in G],float)
    def rates(Mx):
        w, P, D, T, O = (Mx[k].sum(0) for k in "WPDTO")
        ok = w>0
        if not ok.any(): return np.full(4,np.nan)
        r_pain = P/w
        r_dx   = np.divide(D,P,out=np.full_like(D,np.nan),where=P>0)
        r_tx   = np.divide(T,O,out=np.full_like(T,np.nan),where=O>0)
        r_y    = np.divide(T,P,out=np.full_like(T,np.nan),where=P>0)
        def sm(r):
            o = ok & ~np.isnan(r)
            return np.sum(wg[o]*r[o])/np.sum(wg[o]) if o.any() else np.nan
        return np.array([sm(r_pain), sm(r_dx), sm(r_tx), sm(r_y)])
    point = rates(M)
    C = len(clus); reps = np.empty((B,4))
    for b in range(B):
        idx = RNG.integers(0,C,C)
        reps[b] = rates({k:v[idx,:] for k,v in M.items()})
    lo, hi = np.nanpercentile(reps,[2.5,97.5],axis=0)
    n_dx = int(m["D"].sum()); n_obs = int(m["O"].sum())
    return dict(point=point*100, lo=lo*100, hi=hi*100, n=len(m), G=C,
                n_pain=int(m["P"].sum()), n_dx=n_dx, n_obs=n_obs,
                tx_miss=1-n_obs/n_dx if n_dx else np.nan, waves=",".join(str(int(x)) for x in waves_tx))

rows=[]
for tier, cohorts in [("A",["HRS","MHAS","ELSA","SHARE"]), ("B",["CHARLS","KLoSA","LASI","SHARE"])]:
    for c in cohorts:
        d = load(c, SPEC[c], tier)
        r = boot_multi(d)
        if not r: continue
        for i,metric in enumerate(["pain_rate","dx_given_pain","tx_given_dx_pain","tx_per100_pain"]):
            rows.append(dict(tier=tier, cohort=c, metric=metric,
                             est=r["point"][i], lo=r["lo"][i], hi=r["hi"][i],
                             n=r["n"], n_pain=r["n_pain"], n_dx=r["n_dx"], n_tx_obs=r["n_obs"],
                             tx_miss_dx=r["tx_miss"], n_clusters=r["G"], waves_tx=r["waves"]))
        print(f"{tier} {c:7s} [{r['waves']:14s}] pain={r['point'][0]:5.1f} dx={r['point'][1]:5.1f} "
              f"tx={r['point'][2]:5.1f} yield={r['point'][3]:5.1f}  txmiss={r['tx_miss']*100:4.1f}%", flush=True)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT,"42_fig1_funnel.csv"), index=False, encoding="utf-8-sig")
print("\nWROTE 42_fig1_funnel.csv  rows =", len(res))
