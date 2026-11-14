# -*- coding: utf-8 -*-
"""
17: 七库加权 + 年龄标准化的风湿疼痛诊疗级联（第 1 步）
输出：
  17_cascade_tierA.csv / .txt   任何频繁疼痛层
  17_cascade_tierB.csv          关节特异性疼痛层
  17_four_state.csv             四状态 × ADL
  17_missing.csv                三变量联合缺失矩阵
  第七库级联汇总表.md            可读汇总
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os
import numpy as np
import pandas as pd
import pyreadstat

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))

# WHO World Standard Population 2000-2025，50+ 重归一
STD_W = {1:5370, 2:4550, 3:3720, 4:2970, 5:2210, 6:1520, 7:910, 8:630}
AGE_EDGES = [50,55,60,65,70,75,80,85,200]

def agegrp(s):
    return pd.cut(s, bins=AGE_EDGES, right=False, labels=list(range(1,9))).astype(float)

def anypos(df, cols):
    cols=[c for c in cols if c in df.columns]
    if not cols: return pd.Series(np.nan, index=df.index)
    sub=df[cols]
    out=(sub==1).any(axis=1).astype(float)
    out[sub.isna().all(axis=1)]=np.nan
    return out

def wrate(df, y, weight):
    m = df[[y,weight]].dropna()
    if len(m)==0 or m[weight].sum()==0: return np.nan
    return np.average(m[y], weights=m[weight])

def std_rate(df, y, weight, ageg):
    """年龄标准化：组内加权率，再按 WHO 标准人口加权"""
    d = df[[y,weight,ageg]].dropna()
    if len(d)==0: return np.nan, 0
    num=0.0; den=0.0; used=0
    for g,wg in STD_W.items():
        sub=d[d[ageg]==g]
        if len(sub)==0 or sub[weight].sum()==0: continue
        r=np.average(sub[y], weights=sub[weight])
        num += wg*r; den += wg; used+=1
    return (num/den if den>0 else np.nan), used

# ---------------- 载入权重 ----------------
W = pd.read_csv(os.path.join(OUT,"16_weights.csv"), dtype={"ID":str})
W = W[["cohort","ID","wave","weight"]]

rowsA=[]; rowsB=[]; rows4=[]; rowsM=[]

def build(cohort, path, cols, pain_map, dx_col, tx_col, adl_cols, age_col, wave_filter=None, id_col="ID"):
    df,_ = pyreadstat.read_dta(path, usecols=cols)
    if id_col != "ID": df = df.rename(columns={id_col:"ID"})
    df["ID"]=df["ID"].astype(str)
    d = df.copy()
    if "wave" not in d.columns: d["wave"]=1
    d["age"]=d[age_col]
    d = d[(d["age"]>=50) & (d["age"]<=105)]
    if wave_filter is not None:
        d = d[d["wave"].isin(wave_filter)]
    d["ageg"]=agegrp(d["age"])
    if callable(pain_map):
        d["pain"]=pain_map(d)
    else:
        d["pain"]=d[pain_map]
    d["pain"]=d["pain"].where(d["pain"].isin([0,1]))
    d["dx"]=d[dx_col].where(d[dx_col].isin([0,1])) if dx_col in d.columns else np.nan
    d["tx"]=d[tx_col].where(d[tx_col].isin([0,1])) if (tx_col and tx_col in d.columns) else np.nan
    if adl_cols:
        avail=[c for c in adl_cols if c in d.columns]
        d["adl"]=(d[avail]>=1).any(axis=1).astype(float)
        d.loc[d[avail].isna().all(axis=1),"adl"]=np.nan
    else:
        d["adl"]=np.nan
    d["cohort"]=cohort
    # 权重
    d = d.merge(W[W["cohort"]==cohort][["ID","wave","weight"]], on=["ID","wave"], how="left")
    return d[["cohort","ID","wave","age","ageg","pain","dx","tx","adl","weight"]]

# ---------------- 逐库构建 ----------------
frames={}

# CHARLS：w1–w4（权重与 rx 均可用；w5 无权重、无 rx，单独说明）
ch_cols=["ID","wave","age","arthre","rxarthr","rxarthr_c","adlab_c",
         "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10",
         "da042s11","da042s12","da042s13","da042s14","da042s15"]
frames["CHARLS"]=build("CHARLS", os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),
    ch_cols, lambda d: anypos(d,["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9",
                                 "da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"]),
    "arthre","rxarthr_c",["adlab_c"],"age", wave_filter=[1,2,3,4])

# ELSA
frames["ELSA"]=build("ELSA", os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),
    ["idauniqc","wave","agey","arthre","painfr","adltot6"], "painfr","arthre",None,["adltot6"],"agey", id_col="idauniqc")

# HRS
frames["HRS"]=build("HRS", os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),
    ["hhidpn","wave","ragey_b","arthre","rxarthr","painfr","adl6a"], "painfr","arthre","rxarthr",["adl6a"],"ragey_b", id_col="hhidpn")

# KLoSA
frames["KLoSA"]=build("KLoSA", os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),
    ["pid","wave","agey","arthre","rxarthr","pain_2","pain_3","pain_4","pain_5","pain_8","pain_9",
     "pain_10","pain_11","pain_12","pain_13","dressb","bathb","eatb","toiltb","bedb_k","urinb"],
    lambda d: anypos(d,["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"]),
    "arthre","rxarthr",["dressb","bathb","eatb","toiltb","bedb_k","urinb"],"agey", id_col="pid")

# LASI
frames["LASI"]=build("LASI", os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),
    ["prim_key","r1agey","r1arthre","r1rxarthr","r1jointp","r1adltot6"], "r1jointp","r1arthre","r1rxarthr",["r1adltot6"],"r1agey", id_col="prim_key")
frames["LASI"]["wave"]=1

# MHAS
frames["MHAS"]=build("MHAS", os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),
    ["rahhidnp","wave","agey","arthre","rxarthr","painfr","adltot6"], "painfr","arthre","rxarthr",["adltot6"],"agey", id_col="rahhidnp")

# SHARE（任何疼痛层只用 painfr，w5–8）
frames["SHARE"]=build("SHARE", os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),
    ["mergeid","wave","agey","arthre","rxpain","painfr","walkra","dressa","batha","eata","beda","toilta"],
    "painfr","arthre","rxpain",["walkra","dressa","batha","eata","beda","toilta"],"agey", wave_filter=[5,6,7,8], id_col="mergeid")

# ---------------- 表 A：任何频繁疼痛 ----------------
def summarize(cohort, d, tier):
    n_all=len(d); n_w=d["weight"].notna().sum()
    d = d[d["weight"].notna()].copy()
    d["pain_w"]=d["pain"]
    p_cr = wrate(d,"pain","weight")
    p_st,_ = std_rate(d,"pain","weight","ageg")
    # 诊断覆盖：分母 = 有疼痛且 dx 有值
    sp = d[(d["pain"]==1)&d["dx"].notna()].copy()
    dx_cr = wrate(sp,"dx","weight") if len(sp) else np.nan
    dx_st,_ = std_rate(sp,"dx","weight","ageg") if len(sp) else (np.nan,0)
    # 治疗覆盖：分母 = 诊断且 tx 有值
    if d["tx"].notna().any():
        st = d[(d["dx"]==1)&d["tx"].notna()].copy()
        tx_cr = wrate(st,"tx","weight") if len(st) else np.nan
        tx_st,_ = std_rate(st,"tx","weight","ageg") if len(st) else (np.nan,0)
        n_tx_den=len(st)
    else:
        tx_cr=tx_st=np.nan; n_tx_den=0
    if d["tx"].notna().any():
        stp = d[(d["pain"]==1)&(d["dx"]==1)&d["tx"].notna()].copy()
        txp_cr = wrate(stp,"tx","weight") if len(stp) else np.nan
        txp_st,_ = std_rate(stp,"tx","weight","ageg") if len(stp) else (np.nan,0)
        n_txp=len(stp)
    else:
        txp_cr=txp_st=np.nan; n_txp=0
    return dict(cohort=cohort, tier=tier, n_all=n_all, n_weighted=int(n_w),
                n_pain=int((d["pain"]==1).sum()),
                pain_crude=round(p_cr*100,1) if pd.notna(p_cr) else np.nan,
                pain_std=round(p_st*100,1) if pd.notna(p_st) else np.nan,
                n_pain_dx=int(len(sp)),
                dx_crude=round(dx_cr*100,1) if pd.notna(dx_cr) else np.nan,
                dx_std=round(dx_st*100,1) if pd.notna(dx_st) else np.nan,
                n_dx_tx=int(n_tx_den),
                tx_crude=round(tx_cr*100,1) if pd.notna(tx_cr) else np.nan,
                tx_std=round(tx_st*100,1) if pd.notna(tx_st) else np.nan,
                n_tx_pain=int(n_txp),
                txpain_crude=round(txp_cr*100,1) if pd.notna(txp_cr) else np.nan,
                txpain_std=round(txp_st*100,1) if pd.notna(txp_st) else np.nan)

PANEL_A = ["HRS","MHAS","ELSA","SHARE"]
PANEL_B = ["CHARLS","KLoSA","LASI","SHARE"]
for c in PANEL_A:
    rowsA.append(summarize(c,frames[c],"A_frequent_pain"))

# ---------------- 表 B：关节特异性疼痛（CHARLS/KLoSA/LASI/SHARE w1,2,4） ----------------
def build_joint(cohort, path, cols, pain_cols, dx_col, tx_col, age_col, wave_filter=None, pain_single=None, id_col="ID"):
    df,_=pyreadstat.read_dta(path,usecols=cols)
    if id_col != "ID": df = df.rename(columns={id_col:"ID"})
    df["ID"]=df["ID"].astype(str); d=df.copy()
    if "wave" not in d.columns: d["wave"]=1
    d["age"]=d[age_col]; d=d[(d["age"]>=50)&(d["age"]<=105)]
    if wave_filter is not None: d=d[d["wave"].isin(wave_filter)]
    d["ageg"]=agegrp(d["age"])
    d["pain"]= d[pain_single] if pain_single else anypos(d,pain_cols)
    d["pain"]=d["pain"].where(d["pain"].isin([0,1]))
    d["dx"]=d[dx_col].where(d[dx_col].isin([0,1]))
    d["tx"]=d[tx_col].where(d[tx_col].isin([0,1])) if tx_col else np.nan
    d["cohort"]=cohort
    d=d.merge(W[W["cohort"]==cohort][["ID","wave","weight"]],on=["ID","wave"],how="left")
    return d[["cohort","ID","wave","age","ageg","pain","dx","tx","weight"]]

j={}
j["CHARLS"]=build_joint("CHARLS", os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),
    ["ID","wave","age","arthre","rxarthr_c","da042s2","da042s3","da042s4","da042s5","da042s8","da042s9",
     "da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
    ["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s14","da042s15"],
    "arthre","rxarthr_c","age",[1,2,3,4], id_col="ID")
j["KLoSA"]=build_joint("KLoSA", os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),
    ["pid","wave","agey","arthre","rxarthr","pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
    ["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"],
    "arthre","rxarthr","agey", id_col="pid")
j["LASI"]=build_joint("LASI", os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),
    ["prim_key","r1agey","r1arthre","r1rxarthr","r1jointp"], None,"r1arthre","r1rxarthr","r1agey", None, pain_single="r1jointp", id_col="prim_key")
j["LASI"]["wave"]=1
j["SHARE"]=build_joint("SHARE", os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),
    ["mergeid","wave","agey","arthre","rxpain","pain_s"], None,"arthre","rxpain","agey",[1,2,4], pain_single="pain_s", id_col="mergeid")

for c in PANEL_B:
    rowsB.append(summarize(c,j[c],"B_joint_pain"))

# ---------------- 表 C：四状态 × ADL ----------------
for c,d in frames.items():
    dd=d[d["weight"].notna() & d["pain"].notna()].copy()
    if len(dd)==0: continue
    if dd["tx"].notna().any() and c != "SHARE":
        dd["state"]=np.select(
            [dd["pain"]==0, (dd["pain"]==1)&(dd["dx"]==0), (dd["pain"]==1)&(dd["dx"]==1)&(dd["tx"]==0), (dd["pain"]==1)&(dd["dx"]==1)&(dd["tx"]==1)],
            ["0_no_pain","1_pain_undx","2_pain_dx_untx","3_pain_dx_tx"], default=None)
    else:
        dd["state"]=np.select(
            [dd["pain"]==0, (dd["pain"]==1)&(dd["dx"]==0), (dd["pain"]==1)&(dd["dx"]==1)],
            ["0_no_pain","1_pain_undx","2_pain_dx"], default=None)
    for s,g in dd.groupby("state"):
        rows4.append(dict(cohort=c, state=s, n=len(g),
                          pct_of_pain=round(len(g)/max(1,(dd["pain"]==1).sum())*100,1) if s!="0_no_pain" else np.nan,
                          adl_crude=round(wrate(g,"adl","weight")*100,1) if g["adl"].notna().any() else np.nan))

# ---------------- 表 D：缺失矩阵 ----------------
for c,d in frames.items():
    pat = d.assign(p=d["pain"].isna(), x=d["dx"].isna(), t=d["tx"].isna())
    tab = pat.groupby(["p","x","t"]).size().reset_index(name="n")
    tab["cohort"]=c
    rowsM.append(tab)

pd.DataFrame(rowsA).to_csv(os.path.join(OUT,"17_cascade_tierA.csv"),index=False,encoding="utf-8-sig")
pd.DataFrame(rowsB).to_csv(os.path.join(OUT,"17_cascade_tierB.csv"),index=False,encoding="utf-8-sig")
pd.DataFrame(rows4).to_csv(os.path.join(OUT,"17_four_state.csv"),index=False,encoding="utf-8-sig")
pd.concat(rowsM).to_csv(os.path.join(OUT,"17_missing.csv"),index=False,encoding="utf-8-sig")

print(pd.DataFrame(rowsA).to_string(index=False))
print()
print(pd.DataFrame(rowsB).to_string(index=False))
print()
print(pd.DataFrame(rows4).to_string(index=False))
