# -*- coding: utf-8 -*-
"""Descriptive cascade counts (unweighted, feasibility only) for the recommended direction."""
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
rows=[]; lines=[]
def w(s=""):
    lines.append(str(s)); print(s, flush=True)

def add(cohort, wave, n_pain, n_dx_among_pain, n_dx, n_rx_among_dx, n_rx_avail):
    rows.append(dict(cohort=cohort, wave=wave, n_pain=n_pain, dx_among_pain=n_dx_among_pain,
                     n_dx=n_dx, rx_among_dx=n_rx_among_dx, rx_avail=n_rx_avail))

def anypos(df, cols):
    cols=[c for c in cols if c in df.columns]
    if not cols: return pd.Series(np.nan, index=df.index)
    sub=df[cols]
    out=(sub==1).any(axis=1).astype(float)
    out[sub.isna().all(axis=1)]=np.nan
    return out

# ---------- CHARLS ----------
p=os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta")
df,_=pyreadstat.read_dta(p, usecols=["ID","wave","arthre","rxarthr",
    "da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s15"])
df["pain"]=anypos(df,["da042s2","da042s3","da042s4","da042s5","da042s8","da042s9","da042s10","da042s11","da042s12","da042s13","da042s15"])
for wv,g in df.groupby("wave"):
    pain=g[g["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=g[g["arthre"]==1]
    add("CHARLS",wv,len(pain),(dx["arthre"]==1).sum(),len(dxall),(dxall["rxarthr"]==1).sum(),dxall["rxarthr"].notna().sum())
w("CHARLS done")

# ---------- KLoSA ----------
p=os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Working_data/klosa.dta")
df,_=pyreadstat.read_dta(p, usecols=["pid","wave","arthre","rxarthr",
    "pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"])
df["pain"]=anypos(df,["pain_2","pain_3","pain_4","pain_5","pain_8","pain_9","pain_10","pain_11","pain_12","pain_13"])
for wv,g in df.groupby("wave"):
    pain=g[g["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=g[g["arthre"]==1]
    add("KLoSA",wv,len(pain),(dx["arthre"]==1).sum(),len(dxall),(dxall["rxarthr"]==1).sum(),dxall["rxarthr"].notna().sum())
w("KLoSA done")

# ---------- LASI ----------
p=os.path.join(DL, "5. LASI 印度/LASI_印度/Working_data/lasi.dta")
df,_=pyreadstat.read_dta(p, usecols=["prim_key","r1arthre","r1rxarthr","r1jointp"])
df=df.rename(columns={"r1arthre":"arthre","r1rxarthr":"rxarthr","r1jointp":"pain"})
pain=df[df["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=df[df["arthre"]==1]
add("LASI",1,len(pain),len(dx),len(dxall),(dxall["rxarthr"]==1).sum(),dxall["rxarthr"].notna().sum())
w("LASI done")

# ---------- MHAS ----------
p=os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Working_data/mhas.dta")
df,_=pyreadstat.read_dta(p, usecols=["rahhidnp","wave","arthre","rxarthr","painfr"])
df["pain"]=df["painfr"]
for wv,g in df.groupby("wave"):
    pain=g[g["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=g[g["arthre"]==1]
    add("MHAS",wv,len(pain),len(dx),len(dxall),(dxall["rxarthr"]==1).sum(),dxall["rxarthr"].notna().sum())
w("MHAS done")

# ---------- HRS ----------
p=os.path.join(DL, "3. HRS  美国/HRS_美国/Working_data/hrs.dta")
df,_=pyreadstat.read_dta(p, usecols=["hhidpn","wave","arthre","rxarthr","painfr"])
df["pain"]=df["painfr"]
for wv,g in df.groupby("wave"):
    pain=g[g["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=g[g["arthre"]==1]
    add("HRS",wv,len(pain),len(dx),len(dxall),(dxall["rxarthr"]==1).sum(),dxall["rxarthr"].notna().sum())
w("HRS done")

# ---------- SHARE ----------
p=os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Working_data/share.dta")
df,_=pyreadstat.read_dta(p, usecols=["mergeid","wave","country","arthre","pain_s","painfr","rxpain"])
df["pain"]=df["pain_s"].where(df["pain_s"].notna(), df["painfr"])
df.loc[(df["pain_s"].isna())&(df["painfr"].isna()),"pain"]=np.nan
for wv,g in df.groupby("wave"):
    pain=g[g["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=g[g["arthre"]==1]
    add("SHARE",wv,len(pain),len(dx),len(dxall),(dxall["rxpain"]==1).sum(),dxall["rxpain"].notna().sum())
w("SHARE done (rxpain = any analgesic weekly, NOT arthritis-specific)")

# ---------- ELSA ----------
p=os.path.join(DL, "2. ELSA 英国/Working_data/elsa.dta")
df,_=pyreadstat.read_dta(p, usecols=["idauniqc","wave","arthre","painfr"])
df["pain"]=df["painfr"]
for wv,g in df.groupby("wave"):
    pain=g[g["pain"]==1]; dx=pain[pain["arthre"]==1]; dxall=g[g["arthre"]==1]
    add("ELSA",wv,len(pain),len(dx),len(dxall),np.nan,np.nan)
w("ELSA done (no arthritis-medication variable)")

tab=pd.DataFrame(rows)
tab["pct_dx_among_pain"]=(tab["dx_among_pain"]/tab["n_pain"]*100).round(1)
tab["pct_rx_among_dx"]=(tab["rx_among_dx"]/tab["rx_avail"]*100).round(1)
w("\n=== per-wave cascade counts ===")
w(tab.to_string(index=False))

w("\n=== cohort totals (pain stage pooled over waves) ===")
agg=tab.groupby("cohort").agg(n_pain=("n_pain","sum"), dx_among_pain=("dx_among_pain","sum"),
                              n_dx=("n_dx","sum"), rx_among_dx=("rx_among_dx","sum"), rx_avail=("rx_avail","sum"))
agg["pct_dx_among_pain"]=(agg["dx_among_pain"]/agg["n_pain"]*100).round(1)
agg["pct_rx_among_dx"]=(agg["rx_among_dx"]/agg["rx_avail"]*100).round(1)
w(agg.to_string())

tab.to_csv(os.path.join(OUT,"08_cascade_counts.csv"),index=False,encoding="utf-8-sig")
with open(os.path.join(OUT,"08_cascade_counts.txt"),"w",encoding="utf-8") as f:
    f.write("\n".join(lines))
print("WROTE 08_cascade_counts.*")
