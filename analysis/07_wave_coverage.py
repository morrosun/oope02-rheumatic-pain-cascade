# -*- coding: utf-8 -*-
"""Round 5 probe: wave-level coverage of the cascade and of preventive-service variables."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import pyreadstat
import pandas as pd

lines=[]
def w(s=""):
    lines.append(s); print(s, flush=True)

def cov(df, var, by="wave"):
    if var not in df.columns: 
        w(f"  {var}: NOT PRESENT"); return
    t = df.groupby(by)[var].apply(lambda s: pd.Series({"n": s.notna().sum(), "pos": (s==1).sum()}))
    w(f"  {var}:"); w(t.to_string())

# CHARLS
p = os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta")
df,_ = pyreadstat.read_dta(p, usecols=["ID","wave","arthre","rxarthr","rxarthr_c","da042s2","da042s12"])
w("## CHARLS wave coverage"); cov(df,"arthre"); cov(df,"rxarthr"); cov(df,"rxarthr_c")
w("  among arthre==1 by wave: rxarthr non-missing share")
sub=df[df["arthre"]==1]
w(sub.groupby("wave").apply(lambda g: pd.Series({"n_dx":len(g),"rx_nonmiss":g["rxarthr"].notna().sum(),"rx1":(g["rxarthr"]==1).sum()})).to_string())

# HRS
p = os.path.join(DL, "3. HRS  美国/HRS_美国/Working_data/hrs.dta")
df,_ = pyreadstat.read_dta(p, usecols=["hhidpn","wave","arthre","rxarthr","rxpain","painfr","paina",
                                       "flusht","pneushte","shnglshte","cholst","mammog","papsm","prost"])
w("\n## HRS wave coverage"); cov(df,"rxarthr"); cov(df,"rxpain"); cov(df,"paina")
w("  preventive:"); cov(df,"flusht"); cov(df,"pneushte"); cov(df,"shnglshte"); cov(df,"cholst"); cov(df,"mammog"); cov(df,"prost")

# SHARE
p = os.path.join(DL, "7.SHARE 欧洲/SHARE_欧洲/Working_data/share.dta")
df,_ = pyreadstat.read_dta(p, usecols=["mergeid","wave","arthre","pain_s","painfr","painlv","rxpain","rxinflm","rxosteo"])
w("\n## SHARE wave coverage"); cov(df,"pain_s"); cov(df,"painfr"); cov(df,"rxpain"); cov(df,"rxinflm")

# LASI / MHAS / KLoSA preventive
p = os.path.join(DL, "5. LASI 印度/LASI_印度/Working_data/lasi.dta")
df,_ = pyreadstat.read_dta(p, usecols=["prim_key","r1flushte","r1pneushte","r1cholst","r1mammog","r1papsm","r1dentst1y"])
w("\n## LASI coverage (single wave)")
for v in df.columns[1:]: w(f"  {v}: n={df[v].notna().sum()} pos={(df[v]==1).sum()}")

p = os.path.join(DL, "6. MHAS 墨西哥/MHAS_墨西哥/Working_data/mhas.dta")
df,_ = pyreadstat.read_dta(p, usecols=["rahhidnp","wave","flusht","cholst","mammog","papsm","prost","dentst1y"])
w("\n## MHAS coverage"); cov(df,"flusht"); cov(df,"cholst"); cov(df,"mammog"); cov(df,"prost")

p = os.path.join(DL, "4. KLoSA 韩国/KLoSA_韩国/Working_data/klosa.dta")
df,_ = pyreadstat.read_dta(p, usecols=["pid","wave","arthre","rxarthr","painfr","paina","dentst","homcar","trdmed2y"])
w("\n## KLoSA coverage"); cov(df,"rxarthr"); cov(df,"dentst"); cov(df,"homcar"); cov(df,"trdmed2y")

p = os.path.join(DL, "2. ELSA 英国/Working_data/elsa.dta")
df,_ = pyreadstat.read_dta(p, usecols=["idauniqc","wave","arthre","painfr","painlv","mammog","prost","cholst","dentalh","noteeth"])
w("\n## ELSA coverage"); cov(df,"mammog"); cov(df,"prost"); cov(df,"dentalh"); cov(df,"noteeth")

with open(os.path.join(OUT,"07_wave_coverage.txt"),"w",encoding="utf-8") as f:
    f.write("\n".join(lines))
print("WROTE 07_wave_coverage.txt")
