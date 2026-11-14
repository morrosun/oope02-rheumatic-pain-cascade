# -*- coding: utf-8 -*-
"""36: 第二遍复盘 —— 检查新代码引入的问题 + 文件间数字一致性"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
L=[]; w=lambda s:(L.append(str(s)),print(s,flush=True))

w("="*74); w("【1】群标识缺失与唯一性（34/33 的集群单位）"); w("="*74)
CL=[("CHARLS",os.path.join(DL,"1. CHARLS  中国","CHARLS_中国","Working_data","charls.dta"),"communityID","ID"),
    ("ELSA",os.path.join(DL,"2. ELSA 英国","Working_data","elsa.dta"),"hhidc","idauniqc"),
    ("HRS",os.path.join(DL,"3. HRS  美国","HRS_美国","Working_data","hrs.dta"),"hhhid","hhidpn"),
    ("KLoSA",os.path.join(DL,"4. KLoSA 韩国","KLoSA_韩国","Working_data","klosa.dta"),"hhidc","pid"),
    ("LASI",os.path.join(DL,"5. LASI 印度","LASI_印度","Working_data","lasi.dta"),"hhid","prim_key"),
    ("MHAS",os.path.join(DL,"6. MHAS 墨西哥","MHAS_墨西哥","Working_data","mhas.dta"),"hhidc","rahhidnp"),
    ("SHARE",os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),"hhidc","mergeid")]
for c,p,clu,idc in CL:
    df,_=pyreadstat.read_dta(p,usecols=[idc,clu])
    s=df[clu]; miss=int(s.isna().sum())
    # 缺失是否会被 astype(str) 变成同一个 "nan" 群
    w(f"{c:7s} 行={len(df):7d} 群={s.nunique():7d} 群缺失={miss:6d} ({miss/len(df)*100:5.1f}%)  "
      f"每群中位人数={df.groupby(clu).size().median():.1f}")

w("\n"+"="*74); w("【2】SHARE 的 painfr / pain_s 是否按波可用（31 的 A 面板未显式过滤波）"); w("="*74)
sh,mt=pyreadstat.read_dta(os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta"),
                          usecols=["mergeid","wave","painfr","pain_s","painlv","agey"])
g=sh.groupby("wave").agg(n=("mergeid","size"),painfr=("painfr","count"),pain_s=("pain_s","count"),
                         painlv=("painlv","count"))
w(g.to_string())
w("\n→ 若 painfr 在 w1/2/4 有值，则 31 的 A 面板会混入关节痛时代的人，必须加波过滤。")

w("\n"+"="*74); w("【3】关键比值反算（防止文字与数字不一致）"); w("="*74)
ci=pd.read_csv(os.path.join(OUT,"33_std_ci.csv"))
pv=pd.read_csv(os.path.join(OUT,"31_pooled_v2.csv"))
def get(t,m,co):
    r=ci[(ci.tier==t)&(ci.metric==m)&(ci.cohort==co)]
    return (r.est.iloc[0],r.lo.iloc[0],r.hi.iloc[0]) if len(r) else (np.nan,)*3
for t,pairs,label in [("A",[("HRS","MHAS")],"面板A 诊断"),("B",[("CHARLS","LASI")],"面板B 诊断")]:
    for hi,lo in pairs:
        a=get(t,"dx_given_pain",hi); b=get(t,"dx_given_pain",lo)
        w(f"{label}: {hi} {a[0]:.2f}({a[1]:.2f}-{a[2]:.2f}) vs {lo} {b[0]:.2f}({b[1]:.2f}-{b[2]:.2f}) "
          f"= {a[0]/b[0]:.2f} 倍；区间重叠={not (a[1]>b[2] or b[1]>a[2])}")
for t,pairs,label in [("A",[("HRS","MHAS")],"面板A 治疗"),("B",[("KLoSA","CHARLS")],"面板B 治疗")]:
    for hi,lo in pairs:
        a=get(t,"tx_given_dx_pain",hi); b=get(t,"tx_given_dx_pain",lo)
        w(f"{label}: {hi} {a[0]:.2f} vs {lo} {b[0]:.2f} = {a[0]/b[0]:.2f} 倍；区间重叠={not (a[1]>b[2] or b[1]>a[2])}")

w("\n【31 内部】未诊断组承担已诊断组的比例：")
for panel in ["A","B"]:
    s1=pv[(pv.panel==panel)&(pv["mode"]=="three")&(pv.term=="s1")].OR.iloc[0]
    s2=pv[(pv.panel==panel)&(pv["mode"]=="three")&(pv.term=="s2")].OR.iloc[0]
    w(f"  面板{panel}: {s1:.3f}/{s2:.3f} = {s1/s2*100:.1f}%")

w("\n"+"="*74); w("【4】A4 重跑后与旧值的差异"); w("="*74)
inc=pd.read_csv(os.path.join(OUT,"24_incident_pooled.csv"))
w(inc[inc.model=="M_adl"][["term","OR","lo","hi","I2"]].to_string(index=False))

w("\n"+"="*74); w("【5】结局梯度中的 s2 是否已彻底分离三/四状态"); w("="*74)
res=pd.read_csv(os.path.join(OUT,"31_outcome_v2.csv"))
for panel in ["A","B"]:
    for c in (["HRS","MHAS","ELSA","SHARE"] if panel=="A" else ["CHARLS","KLoSA","LASI","SHARE"]):
        modes=res[(res.panel==panel)&(res.cohort==c)]["mode"].unique().tolist()
        w(f"  {panel}/{c:7s} modes={modes}")
w("\n→ 三状态(SHARE/ELSA)与四状态库必须分开合并；检查 pooled 的 cohorts 字段：")
w(pv[["panel","mode","term","k","cohorts"]].to_string(index=False))

open(os.path.join(OUT,"36_verify.txt"),"w",encoding="utf-8").write("\n".join(L))
print("\nWROTE 36_verify.txt")
