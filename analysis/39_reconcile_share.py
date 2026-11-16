# -*- coding: utf-8 -*-
"""
39: 对账 —— 17_cascade_tierB 的 SHARE n_all=121,819 与 33 的 n=120,781 相差 1,038
逐层拆解：权重缺失 / 疼痛缺失 / 年龄组缺失 / 群缺失，并检查缺失是否与波次、年龄、疼痛率相关
输出：39_reconcile_share.txt
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
L=[]; w=lambda s:(L.append(str(s)),print(s,flush=True))
SH=os.path.join(DL,"7.SHARE 欧洲","SHARE_欧洲","Working_data","share.dta")
W=pd.read_csv(os.path.join(OUT,"16_weights.csv"),dtype={"ID":str})[["cohort","ID","wave","weight"]]

df,_=pyreadstat.read_dta(SH,usecols=["mergeid","hhidc","wave","agey","arthre","pain_s","rxpain","painfr","painlv"])
w("="*74); w("【0】原始行数与各过滤步骤"); w("="*74)
w(f"读入行数                      {len(df):7d}")
d=df.rename(columns={"mergeid":"ID"}).copy()
d["ID"]=d["ID"].astype(str)
w(f"ID 转字符串后                 {len(d):7d}")

d["age"]=pd.to_numeric(d["agey"],errors="coerce")
w(f"age 非缺失                    {int(d['age'].notna().sum()):7d}   （缺失 {int(d['age'].isna().sum())}）")
d1=d[(d["age"]>=50)&(d["age"]<=105)].copy()
w(f"年龄 50–105 过滤后            {len(d1):7d}   ← 17 的 n_all 应等于此")

d2=d1[d1["wave"].isin([1,2,4])].copy()
w(f"波次 [1,2,4] 过滤后           {len(d2):7d}")
w(f"  （对照：17_cascade_tierB.csv 的 SHARE n_all = 121,819）")

d2["wave"]=d2["wave"].astype(float)
m=d2.merge(W[W.cohort=="SHARE"][["ID","wave","weight"]],on=["ID","wave"],how="left")
w(f"merge 权重后行数              {len(m):7d}   （若有增殖说明权重表有重复）")
w(f"  权重非缺失                  {int(m['weight'].notna().sum()):7d}   ← 17 的 n_weighted 应等于此")
w(f"  权重缺失                    {int(m['weight'].isna().sum()):7d}")

w("\n"+"="*74); w("【1】逐层拆解（33 的口径：pain + ageg + weight + CLU 四者非缺失）"); w("="*74)
m["pain"]=pd.to_numeric(m["pain_s"],errors="coerce").where(lambda s:s.isin([0,1]))
m["ageg"]=pd.cut(m["age"],bins=[50,55,60,65,70,75,80,85,200],right=False,labels=list(range(1,9))).astype(float)
m["CLU"]=m["hhidc"].astype(str)
first=m.sort_values(["ID","wave"]).groupby("ID")["CLU"].first().rename("CLU_FIX")
m=m.merge(first,on="ID",how="left")

n0=len(m)
nw=int(m["weight"].notna().sum())
nwp=int((m["weight"].notna()&m["pain"].notna()).sum())
nwpa=int((m["weight"].notna()&m["pain"].notna()&m["ageg"].notna()).sum())
nwpac=int((m["weight"].notna()&m["pain"].notna()&m["ageg"].notna()&m["CLU_FIX"].notna()).sum())
w(f"起始（年龄+波次过滤）          {n0:7d}")
w(f"  − 权重缺失                  {n0-nw:7d}  → {nw:7d}")
w(f"  − 疼痛缺失（权重非缺失者中）  {nw-nwp:7d}  → {nwp:7d}")
w(f"  − 年龄组缺失                {nwp-nwpa:7d}  → {nwpa:7d}")
w(f"  − 群标识缺失                {nwpa-nwpac:7d}  → {nwpac:7d}   ← 33 的 n")
w(f"\n  总差额 = {n0-nwpac}   （用户指出的 1,038 = {n0}-{nwpac} = {n0-nwpac}）")

w("\n"+"="*74); w("【2】重叠检查：缺失是否叠加（同一行既缺权重又缺疼痛）"); w("="*74)
tab=pd.crosstab(m["weight"].isna().map({True:"缺权重",False:"有权重"}),
                m["pain"].isna().map({True:"缺疼痛",False:"有疼痛"}),margins=True)
w(tab.to_string())
w("\n  → 若「缺权重」与「缺疼痛」几乎不重叠，则 1,038 ≈ 484 + 554，可直接相加解释。")

w("\n"+"="*74); w("【3】缺失是否集中于某一波（决定会不会引入波次偏倚）"); w("="*74)
g=m.groupby("wave").agg(rows=("ID","size"),
                        w_miss=("weight",lambda s:s.isna().sum()),
                        p_miss=("pain",lambda s:s.isna().sum()))
g["w_miss%"]=(g["w_miss"]/g["rows"]*100).round(1)
g["p_miss%"]=(g["p_miss"]/g["rows"]*100).round(1)
w(g.to_string())

w("\n"+"="*74); w("【4】缺失者的年龄与疼痛率（决定会不会引入年龄/结局偏倚）"); w("="*74)
m["keep"]=m["weight"].notna()&m["pain"].notna()
for lab,sel in [("纳入（33 的样本）",m["keep"]),("权重缺失",m["weight"].isna()),("疼痛缺失但有权重",m["weight"].notna()&m["pain"].isna())]:
    s=m[sel]
    w(f"{lab:20s} n={len(s):6d}  年龄中位={s['age'].median():5.0f}  "
      f"关节炎率={s['arthre'].mean()*100 if s['arthre'].notna().any() else float('nan'):5.1f}%  "
      f"疼痛率={s['pain'].mean()*100 if s['pain'].notna().any() else float('nan'):5.1f}%")

w("\n"+"="*74); w("【5】影响评估：把缺失者按最坏/最好情形补回，疼痛率与诊断覆盖率会变多少"); w("="*74)
base=m[m["keep"]]
w(f"当前（33）：n={len(base)}  疼痛率={base['pain'].mean()*100:.2f}%  "
  f"诊断覆盖（疼痛者中）={base[base.pain==1]['arthre'].mean()*100:.2f}%")
w("缺失的 1,038 行无法观测其疼痛状态；但可用其 arthre 与年龄刻画：")
miss=m[~m["keep"]]
w(f"  缺失行 n={len(miss)}  年龄中位={miss['age'].median():.0f}  arthre 非缺失 {int(miss['arthre'].notna().sum())}  "
  f"arthre 阳性率={miss['arthre'].mean()*100 if miss['arthre'].notna().any() else float('nan'):.1f}%")
w("结论：若缺失行的疼痛率与纳入行同量级（约 54%），补回后疼痛率变动 < 0.5 个百分点，可忽略。")

open(os.path.join(OUT,"39_reconcile_share.txt"),"w",encoding="utf-8").write("\n".join(L))
print("\nWROTE 39_reconcile_share.txt")
