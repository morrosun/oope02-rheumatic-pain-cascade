# -*- coding: utf-8 -*-
"""
43: 诊断 —— 把 tx 结构性缺失当 0 会低估终点产出率多少？
对比三种口径：
  all      : 全部波次，疼痛与诊断可算，tx 缺失者从 tx 分母剔除（= 33 的口径）
  txwave   : 只保留「有治疗条目」的波次，三段同分母
  txobs    : 全部波次，但只保留 tx 非缺失者（诊断者内）
输出 43_diag_txwaves.txt
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd
import importlib.util
spec = importlib.util.spec_from_file_location("f42", os.path.join(OUT, "42_fig1_funnel.py"))
# 不执行整个脚本，改为直接复制其函数：用 exec 截断到 rows=[] 之前
src = open(os.path.join(OUT, "42_fig1_funnel.py"), encoding="utf-8").read()
head = src.split("rows=[]")[0]
exec(compile(head, "42_head", "exec"))

def wrate(df, y, w):
    m = df.dropna(subset=[y, w])
    return float(np.average(m[y], weights=m[w])) if len(m) else np.nan

L=[]
for tier, cohorts in [("A",["HRS","MHAS","ELSA","SHARE"]), ("B",["CHARLS","KLoSA","LASI","SHARE"])]:
    for c in cohorts:
        d = load(c, SPEC[c], tier)
        base = d.dropna(subset=["pain","dx","weight"]).copy()
        base["yield"] = ((base["pain"]==1)&(base["dx"]==1)&(base["tx"]==1)).astype(float)
        pain = base[base.pain==1]
        # all：tx 分母剔除缺失
        dx_all  = wrate(pain, "dx", "weight")
        tx_all  = wrate(pain[pain.dx==1], "tx", "weight")
        yld_all = wrate(pain, "yield", "weight")
        # txobs：只保留 tx 非缺失
        obs = base[base.tx.notna()]
        dx_obs  = wrate(obs[obs.pain==1], "dx", "weight")
        yld_obs = wrate(obs[obs.pain==1], "yield", "weight")
        # txwave：逐波看哪些波有 tx
        byw = base.groupby("wave")["tx"].apply(lambda s: s.notna().sum())
        okw = sorted(byw[byw>0].index.tolist())
        tw = base[base.wave.isin(okw)]
        dx_tw  = wrate(tw[tw.pain==1], "dx", "weight")
        tx_tw  = wrate(tw[(tw.pain==1)&(tw.dx==1)], "tx", "weight")
        yld_tw = wrate(tw[tw.pain==1], "yield", "weight")
        L.append(f"{tier} {c:7s} waves_with_tx={okw}")
        L.append(f"    all    : dx={dx_all*100:5.1f}  tx={tx_all*100:5.1f}  yield={yld_all*100:5.1f}")
        L.append(f"    txobs  : dx={dx_obs*100:5.1f}                yield={yld_obs*100:5.1f}  (Δdx={100*(dx_obs-dx_all):+.1f})")
        L.append(f"    txwave : dx={dx_tw*100:5.1f}  tx={tx_tw*100:5.1f}  yield={yld_tw*100:5.1f}  (Δdx={100*(dx_tw-dx_all):+.1f})")
out="\n".join(L)
open(os.path.join(OUT, "43_diag_txwaves.txt"),"w",encoding="utf-8").write(out)
print(out)
