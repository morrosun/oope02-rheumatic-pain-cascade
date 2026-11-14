# -*- coding: utf-8 -*-
"""23: A3 汇总表"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pandas as pd, numpy as np
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
res=pd.read_csv(os.path.join(OUT,"22_outcome_gradient.csv"))
pool=pd.read_csv(os.path.join(OUT,"22_pooled.csv"))

LBL={"s1":"疼痛·未诊断","s2":"疼痛·已诊断未治疗","s3":"疼痛·已诊断在治疗","dx_vs_undx":"已诊断 vs 未诊断（疼痛者中）"}

def md(df):
    cols=list(df.columns)
    o=["| "+" | ".join(cols)+" |","|"+"|".join(["---"]*len(cols))+"|"]
    for _,r in df.iterrows():
        o.append("| "+" | ".join(str(r[c]) for c in cols)+" |")
    return "\n".join(o)

def fmt(r):
    return f"{r['OR']:.2f}（{r['lo']:.2f}–{r['hi']:.2f}）"

L=[]
L.append("# A3：疼痛、诊断与 ADL 受限\n")
L.append("**日期**：2026-11-14　**脚本**：`out/22_outcome_gradient.py` → `23_outcome_table.py`\n")
L.append("**结局**：ADL 受限（有任意 1 项 ADL 困难/需要帮助）。**模型**：每库 logistic，个体聚类稳健，DL 合并。\n")
L.append("**协变量**：年龄、性别、教育、财富五分位、城乡、共病数、吸烟；M2 另加疼痛严重度 z 分（`painlv`/`painhlv`/部位数/频率）。\n")

L.append("\n## 表 1　四状态 vs 无痛（M1，调整社会人口学与共病，不调整严重度）\n")
rows=[]
for panel in ["A","B"]:
    g=pool[(pool.panel==panel)&(pool.model=="M1_vs_nopain")]
    for t in ["s1","s2","s3"]:
        r=g[g.term==t]
        if len(r)==0: continue
        r=r.iloc[0]
        rows.append({"面板":panel,"状态":LBL[t],"合并 OR":fmt(r),"k":int(r["k"]),"I²%":r["I2"]})
L.append(md(pd.DataFrame(rows)))

L.append("\n## 表 2　疼痛者中「已诊断 vs 未诊断」（M2）\n")
rows=[]
for panel in ["A","B"]:
    for model,lab in [("M2_no_sev","未调整疼痛严重度"),("M2_adj_sev","调整疼痛严重度")]:
        g=pool[(pool.panel==panel)&(pool.model==model)]
        if len(g)==0: continue
        r=g.iloc[0]
        rows.append({"面板":panel,"模型":lab,"合并 OR":fmt(r),"k":int(r["k"]),"I²%":r["I2"]})
L.append(md(pd.DataFrame(rows)))

L.append("\n## 逐库 OR（看一致性）\n")
piv=res[res.model.isin(["M1_vs_nopain","M2_adj_sev"])].copy()
piv["标签"]=piv["model"].str.replace("M1_vs_nopain","M1_").str.replace("M2_adj_sev","M2_")+piv["term"].map(lambda x:"" if x=="dx_vs_undx" else "_"+x)
t3=piv.pivot_table(index=["panel","标签"],columns="cohort",values="OR").round(2).reset_index()
L.append(md(t3.fillna("—")))

L.append("\n## 三条结论（第二条与预期相反，必须正视）\n")
L.append("1. **疼痛本身就是功能受限的强标记，且“未诊断”并不轻**：与无痛者相比，疼痛但从未诊断者的 ADL 受限 OR 为 **3.22（面板 A）/ 1.82（面板 B）**；已诊断者是 4.69/2.33，已诊断在治疗者 6.71/2.94。**未诊断组已经承担了诊断组约 60–70% 的额外功能负担**——冰山不是良性的。")
L.append("2. **但“已诊断”与更高的 ADL 受限相关，而不是更低**：疼痛者中，已诊断 vs 未诊断的 OR 为 **1.61（A）/ 1.46（B）**（已调整疼痛严重度）。这不是“诊断有害”，而是**横断面设计无法把“诊断”与“病更重/更常就医”分开**——ADL 差的人更可能去看医生并被诊断（反向因果），或诊断组确实病变更重（严重度测量不足）。")
L.append("3. **因此 A3 不能回答“漏诊是否造成损害”**。它能回答的是：漏诊人群的功能负担已经很高（结论 1），而这正说明**需要 A4 的纵向设计**——基线未诊断的疼痛者，两年内获得诊断的那部分人，ADL 轨迹是否优于始终未诊断者。")

L.append("\n## 为什么调整严重度后诊断组仍然更差\n")
L.append("- **严重度测量太粗**：`painlv` 只有 0–3 四级、CHARLS 用部位数代替，无法刻画“同一疼痛强度下病变更重”的差别。")
L.append("- **诊断是就医的结果**：ADL 差 → 就医多 → 更容易被诊断，这条通路在横断面里与“诊断改善功能”方向相反，会互相抵消。")
L.append("- **适应证混杂**：治疗组 ADL 最差（M1 的 s3 最高），就是最直接的证据——病重的人才吃药。")

L.append("\n## 口径与局限\n")
L.append("- **KLoSA 是异常点**：ADL 需要帮助者仅 2.9%，四状态梯度几乎平坦（s1 1.41、s2 0.94、s3 1.48），与其他六库不同；建议主分析报六库、KLoSA 单列。")
L.append("- **SHARE 无抗风湿治疗变量**，四状态只有三段；ELSA 同。")
L.append("- **CHARLS 只到 w1–w4**；LASI 单波。")
L.append("- **ADL 口径不完全一致**：KLoSA/SHARE 由条目合成（任一困难/需帮助），其余用总分 ≥1；ELSA 的 `adltot6` 与 HRS 的 `adl6a` 计分方式需在附录列明。")
L.append("- **横断面**：暴露（诊断）与结局（ADL）同时测量，不能做因果解释；所有 OR 只能读作关联。")

with open(os.path.join(OUT,"结局梯度.md"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("WROTE 结局梯度.md")
