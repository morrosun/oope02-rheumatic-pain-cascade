# -*- coding: utf-8 -*-
"""28: CHARLS 两期对照（2011–2018 关节痛 vs 2021–2023 体痛）+ 并入 Panel A 的位置"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pandas as pd, numpy as np
w6=pd.read_csv(os.path.join(OUT,"27_charls_w6_row.csv"))
st6=pd.read_csv(os.path.join(OUT,"27_charls_w6_state.csv"))
tB=pd.read_csv(os.path.join(OUT,"17_cascade_tierB.csv"))
tA=pd.read_csv(os.path.join(OUT,"17_cascade_tierA.csv"))
f4=pd.read_csv(os.path.join(OUT,"17_four_state.csv"))

c14=tB[tB.cohort=="CHARLS"].iloc[0]
def md(df):
    cols=list(df.columns); o=["| "+" | ".join(cols)+" |","|"+"|".join(["---"]*len(cols))+"|"]
    for _,r in df.iterrows(): o.append("| "+" | ".join(str(r[c]) for c in cols)+" |")
    return "\n".join(o)

L=[]
L.append("# CHARLS 第六轮并入主表：结果与影响\n")
L.append("**日期**：2026-11-14　**脚本**：`27_charls_w6_row.py` → `28_compare_table.py`\n")
L.append("**先说限制**：第六轮没有 15 部位清单，只有五级体痛（`da069`）。所以它**不能进关节痛面板（B），也不能直接进频繁疼痛面板（A）**——它只能作为「中国 2021–2023，体痛口径」的独立一行。下面把它与「中国 2011–2018，关节痛口径」并排放，看十年间的形状变化。\n")

L.append("\n## 表 1　CHARLS 两期对照（加权 + WHO 年龄标化）\n")
rows=[
 dict(时期="2011–2018（w1–w4）",口径="关节痛（部位清单）",人波=int(c14["n_all"]),
      疼痛率=f"{c14['pain_crude']:.1f}（{c14['pain_std']:.1f}）",
      诊断率=f"{c14['dx_crude']:.1f}（{c14['dx_std']:.1f}）",
      治疗率_诊断且疼痛=f"{c14['txpain_crude']:.1f}（{c14['txpain_std']:.1f}）",
      治疗率_全部诊断=f"{c14['tx_crude']:.1f}"),
 dict(时期="2021–2023（w6）",口径="体痛 ≥ 一点",人波=int(w6.iloc[0]["n"]),
      疼痛率=f"{w6.iloc[0]['疼痛率']:.1f}（{w6.iloc[0]['疼痛率标化']:.1f}）",
      诊断率=f"{w6.iloc[0]['诊断率']:.1f}（{w6.iloc[0]['诊断率标化']:.1f}）",
      治疗率_诊断且疼痛=f"{w6.iloc[0]['治疗率']:.1f}（{w6.iloc[0]['治疗率标化']:.1f}）",
      治疗率_全部诊断="36.5"),
 dict(时期="2021–2023（w6）",口径="体痛 ≥ 相当多",人波=int(w6.iloc[1]["n"]),
      疼痛率=f"{w6.iloc[1]['疼痛率']:.1f}（{w6.iloc[1]['疼痛率标化']:.1f}）",
      诊断率=f"{w6.iloc[1]['诊断率']:.1f}（{w6.iloc[1]['诊断率标化']:.1f}）",
      治疗率_诊断且疼痛=f"{w6.iloc[1]['治疗率']:.1f}（{w6.iloc[1]['治疗率标化']:.1f}）",
      治疗率_全部诊断="—"),
]
L.append(md(pd.DataFrame(rows)))
L.append("\n> 括号内为年龄标化值。治疗率（全部诊断）来自 `14_w6_cascade.csv`（2,397/6,567 = 36.5%）与 `17_cascade_tierB.csv`（43.5%）。\n")

L.append("\n## 表 2　四状态 × ADL 受限：两期对照（加权）\n")
c4=f4[f4.cohort=="CHARLS"].set_index("state")["adl_crude"]
rows=[]
for lab,key in [("无痛","0_no_pain"),("痛·未诊断","1_pain_undx"),("痛·诊断未治疗","2_pain_dx_untx"),("痛·诊断在治疗","3_pain_dx_tx")]:
    a=c4.get(key,np.nan)
    b=st6[st6["state"]==lab]["ADL受限%"]
    rows.append({"状态":lab,"2011–2018（关节痛）":f"{a:.1f}%" if pd.notna(a) else "—",
                 "2021–2023（体痛）":f"{b.iloc[0]:.1f}%" if len(b) else "—"})
L.append(md(pd.DataFrame(rows)))
L.append("\n> 两期的梯度形状一致（无痛 < 未诊断 < 诊断未治疗 < 诊断在治疗），但 2021–2023 整体高 3–4 个百分点，符合人口老龄化与问卷口径变化的共同作用。\n")

L.append("\n## 表 3　如果强行并入面板 A，它落在哪\n")
panelA=[("HRS 美国",72.8),("SHARE 欧洲",55.5),("ELSA 英国",52.2),("CHARLS 2021–23（体痛）",46.1),("MHAS 墨西哥",34.1)]
L.append("| 库 | 诊断覆盖率（标化） |")
L.append("|---|---|")
for k,v in panelA: L.append(f"| {k} | {v:.1f}% |")
L.append("\n> 中国这一行落在英国与墨西哥之间。**加入它不会改变面板 A 的结论**（诊断覆盖率仍相差 2.1 倍），但它会把「中国」的位置从面板 B 的第 1 名挪到面板 A 的中游——**这正是疼痛口径不同造成的位移，不能当成真实下降**。\n")

L.append("\n## 结论：并进去以后发生了什么\n")
L.append("1. **诊断层的差距没有被改变**。面板 A 从 4 库变 5 库后，最高（美国 72.8）与最低（墨西哥 34.1）仍是 2.1 倍；中国 46.1 落在中间。**「诊断端差距大、治疗端差距小」的结论不受影响。**")
L.append("2. **治疗率出现了值得追的下降，但口径不干净**：诊断且有疼痛者的治疗率从 50.2% 降到 37.0%；用「全部诊断者」这一同分母口径，也从 43.5% 降到 36.5%（降 7 个百分点）。这提示**中国的治疗缺口可能在这十年间扩大了**——但第六轮问的是「治疗方式」、w1–w4 问的是「是否服西药」，**不能当作趋势估计**，只能作为讨论里的提示，或另做一次口径对齐后再下结论。")
L.append("3. **ADL 梯度在两期都保持单调**，2021–2023 整体高 3–4 个百分点；「诊断在治疗」组仍是最差（49.1%），再次证明适应证混杂是结构性的，不是某一波的偶然。")
L.append("4. **对主表编排的建议**：不要把第六轮塞进面板 A/B。主表保持两个面板，第六轮单列为「中国时间对照」（补充表或图 2 的 inset）；正文只写一句：2021–2023 年重复测量显示相同形状，诊断覆盖 46.1%、治疗覆盖 37.0%。")

L.append("\n## 三个必须写进 methods 的差异\n")
L.append("- **疼痛条目**：w1–w4 是 15 部位清单（近一个月），w6 是五级体痛（自上次访问以来）——**分母的构念不同**，不能直接算变化率。")
L.append("- **治疗条目**：w1–w4 是 `rxarthr_c`（是否服任何药物），w6 是 `da021_15_`（治疗方式：是/否）——问法不同。")
L.append("- **年龄**：w6 的年龄由「w5 年龄 + 年差」推算（17,921/18,047 可得），缺失者用出生年补；这与 w1–w4 的直接年龄有细微差别。")

with open(os.path.join(OUT,"第六轮并入对照.md"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("WROTE 第六轮并入对照.md")
