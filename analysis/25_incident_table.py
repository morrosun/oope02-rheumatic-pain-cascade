# -*- coding: utf-8 -*-
"""25: A4 汇总表"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pandas as pd
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
desc=pd.read_csv(os.path.join(OUT,"24_incident_desc.csv"))
res=pd.read_csv(os.path.join(OUT,"24_incident.csv"))
pool=pd.read_csv(os.path.join(OUT,"24_incident_pooled.csv"))
NAME={"CHARLS":"CHARLS 中国","ELSA":"ELSA 英国","HRS":"HRS 美国","KLoSA":"KLoSA 韩国","LASI":"LASI 印度","MHAS":"MHAS 墨西哥","SHARE":"SHARE 欧洲"}
LAB={"age10":"年龄（每 10 岁）","male":"男性（vs 女性）","edu_low":"低教育（vs 高）","edu_mid":"中等教育（vs 高）",
     "wq":"财富五分位（每档）","rural":"农村","comorb":"共病数（每种）","smoke":"现在吸烟",
     "sev_z":"疼痛严重度（z）","adl":"基线 ADL 受限","y_dx":"随访期内获得诊断"}
ORDER=["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke","sev_z","adl","y_dx"]
def fmt(r): return f"{r['OR']:.2f}（{r['lo']:.2f}–{r['hi']:.2f}）"
def md(df):
    cols=list(df.columns); o=["| "+" | ".join(cols)+" |","|"+"|".join(["---"]*len(cols))+"|"]
    for _,r in df.iterrows(): o.append("| "+" | ".join(str(r[c]) for c in cols)+" |")
    return "\n".join(o)

L=[]
L.append("# A4：纵向级联 —— 两年内获得诊断的比例，以及获得诊断后 ADL 的变化\n")
L.append("**日期**：2026-11-14　**脚本**：`out/24_incident.py` → `25_incident_table.py`\n")
L.append("**设计**：把每个库拆成人-波面板，取**相邻两波**（间隔 ≤3 年）。基线样本 = 有疼痛且从未诊断（`pain=1, dx=0`）；结局 = 下一波是否已诊断（`dx=1`）。")
L.append("第二问：以**下一波 ADL** 为结局，调整**基线 ADL** 与疼痛严重度，比较\"获得诊断\"与\"仍未诊断\"。\n")

L.append("\n## 表 1　两年内获得诊断的比例（基线：有疼痛、未诊断）\n")
d=desc.copy(); d["cohort"]=d["cohort"].map(NAME)
d.columns=["库","基线人-波","获得诊断 %","波间隔中位数（年）","可分析 ADL 数"]
d["获得诊断 %"]=d["获得诊断 %"].map(lambda x:f"{x:.1f}")
L.append(md(d))
L.append("\n> KLoSA 只有 2.5%，明显低于其余五库（8.8–18.9%）。原因是韩国 `arthre` 问的是\"风湿病\"、且 w3–w8 之间极少由 0 翻 1；这条要在附录说明，主分析建议报五库、KLoSA 单列。")

L.append("\n## 表 2　谁会获得诊断（合并 OR）\n")
g=pool[pool.model=="M_acquire"].set_index("term")
rows=[]
for t in ORDER:
    if t not in g.index: continue
    r=g.loc[t]; rows.append({"变量":LAB[t],"合并 OR":fmt(r),"k":int(r["k"]),"I²%":r["I2"]})
L.append(md(pd.DataFrame(rows)))

L.append("\n## 表 3　获得诊断与下一波 ADL（调整基线 ADL 与疼痛严重度）\n")
g=pool[pool.model=="M_adl"].set_index("term")
rows=[]
for t in ORDER:
    if t not in g.index: continue
    r=g.loc[t]; rows.append({"变量":LAB[t],"合并 OR":fmt(r),"k":int(r["k"]),"I²%":r["I2"]})
L.append(md(pd.DataFrame(rows)))

L.append("\n## 两条结论\n")
L.append("1. **获得诊断本身不带来功能改善，反而是进展的标记**：在基线功能相同、疼痛严重度相同的前提下，两年内获得诊断者下一波 ADL 受限的 OR 为 **1.73（1.61–1.87）**，且 **I²=0**——六个库方向完全一致。这不能读成\"诊断有害\"，只能读成\"**人被诊断，往往是因为病已经重到去看医生**\"。")
L.append("2. **谁会获得诊断，仍然只跟性别和年龄有关**：男性 OR 0.68（女性更易获得诊断）、年龄 1.12/10 岁、共病 1.07；财富 0.97（每档，虽因样本量大而\"显著\"，但幅度可忽略）、教育、城乡全部跨 1。**这与 A2 的横断面结果完全一致**，说明不是设计造成的假象。")

L.append("\n## 这两问合起来说明什么\n")
L.append("- A3 显示：**未诊断的疼痛者功能负担已经很高**（ADL OR 1.82–3.22），但他们并不比已诊断者更轻；")
L.append("- A4 显示：**拿到诊断并不会让功能变好**，获得诊断反而标志进展；")
L.append("- 结合第 1 步\"诊断者中只有 44–62% 在治疗\"——**级联在\"诊断→治疗\"这一步漏得比\"症状→诊断\"更致命**。")
L.append("\n> 所以论文的主张不能写成\"把诊断率提上去就能改善结局\"。可辩护的写法是：**诊断覆盖率的跨国差异是真实的、可测量的；但在治疗覆盖率只有一半的情况下，单纯提高诊断率不足以改变功能轨迹**。这恰好把\"治疗缺口\"推到前台，也解释了为什么中国的\"高诊断、低治疗\"模式值得单独讨论。")

L.append("\n## 局限\n")
L.append("- **反向因果无法排除**：诊断发生在\"病情重到就医\"的时点，基线 ADL 只能部分校正；`y_dx` 的 OR>1 应读作\"诊断是进展的标记\"，不是\"诊断导致恶化\"。")
L.append("- **随访间隔不齐**：多数库 2 年，MHAS 中位 3 年；已限制 ≤3 年，但未做区间删失处理。")
L.append("- **CHARLS 用到 w1–w5**（w5 无用药但 `arthre` 可用），未接入第六轮 `da011_15_`——后者是「自上次访问以来的新诊断」，可作为独立验证（样本 18,047，见 `第六轮级联小表.md`）。")
L.append("- **LASI 单波**，无法进入 A4。")
L.append("- **ADL 口径**与 A3 相同，KLoSA/SHARE 由条目合成。")

with open(os.path.join(OUT,"纵向级联.md"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("WROTE 纵向级联.md")
