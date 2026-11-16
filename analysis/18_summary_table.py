# -*- coding: utf-8 -*-
"""18: 把 17 的 CSV 组装成可读汇总表"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pandas as pd, numpy as np

A = pd.read_csv(os.path.join(OUT,"17_cascade_tierA.csv"))
B = pd.read_csv(os.path.join(OUT,"17_cascade_tierB.csv"))
F = pd.read_csv(os.path.join(OUT,"17_four_state.csv"))
M = pd.read_csv(os.path.join(OUT,"17_missing.csv"))
W6 = pd.read_csv(os.path.join(OUT,"14_w6_cascade.csv"))

NAME = {"CHARLS":"CHARLS 中国","ELSA":"ELSA 英国","HRS":"HRS 美国","KLoSA":"KLoSA 韩国",
        "LASI":"LASI 印度","MHAS":"MHAS 墨西哥","SHARE":"SHARE 欧洲"}
ORDER_A = ["HRS","SHARE","ELSA","MHAS"]
ORDER_B = ["CHARLS","SHARE","KLoSA","LASI"]

def md_table(df, cols, headers):
    out=["| "+" | ".join(headers)+" |", "|"+"|".join(["---"]*len(headers))+"|"]
    for _,r in df.iterrows():
        vals=[]
        for c in cols:
            v=r[c]
            if isinstance(v,float) and np.isnan(v): vals.append("—")
            else: vals.append(str(v))
        out.append("| "+" | ".join(vals)+" |")
    return "\n".join(out)

L=[]
L.append("# 七库风湿疼痛诊疗级联：加权 + 年龄标准化汇总表\n")
L.append("**日期**：2026-11-14　**脚本**：`out/16_build_weights.py` → `17_cascade_weighted.py` → `18_summary_table.py`\n")
L.append("**样本**：7 库、覆盖 20 余个国家/地区；年龄限制 50–105 岁；只保留有抽样权重的人-波。\n")
L.append("**权重**：CHARLS `r{w}wtresp`（w1–w4）、ELSA `r{w}cwtresp`、HRS `wtresp`、KLoSA `r{w}wtresp`、LASI `r1wtresp`、MHAS `r{w}wtresp`、SHARE `r{w}wtresp`。\n")
L.append("**年龄标准化**：WHO World Standard Population 2000–2025，50 岁以上组重归一（50–54…85+ 权重 5370/4550/3720/2970/2210/1520/910/630）。\n")

L.append("\n## 为什么分成两个面板\n")
L.append("七库的疼痛问法不是一回事，硬拼会造出假梯度：\n")
L.append("| 面板 | 疼痛条目 | 覆盖库 | 含义 |")
L.append("|---|---|---|---|")
L.append("| **A 频繁疼痛** | \"是否经常被疼痛困扰\"（`painfr`） | HRS、MHAS、ELSA、SHARE（w5–8） | 任何部位的频繁疼痛 |")
L.append("| **B 关节疼痛** | 部位清单（肩/腕/手指/背/髋/膝/踝/颈…）或\"关节痛或僵硬\" | CHARLS、KLoSA、LASI、SHARE（w1/2/4） | 肌骨/关节源性疼痛 |")
L.append("\n> SHARE 同时进两个面板，但两种问法**分处不同波次**（`pain_s` 在 w1/2/4，`painfr` 在 w5–8），**从未在同波共存**。所以 63.2%（w5–8，频繁疼痛）与 42.9%（w1/2/4，关节痛）的差异同时包含**问法**与**年代**两个成分，不能全部归因于问法。\n>\n> 但这反而更有力地证明了必须拆面板：**同一个队列在同一份研究里给出两个相差 16–20 个百分点的诊断覆盖率**（标化后 55.5% vs 39.3%，差 16.2 点）。无论差额来自问法还是年代，都说明「关节炎诊断覆盖率」这个数字离开疼痛条目定义就不成立——任何把七库硬拼成一张表的做法，都无法解释 SHARE 自己这两个数。\n")

L.append("\n## 表 A　频繁疼痛面板（加权）\n")
ta = A.set_index("cohort").loc[ORDER_A].reset_index()
ta["cohort"]=ta["cohort"].map(NAME)
L.append(md_table(ta, ["cohort","n_all","n_weighted","n_pain","pain_crude","pain_std","dx_crude","dx_std","txpain_crude","txpain_std"],
                  ["库","人波","有权重","有疼痛","疼痛% 粗","疼痛% 标化","诊断覆盖% 粗","诊断覆盖% 标化","诊断+疼痛者治疗% 粗","治疗% 标化"]))
L.append("\n> ELSA 无关节炎用药变量；SHARE 只有\"每周止痛药\"（`rxpain`），其\"治疗\"列不能与其余库直接比。\n")

L.append("\n## 表 B　关节疼痛面板（加权）\n")
tb = B.set_index("cohort").loc[ORDER_B].reset_index()
tb["cohort"]=tb["cohort"].map(NAME)
L.append(md_table(tb, ["cohort","n_all","n_weighted","n_pain","pain_crude","pain_std","dx_crude","dx_std","txpain_crude","txpain_std"],
                  ["库","人波","有权重","有关节痛","疼痛% 粗","疼痛% 标化","诊断覆盖% 粗","诊断覆盖% 标化","诊断+疼痛者治疗% 粗","治疗% 标化"]))

L.append("\n## 表 C　四状态 × ADL 受限（加权，未标准化）\n")
piv = F.pivot(index="cohort", columns="state", values="adl_crude")
cnt = F.pivot(index="cohort", columns="state", values="n")
rows=[]
for c in ["CHARLS","KLoSA","LASI","MHAS","HRS","SHARE","ELSA"]:
    if c not in piv.index: continue
    r={"库":NAME[c]}
    for st,lab in [("0_no_pain","无痛"),("1_pain_undx","痛·未诊断"),("2_pain_dx_untx","痛·诊断未治疗"),
                   ("2_pain_dx","痛·已诊断"),("3_pain_dx_tx","痛·诊断在治疗")]:
        if st in piv.columns and pd.notna(piv.loc[c,st]):
            r[lab]=f"{piv.loc[c,st]:.1f}%（{int(cnt.loc[c,st]):,}）"
        else:
            r[lab]="—"
    rows.append(r)
t3=pd.DataFrame(rows).fillna("—")
L.append(md_table(t3, list(t3.columns), list(t3.columns)))

L.append("\n> **这张表不能当因果读**。四个库里 ADL 受限都是「无痛 < 未诊断 < 诊断未治疗 < 诊断在治疗」单调上升——治疗组最差，是典型的**适应证混杂**（病重的人才吃药）。A3 只能比较\"未诊断 vs 已诊断\"并调整疼痛严重度，或干脆不做治疗组间比较。\n")

L.append("\n## 表 D　三变量联合缺失（人波）\n")
mm = M.groupby("cohort").apply(lambda g: pd.Series({
    "总人波": f"{int(g['n'].sum()):,}",
    "疼痛缺失%": round(g.loc[g["p"],"n"].sum()/g["n"].sum()*100,2),
    "诊断缺失%": round(g.loc[g["x"],"n"].sum()/g["n"].sum()*100,2),
    "治疗缺失%": round(g.loc[g["t"],"n"].sum()/g["n"].sum()*100,2),
}), include_groups=False).reset_index()
mm["cohort"]=mm["cohort"].map(NAME)
L.append(md_table(mm, ["cohort","总人波","疼痛缺失%","诊断缺失%","治疗缺失%"],
                  ["库","总人波","疼痛缺失%","诊断缺失%","治疗缺失%"]))
L.append("\n> ⚠️ **本表的「治疗缺失%」是按全部人波算的，不能用来判断治疗率的分母是否干净。** 治疗率的分母是「已被诊断者」，按此重算：")
L.append(">")
L.append("> | 库 | 诊断者 n | 有治疗值 | **诊断者中缺失** | 逐波 |")
L.append("> |---|---:|---:|---:|---|")
L.append("> | CHARLS | 31,096 | 22,334 | **28.2%** | w1 0.4%、w2 0.9%、**w3 22.4%**、w4 0%、w5 100% |")
L.append("> | HRS | 115,948 | 73,944 | **36.2%** | w5–11 均 <1%、**w12–15 全部 100%** |")
L.append("> | KLoSA | 10,801 | 10,797 | 0.0% | 各波 ≤0.1% |")
L.append("> | LASI | 6,708 | 6,618 | 1.3% | — |")
L.append("> | MHAS | 11,696 | 11,675 | 0.2% | — |")
L.append("> | SHARE | 105,848 | 105,569 | 0.3% | — |")
L.append(">")
L.append("> **两个直接后果**：① HRS 的 42.8% 治疗率只代表 **2000–2012（w5–11）**，不代表 2014–2020；② CHARLS 的 43.5% 受 w3 那 22.4% 缺失影响。两句话都必须进表注，否则读者会以为覆盖全时段。治疗缺失是结构性的（w12–15 没问、ELSA 无此变量、SHARE 只有止痛药），**不是数据质量差**。\n")

L.append("\n## 四条可以直接写进论文的结果\n")
L.append("1. **诊断层是\"宽带\"，治疗层是\"窄带\"**——但只能在面板内比：")
L.append("   - B 面板（关节痛）：诊断覆盖标化后 **21.1%（印）→ 25.2%（韩）→ 39.3%（欧）→ 55.1%（中）**，相差 2.6 倍；而诊断且有疼痛者的治疗率只有 **50.2–61.7%**，相差 1.2 倍。")
L.append("   - A 面板（频繁疼痛）：诊断覆盖 **34.1%（墨）→ 52.2%（英）→ 55.5%（欧）→ 72.8%（美）**，相差 2.1 倍；治疗率 **44.0–55.2%**。")
L.append("   - 两个面板给出同一个结论：**跨国差异主要在\"有没有被诊断\"，不在\"诊断后治不治\"**。")
L.append("2. **中国的定位只能在 B 面板内说**：关节痛者诊断覆盖 55.1%（四库最高），但诊断且有疼痛者的治疗率 50.2%，低于韩国 61.7%、欧洲 55.6%（止痛药口径）、印度 53.3%。用一句话概括是\"**诊断不差、治疗偏弱**\"。**不能拿中国的 55.1% 去比美国的 72.8%**——两个数来自不同疼痛条目。")
L.append("3. **疼痛定义本身就能改写结论**：SHARE 在「频繁疼痛」口径（w5–8）下诊断覆盖 63.2%，在「关节痛」口径（w1/2/4）下只有 42.9%，差约 20 个百分点（标化后 55.5% vs 39.3%，差 16.2 点）。⚠️ **两种问法从未在同一波共存**，差额同时含问法与年代两个成分，不能写成「纯粹来自问法」。但即便如此，同一队列给出两个相差 16–20 点的数字，已足以证明必须分两个面板而不是硬拼成一张七库表。")
L.append("4. **ADL 梯度在七个库内部都单调（加权粗率）**：无痛·未诊断 → 痛·未诊断 → 痛·诊断未治疗 → 痛·诊断在治疗，如 HRS **3.9→18.4→31.0→42.2%**、CHARLS 9.1→24.2→33.0→45.0%。这不是治疗有害，是适应证混杂。\n> ⚠️ 本表 C 的「无痛」组混着已诊断者（HRS 合并值 7.0%，非 3.9%）；拆分版本见 `37_four_state_v2.csv`。另注：调整协变量后 KLoSA 的「诊断未治疗」OR 为 0.96（<「未诊断」），故「单调」只成立于粗率，不成立于调整后 OR。")

L.append("\n## 口径提醒（写进 methods）\n")
L.append("1. **CHARLS 只用了 w1–w4**：w5 的 harmonized 层没有 `r5wtresp`（无权重），且 w5 无 `rxarthr`，无法进入加权治疗层；w5 的疼痛×诊断可在附录单独报。第六轮（2021–2023）见 `第六轮级联小表.md`，其疼痛是五级体痛、口径不同，不进本表。")
L.append("2. **CHARLS 治疗用 `rxarthr_c`（任何治疗，含中成药）**，不是 `rxarthr`（仅西药）。用西药口径时诊断者治疗率会从 43.5% 降到约 36%。")
L.append("3. **SHARE 的 `rxpain` 是\"每周至少一次止痛药\"**，不是抗风湿治疗，故四状态表里 SHARE 只分三段。")
L.append("4. **LASI 是单波横断面**，不能进纵向分析；其\"有关节痛\"是\"关节痛或僵硬\"，比其他库宽。")
L.append("5. **年龄标准化只解决年龄构成**，不解决抽样设计（分层/整群）。本表已加抽样权重，但未做设计效应校正，标准误需用 bootstrap 或 Taylor 线性化重算。")

with open(os.path.join(OUT,"第七库级联汇总表.md"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("WROTE 第七库级联汇总表.md")
print("\n".join(L[:40]))
