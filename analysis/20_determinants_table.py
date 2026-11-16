# -*- coding: utf-8 -*-
"""20: A2 结果汇总表"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, pandas as pd, numpy as np

det = pd.read_csv(os.path.join(OUT,"19_determinants.csv"))
pool = pd.read_csv(os.path.join(OUT,"19_pooled.csv"))

TERM_LAB = {"age10":"年龄（每 10 岁）","male":"男性（vs 女性）","edu_low":"低教育（vs 高教育）",
            "edu_mid":"中等教育（vs 高教育）","wq":"财富五分位（每升 1 档）",
            "rural":"农村（vs 城市）","comorb":"共病数（每多 1 种）","smoke":"现在吸烟"}
ORDER = ["age10","male","edu_low","edu_mid","wq","rural","comorb","smoke"]

def fmt(r):
    s=f"{r['OR']:.2f}（{r['lo']:.2f}–{r['hi']:.2f}）"
    if r['lo']<=1<=r['hi']: s+="  "
    return s

def tbl(panel, stage):
    g=pool[(pool.panel==panel)&(pool.stage==stage)].set_index("term")
    rows=[]
    for t in ORDER:
        if t not in g.index: continue
        r=g.loc[t]
        rows.append({"变量":TERM_LAB[t],"合并 OR（95% CI）":fmt(r),"k":int(r["k"]),
                     "I²%":r["I2"],"参与库":r["cohorts"]})
    return pd.DataFrame(rows)

def md(df):
    cols=list(df.columns)
    out=["| "+" | ".join(cols)+" |","|"+"|".join(["---"]*len(cols))+"|"]
    for _,r in df.iterrows():
        out.append("| "+" | ".join(str(r[c]) for c in cols)+" |")
    return "\n".join(out)

L=[]
L.append("# 级联流失的决定因素（A2）：谁被漏诊、谁没治疗\n")
L.append("**日期**：2026-11-14　**脚本**：`out/19_determinants.py` → `20_determinants_table.py`\n")
L.append("**模型**：每库 logistic 回归，个体聚类稳健标准误；跨库用 DerSimonian–Laird 随机效应合并 log OR（不做 naïve pooling）。")
L.append("**两个转折点**：① 诊断 = 有疼痛者中是否曾被医生诊断；② 治疗 = 已诊断且有疼痛者中是否在用药。\n")
L.append("**协变量**：年龄（每 10 岁）、性别、教育（低/中 vs 高）、财富五分位（库内分波计算）、城乡（ELSA 无）、共病数（高血压/糖尿病/心脏病/中风/癌症/肺病 6 项）、现在吸烟。\n")

L.append("\n## 表 1　诊断（分母：有疼痛者）\n")
L.append("### 面板 A｜频繁疼痛（HRS、MHAS、ELSA、SHARE）\n")
L.append(md(tbl("A_frequent","D_dx_given_pain")))
L.append("\n### 面板 B｜关节疼痛（CHARLS、KLoSA、LASI、SHARE）\n")
L.append(md(tbl("B_joint","D_dx_given_pain")))

L.append("\n## 表 2　治疗（分母：已诊断且有疼痛者）\n")
L.append("### 面板 A｜仅 HRS、MHAS（ELSA 无该变量、SHARE 只有止痛药）\n")
L.append(md(tbl("A_frequent","T_tx_given_dx_pain")))
L.append("\n### 面板 B｜CHARLS、KLoSA、LASI\n")
L.append(md(tbl("B_joint","T_tx_given_dx_pain")))

L.append("\n## 表 3　关键变量的逐库方向（OR，方便看一致性）\n")
rows=[]
for (p,s),g in det.groupby(["panel","stage"]):
    if s!="D_dx_given_pain": continue
    for t in ["male","edu_low","wq","rural"]:
        sub=g[g.term==t]
        if len(sub)==0: continue
        vals={r.cohort:round(r.OR,2) for r in sub.itertuples()}
        rows.append({"面板":p,"变量":TERM_LAB[t],**{k:vals.get(k,"—") for k in ["HRS","MHAS","ELSA","SHARE","CHARLS","KLoSA","LASI"]}})
t3=pd.DataFrame(rows)
L.append(md(t3))

L.append("\n## 四条结论\n")
L.append("1. **性别是唯一在所有库里都同向的强因素**：男性被诊断的概率只有女性的 **0.54–0.55 倍**，治疗概率是 **0.75–0.83 倍**（7 个模型方向全一致）。这是本文最稳的一个发现。")
L.append("2. **财富梯度不是不存在，而是方向相反、合并后被抵消**：财富五分位每升 1 档，**逐库 OR 全部显著**，但方向分两派——**穷人更易被诊断**：ELSA 0.887、CHARLS 0.922、KLoSA 0.952、SHARE(B) 0.956、MHAS 0.967、HRS 0.968；**富人更易被诊断**：LASI 1.061、SHARE(A) 1.078。合并值 0.97（0.89–1.06）与 0.97（0.92–1.03）之所以跨 1，是**两个相反方向相互抵消**，不是「没有梯度」。")
L.append("> ⚠️ 更值得注意的是 **SHARE 自己就两头占**：A 面板（频繁疼痛，w5–8）OR 1.078、B 面板（关节痛，w1/2/4）OR 0.956。同一队列、不同疼痛定义，财富梯度的方向就翻转——说明该梯度对「谁被算进分母」高度敏感。**因此不能写「财富没有梯度」，应写「财富梯度普遍存在，但方向随卫生系统与测量口径而反转」。**")
L.append("> 合并值接近无效，但逐库方向相反。检验见 `38_wq_linearity.csv`（把五分位改分类变量后，7 库中 5 库的联合 Wald 检验 p<0.05）。")
L.append("3. **年龄与共病是“累积暴露”型因素**：每长 10 岁诊断概率高 20–57%，每多一种共病高 25%——符合“病得越久越可能碰上诊断”的解释，而不是“病人越主动越容易被诊断”。")
L.append("4. **治疗端几乎什么都不是**：治疗层（已诊断且有疼痛）里，除性别外所有变量的合并 OR 都跨 1；点估计方向也不稳定（低教育在 A 面板 0.88、B 面板 1.28）。**一旦确诊，治不治主要不取决于社会人口学特征。**")

L.append("\n## 表 4　敏感性：把疼痛严重度加进诊断模型后，低教育的 OR 变了吗\n")
try:
    sev=pd.read_csv(os.path.join(OUT,"21_severity_sensitivity.csv"))
    sev["变化"]=[f"{(b-a):+.2f}" for a,b in zip(sev["OR_edu_sev"],sev["OR_edu_base"])]
    sv=sev[["panel","cohort","n","sev_var","OR_edu_base","OR_edu_sev","变化"]].copy()
    sv.columns=["面板","库","n","严重度变量","低教育 OR（未加严重度）","低教育 OR（加严重度）","衰减量"]
    L.append(md(sv))
    L.append("\n> 加进疼痛严重度后，低教育的 OR 只小幅下降（HRS 1.19→1.13、ELSA 1.17→1.10、CHARLS 1.57→1.53、KLoSA 2.46→2.26），**方向没有反转**；只有 LASI 从 1.02 掉到 0.90（本来就接近无效）。")
    L.append("> 结论：**「低教育者更容易被诊断」不是由疼痛更重解释的**。剩下的解释只能是报告行为、诊断机会（接触医疗系统的频率）或参照组过小，这三条都要在讨论里写。")
except Exception as e:
    L.append(f"（敏感性表读取失败：{e}）")

L.append("\n## 怎么解释“低教育反而更容易被诊断”\n")
L.append("这与常识相反，但有三种可能，写讨论时必须正面处理：")
L.append("1. **条件于疼痛后的严重度差异**：低教育人群更多从事体力劳动，出现疼痛时病变更重、更“看得见”，因而更可能被医生确诊。若把疼痛严重度放进模型，这个梯度可能消失——**这是下一步必须做的敏感性分析**（各库都有 `painlv`/`painhlv`/`painfrq`）。")
L.append("2. **参照组太小**：CHARLS 高教育只有 2,095 人（2%）、LASI 3,813 人（5%），OR 不稳，这也解释了 I²=87.5%。")
L.append("3. **报告行为**：教育高的人可能更倾向把不适说成“老了”“劳损”，而低教育者更愿意接受“关节炎”这个标签。")
L.append("\n> 因此这条结果**不能写成“低教育者获得更多医疗资源”**。稳妥的写法是：*after conditioning on pain, education was not inversely associated with diagnosis in the expected direction; the association was heterogeneous and should be interpreted cautiously*。")

L.append("\n## 局限（写进 methods）\n")
L.append("- **未加权**：决定因素分析用未加权 logistic（加权会让 CI 失真且不利于跨库比较），因此结果代表“样本内关联”，不外推人群归因。")
L.append("- **财富口径不齐**：CHARLS 用非住房金融财富、HRS 用非住房财富、LASI 用家庭总财富、SHARE 用家庭非住房财富；只能做库内五分位，不能比绝对水平。")
L.append("- **CHARLS 只到 w1–w4**（w5 无权重/无用药），CHARLS 的治疗层样本因此少了 2020 波。")
L.append("- **SHARE 不进治疗层**（`rxpain` 是止痛药，不是抗风湿治疗）；面板 A 的治疗合并只有 2 个库，只能作方向参考。")
L.append("- **横断面**：诊断与疼痛同时测量，「低教育→更重疼痛→更可能诊断」这条通路无法用本设计拆开。")

with open(os.path.join(OUT,"级联流失决定因素.md"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("WROTE 级联流失决定因素.md")
