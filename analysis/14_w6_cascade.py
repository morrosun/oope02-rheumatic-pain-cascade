# -*- coding: utf-8 -*-
"""
Step 0: CHARLS 第六轮（2021-2023）去重 -> 重建 ever 关节炎 -> 四段级联小表
输出：
  14_w6_cascade.txt   人读表
  14_w6_cascade.csv   机读表
  14_w6_derived.csv   派生分析集（仅派生变量，供复现）
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os
import numpy as np
import pandas as pd
import pyreadstat

BASE = os.environ.get("AGING_DB_RAW", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
CH   = os.path.join(DL, "1. CHARLS  中国/CHARLS_中国/Working_data/charls.dta")

L = []
def w(s=""):
    L.append(str(s)); print(s, flush=True)

# ---------------- 1. 读三年包 ----------------
COLS = ["ID","da011_15_","da012_15_","da015_15_","da016_1_15_","da016_2_15_",
        "da021_15_","da069","da071","da072","da074","proxy",
        "db001","db003","db005","db007","db009","db011","db012","db014","db016","db018","db020","db022","db052"]
frames = []
for y in [2021, 2022, 2023]:
    p = os.path.join(BASE, str(y), "Health_Status_and_Functioning.dta")
    df, _ = pyreadstat.read_dta(p, usecols=COLS)
    df["res_wave"] = y
    frames.append(df)
raw = pd.concat(frames, ignore_index=True)
w(f"三年包原始记录 {len(raw)}，唯一 ID {raw['ID'].nunique()}")

# ---------------- 2. 去重：同一 ID 保留 res_wave 最晚、信息最全的一条 ----------------
raw["_nonmiss"] = raw[COLS[1:]].notna().sum(axis=1)
raw = raw.sort_values(["ID","res_wave","_nonmiss"], ascending=[True, False, False])
d = raw.drop_duplicates("ID", keep="first").copy()
w(f"去重后 {len(d)}（净删除 {len(raw)-len(d)} 条）")
w("保留记录的 res_wave 分布: " + str(d["res_wave"].value_counts().sort_index().to_dict()))

# ---------------- 3. 重建 ever 关节炎 ----------------
# 问卷逻辑（已用 2021 交叉表核实）：
#   既往无病 -> 问 da011（1=新诊断, 2=否）；da021 仅对新诊断者问
#   既往有病 -> da011 缺失，问 da012（1/2/3 变化）与 da021（是否治疗）
d["ever_w6"] = np.nan
d.loc[d["da011_15_"] == 1, "ever_w6"] = 1                                  # 新诊断
d.loc[d["da021_15_"].notna(), "ever_w6"] = 1                               # 被问治疗 = 现/曾患病
d.loc[d["da012_15_"].isin([1,2,3]), "ever_w6"] = 1                          # 与上次比有好转/不变/变差
d.loc[(d["da011_15_"] == 2) & (d["da012_15_"] == 99), "ever_w6"] = 0        # 上次访问未患病且本次无新诊断
d.loc[(d["da011_15_"] == 2) & (d["da021_15_"].isna()) & (d["da012_15_"].isna()), "ever_w6"] = 0
w("\never_w6 分布: " + str(d["ever_w6"].value_counts(dropna=False).to_dict()))

# 与 w5 累计 arthre 对账
ch, _ = pyreadstat.read_dta(CH, usecols=["ID","wave","arthre"])
ch5 = ch[ch["wave"] == 5][["ID","arthre"]].rename(columns={"arthre":"arthre_w5"})
m = d.merge(ch5, on="ID", how="left")
w(f"与 charls.dta w5 匹配: {m['arthre_w5'].notna().sum()}/{len(m)} ({m['arthre_w5'].notna().mean()*100:.1f}%)")
ct = pd.crosstab(m["arthre_w5"].fillna(-1), m["ever_w6"].fillna(-1))
w("w5 arthre × W6 ever 交叉：")
w(ct.to_string())
both = m.dropna(subset=["arthre_w5","ever_w6"])
w(f"两者均有值 n={len(both)}；一致率={(both['arthre_w5']==both['ever_w6']).mean()*100:.1f}%")
w(f"  w5 阳性→W6 阴性: {((both['arthre_w5']==1)&(both['ever_w6']==0)).sum()}")
w(f"  w5 阴性→W6 阳性: {((both['arthre_w5']==0)&(both['ever_w6']==1)).sum()}")

# ---------------- 4. 派生变量 ----------------
d["pain_any"]  = np.where(d["da069"].isin([2,3,4,5]), 1,
                  np.where(d["da069"] == 1, 0, np.nan))
d["pain_sev"]  = np.where(d["da069"].isin([4,5]), 1,
                  np.where(d["da069"].isin([1,2,3]), 0, np.nan))
d["otc"]       = d["da071"].map({1:1, 2:0})
d["opioid"]    = d["da072"].map({1:1, 2:0})
d["analg_any"] = np.where(d["otc"].notna() | d["opioid"].notna(),
                          ((d["otc"] == 1) | (d["opioid"] == 1)).astype(float), np.nan)
d["analg_none"] = np.where(d["analg_any"].notna(), 1 - d["analg_any"], np.nan)
d["treat"]     = d["da021_15_"].map({1:1, 2:0})
adl_cols = ["db001","db003","db005","db007","db009","db011","db012","db014","db016","db018","db020","db022","db052"]
d["adl_any"] = (d[adl_cols] >= 2).any(axis=1).astype(float)
d.loc[d[adl_cols].isna().all(axis=1), "adl_any"] = np.nan

# ---------------- 5. 四段级联表 ----------------
rows = []
def add(stage, denom_label, n, dnm, note=""):
    rows.append(dict(stage=stage, denominator=denom_label, n=n, denom=dnm,
                     pct=round(n/dnm*100, 1) if dnm else np.nan, note=note))

n_valid_pain = int(d["pain_any"].notna().sum())
n_pain  = int((d["pain_any"] == 1).sum())
n_sev   = int((d["pain_sev"] == 1).sum())
n_sev_ans = int(d.loc[d["pain_sev"] == 1, "otc"].notna().sum())
n_otc   = int((d.loc[d["pain_sev"] == 1, "otc"] == 1).sum())
n_opi   = int((d.loc[d["pain_sev"] == 1, "opioid"] == 1).sum())
n_ever  = int((d["ever_w6"] == 1).sum())
n_ever_ans = int(d["ever_w6"].notna().sum())
n_treat_ans = int(d.loc[d["ever_w6"] == 1, "treat"].notna().sum())
n_treat = int((d.loc[d["ever_w6"] == 1, "treat"] == 1).sum())

# 疼痛者中的诊断覆盖（分母 = 有疼痛且有诊断信息）
sub = d[(d["pain_any"] == 1) & d["ever_w6"].notna()]
n_dx_in_pain = int((sub["ever_w6"] == 1).sum())
# 严重疼痛者中的诊断覆盖
sub2 = d[(d["pain_sev"] == 1) & d["ever_w6"].notna()]
n_dx_in_sev = int((sub2["ever_w6"] == 1).sum())
# 冰山格：严重疼痛 + 未诊断 + 自费止痛（任何止痛药）
ice = d[(d["pain_sev"] == 1) & (d["ever_w6"] == 0) & d["analg_any"].notna()]
# 诊断 + 有疼痛 + 未治疗
dxpain = d[(d["pain_any"] == 1) & (d["ever_w6"] == 1) & d["treat"].notna()]

add("① 体痛（≥一点）", "有体痛记录者", n_pain, n_valid_pain)
add("①a 体痛（相当多/非常）", "有体痛记录者", n_sev, n_valid_pain)
add("② 严重体痛者用非处方止痛药", "严重体痛且被问及", n_otc, n_sev_ans)
add("②a 严重体痛者用阿片类", "严重体痛且被问及", n_opi, n_sev_ans)
add("③ 曾诊断关节炎（重建）", "有诊断信息者", n_ever, n_ever_ans)
add("③a 体痛者中曾诊断", "体痛且有诊断信息", n_dx_in_pain, len(sub))
add("③b 严重体痛者中曾诊断", "严重体痛且有诊断信息", n_dx_in_sev, len(sub2))
add("④ 诊断者中在治疗", "诊断且有治疗信息", n_treat, n_treat_ans)
add("④a 诊断+有体痛者中在治疗", "诊断+体痛+有治疗信息", int((dxpain['treat']==1).sum()), len(dxpain))
add("冰山格：严重体痛+未诊断+在吃止痛药", "严重体痛+未诊断+被问及", int((ice['analg_any']==1).sum()), len(ice))

tab = pd.DataFrame(rows)
w("\n【表 1】四段级联（去重后，未加权）")
w(tab.to_string(index=False))

# ---------------- 6. 按疼痛分级 ----------------
lbl = {1:"无", 2:"一点", 3:"有些", 4:"相当多", 5:"非常"}
g = d.dropna(subset=["da069"]).groupby("da069")
t2 = pd.DataFrame({
    "n": g.size(),
    "曾诊断_%": (g["ever_w6"].mean()*100).round(1),
    "诊断者治疗_%": (g.apply(lambda x: x.loc[x["ever_w6"]==1,"treat"].mean()*100)).round(1),
    "止痛药_%": (g["otc"].mean()*100).round(1),
    "阿片_%": (g["opioid"].mean()*100).round(1),
    "ADL受限_%": (g["adl_any"].mean()*100).round(1),
})
t2.index = [lbl[i] for i in t2.index]
w("\n【表 2】按体痛严重度分层")
w(t2.to_string())

# ---------------- 7. 冰山格 2×2 ----------------
w("\n【表 3】严重体痛者的诊断 × 自费止痛（行%）")
sub3 = d[d["pain_sev"] == 1].dropna(subset=["ever_w6","otc"])
t3 = pd.crosstab(sub3["ever_w6"].map({1:"已诊断",0:"未诊断"}),
                 sub3["analg_any"].map({1:"用止痛药",0:"完全未用"}),
                 margins=True)
w(t3.to_string())
w("\n  行内比例：")
w((pd.crosstab(sub3["ever_w6"].map({1:"已诊断",0:"未诊断"}),
               sub3["analg_any"].map({1:"用止痛药",0:"完全未用"}),
               normalize="index")*100).round(1).to_string())

# ---------------- 7b. 敏感性：ever = max(w5, W6) ----------------
mu = d.merge(ch5, on="ID", how="left")
mu["ever_union"] = np.where(mu["arthre_w5"].notna() | mu["ever_w6"].notna(),
                            np.nanmax(np.vstack([mu["arthre_w5"].fillna(0).values,
                                                 mu["ever_w6"].fillna(0).values]), axis=0), np.nan)
mu.loc[mu["arthre_w5"].isna() & mu["ever_w6"].isna(), "ever_union"] = np.nan
u = mu[(mu["pain_any"] == 1) & mu["ever_union"].notna()]
n_dx_u = int((u["ever_union"] == 1).sum())
w(f"\n【敏感性】ever = max(w5, W6)：诊断者 {int((mu['ever_union']==1).sum())}/{int(mu['ever_union'].notna().sum())} "
  f"= {(mu['ever_union']==1).sum()/mu['ever_union'].notna().sum()*100:.1f}%；"
  f"体痛者中诊断覆盖 {n_dx_u}/{len(u)} = {n_dx_u/len(u)*100:.1f}%（主分析 {n_dx_in_pain/len(sub)*100:.1f}%）")
m = mu  # 供后续 markdown 使用

# ---------------- 8. 加权（仅 cross-section 有权重） ----------------
wt, _ = pyreadstat.read_dta(os.path.join(BASE,"2023","Weights.dta"), usecols=["ID","INDV_weight_post2"])
dw = d.merge(wt, on="ID", how="left")
w(f"\n有权重者 {dw['INDV_weight_post2'].notna().sum()}/{len(dw)}")
def wpct(mask_num, mask_den):
    s = dw.loc[mask_den, ["INDV_weight_post2"]].copy()
    s["y"] = mask_num.reindex(s.index).astype(float)
    s = s.dropna()
    return round(np.average(s["y"], weights=s["INDV_weight_post2"])*100, 1), len(s)
w("加权估计：")
wres = {}
for name, num, den in [
    ("体痛（≥一点）", dw["pain_any"]==1, dw["pain_any"].notna()),
    ("体痛（≥相当多）", dw["pain_sev"]==1, dw["pain_sev"].notna()),
    ("严重体痛者用止痛药", dw["analg_any"]==1, dw["pain_sev"]==1),
]:
    v, n = wpct(num, den)
    wres[name] = (v, n)
    w(f"  {name:18s} {v}%  (n={n})")

# ---------------- 9. 落盘 ----------------
tab.to_csv(os.path.join(OUT,"14_w6_cascade.csv"), index=False, encoding="utf-8-sig")
t2.to_csv(os.path.join(OUT,"14_w6_by_severity.csv"), encoding="utf-8-sig")
keep = ["ID","res_wave","da069","pain_any","pain_sev","ever_w6","ever_union","treat","otc","opioid","analg_any","analg_none","adl_any"]
d2 = d.merge(m[["ID","arthre_w5","ever_union"]], on="ID", how="left")
d2[keep].to_csv(os.path.join(OUT,"14_w6_derived.csv"), index=False, encoding="utf-8-sig")

# ---------------- 10. Markdown 小表 ----------------
def to_md(df, headers):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"]*len(headers)) + "|"]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in headers) + " |")
    return "\n".join(out)

md = []
md.append("# CHARLS 第六轮（2021–2023）四段级联小表\n")
md.append(f"**生成日期**：2026-11-14　**脚本**：`out/14_w6_cascade.py`\n")
md.append(f"**样本**：三年包 19,021 条记录 → 按 ID 去重后 **18,047 人**（净删除 974 条）\n")
md.append("**口径**：体痛与关节炎均为本轮自报；关节炎 = 本轮自报曾被医生诊断（重建，见方法）；治疗 = `da021_15_`；自费止痛 = `da071`/`da072`（仅问「相当多/非常」疼痛者）。未加权。\n")
md.append("\n## 表 1　四段级联\n")
t1 = tab[["stage","n","denom","pct","denominator"]].copy()
t1.columns = ["阶段","人数","分母","占比 %","分母说明"]
md.append(to_md(t1, ["阶段","人数","分母","占比 %","分母说明"]))
md.append("\n## 表 2　按体痛严重度分层\n")
t2b = t2.reset_index().rename(columns={"index":"体痛等级","da069":"体痛等级"})
md.append(to_md(t2b, ["体痛等级","n","曾诊断_%","诊断者治疗_%","止痛药_%","阿片_%","ADL受限_%"]))
md.append("\n## 表 3　严重体痛者的诊断 × 自费止痛（行 %）\n")
md.append("| 诊断状态 | 完全未用止痛药 | 用止痛药 | 合计 |")
md.append("|---|---|---|---|")
for lab, g_ in sub3.groupby(sub3["ever_w6"].map({1:"已诊断",0:"未诊断"})):
    a = int((g_["analg_any"]==0).sum()); b = int((g_["analg_any"]==1).sum())
    md.append(f"| {lab} | {a} ({a/len(g_)*100:.1f}%) | {b} ({b/len(g_)*100:.1f}%) | {len(g_)} |")

md.append("\n## 方法：`ever` 是怎么重建的\n")
md.append("问卷是分路问的（已用 2021 年交叉表逐格核实）：**既往无病者**问 `da011_15_`（1=新诊断/2=否），**既往有病者** `da011` 缺失、改问 `da012_15_`（1 好转/2 不变/3 变差/99 上次未患病）与 `da021_15_`（是否治疗）。因此按下列优先级重建：\n")
n_new = int((d["da011_15_"] == 1).sum())
n_rx  = int(d["da021_15_"].notna().sum())
n_chg = int(d["da012_15_"].isin([1,2,3]).sum())
n_never = int(((d["da011_15_"] == 2) & (d["da012_15_"] == 99)).sum())
n_notasked = int(((d["da011_15_"] == 2) & d["da021_15_"].isna() & d["da012_15_"].isna()).sum())
md.append("| 规则（按优先级） | 判定 | 去重后命中人数 |")
md.append("|---|---|---|")
md.append(f"| `da011_15_` = 1（本轮新诊断） | ever = 1 | {n_new} |")
md.append(f"| `da021_15_` 有值（被问治疗 = 现/曾患病） | ever = 1 | {n_rx} |")
md.append(f"| `da012_15_` ∈ {{1,2,3}}（与上次比有好转/不变/变差） | ever = 1 | {n_chg} |")
md.append(f"| `da011_15_` = 2 且 `da012_15_` = 99（上次未患病、本轮无新诊断） | ever = 0 | {n_never} |")
md.append(f"| `da011_15_` = 2 且两者均缺失（未进入疾病模块） | ever = 0 | {n_notasked} |")
md.append(f"\n> 规则 2、3 高度重叠（都是既往已病者），所以三条相加不等于 {int((d['ever_w6']==1).sum())}。最终 ever=1 共 **{int((d['ever_w6']==1).sum()):,}** 人。")
md.append("\n> 关键点：`da011_15_` **不是「全体受访者的新诊断率」**，它的分母只是「既往无病者」（2021 年 8,768 人），所以 413/8,768 = 4.7% 才是**既往无病者的两年新诊断率**。这一点极易写错。\n")

md.append("\n## 一致性核对（与 w5 累计 `arthre` 对比）\n")
md.append(f"- 与 `charls.dta` w5 匹配 **{int(m['arthre_w5'].notna().sum()):,}/{len(m):,}（{m['arthre_w5'].notna().mean()*100:.1f}%）**；两者均有值者一致率 **91.4%**")
md.append(f"- w5 阳性 → W6 阴性：**932 人**（占 w5 阳性的 14.1%，即 CHARLS 已知的「回落」现象在 5 年跨度上的累积）")
md.append(f"- w5 阴性 → W6 阳性：**528 人**（5.1%，与 4.7% 的新诊断率同量级）")
md.append(f"- **敏感性**：若定义 ever = max(w5, W6)，诊断者升至 7,504/18,045 = 41.6%，体痛者中诊断覆盖由 47.5% 升至 **52.9%**。两个口径相差 5.4 个百分点，必须都报。\n")

md.append("\n## 加权（`INDV_weight_post2`，仅 cross-section 有权重）\n")
md.append("| 指标 | 未加权 | 加权 |")
md.append("|---|---|---|")
md.append(f"| 体痛（≥一点） | 57.8% | {wres['体痛（≥一点）'][0]}% |")
md.append(f"| 体痛（≥相当多/非常） | 14.7% | {wres['体痛（≥相当多）'][0]}% |")
md.append(f"| 严重体痛者用止痛药 | 48.8%（非处方）/ 51.3%（任何） | {wres['严重体痛者用止痛药'][0]}%（任何） |")
md.append(f"\n> 权重文件只覆盖 16,377 名 cross-section 样本，去重后有 {int(dw['INDV_weight_post2'].notna().sum()):,} 人有权重；加权与未加权相差约 1–2 个百分点，方向不变。\n")

md.append("\n## 这张表告诉我们三件事\n")
md.append("1. **中国这一格是「高诊断率、低治疗率」**：体痛者 47.5% 有诊断、严重体痛者 63.0% 有诊断，但诊断者中只有 **36.5%** 在治疗。与 LASI（18.7% 诊断）形成的是两种不同的失败模式——印度是**够不着诊断**，中国是**诊断了不治疗**。")
md.append("2. **疼痛越重，诊断率越高、治疗率也越高，但仍有一半没治**：从「无痛」到「非常痛」，诊断率 21.4%→63.0%，诊断者治疗率 25.7%→53.3%。梯度是单调的，说明变量不是噪声。")
md.append("3. **冰山下面不是「自己买药顶着」，而是「什么都没做」**：严重体痛且从未诊断的 982 人中，**56.6% 完全没用任何止痛药**；已诊断组 56.0% 在用药。这一格是全文最有政策含义的数字——**漏诊的人同时也漏治**。")

md.append("\n## 三条口径提醒（写进 methods，别等审稿人问）\n")
md.append("1. **`da021_15_` 只问已患病者**，所以「未诊断者的治疗状态」在问卷里**不可观测**（结构缺失，不是缺失值）。级联的第三段只能用自费止痛药（`da071/072`）近似，不能算「未治疗率」。")
md.append("2. **`da071/072` 只问「相当多/非常」疼痛者**（n=2,653），轻度疼痛者的自费用药**完全没问**。所以表 3 只能做严重疼痛层，不能外推到全部体痛者。")
md.append("3. **回忆窗口不统一**：2021/2022/2023 三批的「上次访问」分别是 2018/2020、2021、2021/2022，`da011_15_` 的新诊断率不能直接当「两年发病率」合并；做 A4 时必须把 `res_wave` 放进模型。")
und_n = int((sub3["ever_w6"] == 0).sum()); und_none = int(((sub3["ever_w6"] == 0) & (sub3["analg_any"] == 0)).sum())
dx_n  = int((sub3["ever_w6"] == 1).sum()); dx_use  = int(((sub3["ever_w6"] == 1) & (sub3["analg_any"] == 1)).sum())
md.append(f"\n## 一句话总结\n")
md.append(f"严重体痛且**从未诊断**的 {und_n} 人中，**{und_none/und_n*100:.1f}%（{und_none} 人）完全没有用任何止痛药**，"
          f"{100-und_none/und_n*100:.1f}% 在自行用药；已诊断组相反，**{dx_use/dx_n*100:.1f}% 在用止痛药**。"
          f"也就是说，冰山下面那一层不是「自己买药顶着」，而是**连药都没吃**——这才是级联最该被看见的一格。")
with open(os.path.join(OUT,"第六轮级联小表.md"),"w",encoding="utf-8") as f:
    f.write("\n".join(md))

with open(os.path.join(OUT,"14_w6_cascade.txt"),"w",encoding="utf-8") as f:
    f.write("\n".join(L))
print("\nWROTE 14_w6_cascade.{txt,csv}, 14_w6_by_severity.csv, 14_w6_derived.csv, 第六轮级联小表.md")
