# -*- coding: utf-8 -*-
"""
44: 疼痛缺失的选择偏倚 —— 加权、年龄标准化后的疼痛率上下界
修正 41 号脚本的口径问题：41 的 painrate 是未加权粗率，与表 2（加权+标化）不可比
（KLoSA 粗率 61.6% vs 标化 51.9%）。
本脚本按 33 号的口径（加权 + WHO 2000–2025 标化，50 岁以上重归一）
重算纳入组疼痛率，并把疼痛缺失者分别全部记为「无痛」与「全部记为有痛」，得到上下界。
输出：44_pain_bounds_weighted.csv / 44_pain_bounds_weighted.txt
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os, numpy as np, pandas as pd, pyreadstat
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))

STD_W = {1: 5370, 2: 4550, 3: 3720, 4: 2970, 5: 2210, 6: 1520, 7: 910, 8: 630}
EDGES = [50, 55, 60, 65, 70, 75, 80, 85, 200]
def agegrp(s):
    return pd.cut(s, bins=EDGES, right=False, labels=list(range(1, 9))).astype(float)

W = pd.read_csv(os.path.join(OUT, "16_weights.csv"), dtype={"ID": str})[["cohort", "ID", "wave", "weight"]]

SPEC = {
    "CHARLS": dict(p=os.path.join(DL, "1. CHARLS  中国", "CHARLS_中国", "Working_data", "charls.dta"),
                   idc="ID", base=["ID", "wave", "age", "arthre"],
                   sites=["da042s2", "da042s3", "da042s4", "da042s5", "da042s8", "da042s9",
                          "da042s10", "da042s11", "da042s12", "da042s13", "da042s14", "da042s15"],
                   pain_mode="list", pain=None, dx="arthre", age="age", waves=[1, 2, 3, 4]),
    "ELSA": dict(p=os.path.join(DL, "2. ELSA 英国", "Working_data", "elsa.dta"),
                 idc="idauniqc", base=["idauniqc", "wave", "agey", "arthre", "painfr"],
                 sites=[], pain_mode="single", pain="painfr", dx="arthre", age="agey", waves=None),
    "HRS": dict(p=os.path.join(DL, "3. HRS  美国", "HRS_美国", "Working_data", "hrs.dta"),
                idc="hhidpn", base=["hhidpn", "wave", "ragey_b", "arthre", "painfr"],
                sites=[], pain_mode="single", pain="painfr", dx="arthre", age="ragey_b", waves=None),
    "KLoSA": dict(p=os.path.join(DL, "4. KLoSA 韩国", "KLoSA_韩国", "Working_data", "klosa.dta"),
                  idc="pid", base=["pid", "wave", "agey", "arthre"],
                  sites=["pain_2", "pain_3", "pain_4", "pain_5", "pain_8",
                         "pain_9", "pain_10", "pain_11", "pain_12", "pain_13"],
                  pain_mode="list", pain=None, dx="arthre", age="agey", waves=None),
    "LASI": dict(p=os.path.join(DL, "5. LASI 印度", "LASI_印度", "Working_data", "lasi.dta"),
                 idc="prim_key", base=["prim_key", "r1agey", "r1arthre", "r1jointp"],
                 sites=[], pain_mode="joint", pain="r1jointp", dx="r1arthre", age="r1agey", waves=None),
    "MHAS": dict(p=os.path.join(DL, "6. MHAS 墨西哥", "MHAS_墨西哥", "Working_data", "mhas.dta"),
                 idc="rahhidnp", base=["rahhidnp", "wave", "agey", "arthre", "painfr"],
                 sites=[], pain_mode="single", pain="painfr", dx="arthre", age="agey", waves=None),
    "SHARE": dict(p=os.path.join(DL, "7.SHARE 欧洲", "SHARE_欧洲", "Working_data", "share.dta"),
                  idc="mergeid", base=["mergeid", "wave", "agey", "arthre", "painfr", "pain_s"],
                  sites=[], pain_mode="single", pain="painfr", dx="arthre", age="agey", waves=None),
}

def load(c, spec, tier):
    cols = list(dict.fromkeys(spec["base"] + spec["sites"]))
    df, _ = pyreadstat.read_dta(spec["p"], usecols=cols)
    if spec["idc"] != "ID":
        df = df.rename(columns={spec["idc"]: "ID"})
    d = df.copy()
    if "wave" not in d.columns:
        d["wave"] = 1
    if spec["waves"] is not None:
        d = d[d["wave"].isin(spec["waves"])]
    if c == "SHARE":
        d = d[d["wave"].isin([1, 2, 4] if tier == "B" else [5, 6, 7, 8])]
    d["ID"] = d["ID"].astype(str)
    d["age"] = pd.to_numeric(d[spec["age"]], errors="coerce")
    d = d[(d["age"] >= 50) & (d["age"] <= 105)].copy()
    d["ageg"] = agegrp(d["age"])
    if tier == "B":
        if c == "LASI":
            pv = pd.to_numeric(d["r1jointp"], errors="coerce")
        elif c == "SHARE":
            pv = pd.to_numeric(d["pain_s"], errors="coerce")
        elif spec["pain_mode"] == "list":
            sub = d[spec["sites"]]
            pv = (sub == 1).any(axis=1).astype(float)
            pv[sub.isna().all(axis=1)] = np.nan
        else:
            pv = pd.to_numeric(d["pain_s"], errors="coerce")
    else:
        if spec["pain_mode"] == "list":
            sub = d[spec["sites"]]
            pv = (sub == 1).any(axis=1).astype(float)
            pv[sub.isna().all(axis=1)] = np.nan
        else:
            pv = pd.to_numeric(d[spec["pain"]], errors="coerce")
    d["pain"] = np.where(pv.notna(), np.where(pv > 0, 1.0, 0.0), np.nan)
    d["dx"] = pd.to_numeric(d[spec["dx"]], errors="coerce").where(lambda s: s.isin([0, 1]))
    d["wave"] = d["wave"].astype(float)
    d = d.merge(W[W.cohort == c][["ID", "wave", "weight"]], on=["ID", "wave"], how="left")
    return d

def std_rate(d, pain_col):
    """加权、按 WHO 2000–2025 标化（50+ 重归一）"""
    ok = d[d["ageg"].notna() & d["weight"].notna() & d[pain_col].notna()]
    if ok.empty:
        return np.nan, np.nan
    num = den = 0.0
    for g, grp in ok.groupby("ageg"):
        w = grp["weight"].sum()
        if w <= 0:
            continue
        r = np.average(grp[pain_col], weights=grp["weight"])
        num += STD_W[int(g)] * r
        den += STD_W[int(g)]
    return (num / den * 100, len(ok)) if den > 0 else (np.nan, np.nan)

rows = []
for c, spec in SPEC.items():
    for tier in ["A", "B"]:
        if tier == "A" and c in ["CHARLS", "KLoSA", "LASI"]:
            continue
        if tier == "B" and c in ["ELSA", "HRS", "MHAS"]:
            continue
        d = load(c, spec, tier)
        d = d[d["weight"].notna()]
        inc = d[d["pain"].notna()].copy()
        exc = d[d["pain"].isna()].copy()
        if len(inc) == 0:
            continue
        obs, n_inc = std_rate(inc, "pain")
        # 下界：缺失者全按无痛；上界：缺失者全按有痛
        lo_df = pd.concat([inc, exc.assign(pain=0.0)], ignore_index=True)
        hi_df = pd.concat([inc, exc.assign(pain=1.0)], ignore_index=True)
        lo, _ = std_rate(lo_df, "pain")
        hi, _ = std_rate(hi_df, "pain")
        rows.append(dict(cohort=c, tier=tier, n_inc=n_inc, n_exc=len(exc),
                         pain_std_obs=round(obs, 1),
                         bound=f"{lo:.1f}–{hi:.1f}",
                         dx_inc=round(inc["dx"].mean() * 100, 1) if inc["dx"].notna().any() else np.nan,
                         dx_exc=round(exc["dx"].mean() * 100, 1) if (len(exc) and exc["dx"].notna().any()) else np.nan))
r = pd.DataFrame(rows)
r["gap"] = (r["dx_inc"] - r["dx_exc"]).round(1)
r = r.sort_values("gap", ascending=False)
r.to_csv(os.path.join(OUT, "44_pain_bounds_weighted.csv"), index=False, encoding="utf-8-sig")

L = []
L.append("44 | 疼痛缺失的偏倚方向与疼痛率上下界（加权 + WHO 2000–2025 标化）")
L.append("口径：与 33_std_ci.csv 的 pain_rate 一致；41 号脚本的 painrate 为未加权粗率，不可与表 2 直接比较")
L.append("")
L.append(r.to_string(index=False))
L.append("")
L.append("bound = 疼痛缺失者全部记为无痛（下界）/ 全部记为有痛（上界）时的标化疼痛率")
L.append("gap = 纳入组的关节炎率 − 疼痛缺失组的关节炎率（百分点）；负值表示缺失者关节炎更多")
open(os.path.join(OUT, "44_pain_bounds_weighted.txt"), "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
