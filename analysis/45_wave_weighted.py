# -*- coding: utf-8 -*-
"""
45: 按波次的加权级联率（S3 的加权列）
口径：与 33_std_ci 相同的权重与样本过滤，但不做年龄标准化（波次内样本较小）。
输出：45_wave_weighted.csv
"""
import os, numpy as np, pandas as pd, pyreadstat

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))

W = pd.read_csv(os.path.join(OUT, "16_weights.csv"), dtype={"ID": str})[["cohort", "ID", "wave", "weight"]]

SPEC = {
    "CHARLS": dict(p=os.path.join(DL, "1. CHARLS  中国", "CHARLS_中国", "Working_data", "charls.dta"),
                   idc="ID", base=["ID", "wave", "age", "arthre", "rxarthr_c"],
                   sites=["da042s2", "da042s3", "da042s4", "da042s5", "da042s8", "da042s9",
                          "da042s10", "da042s11", "da042s12", "da042s13", "da042s14", "da042s15"],
                   mode="list", dx="arthre", tx="rxarthr_c", age="age", waves=[1, 2, 3, 4]),
    "ELSA": dict(p=os.path.join(DL, "2. ELSA 英国", "Working_data", "elsa.dta"),
                 idc="idauniqc", base=["idauniqc", "wave", "agey", "arthre", "painfr"],
                 sites=[], mode="single", pain="painfr", dx="arthre", tx=None, age="agey", waves=None),
    "HRS": dict(p=os.path.join(DL, "3. HRS  美国", "HRS_美国", "Working_data", "hrs.dta"),
                idc="hhidpn", base=["hhidpn", "wave", "ragey_b", "arthre", "rxarthr", "painfr"],
                sites=[], mode="single", pain="painfr", dx="arthre", tx="rxarthr", age="ragey_b", waves=None),
    "KLoSA": dict(p=os.path.join(DL, "4. KLoSA 韩国", "KLoSA_韩国", "Working_data", "klosa.dta"),
                  idc="pid", base=["pid", "wave", "agey", "arthre", "rxarthr"],
                  sites=["pain_2", "pain_3", "pain_4", "pain_5", "pain_8",
                         "pain_9", "pain_10", "pain_11", "pain_12", "pain_13"],
                  mode="list", dx="arthre", tx="rxarthr", age="agey", waves=None),
    "LASI": dict(p=os.path.join(DL, "5. LASI 印度", "LASI_印度", "Working_data", "lasi.dta"),
                 idc="prim_key", base=["prim_key", "r1agey", "r1arthre", "r1rxarthr", "r1jointp"],
                 sites=[], mode="joint", pain="r1jointp", dx="r1arthre", tx="r1rxarthr", age="r1agey", waves=None),
    "MHAS": dict(p=os.path.join(DL, "6. MHAS 墨西哥", "MHAS_墨西哥", "Working_data", "mhas.dta"),
                 idc="rahhidnp", base=["rahhidnp", "wave", "agey", "arthre", "rxarthr", "painfr"],
                 sites=[], mode="single", pain="painfr", dx="arthre", tx="rxarthr", age="agey", waves=None),
    "SHARE": dict(p=os.path.join(DL, "7.SHARE 欧洲", "SHARE_欧洲", "Working_data", "share.dta"),
                  idc="mergeid", base=["mergeid", "wave", "agey", "arthre", "rxpain", "painfr", "pain_s"],
                  sites=[], mode="single", pain="painfr", dx="arthre", tx="rxpain", age="agey", waves=None),
}


def load(c, s, tier):
    cols = list(dict.fromkeys(s["base"] + s["sites"]))
    df, _ = pyreadstat.read_dta(s["p"], usecols=cols)
    if s["idc"] != "ID":
        df = df.rename(columns={s["idc"]: "ID"})
    d = df.copy()
    if "wave" not in d.columns:
        d["wave"] = 1
    if s["waves"] is not None:
        d = d[d["wave"].isin(s["waves"])]
    if c == "SHARE":
        d = d[d["wave"].isin([1, 2, 4] if tier == "B" else [5, 6, 7, 8])]
    d["ID"] = d["ID"].astype(str)
    d["age"] = pd.to_numeric(d[s["age"]], errors="coerce")
    d = d[(d["age"] >= 50) & (d["age"] <= 105)].copy()
    if tier == "B":
        if c == "LASI":
            pv = pd.to_numeric(d["r1jointp"], errors="coerce")
        elif c == "SHARE":
            pv = pd.to_numeric(d["pain_s"], errors="coerce")
        else:
            sub = d[s["sites"]]
            pv = (sub == 1).any(axis=1).astype(float)
            pv[sub.isna().all(axis=1)] = np.nan
    else:
        pv = pd.to_numeric(d[s["pain"]], errors="coerce")
    d["pain"] = np.where(pv.notna(), np.where(pv > 0, 1.0, 0.0), np.nan)
    d["dx"] = pd.to_numeric(d[s["dx"]], errors="coerce").where(lambda x: x.isin([0, 1]))
    d["tx"] = (pd.to_numeric(d[s["tx"]], errors="coerce").where(lambda x: x.isin([0, 1]))
               if s["tx"] else np.nan)
    d["wave"] = d["wave"].astype(float)
    d = d.merge(W[W.cohort == c][["ID", "wave", "weight"]], on=["ID", "wave"], how="left")
    return d


def wrate(sub, col):
    s = sub[sub[col].notna() & sub["weight"].notna()]
    if s.empty or s["weight"].sum() <= 0:
        return np.nan, 0
    return float(np.average(s[col], weights=s["weight"]) * 100), len(s)


rows = []
for c, s in SPEC.items():
    tiers = ["A", "B"] if c == "SHARE" else (["B"] if c in ("CHARLS", "KLoSA", "LASI") else ["A"])
    for tier in tiers:
        d = load(c, s, tier)
        d = d[d["weight"].notna()]
        for wv, g in d.groupby("wave"):
            pain, n0 = wrate(g, "pain")
            gp = g[g["pain"] == 1]
            dx, n1 = wrate(gp, "dx")
            gd = gp[gp["dx"] == 1]
            tx, n2 = wrate(gd, "tx") if "tx" in gd and gd["tx"].notna().any() else (np.nan, 0)
            n_dx_valid = int(gp["dx"].notna().sum())
            rows.append(dict(cohort=c, tier="A" if tier == "A" else "B", wave=int(wv),
                             n=len(g), n_pain=len(gp), n_dx=len(gd), n_dx_valid=n_dx_valid,
                             pain_w=round(pain, 1) if np.isfinite(pain) else np.nan,
                             dx_w=round(dx, 1) if np.isfinite(dx) else np.nan,
                             tx_w=round(tx, 1) if np.isfinite(tx) else np.nan))
r = pd.DataFrame(rows)
r.to_csv(os.path.join(OUT, "45_wave_weighted.csv"), index=False, encoding="utf-8-sig")
print(r.to_string(index=False))
print("WROTE 45_wave_weighted.csv")
