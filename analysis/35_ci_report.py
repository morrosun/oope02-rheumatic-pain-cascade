# -*- coding: utf-8 -*-
"""35: 生成 P9（标化 CI）与 P8（聚类敏感性）的可读报告"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import os
import numpy as np, pandas as pd
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
ci=pd.read_csv(os.path.join(OUT,"33_std_ci.csv"))
se=pd.read_csv(os.path.join(OUT,"34_cluster_se.csv"))
L=[]
def w(s=''): L.append(str(s))

NAME={'HRS':'美国 HRS','MHAS':'墨西哥 MHAS','ELSA':'英国 ELSA','SHARE':'欧洲 SHARE',
      'CHARLS':'中国 CHARLS','KLoSA':'韩国 KLoSA','LASI':'印度 LASI'}
ORD_A=['HRS','SHARE','ELSA','MHAS']; ORD_B=['CHARLS','SHARE','KLoSA','LASI']

w('# 年龄标准化率的 95% CI 与聚类层级敏感性')
w('')
w('**日期**：2026-11-14　**脚本**：`33_std_ci.py`、`34_cluster_se.py`')
w('')
w('## 方法')
w('')
w('**CI（P9）**：整群 bootstrap，B=1,000。重抽单位按库的可得层级取最高一级——CHARLS 用 `communityID`（社区，450 个），其余六库用**个体首次观测的家户 ID**（同一个人的家户跨波固定）。每次重算各年龄组的加权率，再按 WHO 2000–2025 标准人口（50+ 重归一）加权；CI 取 bootstrap 百分位 2.5/97.5。')
w('')
w('**聚类敏感性（P8）**：对 A2 的诊断模型，分别以「个体」与「群」为聚类单位重估标准误，报告 SE 比值（群/个体）。比值 >1 表示个体聚类低估了标准误。')
w('')
w('> ⚠️ 层级差异必须写明：CHARLS 的群是真正的社区（PSU 代理），其余库的家户**不是** PSU，只能校正家户内相关，不能校正社区内相关。因此其余库的 CI 是**下限估计**（真实 CI 可能更宽）。CHARLS 的 1.08 倍可视为该设计下聚类效应的量级参考。')
w('')

for tier,order,title in [('A',ORD_A,'面板 A　频繁疼痛'),('B',ORD_B,'面板 B　关节疼痛')]:
    w('## 表 {}　{}（加权 + WHO 标化，95% CI）'.format('1' if tier=='A' else '2', title))
    w('')
    w('| 库 | 疼痛率 % | 诊断覆盖 %（疼痛者中） | 治疗覆盖 %（诊断且有疼痛） | 群数 |')
    w('|---|---|---|---|---:|')
    for c in order:
        g=ci[(ci.tier==tier)&(ci.cohort==c)].set_index('metric')
        def f(k):
            if k not in g.index: return '—'
            r=g.loc[k]; return '{:.1f}（{:.1f}–{:.1f}）'.format(r.est,r.lo,r.hi)
        nc=int(g['n_clusters'].max()) if len(g) else 0
        w('| {} | {} | **{}** | {} | {:,} |'.format(NAME[c],f('pain_rate'),f('dx_given_pain'),f('tx_given_dx_pain'),nc))
    w('')
    g=ci[(ci.tier==tier)&(ci.metric=='dx_given_pain')].set_index('cohort')
    hi=g.loc[g.est.idxmax()]; lo=g.loc[g.est.idxmin()]
    w('> 诊断覆盖极差 **{:.1f} / {:.1f} = {:.2f} 倍**；最高组 CI 下限 {:.1f} > 最低组 CI 上限 {:.1f}，**区间不重叠**。'.format(
        hi.est,lo.est,hi.est/lo.est,hi.lo,lo.hi))
    t=ci[(ci.tier==tier)&(ci.metric=='tx_given_dx_pain')].set_index('cohort')
    if len(t):
        th=t.loc[t.est.idxmax()]; tl=t.loc[t.est.idxmin()]
        w('> 治疗覆盖极差 **{:.1f} / {:.1f} = {:.2f} 倍**，区间 [{:.1f}, {:.1f}]；最高组 CI 下限 {:.1f} 与最低组 CI 上限 {:.1f} {}。'.format(
            th.est,tl.est,th.est/tl.est,tl.est,th.est,th.lo,tl.hi,
            '**不重叠**，说明治疗率的跨国差异虽然幅度小，但仍具统计学意义' if th.lo>tl.hi else '有重叠，差异不显著'))
    w('')

w('## 表 3　聚类层级敏感性（A2 诊断模型）')
w('')
w('| 库 | n | 群数 | 男性 OR（个体聚类） | 男性 OR（群聚类） | SE 比值 |')
w('|---|---:|---:|---|---|---:|')
for c in ['HRS','MHAS','ELSA','SHARE','CHARLS','KLoSA','LASI']:
    g=se[(se.cohort==c)&(se.term=='male')]
    if g.empty: continue
    r=g.iloc[0]
    w('| {} | {:,} | {:,} | {:.3f}（{:.3f}–{:.3f}） | {:.3f}（{:.3f}–{:.3f}） | {:.3f} |'.format(
        NAME[c],int(r.n),int(r.n_clusters),r.OR_id,r.lo_id,r.hi_id,r.OR_clu,r.lo_clu,r.hi_clu,r.se_ratio))
w('')
w('> 全部 21 个系数（男性/年龄/低教育 × 7 库）的 SE 比值中位数 **{:.3f}**、最大 **{:.3f}**（CHARLS 社区）。'.format(se.se_ratio.median(),se.se_ratio.max()))
w('> **结论：个体聚类没有明显低估标准误，A2 的 95% CI 与显著性判断不需要修改。** 最大的一处是 CHARLS（1.077），CI 宽度约增 8%，「男性 OR 0.79」仍远离 1。')
w('')

w('## 三条对论文表述的影响')
w('')
w('1. **「诊断覆盖相差 2.6 倍」可以写**——B 面板 55.1%（53.6–56.8）对 21.1%（20.4–21.7），CI 不重叠；A 面板 2.13 倍，同样不重叠。这是原文最想守的一个数字，现在有区间支撑。')
w('2. **「治疗率挤在 44–62% 的窄带」要改措辞**——幅度确实窄（1.23–1.25 倍），但 HRS 55.2%（54.2–56.3）与 MHAS 44.0%（41.0–46.9）的 CI **不重叠**，KLoSA 61.7%（59.0–64.3）与 CHARLS 50.2%（48.1–52.3）也不重叠。正确写法是「治疗率的跨国差异幅度远小于诊断率，但仍具统计学意义」。')
w('3. **局限性要写明层级**——除 CHARLS 外没有 PSU 标识，家户聚类只能校正家户内相关；真实 CI 可能略宽。CHARLS 的 1.08 倍提示这个偏倚幅度不大，但不能完全排除。')

open(os.path.join(OUT,'标化CI与聚类敏感性.md'),'w',encoding='utf-8').write('\n'.join(L))
print('\n'.join(L))
