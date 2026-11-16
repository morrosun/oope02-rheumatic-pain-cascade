# -*- coding: utf-8 -*-
"""32: Readable tables for the functional-outcome model (31)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import pandas as pd, numpy as np
res=pd.read_csv(os.path.join(OUT,"31_outcome_v2.csv"))
pool=pd.read_csv(os.path.join(OUT,"31_pooled_v2.csv"))
L=[]
def w(s=""): L.append(str(s))

w('# Functional limitation across cascade states')
w('')
w('> Reference group: no pain and no diagnosis.')
w('')
w('**口径**：结局 = ADL 任何受限（0/1）；每库 logistic，个体聚类稳健；跨库 DerSimonian–Laird 合并 log OR；协变量 = 年龄（每 10 岁）、性别、教育（低/中 vs 高）、财富五分位（cohort×wave）、城乡、共病数、吸烟。未加权（见方法说明）。')
w('')

w('## 表 1　三状态：参照 = 无痛且未诊断（七库构念一致）')
w('')
w('| 面板 | 状态 | 合并 OR（95% CI） | k | I² | 参与库 |')
w('|---|---|---|---:|---:|---|')
lab={'s1':'痛 · 未诊断','s2':'痛 · 已诊断','s3':'无痛 · 已诊断'}
for panel,pn in [('A','A 频繁疼痛'),('B','B 关节疼痛')]:
    sub=pool[(pool.panel==panel)&(pool['mode']=='three')].sort_values('term')
    for _,r in sub.iterrows():
        w('| {} | {} | **{:.2f}（{:.2f}–{:.2f}）** | {} | {:.0f} | {} |'.format(pn,lab.get(r['term'],r['term']),r.OR,r.lo,r.hi,int(r.k),r.I2,r.cohorts))
w('')
w('> 读法：未诊断的疼痛者承担了已诊断组约 **60%（A，4.057/6.713）、65%（B，2.040/3.133）** 的额外功能负担——漏诊人群不是轻症人群。')
w('> 「无痛但已诊断」组 ADL 仍显著更差（A 1.96、B 1.87），说明**曾获诊断本身是功能受损的长期标记**，与当前是否疼痛无关。这一格只有在参照组排除「无痛但已诊断」者后才可见。')
w('')

w('## 表 2　四状态（仅 5 个有治疗变量的库，独立报告，不与表 1 合并）')
w('')
w('| 库 | n | 痛·未诊断 | 痛·诊断未治疗 | 痛·诊断在治疗 | 无痛·已诊断 |')
w('|---|---:|---:|---:|---:|---:|')
for c in ['HRS','MHAS','CHARLS','KLoSA','LASI']:
    sub=res[(res.cohort==c)&(res['mode']=='four')].set_index('term')
    if sub.empty: continue
    n=int(sub['n'].iloc[0])
    def g(k):
        return '{:.2f}'.format(sub.loc[k,'OR']) if k in sub.index else '—'
    w('| {} | {:,} | {} | {} | **{}** | {} |'.format(c,n,g('s1'),g('s2'),g('s3'),g('s4')))
w('')
w('> 「诊断在治疗」组在 5 个库中**全部最高**，**适应证混杂是结构性的**——病重的人才吃药。A3 绝不能把它读成「治疗有害」。但 KLoSA 的「诊断未治疗」低于「未诊断」，使梯度非单调。')
w('> **KLoSA 非单调**：原始 ADL 率 2.4→5.1→5.5 单调，但年龄调整后「诊断未治疗」OR 仅 0.96；两组年龄中位数差 5 岁（68 vs 73）。表述必须写「除 KLoSA 外单调」。')
w('')

w('## 表 3　疼痛者中「已诊断 vs 未诊断」')
w('')
w('| 面板 | 未调严重度 | 调整疼痛严重度 | 备注 |')
w('|---|---|---|---|')
for panel,pn in [('A','A'),('B','B')]:
    a=pool[(pool.panel==panel)&(pool['mode']=='M2_no_sev')].iloc[0]
    b=pool[(pool.panel==panel)&(pool['mode']=='M2_adj_sev')]
    if len(b):
        b=b.iloc[0]
        bt='**{:.2f}（{:.2f}–{:.2f}）**，k={}'.format(b.OR,b.lo,b.hi,int(b.k))
        note='显著' if b.lo>1 else '**不显著（CI 跨 1）**'
    else:
        bt='—'; note=''
    if panel=='B': note+='；严重度调整仅 CHARLS/KLoSA/LASI 三库（SHARE 的 `painlv` 在 w1/2/4 不存在）'
    w('| {} | {:.2f}（{:.2f}–{:.2f}），k={} | {} | {} |'.format(pn,a.OR,a.lo,a.hi,int(a.k),bt,note))
w('')
w('> **B 面板调整严重度后不再显著（1.36，1.00–1.84）**，原报告的 1.46（1.19–1.81）依赖了 SHARE 那个**未真正调整**的估计（其 `M2_adj_sev` 与 `M2_no_sev` 数值完全相同，均为 1.804）。')
w('> 这说明「已诊断者功能更差」这个关联的稳健性**依赖于疼痛测量方式**——A 面板（频繁疼痛）稳健，B 面板（关节痛）不稳健。讨论里必须写这一层。')
w('')

w('## 表 4　逐库明细（三状态）')
w('')
w('| 面板 | 库 | n | 痛·未诊断 | 痛·已诊断 | 无痛·已诊断 |')
w('|---|---|---:|---:|---:|---:|')
for panel in ['A','B']:
    for c in (['HRS','MHAS','ELSA','SHARE'] if panel=='A' else ['CHARLS','KLoSA','LASI','SHARE']):
        sub=res[(res.panel==panel)&(res.cohort==c)&(res['mode']=='three')].set_index('term')
        if sub.empty: continue
        n=int(sub['n'].iloc[0])
        def g(k):
            if k not in sub.index: return '—'
            r=sub.loc[k]; return '{:.2f}（{:.2f}–{:.2f}）'.format(r.OR,r.lo,r.hi)
        w('| {} | {} | {:,} | {} | {} | {} |'.format(panel,c,n,g('s1'),g('s2'),g('s3')))
w('')

w('## 三条限制（写进 methods / limitations）')
w('')
w('1. **横断面**：诊断与疼痛同时测量，「已诊断」与「病更重/更常就医」无法分离，OR>1 只能读作**标记**而非效应。')
w('2. **未加权**：描述性统计用抽样权重 + WHO 标化，回归用未加权以估计条件效应——两者分工需在方法里写明。')
w('3. **聚类层级**：标准误仅按个体聚类，未按社区/PSU。CHARLS、LASI、SHARE 为多阶段整群抽样，CI 可能偏窄。')

open(os.path.join(OUT,'31_outcome_table.md'),'w',encoding='utf-8').write('\n'.join(L))
print('\n'.join(L))
