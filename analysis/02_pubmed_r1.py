# -*- coding: utf-8 -*-
"""PubMed novelty scan for AGING_DB x rheumatology directions (round 1)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import urllib.request, urllib.parse, json, time, os, sys, re

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = [
 # --- A. 诊疗缺口 / 未诊断 ---
 ("A1","A","(arthritis[tiab] OR rheumatic[tiab]) AND (undiagnosed[tiab] OR underdiagnos*[tiab] OR \"diagnostic gap\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab] OR cohort[tiab])","未诊断关节炎"),
 ("A2","A","(arthritis[tiab] OR rheumat*[tiab]) AND (\"treatment gap\"[tiab] OR \"care cascade\"[tiab] OR \"treatment cascade\"[tiab] OR \"continuum of care\"[tiab])","风湿治疗级联"),
 ("A3","A","(\"joint pain\"[tiab] OR arthralgia[tiab]) AND (undiagnosed[tiab] OR \"no diagnosis\"[tiab] OR \"without a diagnosis\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","未诊断关节痛"),
 ("A4","A","(CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab]) AND (\"joint pain\"[tiab] OR arthralgia[tiab]) AND (diagnos*[tiab] OR arthrit*[tiab])","多库 关节痛/诊断"),
 ("A5","A","arthritis[tiab] AND (analgesic*[tiab] OR \"pain medication\"[tiab] OR \"pain reliever\"[tiab]) AND (\"disease-modifying\"[tiab] OR \"arthritis medication\"[tiab] OR \"anti-inflammatory\"[tiab])","止痛药 vs 抗风湿药"),
 ("A6","A","(\"unmet need\"[tiab]) AND (rheumat*[tiab] OR arthritis[tiab]) AND (\"older adults\"[tiab] OR aged[tiab]) AND (care[tiab] OR treatment[tiab])","风湿未满足需求"),
 # --- B. 照护供给 ---
 ("B1","B","arthritis[tiab] AND (caregiving[tiab] OR caregiver*[tiab]) AND (\"older adults\"[tiab] OR aged[tiab]) AND (provide*[tiab] OR \"care provision\"[tiab])","关节炎 提供照护"),
 ("B2","B","arthritis[tiab] AND (\"grandchild care\"[tiab] OR grandparent*[tiab] OR \"child care\"[tiab])","关节炎 隔代照护"),
 ("B3","B","(musculoskeletal[tiab] OR arthritis[tiab]) AND (\"informal care\"[tiab]) AND (provid*[tiab] OR hours[tiab]) AND (cohort[tiab] OR longitudinal[tiab])","肌骨病 非正式照护"),
 # --- C. 数字/互联网 ---
 ("C1","C","arthritis[tiab] AND (\"internet use\"[tiab] OR \"digital divide\"[tiab] OR \"digital exclusion\"[tiab] OR \"eHealth\"[tiab] OR telehealth[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","关节炎 数字使用"),
 ("C2","C","(CHARLS[tiab] OR SHARE[tiab] OR HRS[tiab] OR LASI[tiab] OR ELSA[tiab] OR MHAS[tiab]) AND (internet[tiab] OR \"digital\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab])","多库 互联网×风湿"),
 # --- D. 疼痛表型 ---
 ("D1","D","(\"pain interference\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab]) AND (cohort[tiab] OR longitudinal[tiab] OR \"older adults\"[tiab])","疼痛干扰×风湿"),
 ("D2","D","(\"multisite pain\"[tiab] OR \"number of pain sites\"[tiab] OR \"widespread pain\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab]) AND (mortality[tiab] OR disability[tiab] OR frailty[tiab])","多部位疼痛 结局"),
 ("D3","D","(\"pain interference\"[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab])","疼痛干扰 多库"),
 # --- E. 疫苗 ---
 ("E1","E","arthritis[tiab] AND (vaccination[tiab] OR immunization[tiab] OR \"vaccine uptake\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","关节炎 疫苗接种"),
 ("E2","E","(\"rheumatic disease\"[tiab] OR \"inflammatory arthritis\"[tiab] OR \"rheumatoid arthritis\"[tiab]) AND (\"vaccination coverage\"[tiab] OR \"vaccine uptake\"[tiab] OR \"immunization coverage\"[tiab])","风湿病 疫苗覆盖率"),
 ("E3","E","pneumococcal[tiab] AND (arthritis[tiab] OR rheumat*[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","肺炎疫苗×风湿"),
 # --- F. 口腔 ---
 ("F1","F","arthritis[tiab] AND (\"tooth loss\"[tiab] OR edentulism[tiab] OR \"number of teeth\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab] OR cohort[tiab])","关节炎 牙齿脱落"),
 ("F2","F","(CHARLS[tiab] OR SHARE[tiab] OR HRS[tiab] OR LASI[tiab] OR ELSA[tiab]) AND (teeth[tiab] OR dent*[tiab]) AND (arthritis[tiab] OR rheumat*[tiab])","多库 口腔×风湿"),
 # --- G. 尿失禁 ---
 ("G1","G","arthritis[tiab] AND (\"urinary incontinence\"[tiab] OR incontinence[tiab]) AND (\"older adults\"[tiab] OR aged[tiab] OR women[tiab])","关节炎 尿失禁"),
 # --- H. 诊断年龄/病程 ---
 ("H1","H","(\"age at diagnosis\"[tiab] OR \"age at onset\"[tiab] OR \"young-onset\"[tiab] OR \"late-onset\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab]) AND (mortality[tiab] OR disability[tiab] OR \"older adults\"[tiab])","发病年龄 结局"),
 ("H2","H","(\"disease duration\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab]) AND (\"population-based\"[tiab] OR cohort[tiab]) AND (mortality[tiab] OR disability[tiab] OR function*[tiab])","病程 结局"),
 # --- I. 关节置换 ---
 ("I1","I","(\"joint replacement\"[tiab] OR arthroplasty[tiab]) AND (\"unmet need\"[tiab] OR unmet[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","关节置换未满足"),
 ("I2","I","(\"joint replacement\"[tiab] OR arthroplasty[tiab]) AND (\"cross-national\"[tiab] OR \"cross-country\"[tiab] OR international[tiab]) AND (arthritis[tiab] OR osteoarthr*[tiab])","关节置换 跨国"),
 # --- J. 社会隔离/孤独 ---
 ("J1","J","arthritis[tiab] AND (\"social isolation\"[tiab] OR loneliness[tiab] OR lonely[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","关节炎 社会隔离"),
 # --- K. 睡眠 (占用校验) ---
 ("K1","K","arthritis[tiab] AND (sleep[tiab] OR insomnia[tiab]) AND (CHARLS[tiab] OR SHARE[tiab] OR HRS[tiab] OR ELSA[tiab] OR MHAS[tiab])","关节炎 睡眠 多库"),
 # --- L. 感官 ---
 ("L1","L","arthritis[tiab] AND (\"hearing loss\"[tiab] OR \"vision impairment\"[tiab] OR \"dual sensory\"[tiab] OR \"sensory impairment\"[tiab])","关节炎 感官"),
 # --- M. 主观存活预期 ---
 ("M1","M","arthritis[tiab] AND (\"subjective survival\"[tiab] OR \"survival expectations\"[tiab] OR \"life expectancy\"[tiab])","关节炎 存活预期"),
 # --- N. 领域占用校验 ---
 ("N1","N","arthritis[tiab] AND (\"out-of-pocket\"[tiab] OR \"catastrophic health expenditure\"[tiab] OR \"medical expenditure\"[tiab])","关节炎 自付(占用)"),
 ("N2","N","arthritis[tiab] AND (falls[tiab] OR \"fall\"[ti]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab])","关节炎 跌倒(占用)"),
 ("N3","N","arthritis[tiab] AND frailty[tiab] AND (\"older adults\"[tiab] OR cohort[tiab])","关节炎 衰弱(占用)"),
 ("N4","N","(arthritis[ti]) AND (harmonized[tiab] OR harmonised[tiab] OR \"Gateway to Global Aging\"[tiab] OR \"HRS family\"[tiab])","harmonized 关节炎"),
 ("N5","N","arthritis[tiab] AND (CHARLS[tiab] AND HRS[tiab] AND SHARE[tiab])","三大库 关节炎"),
 ("N6","N","(\"rheumatoid arthritis\"[tiab] OR gout[tiab] OR \"systemic lupus\"[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab])","具体风湿病 老龄化库"),
 # --- O. 其他可能空白 ---
 ("O1","O","arthritis[tiab] AND (\"oral health\"[tiab] OR periodont*[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","关节炎 牙周"),
 ("O2","O","(arthritis[tiab] OR rheumat*[tiab]) AND (\"health literacy\"[tiab] OR numeracy[tiab])","风湿 健康素养"),
 ("O3","O","(arthritis[tiab] OR rheumat*[tiab]) AND (\"driving cessation\"[tiab] OR driving[tiab])","风湿 驾驶"),
 ("O4","O","(arthritis[tiab] OR rheumat*[tiab]) AND (\"traditional medicine\"[tiab] OR \"complementary\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","风湿 传统医学"),
 ("O5","O","(arthritis[tiab] OR rheumat*[tiab]) AND (\"social participation\"[tiab]) AND (\"older adults\"[tiab] OR cohort[tiab])","风湿 社会参与"),
 ("O6","O","(arthritis[tiab] OR rheumat*[tiab]) AND (polypharmacy[tiab] OR \"medication burden\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])","风湿 多重用药"),
]

def esearch(term, retmax=6):
    url = BASE + "esearch.fcgi?" + urllib.parse.urlencode({
        "db":"pubmed","term":term,"retmode":"json","retmax":retmax,"sort":"relevance"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                d = json.loads(r.read().decode())
            return d["esearchresult"]
        except Exception as e:
            time.sleep(2+attempt*2)
            err = e
    raise err

def esummary(pmids):
    if not pmids: return {}
    url = BASE + "esummary.fcgi?" + urllib.parse.urlencode({
        "db":"pubmed","id":",".join(pmids),"retmode":"json"})
    with urllib.request.urlopen(url, timeout=40) as r:
        d = json.loads(r.read().decode())
    return d.get("result", {})

lines = []; js = []
for qid, grp, term, note in QUERIES:
    try:
        res = esearch(term)
        cnt = int(res.get("count","0")); ids = res.get("idlist", [])
        summ = esummary(ids)
        recs = []
        for p in ids:
            it = summ.get(p, {})
            recs.append({"pmid":p,"title":it.get("title",""),"journal":it.get("source",""),
                         "year":(it.get("pubdate","") or "")[:4]})
        lines.append(f"\n[{qid}] {note}\n  hits={cnt}\n  query={term}")
        for r_ in recs:
            lines.append(f"    - {r_['pmid']} | {r_['year']} | {r_['journal']} | {r_['title'][:180]}")
        js.append({"id":qid,"group":grp,"note":note,"query":term,"count":cnt,"top":recs})
        print(qid, cnt, flush=True)
    except Exception as e:
        lines.append(f"\n[{qid}] ERROR {e}")
        print(qid, "ERR", e, flush=True)
    time.sleep(0.45)

with open(os.path.join(OUT,"02_pubmed_r1.txt"),"w",encoding="utf-8") as f:
    f.write("PubMed novelty scan round 1 | run 2026-11-14\n"+"\n".join(lines))
with open(os.path.join(OUT,"02_pubmed_r1.json"),"w",encoding="utf-8") as f:
    json.dump(js,f,ensure_ascii=False,indent=1)
print("WROTE", len(js))
