# -*- coding: utf-8 -*-
"""PubMed novelty scan round 2: tighten the surviving candidates."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import urllib.request, urllib.parse, json, time, os

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = [
 # 1. 未诊断 / 症状-诊断缺口
 ("P1","未诊断关节炎 标题级","undiagnosed[ti] AND (arthritis[ti] OR arthralgia[ti] OR \"joint pain\"[ti])",""),
 ("P2","症状未就医/未诊断","(\"joint pain\"[tiab] OR arthralgia[tiab]) AND (\"did not consult\"[tiab] OR \"not consulted\"[tiab] OR \"no medical diagnosis\"[tiab] OR \"unmet need\"[tiab])",""),
 ("P3","关节炎 诊断延迟","arthritis[tiab] AND (\"diagnostic delay\"[tiab] OR \"delay in diagnosis\"[tiab] OR \"time to diagnosis\"[tiab]) AND (population[tiab] OR cohort[tiab])",""),
 ("P4","慢性疼痛 未诊断 老年","(\"chronic pain\"[tiab]) AND (undiagnosed[tiab] OR \"no diagnosis\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab])",""),
 # 2. 诊疗级联
 ("P5","care cascade 肌骨","(\"care cascade\"[tiab] OR \"treatment cascade\"[tiab] OR \"cascade of care\"[tiab]) AND (arthritis[tiab] OR musculoskeletal[tiab] OR rheumat*[tiab])",""),
 ("P6","treatment gap 关节炎 标题","(\"treatment gap\"[ti]) AND (arthritis[ti] OR rheumat*[ti] OR gout[ti] OR osteoporos*[ti])",""),
 ("P7","未治疗关节炎","arthritis[tiab] AND (untreated[tiab] OR \"not treated\"[tiab] OR undertreated[tiab]) AND (\"older adults\"[tiab] OR population[tiab] OR cohort[tiab])",""),
 # 3. 止痛药/抗炎药使用
 ("P8","止痛药使用 关节炎 人群","(analgesic*[tiab] OR \"pain medication\"[tiab]) AND arthritis[tiab] AND (\"population-based\"[tiab] OR \"nationally representative\"[tiab] OR survey[tiab])",""),
 ("P9","SHARE 止痛药","(SHARE[tiab] OR CHARLS[tiab] OR HRS[tiab] OR LASI[tiab]) AND (analgesic*[tiab] OR \"pain medication\"[tiab] OR \"anti-inflammatory\"[tiab])",""),
 # 4. 疼痛干扰 / 多部位疼痛
 ("P10","疼痛干扰 死亡","(\"pain interference\"[tiab]) AND (mortality[tiab] OR death[tiab]) AND (cohort[tiab] OR longitudinal[tiab])",""),
 ("P11","CHARLS 疼痛","CHARLS[tiab] AND (pain[tiab] OR \"painful\"[tiab]) AND (\"chronic pain\"[tiab] OR \"joint pain\"[tiab] OR \"pain sites\"[tiab])",""),
 ("P12","疼痛部位 聚类","(\"pain sites\"[tiab] OR \"pain locations\"[tiab] OR \"multisite\"[tiab]) AND (\"latent class\"[tiab] OR cluster[tiab] OR pattern*[tiab]) AND (\"older adults\"[tiab] OR cohort[tiab])",""),
 # 5. 数字
 ("P13","互联网使用 关节炎","(\"internet use\"[tiab] OR \"internet access\"[tiab] OR \"digital divide\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab] OR \"musculoskeletal\"[tiab])",""),
 # 6. 疫苗
 ("P14","流感疫苗 关节炎 人群","(\"influenza vaccination\"[tiab] OR \"flu shot\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab]) AND (\"older adults\"[tiab] OR population[tiab] OR cohort[tiab])",""),
 # 7. 诊断年龄
 ("P15","诊断年龄 关节炎 人群","(\"age at diagnosis\"[tiab]) AND (arthritis[tiab]) AND (\"older adults\"[tiab] OR population[tiab] OR cohort[tiab]) AND (disability[tiab] OR mortality[tiab] OR function*[tiab])",""),
 # 8. 健康期望寿命
 ("P16","关节炎 无失能期望寿命","arthritis[tiab] AND (\"disability-free life expectancy\"[tiab] OR \"healthy life expectancy\"[tiab] OR \"health expectancy\"[tiab] OR \"Sullivan\"[tiab])",""),
 # 9. 关节置换未满足
 ("P17","关节置换 需求未满足 人群","(\"joint replacement\"[tiab] OR \"hip replacement\"[tiab] OR \"knee replacement\"[tiab]) AND (\"unmet need\"[tiab] OR \"unmet\"[tiab]) AND (survey[tiab] OR cohort[tiab] OR population[tiab])",""),
 # 10. 其他候选
 ("P18","关节炎 尿失禁 队列","(arthritis[tiab]) AND (incontinence[tiab]) AND (cohort[tiab] OR longitudinal[tiab] OR CHARLS[tiab] OR SHARE[tiab])",""),
 ("P19","关节炎 牙齿 队列","(arthritis[tiab]) AND (teeth[tiab] OR edentul*[tiab] OR \"tooth loss\"[tiab]) AND (CHARLS[tiab] OR SHARE[tiab] OR HRS[tiab] OR LASI[tiab] OR ELSA[tiab])",""),
 ("P20","关节炎 感官 损失 队列","(arthritis[tiab]) AND (\"hearing\"[tiab] OR \"vision\"[tiab]) AND (CHARLS[tiab] OR SHARE[tiab] OR HRS[tiab] OR LASI[tiab] OR ELSA[tiab])",""),
 # 11. 我们的 OOPE 论文占位校验
 ("P21","arthritis OOPE 四库","arthritis[tiab] AND (\"out-of-pocket\"[tiab]) AND (MHAS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab])",""),
 # 12. 多库 harmonized 方法学占用
 ("P22","harmonized 多库 肌骨","(harmonized[tiab] OR harmonised[tiab] OR \"Gateway to Global Aging\"[tiab]) AND (musculoskeletal[tiab] OR arthritis[tiab] OR rheumat*[tiab])",""),
 ("P23","CHARLS 关节炎 全部","CHARLS[tiab] AND arthritis[tiab]",""),
 ("P24","SHARE 关节炎 全部","SHARE[tiab] AND arthritis[tiab] AND (wave[tiab] OR Europ*[tiab])",""),
]

def esearch(term, retmax=6):
    url = BASE + "esearch.fcgi?" + urllib.parse.urlencode({
        "db":"pubmed","term":term,"retmode":"json","retmax":retmax,"sort":"relevance"})
    last=None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.loads(r.read().decode())["esearchresult"]
        except Exception as e:
            last=e; time.sleep(2+attempt*2)
    raise last

def esummary(pmids):
    if not pmids: return {}
    url = BASE + "esummary.fcgi?" + urllib.parse.urlencode({"db":"pubmed","id":",".join(pmids),"retmode":"json"})
    with urllib.request.urlopen(url, timeout=40) as r:
        return json.loads(r.read().decode()).get("result", {})

lines=[]; js=[]
for qid, note, term, _ in QUERIES:
    try:
        res=esearch(term); cnt=int(res.get("count","0")); ids=res.get("idlist",[])
        summ=esummary(ids); recs=[]
        for p in ids:
            it=summ.get(p,{})
            recs.append({"pmid":p,"title":it.get("title",""),"journal":it.get("source",""),"year":(it.get("pubdate","") or "")[:4]})
        lines.append(f"\n[{qid}] {note}\n  hits={cnt}\n  query={term}")
        for r_ in recs:
            lines.append(f"    - {r_['pmid']} | {r_['year']} | {r_['journal']} | {r_['title'][:190]}")
        js.append({"id":qid,"note":note,"query":term,"count":cnt,"top":recs})
        print(qid,cnt,flush=True)
    except Exception as e:
        lines.append(f"\n[{qid}] ERROR {e}"); print(qid,"ERR",e,flush=True)
    time.sleep(0.45)

with open(os.path.join(OUT,"03_pubmed_r2.txt"),"w",encoding="utf-8") as f:
    f.write("PubMed novelty scan round 2 | 2026-11-14\n"+"\n".join(lines))
with open(os.path.join(OUT,"03_pubmed_r2.json"),"w",encoding="utf-8") as f:
    json.dump(js,f,ensure_ascii=False,indent=1)
print("WROTE",len(js))
