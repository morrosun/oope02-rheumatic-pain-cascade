# -*- coding: utf-8 -*-
"""Round 4: title-level checks for shortlisted directions + method precedent."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import urllib.request, urllib.parse, json, time, os

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = [
 ("S1","[标题级] 关节炎 尿失禁","arthritis[ti] AND (incontinence[ti] OR incontinent[ti])",""),
 ("S2","[标题级] 关节痛 尿失禁","(\"joint pain\"[ti] OR arthralgia[ti] OR musculoskeletal[ti]) AND (incontinence[ti])",""),
 ("S3","[标题级] 关节炎 互联网/数字","arthritis[ti] AND (internet[ti] OR digital[ti] OR telehealth[ti] OR \"e-health\"[ti])",""),
 ("S4","[标题级] 关节炎 疫苗","arthritis[ti] AND (vaccin*[ti] OR immunis*[ti] OR immuniz*[ti])",""),
 ("S5","[方法先例] 级联 慢病 老龄库","(\"care cascade\"[tiab] OR \"treatment cascade\"[tiab] OR \"cascade of care\"[tiab]) AND (diabetes[tiab] OR hypertension[tiab] OR \"chronic disease\"[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab])",""),
 ("S6","[方法先例] 未诊断 慢病 老龄库","(undiagnosed[tiab]) AND (diabetes[tiab] OR hypertension[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab])",""),
 ("S7","[标题级] 未诊断 风湿 人群","(undiagnosed[ti] OR \"missed diagnosis\"[ti]) AND (rheumatic[ti] OR arthritis[ti] OR \"musculoskeletal\"[ti]) AND (population[tiab] OR cohort[tiab] OR survey[tiab])",""),
 ("S8","[标题级] 关节炎 疼痛 未治疗","arthritis[ti] AND (undertreated[ti] OR untreated[ti] OR \"unmet need\"[ti])",""),
 ("S9","[备选] 关节炎 失能 轨迹 多库","arthritis[tiab] AND (\"activities of daily living\"[tiab] OR disability[tiab]) AND trajectory[tiab] AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab])",""),
 ("S10","[备选] 关节炎 内在能力","arthritis[tiab] AND (\"intrinsic capacity\"[tiab])",""),
 ("S11","[备选] 关节炎 多病 共病 聚集 多库","arthritis[ti] AND (multimorbidity[ti] OR \"comorbidity\"[ti]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab])",""),
 ("S12","[占用校验] 关节炎 自付 已发表","(arthritis[tiab]) AND (\"out-of-pocket\"[tiab]) AND (\"catastrophic\"[tiab] OR \"financial burden\"[tiab] OR \"financial toxicity\"[tiab])",""),
 ("S13","[备选] 关节炎 社会参与 限制","arthritis[tiab] AND (\"social participation\"[tiab] OR \"social activity\"[tiab]) AND (restrict*[tiab] OR limit*[tiab] OR reduc*[tiab])",""),
 ("S14","[备选] 关节炎 照护 接受","arthritis[tiab] AND (\"received help\"[tiab] OR \"unmet need for care\"[tiab] OR \"long-term care\"[tiab]) AND (\"older adults\"[tiab] OR cohort[tiab])",""),
 ("S15","[备选] 关节炎 跌倒恐惧/活动限制","arthritis[tiab] AND (\"fear of falling\"[tiab] OR \"activity restriction\"[tiab])",""),
]

def esearch(term, retmax=8):
    url = BASE + "esearch.fcgi?" + urllib.parse.urlencode({"db":"pubmed","term":term,"retmode":"json","retmax":retmax,"sort":"relevance"})
    last=None
    for a in range(3):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.loads(r.read().decode())["esearchresult"]
        except Exception as e:
            last=e; time.sleep(2+2*a)
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

with open(os.path.join(OUT,"06_pubmed_r4.txt"),"w",encoding="utf-8") as f:
    f.write("PubMed novelty scan round 4 (title-level + method precedent) | 2026-11-14\n"+"\n".join(lines))
with open(os.path.join(OUT,"06_pubmed_r4.json"),"w",encoding="utf-8") as f:
    json.dump(js,f,ensure_ascii=False,indent=1)
print("WROTE",len(js))
