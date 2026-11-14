# -*- coding: utf-8 -*-
"""PubMed novelty scan round 3: decisive queries for the shortlisted directions."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
import urllib.request, urllib.parse, json, time, os

OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERIES = [
 # 主推：诊疗级联 / 未诊断-未治疗
 ("R1","[决定性] 风湿诊疗级联 0命中检验",
  "(arthritis[tiab] OR \"joint pain\"[tiab] OR arthralgia[tiab] OR rheumat*[tiab]) AND (\"care cascade\"[tiab] OR \"treatment cascade\"[tiab] OR \"diagnosis cascade\"[tiab] OR \"diagnostic cascade\"[tiab])",""),
 ("R2","[决定性] 未诊断关节痛 多库",
  "(\"joint pain\"[tiab] OR arthralgia[tiab]) AND (undiagnosed[tiab] OR \"never diagnosed\"[tiab] OR \"no diagnosis\"[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab])",""),
 ("R3","[决定性] 疼痛-诊断-治疗 序贯",
  "(pain[tiab]) AND (diagnos*[tiab]) AND (treatment[tiab] OR medication[tiab]) AND (\"older adults\"[tiab]) AND (\"care gap\"[tiab] OR gap[tiab] OR cascade[tiab]) AND (arthritis[tiab] OR rheumat*[tiab] OR musculoskeletal[tiab])",""),
 ("R4","关节炎 未治疗 症状性",
  "arthritis[tiab] AND (\"symptomatic\"[tiab]) AND (\"no treatment\"[tiab] OR untreated[tiab] OR \"not receiving\"[tiab]) AND (\"older adults\"[tiab] OR cohort[tiab])",""),
 ("R5","关节炎 用药 跨国",
  "(arthritis[tiab] OR rheumat*[tiab]) AND (\"medication use\"[tiab] OR \"drug utilization\"[tiab] OR \"treatment utilisation\"[tiab]) AND (\"cross-national\"[tiab] OR \"cross-country\"[tiab] OR \"across countries\"[tiab] OR international[tiab])",""),
 # 备选：失能恢复/转变
 ("R6","关节炎 失能 恢复",
  "arthritis[tiab] AND (disability[tiab] OR \"functional limitation\"[tiab]) AND (recovery[tiab] OR transition*[tiab] OR revers*[tiab]) AND (\"older adults\"[tiab] OR cohort[tiab])",""),
 ("R7","多状态模型 关节炎",
  "(multistate[tiab] OR \"multi-state\"[tiab] OR \"Markov\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab] OR musculoskeletal[tiab])",""),
 # 备选：预防服务
 ("R8","关节炎 预防服务/筛查",
  "arthritis[tiab] AND (\"preventive services\"[tiab] OR \"cancer screening\"[tiab] OR \"mammography\"[tiab] OR \"vaccination\"[tiab]) AND (\"older adults\"[tiab] OR population[tiab]) AND (\"receipt\"[tiab] OR uptake[tiab] OR utilization[tiab])",""),
 ("R9","关节炎 流感/肺炎疫苗 多库",
  "(arthritis[tiab] OR rheumat*[tiab]) AND (\"influenza\"[tiab] OR \"pneumococcal\"[tiab] OR \"shingles\"[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab])",""),
 # 备选：尿失禁
 ("R10","尿失禁 关节炎 老年人 横断",
  "(incontinence[tiab]) AND (arthritis[tiab] OR \"joint pain\"[tiab] OR \"mobility limitation\"[tiab]) AND (\"older adults\"[tiab] OR women[tiab])",""),
 # 备选：口腔
 ("R11","牙齿脱落 关节炎 多库 结局",
  "(\"tooth loss\"[tiab] OR edentul*[tiab] OR \"number of teeth\"[tiab]) AND (arthritis[tiab]) AND (CHARLS[tiab] OR SHARE[tiab] OR HRS[tiab] OR LASI[tiab] OR ELSA[tiab])",""),
 # 备选：数字
 ("R12","数字排斥 慢病 老年 多库",
  "(\"digital divide\"[tiab] OR \"digital exclusion\"[tiab] OR \"internet non-use\"[tiab]) AND (\"older adults\"[tiab] OR aged[tiab]) AND (chronic[tiab] OR multimorbidity[tiab])",""),
 # 关节置换
 ("R13","关节置换 跨国 不平等",
  "(\"joint replacement\"[tiab] OR \"hip replacement\"[tiab] OR \"knee replacement\"[tiab]) AND (inequalit*[tiab] OR inequit*[tiab] OR disparity[tiab]) AND (\"cross-national\"[tiab] OR \"across countries\"[tiab] OR international[tiab])",""),
 # 其他风湿-老龄化交叉，验证是否已占位
 ("R14","痛风 老龄化队列",
  "gout[tiab] AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab])",""),
 ("R15","高尿酸 无痛风 结局",
  "(\"asymptomatic hyperuricemia\"[tiab] OR \"asymptomatic hyperuricaemia\"[tiab]) AND (\"older adults\"[tiab] OR cohort[tiab]) AND (mortality[tiab] OR disability[tiab] OR \"joint pain\"[tiab])",""),
 ("R16","骨质疏松 治疗缺口 多库",
  "(osteoporosis[tiab]) AND (\"treatment gap\"[tiab] OR undertreat*[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab] OR ELSA[tiab])",""),
 ("R17","关节炎 多重共病 轨迹 多库",
  "arthritis[tiab] AND (multimorbidity[tiab] OR \"comorbidity pattern\"[tiab] OR trajectory[tiab]) AND (CHARLS[tiab] OR HRS[tiab] OR SHARE[tiab] OR LASI[tiab] OR MHAS[tiab])",""),
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

with open(os.path.join(OUT,"04_pubmed_r3.txt"),"w",encoding="utf-8") as f:
    f.write("PubMed novelty scan round 3 (decisive) | 2026-11-14\n"+"\n".join(lines))
with open(os.path.join(OUT,"04_pubmed_r3.json"),"w",encoding="utf-8") as f:
    json.dump(js,f,ensure_ascii=False,indent=1)
print("WROTE",len(js))
