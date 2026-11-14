import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("AGING_DB_RAW_ROOT", os.path.join(ROOT, "data_raw"))
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
DL = os.environ.get("AGING_DB_ROOT", os.path.join(ROOT, "data"))
# -*- coding: utf-8 -*-
import urllib.request, urllib.parse, json, time, os
OUT = os.environ.get("OOPE02_OUT", os.path.join(ROOT, "outputs"))
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
Q=[
 ("T1","传统医学×风湿 多库","(\"traditional medicine\"[tiab] OR \"folk healer\"[tiab] OR \"herbal medicine\"[tiab] OR \"complementary medicine\"[tiab]) AND (arthritis[tiab] OR \"joint pain\"[tiab] OR rheumat*[tiab]) AND (CHARLS[tiab] OR KLoSA[tiab] OR LASI[tiab] OR MHAS[tiab] OR \"cross-national\"[tiab])"),
 ("T2","传统医学 关节炎 老年 人群","(\"traditional medicine\"[tiab] OR \"complementary and alternative\"[tiab]) AND (arthritis[tiab] OR \"musculoskeletal pain\"[tiab]) AND (\"older adults\"[tiab] OR prevalence[tiab])"),
 ("T3","关节炎 失能 多状态 转归","(arthritis[tiab] OR rheumat*[tiab]) AND (\"disability transition\"[tiab] OR \"transitions in disability\"[tiab] OR \"multistate model\"[tiab] OR \"multi-state model\"[tiab])"),
 ("T4","失能 恢复 关节炎 老年","(\"recovery from disability\"[tiab] OR \"disability recovery\"[tiab]) AND (arthritis[tiab] OR rheumat*[tiab] OR musculoskeletal[tiab])"),
]
def es(t):
    u=BASE+"esearch.fcgi?"+urllib.parse.urlencode({"db":"pubmed","term":t,"retmode":"json","retmax":8,"sort":"relevance"})
    with urllib.request.urlopen(u,timeout=40) as r: return json.loads(r.read().decode())["esearchresult"]
def su(ids):
    if not ids: return {}
    u=BASE+"esummary.fcgi?"+urllib.parse.urlencode({"db":"pubmed","id":",".join(ids),"retmode":"json"})
    with urllib.request.urlopen(u,timeout=40) as r: return json.loads(r.read().decode()).get("result",{})
L=["PubMed round 5 | 2026-11-14"]
for qid,note,t in Q:
    r=es(t); ids=r.get("idlist",[]); s=su(ids)
    L+= [f"\n[{qid}] {note}\n  hits={r.get('count')}\n  query={t}"]
    for p in ids:
        it=s.get(p,{}); L.append(f"    - {p} | {(it.get('pubdate','') or '')[:4]} | {it.get('source','')} | {it.get('title','')[:180]}")
    print(qid,r.get("count"),flush=True); time.sleep(0.4)
open(os.path.join(OUT,"09_pubmed_r5.txt"),"w",encoding="utf-8").write("\n".join(L))
