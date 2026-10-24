"""Cross-source mechanistic bridge discovery for V3.1.

Connects:
  primary-gene defect --(source A: target perturbation rescues defect)--> target
  approved medicine --(source B: direct modulation)------------------> target

No patient variants/phenotype are transmitted. Copyrighted abstracts are
processed transiently and are NOT persisted in output; only citation metadata,
hashes and structured paraphrases are retained.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,html,json,os,re,urllib.parse,urllib.request
from v2810.core import norm_name
from . import openai_reasoner

SCHEMA="RDMVA_CROSS_SOURCE_BRIDGE_V31_1"
UA="RareDiseaseMVAPharma-V31/1.0"
EPMC="https://www.ebi.ac.uk/europepmc/webservices/rest/search"

def _clean(s):
    s=re.sub(r"<[^>]+>"," ",html.unescape(str(s or "")))
    return " ".join(s.split())

def _search(query,page_size=100):
    u=EPMC+"?"+urllib.parse.urlencode({
      "query":query,"format":"json","pageSize":int(page_size),"resultType":"core"})
    req=urllib.request.Request(u,headers={"User-Agent":"RareDiseaseMVAPharma-V31/1.0"})
    with urllib.request.urlopen(req,timeout=35) as r:
        raw=json.loads(r.read(8_000_000))
    return (raw.get("resultList") or {}).get("result") or []

def _record(row):
    abstract=_clean(row.get("abstractText"))
    title=_clean(row.get("title"))
    return {
      "pmid":row.get("pmid"),"pmcid":row.get("pmcid"),
      "title":title,"abstract":abstract,
      "abstract_sha256":hashlib.sha256(abstract.encode()).hexdigest() if abstract else None,
      "is_open_access":str(row.get("isOpenAccess") or "").upper()=="Y",
      "journal":row.get("journalTitle"),"year":row.get("pubYear"),
      "license":row.get("license") or "SOURCE_SPECIFIC_REFERENCE_ONLY"}

def _name_map(entities):
    out={}
    for e in entities or []:
        if not bool(e.get("global_approved_as_single") or e.get("solo_approved")
                    or e.get("biologic_licensed") or not e.get("combo_only",False)):
            continue
        for n in [e.get("canonical_label"),*(e.get("member_labels") or [])]:
            n=str(n or "").strip()
            if len(n)>=6:out.setdefault(n.lower(),e)
    return out

def _term(text,term):
    return bool(re.search(r"(?<![A-Za-z0-9])"+re.escape(str(term))+
                          r"(?![A-Za-z0-9])",str(text),re.I))

def _objective(primary_gene,objectives):
    opts=[o for o in objectives.get("objectives",[])
          if o.get("gene")==primary_gene and o.get("drug_search_eligible")]
    for typ in ("bypass_or_compensate_downstream","compensate_downstream_deficit",
                "restore_complex_or_interaction_function","increase_residual_activity"):
        x=next((o for o in opts if o.get("objective_type")==typ),None)
        if x:return x
    return next(iter(opts),None)

def _source_payload(records,include_abstract=True):
    rows=[]
    for i,r in enumerate(records):
        rows.append({
          "index":i,"pmid":r.get("pmid"),"pmcid":r.get("pmcid"),
          "title":r.get("title"),"abstract":r.get("abstract") if include_abstract else None})
    return rows

def _contains_any(text,terms):
    lo=str(text or "").lower()
    return any(str(x).lower() in lo for x in terms if x)

def _rescue_score(r,primary_gene,aliases,primary_gene_aliases=None,case_terms=None):
    text=(r.get("title","")+" "+r.get("abstract","")).lower()
    pgs=[str(primary_gene or "").lower()]+[
      str(x or "").lower() for x in (primary_gene_aliases or []) if x]
    pgs=list(dict.fromkeys(x for x in pgs if x))
    score=0
    if any(x in text for x in pgs):score+=7
    if any(str(x).lower() in text for x in aliases):score+=7
    for k in ("rescu","restore","suppress","compensat","deficien","dysfunction",
              "phenotype","functional","mutant","loss of function","checkpoint","pathway"):
        if k in text:score+=2
    for k in (case_terms or []):
        if str(k).lower() in text:score+=1
    if any(x in r.get("title","").lower() for x in pgs):score+=5
    return score

def _drug_doc_score(r,aliases,candidates):
    text=(r.get("title","")+" "+r.get("abstract","")).lower()
    score=0
    if any(str(x).lower() in text for x in aliases):score+=7
    if "fda-approved" in text or "fda approved" in text:score+=8
    if "approved" in text:score+=3
    for k in ("bind","binding","kinase","inhibit","activity","target","off-target","repurpos"):
        if k in text:score+=2
    if any(str(x.get("candidate") or "").lower() in r.get("title","").lower()
           for x in candidates):score+=6
    return score

def _model_name():
    cfg=openai_reasoner.config()
    return (os.environ.get(cfg.get("environment_literature_model_override",
                                   "RDMVA_OPENAI_LITERATURE_MODEL"))
            or cfg.get("literature_extraction_model") or cfg["default_model"])

def _adjudicate_rescue(primary_gene,target_symbol,aliases,desired_action,records,out,budget_root=None):
    schema={"type":"object","additionalProperties":False,"properties":{
      "reviews":{"type":"array","items":{
        "type":"object","additionalProperties":False,
        "properties":{
          "record_index":{"type":"integer"},"supported":{"type":"boolean"},
          "target_action":{"type":"string","enum":["increase","decrease","unknown"]},
          "target_directly_perturbed":{"type":"boolean"},
          "primary_gene_defect_rescue":{"type":"boolean"},
          "primary_gene_or_disease_context":{"type":"boolean"},
          "reason":{"type":"string","maxLength":600}},
        "required":["record_index","supported","target_action","target_directly_perturbed",
                    "primary_gene_defect_rescue","primary_gene_or_disease_context","reason"]}}
    },"required":["reviews"]}
    payload={
      "primary_gene":primary_gene,"target_symbol":target_symbol,
      "target_aliases":aliases,"required_target_action":desired_action,
      "sources":_source_payload(records),
      "patient_variant_or_phenotype_data_included":False}
    instructions=(
      "Use ONLY the supplied publication titles/abstracts. Accept a source only if the NAMED TARGET ITSELF "
      "is experimentally perturbed in the required direction (for example inhibition/knockdown for decrease or "
      "activation/overexpression for increase) and that perturbation rescues/restores a molecular, cellular, pathway, "
      "physiologic or phenotypic defect attributable to the primary gene or disease context. "
      "Engineering the primary gene to bind the target better, observing target activation downstream of another insult, "
      "pathway membership, correlation, cancer association, or generic target biology is insufficient. "
      "Set target_directly_perturbed=false unless the target itself was manipulated. "
      "If the required direction is not directly supported, return supported=false.")
    body={"model":_model_name(),"input":json.dumps(payload,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],"max_output_tokens":3500,
          "text":{"format":{"type":"json_schema","name":"v31_bridge_rescue",
                           "strict":True,"schema":schema}}}
    return openai_reasoner.request_json(body,out/"rescue_ai.json",budget_root=budget_root or out)

def _adjudicate_drugs(target_symbol,aliases,desired_action,records,candidates,out,budget_root=None):
    schema={"type":"object","additionalProperties":False,"properties":{
      "reviews":{"type":"array","items":{
        "type":"object","additionalProperties":False,
        "properties":{
          "candidate_index":{"type":"integer"},"supported":{"type":"boolean"},
          "target_action":{"type":"string","enum":["increase","decrease","unknown"]},
          "direct_binding_or_activity_assay":{"type":"boolean"},
          "relationship_class":{"type":"string","enum":[
            "DIRECT_BINDING","BIOCHEMICAL_TARGET_ASSAY","INDIRECT_CELLULAR_RESPONSE","UNKNOWN"]},
          "reason":{"type":"string","maxLength":600}},
        "required":["candidate_index","supported","target_action",
                    "direct_binding_or_activity_assay","relationship_class","reason"]}}
    },"required":["reviews"]}
    source_rows=[]
    for i,r in enumerate(records):
        source_rows.append({"record_index":i,"pmid":r.get("pmid"),"pmcid":r.get("pmcid"),
                            "title":r.get("title"),"abstract":r.get("abstract")})
    payload={
      "target_symbol":target_symbol,"target_aliases":aliases,
      "required_target_action":desired_action,
      "sources":source_rows,
      "candidates":[{**x,"candidate_index":i} for i,x in enumerate(candidates)],
      "patient_variant_or_phenotype_data_included":False}
    instructions=(
      "Use ONLY the supplied publication titles/abstracts. Accept a candidate only if the approved medicine "
      "itself DIRECTLY BINDS the named target or a BIOCHEMICAL/TARGET-SPECIFIC activity assay demonstrates modulation "
      "in the required direction. A cellular increase/decrease in target phosphorylation/expression that is merely "
      "downstream of microtubule damage, cell-cycle arrest, stress, toxicity, or another primary target is INDIRECT and "
      "must be rejected. Computational docking alone, class effects, correlation, generic signaling activation, or a "
      "different compound are insufficient. Do not infer clinical efficacy.")
    body={"model":_model_name(),"input":json.dumps(payload,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],"max_output_tokens":3500,
          "text":{"format":{"type":"json_schema","name":"v31_bridge_drug",
                           "strict":True,"schema":schema}}}
    return openai_reasoner.request_json(body,out/"drug_ai.json",budget_root=budget_root or out)

def discover(primary_gene,target_symbol,target_aliases,desired_action,
             entities,objectives,outdir,max_results=100,budget_root=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    aliases=[target_symbol,*list(target_aliases or [])]
    seen=set();aliases=[str(x) for x in aliases if x and
                       not (str(x).lower() in seen or seen.add(str(x).lower()))]
    tq=" OR ".join(f'"{x}"' if " " in x else x for x in aliases)
    pgq=f"({primary_gene} OR BUBR1)" if primary_gene.upper()=="BUB1B" else primary_gene

    rescue_queries=[
      f'(({pgq}) AND ({tq}) AND (rescue OR restore OR suppressor OR compensation OR congression OR checkpoint OR MVA))',
      f'(({pgq}) AND ({tq}) AND (PP2A OR B56 OR congression OR aneuploid*))',
      f'(({pgq}) AND ({tq}) AND (rescue OR restore OR congression)) AND FIRST_PDATE:[2000-01-01 TO 2020-12-31]']
    drug_queries=[
      f'(TITLE_ABS:({tq}) AND "FDA-approved")',
      f'(TITLE_ABS:({tq}) AND (approved OR repurpos* OR reposition* OR "off-target"))',
      f'(({tq}) AND ("FDA approved" OR "FDA-approved"))',
      f'(({tq}) AND (approved OR repurpos* OR reposition* OR "off-target") AND (drug OR inhibitor OR activator OR binding))']

    def collect(queries):
        found=[];seen_ids=set()
        for q in queries:
            try:rows=_search(q,max_results)
            except Exception:continue
            for row in rows:
                r=_record(row)
                sid=r.get("pmcid") or ("PMID:"+str(r.get("pmid")))
                if not r.get("abstract") or sid in seen_ids:continue
                seen_ids.add(sid);r["search_query"]=q;found.append(r)
        return found

    rescue_all=collect(rescue_queries)
    drug_all=collect(drug_queries)

    # Deterministic compression before AI. Important: search 100+ records for
    # recall, but send only biologically dense records to the model.
    rescue=[r for r in rescue_all if _rescue_score(r,primary_gene,aliases)>=14]
    rescue=sorted(rescue,key=lambda r:(-_rescue_score(r,primary_gene,aliases),
                                       int(r.get("year") or 9999),
                                       str(r.get("pmid") or "")))[:8]

    amap=_name_map(entities)
    raw_candidates=[]
    for ri,r in enumerate(drug_all):
        text=r["title"]+" "+r["abstract"];lo=text.lower()
        if not _contains_any(text,aliases):continue
        for alias,e in amap.items():
            if alias in lo and _term(text,alias):
                raw_candidates.append({
                  "record_index":ri,"alias":alias,
                  "candidate":e.get("canonical_label"),
                  "entity_id":e.get("entity_id"),
                  "scientific_group_id":e.get("scientific_group_id"),
                  "member_labels":e.get("member_labels") or [e.get("canonical_label")]})
    dd={}
    for x in raw_candidates:
        dd[(x["record_index"],x["candidate"])]=x
    raw_candidates=list(dd.values())

    # Rank documents with at least one approved medicine mention, then remap
    # record indices into the compact AI payload.
    by_record={}
    for x in raw_candidates:by_record.setdefault(x["record_index"],[]).append(x)
    ranked=[]
    for ri,cands in by_record.items():
        r=drug_all[ri]
        ranked.append((_drug_doc_score(r,aliases,cands),ri,r,cands))
    ranked=sorted(ranked,key=lambda z:(-z[0],int(z[2].get("year") or 9999),
                                       str(z[2].get("pmid") or "")))[:8]
    drug_records=[];drug_candidates=[]
    for new_ri,(score,old_ri,r,cands) in enumerate(ranked):
        drug_records.append(r)
        for x in cands:
            y=dict(x);y["record_index"]=new_ri
            drug_candidates.append(y)
    # Limit candidate rows rather than repeatedly sending salt/alias variants.
    seen_c=set();compact=[]
    for x in drug_candidates:
        key=(x["record_index"],norm_name(x["candidate"]))
        if key in seen_c:continue
        seen_c.add(key);compact.append(x)
    drug_candidates=compact[:16]

    rescue_rec=_adjudicate_rescue(
        primary_gene,target_symbol,aliases,desired_action,rescue,out,
        budget_root=budget_root) if rescue else {"value":{"reviews":[]}}
    rv={int(x["record_index"]):x for x in rescue_rec.get("value",{}).get("reviews",[])
        if isinstance(x.get("record_index"),int)}
    drug_rec=_adjudicate_drugs(
        target_symbol,aliases,desired_action,drug_records,drug_candidates,out,
        budget_root=budget_root) if drug_candidates else {"value":{"reviews":[]}}
    dv={int(x["candidate_index"]):x for x in drug_rec.get("value",{}).get("reviews",[])
        if isinstance(x.get("candidate_index"),int)}

    rescue_ok=[]
    for i,r in enumerate(rescue):
        x=rv.get(i) or {}
        if (x.get("supported") and x.get("target_directly_perturbed") and
            x.get("primary_gene_defect_rescue") and
            x.get("target_action")==desired_action):
            rescue_ok.append({
              "pmid":r.get("pmid"),"pmcid":r.get("pmcid"),"title":r.get("title"),
              "abstract_sha256":r.get("abstract_sha256"),"license":r.get("license"),
              "text_redistributed":False,"reason":x.get("reason")})
    drug_ok=[]
    for ci,x in enumerate(drug_candidates):
        v=dv.get(ci) or {}
        if not (v.get("supported") and v.get("direct_binding_or_activity_assay")
                and v.get("relationship_class") in ("DIRECT_BINDING","BIOCHEMICAL_TARGET_ASSAY")
                and v.get("target_action")==desired_action):continue
        r=drug_records[x["record_index"]]
        drug_ok.append({
          "candidate":x["candidate"],"entity_id":x["entity_id"],
          "scientific_group_id":x.get("scientific_group_id"),
          "member_labels":x.get("member_labels"),
          "pmid":r.get("pmid"),"pmcid":r.get("pmcid"),"title":r.get("title"),
          "abstract_sha256":r.get("abstract_sha256"),"license":r.get("license"),
          "text_redistributed":False,"reason":v.get("reason")})

    # De-duplicate drug/source claims.
    uniq={}
    for d in drug_ok:
        uniq[(norm_name(d["candidate"]),d.get("pmid") or d.get("pmcid"))]=d
    drug_ok=list(uniq.values())

    obj=_objective(primary_gene,objectives)
    bridges=[]
    if obj and rescue_ok:
        for d in drug_ok:
            bridges.append({
              "candidate":d["candidate"],"entity_id":d["entity_id"],
              "scientific_group_id":d.get("scientific_group_id"),
              "member_labels":d.get("member_labels"),
              "primary_gene":primary_gene,"target_symbol":target_symbol,
              "target_aliases":aliases,"desired_action":desired_action,
              "objective_id":obj.get("objective_id"),
              "objective_type":obj.get("objective_type"),
              "objective_priority":obj.get("priority_tier"),
              "rescue_sources":rescue_ok,"drug_target_source":d,
              "bridge_class":"CROSS_SOURCE_MECHANISTIC_RESCUE",
              "same_drug_disease_rescue_demonstrated":False,
              "clinical_recommendation":False})
    result={
      "schema":SCHEMA,"primary_gene":primary_gene,"target_symbol":target_symbol,
      "target_aliases":aliases,"desired_action":desired_action,
      "rescue_search_records":len(rescue_all),"rescue_records_sent_to_ai":len(rescue),
      "drug_search_records":len(drug_all),"drug_records_sent_to_ai":len(drug_records),
      "drug_candidate_rows_sent_to_ai":len(drug_candidates),
      "rescue_source_count":len(rescue_ok),"drug_target_candidate_count":len(drug_ok),
      "bridges":bridges,"bridge_count":len(bridges),
      "abstract_text_persisted":False,"non_cc_text_redistributed":False,
      "patient_data_sent_external":False,"clinical_recommendation":False}
    (out/"CROSS_SOURCE_BRIDGE_V31.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return result

def merge_pairs(pairs,result):
    out=list(pairs or []);routes=[];added=0
    for i,b in enumerate(result.get("bridges",[]),1):
        rid=f"V31XSB{i:04d}"
        route={
          "route_id":rid,"rank":0,"family":"v31_cross_source_bridge",
          "case_gene":b.get("primary_gene"),"target_concept":b.get("target_symbol"),
          "target_aliases":list(b.get("target_aliases") or []),
          "targets":[],"desired_action":b.get("desired_action"),"variable":"activity",
          "mechanism_concept":"approved medicine directly modulates a rescue target supported by an independent primary-gene/MVA rescue source",
          "score":17.0,"state":"conditional","network_hops":1,
          "rescue_objective_id":b.get("objective_id"),
          "rescue_objective_ids":[b.get("objective_id")],
          "rescue_objective_type":b.get("objective_type"),
          "rescue_objective_types":[b.get("objective_type")],
          "rescue_priority_tier":b.get("objective_priority"),
          "source_bound_cross_source_bridge":True,
          "rescue_sources":b.get("rescue_sources"),
          "drug_target_source":b.get("drug_target_source"),
          "scientific_gate":{
            "schema":"V31_CROSS_SOURCE_BRIDGE_GATE_1",
            "primary_gene_target_identity":False,
            "allele_specific_rescue_verified":False,
            "eligible_for_exploratory_screen":True,
            "functional_rescue_validated":False,
            "causal_direction_validated_for_disease":True,
            "requires_independent_endpoint":True,
            "hard_hold":False,
            "issues":["same_drug_disease_model_rescue_not_demonstrated",
                      "bridge_requires_patient_variant_validation"]}}
        routes.append(route)
        out.append({
          "schema":"RDMVA_CROSS_SOURCE_BRIDGE_PAIR_V31_1",
          "ingredient":b.get("candidate"),"entity_id":b.get("entity_id"),
          "scientific_group_id":b.get("scientific_group_id"),
          "member_labels":b.get("member_labels"),
          "route_id":rid,"route":route,"tier":"C","score":14.0,
          "providers":["cross_source_mechanistic_bridge"],
          "signals":[{
            "route_id":rid,"provider":"cross_source_mechanistic_bridge",
            "direct_target_assay_verified":True,"functional_rescue":False,
            "target_action":b.get("desired_action"),
            "drug_target_source":b.get("drug_target_source"),
            "rescue_sources":b.get("rescue_sources")}],
          "combo_only":False,"global_approved_as_single":True,
          "publication_strength":"independent_target_rescue_plus_direct_drug_target_assay",
          "evidence_blockers":["same_drug_disease_model_rescue_not_demonstrated",
                               "patient_variant_specific_rescue_not_demonstrated"],
          "shortest_network_hops":1,"drug_to_route_hops":0,"route_to_case_hops":1,
          "source_bound_cross_source_bridge":True,
          "scientific_gate":route["scientific_gate"],
          "promotion_eligible_validated_lead":False,
          "clinical_recommendation":False})
        added+=1
    return out,{"added_pairs":added,"routes":routes}


def _route_target_candidates(routes,max_targets=60):
    pri={"O1":3,"O2":2,"O3":1,None:0}
    best={}
    for r in routes or []:
        if r.get("state")=="hold":continue
        targets=r.get("targets") or []
        if len(targets)!=1:continue
        tg=str(targets[0].get("gene") or "").upper().strip()
        pg=str(r.get("case_gene") or "").upper().strip()
        if not tg or not pg or tg==pg:continue
        if not re.fullmatch(r"[A-Z][A-Z0-9-]{1,19}",tg):continue
        action=r.get("desired_action")
        if action not in ("increase","decrease"):continue
        hops=int(r.get("network_hops") or 0)
        score=float(r.get("score") or 0)
        score+=3*pri.get(r.get("rescue_priority_tier"),1)
        score-=1.5*hops
        if r.get("family")=="v31_source_bound_modifier":score+=6
        if (r.get("scientific_gate") or {}).get("causal_direction_validated_for_disease"):
            score+=6
        key=(pg,tg,action)
        rec={
          "primary_gene":pg,"target_symbol":tg,
          "target_aliases":[tg,str(r.get("target_concept") or tg)],
          "desired_action":action,"route_id":r.get("route_id"),
          "route_family":r.get("family"),"route_score":score,
          "rescue_objective_ids":list(r.get("rescue_objective_ids") or []),
          "network_hops":hops}
        if key not in best or score>best[key]["route_score"]:best[key]=rec
    rows=sorted(best.values(),key=lambda x:(
        -x["route_score"],x["network_hops"],x["primary_gene"],
        x["target_symbol"],x["desired_action"]))
    return rows[:max_targets]

def _scout(candidate,entities):
    pg=candidate["primary_gene"];tg=candidate["target_symbol"]
    aliases=list(dict.fromkeys([tg,*candidate.get("target_aliases",[])]))
    tq=" OR ".join(f'"{x}"' if " " in x else x for x in aliases)
    pgq=f"({pg} OR BUBR1)" if pg=="BUB1B" else pg
    rescue_q=f'(({pgq}) AND ({tq}) AND (rescue OR restore OR suppressor OR congression OR checkpoint OR MVA))'
    drug_q=f'(TITLE_ABS:({tq}) AND (approved OR repurpos* OR reposition* OR "off-target") AND (drug OR inhibitor OR activator OR binding))'
    try:rescue=[_record(x) for x in _search(rescue_q,35)]
    except Exception:rescue=[]
    try:drug=[_record(x) for x in _search(drug_q,35)]
    except Exception:drug=[]
    rescue=[r for r in rescue if r.get("abstract") and
            _rescue_score(r,pg,aliases)>=12]
    amap=_name_map(entities)
    named=0
    for r in drug:
        text=(r.get("title","")+" "+r.get("abstract","")).lower()
        if not _contains_any(text,aliases):continue
        if any(alias in text and _term(text,alias) for alias in amap):
            named+=1
            if named>=1:break
    return {
      **candidate,
      "rescue_scout_hits":len(rescue),
      "approved_drug_scout_hits":named,
      "eligible_for_full_bridge":bool(rescue and named)}

def discover_from_routes(routes,entities,objectives,outdir,
                         max_targets=60,max_full_bridges=12):
    """High-recall target scouting followed by sparse full adjudication.

    Public Europe PMC searches are cheap and carry no patient-level data.
    OpenAI is called only for targets with both rescue-like literature and at
    least one approved medicine mention.
    """
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    candidates=_route_target_candidates(routes,max_targets=max_targets)
    scouted=[];full=[];all_bridges=[];errors=[]
    for c0 in candidates:
        c=_enrich_candidate_aliases(c0)
        s=_scout(c,entities);scouted.append(s)
        if not s["eligible_for_full_bridge"]:continue
        if len(full)>=max_full_bridges:continue
        sub=out/f'{len(full)+1:02d}_{c["primary_gene"]}_{c["target_symbol"]}_{c["desired_action"]}'
        try:
            result=discover(
              c["primary_gene"],c["target_symbol"],c.get("target_aliases",[]),
              c["desired_action"],entities,objectives,sub,max_results=100)
            full.append({
              "primary_gene":c["primary_gene"],"target_symbol":c["target_symbol"],
              "desired_action":c["desired_action"],"route_score":c["route_score"],
              "bridge_count":result.get("bridge_count",0),
              "result_file":str((sub/"CROSS_SOURCE_BRIDGE_V31.json").relative_to(out))})
            all_bridges.extend(result.get("bridges") or [])
        except Exception as exc:
            errors.append({
              "primary_gene":c["primary_gene"],"target_symbol":c["target_symbol"],
              "desired_action":c["desired_action"],
              "error":type(exc).__name__+":"+str(exc)[:500]})
    # De-duplicate identical candidate/target/evidence bridge.
    uniq={}
    for b in all_bridges:
        d=b.get("drug_target_source") or {}
        key=(norm_name(b.get("candidate")),b.get("primary_gene"),
             b.get("target_symbol"),b.get("desired_action"),
             d.get("pmid") or d.get("pmcid"))
        uniq[key]=b
    all_bridges=list(uniq.values())
    summary={
      "schema":"RDMVA_CROSS_SOURCE_BRIDGE_BATCH_V31_1",
      "targets_considered":len(candidates),
      "targets_scouted":len(scouted),
      "targets_full_adjudicated":len(full),
      "scout":scouted,"full_adjudications":full,
      "bridges":all_bridges,"bridge_count":len(all_bridges),
      "errors":errors,
      "patient_data_sent_external":False,
      "abstract_text_persisted_in_summary":False,
      "non_cc_text_redistributed":False,
      "clinical_recommendation":False}
    (out/"CROSS_SOURCE_BRIDGE_BATCH_V31.json").write_text(
      json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return summary


def _name_variants(name):
    s=" ".join(str(name or "").split()).strip()
    if not s:return []
    vals=[s]
    # Literature often drops the generic word "kinase" from established
    # family names (e.g. Aurora kinase B -> Aurora B).
    m=re.match(r"^(.+?)\s+kinase\s+([A-Za-z0-9-]+)$",s,re.I)
    if m:
        vals.append((m.group(1)+" "+m.group(2)).strip())
    if "-" in s:
        vals.append(s.replace("-"," "))
    return [x for x in vals if len(x)>=3]

def _uniprot_target_aliases(accession,gene):
    vals=[str(gene or "").upper()]
    if not accession:return vals
    try:
        url=f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(str(accession),safe='')}.json"
        req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
        with urllib.request.urlopen(req,timeout=25) as r:
            j=json.loads(r.read(2_000_000))
        pd=j.get("proteinDescription") or {}
        rn=((pd.get("recommendedName") or {}).get("fullName") or {}).get("value")
        if rn:vals.extend(_name_variants(rn))
        for a in pd.get("alternativeNames") or []:
            fn=(a.get("fullName") or {}).get("value")
            if fn:vals.extend(_name_variants(fn))
            for sn in a.get("shortNames") or []:
                v=sn.get("value")
                if v:vals.append(v)
        for g in j.get("genes") or []:
            gn=(g.get("geneName") or {}).get("value")
            if gn:vals.append(gn)
            for s in g.get("synonyms") or []:
                v=s.get("value")
                if v:vals.append(v)
    except Exception:
        pass
    clean=[];seen=set()
    for x in vals:
        x=" ".join(str(x or "").split()).strip()
        if len(x)<3:continue
        k=x.lower()
        if k in seen:continue
        seen.add(k);clean.append(x)
        if len(clean)>=16:break
    return clean

def _enrich_candidate_aliases(candidate):
    c=dict(candidate)
    aliases=list(c.get("target_aliases") or [])
    aliases.extend(_uniprot_target_aliases(c.get("target_accession"),
                                            c.get("target_symbol")))
    c["target_aliases"]=list(dict.fromkeys(x for x in aliases if x))
    return c
