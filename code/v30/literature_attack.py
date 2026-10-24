"""Source-bound public literature attack for V3.0 research candidates.

Search uses Europe PMC open-access full text through the V2.8.10 retriever.
OpenAI only extracts/criticizes evidence from those public passages. Every
accepted quote is checked byte-for-byte against the retrieved passage.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import json,os,re
from v2810 import ai as v28_ai
from . import openai_reasoner

SCHEMA="RDMVA_SOURCE_BOUND_LITERATURE_ATTACK_V30_1"

def _norm(x):return re.sub(r"[^a-z0-9]+","",str(x or "").lower())

def select_candidates(falsification,graph,limit=8):
    pidx=defaultdict(list)
    for p in graph.get("candidate_paths",[]):
        pidx[p.get("candidate")].append(p)
    rows=[];seen_axes=set()
    ordered=[x for x in falsification.get("candidates",[]) if x.get("decision")=="HOLD"]
    for c in ordered:
        paths=pidx.get(c.get("candidate"),[])
        path=next((p for p in paths if p.get("route_id")==c.get("route_id")
                   and p.get("objective_id")==c.get("objective_id")), paths[0] if paths else None)
        if not path:continue
        axis=(path.get("gene"),path.get("objective_type"))
        if axis in seen_axes:continue
        seen_axes.add(axis);rows.append((c,path))
        if len(rows)>=limit:return rows
    for c in ordered:
        if any(x[0].get("candidate")==c.get("candidate") for x in rows):continue
        paths=pidx.get(c.get("candidate"),[])
        path=paths[0] if paths else None
        if path:rows.append((c,path))
        if len(rows)>=limit:break
    return rows

def _passage_payload(docs,candidate,path,max_passages=18):
    aliases=[x for x in path.get("member_labels",[]) if x]+[candidate]
    terms=[*aliases,path.get("gene"),path.get("route_target"),
           "rescue","restore","stabil","activity","function","inhibit","activate","degrad"]
    scored=[]
    for d in docs:
        for p in d.get("passages",[]):
            text=str(p.get("text") or "")
            lo=text.lower()
            score=sum(3 for t in terms if t and str(t).lower() in lo)
            score+=sum(1 for t in ("dose","concentration","rescue","restored","inhibition",
                                   "activation","knockdown","overexpression","degradation") if t in lo)
            scored.append((score,d,p))
    scored.sort(key=lambda x:(-x[0],str(x[1].get("source_id")),str(x[2].get("id"))))
    out=[];seen=set()
    for score,d,p in scored:
        key=(d.get("source_id"),p.get("id"))
        if key in seen:continue
        seen.add(key)
        out.append({"source_id":d.get("source_id"),"citation":d.get("citation"),
          "license":d.get("license"),"passage_id":p.get("id"),
          "text":str(p.get("text") or "")[:2200]})
        if len(out)>=max_passages:break
    return out

def _call_extractor(candidate,path,passages,outfile,budget_root):
    cfg=openai_reasoner.config()
    model=os.environ.get(cfg.get("environment_literature_model_override","RDMVA_OPENAI_LITERATURE_MODEL"),
                         cfg.get("literature_extraction_model") or cfg["default_model"])
    aliases=[x for x in path.get("member_labels",[]) if x] or [candidate]
    payload={
      "candidate":candidate,"candidate_aliases":aliases,
      "gene_hypothesis":path.get("gene"),"route_target":path.get("route_target"),
      "route_desired_action":path.get("route_desired_action"),
      "rescue_objective":path.get("objective_type"),
      "public_passages":passages,
      "patient_specific_variants_or_phenotype_included":False,
    }
    schema={"type":"object","additionalProperties":False,"properties":{
      "evidence_items":{"type":"array","items":{"type":"object","additionalProperties":False,
        "properties":{
          "source_id":{"type":"string"},"passage_id":{"type":"string"},
          "exact_quote":{"type":"string","maxLength":900},
          "evidence_type":{"type":"string","enum":[
             "DIRECT_TARGET_MODULATION","SAME_DRUG_FUNCTIONAL_RESCUE",
             "MECHANISTIC_COUNTEREVIDENCE","EXPOSURE_CONTEXT","OTHER"]},
          "direction":{"type":"string","enum":["increase","decrease","restore","unknown"]},
          "mechanistically_relevant":{"type":"boolean"},
          "reason":{"type":"string","maxLength":500}
        },
        "required":["source_id","passage_id","exact_quote","evidence_type",
                    "direction","mechanistically_relevant","reason"]}},
      "overall_assessment":{"type":"string","enum":[
         "DIRECT_RESCUE_EVIDENCE_FOUND","PARTIAL_MECHANISTIC_SUPPORT",
         "COUNTEREVIDENCE_DOMINANT","INSUFFICIENT"]},
      "missing_evidence":{"type":"array","items":{"type":"string"}}
    },"required":["evidence_items","overall_assessment","missing_evidence"]}
    instructions=(
      "Use only the supplied open-access public passages. Extract evidence for the named compound and the specified "
      "gene/route/rescue objective. Do not use outside knowledge. DIRECT_TARGET_MODULATION requires an experimental "
      "statement that the compound itself modulates the relevant target or route in the required direction. "
      "SAME_DRUG_FUNCTIONAL_RESCUE requires the same compound to improve/restore a function that is mechanistically "
      "relevant to the stated gene/objective; generic anticancer cytotoxicity, correlation, expression signatures, "
      "or pathway membership do not qualify. MECHANISTIC_COUNTEREVIDENCE includes opposite-direction effects or evidence "
      "that would worsen the stated objective. Copy every supporting quote exactly and contiguously from one supplied passage. "
      "If exact evidence is absent return INSUFFICIENT. This is research evidence extraction, not treatment advice."
    )
    body={"model":model,"input":json.dumps(payload,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],
          "max_output_tokens":min(2200,int(cfg.get("max_output_tokens",4500))),
          "text":{"format":{"type":"json_schema","name":"v30_source_bound_evidence",
                            "strict":True,"schema":schema}}}
    return openai_reasoner.request_json(body,outfile,budget_root=budget_root)


def _deterministic_counter_items(passages,path,candidate=None):
    """Catch explicit adverse route biology even if the LLM omits it.

    Only high-specificity sentences are admitted: the same sentence must name
    the primary gene and route target (or target stem), describe loss/defect of
    the gene, and describe resistance/escape/reactivation/worsening.
    """
    gene=str(path.get("gene") or "").strip()
    target=str(path.get("route_target") or "").strip()
    if not gene or not target:return []
    gene_rx=re.compile(r"\b"+re.escape(gene)+r"\b",re.I)
    # For concepts like "KRAS pathway", use the first alphanumeric token.
    target_token=(re.findall(r"[A-Za-z0-9_-]{3,}",target) or [target])[0]
    target_rx=re.compile(r"\b"+re.escape(target_token)+r"\b",re.I)
    defect_rx=re.compile(
      r"\b(?:loss|loss-of-function|inactivat(?:e|ed|ion)|deficien(?:t|cy)|"
      r"knock(?:out|down)|delet(?:e|ed|ion)|mutant|mutation)\b",re.I)
    adverse_rx=re.compile(
      r"\b(?:resistan(?:ce|t)|reactivat(?:e|ed|ion)|escape|rebound|"
      r"worsen(?:s|ed|ing)?|fails?\s+to\s+respond|refractor)\b",re.I)
    out=[];seen=set()
    for passage in passages:
        text=str(passage.get("text") or "")
        for sentence in [x.strip() for x in re.split(r"(?<=[.!?])\s+",text) if x.strip()]:
            if not (gene_rx.search(sentence) and target_rx.search(sentence)
                    and defect_rx.search(sentence) and adverse_rx.search(sentence)):
                continue
            # Direct source-bound drug-rescue routes already have compound-
            # specific positive evidence. Generic resistance statements about
            # the disease gene are not counterevidence to that compound.
            # The LLM extractor may still admit class-level counterevidence
            # when the passage actually identifies the relevant drug class.
            if str(path.get("route_family") or "")=="v31_direct_source_drug_rescue":
                aliases=[str(candidate or ""),*(path.get("member_labels") or [])]
                if not any(a and re.search(r"\b"+re.escape(a)+r"\b",sentence,re.I)
                           for a in aliases):
                    continue
            key=(passage.get("source_id"),passage.get("passage_id"),sentence)
            if key in seen:continue
            seen.add(key)
            out.append({
              "source_id":passage.get("source_id"),
              "passage_id":passage.get("passage_id"),
              "exact_quote":sentence,
              "evidence_type":"MECHANISTIC_COUNTEREVIDENCE",
              "direction":"unknown",
              "mechanistically_relevant":True,
              "reason":"Deterministic high-specificity primary-gene/route resistance or reactivation sentence.",
              "deterministic_counter_scan":True,
            })
    return out

def _validate_items(candidate,path,passages,model_result):
    pmap={(x["source_id"],x["passage_id"]):x["text"] for x in passages}
    aliases=[_norm(x) for x in path.get("member_labels",[]) if _norm(x)]
    if _norm(candidate):aliases.append(_norm(candidate))
    gene=_norm(path.get("gene"));target=_norm(path.get("route_target"))
    accepted=[];rejected=[]
    for item in model_result.get("value",{}).get("evidence_items",[]):
        key=(item.get("source_id"),item.get("passage_id"));quote=item.get("exact_quote") or ""
        source=pmap.get(key)
        reason=None
        if not source:reason="unknown_source_or_passage"
        elif not quote or quote not in source:reason="quote_not_exact_contiguous_source_text"
        qn=_norm(quote);etype=item.get("evidence_type")
        candidate_named=any(a and a in qn for a in aliases)
        gene_named=bool(gene and gene in qn)
        target_named=bool(target and target in qn)

        # Positive evidence must be compound-specific. Counterevidence may be
        # class/mechanism-level (e.g. resistance to G12C inhibitors) and thus
        # need not name the individual candidate, but it must name the primary
        # gene or route target.
        if not reason and etype!="MECHANISTIC_COUNTEREVIDENCE" and not candidate_named:
            reason="candidate_name_or_alias_absent_from_exact_quote"
        if not reason and etype=="DIRECT_TARGET_MODULATION" and not (gene_named or target_named):
            reason="gene_or_route_target_absent_from_exact_quote"
        if not reason and etype=="SAME_DRUG_FUNCTIONAL_RESCUE":
            # A generic downstream/cancer rescue around the route target is
            # not patient-mechanism rescue. Require the primary disease-gene
            # context in the exact evidence sentence.
            if not gene_named:
                reason="primary_gene_absent_from_functional_rescue_quote"
        if not reason and etype=="MECHANISTIC_COUNTEREVIDENCE" and not (gene_named or target_named):
            reason="gene_or_route_target_absent_from_counterevidence_quote"
        if reason:
            rejected.append({**item,"rejection_reason":reason});continue
        accepted.append(item)
    return accepted,rejected

def attack_one(candidate_eval,path,outdir,max_docs=8,budget_root=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    candidate=candidate_eval["candidate"]
    aliases=list(dict.fromkeys([candidate,*path.get("member_labels",[])]))
    docs=v28_ai.literature(path.get("gene"),out/"literature",max_docs=max_docs,
                           candidate=aliases,targets=[path.get("route_target")])
    if not docs:
        result={"schema":SCHEMA,"candidate":candidate,"gene":path.get("gene"),
          "route_id":path.get("route_id"),"documents":0,"accepted_evidence":[],
          "rejected_model_items":[],"direct_sources":[],"rescue_sources":[],
          "counterevidence_sources":[],"independent_support_sources":[],
          "evidence_level":"P0_NO_SOURCE_BOUND_SUPPORT","overall_assessment":"INSUFFICIENT",
          "patient_data_sent_external":False}
        (out/"SOURCE_BOUND_ATTACK_RESULT.json").write_text(json.dumps(result,indent=2)+"\n")
        return result
    passages=_passage_payload(docs,candidate,path)
    model=_call_extractor(candidate,path,passages,out/"openai_extraction.json",
                          budget_root or out.parent)
    accepted,rejected=_validate_items(candidate,path,passages,model)
    # Deterministic counter-scan is additive: it can only surface adverse
    # evidence the extractor missed; it never creates positive support.
    for item in _deterministic_counter_items(passages,path,candidate=candidate):
        key=(item.get("source_id"),item.get("passage_id"),item.get("exact_quote"))
        if not any((x.get("source_id"),x.get("passage_id"),x.get("exact_quote"))==key for x in accepted):
            accepted.append(item)
    direct={x["source_id"] for x in accepted if x["evidence_type"]=="DIRECT_TARGET_MODULATION"
            and x.get("mechanistically_relevant")}
    rescue={x["source_id"] for x in accepted if x["evidence_type"]=="SAME_DRUG_FUNCTIONAL_RESCUE"
            and x.get("mechanistically_relevant")}
    counter={x["source_id"] for x in accepted if x["evidence_type"]=="MECHANISTIC_COUNTEREVIDENCE"
             and x.get("mechanistically_relevant")}
    supportive=direct|rescue
    if counter and supportive:level="P0_CONFLICTED_SOURCE_BOUND_EVIDENCE"
    elif direct and rescue and len(supportive)>=2:level="P2_INDEPENDENT_DIRECT_AND_RESCUE"
    elif direct and rescue:level="P1_SINGLE_SOURCE_DIRECT_AND_RESCUE"
    elif supportive:level="P1_PARTIAL_SOURCE_BOUND_SUPPORT"
    elif counter:level="P0_COUNTEREVIDENCE_ONLY"
    else:level="P0_NO_SOURCE_BOUND_SUPPORT"
    result={"schema":SCHEMA,"candidate":candidate,"gene":path.get("gene"),
      "route_id":path.get("route_id"),"objective_id":path.get("objective_id"),
      "documents":len(docs),"passages_supplied_to_model":len(passages),
      "accepted_evidence":accepted,"rejected_model_items":rejected,
      "direct_sources":sorted(direct),"rescue_sources":sorted(rescue),
      "counterevidence_sources":sorted(counter),"independent_support_sources":sorted(supportive),
      "evidence_level":level,"overall_assessment":model.get("value",{}).get("overall_assessment"),
      "missing_evidence":model.get("value",{}).get("missing_evidence",[]),
      "quotes_verified_exactly":True,"patient_data_sent_external":False,
      "clinical_recommendation":False}
    (out/"SOURCE_BOUND_ATTACK_RESULT.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return result

def run(falsification,graph,outdir,count=8,max_docs=8):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    budget_root=out.parent
    selected=select_candidates(falsification,graph,count)
    results=[]
    for i,(candidate,path) in enumerate(selected,1):
        safe=re.sub(r"[^A-Za-z0-9._-]+","_",candidate["candidate"])[:80]
        try:
            results.append(attack_one(candidate,path,out/f"{i:02d}_{safe}",
                                      max_docs=max_docs,budget_root=budget_root))
        except Exception as exc:
            results.append({"schema":SCHEMA,"candidate":candidate.get("candidate"),
              "gene":path.get("gene"),"route_id":path.get("route_id"),
              "error":type(exc).__name__+":"+str(exc)[:600],
              "patient_data_sent_external":False})
    from . import api_budget
    summary={"schema":SCHEMA,"requested":count,"completed":len(results),
      "results":results,"patient_data_sent_external":False,
      "api_budget":api_budget.status(budget_root)}
    (out/"SOURCE_BOUND_ATTACK_SUMMARY.json").write_text(
      json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return summary

def apply(falsification,attack_summary):
    by={x.get("candidate"):x for x in attack_summary.get("results",[]) if x.get("candidate")}
    updated=[]
    for c0 in falsification.get("candidates",[]):
        c=dict(c0);ev=by.get(c.get("candidate"))
        if not ev:
            updated.append(c);continue
        c["source_bound_literature_evidence"]=ev
        if ev.get("counterevidence_sources"):
            c["hold_reasons"]=sorted(set(c.get("hold_reasons",[]))|{"source_bound_counterevidence_requires_resolution"})
        if ev.get("evidence_level")=="P2_INDEPENDENT_DIRECT_AND_RESCUE" and not c.get("kill_reasons"):
            c["supporting_factors"]=sorted(set(c.get("supporting_factors",[]))|{
              "independent_source_bound_direct_target_modulation",
              "independent_source_bound_same_drug_functional_rescue"})
            if c.get("decision")=="HOLD":
                c["decision"]="CONDITIONAL_RESEARCH_LEAD"
                c["hold_reasons"]=sorted(set(c.get("hold_reasons",[]))-{
                  "published_or_assay_evidence_not_yet_A_or_B",
                  "direct_directional_target_assay_missing",
                  "same_drug_disease_relevant_functional_rescue_missing"})
        updated.append(c)
    updated.sort(key=lambda x:(
      {"CONDITIONAL_RESEARCH_LEAD":0,"HOLD":1,"KILLED":2}.get(x.get("decision"),3),
      -x.get("research_priority_score",0),str(x.get("candidate") or "")))
    out=dict(falsification);out["candidates"]=updated
    out["decision_counts"]={d:sum(x.get("decision")==d for x in updated)
      for d in ("CONDITIONAL_RESEARCH_LEAD","HOLD","KILLED")}
    out["source_bound_literature_attack_applied"]=True
    return out


def revalidate_cached_attack(folder,path):
    """Re-apply deterministic evidence gates to an existing attack.

    This performs NO network or OpenAI request. It reuses the exact public
    passages and structured model extraction already stored in the run.
    """
    out=Path(folder)
    candidate=(json.loads((out/"SOURCE_BOUND_ATTACK_RESULT.json").read_text(encoding="utf-8"))
               .get("candidate"))
    bundle=json.loads((out/"literature/literature_bundle.json").read_text(encoding="utf-8"))
    docs=bundle.get("documents") or []
    model_file=out/"openai_extraction.json"
    if not model_file.is_file():
        # No admitted public passages meant the model was never called.
        # Preserve the deterministic no-support result and mark it revalidated.
        previous=json.loads((out/"SOURCE_BOUND_ATTACK_RESULT.json").read_text(encoding="utf-8"))
        previous["gene"]=path.get("gene")
        previous["route_id"]=path.get("route_id")
        previous["objective_id"]=path.get("objective_id")
        previous["revalidated_without_external_request"]=True
        previous["evidence_level"]="P0_NO_SOURCE_BOUND_SUPPORT"
        previous["accepted_evidence"]=[]
        previous["direct_sources"]=[];previous["rescue_sources"]=[]
        previous["counterevidence_sources"]=[];previous["independent_support_sources"]=[]
        (out/"SOURCE_BOUND_ATTACK_RESULT.json").write_text(
          json.dumps(previous,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        return previous
    model=json.loads(model_file.read_text(encoding="utf-8"))
    passages=_passage_payload(docs,candidate,path)
    accepted,rejected=_validate_items(candidate,path,passages,model)
    # Deterministic counter-scan is additive: it can only surface adverse
    # evidence the extractor missed; it never creates positive support.
    for item in _deterministic_counter_items(passages,path,candidate=candidate):
        key=(item.get("source_id"),item.get("passage_id"),item.get("exact_quote"))
        if not any((x.get("source_id"),x.get("passage_id"),x.get("exact_quote"))==key for x in accepted):
            accepted.append(item)
    direct={x["source_id"] for x in accepted if x["evidence_type"]=="DIRECT_TARGET_MODULATION"
            and x.get("mechanistically_relevant")}
    rescue={x["source_id"] for x in accepted if x["evidence_type"]=="SAME_DRUG_FUNCTIONAL_RESCUE"
            and x.get("mechanistically_relevant")}
    counter={x["source_id"] for x in accepted if x["evidence_type"]=="MECHANISTIC_COUNTEREVIDENCE"
             and x.get("mechanistically_relevant")}
    supportive=direct|rescue
    if counter and supportive:level="P0_CONFLICTED_SOURCE_BOUND_EVIDENCE"
    elif direct and rescue and len(supportive)>=2:level="P2_INDEPENDENT_DIRECT_AND_RESCUE"
    elif direct and rescue:level="P1_SINGLE_SOURCE_DIRECT_AND_RESCUE"
    elif supportive:level="P1_PARTIAL_SOURCE_BOUND_SUPPORT"
    elif counter:level="P0_COUNTEREVIDENCE_ONLY"
    else:level="P0_NO_SOURCE_BOUND_SUPPORT"
    result={"schema":SCHEMA,"candidate":candidate,"gene":path.get("gene"),
      "route_id":path.get("route_id"),"objective_id":path.get("objective_id"),
      "documents":len(docs),"passages_supplied_to_model":len(passages),
      "accepted_evidence":accepted,"rejected_model_items":rejected,
      "direct_sources":sorted(direct),"rescue_sources":sorted(rescue),
      "counterevidence_sources":sorted(counter),"independent_support_sources":sorted(supportive),
      "evidence_level":level,"overall_assessment":model.get("value",{}).get("overall_assessment"),
      "missing_evidence":model.get("value",{}).get("missing_evidence",[]),
      "quotes_verified_exactly":True,"revalidated_without_external_request":True,
      "patient_data_sent_external":False,"clinical_recommendation":False}
    (out/"SOURCE_BOUND_ATTACK_RESULT.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return result
