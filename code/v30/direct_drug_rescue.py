"""Direct source-bound approved-drug rescue discovery for V3.1.

This channel catches approved medicines explicitly observed in the drug-blind
mechanistic literature to restore/stabilize a primary gene product. It does
not infer efficacy from class effects or from another compound.
"""
from __future__ import annotations
from pathlib import Path
import json,re,os
from v2810.core import norm_name
from . import openai_reasoner

SCHEMA="RDMVA_DIRECT_DRUG_RESCUE_V31_1"
RESCUE_RX=re.compile(
    r"\b(rescu\w*|recover\w*|restor\w*|stabili[sz]\w*|protect\w*|"
    r"prevent\w*.{0,35}degrad\w*|increase\w*.{0,35}(?:protein|expression|level)|"
    r"half[- ]life|degrad\w*.{0,35}(?:blocked|inhibit\w*))\b",re.I)

def _sentences(text):
    # Keep exact source text; return one- or two-sentence spans.
    parts=re.split(r"(?<=[.!?])\s+",str(text or ""))
    spans=[]
    for i,p in enumerate(parts):
        if p.strip():spans.append(p.strip())
        if i+1<len(parts) and p.strip() and parts[i+1].strip():
            spans.append((p+" "+parts[i+1]).strip())
    return spans

def _term_in(text,term):
    return bool(re.search(r"(?<![A-Za-z0-9])"+re.escape(str(term))+
                          r"(?![A-Za-z0-9])",str(text),re.I))

def _gene_in(text,gene):
    # Exact symbol, or symbol immediately followed by a protein-variant token
    # such as Lztr1R409C/+; reject other gene names such as POU6F2-AS1.
    g=re.escape(str(gene))
    return bool(re.search(r"(?<![A-Za-z0-9])"+g+
                          r"(?=$|[^A-Za-z0-9_-]|[A-Z][0-9])",str(text),re.I))

def _name_map(entities):
    aliases={}
    for e in entities or []:
        names=[e.get("canonical_label"),*(e.get("member_labels") or [])]
        for n in names:
            s=str(n or "").strip()
            if len(s)<6:continue
            aliases.setdefault(s.lower(),e)
    return aliases

def _approved_single(e):
    return bool(e.get("global_approved_as_single") or e.get("solo_approved") or
                e.get("biologic_licensed") or not e.get("combo_only",False))

def _objective_for(gene,objectives,span):
    opts=[o for o in objectives.get("objectives",[])
          if o.get("gene")==gene and o.get("drug_search_eligible")]
    lo=span.lower()
    preferred=[]
    if any(k in lo for k in ("degrad","ubiquit","proteasom","stabili","half-life","half life")):
        preferred=["stabilize_residual_protein","reduce_causal_degradation",
                   "restore_functional_gene_product"]
    elif any(k in lo for k in ("protein","expression","level","recover","restore")):
        preferred=["restore_functional_gene_product","stabilize_residual_protein"]
    for typ in preferred:
        x=next((o for o in opts if o.get("objective_type")==typ),None)
        if x:return x
    return next(iter(opts),None)

def _bundles(run):
    run=Path(run);seen=set()
    for base in ("mechanism_discovery","modifier_discovery"):
        d=run/base
        if not d.is_dir():continue
        for p in d.rglob("literature_bundle.json"):
            rp=str(p.resolve()).lower()
            if rp in seen:continue
            seen.add(rp)
            try:j=json.loads(p.read_text(encoding="utf-8"))
            except Exception:continue
            yield p,j

def discover(run_dir,entities,primary_genes,objectives,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    amap=_name_map(entities)
    genes=[str(g).upper() for g in primary_genes]
    findings=[];seen=set()
    for path,bundle in _bundles(run_dir):
        for doc in bundle.get("documents",[]):
            lic=doc.get("license")
            if not lic:continue
            sid=doc.get("source_id")
            for passage in doc.get("passages",[]):
                text=str(passage.get("text") or "")
                lo=text.lower()
                present_genes=[g for g in genes if _gene_in(text,g)]
                if not present_genes or not RESCUE_RX.search(text):continue
                # Candidate aliases mentioned anywhere in the passage.
                mentioned=[]
                for alias,e in amap.items():
                    # Fast substring prefilter; regex boundaries are evaluated
                    # only for the tiny subset actually present in this passage.
                    if alias in lo and _term_in(text,alias) and _approved_single(e):
                        mentioned.append((alias,e))
                if not mentioned:continue
                for span in _sentences(text):
                    slo=span.lower()
                    if not RESCUE_RX.search(span):continue
                    for gene in present_genes:
                        if not _gene_in(span,gene):
                            continue
                        for alias,e in mentioned:
                            if alias not in slo or not _term_in(span,alias):continue
                            # Require gene and approved medicine in the same exact span.
                            obj=_objective_for(gene,objectives,span)
                            if not obj:continue
                            drug=e.get("canonical_label")
                            key=(drug,gene,sid,passage.get("id"),obj.get("objective_id"))
                            if key in seen:continue
                            seen.add(key)
                            findings.append({
                              "candidate":drug,"entity_id":e.get("entity_id"),
                              "scientific_group_id":e.get("scientific_group_id"),
                              "member_labels":e.get("member_labels") or [drug],
                              "gene":gene,"objective_id":obj.get("objective_id"),
                              "objective_type":obj.get("objective_type"),
                              "objective_priority":obj.get("priority_tier"),
                              "source_id":sid,"passage_id":passage.get("id"),
                              "citation":doc.get("citation"),"license":lic,
                              "exact_quote":span,
                              "observed_effect":"primary_gene_product_rescue_or_stabilization",
                              "disease_functional_rescue_demonstrated":False,
                              "approved_single_agent":True,
                              "patient_context_sent_external":False,
                              "clinical_recommendation":False})
    payload={
      "schema":SCHEMA,"findings":findings,"finding_count":len(findings),
      "approved_drug_and_primary_gene_required_same_exact_span":True,
      "class_effect_extrapolation_allowed":False,
      "disease_functional_rescue_inferred":False,
      "clinical_recommendation":False}
    (out/"DIRECT_DRUG_RESCUE_V31.json").write_text(
        json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload

def merge_pairs(pairs,payload,objectives):
    out=list(pairs or []);existing={(p.get("ingredient"),p.get("route_id")) for p in out}
    added=0;routes=[]
    for i,x in enumerate(payload.get("findings",[]),1):
        rid=f"V31DDR{i:04d}"
        route={
          "route_id":rid,"rank":0,"family":"v31_direct_source_drug_rescue",
          "case_gene":x.get("gene"),
          "target_concept":f"{x.get('gene')} protein abundance",
          "targets":[],"desired_action":"increase","variable":"functional_abundance",
          "mechanism_concept":"approved compound directly observed to restore/stabilize primary-gene product in source experiment",
          "score":18.0,"state":"conditional","network_hops":0,
          "rescue_objective_id":x.get("objective_id"),
          "rescue_objective_ids":[x.get("objective_id")],
          "rescue_objective_type":x.get("objective_type"),
          "rescue_objective_types":[x.get("objective_type")],
          "rescue_priority_tier":x.get("objective_priority"),
          "source":x.get("source_id"),
          "evidence":[{"source_id":x.get("source_id"),"passage_id":x.get("passage_id"),
                       "quote":x.get("exact_quote"),"license":x.get("license")}],
          "source_bound_direct_drug_rescue":True,
          "scientific_gate":{
            "schema":"V31_DIRECT_DRUG_RESCUE_GATE_1",
            "primary_gene_target_identity":True,
            "allele_specific_rescue_verified":False,
            "eligible_for_exploratory_screen":True,
            "functional_rescue_validated":False,
            "causal_direction_validated_for_disease":False,
            "requires_independent_endpoint":True,
            "hard_hold":False,
            "issues":["primary_gene_product_rescue_observed_but_disease_functional_rescue_not_demonstrated"]}}
        routes.append(route)
        key=(x.get("candidate"),rid)
        if key in existing:continue
        out.append({
          "schema":"RDMVA_DIRECT_DRUG_RESCUE_PAIR_V31_1",
          "ingredient":x.get("candidate"),"entity_id":x.get("entity_id"),
          "scientific_group_id":x.get("scientific_group_id"),
          "member_labels":x.get("member_labels"),
          "route_id":rid,"route":route,"tier":"C","score":12.0,
          "providers":["source_bound_direct_drug_rescue"],
          "signals":[{
            "route_id":rid,"provider":"source_bound_direct_drug_rescue",
            "source_id":x.get("source_id"),"passage_id":x.get("passage_id"),
            "exact_quote":x.get("exact_quote"),"license":x.get("license"),
            "direct_target_assay_verified":False,
            "same_drug_primary_gene_product_rescue":True,
            "functional_rescue":False}],
          "combo_only":False,"global_approved_as_single":True,
          "publication_strength":"source_bound_primary_gene_product_rescue",
          "evidence_blockers":["disease_functional_rescue_not_demonstrated",
                               "patient_variant_specific_rescue_not_demonstrated"],
          "shortest_network_hops":0,"drug_to_route_hops":0,"route_to_case_hops":0,
          "source_bound_primary_gene_rescue":True,
          "scientific_gate":route["scientific_gate"],
          "promotion_eligible_validated_lead":False,
          "clinical_recommendation":False})
        existing.add(key);added+=1
    return out,{"added_pairs":added,"routes":routes,
                "validated_disease_rescue_claimed":False}


def _objective_for_class(gene,objectives,evidence_class,quote):
    opts=[o for o in objectives.get("objectives",[])
          if o.get("gene")==gene and o.get("drug_search_eligible")]
    if evidence_class=="MUTANT_MODEL_FUNCTIONAL_RESCUE":
        pref=["bypass_or_compensate_downstream","compensate_downstream_deficit",
              "restore_functional_gene_product"]
    elif evidence_class=="PRIMARY_GENE_PRODUCT_RESCUE":
        lo=str(quote or "").lower()
        pref=(["stabilize_residual_protein","reduce_causal_degradation",
               "restore_functional_gene_product"]
              if any(k in lo for k in ("degrad","ubiquit","proteasom","stabili","half-life","half life"))
              else ["restore_functional_gene_product","stabilize_residual_protein"])
    else:
        pref=[]
    for typ in pref:
        x=next((o for o in opts if o.get("objective_type")==typ),None)
        if x:return x
    return next(iter(opts),None)

def adjudicate(payload,objectives,outfile,budget_root):
    """AI-assisted classification with deterministic source-bound validation.

    The model only classifies already-extracted exact public passages. It cannot
    add a drug, gene, quote or evidence that was not found deterministically.
    """
    findings=list(payload.get("findings") or [])
    if not findings:
        out={**payload,"accepted_findings":[],"rejected_findings":[],
             "adjudication_completed":True}
        Path(outfile).write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        return out

    cfg=openai_reasoner.config()
    model=os.environ.get(cfg.get("environment_literature_model_override","RDMVA_OPENAI_LITERATURE_MODEL"),
                         cfg.get("literature_extraction_model") or cfg["default_model"])
    schema={
      "type":"object","additionalProperties":False,
      "properties":{"reviews":{"type":"array","items":{
        "type":"object","additionalProperties":False,
        "properties":{
          "finding_index":{"type":"integer","minimum":0},
          "evidence_class":{"type":"string","enum":[
            "MUTANT_MODEL_FUNCTIONAL_RESCUE",
            "PRIMARY_GENE_PRODUCT_RESCUE",
            "INDIRECT_DOWNSTREAM_EFFECT",
            "NON_CAUSAL_COOCCURRENCE",
            "ASSAY_OR_INDUCTION_REAGENT_ONLY",
            "OPPOSITE_OR_WORSENING",
            "INSUFFICIENT"]},
          "causal_intervention":{"type":"boolean"},
          "primary_gene_product_rescue":{"type":"boolean"},
          "mutant_model_functional_rescue":{"type":"boolean"},
          "nonhuman_model":{"type":"boolean"},
          "reason":{"type":"string","maxLength":700}
        },
        "required":["finding_index","evidence_class","causal_intervention",
                    "primary_gene_product_rescue","mutant_model_functional_rescue",
                    "nonhuman_model","reason"]}}},
      "required":["reviews"]}
    inp=[]
    for i,x in enumerate(findings):
        inp.append({
          "finding_index":i,"candidate":x.get("candidate"),"gene":x.get("gene"),
          "citation":x.get("citation"),"source_id":x.get("source_id"),
          "exact_quote":x.get("exact_quote")})
    instructions=(
      "Classify whether the named approved medicine CAUSALLY produced a rescue in the exact quoted experiment. "
      "MUTANT_MODEL_FUNCTIONAL_RESCUE requires the medicine to improve a phenotype/function in a model carrying "
      "a mutation/loss/deficiency of the named gene. PRIMARY_GENE_PRODUCT_RESCUE requires the medicine itself to "
      "restore, recover, stabilize or increase the named gene's protein/product after a deficit/degradation context. "
      "If the medicine is merely a comparison treatment, resistance agent, assay reagent, induction reagent, "
      "co-occurring exposure, or the effect concerns a different gene/transcript, reject it. "
      "Do not infer class effects or mechanisms beyond the exact quote."
    )
    body={"model":model,"input":json.dumps(inp,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],
          "max_output_tokens":2400,
          "text":{"format":{"type":"json_schema","name":"v31_direct_drug_rescue_adjudication",
                           "strict":True,"schema":schema}}}
    rec=openai_reasoner.request_json(body,Path(outfile).with_suffix(".ai.json"),
                                     budget_root=budget_root)
    reviews={int(x["finding_index"]):x for x in rec.get("value",{}).get("reviews",[])
             if isinstance(x.get("finding_index"),int)}
    accepted=[];rejected=[]
    for i,x0 in enumerate(findings):
        x=dict(x0);rv=reviews.get(i) or {
          "evidence_class":"INSUFFICIENT","causal_intervention":False,
          "primary_gene_product_rescue":False,"mutant_model_functional_rescue":False,
          "nonhuman_model":False,"reason":"No adjudication returned."}
        cls=rv.get("evidence_class")
        ok=bool(rv.get("causal_intervention") and cls in (
            "MUTANT_MODEL_FUNCTIONAL_RESCUE","PRIMARY_GENE_PRODUCT_RESCUE"))
        # Deterministic invariants: candidate and primary gene remain in exact source span.
        if ok and not _term_in(x.get("exact_quote"),x.get("candidate")):ok=False
        if ok and not _gene_in(x.get("exact_quote"),x.get("gene")):ok=False
        x["adjudication"]=rv
        if ok:
            obj=_objective_for_class(x.get("gene"),objectives,cls,x.get("exact_quote"))
            if not obj:
                x["rejection_reason"]="no_matching_rescue_objective";rejected.append(x);continue
            x["objective_id"]=obj.get("objective_id")
            x["objective_type"]=obj.get("objective_type")
            x["objective_priority"]=obj.get("priority_tier")
            x["evidence_class"]=cls
            x["translational_context"]="NONHUMAN_ONLY" if rv.get("nonhuman_model") else "HUMAN_OR_UNRESOLVED"
            x["source_bound_primary_gene_rescue"]=cls=="PRIMARY_GENE_PRODUCT_RESCUE"
            x["source_bound_mutant_model_functional_rescue"]=cls=="MUTANT_MODEL_FUNCTIONAL_RESCUE"
            accepted.append(x)
        else:
            x["rejection_reason"]=cls
            rejected.append(x)
    # Collapse duplicate passages from the same experiment/source into one
    # drug-gene-evidence hypothesis while preserving every exact supporting span.
    grouped={}
    for x in accepted:
        key=(x.get("candidate"),x.get("gene"),x.get("evidence_class"),
             x.get("objective_id"),x.get("source_id"))
        if key not in grouped:
            y=dict(x)
            y["supporting_passages"]=[{
              "passage_id":x.get("passage_id"),"exact_quote":x.get("exact_quote"),
              "citation":x.get("citation"),"license":x.get("license")}]
            grouped[key]=y
        else:
            grouped[key]["supporting_passages"].append({
              "passage_id":x.get("passage_id"),"exact_quote":x.get("exact_quote"),
              "citation":x.get("citation"),"license":x.get("license")})
    accepted=list(grouped.values())
    out={**payload,"accepted_findings":accepted,"rejected_findings":rejected,
         "accepted_count":len(accepted),"rejected_count":len(rejected),
         "adjudication_completed":True,
         "model_cannot_add_source_evidence":True,
         "clinical_recommendation":False}
    Path(outfile).write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return out


def merge_pairs(pairs,payload,objectives):
    """Merge only adjudicated direct-rescue findings into the candidate graph."""
    findings=(payload.get("accepted_findings") if payload.get("adjudication_completed")
              else [])
    out=list(pairs or []);existing={(p.get("ingredient"),p.get("route_id")) for p in out}
    added=0;routes=[]
    for i,x in enumerate(findings,1):
        rid=f"V31DDR{i:04d}"
        cls=x.get("evidence_class")
        mutant=bool(x.get("source_bound_mutant_model_functional_rescue"))
        primary=bool(x.get("source_bound_primary_gene_rescue"))
        nonhuman=x.get("translational_context")=="NONHUMAN_ONLY"
        issue=("mutant_model_rescue_requires_human_and_patient_variant_validation"
               if mutant else "primary_gene_product_rescue_requires_disease_function_validation")
        route={
          "route_id":rid,"rank":0,"family":"v31_direct_source_drug_rescue",
          "case_gene":x.get("gene"),
          "target_concept":(
              f"{x.get('gene')} mutant phenotype" if mutant
              else f"{x.get('gene')} protein abundance"),
          "targets":[],"desired_action":"restore",
          "variable":"disease_relevant_function" if mutant else "functional_abundance",
          "mechanism_concept":(
              "approved compound directly observed to improve a phenotype/function in a primary-gene mutant model"
              if mutant else
              "approved compound directly observed to restore/stabilize primary-gene product"),
          "score":19.0 if mutant else 18.0,"state":"conditional","network_hops":0,
          "rescue_objective_id":x.get("objective_id"),
          "rescue_objective_ids":[x.get("objective_id")],
          "rescue_objective_type":x.get("objective_type"),
          "rescue_objective_types":[x.get("objective_type")],
          "rescue_priority_tier":x.get("objective_priority"),
          "source":x.get("source_id"),
          "evidence":[{
             "source_id":x.get("source_id"),"passage_id":sp.get("passage_id"),
             "quote":sp.get("exact_quote"),"license":sp.get("license") or x.get("license")}
             for sp in (x.get("supporting_passages") or [{
               "passage_id":x.get("passage_id"),"exact_quote":x.get("exact_quote"),
               "license":x.get("license")}])],
          "source_bound_direct_drug_rescue":True,
          "source_bound_primary_gene_rescue":primary,
          "source_bound_mutant_model_functional_rescue":mutant,
          "translational_context":x.get("translational_context"),
          "scientific_gate":{
            "schema":"V31_DIRECT_DRUG_RESCUE_GATE_2",
            "primary_gene_target_identity":True,
            "allele_specific_rescue_verified":False,
            "eligible_for_exploratory_screen":True,
            "functional_rescue_validated":False,
            "causal_direction_validated_for_disease":bool(mutant and not nonhuman),
            "requires_independent_endpoint":True,
            "hard_hold":False,
            "issues":[issue]+(["nonhuman_model_requires_human_validation"] if nonhuman else [])}}
        routes.append(route)
        key=(x.get("candidate"),rid)
        if key in existing:continue
        blockers=["patient_variant_specific_rescue_not_demonstrated"]
        if mutant:
            blockers+=["independent_replication_required",
                       "human_validation_required" if nonhuman else "independent_human_validation_required"]
        else:
            blockers+=["disease_functional_rescue_not_demonstrated",
                       "independent_replication_required"]
        out.append({
          "schema":"RDMVA_DIRECT_DRUG_RESCUE_PAIR_V31_2",
          "ingredient":x.get("candidate"),"entity_id":x.get("entity_id"),
          "scientific_group_id":x.get("scientific_group_id"),
          "member_labels":x.get("member_labels"),
          "route_id":rid,"route":route,"tier":"C",
          "score":15.0 if mutant else 12.0,
          "providers":["source_bound_direct_drug_rescue"],
          "signals":[{
            "route_id":rid,"provider":"source_bound_direct_drug_rescue",
            "source_id":x.get("source_id"),"passage_id":x.get("passage_id"),
            "exact_quote":x.get("exact_quote"),"license":x.get("license"),
            "evidence_class":cls,
            "direct_target_assay_verified":False,
            "same_drug_primary_gene_product_rescue":primary,
            "same_drug_mutant_model_functional_rescue":mutant,
            "functional_rescue":False}],
          "combo_only":False,"global_approved_as_single":True,
          "publication_strength":(
              "source_bound_mutant_model_functional_rescue"
              if mutant else "source_bound_primary_gene_product_rescue"),
          "evidence_blockers":blockers,
          "shortest_network_hops":0,"drug_to_route_hops":0,"route_to_case_hops":0,
          "source_bound_primary_gene_rescue":primary,
          "source_bound_mutant_model_functional_rescue":mutant,
          "direct_source_rescue_evidence_class":cls,
          "translational_context":x.get("translational_context"),
          "scientific_gate":route["scientific_gate"],
          "promotion_eligible_validated_lead":False,
          "clinical_recommendation":False})
        existing.add(key);added+=1
    out.sort(key=lambda x:({"A":0,"B":1,"C":2,"D":3}.get(x.get("tier"),4),
                           -float(x.get("score") or 0),str(x.get("ingredient")),str(x.get("route_id"))))
    return out,{"added_pairs":added,"routes":routes,
                "accepted_findings":len(findings),
                "validated_disease_rescue_claimed":False}
