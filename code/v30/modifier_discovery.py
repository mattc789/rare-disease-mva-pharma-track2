"""Drug-blind source-bound modifier discovery for V3.1.

Looks for genetic/biological suppressors, bypasses, proteostasis regulators and
other modifiers before drug names enter the discovery logic. Every accepted
route is bound to an exact public quote. AI cannot invent evidence.
"""
from __future__ import annotations
from pathlib import Path
import json,os,re
from v2810 import ai as v28_ai,routes as v28_routes,science as v28_science
from . import openai_reasoner

SCHEMA="RDMVA_MODIFIER_DISCOVERY_V31_1"
GENE_RX=re.compile(r"^[A-Z][A-Z0-9-]{1,19}$")

def _gene_exact(text,gene):
    """Exact gene-symbol token match with narrow experimental suffix support.

    BUB1B-OE/KO/KD are accepted as BUB1B experimental states, while biological
    names such as POU6F2-AS1 are distinct symbols and must not match POU6F2.
    """
    if not text or not gene:return False
    g=re.escape(str(gene).upper())
    pat=(r"(?<![A-Z0-9_-])"+g+
         r"(?:(?=[^A-Z0-9_-]|$)|(?=-(?:OE|KO|KD|WT|MUT|NULL)\b))")
    return bool(re.search(pat,str(text).upper()))

def _passages(docs,max_passages=36):
    out=[]
    for d in docs:
        for p in d.get("passages") or []:
            out.append({"source_id":d.get("source_id"),"citation":d.get("citation"),
                        "license":d.get("license"),"passage_id":p.get("id"),
                        "text":str(p.get("text") or "")[:2600]})
    return out[:max_passages]

def _extract(primary_gene,passages,outfile,budget_root,rescue_context=None):
    cfg=openai_reasoner.config()
    model=os.environ.get(cfg.get("environment_literature_model_override","RDMVA_OPENAI_LITERATURE_MODEL"),
                         cfg.get("literature_extraction_model") or cfg["default_model"])
    schema={"type":"object","additionalProperties":False,"properties":{
      "modifiers":{"type":"array","items":{"type":"object","additionalProperties":False,
        "properties":{
          "modifier_gene":{"type":"string"},
          "intervention_action":{"type":"string","enum":["increase","decrease","restore","unknown"]},
          "effect_on_primary_defect":{"type":"string","enum":["improves","worsens","unknown"]},
          "mechanism_class":{"type":"string","enum":[
            "GENETIC_SUPPRESSOR","GENETIC_ENHANCER","BYPASS","PROTEOSTASIS",
            "DEGRADATION","SPICING","TRAFFICKING","COMPLEX_ASSEMBLY",
            "METABOLIC_COMPENSATION","SIGNALING_COMPENSATION","OTHER"]},
          "functional_rescue_observed":{"type":"boolean"},
          "experiment_context":{"type":"string"},
          "source_id":{"type":"string"},"passage_id":{"type":"string"},
          "exact_quote":{"type":"string","maxLength":1200},
          "reason":{"type":"string","maxLength":700}
        },
        "required":["modifier_gene","intervention_action","effect_on_primary_defect",
                    "mechanism_class","functional_rescue_observed","experiment_context",
                    "source_id","passage_id","exact_quote","reason"]}},
      "unresolved_questions":{"type":"array","items":{"type":"string"}}
    },"required":["modifiers","unresolved_questions"]}
    payload={"primary_gene":primary_gene,
             "primary_gene_rescue_context":rescue_context or {},
             "public_passages":passages,
             "patient_variant_or_phenotype_data_included":False,
             "drug_names_included":False}
    instructions=(
      "Find experimentally supported biological modifiers of the primary gene defect using ONLY the supplied passages. "
      "A modifier is useful only if perturbing a named human gene/protein (increase or decrease) suppresses, rescues, "
      "compensates for, stabilizes, bypasses, or worsens a functional defect connected to the primary gene. "
      "Do not infer direction from mere interaction, coexpression, pathway membership, cancer correlation or binding. "
      "Respect primary_gene_rescue_context: if the patient objective is to restore/increase the primary gene, experiments "
      "performed only in primary-gene overexpression states are opposite-context and must not be reported as rescue; likewise "
      "loss/knockout-only experiments are opposite-context when the patient objective is to decrease excessive primary-gene activity. "
      "GENETIC_SUPPRESSOR requires an experiment in which modifier perturbation improves a relevant phenotype/function in a compatible state. "
      "Copy the exact quote contiguously. If direction/action is not explicit use unknown. Do not mention drugs or treatment."
    )
    body={"model":model,"input":json.dumps(payload,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],
          "max_output_tokens":2600,
          "text":{"format":{"type":"json_schema","name":"v31_modifier_discovery",
                           "strict":True,"schema":schema}}}
    return openai_reasoner.request_json(body,outfile,budget_root=budget_root)

def _source_exact_modifier_quote(text,primary_gene,modifier_gene):
    """Recover an exact contiguous source sentence/span when an AI quote is
    substantively correct but omitted inline figure/citation text.

    The recovered span must itself contain both gene symbols and a mechanistic
    modifier keyword; no paraphrase is ever manufactured.
    """
    if not text:return None
    parts=re.split(r"(?<=[.!?])\s+",str(text))
    keys=("rescu","suppress","stabili","half-life","half life","degrad",
          "ubiquit","compens","restore","knockdown","overexpress","silenc",
          "increase","decrease","prevent")
    for width in (1,2):
        for i in range(0,max(0,len(parts)-width+1)):
            span=" ".join(parts[i:i+width]).strip()
            lo=span.lower()
            if (_gene_exact(span,primary_gene) and _gene_exact(span,modifier_gene)
                    and any(k in lo for k in keys)):
                return span
    return None

def _validate(primary_gene,passages,rec,rescue_context=None):
    pmap={(x["source_id"],x["passage_id"]):x["text"] for x in passages}
    accepted=[];rejected=[]
    for x0 in rec.get("value",{}).get("modifiers",[]):
        x=dict(x0);reason=None
        quote=x.get("exact_quote") or ""
        text=pmap.get((x.get("source_id"),x.get("passage_id")))
        gene=str(x.get("modifier_gene") or "").upper().strip()
        repaired=False
        if text and (not quote or quote not in text):
            fixed=_source_exact_modifier_quote(text,primary_gene,gene)
            if fixed:
                quote=fixed;x["exact_quote"]=fixed
                x["quote_repaired_from_source"]=True
                repaired=True
        qlo=quote.lower()
        if not text:reason="unknown_source_or_passage"
        elif not quote or quote not in text:reason="quote_not_exact"
        elif not GENE_RX.fullmatch(gene):reason="modifier_not_gene_symbol"
        elif gene==str(primary_gene).upper():reason="modifier_same_as_primary_gene"
        elif not _gene_exact(quote,gene):reason="modifier_gene_absent_from_exact_quote"
        elif not _gene_exact(quote,primary_gene):
            # A functional rescue can legitimately name an experimental
            # ortholog rather than the human primary symbol. Accept only when
            # the SAME verified passage explicitly links the human primary
            # gene and the modifier/rescue context.
            if (x.get("functional_rescue_observed") and
                    _gene_exact(text,primary_gene) and _gene_exact(text,gene)):
                x["primary_gene_linkage"]="same_verified_passage"
            else:
                reason="primary_gene_absent_from_exact_quote_and_passage_linkage_insufficient"
        else:
            x["primary_gene_linkage"]="same_exact_quote"
        if not reason and x.get("intervention_action") not in ("increase","decrease"):
            reason="modifier_direction_not_explicit"
        elif not reason and x.get("effect_on_primary_defect")!="improves":
            reason="not_rescue_or_compensation"
        if not reason:
            desired=set((rescue_context or {}).get("desired_actions") or [])
            # Opposite-context protection should inspect the source passage,
            # not only the short extracted quote.
            slo=(str(text)+" "+str(x.get("experiment_context") or "")).lower()
            scompact=re.sub(r"[^a-z0-9]+"," ",slo)
            pg=primary_gene.lower()
            overexpr=(f"{pg} oe" in scompact or
                      re.search(rf"\b{re.escape(pg)}\b.{{0,45}}overexpress",slo) or
                      re.search(rf"overexpress.{{0,45}}\b{re.escape(pg)}\b",slo))
            lossctx=(re.search(rf"\b(loss|deficien|knockout|knockdown|silenc)\b.{{0,45}}\b{re.escape(pg)}\b",slo) or
                     re.search(rf"\b{re.escape(pg)}\b.{{0,45}}\b(loss|deficien|knockout|knockdown|silenc)\b",slo))
            # Only reject when the candidate rescue itself is explicitly tied
            # to the opposite primary-gene state; broad papers often discuss
            # both gain and loss contexts in the same paragraph.
            qstate=(quote+" "+str(x.get("experiment_context") or "")).lower()
            # Primary-gene state must be syntactically attached to the primary
            # gene; "KIFC1 overexpression ... BUB1B" is NOT BUB1B overexpression.
            qover=(
              re.search(rf"\b{re.escape(pg)}(?:-|\s+)(?:oe|overexpress(?:ed|ing|ion)?)\b",qstate) or
              re.search(rf"\boverexpress(?:ion|ed|ing)?\s+(?:of\s+)?(?:the\s+)?\b{re.escape(pg)}\b",qstate)
            )
            qloss=(
              re.search(rf"\b{re.escape(pg)}(?:-|\s+)(?:ko|kd|null|deficien(?:t|cy)?|knockout|knockdown|silenc(?:ed|ing)?)\b",qstate) or
              re.search(rf"\b(?:loss|deficien(?:cy|t)?|knockout|knockdown|silenc(?:ing|ed)?)\s+(?:of\s+)?(?:the\s+)?\b{re.escape(pg)}\b",qstate)
            )
            if desired & {"increase","restore"} and qover:
                reason="opposite_primary_gene_overexpression_context"
            elif desired=={"decrease"} and qloss:
                reason="opposite_primary_gene_loss_context"
        if reason:
            rejected.append({**x,"rejection_reason":reason})
        else:
            context=(str(x.get("experiment_context") or "")+" "+str(text)).lower()
            nonhuman_terms=("drosophila","fruit fly","mouse","mice","murine","zebrafish",
                            "c. elegans","caenorhabditis","xenopus","yeast","rat model"," rats ")
            context_class="NONHUMAN_ONLY" if any(term in context for term in nonhuman_terms) else "HUMAN_OR_UNRESOLVED"
            accepted.append({**x,
                             "quote_repaired_from_source":bool(x.get("quote_repaired_from_source",False)),
                             "source_passage_contains_primary_gene":_gene_exact(text,primary_gene),
                             "translational_context":context_class,
                             "nonhuman_rescue_is_human_validation":False})
    return accepted,rejected

def discover(primary_genes,network_expansion,outdir,max_docs=18,max_network_targets=18,objectives=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    rescue_by_gene={}
    for o in (objectives or {}).get("objectives",[]):
        g=o.get("gene")
        if not g:continue
        rec=rescue_by_gene.setdefault(g,{"desired_actions":set(),"objective_types":[]})
        if o.get("desired_action"):rec["desired_actions"].add(o.get("desired_action"))
        if o.get("objective_type"):rec["objective_types"].append(o.get("objective_type"))
    for g,rec in rescue_by_gene.items():
        rec["desired_actions"]=sorted(rec["desired_actions"])
        rec["objective_types"]=sorted(set(rec["objective_types"]))
    candidates=[x.get("gene") for x in network_expansion.get("modifier_candidates",[])
                if x.get("gene") and not x.get("seed_gene")]
    all_acc=[];all_rej=[]
    for gene in primary_genes:
        # Run three diversified target batches so the unsigned network can
        # improve recall without supplying pharmacologic direction.
        batches=[[]]
        top=[x for x in candidates if x!=gene][:max_network_targets]
        for i in range(0,len(top),6):batches.append(top[i:i+6])
        seen=set()
        for bi,batch in enumerate(batches):
            docs=v28_ai.literature(gene,out/gene/f"batch_{bi:02d}",
                                   max_docs=max_docs,targets=batch)
            passages=_passages(docs)
            if not passages:continue
            context=rescue_by_gene.get(gene,{})
            rec=_extract(gene,passages,out/gene/f"batch_{bi:02d}"/"modifier_ai.json",out,context)
            acc,rej=_validate(gene,passages,rec,context)
            for x in acc:
                key=(gene,x["modifier_gene"],x["intervention_action"],
                     x["source_id"],x["passage_id"])
                if key in seen:continue
                seen.add(key);all_acc.append({"primary_gene":gene,**x})
            all_rej.extend({"primary_gene":gene,**x} for x in rej)
    payload={"schema":SCHEMA,"accepted_modifiers":all_acc,
             "rejected_items":all_rej,
             "accepted_count":len(all_acc),
             "drug_names_used_in_discovery":False,
             "patient_raw_data_sent_to_ai":False,
             "direction_requires_exact_source_quote":True,
             "clinical_recommendation":False}
    (out/"MODIFIER_DISCOVERY_V31.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload

def to_routes(modifiers,cache_dir,objectives=None):
    by_gene={}
    for o in (objectives or {}).get("objectives",[]):
        by_gene.setdefault(o.get("gene"),[]).append(o)
    pref={
      "PROTEOSTASIS":["stabilize_residual_protein","reduce_causal_degradation",
                       "restore_functional_gene_product"],
      "DEGRADATION":["reduce_causal_degradation","stabilize_residual_protein"],
      "SPICING":["restore_productive_splicing"],
      "COMPLEX_ASSEMBLY":["restore_complex_or_interaction_function"],
      "BYPASS":["bypass_or_compensate_downstream","compensate_downstream_deficit"],
      "SIGNALING_COMPENSATION":["bypass_or_compensate_downstream","compensate_downstream_deficit"],
      "METABOLIC_COMPENSATION":["bypass_or_compensate_downstream","compensate_downstream_deficit"],
      "GENETIC_SUPPRESSOR":["bypass_or_compensate_downstream","compensate_downstream_deficit",
                            "restore_functional_gene_product"],
      "GENETIC_ENHANCER":["bypass_or_compensate_downstream","compensate_downstream_deficit"],
      "TRAFFICKING":["restore_functional_gene_product","restore_complex_or_interaction_function"],
      "OTHER":["bypass_or_compensate_downstream","compensate_downstream_deficit"],
    }
    routes=[]
    for x in modifiers.get("accepted_modifiers",[]):
        gene=x.get("modifier_gene");primary=x.get("primary_gene")
        try:acc=v28_routes.accession_for_gene(gene,cache_dir)
        except Exception:acc=None
        if not acc:continue
        action=x.get("intervention_action")
        if action not in ("increase","decrease"):continue
        opts=by_gene.get(primary,[])
        wanted=pref.get(x.get("mechanism_class"),pref["OTHER"])
        chosen=next((o for typ in wanted for o in opts
                     if o.get("objective_type")==typ and o.get("drug_search_eligible")),None)
        if chosen is None:
            chosen=next((o for o in opts if o.get("drug_search_eligible")),None)
        nonhuman=x.get("translational_context")=="NONHUMAN_ONLY"
        score=(12.0 if nonhuman else 16.0) if x.get("functional_rescue_observed") else 13.0
        route={
          "family":"v31_source_bound_modifier",
          "case_gene":primary,"target_concept":gene,
          "targets":[{"gene":gene,"accession":acc}],
          "desired_action":action,"variable":"activity",
          "mechanism_concept":f"{action} {gene} as source-bound modifier of {primary}",
          "score":score,"state":"conditional",
          "source":x.get("source_id"),
          "evidence_context":x.get("experiment_context"),
          "evidence":[{
             "source_id":x.get("source_id"),"passage_id":x.get("passage_id"),
             "quote":x.get("exact_quote"),"license":"open_access_source_bound"}],
          "modifier_mechanism_class":x.get("mechanism_class"),
          "functional_rescue_observed_in_source":bool(x.get("functional_rescue_observed")),
          "translational_context":x.get("translational_context","HUMAN_OR_UNRESOLVED"),
          "source_bound_direction_verified":True,
          "network_hops":0,
        }
        if chosen:
            route["rescue_objective_id"]=chosen.get("objective_id")
            route["rescue_objective_ids"]=[chosen.get("objective_id")]
            route["rescue_objective_type"]=chosen.get("objective_type")
            route["rescue_objective_types"]=[chosen.get("objective_type")]
            route["rescue_priority_tier"]=chosen.get("priority_tier")
        routes.append(route)
    for i,r in enumerate(routes,1):
        r["route_id"]=f"V31M{i:04d}";r["rank"]=i
        nonhuman=r.get("translational_context")=="NONHUMAN_ONLY"
        source_rescue=bool(r["functional_rescue_observed_in_source"])
        validated=bool(source_rescue and not nonhuman)
        issues=[]
        if nonhuman:issues.append("nonhuman_model_requires_human_validation")
        if source_rescue and not validated:issues.append("source_rescue_not_human_validated")
        if not source_rescue:issues.append("modifier_direction_source_bound_but_disease_rescue_not_directly_validated")
        r["scientific_gate"]={
          "schema":"V31_SOURCE_BOUND_MODIFIER_GATE_1",
          "primary_gene_target_identity":False,
          "allele_specific_rescue_verified":False,
          "eligible_for_exploratory_screen":True,
          "functional_rescue_validated":validated,
          "causal_direction_validated_for_disease":validated,
          "requires_independent_endpoint":True,
          "hard_hold":False,
          "issues":issues,
        }
    return routes


def revalidate_existing(primary_genes,outdir,objectives=None):
    """Re-run deterministic source validation over persisted literature/AI outputs.

    No network or OpenAI call is made. This supports stricter/safer validator
    upgrades without repeating expensive discovery.
    """
    out=Path(outdir)
    rescue_by_gene={}
    for o in (objectives or {}).get("objectives",[]):
        g=o.get("gene")
        if not g:continue
        rec=rescue_by_gene.setdefault(g,{"desired_actions":set(),"objective_types":[]})
        if o.get("desired_action"):rec["desired_actions"].add(o.get("desired_action"))
        if o.get("objective_type"):rec["objective_types"].append(o.get("objective_type"))
    for rec in rescue_by_gene.values():
        rec["desired_actions"]=sorted(rec["desired_actions"])
        rec["objective_types"]=sorted(set(rec["objective_types"]))
    accepted=[];rejected=[];seen=set()
    for gene in primary_genes:
        gd=out/gene
        if not gd.is_dir():continue
        for bd in sorted(gd.glob("batch_*")):
            lp=bd/"literature_bundle.json";ap=bd/"modifier_ai.json"
            if not lp.is_file() or not ap.is_file():continue
            bundle=json.loads(lp.read_text(encoding="utf-8"))
            ai=json.loads(ap.read_text(encoding="utf-8"))
            passages=_passages(bundle.get("documents") or [])
            context=rescue_by_gene.get(gene,{})
            acc,rej=_validate(gene,passages,ai,context)
            for x in acc:
                key=(gene,x.get("modifier_gene"),x.get("intervention_action"),
                     x.get("source_id"),x.get("passage_id"))
                if key in seen:continue
                seen.add(key);accepted.append({"primary_gene":gene,**x})
            rejected.extend({"primary_gene":gene,**x} for x in rej)
    payload={"schema":SCHEMA,"accepted_modifiers":accepted,
             "rejected_items":rejected,"accepted_count":len(accepted),
             "revalidated_from_persisted_source_bundles":True,
             "network_or_ai_calls_performed":False,
             "drug_names_used_in_discovery":False,
             "patient_raw_data_sent_to_ai":False,
             "direction_requires_exact_source_quote":True,
             "clinical_recommendation":False}
    (out/"MODIFIER_DISCOVERY_V31_REVALIDATED.json").write_text(
        json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload
