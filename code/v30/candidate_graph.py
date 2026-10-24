"""Causal candidate graph V3.0.

Connects local patient-state mechanisms to rescue objectives, signed routes and
drug evidence. The graph records support; it does not infer clinical efficacy.
"""
from __future__ import annotations
from collections import defaultdict

SCHEMA="RDMVA_CANDIDATE_GRAPH_V30_1"

def _route_gene(route,known_genes):
    gene=str(route.get("case_gene") or "").upper()
    if gene:return gene
    text=" ".join(str(route.get(k) or "") for k in ("target_concept","mechanism_concept")).upper()
    hits=[g for g in known_genes if g and g in text]
    return hits[0] if len(hits)==1 else None

def _objective_match(obj,route):
    if obj.get("gene")!=route.get("_case_gene"):return 0,[]
    # Native V3 routes retain the exact rescue objective(s) that generated
    # them. Preserve this provenance instead of heuristically letting one
    # route appear to solve unrelated objectives from the same gene.
    native_ids=set(route.get("rescue_objective_ids") or [])
    if route.get("rescue_objective_id"):native_ids.add(route["rescue_objective_id"])
    if native_ids:
        if obj.get("objective_id") not in native_ids:return 0,[]
        score=8;why=["native_route_derived_from_objective"]
        ov=str(obj.get("variable") or "").lower()
        rv=str(route.get("variable") or "").lower()
        if ov and rv and (ov==rv or ov in rv or rv in ov):
            score+=2;why.append("same_or_related_variable")
        oa=obj.get("desired_action");ra=route.get("desired_action")
        if oa in ("increase","decrease") and ra==oa:
            score+=1;why.append("direction_matches")
        elif oa=="restore" and ra in ("increase","decrease"):
            score+=1;why.append("directional_compensation_for_restore_objective")
        return score,why

    family=str(route.get("family") or "")
    if family.startswith(("v30_","v31_")) or route.get("objective_binding_state"):
        # Native routes without an explicit rescue-objective binding are
        # exploratory network/literature observations only and must not be
        # heuristically attached to another objective.
        return 0,[]

    # Legacy-route fallback used only for compatibility/replay.
    score=2;why=["same_case_gene"]
    ov=str(obj.get("variable") or "").lower()
    rv=str(route.get("variable") or "").lower()
    if ov and rv and (ov==rv or ov in rv or rv in ov):
        score+=3;why.append("same_or_related_variable")
    oa=obj.get("desired_action");ra=route.get("desired_action")
    if oa in ("increase","decrease") and ra==oa:
        score+=2;why.append("direction_matches")
    family=str(route.get("family") or "")
    if oa=="restore" and family in ("signed_upstream","signed_downstream","literature_discovered"):
        score+=2;why.append("signed_network_compensation")
    if obj.get("target_layer")=="primary_gene" and family=="direct_gene":
        score+=3;why.append("direct_primary_gene_route")
    if obj.get("objective_type")=="restore_productive_splicing" and family=="splicing":
        score+=4;why.append("splicing_route")
    return (score,why) if score>=4 else (0,[])

def build(objectives,routes,pairs,predictions=None):
    known={x.get("gene") for x in objectives.get("objectives",[]) if x.get("gene")}
    route_by={}
    for r0 in routes or []:
        r=dict(r0)
        r["_case_gene"]=_route_gene(r,known)
        route_by[r.get("route_id")]=r
    pred_idx={(x.get("candidate"),x.get("route_id")):x for x in (predictions or [])}
    obj_by_gene=defaultdict(list)
    for o in objectives.get("objectives",[]):
        if o.get("drug_search_eligible"):obj_by_gene[o.get("gene")].append(o)
    paths=[];nodes={};edges=[]
    for p in pairs or []:
        rid=p.get("route_id");route=dict(p.get("route") or route_by.get(rid) or {})
        route.setdefault("route_id",rid)
        route["_case_gene"]=_route_gene(route,known)
        gene=route.get("_case_gene")
        if not gene:continue
        scored=[]
        for o in obj_by_gene.get(gene,[]):
            sc,why=_objective_match(o,route)
            if sc:scored.append((sc,o,why))
        if not scored:continue
        scored.sort(key=lambda x:(-x[0],x[1]["priority_tier"],x[1]["objective_id"]))
        prediction=pred_idx.get((p.get("ingredient"),rid)) or {}
        pf=p.get("patient_first_prediction") or {}
        expert=p.get("expert_triage") or {}
        functional=p.get("functional_evidence") or {}
        signals=p.get("signals") or []
        direct=sum(bool(s.get("direct_target_assay_verified")) for s in signals)
        providers=sorted({str(s.get("provider")) for s in signals if s.get("provider")})
        for match_score,o,why in scored[:3]:
            path={
              "candidate":p.get("ingredient"),
              "member_labels":list(p.get("member_labels") or prediction.get("member_labels") or [p.get("ingredient")]),
              "scientific_group_id":p.get("scientific_group_id"),
              "gene":gene,"route_id":rid,
              "route_family":route.get("family"),
              "route_target":route.get("target_concept"),
              "route_desired_action":route.get("desired_action"),
              "route_variable":route.get("variable"),
              "objective_id":o["objective_id"],
              "objective_type":o["objective_type"],
              "objective_axis":o.get("objective_axis") or o["objective_type"],
              "objective_priority":o["priority_tier"],
              "objective_match_score":match_score,
              "objective_match_reasons":why,
              "evidence_tier":p.get("tier") or pf.get("research_evidence_tier"),
              "shortest_network_hops":p.get("shortest_network_hops"),
              "drug_to_route_hops":p.get("drug_to_route_hops"),
              "route_to_case_hops":p.get("route_to_case_hops"),
              "source_bound_modifier_rescue":bool(
                  (route.get("scientific_gate") or {}).get("functional_rescue_validated")),
              "source_bound_direction_verified":bool(route.get("source_bound_direction_verified")),
              "source_bound_primary_gene_rescue":bool(
                  p.get("source_bound_primary_gene_rescue") or route.get("source_bound_primary_gene_rescue")),
              "source_bound_mutant_model_functional_rescue":bool(
                  p.get("source_bound_mutant_model_functional_rescue") or
                  route.get("source_bound_mutant_model_functional_rescue")),
              "source_bound_cross_source_bridge":bool(
                  p.get("source_bound_cross_source_bridge") or
                  route.get("source_bound_cross_source_bridge")),
              "direct_source_rescue_evidence_class":p.get("direct_source_rescue_evidence_class"),
              "source_bound_cross_source_bridge":bool(
                  p.get("source_bound_cross_source_bridge") or route.get("source_bound_cross_source_bridge")),
              "cross_source_bridge_rescue_sources":list(route.get("rescue_sources") or []),
              "cross_source_bridge_drug_target_source":route.get("drug_target_source"),
              "cross_source_bridge_target_aliases":list(route.get("target_aliases") or []),
              "translational_context":p.get("translational_context") or route.get("translational_context"),
              "polypharmacology":p.get("polypharmacology") or {},
              "investigation_bucket":p.get("expert_investigation_bucket") or expert.get("investigation_bucket"),
              "prediction_class":prediction.get("prediction_class") or pf.get("prediction_class"),
              "mechanistic_priority_score_audit_only":prediction.get("mechanistic_priority_score") or pf.get("mechanistic_priority_score"),
              "source_channels":providers,
              "direct_target_assay_count":direct or int(expert.get("direct_functional_target_assay_count") or 0),
              "same_drug_functional_rescue":bool(functional.get("same_drug_functional_rescue") or
                                                   functional.get("same_drug_disease_functional_rescue")),
              "functional_rescue_validated":bool((route.get("scientific_gate") or {}).get("functional_rescue_validated")),
              "scientific_hard_hold":bool((route.get("scientific_gate") or {}).get("hard_hold")),
              "scientific_gate_issues":list((route.get("scientific_gate") or {}).get("issues") or []),
              "evidence_blockers":list(p.get("evidence_blockers") or expert.get("hypothesis_blockers") or []),
              "polypharmacology":dict(p.get("polypharmacology") or {}),
              "variant_specific_computational_context":dict(o.get("variant_specific_computational_context") or {}),
              "public_label_safety_state":prediction.get("public_label_safety_state") or expert.get("label_state") or "NOT_ASSESSED",
              "public_label_safety_flags":list(prediction.get("public_label_safety_flags") or []),
              "market_approved_single_agent":bool(p.get("global_approved_as_single",not bool(p.get("combo_only",False)))),
              "registry_eligible":bool((prediction.get("registry") or {}).get(
                  "eligible_for_research_search",not bool(p.get("combo_only",False)))),
              "clinical_recommendation":False,
            }
            paths.append(path)
            cid="DRUG:"+str(path["candidate"]);ridn="ROUTE:"+str(rid);oid="OBJECTIVE:"+o["objective_id"]
            nodes[cid]={"id":cid,"type":"drug","label":path["candidate"]}
            nodes[ridn]={"id":ridn,"type":"route","label":route.get("target_concept"),"gene":gene}
            nodes[oid]={"id":oid,"type":"objective","label":o["objective_type"],"gene":gene}
            edges += [
              {"source":cid,"target":ridn,"type":"supports_route","evidence_tier":path["evidence_tier"]},
              {"source":ridn,"target":oid,"type":"maps_to_rescue_objective","match_score":match_score},
            ]
    paths.sort(key=lambda x:(
      {"A":0,"B":1,"C":2,"D":3,"U":4,None:5}.get(x.get("evidence_tier"),5),
      -x.get("direct_target_assay_count",0),
      -x.get("objective_match_score",0),
      str(x.get("candidate") or "")))
    return {"schema":SCHEMA,"nodes":list(nodes.values()),"edges":edges,
            "candidate_paths":paths,"clinical_efficacy_inferred":False}
