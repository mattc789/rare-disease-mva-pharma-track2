"""V3 evidence ranker for the universal approved-drug screen."""
from __future__ import annotations
from collections import defaultdict

SCHEMA="RDMVA_EVIDENCE_RANKER_V30_1"

def prioritize(ledger,routes):
    rmap={r["route_id"]:r for r in routes}
    pairs=[]
    for entity in ledger:
        if entity.get("status")!="SIGNAL":continue
        groups=defaultdict(list)
        for signal in entity.get("signals",[]):
            rid=signal.get("route_id");route=rmap.get(rid)
            if not route or route.get("state")=="hold":continue
            groups[rid].append(signal)
        for rid,sigs in groups.items():
            route=rmap[rid];sci=route.get("scientific_gate") or {}
            net=[int(s.get("depth",99)) for s in sigs if s.get("provider")=="ChEMBL_MoA_SIGNOR"]
            shortest=min(net or [99])
            route_hops=int(route.get("network_hops") or 0)
            effective_hops=(shortest+route_hops) if shortest<99 else None
            direct_moa=(shortest==0 and route_hops==0)
            assay_direct=any(
              s.get("provider")=="ChEMBL_measured_activity" and
              s.get("value_nM") is not None and float(s["value_nM"])<=1000 and
              s.get("direct_target_assay_verified") is True for s in sigs)
            real_rescue=any(
              s.get("provider")=="patient_disease_functional_rescue" and
              all(s.get(k) is True for k in ("same_drug","functional_endpoint_verified",
                  "disease_relevant_model","exposure_supported","source_independently_verified"))
              for s in sigs)
            family=str(route.get("family") or "")
            direct_primary=(family in ("v30_direct_gene","direct_gene") and
                            str(route.get("target_concept") or "").upper()==
                            str(route.get("case_gene") or "").upper())
            route_direction_ok=route.get("desired_action") in ("increase","decrease")
            context_verified=(sci.get("causal_direction_validated_for_disease") is True or
                              (direct_primary and route_direction_ok and assay_direct))
            if real_rescue and context_verified:tier="A"
            elif assay_direct and (context_verified or direct_primary):tier="B"
            elif direct_moa or (effective_hops==1 and route_direction_ok):tier="C"
            else:tier="D"
            route_source_rescue=bool(route.get("functional_rescue_observed_in_source"))
            score=(
              (8 if real_rescue and context_verified else 0)+
              (4 if assay_direct else 0)+
              (3 if route_source_rescue else 0)+
              ({0:2.0,1:0.8,2:0.35,3:0.1}.get(effective_hops,0))+
              min(2.0,float(route.get("score") or 0)/10.0))
            providers=sorted({s.get("provider") for s in sigs if s.get("provider")})
            reasons=[]
            if not real_rescue:reasons.append("no_same_drug_disease_functional_rescue")
            if not assay_direct:reasons.append("no_verified_direct_target_modulation")
            if not direct_primary and not sci.get("causal_direction_validated_for_disease"):
                reasons.append("target_to_disease_direction_not_independently_verified")
            if effective_hops is not None and effective_hops>=2:
                reasons.append("multi_hop_context_sensitive_inference")
            if route_source_rescue:
                reasons.append("source_bound_modifier_rescue_support")
            if route.get("case_gate_reasons"):reasons.extend(route["case_gate_reasons"])
            pairs.append({
              "schema":SCHEMA,
              "ingredient":entity["canonical_label"],"entity_id":entity["entity_id"],
              "scientific_group_id":entity.get("scientific_group_id"),
              "member_labels":entity.get("member_labels",[entity["canonical_label"]]),
              "route_id":rid,"route":route,"tier":tier,"score":round(score,4),
              "providers":providers,"signals":sigs,
              "combo_only":bool(entity.get("combo_only")),
              "global_approved_as_single":bool(entity.get("global_approved_as_single") or
                                                entity.get("solo_approved") or
                                                entity.get("biologic_licensed")),
              "modality":entity.get("modality"),
              "publication_strength":"experimental" if tier in ("A","B") else "hypothesis_only",
              "evidence_blockers":sorted(set(reasons)),
              "shortest_network_hops":effective_hops,
              "drug_to_route_hops":shortest if shortest<99 else None,
              "route_to_case_hops":route_hops,
              "scientific_gate":sci,
              "rescue_objective_ids":route.get("rescue_objective_ids") or
                                     ([route.get("rescue_objective_id")] if route.get("rescue_objective_id") else []),
              "rescue_objective_types":route.get("rescue_objective_types") or
                                       ([route.get("rescue_objective_type")] if route.get("rescue_objective_type") else []),
              "promotion_eligible_validated_lead":tier=="A",
              "clinical_recommendation":False,
            })
    pairs.sort(key=lambda x:({"A":0,"B":1,"C":2,"D":3}[x["tier"]],
                              -x["score"],x["ingredient"],x["route_id"]))
    return pairs
