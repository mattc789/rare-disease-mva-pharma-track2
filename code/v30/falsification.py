"""Falsification/Kill-test engine V3.0.

The engine is asymmetric: weak or missing evidence can HOLD a hypothesis but
never becomes evidence of benefit. Only explicit contradictions create KILLs.
"""
from __future__ import annotations
from collections import defaultdict

SCHEMA="RDMVA_FALSIFICATION_V30_1"
TIER={"A":4,"B":3,"C":2,"D":1,"U":0,None:0}

def _tests(path):
    gene=path.get("gene") or "implicated gene"
    target=path.get("route_target") or "route target"
    return [
      f"Confirm that {gene} is causally relevant to the observed disease state and resolve inheritance/phase where required.",
      f"Measure the proposed molecular defect and show that changing {target} moves the disease-relevant endpoint in the predicted direction.",
      "Verify direct compound-target action and concentration-response in a relevant assay; pathway proximity alone is insufficient.",
      "Demonstrate same-compound rescue in a disease-relevant cellular/functional model with appropriate controls.",
      "Establish whether target-relevant unbound exposure is plausible before interpreting an in-vitro effect.",
      "Actively test counterevidence, genotype-dependent toxicity and mechanism-specific safety liabilities."
    ]

def _evaluate_path(path):
    kill=[];hold=[];support=[]
    pc=str(path.get("prediction_class") or "")
    if not path.get("registry_eligible",True):kill.append("not_verified_as_registry_eligible_single_ingredient_research_candidate")
    if pc in ("CONFLICTING_DIRECTION","NOT_ELIGIBLE_REGISTRY"):kill.append("deterministic_direction_or_registry_conflict")
    if path.get("scientific_hard_hold"):kill.append("scientific_route_hard_hold")
    if str(path.get("public_label_safety_state"))=="CLOSED":kill.append("public_label_safety_closed")
    if path.get("objective_match_score",0)<=0:kill.append("no_rescue_objective_mapping")
    tier=path.get("evidence_tier")
    if tier in ("A","B"):support.append("higher_research_evidence_tier")
    else:hold.append("published_or_assay_evidence_not_yet_A_or_B")
    if path.get("direct_target_assay_count",0)>0:support.append("direct_directional_target_assay")
    else:hold.append("direct_directional_target_assay_missing")
    if path.get("same_drug_functional_rescue") or path.get("functional_rescue_validated"):
        support.append("same_drug_or_route_functional_rescue")
    else:hold.append("same_drug_disease_relevant_functional_rescue_missing")
    if path.get("source_bound_mutant_model_functional_rescue"):
        support.append("source_bound_mutant_model_functional_rescue")
        hold.append("patient_variant_specific_rescue_not_demonstrated")
        if path.get("translational_context")=="NONHUMAN_ONLY":
            hold.append("human_translation_required")
    elif path.get("source_bound_primary_gene_rescue"):
        support.append("source_bound_primary_gene_product_rescue")
        hold.append("disease_functional_rescue_not_demonstrated")
    if path.get("source_bound_cross_source_bridge"):
        support.append("independent_cross_source_mechanistic_bridge")
        hold.append("same_drug_disease_model_rescue_not_demonstrated")
    if pc in ("HYPOTHESIS_ONLY","CASE_LINK_UNRESOLVED","",None):
        hold.append("patient_mechanistic_link_not_established")
    if path.get("public_label_safety_state") in ("HOLD","NOT_ASSESSED","UNASSESSED"):
        hold.append("safety_not_cleared")
    hold += list(path.get("scientific_gate_issues") or [])
    hold += list(path.get("evidence_blockers") or [])
    if kill:decision="KILLED"
    elif ("same_drug_or_route_functional_rescue" in support and
          ("higher_research_evidence_tier" in support or "direct_directional_target_assay" in support)):
        decision="CONDITIONAL_RESEARCH_LEAD"
    else:decision="HOLD"
    score=TIER.get(tier,0)*10 + path.get("direct_target_assay_count",0)*5
    if path.get("same_drug_functional_rescue") or path.get("functional_rescue_validated"):score+=15
    if path.get("source_bound_mutant_model_functional_rescue"):score+=12
    elif path.get("source_bound_primary_gene_rescue"):score+=8
    if path.get("source_bound_cross_source_bridge"):score+=9
    score+=min(10,path.get("objective_match_score",0))
    score-=min(20,len(set(hold))*2)
    if kill:score=-100
    return {
      "candidate":path.get("candidate"),"gene":path.get("gene"),
      "route_id":path.get("route_id"),"route_target":path.get("route_target"),
      "objective_id":path.get("objective_id"),"objective_type":path.get("objective_type"),
      "decision":decision,"research_priority_score":score,
      "research_priority_score_is_efficacy_probability":False,
      "supporting_factors":sorted(set(support)),
      "hold_reasons":sorted(set(str(x) for x in hold if x)),
      "kill_reasons":sorted(set(str(x) for x in kill if x)),
      "falsification_tests":_tests(path),
      "clinical_recommendation":False,
    }

def evaluate(graph):
    by=defaultdict(list)
    for p in graph.get("candidate_paths",[]):
        by[p.get("candidate")].append(_evaluate_path(p))
    candidates=[]
    for candidate,rows in by.items():
        rows.sort(key=lambda x:(
          {"CONDITIONAL_RESEARCH_LEAD":0,"HOLD":1,"KILLED":2}[x["decision"]],
          -x["research_priority_score"],x.get("route_id") or ""))
        best=rows[0]
        candidates.append({**best,"alternative_paths":rows[1:5]})
    candidates.sort(key=lambda x:(
      {"CONDITIONAL_RESEARCH_LEAD":0,"HOLD":1,"KILLED":2}[x["decision"]],
      -x["research_priority_score"],str(x.get("candidate") or "")))
    return {
      "schema":SCHEMA,"candidates":candidates,
      "decision_counts":{d:sum(x["decision"]==d for x in candidates)
                         for d in ("CONDITIONAL_RESEARCH_LEAD","HOLD","KILLED")},
      "absence_of_evidence_never_treated_as_evidence_of_absence":True,
      "clinical_recommendation":False,
    }
