"""Track-2 hypothesis classifier.

Maps mechanistic research evidence to the challenge's actual scope:
market-approved medications that are plausible candidates for further
investigation. This is deliberately distinct from a validated rescue lead.
"""
from __future__ import annotations
from collections import defaultdict

SCHEMA="RDMVA_CHALLENGE_HYPOTHESES_V30_1"
TIER={"A":4,"B":3,"C":2,"D":1,"U":0,None:0}
PRIORITY={"O1":3,"O2":2,"O3":1,"O0":0}

def _path_strength(p):
    if not p.get("market_approved_single_agent",False):return -999,[]
    if p.get("scientific_hard_hold"):return -999,[]
    reasons=[]
    score=0.0
    tier=p.get("evidence_tier")
    score+=TIER.get(tier,0)*3
    if tier in ("A","B"):reasons.append("higher_research_evidence")
    elif tier=="C":reasons.append("direct_moa_or_close_mechanistic_evidence")
    elif tier=="D":reasons.append("indirect_mechanistic_evidence_only")
    match=float(p.get("objective_match_score") or 0)
    score+=match
    if match>=5:reasons.append("strong_rescue_objective_mapping")
    direct=int(p.get("direct_target_assay_count") or 0)
    if direct:
        score+=5+min(3,direct);reasons.append("direct_directional_target_assay")
    hops=p.get("shortest_network_hops")
    if hops==0:
        score+=4;reasons.append("direct_drug_target_route")
    elif hops==1:
        score+=2;reasons.append("one_hop_signed_mechanistic_route")
    elif isinstance(hops,int) and hops>=2:
        score-=min(4,hops);reasons.append("multi_hop_route_uncertainty")
    if p.get("same_drug_functional_rescue") or p.get("functional_rescue_validated"):
        score+=8;reasons.append("functional_rescue_evidence")
    source_primary=bool(p.get("source_bound_primary_gene_rescue"))
    source_mutant=bool(p.get("source_bound_mutant_model_functional_rescue"))
    if source_mutant:
        score+=15;reasons.append("source_bound_mutant_model_functional_rescue")
        if p.get("translational_context")=="NONHUMAN_ONLY":
            score-=2;reasons.append("nonhuman_model_requires_translation")
    elif source_primary:
        score+=6;reasons.append("source_bound_primary_gene_product_rescue")
    cross_bridge=bool(p.get("source_bound_cross_source_bridge"))
    if cross_bridge:
        score+=9;reasons.append("independent_rescue_target_plus_direct_approved_drug_modulation")
    priority=PRIORITY.get(p.get("objective_priority"),1)
    score+=priority*2
    blockers=set(p.get("evidence_blockers") or [])
    score-=min(8,len(blockers))
    if p.get("prediction_class")=="CONFLICTING_DIRECTION":return -999,[]
    return score,reasons

def classify_singles(graph,max_results=30):
    by=defaultdict(list)
    for p in graph.get("candidate_paths",[]):
        s,reasons=_path_strength(p)
        if s<=-900:continue
        by[p.get("candidate")].append((s,p,reasons))
    rows=[]
    for name,paths in by.items():
        paths.sort(key=lambda x:(-x[0],str(x[1].get("route_id") or "")))
        best=paths[0]
        score,p,reasons=best
        def _axis_supported(path):
            return bool(
              path.get("source_bound_mutant_model_functional_rescue") or
              path.get("source_bound_primary_gene_rescue") or
              path.get("source_bound_cross_source_bridge") or
              path.get("same_drug_functional_rescue") or
              path.get("functional_rescue_validated") or
              int(path.get("direct_target_assay_count") or 0)>0 or
              path.get("evidence_tier") in ("A","B"))
        supportive_paths=[x for x in paths if x[0]>0 and _axis_supported(x[1])]
        axes=sorted({f"{x[1].get('gene')}:{x[1].get('objective_axis') or x[1].get('objective_type')}"
                     for x in supportive_paths})
        coverage_bonus=min(6.0,max(0,len(axes)-1)*1.5)
        score+=coverage_bonus
        poly=dict(p.get("polypharmacology") or {})
        if poly.get("polypharmacology_present"):
            score+=0.5;reasons.append("known_multi_target_pharmacology")
        if int(poly.get("compatible_route_target_mechanism_count") or 0)>0:
            score+=1.0;reasons.append("chembl_signed_mechanism_matches_route")
        if len(axes)>1:reasons.append("multiple_rescue_axes_covered")
        tier=p.get("evidence_tier")
        direct=int(p.get("direct_target_assay_count") or 0)
        hops=p.get("shortest_network_hops")
        drug_hops=p.get("drug_to_route_hops")
        route_hops=int(p.get("route_to_case_hops") or 0)
        match=int(p.get("objective_match_score") or 0)
        # Track-2 hypothesis does NOT require demonstrated disease rescue.
        direct_family=p.get("route_family") in ("v30_direct_gene","direct_gene")
        has_rescue=bool(p.get("same_drug_functional_rescue") or p.get("functional_rescue_validated"))
        direct_drug_route=(drug_hops==0)
        total_close=(isinstance(hops,int) and hops<=1)
        # STRONG needs high evidence and a very short chain. A C-tier network
        # hypothesis is strong only when it is truly 0-hop to the patient
        # rescue route; a one-hop route remains PLAUSIBLE pending validation.
        strong=(match>=5 and (
                  tier in ("A","B") and total_close and
                    (direct>0 or has_rescue or direct_drug_route or direct_family)
                  or tier=="C" and hops==0 and (direct_family or direct>0)
                ))
        # PLAUSIBLE may include one causal hop if the drug directly modulates
        # the route target. Total causal distance >=2 remains exploratory.
        source_direct=bool(p.get("source_bound_primary_gene_rescue") or
                           p.get("source_bound_mutant_model_functional_rescue"))
        cross_bridge=bool(p.get("source_bound_cross_source_bridge"))
        plausible=((match>=4 and tier in ("A","B","C") and total_close and
                    (direct>0 or has_rescue or direct_drug_route or direct_family))
                   or (match>=5 and tier in ("B","C") and source_direct)
                   or (match>=5 and tier=="C" and cross_bridge and direct>0))
        if strong:status="STRONG_CHALLENGE_HYPOTHESIS"
        elif plausible:status="PLAUSIBLE_CHALLENGE_HYPOTHESIS"
        else:status="EXPLORATORY_ONLY"
        rows.append({
          "candidate":name,"status":status,"challenge_priority_score":round(score,3),
          "score_is_clinical_efficacy_probability":False,
          "gene":p.get("gene"),"objective_id":p.get("objective_id"),
          "objective_type":p.get("objective_type"),"objective_priority":p.get("objective_priority"),
          "route_id":p.get("route_id"),"route_target":p.get("route_target"),
          "route_desired_action":p.get("route_desired_action"),
          "evidence_tier":tier,"shortest_network_hops":hops,
          "drug_to_route_hops":drug_hops,"route_to_case_hops":route_hops,
          "direct_target_assay_count":direct,
          "same_drug_functional_rescue":bool(p.get("same_drug_functional_rescue") or p.get("functional_rescue_validated")),
          "source_bound_primary_gene_rescue":bool(p.get("source_bound_primary_gene_rescue")),
          "source_bound_mutant_model_functional_rescue":bool(p.get("source_bound_mutant_model_functional_rescue")),
          "source_bound_cross_source_bridge":bool(p.get("source_bound_cross_source_bridge")),
          "direct_source_rescue_evidence_class":p.get("direct_source_rescue_evidence_class"),
          "cross_source_bridge_rescue_sources":list(p.get("cross_source_bridge_rescue_sources") or []),
          "cross_source_bridge_drug_target_source":p.get("cross_source_bridge_drug_target_source"),
          "cross_source_bridge_target_aliases":list(p.get("cross_source_bridge_target_aliases") or []),
          "translational_context":p.get("translational_context"),
          "covered_rescue_axes":axes,
          "covered_rescue_axis_count":len(axes),
          "polypharmacology":poly,
          "why_plausible":reasons,
          "main_blockers":list(p.get("evidence_blockers") or p.get("scientific_gate_issues") or []),
          "market_approved_single_agent":True,
          "hypothesis_for_follow_up_not_evidence_it_works":True,
          "clinical_recommendation":False,
        })
    order={"STRONG_CHALLENGE_HYPOTHESIS":0,"PLAUSIBLE_CHALLENGE_HYPOTHESIS":1,"EXPLORATORY_ONLY":2}
    rows.sort(key=lambda x:(order[x["status"]],-x["challenge_priority_score"],x["candidate"]))
    return rows[:max_results]

def classify_combinations(combination_result,single_rows,max_results=20):
    smap={x["candidate"]:x for x in single_rows}
    out=[]
    for c in combination_result.get("combinations",[]):
        components=c.get("components") or []
        comps=[smap.get(x) for x in components]
        if any(x is None for x in comps):continue
        if any(x["status"]=="EXPLORATORY_ONLY" for x in comps):continue
        strong_count=sum(x["status"]=="STRONG_CHALLENGE_HYPOTHESIS" for x in comps)
        status="STRONG_COMBINATION_HYPOTHESIS" if strong_count==len(comps) else "PLAUSIBLE_COMBINATION_HYPOTHESIS"
        out.append({
          **c,"status":status,
          "component_hypotheses":[{"candidate":x["candidate"],"status":x["status"],
                                   "gene":x["gene"],"objective_type":x["objective_type"],
                                   "evidence_tier":x["evidence_tier"]} for x in comps],
          "challenge_rationale":"Each approved component contributes unique mechanistic rescue-objective coverage.",
          "combination_effect_demonstrated":False,
          "synergy_claimed":False,
          "requires_combination_specific_safety_interaction_review":True,
          "hypothesis_for_follow_up_not_evidence_it_works":True,
          "clinical_recommendation":False,
        })
    order={"STRONG_COMBINATION_HYPOTHESIS":0,"PLAUSIBLE_COMBINATION_HYPOTHESIS":1}
    out.sort(key=lambda x:(order[x["status"]],-x["covered_objective_count"],
                           -x["combination_priority_score"],x["component_count"]))
    return out[:max_results]

def build(graph,objectives,combination_result):
    # Keep a broad internal pool so a useful combination component is not lost
    # merely because it sits outside the top standalone shortlist.
    singles_all=classify_singles(graph,max_results=500)
    combos=classify_combinations(combination_result,singles_all)
    singles=singles_all[:100]
    submission_singles=[x for x in singles_all if x["status"]!="EXPLORATORY_ONLY"][:30]
    return {
      "schema":SCHEMA,
      "track2_scope":"market-approved medications plausible for further investigation",
      "single_agent_hypotheses":singles,
      "submission_level_single_agents":submission_singles,
      "combination_hypotheses":combos,
      "has_submission_level_hypothesis":bool(submission_singles or combos),
      "validated_disease_rescue_required_for_hypothesis":False,
      "approved_medication_required":True,
      "clinical_recommendation":False,
    }


def apply_source_bound_evidence(track2,attack_summary):
    """Annotate/re-rank Track-2 hypotheses from exact-quote public evidence.

    Lack of source-bound support does not kill a mechanistic hypothesis.
    Explicit counterevidence creates a submission hold until resolved.
    """
    attacks={x.get("candidate"):x for x in attack_summary.get("results",[])
             if x.get("candidate")}
    bonus_map={
      "P2_INDEPENDENT_DIRECT_AND_RESCUE":20.0,
      "P1_SINGLE_SOURCE_DIRECT_AND_RESCUE":10.0,
      "P1_PARTIAL_SOURCE_BOUND_SUPPORT":5.0,
      "P0_NO_SOURCE_BOUND_SUPPORT":0.0,
      "P0_CONFLICTED_SOURCE_BOUND_EVIDENCE":-25.0,
      "P0_COUNTEREVIDENCE_ONLY":-25.0,
    }
    def annotate(x0):
        x=dict(x0);a=attacks.get(x.get("candidate"))
        if not a:
            if x.get("source_bound_cross_source_bridge"):
                rescue=list(x.get("cross_source_bridge_rescue_sources") or [])
                drug=x.get("cross_source_bridge_drug_target_source")
                x["source_bound_evidence_state"]="P2_CROSS_SOURCE_MECHANISTIC_SUPPORT"
                x["source_bound_public_evidence"]={
                  "overall_assessment":"CROSS_SOURCE_MECHANISTIC_SUPPORT",
                  "direct_source_count":1 if drug else 0,
                  "rescue_source_count":len(rescue),
                  "independent_support_source_count":(1 if drug else 0)+len(rescue),
                  "counterevidence_source_count":0,
                  "missing_evidence":["same_drug_disease_model_rescue_not_demonstrated"],
                  "source_metadata_verified":True,
                  "raw_source_text_redistributed":False,
                  "same_article_contains_full_bridge":False}
            else:
                x["source_bound_evidence_state"]="NOT_ATTACKED"
            x["challenge_priority_score_source_adjusted"]=x.get("challenge_priority_score")
            x["submission_hold_from_source_evidence"]=False
            return x
        level=a.get("evidence_level","UNKNOWN")
        has_counter=bool(a.get("counterevidence_sources"))
        direct_rescue_context=bool(
            a.get("rescue_sources") and
            a.get("overall_assessment")=="DIRECT_RESCUE_EVIDENCE_FOUND")
        bridge_context=bool(x.get("source_bound_cross_source_bridge") and has_counter and
                            a.get("overall_assessment")!="COUNTEREVIDENCE_DOMINANT")
        context_conflict=bool(has_counter and (direct_rescue_context or bridge_context))
        conflict=bool(
            a.get("overall_assessment")=="COUNTEREVIDENCE_DOMINANT" or
            (has_counter and not (direct_rescue_context or bridge_context)))
        bridge_verified=bool(
            x.get("source_bound_cross_source_bridge") and
            level=="P0_NO_SOURCE_BOUND_SUPPORT" and not has_counter)
        if bridge_verified:
            rescue=list(x.get("cross_source_bridge_rescue_sources") or [])
            drug=x.get("cross_source_bridge_drug_target_source")
            x["open_access_attack_evidence_state"]=level
            x["source_bound_evidence_state"]="P2_CROSS_SOURCE_MECHANISTIC_SUPPORT"
            x["source_bound_public_evidence"]={
              "overall_assessment":"CROSS_SOURCE_MECHANISTIC_SUPPORT",
              "direct_source_count":1 if drug else 0,
              "rescue_source_count":len(rescue),
              "independent_support_source_count":(1 if drug else 0)+len(rescue),
              "counterevidence_source_count":0,
              "missing_evidence":["same_drug_disease_model_rescue_not_demonstrated"],
              "source_metadata_verified":True,
              "raw_source_text_redistributed":False,
              "same_article_contains_full_bridge":False}
        else:
            x["source_bound_evidence_state"]=level
            x["source_bound_public_evidence"]={
              "overall_assessment":a.get("overall_assessment"),
              "direct_source_count":len(a.get("direct_sources") or []),
              "rescue_source_count":len(a.get("rescue_sources") or []),
              "independent_support_source_count":len(a.get("independent_support_sources") or []),
              "counterevidence_source_count":len(a.get("counterevidence_sources") or []),
              "missing_evidence":a.get("missing_evidence") or [],
              "quotes_verified_exactly":bool(a.get("quotes_verified_exactly")),
            }
        x["submission_hold_from_source_evidence"]=conflict
        x["context_conflict_requires_resolution"]=context_conflict
        base=float(x.get("challenge_priority_score") or 0)
        bonus=bonus_map.get(level,0.0)
        # A direct-source hypothesis already received credit for its discovery
        # source in the base score. Do not double-count the same source merely
        # because the literature attack rediscovered it.
        if x.get("direct_source_rescue_evidence_class") and bonus>0:
            bonus=0.0
        if context_conflict:
            # Direct same-drug mutant-model rescue survives as a Track-2
            # hypothesis, but competing model/genotype evidence is a material
            # uncertainty and receives a penalty rather than being ignored.
            bonus=-5.0
        x["challenge_priority_score_source_adjusted"]=round(base+bonus,3)
        return x

    singles=[annotate(x) for x in track2.get("single_agent_hypotheses",[])]
    submission=[annotate(x) for x in track2.get("submission_level_single_agents",[])]
    submission=[x for x in submission if not x.get("submission_hold_from_source_evidence")]
    submission.sort(key=lambda x:(
        -float(x.get("challenge_priority_score_source_adjusted") or 0),
        -float(x.get("challenge_priority_score") or 0),
        x.get("candidate") or ""))
    combos=[]
    single_by={x["candidate"]:x for x in singles}
    for c0 in track2.get("combination_hypotheses",[]):
        c=dict(c0);comp_evidence={}
        hold=False;bonus=0.0
        for name in c.get("components") or []:
            sx=single_by.get(name) or {}
            comp_evidence[name]={
              "source_bound_evidence_state":sx.get("source_bound_evidence_state","NOT_ATTACKED"),
              "submission_hold":bool(sx.get("submission_hold_from_source_evidence")),
            }
            hold=hold or bool(sx.get("submission_hold_from_source_evidence"))
            if name in attacks:
                bonus+=bonus_map.get(attacks[name].get("evidence_level"),0.0)
        c["component_source_bound_evidence"]=comp_evidence
        c["submission_hold_from_source_evidence"]=hold
        c["combination_priority_score_source_adjusted"]=round(
            float(c.get("combination_priority_score") or 0)+bonus,3)
        combos.append(c)
    combos.sort(key=lambda x:(bool(x.get("submission_hold_from_source_evidence")),
                              -float(x.get("combination_priority_score_source_adjusted") or 0),
                              -int(x.get("covered_rescue_axis_count") or x.get("covered_objective_count") or 0),
                              x.get("component_count") or 9))
    result=dict(track2)
    result["single_agent_hypotheses"]=singles
    result["submission_level_single_agents"]=submission
    result["combination_hypotheses"]=combos
    result["source_bound_evidence_applied"]=True
    result["has_submission_level_hypothesis"]=bool(
        submission or any(not c.get("submission_hold_from_source_evidence") for c in combos))
    return result
