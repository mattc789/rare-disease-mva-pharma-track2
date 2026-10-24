"""Mechanistic combination discovery for approved medicines.

Searches pairs/triples only after full single-drug screening. A combination is
proposed for complementary rescue-objective coverage, never assumed synergy.
"""
from __future__ import annotations
from collections import defaultdict
from itertools import combinations

SCHEMA="RDMVA_COMBINATION_ENGINE_V30_1"
OBJ_WEIGHT={"O1":5.0,"O2":3.0,"O3":1.5,"O0":0.0}
TIER_WEIGHT={"A":5.0,"B":4.0,"C":2.0,"D":1.0,"U":0.0,None:0.0}

def _candidate_profiles(graph,objectives):
    objmeta={x["objective_id"]:x for x in objectives.get("objectives",[])}
    by=defaultdict(list)
    for p in graph.get("candidate_paths",[]):
        if p.get("scientific_hard_hold"):continue
        if p.get("registry_eligible") is False:continue
        by[p.get("candidate")].append(p)
    profiles={}
    for name,paths in by.items():
        coverage={}
        for p in paths:
            oid=p.get("objective_id");meta=objmeta.get(oid) or {}
            gene=p.get("gene")
            axis=p.get("objective_axis") or meta.get("objective_axis") or p.get("objective_type")
            axis_id=f"{gene}:{axis}"
            quality=TIER_WEIGHT.get(p.get("evidence_tier"),0)
            quality+=min(2.0,float(p.get("direct_target_assay_count") or 0))
            if p.get("same_drug_functional_rescue") or p.get("functional_rescue_validated"):quality+=4
            quality+=min(2.0,float(p.get("objective_match_score") or 0)/3.0)
            rec={"objective_id":oid,"objective_type":p.get("objective_type"),
                 "objective_axis":axis,"axis_id":axis_id,
                 "gene":gene,"priority":meta.get("priority_tier","O3"),
                 "quality":quality,"route_id":p.get("route_id"),
                 "route_target":p.get("route_target"),
                 "desired_action":p.get("route_desired_action"),
                 "evidence_tier":p.get("evidence_tier"),
                 "source_channels":p.get("source_channels") or [],
                 "public_label_safety_state":p.get("public_label_safety_state")}
            weighted=OBJ_WEIGHT.get(rec["priority"],1.0)*quality
            previous=coverage.get(axis_id)
            previous_weighted=(OBJ_WEIGHT.get(previous["priority"],1.0)*previous["quality"]) if previous else -1
            if oid and (previous is None or weighted>previous_weighted):coverage[axis_id]=rec
        if coverage:
            profiles[name]={"candidate":name,"coverage":coverage,
                            "rescue_axes":set(coverage),"paths":paths}
    return profiles,objmeta

def _conflicts(names,profiles):
    actions=defaultdict(set);reasons=[]
    for name in names:
        for p in profiles[name]["paths"]:
            target=str(p.get("route_target") or "").upper()
            act=str(p.get("route_desired_action") or "")
            if target and act in ("increase","decrease"):actions[target].add(act)
            if p.get("public_label_safety_state")=="CLOSED":
                reasons.append(f"{name}:safety_closed")
    for target,acts in actions.items():
        if "increase" in acts and "decrease" in acts:
            reasons.append("opposing_actions_same_target:"+target)
    return sorted(set(reasons))

def _score(names,profiles,objmeta):
    covered={}
    component_scores={}
    for name in names:
        subtotal=0.0
        for axis_id,rec in profiles[name]["coverage"].items():
            w=OBJ_WEIGHT.get(rec.get("priority","O3"),1.0)
            val=w*max(0.25,rec["quality"])
            subtotal+=val
            if axis_id not in covered or val>covered[axis_id][0]:
                covered[axis_id]=(val,name,rec)
        component_scores[name]=round(subtotal,3)
    total=sum(v[0] for v in covered.values())
    # Penalize redundancy: adding a drug that contributes no unique rescue
    # axis should not produce a high-ranked combination.
    unique={name:0 for name in names}
    for axis_id,(v,name,rec) in covered.items():unique[name]+=1
    redundancy=sum(unique[n]==0 for n in names)*4.0
    return total-redundancy,covered,component_scores,unique

def discover(graph,objectives,max_components=3,shortlist=50,max_results=30):
    profiles,objmeta=_candidate_profiles(graph,objectives)
    singles=[]
    for name,p in profiles.items():
        s,c,cs,u=_score((name,),profiles,objmeta)
        singles.append((s,name))
    singles.sort(reverse=True)
    pool=[name for _,name in singles[:shortlist]]
    combos=[]
    for size in range(2,max_components+1):
        for names in combinations(pool,size):
            conflicts=_conflicts(names,profiles)
            if conflicts:continue
            score,covered,component_scores,unique=_score(names,profiles,objmeta)
            # Every component must add at least one objective not better covered
            # by the other components.
            if any(unique[n]==0 for n in names):continue
            priorities=[rec.get("priority") for _,_,rec in covered.values()]
            if "O1" not in priorities:continue
            rows=[]
            for axis_id,(value,name,rec) in covered.items():
                rows.append({**rec,"best_component":name,"weighted_value":round(value,3)})
            rows.sort(key=lambda x:(x.get("priority"),x.get("gene"),x.get("objective_type")))
            combos.append({
              "components":list(names),"component_count":size,
              "combination_priority_score":round(score,3),
              "score_is_clinical_benefit_probability":False,
              "covered_objective_count":len(covered),
              "covered_rescue_axis_count":len(covered),
              "covered_objectives":rows,
              "component_objective_contributions":unique,
              "component_rescue_axis_contributions":unique,
              "component_scores":component_scores,
              "combination_logic":"complementary_mechanistic_rescue_axis_coverage",
              "synergy_claimed":False,
              "requires_combination_specific_interaction_and_safety_review":True,
              "clinical_recommendation":False,
            })
    combos.sort(key=lambda x:(-x["covered_objective_count"],
                              -x["combination_priority_score"],
                              x["component_count"],x["components"]))
    # Remove dominated combinations with identical coverage where a smaller set exists.
    kept=[]
    for c in combos:
        cov=frozenset(x["axis_id"] for x in c["covered_objectives"])
        if any(frozenset(x["axis_id"] for x in k["covered_objectives"])==cov
               and k["component_count"]<c["component_count"] for k in kept):continue
        kept.append(c)
        if len(kept)>=max_results:break
    return {"schema":SCHEMA,"candidate_pool_size":len(pool),
            "pairs_and_triples_evaluated":True,"max_components":max_components,
            "combinations":kept,"combination_count":len(kept),
            "market_approved_single_agent_filter_required":True,
            "synergy_not_assumed":True,"clinical_recommendation":False}
