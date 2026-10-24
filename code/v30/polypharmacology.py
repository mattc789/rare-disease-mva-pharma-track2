"""Whole-drug polypharmacology reconstruction for V3.1.

Uses local ChEMBL mechanism records for every approved entity and annotates
route-drug pairs with the drug's known multi-target profile. ChEMBL bulk data
remain external and are not redistributed in the source package.
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import json
from v2810 import screen

SCHEMA="RDMVA_POLYPHARMACOLOGY_V31_1"

def _accession_gene_map():
    out={}
    for t in screen.load(screen.CHEMBL_TARGETS):
        if t.get("organism")!="Homo sapiens":continue
        comps=t.get("target_components") or []
        accs={x.get("accession") for x in comps if x.get("accession")}
        if len(accs)!=1:continue
        acc=next(iter(accs));genes=[]
        for c in comps:
            for s in c.get("target_component_synonyms") or []:
                if s.get("syn_type") in ("GENE_SYMBOL","GENE_SYMBOL_OTHER"):
                    g=str(s.get("component_synonym") or "").strip().upper()
                    if g and " " not in g and "=" not in g and "{" not in g:
                        genes.append(g)
        out[acc]=genes[0] if genes else t.get("pref_name")
    return out

def _profile(entities,routes):
    mech=screen.mechanisms();a2g=_accession_gene_map()
    route_targets=defaultdict(list)
    for r in routes or []:
        for t in r.get("targets") or []:
            acc=t.get("accession")
            if acc:
                route_targets[acc].append({
                  "route_id":r.get("route_id"),"case_gene":r.get("case_gene"),
                  "target_concept":r.get("target_concept"),
                  "desired_action":r.get("desired_action"),
                  "rescue_objective_ids":r.get("rescue_objective_ids") or
                    ([r.get("rescue_objective_id")] if r.get("rescue_objective_id") else []),
                })
    profiles={}
    for e in entities or []:
        mrecs=[]
        for pid in e.get("chembl_parent_ids") or []:mrecs.extend(mech.get(pid,[]))
        targets={};useful=[]
        for m in mrecs:
            acc=m.get("accession");g=a2g.get(acc)
            key=(acc,m.get("action_sign"),m.get("action_type"))
            targets[key]={
              "accession":acc,"gene":g,"action_sign":m.get("action_sign"),
              "action_type":m.get("action_type"),"mechanism":m.get("mechanism"),
              "direct":m.get("direct"),"mechanism_id":m.get("mechanism_id")}
            for r in route_targets.get(acc,[]):
                desired=1 if r.get("desired_action")=="increase" else -1 if r.get("desired_action")=="decrease" else 0
                useful.append({**r,"target_gene":g,"accession":acc,
                               "drug_action_sign":m.get("action_sign"),
                               "drug_action_type":m.get("action_type"),
                               "compatible_direction":bool(desired and desired==m.get("action_sign")),
                               "mechanism_id":m.get("mechanism_id")})
        label=e.get("canonical_label")
        profiles[label]={
          "candidate":label,"entity_id":e.get("entity_id"),
          "approved_single":bool(e.get("global_approved_as_single") or
                                 e.get("solo_approved") or e.get("biologic_licensed") or
                                 not e.get("combo_only",False)),
          "known_mechanism_target_count":len(targets),
          "mechanism_targets":list(targets.values()),
          "route_target_matches":useful,
          "compatible_route_target_matches":sum(x["compatible_direction"] for x in useful),
          "polypharmacology_present":len(targets)>1}

    return profiles

def annotate_pairs(pairs,profiles):
    out=[]
    for p0 in pairs or []:
        p=dict(p0);prof=profiles.get(p.get("ingredient")) or {}
        matches=[x for x in prof.get("route_target_matches",[])
                 if x.get("route_id")==p.get("route_id")]
        p["polypharmacology"]={
          "known_mechanism_target_count":prof.get("known_mechanism_target_count",0),
          "polypharmacology_present":prof.get("polypharmacology_present",False),
          "route_target_mechanism_matches":matches,
          "compatible_route_target_mechanism_count":sum(bool(x.get("compatible_direction")) for x in matches),
          "all_known_targets":prof.get("mechanism_targets",[])[:40],
        }
        out.append(p)
    return out

def build(entities,pairs,outdir):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    route_map={}
    for p in pairs or []:
        r=p.get("route") or {}
        if r.get("route_id"):route_map[r["route_id"]]=r
    profiles=_profile(entities,list(route_map.values()))
    annotated=annotate_pairs(pairs,profiles)
    summary={
      "schema":SCHEMA,
      "profiled_entity_count":len(profiles),
      "entities_with_known_mechanisms":sum(bool(x["known_mechanism_target_count"]) for x in profiles.values()),
      "entities_with_multiple_known_targets":sum(bool(x["polypharmacology_present"]) for x in profiles.values()),
      "entities_with_compatible_route_mechanism":sum(bool(x["compatible_route_target_matches"]) for x in profiles.values()),
      "source":"ChEMBL mechanism records","source_license":"CC BY-SA 3.0",
      "bulk_data_redistributed_in_source_release":False,
      "clinical_recommendation":False}
    (outdir/"POLYPHARMACOLOGY_SUMMARY_V31.json").write_text(
      json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    # Full profiles stay local research output, not source release.
    (outdir/"POLYPHARMACOLOGY_PROFILES_V31.json").write_text(
      json.dumps({"schema":SCHEMA,"profiles":profiles},indent=2,ensure_ascii=False)+"\n",
      encoding="utf-8")
    return annotated,summary
