"""Public-label safety review for V3 single and combination hypotheses.

Uses public FDA label data through the existing V2 adapter. It can flag
warnings/contraindications and explicit pair mentions, but absence of a label
interaction is NEVER treated as safety clearance.
"""
from __future__ import annotations
from pathlib import Path
import json,re
from v2810 import deep

SCHEMA="RDMVA_SAFETY_REVIEW_V30_1"

def _norm(x):
    return re.sub(r"[^a-z0-9]+","",str(x or "").lower())

def _interaction_text(result):
    if not result:return ""
    vals=[]
    for key in ("drug_interactions","contraindications","boxed_warning","warnings","warnings_and_cautions"):
        v=result.get(key,[])
        if isinstance(v,list):vals.extend(str(x) for x in v)
        elif isinstance(v,str):vals.append(v)
    return " ".join(vals)

def review_candidates(names,entities,case_private,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    emap={x.get("canonical_label"):x for x in entities}
    pref=deep.chembl_pref_names(entities)
    results={}
    raw_labels={}
    for name in sorted(set(names)):
        entity=emap.get(name)
        if not entity:
            results[name]={"state":"HOLD","flags":["entity_not_resolved_in_approved_universe"],
                           "clinical_clearance":False}
            continue
        gate=deep.safety_gate(name,entity,case_private,pref,out/"labels")
        lab=deep.label(name,entity,pref,out/"labels")
        raw_labels[name]=lab.get("result")
        results[name]={**gate,"clinical_clearance":False,
                       "absence_of_flags_is_not_safety_clearance":True}
    payload={"schema":SCHEMA,"candidates":results,
             "clinical_clearance":False,
             "patient_identifiers_sent_external":False}
    (out/"SINGLE_AGENT_PUBLIC_LABEL_REVIEW.json").write_text(
        json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload,raw_labels

def annotate_hypotheses(track2,safety_payload,raw_labels):
    safety=safety_payload.get("candidates") or {}
    singles=[]
    for x0 in track2.get("single_agent_hypotheses",[]):
        x=dict(x0);x["public_label_safety_review"]=safety.get(x["candidate"],{"state":"NOT_REVIEWED"})
        singles.append(x)
    submission=[]
    for x0 in track2.get("submission_level_single_agents",[]):
        x=dict(x0);x["public_label_safety_review"]=safety.get(x["candidate"],{"state":"NOT_REVIEWED"})
        submission.append(x)
    combos=[]
    for c0 in track2.get("combination_hypotheses",[]):
        c=dict(c0);names=c.get("components") or []
        flags=[];explicit=[]
        for name in names:
            s=safety.get(name) or {"state":"NOT_REVIEWED"}
            if s.get("state")!="PASS":
                flags.append(name+":"+str(s.get("state")))
        for a in names:
            txt=_interaction_text(raw_labels.get(a))
            ntxt=_norm(txt)
            for b in names:
                if a==b:continue
                nb=_norm(b)
                if len(nb)>=5 and nb in ntxt:
                    explicit.append({"label_of":a,"mentions":b,
                                     "interpretation":"explicit_name_appears_in_public_label_interaction_or_warning_sections"})
        if explicit:
            state="EXPLICIT_PAIR_MENTION_REQUIRES_EXPERT_REVIEW"
        elif flags:
            state="COMPONENT_LABEL_FLAGS_REQUIRE_REVIEW"
        else:
            state="NO_EXPLICIT_PAIR_MENTION_NOT_SAFETY_CLEARANCE"
        c["combination_safety_review"]={
          "state":state,
          "component_states":{n:(safety.get(n) or {}).get("state","NOT_REVIEWED") for n in names},
          "component_review_flags":flags,
          "explicit_pair_mentions":explicit,
          "no_interaction_found_does_not_mean_safe":True,
          "patient_specific_medication_reconciliation_required":True,
          "clinical_clearance":False,
        }
        combos.append(c)
    result=dict(track2)
    result["single_agent_hypotheses"]=singles
    result["submission_level_single_agents"]=submission
    result["combination_hypotheses"]=combos
    result["public_label_safety_review_applied"]=True
    result["clinical_safety_clearance"]=False
    return result
