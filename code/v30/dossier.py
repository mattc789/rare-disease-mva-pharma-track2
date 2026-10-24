"""Final research dossier writer V3.0."""
from __future__ import annotations
from pathlib import Path
import json

SCHEMA="RDMVA_RESEARCH_DOSSIER_V30_1"

def build(patient_state,mechanisms,objectives,falsification,ai_review=None):
    candidates=falsification.get("candidates",[])
    leads=[x for x in candidates if x.get("decision")=="CONDITIONAL_RESEARCH_LEAD"]
    holds=[x for x in candidates if x.get("decision")=="HOLD"]
    primary=leads[0] if leads else None
    unresolved=holds[0] if holds else None
    return {
      "schema":SCHEMA,
      "patient_state_id":patient_state.get("state_id"),
      "research_conclusion":"CONDITIONAL_RESEARCH_LEAD" if primary else "NO_VALIDATED_LEAD",
      "proposed_research_lead":primary,
      "highest_priority_unresolved_hypothesis":unresolved,
      "candidate_evaluations":candidates,
      "mechanism_hypothesis_count":len(mechanisms.get("hypotheses",[])),
      "rescue_objective_count":len(objectives.get("objectives",[])),
      "ai_critic":ai_review,
      "ai_can_upgrade_evidence":False,
      "clinical_safety_clearance":False,
      "clinical_efficacy_claim":False,
      "dose_recommendation":False,
      "clinical_recommendation":False,
      "interpretation":"A conditional research lead is an experimentally testable hypothesis, not a treatment recommendation."
    }

def write(outdir,dossier):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/"V30_RESEARCH_DOSSIER.json").write_text(
      json.dumps(dossier,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lead=dossier.get("proposed_research_lead")
    unresolved=dossier.get("highest_priority_unresolved_hypothesis")
    lines=[
      "# Rare Disease MVA Pharma — V3.0 Research Dossier","",
      "**Research-only output. This document is not a prescription, treatment recommendation, safety clearance, or efficacy estimate.**","",
      "## Conclusion",dossier["research_conclusion"],"",
    ]
    if lead:
        lines += ["## Conditional research lead",f"**Candidate:** {lead.get('candidate')}",
          f"**Gene hypothesis:** {lead.get('gene')}",
          f"**Rescue objective:** {lead.get('objective_type')}",
          f"**Route target:** {lead.get('route_target')}",
          "","### Why it remains conditional"]
        lines += ["- "+x.replace("_"," ") for x in lead.get("hold_reasons",[])[:12]]
        lines += ["","### Highest-information falsification tests"]
        lines += ["- "+x for x in lead.get("falsification_tests",[])[:8]]
    elif unresolved:
        lines += ["## No conditional lead passed the current evidence gates",
          f"Highest-priority unresolved research hypothesis: **{unresolved.get('candidate')}**.",
          "It remains on HOLD and must not be interpreted as a proposed treatment.","",
          "### Main unresolved blockers"]
        lines += ["- "+x.replace("_"," ") for x in unresolved.get("hold_reasons",[])[:12]]
        lines += ["","### Experiments that could resolve it"]
        lines += ["- "+x for x in unresolved.get("falsification_tests",[])[:8]]
    else:
        lines += ["No candidate currently maps to a sufficiently supported patient-specific rescue path."]
    lines += ["","## Methodological guardrails",
      "- Raw VCF alleles, coordinates, genotype strings and full phenotype profiles are not part of the AI critique envelope.",
      "- Missing evidence creates HOLD, not benefit.",
      "- Contradictory direction or hard scientific/safety gates create KILL.",
      "- AI critique cannot change deterministic evidence tiers.",
      "- Every candidate must remain experimentally falsifiable."
    ]
    (out/"V30_RESEARCH_DOSSIER.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return out/"V30_RESEARCH_DOSSIER.json"
