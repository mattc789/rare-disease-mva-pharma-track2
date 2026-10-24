"""Human- and machine-readable A/B/C/D/U evidence and limitations reports.

Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""
from __future__ import annotations
import csv,json
from pathlib import Path
from collections import Counter
from . import evidence,registry_output,case_science
from .core import save,jsonl

HEADING="Rare Disease MVA Pharma — V2.8.6"
BYLINE="Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch"

def baseline_safety(private,public_label=None):
    x=case_science.case_suitability_constraints(private)
    if isinstance(public_label,dict):
        x["public_label_flags"]=public_label.get("flags",[])
        x["public_label_review_state"]=public_label.get("state","UNKNOWN")
        x["public_label_query"]=public_label.get("label_query")
    else:
        x["public_label_flags"]=[]
        x["public_label_review_state"]="not_attempted_or_not_available"
    x["can_recommend_for_patient"]=False
    x["clinical_safety_clearance"]=False
    return x

def report_entries(grouped,pairs,private,registry_index,safety_by_name=None):
    safety_by_name=safety_by_name or {}
    best={}
    alternative=Counter()
    for p in pairs:
        key=p.get("scientific_group_id") or p["ingredient"]
        if key not in best:best[key]=p
        alternative[key]+=1
    rows=[]
    for entry in grouped:
        gid=entry.get("scientific_group_id") or entry["canonical_label"]
        pair=best.get(gid)
        name=entry["canonical_label"]
        src=registry_output.for_candidate(name,registry_index)
        tier=pair["tier"] if pair else "U"
        details=evidence.evidence_explanation(pair) if pair else {
          "tier":"U","definition":evidence.TIER_DEFINITIONS["U"],
          "reasons":["No qualified mechanism/assay signal in the indexed public sources; not proof of therapeutic inactivity"],
          "limitations":["Incomplete experimental coverage or unlinked compound identity may explain absence",
                         "Patient-specific suitability is unverified","Missing data are not a negative result"],
          "scientific_strength":"not_assessed_or_no_qualifying_evidence",
          "case_specific_treatment_recommendation":False}
        row={"candidate":name,"scientific_group_id":gid,
             "related_product_forms":entry.get("member_labels",[name]),
             "evidence_tier":tier,"research_score":pair.get("score",0) if pair else 0,
             "evidence_definition":evidence.TIER_DEFINITIONS[tier],
             "mechanism_route":pair["route_id"] if pair else None,
             "route_target":pair["route"].get("target_concept") if pair else None,
             "mechanism_action":pair["route"].get("desired_action") if pair else None,
             "route_family":pair["route"].get("family") if pair else None,
             "number_of_explored_routes":alternative.get(gid,0),
             "supporting_signal_providers":pair.get("providers",[]) if pair else [],
             "reasons_for_tier":details["reasons"],
             "limitations":details["limitations"],
             "reference_studies":(pair.get("functional_evidence") or {}).get("validated_source_ids",[]) if pair else [],
             "regulatory":src,
             "safety_and_suitability":baseline_safety(private,safety_by_name.get(name)),
             "clinical_recommendation":False,
             "implies_approval_for_mva":False}
        rows.append(row)
    rows.sort(key=lambda r:({"A":0,"B":1,"C":2,"D":3,"U":4}[r["evidence_tier"]],
                  -r["research_score"],r["candidate"]))
    return rows

def create(outdir,grouped,pairs,private,case_model_summary,registry_index,
           shortlist=None,ai_review_queue=None,safety_by_name=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    rows=report_entries(grouped,pairs,private,registry_index,safety_by_name)
    summary={"version":"2.8.5","report_schema":"RESEARCH_EVIDENCE_REPORT_1",
        "title":HEADING,"developer":BYLINE,
        "registered_label_count":sum(len(x.get("member_labels",[])) for x in grouped),
        "distinct_research_groups":len(grouped),
        "route_drug_pairs":len(pairs),
        "evidence_tier_counts":dict(Counter(row["evidence_tier"] for row in rows)),
        "evidence_definitions":evidence.TIER_DEFINITIONS,
        "genetic_summary":case_model_summary,
        "clinical_recommendation":False,"candidate_suitability_assessed":False,
        "source_snapshot_date":"2026-10-10",
        "limitation":"Research evidence categories do not certify clinical effectiveness. Scores depend on public evidence coverage; U is not inactivity.",
        "literature_reviewed":len(ai_review_queue or []),
        "public_label_safety_lookups":len(safety_by_name or {}),
    }
    save(out/"CANDIDATE_REPORT_SUMMARY.json",summary)
    jsonl(out/"CANDIDATE_REPORT_ALL.jsonl",rows)
    fields=["candidate","scientific_group_id","evidence_tier","research_score","route_target",
            "mechanism_route","number_of_explored_routes","registries",
            "registered_forms","reasons","limitations","safety_state"]
    with (out/"CANDIDATES_ALL.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for x in rows:
            writer.writerow({
                "candidate":x["candidate"],
                "scientific_group_id":x["scientific_group_id"],
                "evidence_tier":x["evidence_tier"],
                "research_score":x["research_score"],
                "route_target":x["route_target"],
                "mechanism_route":x["mechanism_route"],
                "number_of_explored_routes":x["number_of_explored_routes"],
                "registries":"; ".join(x["regulatory"].get("registry_names",[])),
                "registered_forms":"; ".join(x["related_product_forms"]),
                "reasons":"; ".join(x["reasons_for_tier"]),
                "limitations":"; ".join(x["limitations"]),
                "safety_state":x["safety_and_suitability"]["public_label_review_state"],
            })
    valid=[x for x in rows if x["evidence_tier"]!="U"]
    candidates=valid[:min(30,len(valid))]
    # Distinct mechanism diversity is preserved in full JSONL, top ranking is
    # evidence-based with no clinical-suitability claim.
    save(out/"CANDIDATES_TOP.json",candidates)
    lines=["# "+HEADING,BYLINE,"",
      "## Research results",
      "**These are research hypotheses, not prescriptions or clinical treatment recommendations.**",
      "",
      f"Regulatory-labelled research groups reviewed: {len(rows)}. Evidence-eligible groups: {len(valid)}.",
      "",
      "## Evidence tiers"]
    for tier in ("A","B","C","D","U"):
        lines.append(f"- **{tier}** — {evidence.TIER_DEFINITIONS[tier]} "
                     f"Count: {summary['evidence_tier_counts'].get(tier,0)}.")
    lines += ["",
      "## What was observed from the local VCF (public-safe summary)",
      f"- Primary gene hypotheses: {len(case_model_summary.get('primary_gene_hypotheses',[]))}",
      f"- Clinically reported phenotype features: {case_model_summary.get('phenotype_counts',{}).get('present',0)}",
      "- The software does **not** establish allele-specific residual protein function, patient drug exposure or clinical safety.",
      "",
      "## Ranked research candidates — why this tier, and what is missing"]
    for index,row in enumerate(candidates,1):
        lines+=["",f"### {index}. {row['candidate']} — tier {row['evidence_tier']}",
          f"**Target/route:** {row['route_target'] or 'not established'}",
          f"**Related forms:** {', '.join(row['related_product_forms']) or 'none'}",
          "**Why assigned this grade:**"]
        lines += ["- "+reason for reason in row["reasons_for_tier"]]
        lines += ["**Limitations:**"]
        lines += ["- "+limitation for limitation in row["limitations"]]
        lines += ["**Official product records:**"]
        lines += registry_output.text_lines(row["regulatory"])
        lines+=["**Patient suitability:** not established; clinical review required.",
               f"**Safety information:** {row['safety_and_suitability']['public_label_review_state']}."]
    lines += ["","## Key cautions",
      "- Gene-expression perturbations and signed network edges do not demonstrate rescue of the genetic defect.",
      "- Authorisation of a medicine for another indication is not an authorisation for MVA.",
      "- Drug exposure, paediatric safety, organ function, other medicines, allele phase and functional assays must be independently assessed.",
      "- Missing evidence or an unlinked ChEMBL identity must not be interpreted as proof a drug is ineffective.",
      "- Any provisional A/B research grade requires expert confirmation and does not imply a treatment is safe or appropriate."]
    (out/"CANDIDATE_REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return summary,candidates
