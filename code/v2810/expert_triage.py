"""Mechanism direction, source independence, gene confidence and safety-aware expert triage.

V2.8.8 is a research workflow; categorization is NOT an efficacy prediction.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from copy import deepcopy
from pathlib import Path
import csv,json,re

SCHEMA="V288_EXPERT_RESEARCH_TRIAGE_1"
# Small ordinal numbers are investigation-utility buckets, not likelihoods.
BUCKETS={
 "T1": "Investigate verified direct target pharmacology or source-checked functional model",
 "T2": "Investigate directly annotated drug mechanism with a plausible gene-linked route",
 "T3": "Indirect evidence only; generate experiments, do NOT rank by presumed efficacy",
 "T4": "Link to the implicated defect cannot be justified from available data",
 "HOLD": "Conflicting direction, nontransferable context, or genetic causal premise contradicted",
}
BUCKET_ORDER={"T1":0,"T2":1,"T3":2,"T4":3,"HOLD":4}
NEGATIVE=r"\b(?:suppress(?:ion|ed|es)?|down.?regulat(?:e|ion|ed|es)|decreas(?:e|es|ed|ing)|reduc(?:e|es|ed|tion)|knockdown|deplet(?:ion|ed|e))\b"
POSITIVE=r"\b(?:increas(?:e|es|ed|ing)|restor(?:e|es|ed)|up.?regulat(?:e|ion|ed|es)|stabiliz(?:e|es|ed)|enhanc(?:e|ed|es))\b"
GENE_TOKEN=r"(?<![A-Z0-9]){}(?![A-Z0-9])"

def is_gene_suppression_mechanism(text,gene):
    """Recognize suppression *of the patient's gene*, not nearby antitumor effects.

    This is a conservative flag; natural-language mechanism direction still
    needs published causal experiments for validation.
    """
    if not gene or not text:return False
    plain=re.sub(r"\s+"," ",str(text)).upper()
    token=GENE_TOKEN.format(re.escape(gene.upper()))
    negative=r"(?:DOWN.?REGULAT(?:E|ED|ES|ION)|SUPPRESS(?:ED|ES|ION)?|DECREAS(?:E|ED|ES)|REDUC(?:E|ED|ES|TION)|DEPLET(?:E|ED|ION)|INHIBIT(?:S|ED|ION)|SILENC(?:E|ED|ING))"
    # e.g. "ERCC2 downregulation", "ERCC2 expression was reduced".
    after=rf"{token}(?:\s+(?:GENE|PROTEIN|MRNA|TRANSCRIPT|EXPRESSION|ABUNDANCE|LEVELS?|WAS|IS|ARE|WERE|HAS|BEEN|DIRECTLY|SIGNIFICANTLY)){{0,5}}\s+{negative}"
    # e.g. "downregulates ERCC2", "decreased ERCC2 expression".
    before=rf"\b{negative}(?:\s+(?:THE|EXPRESSION|OF|PROTEIN|MRNA|TRANSCRIPT|LEVELS?|DIRECTLY|SIGNIFICANTLY)){{0,4}}\s+{token}"
    for pat in (before,after):
        for m in re.finditer(pat,plain,re.I):
            prefix=plain[max(0,m.start()-17):m.start()]
            if re.search(r"\b(?:NOT|NO|WITHOUT)\s*$",prefix):
                continue
            return True
    return False

def _direct_assay(signals):
    return [s for s in signals if s.get("provider")=="ChEMBL_measured_activity"
        and s.get("direct_target_assay_verified") is True and
        s.get("value_nM") is not None]

def _moa(signals):
    return [s for s in signals if s.get("provider")=="ChEMBL_MoA_SIGNOR"
            and s.get("depth")==0]

def _indirect(signals):
    return [s for s in signals if s.get("provider") in (
      "perturbational_transcriptomics","rdkit_morgan_similarity")]

def assessment(pair,genetic_model,patient_model=None,safety=None):
    """Route assessment with explicit unresolvable causality, not artificial patient efficacy."""
    route=pair.get("route") or {}
    scientific=route.get("scientific_gate") or {}
    predictor=pair.get("patient_first_prediction") or {}
    gene=str(predictor.get("primary_gene_hypothesis") or "").upper()
    allgenes=(genetic_model or {}).get("genes") or {}
    valid_gene=gene in allgenes
    gmodel=allgenes.get(gene) or {}
    signals=pair.get("signals") or []
    direct=_direct_assay(signals);moa=_moa(signals);indirect=_indirect(signals)
    network=[x for x in signals if x.get("provider")=="ChEMBL_MoA_SIGNOR"]
    hops=min((x.get("depth",99) for x in network),default=None)
    state=str(route.get("state") or "")
    link=str(predictor.get("gene_association_method") or pair.get("patient_first_hypothesis_status") or "")
    # v286 safe summaries store prediction_class, not full gene link.
    if not link and route.get("family")=="direct_gene" and gene==str(route.get("target_concept") or "").upper():
        link="direct_gene_target"
    issues=list(scientific.get("issues") or [])
    problems=[];needs=[]
    sources={x.get("provider") for x in signals if x.get("provider")}
    if not valid_gene:problems.append("drug_route_not_attributable_to_a_selected_patient_gene")
    if gmodel.get("gene_mechanism_confidence")=="ALLELES_CONTRADICTED_BY_CONCORDANT_BENIGN_EVIDENCE":
        problems.append("gene_variants_currently_supported_as_benign")
    if scientific.get("hard_hold") or state=="hold":
        problems.append("biological_route_explicitly_held_by_source_adjudication")
    if route.get("family")=="literature_discovered" and valid_gene:
        text=" ".join((str(route.get("target_concept") or ""),
                      str(route.get("mechanism_concept") or "")))
        if is_gene_suppression_mechanism(text,gene):
            problems.append("source_describes_suppression_of_the_gene_we_aim_to_restore")
    risk_targets=sorted({str(x.get("gene") or "").upper() for x in route.get("targets",[])}
           & {"BUB1B","CDC20","CDK1","CDK2","PLK1","AURKB","CHEK1"})
    if risk_targets:
        needs.append("cell_cycle_and_cancer_genome_stability_risk_critical")
    if not gmodel.get("protein_rescue_endpoint_measured",False):
        needs.append("patient_relevant_functional_assay_not_performed")
    if gmodel.get("phase_state") not in ("homozygous_observed","trans_supported_by_vcf"):
        needs.append("trans_phase_or_gene_inheritance_evidence_not_established")
    if not direct:needs.append("no_direct_directional_functional_target_assay")
    if not (pair.get("functional_evidence") or {}).get("distinct_publications_with_qualified_quotes",0):
        needs.append("no_verified_same_compound_disease_rescue_study")
    if indirect and not direct and not moa:
        needs.append("RNA_perturbation_or_similarity_does_not_establish_direction_of_protein_effect")
    if not (safety or {}).get("state") in ("PASS","HOLD","CLOSED"):
        needs.append("public_product_label_review_missing")
    if (safety or {}).get("state") in ("HOLD","CLOSED"):
        needs.append("known_product_label_warning_requires_expert_risk_review")
    if (pair.get("patient_first_prediction") or {}).get("prediction_class")=="NOT_ELIGIBLE_REGISTRY":
        problems.append("single_ingredient_regulatory_eligibility_absent")
    if problems:
        bucket="HOLD" if any(z.startswith(("gene_variants_","biological_route_",
                "source_describes_","single_ingredient_")) for z in problems) else "T4"
    elif direct or ((pair.get("functional_evidence") or {}).get("distinct_publications_with_qualified_quotes",0)>0
                    and pair.get("tier") in ("A","B")):
        bucket="T1"
    elif moa and route.get("family")=="direct_gene":
        bucket="T2"
    elif (network and hops is not None and hops<=1
          and str(route.get("family") or "") in ("direct_gene","signed_upstream","proteostasis")):
        bucket="T3"  # one-hop remains conditional, not a direct drug-target proof
    elif signals:
        bucket="T3"
    else:bucket="T4"
    direction="not_validated"
    if direct:direction="verified_for_direct_target_assay_not_patient_rescue"
    elif moa:direction="ChEMBL_MoA_context_not_disease_rescue"
    elif indirect:direction="unknown_from_expression_or_similarity"
    elif network:direction="signed_network_path_context_specific"
    reasons=[]
    if direct:reasons.append("verified_directional_functional_target_assay="+str(len(direct)))
    if moa:reasons.append("direct_drug_mechanism_annotated")
    if network:reasons.append("drug_to_route_network_shortest_hops="+str(hops))
    if indirect:reasons.append("indirect_expression_or_similarity_signals="+str(len(indirect)))
    if not reasons:reasons=["no_qualifying_drug_mechanistic_signal"]
    why_not=[]
    if pair.get("tier") not in ("A","B"):why_not.append("no_credible_high_tier_disease_specific_rescue")
    why_not += needs + problems
    if route.get("family")=="direct_gene" and gene:
        reasons.append("route_targets_patient_hypothesis_gene="+gene)
    # Literature evidence does not certify independent replication, even with distinct PMCs.
    functional=(pair.get("functional_evidence") or {})
    sources_independence="NOT_VERIFIED" if functional.get("independent_experimental_replication_not_proven",True) else "MANUALLY_REVIEWED"
    return {
      "schema":SCHEMA,"candidate":pair.get("ingredient"),"scientific_group_id":pair.get("scientific_group_id"),
      "route_id":pair.get("route_id"),"gene_hypothesis":gene or None,
      "tier_evidence_A_B_C_D":pair.get("tier","D"),
      "investigation_bucket":bucket,"bucket_meaning":BUCKETS[bucket],
      "functionally_predicted_patient_benefit":None,
      "calibrated_efficacy_probability":None,
      "mechanistic_compatibility_not_computable":bucket not in ("T1","T2"),
      "drug_target_action_direction_evidence":direction,
      "gene_causality_state":gmodel.get("gene_mechanism_confidence","NOT_ASSESSED"),
      "genetic_phase_state":gmodel.get("phase_state","unknown"),
      "clinical_safety_clearance":False,
      "label_state":(safety or {}).get("state","NOT_ASSESSED"),
      "source_study_independence":sources_independence,
      "direct_functional_target_assay_count":len(direct),
      "direct_moa_count":len(moa),
      "shortest_network_hops":hops,
      "signal_provider_types":sorted(sources),
      "positive_experimental_disease_rescue_not_established":not bool(functional.get("distinct_publications_with_qualified_quotes")),
      "rationale":reasons,
      "hypothesis_blockers":list(dict.fromkeys(issues+problems+needs)),
      "falsification_plan":design_experiments(pair,gmodel),
      "treatment_recommendation":False,
    }

def design_experiments(pair,gene):
    route=pair.get("route") or {}
    cls={x.get("consequence") for x in gene.get("variants",[])}
    target=str(route.get("target_concept") or "proposed target")
    tests=[
       "Reconcile ClinVar condition-specific assertions and variant provenance; establish mode of inheritance.",
       "Measure allele-specific RNA/protein in an appropriate disease model, with pathogenic/benign controls.",
       "Confirm direct drug action on "+target+" and assay concentration-response under appropriate conditions.",
       "Assess drug concentration in relevant tissue or cells and compare with human unbound exposure where available.",
       "Replicate disease-specific functional rescue independently before interpreting as efficacy.",
       "Assess genotype-dependent toxicity, off-target actions and safety in the intended age group.",
    ]
    if "stop_gained" in cls or "frameshift_variant" in cls:
        tests.insert(2,"Distinguish transcript decay, truncated protein and retained functional activity.")
    if any(x and "splice" in x for x in cls):
        tests.insert(2,"Use RNA sequencing or a validated splice assay to quantify productive isoforms.")
    if "missense_variant" in cls:
        tests.insert(2,"Quantify mutant protein folding, localization and interactions; test residual function.")
    return tests

def order_pairs(pairs,genetic_model,safety=None):
    """The old 0–100 heuristic is retained for audit but NOT used to order weak discoveries."""
    merged=[]
    safety=safety or {}
    counts=Counter()
    for p in pairs:
        row=deepcopy(p)
        assessment_row=assessment(row,genetic_model,safety=safety.get(p.get("ingredient")))
        row["expert_triage"]=assessment_row
        row["expert_investigation_bucket"]=assessment_row["investigation_bucket"]
        counts[assessment_row["investigation_bucket"]]+=1
        merged.append(row)
    # Within T3, order by genetic model relevance and network path,
    # not by an invented patient efficacy estimate; identical evidence is tied.
    merged.sort(key=lambda r:(
       BUCKET_ORDER[r["expert_investigation_bucket"]],
       {"A":0,"B":1,"C":2,"D":3}.get(r.get("tier"),4),
       0 if r["expert_triage"].get("gene_causality_state","").startswith("SUPPORTED_VARIANT") else 1,
       r["expert_triage"].get("shortest_network_hops") if r["expert_triage"].get("shortest_network_hops") is not None else 99,
       str(r.get("route_id") or ""), str(r.get("ingredient") or "")))
    return merged,{"schema":SCHEMA,"route_drug_pairs":len(merged),
      "bucket_counts":dict(counts),
      "no_calibrated_clinical_efficacy":True,
      "weak_scores_are_not_used_to_assign_clinical_priorities":True,
      "clinical_suitability_assessed":False}

def report(outdir,pairs,genetics,safety=None,top=30):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    safety=safety or {}
    rows,_=order_pairs(pairs,genetics,safety)
    best={}
    alternatives=defaultdict(list)
    for p in rows:
        x=p["expert_triage"];key=x.get("scientific_group_id") or x.get("candidate")
        alternatives[key].append(x)
        if key not in best:best[key]=deepcopy(x)
    catalogue=list(best.values())
    catalogue.sort(key=lambda x:(
       BUCKET_ORDER[x["investigation_bucket"]],
       {"A":0,"B":1,"C":2,"D":3}.get(x["tier_evidence_A_B_C_D"],4),
       x["gene_hypothesis"] or "",x["route_id"] or "",x["candidate"] or ""))
    for r in catalogue:
        k=r.get("scientific_group_id") or r.get("candidate")
        r["mechanism_hypotheses_for_same_molecule"]=len(alternatives[k])
        r["conflicting_route_directions_not_resolved"]=False
    with (out/"EXPERT_ALL_ROUTE_EVIDENCE.jsonl").open("w",encoding="utf-8") as f:
        for p in rows:f.write(json.dumps(p["expert_triage"],ensure_ascii=False,sort_keys=True)+"\n")
    with (out/"EXPERT_DRUG_CATALOGUE.jsonl").open("w",encoding="utf-8") as f:
        for item in catalogue:f.write(json.dumps(item,ensure_ascii=False,sort_keys=True)+"\n")
    counts=Counter(x["investigation_bucket"] for x in catalogue)
    summary={
       "schema":SCHEMA,"distinct_approved_drug_groups_with_signals":len(catalogue),
       "evaluated_drug_route_pairs":len(rows),"expert_bucket_counts":dict(counts),
       "cohort_variant_diagnosis_confirmed":False,
       "experimental_disease_rescue_confirmed":False,
       "clinical_treatment_recommendations":0,
       "report":"EXPERT_RESEARCH_REPORT.md",
       "gene_hypotheses":list(genetics.get("genes",{})),
       "annotation_is_not_qualified_clinical_diagnosis":True,
    }
    (out/"EXPERT_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    with (out/"EXPERT_RESEARCH_RANKING.csv").open("w",encoding="utf-8",newline="") as f:
        fields=["candidate","gene","bucket","research_tier","route","direct_assays","mechanism_hops",
                "safety_label","independent_study_proof","rationale","blockers"]
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in catalogue:
            writer.writerow({
              "candidate":r["candidate"],"gene":r["gene_hypothesis"],
              "bucket":r["investigation_bucket"],"research_tier":r["tier_evidence_A_B_C_D"],
              "route":r["route_id"],"direct_assays":r["direct_functional_target_assay_count"],
              "mechanism_hops":r["shortest_network_hops"],"safety_label":r["label_state"],
              "independent_study_proof":r["source_study_independence"],
              "rationale":"; ".join(r["rationale"]), "blockers":"; ".join(r["hypothesis_blockers"]),
            })
    lines=["# Rare Disease MVA Pharma — Expert computational research report (V2.8.8)",
      "Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch",
      "",
      "**No patient-specific medicine has been shown safe or effective by these computations.**",
      "Research-evidence categories A/B/C/D remain separate from the mechanistic investigation buckets T1–T4/HOLD.",
      "The previous 0–100 heuristic is NOT a validated efficacy estimate and is not used to calculate clinical benefit.",
      "",
      f"Drug groups: {len(catalogue)}; drug–mechanism routes: {len(rows)}.",
      "## What can and cannot be concluded"]
    for b,detail in BUCKETS.items():
        lines.append(f"- **{b}** ({counts[b]}): {detail}.")
    lines += ["","## Ranked by evidence utility, not presumed therapeutic benefit"]
    for idx,r in enumerate(catalogue[:top],1):
        lines += ["",f"### {idx}. {r['candidate']} — {r['investigation_bucket']} — evidence {r['tier_evidence_A_B_C_D']}",
          f"**Disease-gene hypothesis:** {r['gene_hypothesis'] or 'unresolved'}; genotype causality: {r['gene_causality_state']}",
          f"**Mechanism:** {r['route_id']}; drug target direction: {r['drug_target_action_direction_evidence']}",
          f"**Gene–drug evidence:** {', '.join(r['rationale'])}",
          f"**Safety source status:** {r['label_state']} (not a patient clearance)",
          "**Critical limitations:**"]
        lines += ["- "+x for x in r["hypothesis_blockers"][:7]]
        lines += ["**Next disprovable research tests:**"]
        lines += ["- "+x for x in r["falsification_plan"][:4]]
    lines += ["","## Scope limits",
      "No gnomAD population annotation is inferred from generic VCF FORMAT/INFO AF.",
      "Computational estimates do not replace genetic expert diagnosis, family segregation, functional assays, tissue exposure or toxicology.",
      "Cellular RNA expression, chemical similarity or a network path does not imply rescue of an allele-specific phenotype.",
      "A source-checked paper is not experimental replication or guaranteed therapeutic benefit."]
    (out/"EXPERT_RESEARCH_REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return summary,catalogue[:top]


def safety_review_queue(pairs,already_reviewed,count=20):
    """Broaden official-label lookups over distinct molecules, genes and mechanisms.

    Cost-bounded. A missing label still means UNKNOWN, never 'safe'.
    """
    if count<=0:return []
    previous=set(already_reviewed or ())
    accepted=[];seen=set(previous);touched=set()
    candidates=sorted(pairs,key=lambda p:(
       BUCKET_ORDER.get(p.get("expert_investigation_bucket"),5),
       {"A":0,"B":1,"C":2,"D":3}.get(p.get("tier"),4),
       p.get("route",{}).get("case_gene") or "",
       p.get("route_id") or "",p.get("ingredient") or ""))
    for p in candidates:
        x=p.get("expert_triage") or {}
        if x.get("investigation_bucket") in ("HOLD","T4"):continue
        k=p.get("scientific_group_id") or p.get("ingredient")
        drug=p.get("ingredient")
        bucket=(p.get("route",{}).get("case_gene"),p.get("route_id"))
        if not drug or drug in previous or k in seen or bucket in touched:continue
        accepted.append(drug);seen.add(k);touched.add(bucket)
        if len(accepted)>=count:return accepted
    for p in candidates:
        x=p.get("expert_triage") or {}
        if x.get("investigation_bucket") in ("HOLD","T4"):continue
        k=p.get("scientific_group_id") or p.get("ingredient")
        drug=p.get("ingredient")
        if not drug or drug in previous or k in seen:continue
        accepted.append(drug);seen.add(k)
        if len(accepted)>=count:return accepted
    return accepted
