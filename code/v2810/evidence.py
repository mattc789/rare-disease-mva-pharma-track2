"""Research evidence upgrades from actually quoted, externally sourced experiments.

Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch

AI may suggest candidates and cite paper passages; it never supplies
"independently reviewed by a clinician". Tier A/B means research-evidence
quality (with disclosed verification limits), not suitability for treatment.
"""
from __future__ import annotations
import re
from collections import Counter,defaultdict
from copy import deepcopy

SCHEMA="V285_EVIDENCE_1"
TIER_DEFINITIONS={
 "A":"Strong source-checked disease-model functional rescue evidence in at least two separate publications, with exposure; not independent clinical validation.",
 "B":"Direct, quantitatively measured modulation of a relevant primary target or one source-checked functional study; patient-specific rescue unproven.",
 "C":"Direct recorded pharmacological mechanism without verified disease functional rescue.",
 "D":"Indirect network, perturbational-expression or structural hypothesis, no direct disease rescue.",
 "U":"No qualifying discovery signal in the current public snapshot; NOT proof of inactivity."
}

def _name_tokens(form_names):
    names=set()
    for name in form_names:
        low=re.sub(r"[^a-z0-9 ]"," ",str(name).lower())
        low=" ".join(low.split())
        if len(low)>=4:names.add(low)
        if len(low.split())>1 and len(low.split()[0])>=6:names.add(low.split()[0])
    return names

def _gene_aliases(primary_gene):
    g=str(primary_gene).upper()
    return {g.lower()}|({"bubr1","bub1b"} if g=="BUB1B" else set())

def _source_qualified_quote(q,form_names,gene):
    source=str(q.get("source_id") or "").strip()
    quote=str(q.get("quote") or "").strip()
    if not re.fullmatch(r"PMC\d{4,15}",source,re.I):return False
    low=" ".join(re.sub(r"[^a-z0-9 ]"," ",quote.lower()).split())
    names=_name_tokens(form_names)
    if not any(n in low for n in names):return False
    if not any(g in low for g in _gene_aliases(gene)):return False
    functional=(r"\b(?:aneuploidy|chromosome segregation|chromosomal instability|mitotic checkpoint|"
                r"spindle checkpoint|micronuclei|kinetochore|protein function|functional activity)\b")
    improvement=r"\b(?:restore|restores|restored|rescue|rescues|rescued|improv(?:e|es|ed)|recover(?:s|ed)?)\b"
    if not re.search(functional,low) or not re.search(improvement,low):return False
    if re.search(r"\b(?:failed to|did not|does not|doesn't|not restored|no rescue|lack of)\b",low):
        return False
    # Exact quotes are already checked by v285.ai.deep_attack against source passages.
    return True

def enrich_pairs_from_literature(pairs,review_rows,primary_genes):
    """Research-only upgrade; one candidate->route review cannot upgrade others."""
    gene=next(iter(sorted(primary_genes)),"")
    bykey={(x.get("scientific_group_id"),x.get("source_route_id")):x for x in review_rows}
    out=[];events=[]
    for p in pairs:
        row=deepcopy(p)
        matching=bykey.get((p.get("scientific_group_id"),p.get("route_id")))
        review=(matching or {}).get("ai_evaluation") or {}
        forms=row.get("member_labels") or [row["ingredient"]]
        gene_for_route=str((row.get('route') or {}).get('case_gene') or gene)
        valid_quotes=[q for q in review.get("support_quotes",[])
                      if _source_qualified_quote(q,forms,gene_for_route)]
        pubs={q["source_id"] for q in valid_quotes}
        no_conflict=not review.get("counterevidence_present",False) and review.get("adjudication")=="SUPPORTED"
        exposure=review.get("exposure_reported") is True
        direct=review.get("direct_route_modulation") is True
        same=review.get("same_drug_functional_rescue") is True
        original=row["tier"]
        if original not in ("A","B"):
            if no_conflict and exposure and direct and same and len(pubs)>=2:
                row["tier"]="A"
            elif no_conflict and direct and same and len(pubs)>=1:
                row["tier"]="B"
        row["publication_strength"]="source_checked_provisional" if row["tier"] in ("A","B") else "hypothesis_only"
        row["functional_evidence"]={
           "source_quoted_candidate_gene_functional_passages":len(valid_quotes),
           "distinct_publications_with_qualified_quotes":len(pubs),
           "validated_source_ids":sorted(pubs),
           "independent_experimental_replication_not_proven":True,
           "no_clinical_benefit_inferred":True,
           "AI_quote_extraction_alone_is_not_expert_curation":True,
           "AI_evidence_claim_present":bool(review),
           "target_modulation_reported":direct,
           "exposure_reported":exposure,
           "counterevidence_reported":bool(review.get("counterevidence_present")),
        }
        row["promotion_eligible"]=False # Scientific A != treatment or clinical lead.
        if row["tier"]!=original:events.append({
          "compound":row["ingredient"],"route":row["route_id"],
          "old_tier":original,"research_evidence_tier":row["tier"],
          "grounding":"quoted independent PMC IDs and structured requirements; experimental replication unverified"})
        out.append(row)
    out.sort(key=lambda x:({"A":0,"B":1,"C":2,"D":3}[x["tier"]],
      -x["score"],x["ingredient"],x["route_id"]))
    return out,events

def evidence_explanation(pair):
    level=pair.get("tier","D")
    signals=pair.get("signals") or []
    network=[s for s in signals if s.get("provider")=="ChEMBL_MoA_SIGNOR"]
    transcript=[s for s in signals if s.get("provider")=="perturbational_transcriptomics"]
    measured=[s for s in signals if s.get("provider")=="ChEMBL_measured_activity"]
    why=[]
    if network:
        steps=min(int(s.get("depth",99)) for s in network)
        why.append("Drug mechanism connected to route via SIGNOR network at "+str(steps)+" edge(s); these edges are context-sensitive.")
    if transcript:
        why.append("Matching perturbation-related RNA expression signals; not independent protein/function rescue.")
    if measured:
        why.append("ChEMBL target assay records exist and require quantitative target, assay direction and exposure review.")
    if any(s.get("provider")=="rdkit_morgan_similarity" for s in signals):
        why.append("Structural similarity only; not proof of similar functional effects.")
    ev=pair.get("functional_evidence") or {}
    if ev.get("distinct_publications_with_qualified_quotes"):
        why.append("Public full-text quotations mention the compound, primary gene and functional phenotype; source-bound AI check, not independent validation.")
    limits=[]
    if level=="A":
        limits.extend(["Independent cellular replication and clinical validation not established",
          "Assay doses and formulation may not be achievable/safe in this patient"])
    elif level=="B":
        limits.extend(["Insufficient independent disease-model replication",
          "No demonstrated restoration of this patient's allele-specific protein function"])
    elif level=="C":
        limits.extend(["Pharmacological mechanism does not prove cellular disease rescue",
          "Specific downstream functional effect and exposure remain unknown"])
    else:
        limits.extend(["Only indirect or exploratory observations",
          "No direct drug-to-defect rescue experiment established"])
    if pair.get("scientific_gate",{}).get("issues"):
        limits+=["Route uncertainties: "+", ".join(pair["scientific_gate"]["issues"][:5])]
    limits.extend(["Clinical safety and suitability remain unverified",
      "Allelic phase and residual functional protein may remain unknown",
      "Regulatory approval is not approval for this rare disease"])
    return {"tier":level,"definition":TIER_DEFINITIONS.get(level),
      "reasons":why or ["No qualifying positive pharmacological evidence"],
      "limitations":list(dict.fromkeys(limits)),
      "scientific_strength":"research_only_not_clinical_validation",
      "case_specific_treatment_recommendation":False}
