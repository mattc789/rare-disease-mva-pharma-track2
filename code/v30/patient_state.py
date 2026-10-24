"""Patient State Engine V3.0.

Builds a local mechanistic state from VCF-derived evidence without turning
uncertainty into diagnosis. The local state may contain hashed allele refs;
the API-safe envelope produced elsewhere never contains them.
"""
from __future__ import annotations
import hashlib,json
from collections import Counter

SCHEMA="RDMVA_PATIENT_STATE_V30_1"

def _sid(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,ensure_ascii=True).encode()).hexdigest()[:16]

def _by_gene(rows):
    return {str(x.get("gene") or "").upper():x for x in rows or [] if x.get("gene")}

def build(private,patient_mechanisms,genomic_expert=None,protein_context=None,
          specialist_summary=None):
    genomic_expert=genomic_expert or {}
    protein_context=protein_context or {}
    specialist_summary=specialist_summary or {}
    mech=_by_gene(patient_mechanisms.get("primary_gene_hypotheses"))
    genes=[]
    for gene,variants in sorted((private.get("genes") or {}).items()):
        m=mech.get(gene,{})
        effects=list(m.get("variant_effect_models") or [])
        classes=sorted({str(x.get("effect_class") or "uncertain_effect") for x in effects})
        assertions=Counter(str(x.get("clinical_assertion") or "not_available") for x in effects)
        phase=m.get("phase") or {"state":"unresolved"}
        gx=(genomic_expert.get("genes") or {}).get(gene,{})
        px=(protein_context.get("genes") or {}).get(gene,{})
        evidence_flags={
          "variant_class_modelled":bool(effects),
          "phase_evidence_available":phase.get("state") not in (None,"unresolved","not_assessable"),
          "clinvar_reviewed_pathogenic_present":assertions.get("reported_pathogenic",0)>0,
          "clinvar_conflict_present":assertions.get("conflicting",0)+assertions.get("conflicting_or_mixed",0)>0,
          "protein_context_available":bool(px),
          "population_frequency_source_available":bool(
              genomic_expert.get("population_frequency_source",{}).get("state") not in
              (None,"NOT_PROVIDED","not_provided","unavailable")),
          "specialist_genomics_integrated":bool(specialist_summary),
        }
        completeness=sum(bool(v) for v in evidence_flags.values())
        tier="E3" if completeness>=6 else "E2" if completeness>=4 else "E1" if completeness>=2 else "E0"
        gene_state={
          "gene":gene,
          "variant_count":len(variants),
          "effect_classes":classes,
          "phase_state":phase.get("state","unresolved"),
          "phase_trans_established":bool(phase.get("trans_alleles_established",False)),
          "clinical_assertion_counts":dict(assertions),
          "residual_product_state":private.get("residual_protein_status","unknown"),
          "functional_measurement_available":bool(private.get("expression_measured",False)),
          "evidence_flags":evidence_flags,
          "evidence_completeness_tier":tier,
          "evidence_completeness_is_probability":False,
          "mechanism_confirmed":bool(m.get("actual_disease_mechanism_confirmed",False)),
          "variant_refs":[x.get("variant_ref") for x in effects if x.get("variant_ref")],
          "protein_context_source":px.get("source") if isinstance(px,dict) else None,
          "expert_genomic_state":gx.get("state") or gx.get("gene_state") or "available" if gx else "not_available",
        }
        gene_state["state_id"]="GENE_STATE:"+_sid({
          "gene":gene,"classes":classes,"phase":gene_state["phase_state"],
          "assertions":dict(assertions)})
        genes.append(gene_state)
    present=[x for x in private.get("phenotype_features",[])
             if x.get("subject")=="proband" and x.get("presence")=="reported_present"]
    absent=[x for x in private.get("phenotype_features",[])
            if x.get("subject")=="proband" and x.get("presence")=="reported_absent"]
    state={
      "schema":SCHEMA,
      "case_id_hash":hashlib.sha256(str(private.get("case_id","")).encode()).hexdigest()[:16],
      "genes":genes,
      "phenotype_summary":{"reported_present":len(present),"reported_absent":len(absent),
                           "other_or_unknown":max(0,len(private.get("phenotype_features",[]))-len(present)-len(absent))},
      "diagnostic_certainty":private.get("diagnostic_certainty","unknown"),
      "residual_protein_status":private.get("residual_protein_status","unknown"),
      "compound_heterozygosity_confirmed":bool(private.get("compound_heterozygosity_confirmed",False)),
      "raw_variant_strings_exported":False,
      "clinical_diagnosis_asserted":False,
      "clinical_recommendation":False,
    }
    state["state_id"]="PATIENT_STATE:"+_sid({
      "case":state["case_id_hash"],
      "genes":[{"g":x["gene"],"c":x["effect_classes"],"p":x["phase_state"]} for x in genes]})
    return state
