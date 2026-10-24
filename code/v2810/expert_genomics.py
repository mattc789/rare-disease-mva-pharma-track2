"""Expert genetic-evidence dossier from the locally interpreted VCF.

Rare Disease MVA Pharma V2.8.8 | ExaMedix Technologies.
No ACMG/AMP classification is assigned by this lightweight pipeline.
No allele/genome strings are sent to an external AI.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import hashlib,json,math

SCHEMA="V288_EXPERT_GENOMICS_1"
# Quality screening markers; thresholds are heuristics, not diagnostic rules.
HET={"0/1","1/0","0|1","1|0"}
HOM={"1/1","1|1"}
LOSS_CLASSES={"stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant"}

def safe_hash(v):
    raw="|".join(str(v.get(k) or "") for k in ("gene","chrom","pos","ref","alt","hgvsp"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]

def _value(v):
    try:
        out=float(v)
        return out if math.isfinite(out) else None
    except (ValueError,TypeError,OverflowError):return None

def load_local_population(path):
    """Optional sourced, locally curated GRCh38 AF table; NEVER use caller VCF INFO/AF as gnomAD."""
    if not path or not Path(path).exists():return {},{"state":"not_available","entries":0}
    j=json.loads(Path(path).read_text(encoding="utf-8"))
    if j.get("schema")!="V288_CURATED_POPULATION_AF_1" or j.get("assembly")!="GRCh38":
        raise ValueError("V288_POPULATION_AF_SCHEMA_OR_BUILD_MISMATCH")
    data={}
    for row in j.get("entries",[]):
        af=_value(row.get("allele_frequency"))
        key=str(row.get("variant_key") or "")
        if (not key or af is None or not 0<=af<=1 or not row.get("source_reference")
            or not row.get("population") or row.get("source_checked") is not True
            or _value(row.get("allele_number")) is None or float(row["allele_number"])<=0):
            continue
        data[key]={
         "af":af,"an":int(float(row["allele_number"])),
         "population":row["population"],"source_reference":row["source_reference"],
         "population_context":"sourced_and_allele_specific",
        }
    return data,{"state":"local_curated_reference","entries":len(data),
                 "not_a_clinical_diagnosis":True}

def quality_flags(variant):
    gt=str(variant.get("GT") or "")
    gq=_value(variant.get("GQ"))
    dp=_value(variant.get("DP"))
    ab=_value(variant.get("allele_balance"))
    flt=str(variant.get("FILTER") or variant.get("filter") or "")
    out=[]
    if flt not in ("PASS","."):out.append("variant_filter_not_pass")
    if gq is None:out.append("genotype_quality_missing")
    elif gq<20:out.append("genotype_quality_below_20_needs_review")
    if dp is None:out.append("read_depth_missing")
    elif dp<10:out.append("read_depth_below_10_needs_review")
    if gt in HET:
        if ab is None:out.append("heterozygote_allele_balance_missing")
        elif not (0.20<=ab<=0.80):out.append("heterozygote_allele_balance_outside_20_80")
    if gt in HOM:
        if ab is None:out.append("homozygote_allele_balance_missing")
        elif ab<0.90:out.append("homozygote_allele_balance_below_90")
    if gt not in HET|HOM:out.append("genotype_not_resolved_as_standard_diploid")
    return sorted(set(out))

def inspect_allele(v,pop_index=None):
    from . import case_science
    cv=case_science.clinvar_interpretation(v.get("clinvar"))
    annotations={
      "consequence":v.get("consequence") or (v.get("best_annotation") or {}).get("consequence"),
      "impact_score":v.get("impact_score") or (v.get("best_annotation") or {}).get("severity"),
    }
    gt=str(v.get("GT") or "")
    allele_key=":".join(str(v.get(k) or "") for k in ("chrom","pos","ref","alt"))
    population=(pop_index or {}).get(allele_key)
    problems=quality_flags(v)
    if not population:problems.append("source_verified_population_af_unavailable")
    elif population["af"]>=0.01:problems.append("allele_frequency_above_1pct_requires_disease_model_review")
    if cv["state"] in ("reported_benign","conflicting","conflicting_or_mixed"):
        problems.append("clinical_assertion_not_supportive_of_pathogenicity")
    if annotations["consequence"] in LOSS_CLASSES:
        problems.append("predicted_loss_of_function_not_equivalent_to_validated_PVS1")
    clndn=(v.get("clinvar") or {}).get("CLNDN")
    return {
      "variant_ref":"SHA256:"+safe_hash(v),
      "genotype_group":"homozygous" if gt in HOM else "heterozygous" if gt in HET else "other_or_unknown",
      "phase_tag_present":bool((v.get("phase_tags") or {}).get("PS") or (v.get("phase_tags") or {}).get("PID")),
      "consequence":annotations["consequence"],
      "clinvar_assertion":cv["state"],"clinvar_review_status":cv["review"],
      "pathogenicity_weight_heuristic":cv["weight"],
      "clinvar_condition_specificity":"unverified_for_present_phenotype" if clndn else "not_provided",
      "population_af":(population or {}).get("af"),
      "population_af_source":(population or {}).get("source_reference"),
      "qc_pass_for_exploration":not any(x.startswith(("genotype_quality_below","read_depth_below",
        "variant_filter_not_pass","heterozygote_allele_balance_outside","homozygote_allele_balance_below"))
        for x in problems),
      "gaps":sorted(set(problems)),
      "acmg_amp_classification":"NOT_ESTABLISHED_BY_AUTOMATION",
      "functional_effect_observed_in_patient":False,
      "variant_coordinates_in_public_report":False,
    }

def build(private,population_path=None):
    from . import case_science
    pop,popstatus=load_local_population(population_path)
    ranks={str(g.get("gene") or "").upper():g for g in
           (private.get("genome_wide_candidate_genes") or [])}
    genes={}
    for gene,vs in sorted((private.get("genes") or {}).items()):
        g=gene.upper()
        p=(private.get("phase_evidence") or {})
        phase=p.get(g) if isinstance(p,dict) else None
        phase=phase or case_science.phase_assessment(vs)
        items=[inspect_allele({**v,"gene":g},pop) for v in vs]
        states=Counter(x["clinvar_assertion"] for x in items)
        medically_supported=states.get("reported_pathogenic",0)>0
        genetics_missing=(phase.get("state") not in ("trans_supported_by_vcf","homozygous_observed")
                    and len(vs)>1)
        if all(x["clinvar_assertion"]=="reported_benign" for x in items):
            confidence="ALLELES_CONTRADICTED_BY_CONCORDANT_BENIGN_EVIDENCE"
        elif medically_supported and not genetics_missing:
            confidence="REVIEW_STRONGER_GENETIC_EVIDENCE_NOT_DIAGNOSIS"
        elif medically_supported:
            confidence="SUPPORTED_VARIANT_BUT_PHASE_OR_SECOND_ALLELE_UNRESOLVED"
        else:
            confidence="GENETIC_CAUSALITY_UNCERTAIN"
        r=ranks.get(g) or {}
        genes[g]={
          "gene":g,"rank_score_not_diagnostic":r.get("score"),
          "phenotype_similarity_score":r.get("phenotype_score"),
          "phase_state":phase.get("state","unresolved"),
          "phase_requires_orthogonal_validation":True,
          "variants":items,
          "clinical_allele_states":dict(states),
          "gene_mechanism_confidence":confidence,
          "clinvar_disease_condition_match_unverified":True,
          "gene_disease_inheritance_verified":False,
          "protein_rescue_endpoint_measured":False,
          "missing_evidence":[
            "expert_ACMG_AMP_with_condition_and_variant_review",
            "gene_disease_inheritance_and_phenotypic_specificity",
            "orthogonal_phasing_or_segregation_when_needed",
            "protein_activity_localization_and_splicing_function",
            "read_level_validation_of_candidate_calls",
          ]+([] if popstatus["entries"] else ["population_frequency_reference_missing"]),
        }
    phen=private.get("phenotype_features") or []
    result={
      "schema":SCHEMA,"local_analysis_only":True,
      "population_frequency_source":popstatus,
      "genes":genes,"primary_hypotheses":list(genes),
      "reported_present_proband_features":sum(
         x.get("presence")=="reported_present" and x.get("subject")=="proband" for x in phen),
      "reported_absent_proband_features":sum(
         x.get("presence")=="reported_absent" and x.get("subject")=="proband" for x in phen),
      "not_full_acmg_amp":True,
      "clinical_diagnosis_confirmed":False,
      "variant_strings_sent_to_external_models":False,
      "disease_mechanism_experimentally_resolved":False,
    }
    return result

def report(model,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/"GENETIC_EXPERT_DOSSIER.json").write_text(
       json.dumps(model,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=["# Genetic expert audit — research-only",
      "Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch",
      "This is NOT an ACMG/AMP diagnosis, and gene ranks are hypotheses.",
      "",
      "Population frequency: "+model["population_frequency_source"]["state"],""]
    for gene,r in model["genes"].items():
        lines.extend([f"## {gene} — {r['gene_mechanism_confidence']}",
          f"Variants examined: {len(r['variants'])}; phase: {r['phase_state']}",
          "Variant assertions: "+", ".join(f"{k}={v}" for k,v in r["clinical_allele_states"].items()),
          "Read-level validation, disease-specific inheritance, population AF and functional rescue must be checked.",
          ""])
    (out/"GENETIC_EXPERT_SUMMARY.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return model
