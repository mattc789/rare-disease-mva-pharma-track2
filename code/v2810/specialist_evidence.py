"""Merge independent specialist evidence into a patient-private research dossier.

Annotative only: no source absence converts to false clinical normality,
and no externally computed score upgrades A/B or creates a drug recommendation.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from pathlib import Path
import json

SCHEMA="V289_SPECIALIST_EVIDENCE_INTEGRATION_1"

def attach(genomics,private,vep_records,splice_results,gnomad_results,pangolin_records=None):
    from .expert_genomics import safe_hash
    by_hash_sp={x.get("variant_ref"):x for x in splice_results}
    by_hash_gn={x.get("variant_ref"):x for x in gnomad_results}
    by_hash_open={x.get("variant_ref"):x for x in (pangolin_records or [])}
    mapped={(str(x.get("chrom") or "").removeprefix("chr"),x.get("pos"),
             str(x.get("ref") or "").upper(),str(x.get("alt") or "").upper()):x
            for x in vep_records}
    stats=Counter()
    for gene,record in (genomics.get("genes") or {}).items():
        for raw,item in zip(private.get("genes",{}).get(gene,[]),record.get("variants",[])):
            reference="SHA256:"+safe_hash({**raw,"gene":gene})
            vep=mapped.get((str(raw.get("chrom") or "").removeprefix("chr"),raw.get("pos"),
                 str(raw.get("ref") or "").upper(),str(raw.get("alt") or "").upper()))
            sp=by_hash_sp.get(reference)
            gn=by_hash_gn.get(reference)
            open_splice=by_hash_open.get(reference)
            q={
             "vep_annotation_state":"FOUND_LOCAL_GRCH38" if vep else "NOT_FOUND_OR_NOT_RUN",
             "vep_most_severe_consequence":(vep or {}).get("most_severe_consequence"),
             "vep_vs_local_class_discordance":bool(vep and item.get("consequence") and
                 vep.get("most_severe_consequence") and
                 vep["most_severe_consequence"]!=item["consequence"]),
             "vep_transcript_consequences_count":len((vep or {}).get("transcript_annotations") or []),
             "vep_mane_transcript_effect_independently_reviewed":False,
             "spliceai_state":(sp or {}).get("state","NOT_AVAILABLE"),
             "spliceai_max_delta_score":(sp or {}).get("maximum_delta_score"),
             "spliceai_is_functional_splicing_proof":False,
             "open_source_pangolin_state":(open_splice or {}).get("state","NOT_AVAILABLE"),
             "open_source_pangolin_max_absolute_score_change":(open_splice or {}).get("max_absolute_splice_score_change"),
             "open_source_pangolin_is_rna_experiment":False,
             "gnomad_state":(gn or {}).get("state","NOT_AVAILABLE"),
             "gnomad_max_global_af":(gn or {}).get("max_valid_global_af"),
             "gnomad_exome_genome_context":(gn or {}).get("measurements",[]),
             "gnomad_frequency_is_causality_proof":False,
             "experimentally_confirmed_variant_impact":False,
            }
            if q["vep_vs_local_class_discordance"]:
                item.setdefault("gaps",[]).append("VEP_GTF_vs_local_MANE_consequence_discordance_requires_isoform_review")
                stats["vep_local_annotation_class_discordances"]+=1
            if q["gnomad_max_global_af"] is not None and q["gnomad_max_global_af"]>=0.01:
                item.setdefault("gaps",[]).append("high_observed_population_af_review_in_disease_inheritance_context")
            if q["open_source_pangolin_max_absolute_score_change"] is not None:
                stats["pangolin_genetic_variants_with_scores"]+=1
                item.setdefault("gaps",[]).append("Pangolin_splicing_effect_requires_orthogonal_RNA_validation")
            if q["spliceai_max_delta_score"] is not None and q["spliceai_max_delta_score"]>=0.5:
                item.setdefault("gaps",[]).append("computational_splice_effect_rna_confirmation_needed")
            item["independent_specialist_annotation"]=q
            if vep:stats["vep_local_records_matched"]+=1
            if q["spliceai_max_delta_score"] is not None:stats["spliceai_variants_matched"]+=1
            if q["gnomad_max_global_af"] is not None:stats["gnomad_variants_matched"]+=1
    summary={"schema":SCHEMA,**dict(stats),
        "patient_variant_count":sum(len(v.get("variants",[])) for v in genomics["genes"].values()),
        "clinical_acmg_classification_not_created":True,
        "disease_causality_not_promoted_by_database_presence":True,
        "patient_allele_strings_not_exposed_to_external_services":True}
    return summary

def make_report(outdir,sources,pk_summary,merge_summary):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    rows=[
      "# Rare Disease MVA Pharma — Specialist Evidence Status (V2.8.9)",
      "Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch",
      "",
      "**These are data-source and model qualification results, not prescriptions or evidence of clinical drug effectiveness.**",
      "",
      "## Integrated specialist tools",
    ]
    display={
      "vep":("Ensembl VEP offline annotations","tool_state"),
      "spliceai":("SpliceAI local, licensed predictions","state"),
      "gnomad":("gnomAD v4.1 local exomes and genomes","state"),
      "pangolin":("Pangolin GPLv3 independent local splice prediction","prediction_state"),
    }
    for key,(label,field) in display.items():
        obj=sources.get(key) or {}
        rows.append(f"- **{label}:** {obj.get(field,'not run')}.")
        if key=="gnomad":
            rows.append("  - Complete coverage of 24 nuclear chromosomes in both cohorts: "+
               str(obj.get("v4_1_full_coverage",False))+".")
        if key=="spliceai":
            rows.append("  - Commercial license asserted in supplied manifest (not independently legally verified): "+
               str(obj.get("commercial_use_automatically_permitted",False))+".")
    rows += [
      f"- **Quantitative pharmacokinetic model assessment:** {pk_summary.get('state','unknown')}.",
      f"  - Official OSP midazolam technical control exported and parsed: {((pk_summary.get('osp_midazolam_public_engineering_control') or {}).get('simulations_with_parsable_concentrations',0))} simulations. Not a case drug candidate.",
      f"  - Source-qualified model profiles: {pk_summary.get('models_source_qualified',0)}.",
      "",
      "## Interpretation",
      f"- Specialist annotations matched to selected variants: {merge_summary.get('patient_variant_count',0)} input variants.",
      "- No genomic allele strings were sent to web services.",
      "- No population frequency is assumed to be zero when the source is absent.",
      "- Model absence is NOT lack of drug activity.",
      "- VEP with local GENCODE GTF+GRCh38 FASTA is not equivalent to a full Ensembl VEP variation cache.",
      "- Source-qualified PK plasma profiles do not establish tissue exposure, pediatric safety or drug rescue.",
      "- Clinical ACMG/AMP/ClinGen assessment, drug-specific functional rescue and independent experimental replication remain necessary.",
    ]
    (out/"SPECIALIST_EVIDENCE_STATUS.md").write_text("\n".join(rows)+"\n",encoding="utf-8")
    summary={"schema":SCHEMA,"source_status":sources,"pk_model_summary":pk_summary,
       "merged_variant_evidence":merge_summary,
       "none_of_these_outputs_are_therapy_recommendations":True}
    (out/"SPECIALIST_EVIDENCE_STATUS.json").write_text(
       json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return summary
