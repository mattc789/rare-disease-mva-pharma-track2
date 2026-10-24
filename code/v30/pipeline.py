"""Rare Disease MVA Pharma V3.0 coordinator.

V2.8.10 remains the validated evidence-acquisition core. V3.0 adds a
mechanism-first decision layer and can also upgrade an existing completed
V2.8.10 run without repeating downloads or API work.
"""
from __future__ import annotations
from pathlib import Path
import argparse,datetime,json,traceback
from v2810 import pipeline as v2810_pipeline
from v2810.core import ROOT,RUN_ROOT,load,save,read_jsonl,now
from . import license_gate,patient_state,mechanism_engine,rescue_objectives
from . import candidate_graph,falsification,dossier,openai_reasoner,literature_attack,public_enrichment,native_pipeline

VERSION="3.1.0-dev"
PROJECT_TITLE="Rare Disease MVA Pharma — V3.1 Expanded Patient-First Rescue Engine"

def _optional_json(path,default=None):
    p=Path(path)
    return load(p) if p.is_file() else ({} if default is None else default)

def _optional_jsonl(path):
    p=Path(path)
    return list(read_jsonl(p)) if p.is_file() else []

def _validate_core(core):
    required=[
      "case_private.json","PATIENT_VARIANT_MECHANISM_MODEL.json",
      "routes_final.json","priority_pairs.json","FINAL_RESULT.json"
    ]
    missing=[x for x in required if not (core/x).is_file()]
    if missing:raise RuntimeError("V30_CORE_RUN_INCOMPLETE:"+",".join(missing))
    final=load(core/"FINAL_RESULT.json")
    if not (core/"STATUS.json").is_file() or not load(core/"STATUS.json").get("completed"):
        raise RuntimeError("V30_CORE_RUN_NOT_COMPLETED")
    return final

def enhance(core_run,outdir,enable_ai_critic=False,literature_attack_count=0,public_gene_enrichment=False):
    core=Path(core_run);out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    core_final=_validate_core(core)
    policy=license_gate.assert_default_profile()
    save(out/"V30_LICENSE_GATE.json",policy)

    private=load(core/"case_private.json")
    v28_mechanisms=load(core/"PATIENT_VARIANT_MECHANISM_MODEL.json")
    genomic=_optional_json(core/"expert_genetics/GENETIC_EXPERT_DOSSIER.json")
    protein=_optional_json(core/"expert_genetics/UNIPROT_PROTEIN_CONTEXT.json")
    specialist=_optional_json(core/"patient_first_predictions/PREDICTOR_SUMMARY.json")
    routes=load(core/"routes_final.json")
    pairs=load(core/"priority_pairs.json")
    predictions=_optional_json(core/"patient_first_predictions/PATIENT_FIRST_TOP.json",[])
    if not isinstance(predictions,list):predictions=[]

    state=patient_state.build(private,v28_mechanisms,genomic,protein,specialist)
    save(out/"PATIENT_STATE_V30.json",state)
    mechanisms=mechanism_engine.build(state)
    save(out/"MECHANISM_HYPOTHESES_V30.json",mechanisms)
    objectives=rescue_objectives.build(mechanisms)
    save(out/"RESCUE_OBJECTIVES_V30.json",objectives)

    public_gene_context=None
    if public_gene_enrichment:
        uniprot_by_gene={}
        for route in routes:
            gene=str(route.get("case_gene") or "").upper()
            if gene not in private.get("genes",{}):continue
            for target in route.get("targets",[]) or []:
                if str(target.get("gene") or "").upper()==gene and target.get("accession"):
                    uniprot_by_gene.setdefault(gene,target["accession"])
        public_gene_context=public_enrichment.enrich_genes(
          sorted(private.get("genes",{})),uniprot_by_gene,
          out/"public_gene_enrichment",online=True)

    graph=candidate_graph.build(objectives,routes,pairs,predictions)
    save(out/"CAUSAL_CANDIDATE_GRAPH_V30.json",graph)
    falsified=falsification.evaluate(graph)
    save(out/"FALSIFICATION_RESULTS_V30.json",falsified)

    attack_summary=None
    if literature_attack_count and falsified.get("candidates"):
        attack_summary=literature_attack.run(
          falsified,graph,out/"source_bound_literature_attack",
          count=int(literature_attack_count),max_docs=8)
        falsified=literature_attack.apply(falsified,attack_summary)
        save(out/"FALSIFICATION_RESULTS_V30.json",falsified)

    ai_review=None
    if enable_ai_critic and falsified.get("candidates"):
        try:
            ai_review=openai_reasoner.review(
              state,mechanisms,objectives,falsified,out/"AI_RESEARCH_CRITIC_V30.json")
        except Exception as exc:
            save(out/"AI_RESEARCH_CRITIC_ERROR_V30.json",{
              "error":type(exc).__name__,"message":str(exc)[:1000],
              "pipeline_continued":True,"deterministic_results_modified":False})
    research=dossier.build(state,mechanisms,objectives,falsified,ai_review)
    dossier.write(out,research)
    result={
      "version":VERSION,"project_title":PROJECT_TITLE,"completed_utc":now(),
      "core_v2810_run":str(core.resolve()),
      "core_v2810_version":core_final.get("version"),
      "patient_state":"PATIENT_STATE_V30.json",
      "mechanism_hypotheses":"MECHANISM_HYPOTHESES_V30.json",
      "rescue_objectives":"RESCUE_OBJECTIVES_V30.json",
      "causal_candidate_graph":"CAUSAL_CANDIDATE_GRAPH_V30.json",
      "falsification_results":"FALSIFICATION_RESULTS_V30.json",
      "research_dossier":"V30_RESEARCH_DOSSIER.json",
      "research_conclusion":research["research_conclusion"],
      "proposed_research_lead":research.get("proposed_research_lead"),
      "highest_priority_unresolved_hypothesis":research.get("highest_priority_unresolved_hypothesis"),
      "ai_critic_requested":bool(enable_ai_critic),
      "ai_critic_completed":bool(ai_review),
      "source_bound_literature_attack_requested":int(literature_attack_count or 0),
      "source_bound_literature_attack_completed":len((attack_summary or {}).get("results",[])),
      "public_gene_enrichment_requested":bool(public_gene_enrichment),
      "public_gene_enrichment_file":"public_gene_enrichment/PUBLIC_GENE_ENRICHMENT_V30.json" if public_gene_context else None,
      "public_gene_queries_exposed_variants_or_genotypes":False,
      "raw_patient_data_sent_to_v30_ai":False,
      "ai_can_upgrade_deterministic_evidence":False,
      "clinical_recommendation":False,
      "clinical_efficacy_probability_computed":False,
      "dose_recommendation":False,
      "source_policy_valid":policy["default_profile_valid"],
    }
    save(out/"FINAL_RESULT_V30.json",result)
    save(out/"STATUS_V30.json",{"completed":True,"version":VERSION,
                               "current_stage":"COMPLETED","final_result":"FINAL_RESULT_V30.json"})
    return result

def run(vcf_path=None,clinical_path=None,case_path=None,outdir=None,from_v2810_run=None,
        final_count=5,use_ai=False,functional_review_count=12,population_af_path=None,
        safety_review_count=20,gnomad_manifest=None,pk_manifest=None,vep_on_case_replay=False,
        auto_resources=True,resource_download=True,resource_max_gib=240,
        enable_v30_ai_critic=None,literature_attack_count=0,public_gene_enrichment=True,
        public_safety_review=True,network_expansion=True,variant_structure_enabled=True,
        alphamissense=False,alphamissense_download=False,modifier_discovery=True,dgidb=True):
    run_dir=Path(outdir) if outdir else RUN_ROOT/(
      "E2E_RDMVA_V30_"+datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    run_dir.mkdir(parents=True,exist_ok=True)
    try:
        if from_v2810_run:
            core=Path(from_v2810_run)
            critic=use_ai if enable_v30_ai_critic is None else enable_v30_ai_critic
            return enhance(core,run_dir/"v30",enable_ai_critic=critic,
                           literature_attack_count=literature_attack_count,
                           public_gene_enrichment=public_gene_enrichment)
        # New V3 runs are truly native: rescue objectives and V3 routes are
        # constructed before the complete approved-drug universe is screened.
        return native_pipeline.run(
          vcf_path=vcf_path,clinical_path=clinical_path,case_path=case_path,
          outdir=run_dir,enable_online_ai=use_ai,
          enable_public_gene_enrichment=public_gene_enrichment,
          literature_route_discovery=use_ai,
          literature_attack_count=literature_attack_count if use_ai else 0,
          population_af_path=population_af_path,gnomad_manifest=gnomad_manifest,
          pk_manifest=pk_manifest,auto_resources=auto_resources,
          resource_download=resource_download,resource_max_gib=resource_max_gib,
          enable_public_safety_review=public_safety_review,
          enable_network_expansion=network_expansion,
          enable_variant_structure=variant_structure_enabled,
          enable_alphamissense=alphamissense,
          alphamissense_download=alphamissense_download,
          enable_modifier_discovery=modifier_discovery,
          enable_dgidb=dgidb)
    except Exception as exc:
        save(run_dir/"FAILURE_V30.json",{"utc":now(),"error":type(exc).__name__,
          "message":str(exc)[:1200],"traceback":traceback.format_exc()[-6000:]})
        raise

def main():
    ap=argparse.ArgumentParser(description=PROJECT_TITLE)
    ap.add_argument("--vcf")
    ap.add_argument("--clinical")
    ap.add_argument("--case")
    ap.add_argument("--out")
    ap.add_argument("--from-v2810-run",help="Upgrade an existing completed V2.8.10 run without rerunning the core.")
    ap.add_argument("--final-count",type=int,default=5)
    ap.add_argument("--enable-online-ai",action="store_true",
      help="Enable V2.8.10 source-bound AI plus the V3 sanitized OpenAI critic.")
    ap.add_argument("--no-v30-ai-critic",action="store_true",
      help="Run core AI if enabled but skip the V3 sanitized critic.")
    ap.add_argument("--literature-attack-count",type=int,default=0,
      help="Source-bound Europe PMC + OpenAI evidence attacks for top HOLD candidates (0 disables).")
    ap.add_argument("--no-public-gene-enrichment",action="store_true",
      help="Disable Open Targets/Reactome gene/network queries. Enabled by default in native V3.1.")
    ap.add_argument("--no-network-expansion",action="store_true",
      help="Disable genome-wide Open Targets/Reactome network expansion.")
    ap.add_argument("--no-variant-structure",action="store_true",
      help="Disable AlphaFold DB structural context.")
    ap.add_argument("--enable-alphamissense-restricted",action="store_true",
      help="Explicitly opt in to AlphaMissense external research data despite conflicting official license notices. Disabled by default.")
    ap.add_argument("--no-alphamissense",action="store_true",
      help="Deprecated compatibility switch; AlphaMissense is already disabled by default.")
    ap.add_argument("--download-alphamissense",action="store_true",
      help="With --enable-alphamissense-restricted only: download the external GRCh38 dataset if absent.")
    ap.add_argument("--no-modifier-discovery",action="store_true",
      help="Disable drug-blind source-bound modifier discovery.")
    ap.add_argument("--no-dgidb",action="store_true",
      help="Disable DGIdb discovery-index channel.")
    ap.add_argument("--no-public-safety-review",action="store_true",
      help="Disable public FDA-label review of Track-2 single/combination hypotheses.")
    ap.add_argument("--functional-review-count",type=int,default=12)
    ap.add_argument("--safety-review-count",type=int,default=20)
    ap.add_argument("--population-af")
    ap.add_argument("--gnomad-manifest")
    ap.add_argument("--pk-manifest")
    ap.add_argument("--vep-on-replay",action="store_true")
    ap.add_argument("--no-auto-resources",action="store_true")
    ap.add_argument("--no-resource-download",action="store_true")
    ap.add_argument("--resource-max-gib",type=int,default=240)
    a=ap.parse_args()
    if not a.from_v2810_run and not a.vcf and not a.case:
        ap.error("Supply --vcf/--case, or --from-v2810-run.")
    result=run(
      vcf_path=a.vcf,clinical_path=a.clinical,case_path=a.case,outdir=a.out,
      from_v2810_run=a.from_v2810_run,final_count=a.final_count,
      use_ai=a.enable_online_ai,functional_review_count=a.functional_review_count,
      population_af_path=a.population_af,safety_review_count=a.safety_review_count,
      gnomad_manifest=a.gnomad_manifest,pk_manifest=a.pk_manifest,
      vep_on_case_replay=a.vep_on_replay,auto_resources=not a.no_auto_resources,
      resource_download=not a.no_resource_download,resource_max_gib=a.resource_max_gib,
      enable_v30_ai_critic=(a.enable_online_ai and not a.no_v30_ai_critic),
      literature_attack_count=a.literature_attack_count,
      public_gene_enrichment=not a.no_public_gene_enrichment,
      public_safety_review=not a.no_public_safety_review,
      network_expansion=not a.no_network_expansion,
      variant_structure_enabled=not a.no_variant_structure,
      alphamissense=(a.enable_alphamissense_restricted and not a.no_alphamissense),
      alphamissense_download=(a.download_alphamissense and a.enable_alphamissense_restricted),
      modifier_discovery=not a.no_modifier_discovery,
      dgidb=not a.no_dgidb)
    conclusion=(result.get("research_conclusion") or
                result.get("strict_research_conclusion") or
                ("TRACK2_HYPOTHESIS_AVAILABLE" if result.get("track2_has_hypothesis") else "NO_LEAD"))
    print(json.dumps({"completed":True,"research_conclusion":conclusion,
                      "track2_has_hypothesis":bool(result.get("track2_has_hypothesis")),
                      "clinical_recommendation":False},sort_keys=True))

if __name__=="__main__":main()
