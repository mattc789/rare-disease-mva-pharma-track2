"""True V3.0 end-to-end patient-first drug-rescue pipeline.

Key invariant: Patient State -> Mechanism -> Rescue Objectives -> V3 routes
are built BEFORE the approved-drug universe is screened.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import datetime,hashlib,json,os,re,traceback,time

from v2810.core import ROOT,RUN_ROOT,RAW,load,save,jsonl,read_jsonl,now,normalize_case
from v2810 import (intake,case_science,expert_genomics,expert_protein,
                   resource_bootstrap,vep_local,gnomad_local,open_splicing,
                   specialist_evidence,master,registry_output,screen,science,
                   ai as v28_ai,routes as v28_routes,quant_pk)
from . import (patient_state,mechanism_engine,rescue_objectives,native_routes,
               evidence_ranker,candidate_graph,falsification,combination_engine,
               challenge_hypotheses,public_enrichment,literature_attack,
               openai_reasoner,dossier,license_gate,safety_review,
               disease_state_graph,network_expansion,variant_structure,
               modifier_discovery,dgidb_discovery,polypharmacology,
               signed_network_expansion,fast_screen,direct_drug_rescue,
               cross_source_bridge,exposure_feasibility,challenge_dossier,
               rights_audit)

VERSION="3.1.0-dev-native"
SCHEMA="RDMVA_NATIVE_PIPELINE_V31_1"

def _safe_save(path,obj,retries=12,delay=0.5):
    """Extra retry layer for transient Windows/AV/indexer file locks."""
    last=None
    for attempt in range(retries):
        try:
            save(path,obj);return
        except PermissionError as exc:
            last=exc;time.sleep(delay*(1+attempt/4))
    raise last


def _stage(run,name,state,detail=None):
    p=Path(run)/"STATUS_V30_NATIVE.json"
    s=load(p) if p.exists() else {"version":VERSION,"stages":{},"completed":False}
    s["stages"][name]={"state":state,"utc":now(),"detail":detail}
    s["current_stage"]=name
    _safe_save(p,s)
    if state=="running":
        print("V31_NATIVE_STAGE_START",name,flush=True)
    elif state=="completed":
        print("V31_NATIVE_STAGE_DONE",name,flush=True)

def verify_required_genomic_evidence(private,vep_setup,vep_status,
                                     gnomad_audit,gnomad_status,gnomad_rows):
    # Case-targeted local gnomAD shards are not a claim of *whole-genome*
    # gnomAD coverage. Validate each cohort for every needed case chromosome.
    from v2810.resource_bootstrap import needed_chromosomes
    chroms,unsupported=needed_chromosomes(private)
    coverage=gnomad_status.get("cohort_chromosome_completeness") or {}
    check={
      "vep_setup_ready":vep_setup.get("state")=="READY",
      "vep_executed_offline":vep_status.get("tool_state")=="COMPLETED" and
           vep_status.get("tool_provenance",{}).get("run_command_kind")=="wsl_local_perl_offline",
      "vep_annotation_output_present":int(vep_status.get("annotation_records_returned") or 0)>0,
      "gnomad_all_case_shards_downloaded":gnomad_audit.get("completion")=="READY" and
           int(gnomad_audit.get("required_jobs") or 0)>0 and
           int(gnomad_audit.get("installed_jobs") or 0)==int(gnomad_audit.get("required_jobs") or 0),
      "gnomad_selected_chromosomes_supported":bool(chroms) and not unsupported,
      "gnomad_case_targeted_exomes_and_genomes_present":bool(chroms) and all(
           set(chroms)<=set(coverage.get(cohort,[]))
           for cohort in ("exomes","genomes")),
      "gnomad_allele_lookups_have_no_index_failures":bool(gnomad_rows) and all(
           all(item.get("status") not in
               ("INDEX_QUERY_FAILED","DATASET_OR_CHROMOSOME_MISSING")
               for item in row.get("measurements",[]))
           for row in gnomad_rows),
      "gnomad_full_24_chromosome_both_cohorts":bool(gnomad_status.get("v4_1_full_coverage"))
    }
    return {
      "schema":"RDMVA_CASE_GENOMIC_EVIDENCE_AUDIT_V31_1",
      "required_science_pass":all(v for k,v in check.items()
                                  if k!="gnomad_full_24_chromosome_both_cohorts"),
      "checks":check,
      "case_chromosome_cohort_count":len(chroms)*2,
      "case_focused_shards_not_genome_wide_full_coverage":True,
      "patient_coordinates_exported":False,
      "vep_offline_required":True}

def audit_supplied_gnomad_manifest(private,manifest_path):
    """Case-scoped audit of a predownloaded local gnomAD v4.1 VCF manifest.

    An explicitly supplied local manifest bypasses remote gnomAD download,
    but is never treated as verified merely because it exists. Validate
    all required cohort/chromosome pairs, exact public source VCF SHA-256
    and local index SHA-256. Never claim 48-shard genome-wide coverage
    when only the current case's chromosome shards are present.
    """
    chroms,unsupported=resource_bootstrap.needed_chromosomes(private)
    if not chroms or unsupported:
        raise RuntimeError("V31_SUPPLIED_GNOMAD_CASE_CHROMOSOMES_INVALID")
    manifest_path=Path(manifest_path)
    shards,summary=gnomad_local.validate_manifest(manifest_path)
    required={(cohort,chrom) for chrom in chroms for cohort in ("exomes","genomes")}
    missing=required-set(shards)
    if missing:
        raise RuntimeError(
            "V31_SUPPLIED_GNOMAD_MISSING_CASE_SHARDS:"+str(len(missing))+
            "_OF_"+str(len(required))+
            ":Provide indexed exome AND genome data for every case chromosome")
    original=json.loads(manifest_path.read_text(encoding="utf-8"))
    entries={(str(item.get("cohort")),str(item.get("chrom")).removeprefix("chr").upper()):item
             for item in original.get("shards",[])}
    rows=[]
    for cohort,chrom in sorted(required):
        item=entries.get((cohort,chrom)) or {}
        shard=shards[(cohort,chrom)]
        public_file=Path(shard["file"])
        index_file=Path(str(public_file)+".tbi")
        if not index_file.is_file():
            index_file=Path(str(public_file)+".csi")
        expected=item.get("source_sha256")
        index_expected=item.get("index_sha256")
        if (not isinstance(expected,str) or
                not re.fullmatch("[a-fA-F0-9]{64}",expected) or
                not isinstance(index_expected,str) or
                not re.fullmatch("[a-fA-F0-9]{64}",index_expected)):
            raise RuntimeError(
                "V31_SUPPLIED_GNOMAD_NEEDS_SOURCE_AND_INDEX_SHA256_PINS:"+
                cohort+":"+chrom)
        if not public_file.is_file() or not index_file.is_file():
            raise RuntimeError("V31_SUPPLIED_GNOMAD_SHARD_OR_INDEX_MISSING")
        actual=resource_bootstrap._sha256_native_file(
            public_file,"gnomAD_"+cohort+"_chr"+chrom)
        if actual.lower()!=expected.lower():
            raise RuntimeError(
                "V31_SUPPLIED_GNOMAD_SOURCE_SHA256_MISMATCH:"+cohort+":"+chrom)
        with index_file.open("rb") as inp:
            digest=hashlib.sha256()
            for block in iter(lambda:inp.read(1<<20),b""):
                digest.update(block)
        if digest.hexdigest()!=index_expected.lower():
            raise RuntimeError(
                "V31_SUPPLIED_GNOMAD_INDEX_SHA256_MISMATCH:"+cohort+":"+chrom)
        rows.append({"cohort":cohort,"chrom":chrom,
                     "status":"EXPLICIT_LOCAL_SOURCE_AND_INDEX_SHA256_VERIFIED",
                     "source_sha256":actual.lower(),"index_sha256":digest.hexdigest(),
                     "bytes":public_file.stat().st_size})
        print("V31_GNOMAD_LOCAL_CACHE_VERIFIED",cohort,
              "chr"+chrom,"source_and_index_SHA256_PASS",flush=True)
    return {
        "schema":"V2810_PREVERIFIED_LOCAL_GNOMAD_AUDIT_V31_1",
        "mode":"EXPLICIT_LOCAL_MANIFEST_NO_DOWNLOAD",
        "completion":"READY",
        "required_jobs":len(required),"installed_jobs":len(rows),
        "selected_case_chromosome_count":len(chroms),
        "local_indexed_shards_available":summary.get("shards_found"),
        "v4_1_full_coverage":summary.get("v4_1_full_coverage",False),
        "no_remote_variant_lookups":True,"downloaded_gib_this_run":0,
        "patient_variant_coordinates_sent_to_public_hosts":False,
        "manifest_source_path":str(manifest_path),
        "records":rows}


def _discover_literature_routes(private,run,max_docs=15):
    discovered=[]
    for gene in sorted(private.get("genes") or {}):
        aliases=v28_routes.search_aliases(gene,run/"sources")
        docs=v28_ai.literature(gene,run/"mechanism_discovery"/gene,
                               max_docs=max_docs,targets=aliases[1:])
        resolver=lambda g:v28_routes.accession_for_gene(g,run/"sources")
        wave1=v28_ai.discover_routes(
            gene,docs,run/"mechanism_discovery"/gene,resolver,
            tag="v30_route_discovery_wave1")
        exclude=[x.get("target_gene") for x in wave1 if x.get("target_gene")]
        wave2=v28_ai.discover_routes(
            gene,docs,run/"mechanism_discovery"/gene,resolver,
            exclude_targets=exclude,tag="v30_route_discovery_wave2")
        seen=set();rows=[]
        for x in wave1+wave2:
            key=(x.get("target_accession"),x.get("desired_action"),
                 x.get("variable"),x.get("source_id"))
            if key in seen:continue
            seen.add(key);rows.append(x)
        rows=v28_ai.verify_routes(gene,rows,run/"mechanism_discovery"/gene)
        for x in rows:
            x["case_gene"]=gene
        discovered.extend(rows)
    return discovered

def run(*,vcf_path=None,clinical_path=None,case_path=None,outdir=None,
        enable_online_ai=False,enable_public_gene_enrichment=True,
        literature_route_discovery=True,literature_attack_count=8,
        population_af_path=None,gnomad_manifest=None,pk_manifest=None,
        auto_resources=True,resource_download=True,resource_max_gib=240,
        enable_public_safety_review=True,enable_network_expansion=True,
        enable_variant_structure=True,enable_alphamissense=False,
        alphamissense_download=False,enable_modifier_discovery=True,
        enable_dgidb=True,enable_cross_source_bridge=True):
    run=Path(outdir) if outdir else RUN_ROOT/(
      "E2E_RDMVA_V30_NATIVE_"+datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    run.mkdir(parents=True,exist_ok=True)
    old_enrich=os.environ.get("RDMVA_ENRICHR_MODE")
    os.environ["RDMVA_ENRICHR_MODE"]="disabled"
    try:
        _stage(run,"LICENSE_GATE","running")
        policy=license_gate.assert_default_profile()
        save(run/"V30_LICENSE_GATE.json",policy)
        rights=rights_audit.write(
            run/"V31_THIRD_PARTY_RIGHTS_AUDIT.json",
            run/"THIRD_PARTY_NOTICES_V31.md")
        _stage(run,"LICENSE_GATE","completed",{
          **policy,
          "rights_audit_pass":rights.get("challenge_public_profile_pass"),
          "rights_audit_file":"V31_THIRD_PARTY_RIGHTS_AUDIT.json"})

        if auto_resources:
            _stage(run,"PUBLIC_BASELINE","running")
            baseline=resource_bootstrap.ensure_reproducible_public_baseline(
                ROOT,run/"resource_bootstrap/baseline",install=resource_download)
            _stage(run,"PUBLIC_BASELINE","completed",{"state":baseline.get("state")})
            if baseline.get("state")!="READY":
                raise RuntimeError("V30_PUBLIC_BASELINE_NOT_READY")

        _stage(run,"PATIENT_INTAKE","running")
        if case_path:
            case_path=Path(case_path)
        else:
            if not vcf_path:raise ValueError("V30_VCF_OR_CASE_REQUIRED")
            intake.build_case_model(Path(vcf_path),
                                    Path(clinical_path) if clinical_path else None,
                                    run/"intake")
            case_path=run/"intake/case_model_v2.json"
        private,public=normalize_case(case_path)
        save(run/"case_private.json",private);save(run/"case_public_safe.json",public)
        if not private.get("genes"):raise RuntimeError("V30_NO_GENE_HYPOTHESIS_FROM_INTAKE")
        _stage(run,"PATIENT_INTAKE","completed",
               {"gene_hypotheses":len(private["genes"]),"raw_data_external":False})

        _stage(run,"GENOMIC_EVIDENCE","running")
        v28_mechanisms=case_science.build_case_science(private)
        save(run/"PATIENT_VARIANT_MECHANISM_MODEL_V28_COMPAT.json",v28_mechanisms)
        genomic=expert_genomics.build(private,population_af_path)
        protein=expert_protein.build_from_private(private,genomic,RAW)

        local_gnomad=gnomad_manifest
        gnomad_audit={"completion":"NOT_ATTEMPTED","required_jobs":0,"installed_jobs":0}
        vep_setup={"state":"UNASSESSED","vep_config":{}}
        if auto_resources:
            vep_setup=resource_bootstrap.ensure_vep_prerequisites(
                ROOT,run/"resource_bootstrap/vep",install=resource_download)
            if resource_download and vep_setup.get("state")!="READY":
                raise RuntimeError("V31_REQUIRED_VEP_SETUP_INCOMPLETE: "+
                                   str(vep_setup.get("failure_reason") or "no usable VEP")[:300])
            if local_gnomad:
                gnomad_audit=audit_supplied_gnomad_manifest(private,local_gnomad)
                _safe_save(run/"resource_bootstrap/gnomad/GNOMAD_AUTO_RESOURCE_AUDIT.json",
                           gnomad_audit)
                print("V31_GNOMAD_REUSED_VERIFIED_LOCAL_CACHE_NO_DOWNLOAD",
                      gnomad_audit["installed_jobs"],flush=True)
            else:
                local_gnomad,gnomad_audit=resource_bootstrap.prepare(
                    private,run/"resource_bootstrap/gnomad",
                    download=resource_download,max_gib=resource_max_gib)
                if resource_download and gnomad_audit.get("completion")!="READY":
                    raise RuntimeError("V31_REQUIRED_GNOMAD_INCOMPLETE: "+
                        str(gnomad_audit.get("error_summary") or
                            gnomad_audit.get("completion"))[:260])
        vep_status,vep_rows=vep_local.annotation_stage(
            private,run/"specialist_genomics/vep",
            config=vep_setup.get("vep_config") if vep_setup.get("state")=="READY" else None)
        gnomad_status,gnomad_rows=gnomad_local.run(
            private,run/"specialist_genomics/gnomad",manifest_path=local_gnomad)
        if auto_resources and resource_download:
            genomic_audit=verify_required_genomic_evidence(
                private,vep_setup,vep_status,gnomad_audit,gnomad_status,gnomad_rows)
            _safe_save(run/"GENOMIC_REQUIRED_EVIDENCE_V31.json",genomic_audit)
            if not genomic_audit["required_science_pass"]:
                failed=sorted(k for k,v in genomic_audit["checks"].items()
                              if not v and k!="gnomad_full_24_chromosome_both_cohorts")
                raise RuntimeError("V31_REQUIRED_GENOMIC_EVIDENCE_FAILED:"+",".join(failed))
        if auto_resources:
            pangolin_status=open_splicing.score_case(
                private,run/"specialist_genomics/pangolin_gpl3",
                install=resource_download,run_if_ready=resource_download)
        else:
            pangolin_status={"prediction_state":"SKIPPED"}
        pangolin_rows=[]
        pf=run/"specialist_genomics/pangolin_gpl3/PANGOLIN_PRIVATE_HASHED_SCORES.jsonl"
        if pf.is_file():
            pangolin_rows=list(__import__("v2810.core",fromlist=["read_jsonl"]).read_jsonl(pf))
        integrated=specialist_evidence.attach(
            genomic,private,vep_rows,[],gnomad_rows,pangolin_records=pangolin_rows)
        expert_protein.save_report(protein,run/"expert_genetics")
        expert_genomics.report(genomic,run/"expert_genetics")
        _stage(run,"GENOMIC_EVIDENCE","completed",{
          "vep":vep_status.get("tool_state"),"gnomad":gnomad_status.get("state"),
          "pangolin":pangolin_status.get("prediction_state"),
          "patient_variants_external":False})

        # The defining V3 stage: derive biology BEFORE seeing drug names.
        _stage(run,"PATIENT_STATE_MECHANISM_OBJECTIVES","running")
        state=patient_state.build(private,v28_mechanisms,genomic,protein,integrated)
        save(run/"PATIENT_STATE_V30.json",state)
        mechanisms=mechanism_engine.build(state)
        save(run/"MECHANISM_HYPOTHESES_V30.json",mechanisms)
        objectives=rescue_objectives.build(mechanisms)

        # V3.1: preserve primary directional hypotheses while expanding the
        # genome-wide causal state used for network/modifier discovery.
        dsg=disease_state_graph.build(
            private,state,mechanisms,objectives,signor_edges=v28_routes.signor_edges())
        _safe_save(run/"DISEASE_STATE_GRAPH_V31.json",dsg)

        uniprot={}
        for gene in dsg.get("network_seed_genes",[]):
            try:uniprot[gene]=v28_routes.accession_for_gene(gene,run/"sources")
            except Exception:pass

        structure_context={"schema":"RDMVA_VARIANT_STRUCTURE_V31_1","genes":{}}
        alphamissense={"schema":"RDMVA_ALPHAMISSENSE_V31_1","variants":[],"matched":0}
        if enable_variant_structure:
            _stage(run,"VARIANT_STRUCTURE_CONTEXT","running")
            structure_context=variant_structure.structural_variant_context(
                private,protein,{g:uniprot.get(g) for g in dsg.get("primary_genes",[])},
                run/"variant_structure",online=bool(enable_public_gene_enrichment))
            if enable_alphamissense:
                am_resource=variant_structure.ensure_alphamissense(
                    ROOT,download=bool(alphamissense_download))
                alphamissense=variant_structure.lookup_alphamissense(
                    private,am_resource,run/"variant_structure")
            objectives=variant_structure.annotate_objectives(
                objectives,structure_context,alphamissense)
            _stage(run,"VARIANT_STRUCTURE_CONTEXT","completed",{
              "afdb_genes":sum((x.get("afdb_state")=="RESOLVED")
                              for x in structure_context.get("genes",{}).values()),
              "alphamissense_matched":alphamissense.get("matched",0),
              "clinical_validation":False})

        _safe_save(run/"RESCUE_OBJECTIVES_V30.json",objectives)
        _safe_save(run/"RESCUE_OBJECTIVES_V31.json",objectives)
        _stage(run,"PATIENT_STATE_MECHANISM_OBJECTIVES","completed",{
          "genes":len(state.get("genes",[])),
          "genome_wide_network_seed_genes":len(dsg.get("network_seed_genes",[])),
          "mechanisms":len(mechanisms.get("hypotheses",[])),
          "objectives":len(objectives.get("objectives",[])),
          "drug_names_used_to_build_objectives":False})

        network={"schema":"RDMVA_NETWORK_EXPANSION_V31_1","nodes":[],"edges":[],
                 "modifier_candidates":[],"modifier_candidate_count":0}
        if enable_network_expansion:
            _stage(run,"NETWORK_PATHWAY_EXPANSION","running")
            network=network_expansion.build(
                dsg,uniprot,run/"network_expansion",
                online=bool(enable_public_gene_enrichment))
            _stage(run,"NETWORK_PATHWAY_EXPANSION","completed",{
              "seed_genes":len(dsg.get("network_seed_genes",[])),
              "nodes":len(network.get("nodes",[])),
              "edges":len(network.get("edges",[])),
              "modifier_candidates":network.get("modifier_candidate_count",0),
              "direction_inferred_from_unsigned_network":False})

        _stage(run,"V31_NATIVE_ROUTE_GENERATION","running")
        route_bundle=native_routes.build(private,objectives,run/"sources")
        discovered=[]
        if enable_online_ai and literature_route_discovery:
            discovered=_discover_literature_routes(private,run)
            route_bundle=native_routes.merge_literature_discovered(route_bundle,discovered,objectives=objectives)

        signed_bundle=signed_network_expansion.build(
            private,objectives,run/"sources",max_hops=2)
        route_bundle=native_routes.merge_objective_routes(
            route_bundle,signed_bundle.get("routes",[]),label="v31_signed_network")

        modifier_bundle={"routes":[],"route_count":0,"summaries":[]}
        if enable_online_ai and enable_modifier_discovery and enable_network_expansion:
            modifier_findings=modifier_discovery.discover(
                dsg.get("primary_genes",[]),network,
                run/"modifier_discovery",max_docs=30,max_network_targets=20,
                objectives=objectives)
            modifier_routes=modifier_discovery.to_routes(
                modifier_findings,run/"sources",objectives)
            modifier_bundle={"routes":modifier_routes,
                             "route_count":len(modifier_routes),
                             "summaries":modifier_findings.get("accepted_modifiers",[])}
            route_bundle=native_routes.merge_objective_routes(
                route_bundle,modifier_routes,label="v31_modifier_discovery")

        routes=route_bundle["routes"]
        _safe_save(run/"NATIVE_ROUTES_V30.json",route_bundle)
        _safe_save(run/"NATIVE_ROUTES_V31.json",route_bundle)
        _stage(run,"V31_NATIVE_ROUTE_GENERATION","completed",{
          "routes":len(routes),"legacy_literature_discovered":len(discovered),
          "signed_network_routes":signed_bundle.get("route_count",0),
          "drug_blind_modifier_routes":modifier_bundle.get("route_count",0),
          "derived_before_drug_screen":True})

        if enable_public_gene_enrichment and not enable_network_expansion:
            public_enrichment.enrich_genes(
                sorted(private["genes"]),uniprot,run/"public_gene_enrichment",online=True)

        _stage(run,"APPROVED_DRUG_UNIVERSE","running")
        entities,products,master_summary=master.build(run/"master")
        registry_index=registry_output.build_index(products)
        if not entities or len(entities)!=master_summary.get("entities",len(entities)):
            raise RuntimeError("V30_APPROVED_DRUG_UNIVERSE_COVERAGE_MISMATCH")
        _stage(run,"APPROVED_DRUG_UNIVERSE","completed",{
          "entities":len(entities),"products":len(products),
          "single_approved":sum(bool(x.get("global_approved_as_single") or x.get("solo_approved") or x.get("biologic_licensed")) for x in entities)})

        # Public ChEMBL target responses can be recovered only from a
        # previously rebuilt, byte-identical regulatory/ChEMBL universe.
        # Never reuse original genomic annotations or clinical conclusions.
        previous_public_cache=os.environ.get("RDMVA_RECOVERY_CHEMBL_FROM")
        if previous_public_cache:
            fast_screen.restore_public_chembl_cache(previous_public_cache,run)
        # ALL approved-universe entities pass through the V3 route screen.
        _stage(run,"V30_UNIVERSAL_SCREEN","running")
        ledger,screen_summary=fast_screen.universal_screen(entities,routes,run/"screen",max_workers=8)
        if len(ledger)!=len(entities):raise RuntimeError("V30_UNIVERSAL_SCREEN_INCOMPLETE")
        grouped,group_summary=science.canonicalize_ledger(ledger)
        registry_index=science.merge_registry_by_family(grouped,registry_index)
        jsonl(run/"screen/scientific_parent_ledger.jsonl",grouped)
        save(run/"screen/SCIENTIFIC_GROUP_SUMMARY.json",group_summary)
        pairs=evidence_ranker.prioritize(grouped,routes)

        dgidb_payload={"hints":[],"hint_count":0};dgidb_audit={"added_pairs":0}
        if enable_dgidb:
            try:
                dgidb_payload=dgidb_discovery.discover(routes,grouped,run/"dgidb_discovery")
                pairs,dgidb_audit=dgidb_discovery.merge_pairs(pairs,dgidb_payload,routes,grouped)
            except Exception as exc:
                _safe_save(run/"dgidb_discovery/DGIDB_ERROR_V31.json",{
                  "error":type(exc).__name__,"message":str(exc)[:400],
                  "aggregator_claims_are_submission_evidence":False})

        pairs,poly_summary=polypharmacology.build(grouped,pairs,run/"polypharmacology")

        direct_payload=direct_drug_rescue.discover(
            run,entities,sorted(private.get("genes") or {}),objectives,
            run/"direct_drug_rescue")
        if enable_online_ai and direct_payload.get("findings"):
            direct_payload=direct_drug_rescue.adjudicate(
                direct_payload,objectives,
                run/"direct_drug_rescue/DIRECT_DRUG_RESCUE_ADJUDICATED_V31.json",
                run/"direct_drug_rescue")
        else:
            direct_payload={**direct_payload,"accepted_findings":[],
                            "rejected_findings":direct_payload.get("findings",[]),
                            "adjudication_completed":False}
        pairs,direct_audit=direct_drug_rescue.merge_pairs(
            pairs,direct_payload,objectives)

        bridge_payload={"bridges":[],"bridge_count":0,"targets_scouted":0}
        bridge_audit={"added_pairs":0,"routes":[]}
        if enable_online_ai and enable_cross_source_bridge:
            bridge_payload=cross_source_bridge.discover_from_routes(
                routes,entities,objectives,run/"cross_source_bridge_strict_batch",
                max_targets=60,max_full_bridges=30)
            pairs,bridge_audit=cross_source_bridge.merge_pairs(
                pairs,bridge_payload)

        graph_routes=(list(routes)+list(direct_audit.get("routes") or [])+
                      list(bridge_audit.get("routes") or []))
        if direct_audit.get("routes") or bridge_audit.get("routes"):
            discovery_route_bundle=dict(route_bundle)
            discovery_route_bundle["routes"]=graph_routes
            discovery_route_bundle["direct_source_drug_rescue_route_count"]=len(
                direct_audit.get("routes") or [])
            discovery_route_bundle["cross_source_bridge_route_count"]=len(
                bridge_audit.get("routes") or [])
            _safe_save(run/"NATIVE_ROUTES_V31_WITH_DISCOVERY_BRIDGES.json",
                       discovery_route_bundle)

        _safe_save(run/"V30_ROUTE_DRUG_PAIRS.json",pairs)
        _safe_save(run/"V31_ROUTE_DRUG_PAIRS.json",pairs)
        _stage(run,"V30_UNIVERSAL_SCREEN","completed",{
          "all_entities_screened":len(ledger),
          "signal_entities":screen_summary.get("signal_entities"),
          "route_drug_pairs":len(pairs),
          "dgidb_discovery_hints":dgidb_payload.get("hint_count",0),
          "dgidb_added_direction_matched_pairs":dgidb_audit.get("added_pairs",0),
          "polypharmacology_profiled_entities":poly_summary.get("profiled_entity_count",0),
          "direct_source_drug_rescue_accepted":direct_payload.get("accepted_count",0),
          "cross_source_targets_scouted":bridge_payload.get("targets_scouted",0),
          "cross_source_bridges":bridge_payload.get("bridge_count",0),
          "cross_source_pairs_added":bridge_audit.get("added_pairs",0),
          "screen_status_counts":screen_summary.get("status_counts")})

        _stage(run,"CAUSAL_GRAPH_AND_COMBINATIONS","running")
        graph=candidate_graph.build(objectives,graph_routes,pairs,[])
        save(run/"CAUSAL_CANDIDATE_GRAPH_V30.json",graph)
        falsified=falsification.evaluate(graph)
        save(run/"FALSIFICATION_RESULTS_V30.json",falsified)
        combos=combination_engine.discover(graph,objectives,max_components=3,shortlist=100,max_results=80)
        save(run/"COMBINATION_HYPOTHESES_RAW_V30.json",combos)
        track2=challenge_hypotheses.build(graph,objectives,combos)
        safety_reviewed=0
        raw_labels={}
        if enable_public_safety_review:
            names={x["candidate"] for x in track2.get("submission_level_single_agents",[])}
            for c in track2.get("combination_hypotheses",[]):
                names.update(c.get("components") or [])
            if names:
                safety_payload,raw_labels=safety_review.review_candidates(
                    sorted(names),grouped,private,run/"safety_review")
                track2=safety_review.annotate_hypotheses(track2,safety_payload,raw_labels)
                safety_reviewed=len(names)
        exposure_assessments=[]
        if enable_online_ai and raw_labels:
            track2,exposure_assessments=exposure_feasibility.assess_track2(
                track2,grouped,raw_labels,run/"exposure_feasibility")
        _safe_save(run/"TRACK2_HYPOTHESES_V30.json",track2)
        _stage(run,"CAUSAL_GRAPH_AND_COMBINATIONS","completed",{
          "candidate_paths":len(graph.get("candidate_paths",[])),
          "single_submission_hypotheses":len(track2.get("submission_level_single_agents",[])),
          "combination_hypotheses":len(track2.get("combination_hypotheses",[])),
          "public_label_safety_candidates_reviewed":safety_reviewed,
          "exposure_feasibility_candidates_assessed":len(exposure_assessments)})

        attack_summary=None
        if enable_online_ai and literature_attack_count>0:
            _stage(run,"SOURCE_BOUND_ATTACK","running")
            # Attack the candidates that actually matter for Track 2 first,
            # including components needed by the best combinations.
            priority=[]
            for x in track2.get("submission_level_single_agents",[]):
                if x.get("candidate") not in priority:priority.append(x.get("candidate"))
            for c in track2.get("combination_hypotheses",[]):
                for name in c.get("components") or []:
                    if name not in priority:priority.append(name)
            rank={name:i for i,name in enumerate(priority)}
            attack_input=dict(falsified)
            attack_input["candidates"]=sorted(
                list(falsified.get("candidates",[])),
                key=lambda x:(rank.get(x.get("candidate"),10**6),
                              -float(x.get("research_priority_score") or -999),
                              str(x.get("candidate") or "")))
            attack_summary=literature_attack.run(
                attack_input,graph,run/"source_bound_literature_attack",
                count=literature_attack_count,max_docs=10)
            falsified=literature_attack.apply(falsified,attack_summary)
            save(run/"FALSIFICATION_RESULTS_V30.json",falsified)
            track2=challenge_hypotheses.apply_source_bound_evidence(track2,attack_summary)
            _safe_save(run/"TRACK2_HYPOTHESES_V30_FINAL.json",track2)
            _stage(run,"SOURCE_BOUND_ATTACK","completed",{
              "attacked":len(attack_summary.get("results",[])),
              "track2_priority_candidates_requested":min(len(priority),literature_attack_count)})

        ai_review=None
        if enable_online_ai and falsified.get("candidates"):
            try:
                ai_review=openai_reasoner.review(
                    state,mechanisms,objectives,falsified,run/"AI_RESEARCH_CRITIC_V30.json")
            except Exception as exc:
                save(run/"AI_RESEARCH_CRITIC_ERROR_V30.json",{
                  "error":type(exc).__name__,"message":str(exc)[:800],
                  "deterministic_results_modified":False})

        # Optional PK models are evaluated separately and never invented.
        pk_status,_=quant_pk.run(
            run/"specialist_pharmacology",manifest_path=pk_manifest,
            approved_drug_names=[e["canonical_label"] for e in entities])

        strict_dossier=dossier.build(state,mechanisms,objectives,falsified,ai_review)
        dossier.write(run/"research_dossier",strict_dossier)
        track2_dossier_meta=_write_track2_dossier(run,track2,state)
        result={
          "schema":SCHEMA,"version":VERSION,"completed_utc":now(),
          "track2_research_dossier":track2_dossier_meta["json"],
          "track2_research_dossier_markdown":track2_dossier_meta["markdown"],
          "primary_submission_candidate":track2_dossier_meta["primary_candidate"],
          "primary_submission_evidence_class":track2_dossier_meta["primary_evidence_class"],
          "patient_state":"PATIENT_STATE_V30.json",
          "mechanism_hypotheses":"MECHANISM_HYPOTHESES_V30.json",
          "rescue_objectives":"RESCUE_OBJECTIVES_V31.json",
          "disease_state_graph":"DISEASE_STATE_GRAPH_V31.json",
          "network_expansion":"network_expansion/NETWORK_EXPANSION_V31.json" if enable_network_expansion else None,
          "variant_structure_context":"variant_structure/VARIANT_STRUCTURE_CONTEXT_V31.json" if enable_variant_structure else None,
          "alphamissense_context":"variant_structure/ALPHAMISSENSE_V31.json" if (enable_variant_structure and enable_alphamissense) else None,
          "native_routes":"NATIVE_ROUTES_V31.json",
          "drug_blind_modifier_route_count":modifier_bundle.get("route_count",0),
          "network_modifier_candidate_count":network.get("modifier_candidate_count",0),
          "dgidb_discovery_hint_count":dgidb_payload.get("hint_count",0),
          "dgidb_added_pair_count":dgidb_audit.get("added_pairs",0),
          "polypharmacology_profiled_entities":poly_summary.get("profiled_entity_count",0),
          "all_approved_entities_screened":len(ledger),
          "approved_product_records":len(products),
          "route_drug_pairs":len(pairs),
          "candidate_paths":len(graph.get("candidate_paths",[])),
          "strict_research_conclusion":strict_dossier.get("research_conclusion"),
          "validated_research_lead":strict_dossier.get("proposed_research_lead"),
          "track2_submission_hypotheses":"TRACK2_HYPOTHESES_V30_FINAL.json" if (run/"TRACK2_HYPOTHESES_V30_FINAL.json").is_file() else "TRACK2_HYPOTHESES_V30.json",
          "track2_single_agent_hypothesis_count":len(track2.get("submission_level_single_agents",[])),
          "track2_combination_hypothesis_count":len(track2.get("combination_hypotheses",[])),
          "track2_has_hypothesis":track2.get("has_submission_level_hypothesis"),
          "combination_search_max_components":3,
          "combination_synergy_assumed":False,
          "public_label_safety_review_enabled":bool(enable_public_safety_review),
          "clinical_safety_clearance":False,
          "market_approved_medications_only_for_track2":True,
          "pk_models_source_qualified":pk_status.get("models_source_qualified",0),
          "openai_enabled":bool(enable_online_ai),
          "patient_raw_data_sent_to_ai":False,
          "clinical_recommendation":False,
          "clinical_efficacy_probability_computed":False,
          "dose_recommendation":False,
        }
        save(run/"FINAL_RESULT_V30_NATIVE.json",result)
        s=load(run/"STATUS_V30_NATIVE.json");s["completed"]=True
        s["current_stage"]="COMPLETED";s["final_result"]="FINAL_RESULT_V30_NATIVE.json"
        save(run/"STATUS_V30_NATIVE.json",s)
        return result
    except Exception as exc:
        save(run/"FAILURE_V30_NATIVE.json",{
          "utc":now(),"error":type(exc).__name__,"message":str(exc)[:1200],
          "traceback":traceback.format_exc()[-7000:]})
        raise
    finally:
        if old_enrich is None:os.environ.pop("RDMVA_ENRICHR_MODE",None)
        else:os.environ["RDMVA_ENRICHR_MODE"]=old_enrich


def resume_post_attack(run_dir,*,pk_manifest=None,enable_online_ai=True):
    """Resume a V3 native run after the expensive universal screen/attack stages.

    Intended for transient filesystem/AV locks or other late non-scientific
    failures. It reuses persisted intermediates and does not repeat screening.
    """
    run=Path(run_dir)
    required=[
      "PATIENT_STATE_V30.json","MECHANISM_HYPOTHESES_V30.json",
      "RESCUE_OBJECTIVES_V30.json","CAUSAL_CANDIDATE_GRAPH_V30.json",
      "FALSIFICATION_RESULTS_V30.json","TRACK2_HYPOTHESES_V30.json",
      "source_bound_literature_attack/SOURCE_BOUND_ATTACK_SUMMARY.json",
      "master/entities.jsonl","master/summary.json","V30_ROUTE_DRUG_PAIRS.json",
    ]
    missing=[x for x in required if not (run/x).is_file()]
    if missing:raise RuntimeError("V30_RESUME_MISSING:"+",".join(missing))

    state=load(run/"PATIENT_STATE_V30.json")
    mechanisms=load(run/"MECHANISM_HYPOTHESES_V30.json")
    objectives=load(run/"RESCUE_OBJECTIVES_V30.json")
    graph=load(run/"CAUSAL_CANDIDATE_GRAPH_V30.json")
    falsified=load(run/"FALSIFICATION_RESULTS_V30.json")
    track2=load(run/"TRACK2_HYPOTHESES_V30.json")
    attack_summary=load(run/"source_bound_literature_attack/SOURCE_BOUND_ATTACK_SUMMARY.json")

    falsified=literature_attack.apply(falsified,attack_summary)
    _safe_save(run/"FALSIFICATION_RESULTS_V30.json",falsified)
    track2=challenge_hypotheses.apply_source_bound_evidence(track2,attack_summary)
    _safe_save(run/"TRACK2_HYPOTHESES_V30_FINAL.json",track2)
    _stage(run,"SOURCE_BOUND_ATTACK","completed",{
      "attacked":len(attack_summary.get("results",[])),
      "resumed_from_persisted_intermediates":True})

    ai_review=None
    if enable_online_ai and falsified.get("candidates"):
        try:
            ai_review=openai_reasoner.review(
                state,mechanisms,objectives,falsified,run/"AI_RESEARCH_CRITIC_V30.json")
        except Exception as exc:
            _safe_save(run/"AI_RESEARCH_CRITIC_ERROR_V30.json",{
              "error":type(exc).__name__,"message":str(exc)[:800],
              "deterministic_results_modified":False})

    entities=list(read_jsonl(run/"master/entities.jsonl"))
    master_summary=load(run/"master/summary.json")
    pairs=load(run/"V30_ROUTE_DRUG_PAIRS.json")
    pk_status,_=quant_pk.run(
        run/"specialist_pharmacology",manifest_path=pk_manifest,
        approved_drug_names=[e["canonical_label"] for e in entities])

    strict_dossier=dossier.build(state,mechanisms,objectives,falsified,ai_review)
    dossier.write(run/"research_dossier",strict_dossier)
    track2_dossier_meta=_write_track2_dossier(run,track2,state)
    result={
      "schema":SCHEMA,"version":VERSION,"completed_utc":now(),
      "patient_state":"PATIENT_STATE_V30.json",
      "mechanism_hypotheses":"MECHANISM_HYPOTHESES_V30.json",
      "rescue_objectives":"RESCUE_OBJECTIVES_V30.json",
      "native_routes":"NATIVE_ROUTES_V30.json",
      "all_approved_entities_screened":len(entities),
      "approved_product_records":master_summary.get("products"),
      "route_drug_pairs":len(pairs),
      "candidate_paths":len(graph.get("candidate_paths",[])),
      "strict_research_conclusion":strict_dossier.get("research_conclusion"),
      "validated_research_lead":strict_dossier.get("proposed_research_lead"),
      "track2_submission_hypotheses":"TRACK2_HYPOTHESES_V30_FINAL.json",
      "track2_single_agent_hypothesis_count":len(track2.get("submission_level_single_agents",[])),
      "track2_combination_hypothesis_count":len(track2.get("combination_hypotheses",[])),
      "track2_has_hypothesis":track2.get("has_submission_level_hypothesis"),
      "combination_search_max_components":3,
      "combination_synergy_assumed":False,
      "market_approved_medications_only_for_track2":True,
      "pk_models_source_qualified":pk_status.get("models_source_qualified",0),
      "openai_enabled":bool(enable_online_ai),
      "patient_raw_data_sent_to_ai":False,
      "clinical_recommendation":False,
      "clinical_efficacy_probability_computed":False,
      "dose_recommendation":False,
      "resumed_after_late_failure":True,
    }
    _safe_save(run/"FINAL_RESULT_V30_NATIVE.json",result)
    s=load(run/"STATUS_V30_NATIVE.json") if (run/"STATUS_V30_NATIVE.json").is_file() else {"stages":{}}
    s["completed"]=True;s["current_stage"]="COMPLETED";s["final_result"]="FINAL_RESULT_V30_NATIVE.json"
    _safe_save(run/"STATUS_V30_NATIVE.json",s)
    return result


def resume_from_screen_stage(run_dir,*,literature_attack_count=20,
                             enable_online_ai=True,
                             enable_public_safety_review=True,
                             enable_dgidb=True,enable_cross_source_bridge=True,
                             pk_manifest=None,revalidate_modifiers=True):
    """Resume V3.1 from persisted pre-drug discovery intermediates.

    Reuses the approved-drug master universe, Europe-PMC/OpenAI modifier
    bundles, AlphaFold context, network expansion and per-target ChEMBL caches.
    It never repeats patient intake or earlier AI discovery unless explicitly
    requested elsewhere.
    """
    run=Path(run_dir)
    required=[
      "case_private.json","PATIENT_STATE_V30.json","MECHANISM_HYPOTHESES_V30.json",
      "RESCUE_OBJECTIVES_V31.json","NATIVE_ROUTES_V31.json",
      "master/entities.jsonl","master/products.jsonl","master/summary.json"]
    missing=[x for x in required if not (run/x).is_file()]
    if missing:
        raise RuntimeError("V31_RESUME_SCREEN_MISSING:"+",".join(missing))

    private=load(run/"case_private.json")
    state=load(run/"PATIENT_STATE_V30.json")
    mechanisms=load(run/"MECHANISM_HYPOTHESES_V30.json")
    objectives=load(run/"RESCUE_OBJECTIVES_V31.json")
    route_bundle=load(run/"NATIVE_ROUTES_V31.json")

    modifier_routes=[]
    if revalidate_modifiers and (run/"modifier_discovery").is_dir():
        primary=sorted(private.get("genes") or {})
        findings=modifier_discovery.revalidate_existing(
            primary,run/"modifier_discovery",objectives)
        modifier_routes=modifier_discovery.to_routes(
            findings,run/"sources",objectives)
        if modifier_routes:
            route_bundle=native_routes.merge_objective_routes(
                route_bundle,modifier_routes,label="v31_modifier_revalidated")

    routes=route_bundle.get("routes") or []
    final_routes_path=run/"NATIVE_ROUTES_V31_FINAL.json"
    previous_bundle=load(final_routes_path) if final_routes_path.is_file() else None
    def _route_signature(bundle):
        rows=[]
        for r in (bundle or {}).get("routes",[]):
            rows.append((
              r.get("route_id"),r.get("case_gene"),str(r.get("target_concept") or "").upper(),
              r.get("desired_action"),r.get("variable"),
              tuple(sorted(str(x.get("accession") or "") for x in (r.get("targets") or []))),
              tuple(sorted(str(x) for x in (r.get("rescue_objective_ids") or [])))))
        return rows
    same_route_signature=bool(
      previous_bundle and _route_signature(previous_bundle)==_route_signature(route_bundle))
    _safe_save(final_routes_path,route_bundle)

    # Invalidate route-dependent caches only when the route signature changed.
    # Per-target act_CHEMBL*.json files are always preserved and reused.
    source_cache=run/"screen/sources"
    source_cache.mkdir(parents=True,exist_ok=True)
    if not same_route_signature:
        for name in ("target_index.json","route_activities.json"):
            (source_cache/name).unlink(missing_ok=True)

    entities=list(read_jsonl(run/"master/entities.jsonl"))
    products=list(read_jsonl(run/"master/products.jsonl"))
    master_summary=load(run/"master/summary.json")
    registry_index=registry_output.build_index(products)

    _stage(run,"V31_FAST_UNIVERSAL_SCREEN","running",{
      "routes":len(routes),"modifier_routes_revalidated":len(modifier_routes),
      "per_target_chembl_cache_preserved":True,
      "route_dependent_cache_reused":same_route_signature})
    ledger,screen_summary=fast_screen.universal_screen(
        entities,routes,run/"screen",max_workers=8)
    if len(ledger)!=len(entities):
        raise RuntimeError("V31_FAST_SCREEN_COVERAGE_FAILURE")

    grouped,group_summary=science.canonicalize_ledger(ledger)
    registry_index=science.merge_registry_by_family(grouped,registry_index)
    jsonl(run/"screen/scientific_parent_ledger.jsonl",grouped)
    _safe_save(run/"screen/SCIENTIFIC_GROUP_SUMMARY.json",group_summary)

    pairs=evidence_ranker.prioritize(grouped,routes)

    dgidb_payload={"hints":[],"hint_count":0}
    dgidb_audit={"added_pairs":0}
    if enable_dgidb:
        try:
            dgidb_payload=dgidb_discovery.discover(
                routes,grouped,run/"dgidb_discovery")
            pairs,dgidb_audit=dgidb_discovery.merge_pairs(
                pairs,dgidb_payload,routes,grouped)
        except Exception as exc:
            err=run/"dgidb_discovery/DGIDB_ERROR_V31.json"
            err.parent.mkdir(parents=True,exist_ok=True)
            _safe_save(err,{
              "error":type(exc).__name__,"message":str(exc)[:500],
              "aggregator_claims_are_submission_evidence":False})

    pairs,poly_summary=polypharmacology.build(
        grouped,pairs,run/"polypharmacology")

    # V3.1 direct source-bound approved-drug rescue: catches a medicine that
    # was explicitly observed to restore/stabilize a primary-gene product or
    # improve a primary-gene mutant phenotype in the drug-blind literature.
    direct_payload=direct_drug_rescue.discover(
        run,entities,sorted(private.get("genes") or {}),objectives,
        run/"direct_drug_rescue")
    if enable_online_ai and direct_payload.get("findings"):
        direct_payload=direct_drug_rescue.adjudicate(
            direct_payload,objectives,
            run/"direct_drug_rescue/DIRECT_DRUG_RESCUE_ADJUDICATED_V31.json",
            run/"direct_drug_rescue")
    else:
        direct_payload={**direct_payload,"accepted_findings":[],
                        "rejected_findings":direct_payload.get("findings",[]),
                        "adjudication_completed":False}
    pairs,direct_audit=direct_drug_rescue.merge_pairs(
        pairs,direct_payload,objectives)

    bridge_payload={"bridges":[],"bridge_count":0,"targets_scouted":0}
    bridge_audit={"added_pairs":0,"routes":[]}
    if enable_online_ai and enable_cross_source_bridge:
        bridge_payload=cross_source_bridge.discover_from_routes(
            routes,entities,objectives,run/"cross_source_bridge_strict_batch",
            max_targets=60,max_full_bridges=30)
        pairs,bridge_audit=cross_source_bridge.merge_pairs(
            pairs,bridge_payload)

    graph_routes=(list(routes)+list(direct_audit.get("routes") or [])+
                  list(bridge_audit.get("routes") or []))
    _safe_save(run/"V31_ROUTE_DRUG_PAIRS_FINAL.json",pairs)
    if direct_audit.get("routes") or bridge_audit.get("routes"):
        discovery_route_bundle=dict(route_bundle)
        discovery_route_bundle["routes"]=graph_routes
        discovery_route_bundle["direct_source_drug_rescue_route_count"]=len(
            direct_audit.get("routes") or [])
        discovery_route_bundle["cross_source_bridge_route_count"]=len(
            bridge_audit.get("routes") or [])
        _safe_save(run/"NATIVE_ROUTES_V31_WITH_DISCOVERY_BRIDGES.json",
                   discovery_route_bundle)
    _stage(run,"V31_FAST_UNIVERSAL_SCREEN","completed",{
      "all_entities_screened":len(ledger),
      "signal_entities":screen_summary.get("signal_entities"),
      "route_drug_pairs":len(pairs),
      "routes":len(routes),
      "modifier_routes_revalidated":len(modifier_routes),
      "dgidb_hints":dgidb_payload.get("hint_count",0),
      "dgidb_added_pairs":dgidb_audit.get("added_pairs",0),
      "polypharmacology_profiled_entities":poly_summary.get("profiled_entity_count",0),
      "direct_source_drug_rescue_raw_findings":direct_payload.get("finding_count",0),
      "direct_source_drug_rescue_accepted":direct_payload.get("accepted_count",0),
      "direct_source_drug_rescue_pairs":direct_audit.get("added_pairs",0),
      "cross_source_targets_scouted":bridge_payload.get("targets_scouted",0),
      "cross_source_bridges":bridge_payload.get("bridge_count",0),
      "cross_source_pairs_added":bridge_audit.get("added_pairs",0),
      "screen_status_counts":screen_summary.get("status_counts")})

    _stage(run,"V31_CAUSAL_GRAPH_AND_COMBINATIONS","running")
    graph=candidate_graph.build(objectives,graph_routes,pairs,[])
    _safe_save(run/"CAUSAL_CANDIDATE_GRAPH_V31.json",graph)
    falsified=falsification.evaluate(graph)
    _safe_save(run/"FALSIFICATION_RESULTS_V31.json",falsified)
    combos=combination_engine.discover(
        graph,objectives,max_components=3,shortlist=100,max_results=80)
    _safe_save(run/"COMBINATION_HYPOTHESES_RAW_V31.json",combos)
    track2=challenge_hypotheses.build(graph,objectives,combos)

    safety_reviewed=0
    raw_labels={}
    if enable_public_safety_review:
        names={x["candidate"] for x in track2.get("submission_level_single_agents",[])}
        for c in track2.get("combination_hypotheses",[]):
            names.update(c.get("components") or [])
        if names:
            safety_payload,raw_labels=safety_review.review_candidates(
                sorted(names),grouped,private,run/"safety_review_v31")
            track2=safety_review.annotate_hypotheses(
                track2,safety_payload,raw_labels)
            safety_reviewed=len(names)
    exposure_assessments=[]
    if enable_online_ai and raw_labels:
        track2,exposure_assessments=exposure_feasibility.assess_track2(
            track2,grouped,raw_labels,run/"exposure_feasibility_v31")
    _safe_save(run/"TRACK2_HYPOTHESES_V31_PREATTACK.json",track2)
    _stage(run,"V31_CAUSAL_GRAPH_AND_COMBINATIONS","completed",{
      "candidate_paths":len(graph.get("candidate_paths",[])),
      "single_submission_hypotheses":len(track2.get("submission_level_single_agents",[])),
      "combination_hypotheses":len(track2.get("combination_hypotheses",[])),
      "safety_candidates_reviewed":safety_reviewed})

    attack_summary=None
    if enable_online_ai and literature_attack_count>0 and falsified.get("candidates"):
        _stage(run,"V31_SOURCE_BOUND_ATTACK","running")
        priority=[]
        for x in track2.get("submission_level_single_agents",[]):
            if x.get("candidate") not in priority:
                priority.append(x.get("candidate"))
        for c in track2.get("combination_hypotheses",[]):
            for name in c.get("components") or []:
                if name not in priority:
                    priority.append(name)
        rank={name:i for i,name in enumerate(priority)}
        attack_input=dict(falsified)
        attack_input["candidates"]=sorted(
            list(falsified.get("candidates",[])),
            key=lambda x:(rank.get(x.get("candidate"),10**6),
                          -float(x.get("research_priority_score") or -999),
                          str(x.get("candidate") or "")))
        attack_summary=literature_attack.run(
            attack_input,graph,run/"source_bound_literature_attack_v31",
            count=int(literature_attack_count),max_docs=12)
        falsified=literature_attack.apply(falsified,attack_summary)
        _safe_save(run/"FALSIFICATION_RESULTS_V31_FINAL.json",falsified)
        track2=challenge_hypotheses.apply_source_bound_evidence(
            track2,attack_summary)
        _safe_save(run/"TRACK2_HYPOTHESES_V31_FINAL.json",track2)
        _stage(run,"V31_SOURCE_BOUND_ATTACK","completed",{
          "attacked":len(attack_summary.get("results",[])),
          "priority_candidates":len(priority)})
    else:
        _safe_save(run/"TRACK2_HYPOTHESES_V31_FINAL.json",track2)

    ai_review=None
    if enable_online_ai and falsified.get("candidates"):
        try:
            ai_review=openai_reasoner.review(
                state,mechanisms,objectives,falsified,
                run/"AI_RESEARCH_CRITIC_V31.json")
        except Exception as exc:
            _safe_save(run/"AI_RESEARCH_CRITIC_ERROR_V31.json",{
              "error":type(exc).__name__,"message":str(exc)[:800],
              "deterministic_results_modified":False})

    pk_status,_=quant_pk.run(
        run/"specialist_pharmacology_v31",manifest_path=pk_manifest,
        approved_drug_names=[e["canonical_label"] for e in entities])

    strict_dossier=dossier.build(
        state,mechanisms,objectives,falsified,ai_review)
    dossier.write(run/"research_dossier_v31",strict_dossier)
    track2_dossier_meta=_write_track2_dossier(run,track2,state)

    result={
      "schema":SCHEMA,"version":VERSION,"completed_utc":now(),
      "resumed_from_screen_stage":True,
      "patient_state":"PATIENT_STATE_V30.json",
      "disease_state_graph":"DISEASE_STATE_GRAPH_V31.json",
      "rescue_objectives":"RESCUE_OBJECTIVES_V31.json",
      "native_routes":("NATIVE_ROUTES_V31_WITH_DIRECT_DRUG_RESCUE.json" if (run/"NATIVE_ROUTES_V31_WITH_DIRECT_DRUG_RESCUE.json").is_file() else "NATIVE_ROUTES_V31_FINAL.json"),
      "route_count":len(graph_routes),
      "revalidated_modifier_route_count":len(modifier_routes),
      "direct_source_drug_rescue_raw_findings":direct_payload.get("finding_count",0),
      "direct_source_drug_rescue_accepted":direct_payload.get("accepted_count",0),
      "direct_source_drug_rescue_pair_count":direct_audit.get("added_pairs",0),
      "all_approved_entities_screened":len(ledger),
      "approved_product_records":master_summary.get("products") or len(products),
      "route_drug_pairs":len(pairs),
      "candidate_paths":len(graph.get("candidate_paths",[])),
      "polypharmacology_profiled_entities":poly_summary.get("profiled_entity_count",0),
      "dgidb_discovery_hint_count":dgidb_payload.get("hint_count",0),
      "strict_research_conclusion":strict_dossier.get("research_conclusion"),
      "validated_research_lead":strict_dossier.get("proposed_research_lead"),
      "track2_submission_hypotheses":"TRACK2_HYPOTHESES_V31_FINAL.json",
      "track2_single_agent_hypothesis_count":len(track2.get("submission_level_single_agents",[])),
      "track2_combination_hypothesis_count":len(track2.get("combination_hypotheses",[])),
      "track2_has_hypothesis":bool(track2.get("has_submission_level_hypothesis")),
      "combination_search_max_components":3,
      "combination_synergy_assumed":False,
      "alphamissense_default_profile_used":False,
      "alphafold_context_used":(run/"variant_structure/VARIANT_STRUCTURE_CONTEXT_V31.json").is_file(),
      "market_approved_medications_only_for_track2":True,
      "pk_models_source_qualified":pk_status.get("models_source_qualified",0),
      "openai_enabled":bool(enable_online_ai),
      "patient_raw_data_sent_to_ai":False,
      "clinical_recommendation":False,
      "clinical_efficacy_probability_computed":False,
      "dose_recommendation":False}
    _safe_save(run/"FINAL_RESULT_V31_NATIVE.json",result)

    s=load(run/"STATUS_V30_NATIVE.json") if (run/"STATUS_V30_NATIVE.json").is_file() else {"stages":{}}
    s["completed"]=True
    s["current_stage"]="COMPLETED"
    s["final_result"]="FINAL_RESULT_V31_NATIVE.json"
    s["version"]=VERSION
    _safe_save(run/"STATUS_V30_NATIVE.json",s)
    return result


def resume_with_direct_drug_rescue(run_dir,*,literature_attack_count=12,
                                   enable_online_ai=True,
                                   enable_public_safety_review=True,
                                   enable_cross_source_bridge=True,
                                   pk_manifest=None):
    """Fast V3.1 resume after universal screening.

    Reuses persisted 7,523-entity screening and adds only approved compounds
    explicitly observed in the already-retrieved drug-blind literature to
    rescue/stabilize a primary-gene product or mutant-model phenotype.
    """
    run=Path(run_dir)
    required=[
      "case_private.json","PATIENT_STATE_V30.json","MECHANISM_HYPOTHESES_V30.json",
      "RESCUE_OBJECTIVES_V31.json","NATIVE_ROUTES_V31_FINAL.json",
      "V31_ROUTE_DRUG_PAIRS_FINAL.json","screen/scientific_parent_ledger.jsonl",
      "master/entities.jsonl","master/products.jsonl","master/summary.json"]
    missing=[x for x in required if not (run/x).is_file()]
    if missing:
        raise RuntimeError("V31_DIRECT_RESUME_MISSING:"+",".join(missing))

    private=load(run/"case_private.json")
    state=load(run/"PATIENT_STATE_V30.json")
    mechanisms=load(run/"MECHANISM_HYPOTHESES_V30.json")
    objectives=load(run/"RESCUE_OBJECTIVES_V31.json")
    route_bundle=load(run/"NATIVE_ROUTES_V31_FINAL.json")
    routes=list(route_bundle.get("routes") or [])
    pairs=load(run/"V31_ROUTE_DRUG_PAIRS_FINAL.json")
    entities=list(read_jsonl(run/"master/entities.jsonl"))
    products=list(read_jsonl(run/"master/products.jsonl"))
    grouped=list(read_jsonl(run/"screen/scientific_parent_ledger.jsonl"))
    master_summary=load(run/"master/summary.json")

    _stage(run,"V31_DIRECT_SOURCE_DRUG_RESCUE","running")
    direct_payload=direct_drug_rescue.discover(
        run,entities,sorted(private.get("genes") or {}),objectives,
        run/"direct_drug_rescue_final")
    if enable_online_ai and direct_payload.get("findings"):
        direct_payload=direct_drug_rescue.adjudicate(
            direct_payload,objectives,
            run/"direct_drug_rescue_final/DIRECT_DRUG_RESCUE_ADJUDICATED_V31.json",
            run/"direct_drug_rescue_final")
    else:
        direct_payload={**direct_payload,"accepted_findings":[],
                        "rejected_findings":direct_payload.get("findings",[]),
                        "accepted_count":0,"adjudication_completed":False}
    pairs2,direct_audit=direct_drug_rescue.merge_pairs(
        pairs,direct_payload,objectives)

    bridge_payload={"bridges":[],"bridge_count":0,"targets_scouted":0}
    bridge_audit={"added_pairs":0,"routes":[]}
    if enable_online_ai and enable_cross_source_bridge:
        bridge_payload=cross_source_bridge.discover_from_routes(
            routes,entities,objectives,run/"cross_source_bridge_strict_batch",
            max_targets=60,max_full_bridges=30)
        pairs2,bridge_audit=cross_source_bridge.merge_pairs(
            pairs2,bridge_payload)

    graph_routes=(routes+list(direct_audit.get("routes") or [])+
                  list(bridge_audit.get("routes") or []))
    discovery_bundle=dict(route_bundle)
    discovery_bundle["routes"]=graph_routes
    discovery_bundle["route_count"]=len(graph_routes)
    discovery_bundle["direct_source_drug_rescue_route_count"]=len(
        direct_audit.get("routes") or [])
    discovery_bundle["cross_source_bridge_route_count"]=len(
        bridge_audit.get("routes") or [])
    _safe_save(run/"NATIVE_ROUTES_V31_WITH_DISCOVERY_BRIDGES.json",discovery_bundle)
    _safe_save(run/"V31_ROUTE_DRUG_PAIRS_WITH_DISCOVERY_BRIDGES.json",pairs2)
    _stage(run,"V31_DIRECT_SOURCE_DRUG_RESCUE","completed",{
      "raw_findings":direct_payload.get("finding_count",0),
      "accepted_hypotheses":direct_payload.get("accepted_count",0),
      "added_direct_pairs":direct_audit.get("added_pairs",0),
      "cross_source_targets_scouted":bridge_payload.get("targets_scouted",0),
      "cross_source_bridges":bridge_payload.get("bridge_count",0),
      "cross_source_pairs_added":bridge_audit.get("added_pairs",0),
      "universal_screen_repeated":False})

    _stage(run,"V31_DIRECT_RESCUE_GRAPH","running")
    graph=candidate_graph.build(objectives,graph_routes,pairs2,[])
    _safe_save(run/"CAUSAL_CANDIDATE_GRAPH_V31_DIRECT.json",graph)
    falsified=falsification.evaluate(graph)
    _safe_save(run/"FALSIFICATION_RESULTS_V31_DIRECT.json",falsified)
    combos=combination_engine.discover(
        graph,objectives,max_components=3,shortlist=150,max_results=100)
    _safe_save(run/"COMBINATION_HYPOTHESES_V31_DIRECT.json",combos)
    track2=challenge_hypotheses.build(graph,objectives,combos)

    safety_reviewed=0
    raw_labels={}
    if enable_public_safety_review:
        names={x["candidate"] for x in track2.get("submission_level_single_agents",[])}
        for c in track2.get("combination_hypotheses",[]):
            names.update(c.get("components") or [])
        if names:
            safety_payload,raw_labels=safety_review.review_candidates(
                sorted(names),grouped,private,run/"safety_review_v31_direct")
            track2=safety_review.annotate_hypotheses(
                track2,safety_payload,raw_labels)
            safety_reviewed=len(names)
    exposure_assessments=[]
    if enable_online_ai and raw_labels:
        track2,exposure_assessments=exposure_feasibility.assess_track2(
            track2,grouped,raw_labels,run/"exposure_feasibility_v31_direct")
    _safe_save(run/"TRACK2_HYPOTHESES_V31_DIRECT_PREATTACK.json",track2)
    _stage(run,"V31_DIRECT_RESCUE_GRAPH","completed",{
      "candidate_paths":len(graph.get("candidate_paths",[])),
      "single_submission_hypotheses":len(track2.get("submission_level_single_agents",[])),
      "combination_hypotheses":len(track2.get("combination_hypotheses",[])),
      "safety_candidates_reviewed":safety_reviewed})

    attack_summary=None
    if enable_online_ai and literature_attack_count>0 and falsified.get("candidates"):
        priority=[]
        for x in track2.get("submission_level_single_agents",[]):
            if x.get("candidate") not in priority:priority.append(x.get("candidate"))
        for c in track2.get("combination_hypotheses",[]):
            for name in c.get("components") or []:
                if name not in priority:priority.append(name)
        rank={name:i for i,name in enumerate(priority)}
        attack_input=dict(falsified)
        attack_input["candidates"]=sorted(
            list(falsified.get("candidates",[])),
            key=lambda x:(rank.get(x.get("candidate"),10**6),
                          -float(x.get("research_priority_score") or -999),
                          str(x.get("candidate") or "")))
        _stage(run,"V31_DIRECT_RESCUE_COUNTEREVIDENCE","running",{
          "priority_candidates":priority[:20]})
        attack_summary=literature_attack.run(
            attack_input,graph,run/"source_bound_literature_attack_v31_direct",
            count=int(literature_attack_count),max_docs=14)
        falsified=literature_attack.apply(falsified,attack_summary)
        _safe_save(run/"FALSIFICATION_RESULTS_V31_DIRECT_FINAL.json",falsified)
        track2=challenge_hypotheses.apply_source_bound_evidence(
            track2,attack_summary)
        _stage(run,"V31_DIRECT_RESCUE_COUNTEREVIDENCE","completed",{
          "attacked":len(attack_summary.get("results",[]))})
    _safe_save(run/"TRACK2_HYPOTHESES_V31_DIRECT_FINAL.json",track2)

    ai_review=None
    if enable_online_ai and falsified.get("candidates"):
        try:
            ai_review=openai_reasoner.review(
                state,mechanisms,objectives,falsified,
                run/"AI_RESEARCH_CRITIC_V31_DIRECT.json")
        except Exception as exc:
            _safe_save(run/"AI_RESEARCH_CRITIC_ERROR_V31_DIRECT.json",{
              "error":type(exc).__name__,"message":str(exc)[:800],
              "deterministic_results_modified":False})

    pk_status,_=quant_pk.run(
        run/"specialist_pharmacology_v31_direct",manifest_path=pk_manifest,
        approved_drug_names=[e["canonical_label"] for e in entities])
    strict_dossier=dossier.build(state,mechanisms,objectives,falsified,ai_review)
    dossier.write(run/"research_dossier_v31_direct",strict_dossier)
    track2_dossier_meta=_write_track2_dossier(run,track2,state)

    result={
      "schema":SCHEMA,"version":VERSION,"completed_utc":now(),
      "track2_research_dossier":track2_dossier_meta["json"],
      "track2_research_dossier_markdown":track2_dossier_meta["markdown"],
      "primary_submission_candidate":track2_dossier_meta["primary_candidate"],
      "primary_submission_evidence_class":track2_dossier_meta["primary_evidence_class"],
      "resume_mode":"post_screen_direct_source_drug_rescue",
      "universal_screen_repeated":False,
      "all_approved_entities_screened":len(entities),
      "approved_product_records":master_summary.get("products") or len(products),
      "base_route_drug_pairs":len(pairs),
      "route_drug_pairs_with_direct_rescue":len(pairs2),
      "route_count":len(graph_routes),
      "candidate_paths":len(graph.get("candidate_paths",[])),
      "direct_source_drug_rescue_raw_findings":direct_payload.get("finding_count",0),
      "direct_source_drug_rescue_accepted":direct_payload.get("accepted_count",0),
      "direct_source_drug_rescue_pair_count":direct_audit.get("added_pairs",0),
      "strict_research_conclusion":strict_dossier.get("research_conclusion"),
      "validated_research_lead":strict_dossier.get("proposed_research_lead"),
      "track2_submission_hypotheses":"TRACK2_HYPOTHESES_V31_DIRECT_FINAL.json",
      "track2_single_agent_hypothesis_count":len(track2.get("submission_level_single_agents",[])),
      "track2_combination_hypothesis_count":len(track2.get("combination_hypotheses",[])),
      "track2_has_hypothesis":bool(track2.get("has_submission_level_hypothesis")),
      "alphafold_context_used":(run/"variant_structure/VARIANT_STRUCTURE_CONTEXT_V31.json").is_file(),
      "alphamissense_default_profile_used":False,
      "market_approved_medications_only_for_track2":True,
      "pk_models_source_qualified":pk_status.get("models_source_qualified",0),
      "openai_enabled":bool(enable_online_ai),
      "patient_raw_data_sent_to_ai":False,
      "clinical_recommendation":False,
      "clinical_efficacy_probability_computed":False,
      "dose_recommendation":False}
    _safe_save(run/"FINAL_RESULT_V31_DIRECT_RESCUE.json",result)
    s=load(run/"STATUS_V30_NATIVE.json") if (run/"STATUS_V30_NATIVE.json").is_file() else {"stages":{}}
    s["completed"]=True
    s["current_stage"]="COMPLETED"
    s["final_result"]="FINAL_RESULT_V31_DIRECT_RESCUE.json"
    s["version"]=VERSION
    _safe_save(run/"STATUS_V30_NATIVE.json",s)
    return result


def _write_track2_dossier(run,track2,patient_state):
    try:
        dsg=load(Path(run)/"DISEASE_STATE_GRAPH_V31.json")
    except Exception:
        dsg={"nodes":[]}
    try:
        policy=load(ROOT/"config/v30_source_policy.json")
    except Exception:
        policy={}
    d=challenge_dossier.build(track2,patient_state,dsg,policy)
    out=Path(run)/"challenge_submission_v31"
    challenge_dossier.write(out,d)
    return {
      "json":"challenge_submission_v31/TRACK2_RESEARCH_DOSSIER_V31.json",
      "markdown":"challenge_submission_v31/TRACK2_RESEARCH_DOSSIER_V31.md",
      "primary_candidate":(d.get("primary_hypothesis") or {}).get("candidate"),
      "primary_evidence_class":(d.get("primary_hypothesis") or {}).get("evidence_class")}
