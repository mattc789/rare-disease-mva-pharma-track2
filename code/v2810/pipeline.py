"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
from collections import Counter
import argparse,json,traceback
from .core import *
from . import intake,master,routes,ai,screen,deep,registry_output,science,case_science,evidence,reporting,predictive,expert_genomics,expert_triage,expert_protein,vep_local,gnomad_local,quant_pk,specialist_evidence,pksim_bridge,resource_bootstrap,open_splicing

DEFAULT_VCF=None
DEFAULT_CLINICAL=None

def stage(run,name,state,detail=None):
 p=Path(run)/'STATUS.json'
 s=load(p) if p.exists() else {'version':VERSION,'stages':{},'completed':False}
 s['stages'][name]={'state':state,'utc':now(),'detail':detail};s['current_stage']=name;save(p,s)

def run(vcf_path=DEFAULT_VCF,clinical_path=DEFAULT_CLINICAL,case_path=None,outdir=None,final_count=5,use_ai=False,functional_review_count=12,population_af_path=None,safety_review_count=20,gnomad_manifest=None,spliceai_manifest=None,pk_manifest=None,vep_on_case_replay=False,auto_resources=False,resource_download=True,resource_max_gib=240):
 raw_intake_requested=not bool(case_path)
 run=Path(outdir) if outdir else RUN_ROOT/('E2E_OPEN_RESEARCH_V2810_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
 run.mkdir(parents=True,exist_ok=True)
 try:
  if auto_resources:
   stage(run,'PUBLIC_PINNED_BASELINE_BOOTSTRAP','running')
   base_status=resource_bootstrap.ensure_reproducible_public_baseline(
      ROOT,run/'resource_bootstrap/baseline',install=resource_download)
   stage(run,'PUBLIC_PINNED_BASELINE_BOOTSTRAP',
      'completed' if base_status.get('state')=='READY' else 'failed',
      {'state':base_status.get('state'),
       'gated_private_genomes_downloaded':False,
       'scientific_source_snapshot_pinned':True})
   if base_status.get('state')!='READY':
      raise RuntimeError('V2810_PUBLIC_SOURCE_BASELINE_UNAVAILABLE_OR_DRIFTED')
  if case_path:
   case_path=Path(case_path)
   stage(run,'RAW_INTAKE','skipped',{'reason':'explicit_prebuilt_case_supplied','case_path':str(case_path)})
  else:
   stage(run,'RAW_INTAKE','running')
   raw_case,im=intake.build_case_model(Path(vcf_path),Path(clinical_path) if clinical_path else None,run/'intake')
   case_path=run/'intake/case_model_v2.json'
   stage(run,'RAW_INTAKE','completed',{'case_model':str(case_path),'primary_genes':raw_case.get('primary_genes',[]),
     'genome_wide_candidates':len(raw_case.get('genome_wide_candidate_genes',[])),'patient_data_external':False})
  stage(run,'CASE_MODEL','running')
  private,public=normalize_case(case_path);save(run/'case_private.json',private);save(run/'case_public_safe.json',public)
  if not private['genes']:raise RuntimeError('V28_INTAKE_NO_PRIMARY_GENE_HYPOTHESIS')
  stage(run,'CASE_MODEL','completed',{'genes':sorted(private['genes']),'variant_classes':private['variant_classes'],
    'diagnostic_certainty':private.get('diagnostic_certainty'),'patient_api_payload':False})

  # License-safe public source downloads are chromosome-aware, never gene/disease hardcoded.
  # This is enabled by default in the CLI; direct programmatic API callers opt in.
  resource_bootstrap.write_license_audit(run/'license_and_redistribution')
  bootstrap_report={'schema':resource_bootstrap.VERSION,'mode':'DISABLED_BY_PROGRAMMATIC_CALL'}
  vep_setup={'state':'UNASSESSED','vep_config':{}}
  local_gnomad_manifest=gnomad_manifest
  if auto_resources:
   stage(run,'PUBLIC_RESOURCE_BOOTSTRAP','running')
   try:
    vep_setup=resource_bootstrap.ensure_vep_prerequisites(
       ROOT,run/'resource_bootstrap/vep',install=resource_download)
    if local_gnomad_manifest:
     bootstrap_report={'schema':resource_bootstrap.VERSION,
       'mode':'EXPLICIT_LOCAL_GNOMAD_MANIFEST_SUPPLIED',
       'completion':'USER_PINNED_REFERENCE'}
    else:
     local_gnomad_manifest,bootstrap_report=resource_bootstrap.prepare(
       private,run/'resource_bootstrap/gnomad',
       download=resource_download,max_gib=resource_max_gib)
    stage(run,'PUBLIC_RESOURCE_BOOTSTRAP','completed',{
      'vep':vep_setup.get('state'),'gnomad':bootstrap_report.get('completion'),
      'public_chromosomes':resource_bootstrap.needed_chromosomes(private)[0],
      'patient_genetic_coordinates_sent_to_remote_hosts':False,
      'downloaded_data_in_redistributable_code':False})
   except Exception as exc:
    bootstrap_report={'schema':resource_bootstrap.VERSION,
      'completion':'FAILED_NONFATAL','error':type(exc).__name__+':'+str(exc)[-240:],
      'no_private_genome_uploaded':True}
    stage(run,'PUBLIC_RESOURCE_BOOTSTRAP','completed',bootstrap_report)
  stage(run,'PATIENT_VARIANT_MECHANISM_MODEL','running')
  patient_mechanisms=case_science.build_case_science(private)
  save(run/'PATIENT_VARIANT_MECHANISM_MODEL.json',patient_mechanisms)
  stage(run,'PATIENT_VARIANT_MECHANISM_MODEL','completed',{
    'primary_gene_hypotheses':len(patient_mechanisms['primary_gene_hypotheses']),
    'variant_specific_function_requires_validation':True,
    'allelic_phase':private.get('phase')})
  stage(run,'EXPERT_VARIANT_EVIDENCE_AUDIT','running')
  genomic_expert=expert_genomics.build(private,population_af_path)
  protein_context=expert_protein.build_from_private(private,genomic_expert,RAW)
  expert_protein.save_report(protein_context,run/'expert_genetics')
  for gene,rec in genomic_expert['genes'].items():
   rec['uniprot_protein_context']=(protein_context['genes'].get(gene) or {}).get('source')
   rec['protein_region_annotations']=(protein_context['genes'].get(gene) or {}).get('annotated_protein_positions',[])
  stage(run,'SPECIALIST_VEP_SPLICEAI_GNOMAD','running')
  if raw_intake_requested or vep_on_case_replay:
   vep_status,vep_records=vep_local.annotation_stage(
     private,run/'specialist_genomics/vep',
     config=vep_setup.get('vep_config') if vep_setup.get('state')=='READY' else None)
  else:
   vep_status={'schema':vep_local.SCHEMA,
     'tool_state':'SKIPPED_EXPLICIT_PREBUILT_CASE_NO_VEP_ON_REPLAY',
     'no_patient_variants_sent_external':True}
   vep_records=[]
  # Illumina SpliceAI PolyForm-Strict/NC weights are excluded from a
  # free, redistributable hackathon build. A GPL-3 Pangolin adapter is separate.
  if spliceai_manifest:raise ValueError('V2810_ILLUMINA_SPLICEAI_DISALLOWED_BY_OPEN_RESEARCH_LICENSE_PROFILE')
  spliceai_status={'state':'EXCLUDED_RESTRICTED_UPSTREAM_LICENSE',
    'commercial_use_automatically_permitted':False,
    'source':'Illumina models not downloaded or included in CC-BY submission'}
  splice_records=[]
  gnomad_status,gnomad_records=gnomad_local.run(private,
      run/'specialist_genomics/gnomad',manifest_path=local_gnomad_manifest)
  if auto_resources:
   stage(run,'FREE_OPEN_SOURCE_SPLICE_PREDICTION','running')
   pangolin_status=open_splicing.score_case(private,
     run/'specialist_genomics/pangolin_gpl3',
     install=resource_download,run_if_ready=resource_download)
   stage(run,'FREE_OPEN_SOURCE_SPLICE_PREDICTION','completed',{
     'state':pangolin_status.get('prediction_state'),
     'license':'GPL-3.0 independent software',
     'illumina_models_used':False,
     'patient_genotypes_sent_external':False})
  else:
   pangolin_status={'state':'SOURCE_CHECK_DISABLED_BY_PROGRAMMATIC_API',
     'prediction_state':'SKIPPED','license':'GPL-3.0 standalone'}
  pangolin_rows=[]
  pangolin_file=run/'specialist_genomics/pangolin_gpl3/PANGOLIN_PRIVATE_HASHED_SCORES.jsonl'
  if pangolin_file.is_file():
   with pangolin_file.open(encoding='utf-8') as f:
    pangolin_rows=[json.loads(line) for line in f if line.strip()]
  integrated=specialist_evidence.attach(genomic_expert,private,
      vep_records,splice_records,gnomad_records,pangolin_records=pangolin_rows)
  expert_genomics.report(genomic_expert,run/'expert_genetics')
  stage(run,'SPECIALIST_VEP_SPLICEAI_GNOMAD','completed',{
    'vep_state':vep_status['tool_state'],
    'spliceai_state':spliceai_status['state'],
    'gnomad_state':gnomad_status['state'],
    'source_verified_genomic_variants':integrated.get('gnomad_variants_matched',0),
    'patient_variants_sent_external':False})
  stage(run,'EXPERT_VARIANT_EVIDENCE_AUDIT','completed',{
    'genes':len(genomic_expert['genes']),
    'population_af_status':genomic_expert['population_frequency_source']['state'],
    'no_full_acmg_claim':True,
    'clinical_diagnosis_confirmed':False})
  stage(run,'MASTER_UNIVERSE','running')
  entities,products,ms=master.build(run/'master')
  if not entities or len(entities)!=ms.get('entities',len(entities)):
   raise RuntimeError('V283_GLOBAL_UNIVERSE_COVERAGE_MISMATCH')
  registry_index=registry_output.build_index(products)
  stage(run,'MASTER_UNIVERSE','completed',{**ms,'registry_snapshot':registry_output.SNAPSHOT,'indexed_ingredients':len(registry_index)})

  stage(run,'QUANTITATIVE_PK_QUALIFICATION','running')
  approved_names=[e['canonical_label'] for e in entities]
  quant_pk_status,quant_pk_models=quant_pk.run(
      run/'specialist_pharmacology',manifest_path=pk_manifest,
      approved_drug_names=approved_names)
  control_path=ROOT/'audits/PKSIM_V289_MIDAZOLAM_REFERENCE_CONTROL'
  if control_path.is_dir():
   pksim_control=pksim_bridge.audit_control(control_path,
      run/'specialist_pharmacology/PKSIM_PUBLIC_REFERENCE_CONTROL.json')
   quant_pk_status['osp_midazolam_public_engineering_control']=pksim_control
   quant_pk_status['midazolam_control_is_not_a_case_research_candidate']=True
  specialist_sources={'vep':vep_status,'spliceai':spliceai_status,'gnomad':gnomad_status,'pangolin':pangolin_status}
  specialist_report=specialist_evidence.make_report(
      run/'specialist_reports',specialist_sources,quant_pk_status,integrated)
  stage(run,'QUANTITATIVE_PK_QUALIFICATION','completed',{
    'qualified_populations':quant_pk_status.get('models_source_qualified',0),
    'no_patient_efficacy_estimated':True,
    'report':'specialist_reports/SPECIALIST_EVIDENCE_STATUS.md'})
  stage(run,'ROUTE_GENERATION','running')
  base=routes.base_routes(private,run/'sources')
  discovered=[]
  if use_ai:
   for gene in sorted(private['genes']):
    aliases=routes.search_aliases(gene,run/'sources')
    docs=ai.literature(gene,run/'mechanism_discovery'/gene,max_docs=15,targets=aliases[1:])
    save(run/'mechanism_discovery'/gene/'search_aliases.json',aliases)
    try:
     resolver=lambda g:routes.accession_for_gene(g,run/'sources')
     dr1=ai.discover_routes(gene,docs,run/'mechanism_discovery'/gene,resolver,tag='route_discovery_wave1')
     exclude=[x.get('target_gene') for x in dr1 if x.get('target_gene')]
     dr2=ai.discover_routes(gene,docs,run/'mechanism_discovery'/gene,resolver,exclude_targets=exclude,tag='route_discovery_wave2')
     seen=set();dr=[]
     for x in dr1+dr2:
      k=(x.get('target_accession'),x.get('desired_action'),x.get('variable'),x.get('source_id'))
      if k in seen:continue
      seen.add(k);dr.append(x)
     dr=ai.verify_routes(gene,dr,run/'mechanism_discovery'/gene)
     for x in dr:
      x['case_gene']=gene
      try:
       cg=routes.canonical_gene_for_name(x.get('target_gene'),run/'sources')
       x['target_gene']=cg
       if x.get('targets'):x['targets'][0]['gene']=cg
      except Exception:pass
     discovered+=dr
    except Exception as exc:
     save(run/'mechanism_discovery'/gene/'AI_ERROR.json',{'error':type(exc).__name__,'message':str(exc)[:300]})
     raise RuntimeError('V28_ROUTE_DISCOVERY_REQUIRED_FAILED:'+gene+':'+type(exc).__name__) from exc
  merged=routes.merge_discovered(base,discovered)
  merged=routes.apply_case_gates(merged,private)
  merged,route_gate_audit=science.adjudicate_routes(merged)
  for r in merged:
   r['allele_specific_compatibility']='unknown_requires_variant_functional_evidence'
   r['patient_variant_mechanism_model_ref']='PATIENT_VARIANT_MECHANISM_MODEL.json'
  save(run/'SCIENTIFIC_ROUTE_GATES.json',route_gate_audit)
  save(run/'routes_base.json',base);save(run/'routes_discovered.json',discovered);save(run/'routes_final.json',merged)
  stage(run,'ROUTE_GENERATION','completed',{'base':len(base),'discovered':len(discovered),'final':len(merged),'ai_enabled':use_ai})

  stage(run,'UNIVERSAL_SCREEN','running')
  ledger,ss=screen.universal_screen(entities,merged,run/'screen')
  if len(ledger)!=len(entities):raise RuntimeError('V28_NOT_ALL_ENTITIES_SCREENED')
  grouped,group_summary=science.canonicalize_ledger(ledger)
  txcache=run/'screen/sources/transcriptomics.json'
  txrows=load(txcache).get('rows',[]) if txcache.is_file() else []
  science.attach_transcript_context(grouped,txrows)
  registry_index=science.merge_registry_by_family(grouped,registry_index)
  jsonl(run/'screen/scientific_parent_ledger.jsonl',grouped)
  save(run/'screen/SCIENTIFIC_GROUP_SUMMARY.json',group_summary)
  stage(run,'UNIVERSAL_SCREEN','completed',{**ss,'scientific_groups':len(grouped),
    'parent_signal_groups':group_summary['signal_groups'],'source_identity_warnings':True})

  stage(run,'INITIAL_EVIDENCE_DISCOVERY','running')
  initial_pairs=deep.prioritize(grouped,merged)
  # Locally score patient-allele-compatible hypotheses BEFORE AI literature search.
  # Require genuine single-ingredient registry evidence; no raw variants go to AI.
  v286_pk_reference=Path(__file__).resolve().parents[2]/'config'/'curated_human_pk_v286.json'
  initial_pairs,patient_first_early_audit=predictive.rank_pairs(
     initial_pairs,patient_mechanisms,registry_index,private,
     curated_pk_path=v286_pk_reference,exclude_unregistered=True)
  initial_pairs,expert_early_audit=expert_triage.order_pairs(initial_pairs,genomic_expert)
  save(run/'EXPERT_DISCOVERY_TRIAGE.json',expert_early_audit)
  save(run/'PATIENT_FIRST_DISCOVERY_PRIORITIES.json',patient_first_early_audit)
  save(run/'INITIAL_EVIDENCE_PAIRS.json',initial_pairs)
  stage(run,'INITIAL_EVIDENCE_DISCOVERY','completed',{'pairs':len(initial_pairs),
    'evidence_categories':dict(Counter(p['tier'] for p in initial_pairs))})
  stage(run,'EXPLORATORY_FUNCTIONAL_SEARCH','running')
  exploration=science.exploratory_functional_evidence_search(
     initial_pairs,next(iter(sorted(private['genes'])),''),
     run/'exploratory_functional_search',registry_index=registry_index,
     count=functional_review_count,use_ai=use_ai)
  stage(run,'EXPLORATORY_FUNCTIONAL_SEARCH','completed',
     {'candidate_groups_reviewed':len(exploration),
      'states':dict(Counter(x['evidence_state'] for x in exploration)),
      'automatic_clinical_promotion_allowed':False})
  stage(run,'EVIDENCE_RECLASSIFICATION','running')
  pairs,upgrades=evidence.enrich_pairs_from_literature(initial_pairs,exploration,private['genes'])
  # Recalculate after source verification without replacing evidence A/B/C/D grades.
  pairs,patient_first_post_evidence=predictive.rank_pairs(
     pairs,patient_mechanisms,registry_index,private,
     curated_pk_path=v286_pk_reference,exclude_unregistered=True)
  pairs,expert_post_evidence=expert_triage.order_pairs(pairs,genomic_expert)
  save(run/'EXPERT_POST_EVIDENCE_TRIAGE.json',expert_post_evidence)
  save(run/'PATIENT_FIRST_POST_EVIDENCE_PRIORITIES.json',patient_first_post_evidence)
  save(run/'priority_pairs.json',pairs)
  save(run/'RESEARCH_EVIDENCE_UPGRADES.json',upgrades)
  stage(run,'EVIDENCE_RECLASSIFICATION','completed',{'research_evidence_upgrades':len(upgrades),
     'evidence_categories':dict(Counter(p['tier'] for p in pairs))})
  stage(run,'PRIORITY_AND_GATES','running')
  gated=deep.deterministic_gates(pairs,grouped,private,run/'deep_gate/labels')
  save(run/'deep_gate/AB_gate_results.json',gated)
  finalists=deep.choose_finalists(gated,final_count)
  for f in finalists:f['regulatory_approval']=registry_output.for_candidate(f['ingredient'],registry_index)
  save(run/'finalists.json',finalists)
  stage(run,'PRIORITY_AND_GATES','completed',{'hypothesis_pairs':len(pairs),
    'A_B_research_candidates':len(gated),'finalists':len(finalists),
    'new_finalists':[x['ingredient'] for x in finalists]})
  stage(run,'FINALIST_DEEP_ATTACK','running')
  attack=deep.attack_finalists(finalists,sorted(private['genes']),run/'finalist_attack',
    registry_index=registry_index) if use_ai else []
  stage(run,'FINALIST_DEEP_ATTACK','completed',{'attacked':len(attack),
    'verdicts':dict(Counter(x['verdict'] for x in attack))})
  stage(run,'RESEARCH_SHORTLIST','running')
  research_shortlist,shortlist_summary=deep.build_research_shortlist(
    pairs,gated,finalists,attack,run/'research_shortlist',count=10,
    use_ai=use_ai,entities=grouped,registry_index=registry_index)
  stage(run,'RESEARCH_SHORTLIST','completed',shortlist_summary)
  stage(run,'BASIC_SAFETY_AND_EVIDENCE_REPORTS','running')
  safety={}
  public_labels_retrieved=0
  lookup={e['canonical_label']:e for e in grouped}
  prefmap=deep.chembl_pref_names(grouped) if use_ai else {}
  for row in research_shortlist:
   name=row['candidate']
   if use_ai:
    try:
     safety[name]=deep.safety_gate(name,lookup[name],private,prefmap,run/'safety_label_sources')
     public_labels_retrieved+=1
    except Exception as exc:
     safety[name]={'state':'HOLD','flags':['label_lookup_unavailable'],
       'reason':type(exc).__name__,'label_query':name}
   else:
    safety[name]={'state':'UNASSESSED','flags':['no_label_lookup_in_offline_mode'],
       'label_query':name}
   row['safety_and_suitability']=reporting.baseline_safety(private,safety[name])
   pair=next((x for x in pairs if x['ingredient']==name and x['route_id']==row['route_id']),None)
   row['grade_explanation']=evidence.evidence_explanation(pair or row)
   row['clinical_recommendation']=False
  save(run/'research_shortlist/RESEARCH_SHORTLIST.json',research_shortlist)
  rep_summary,top=reporting.create(run/'research_reports',grouped,pairs,private,
      patient_mechanisms,registry_index,shortlist=research_shortlist,
      ai_review_queue=exploration,safety_by_name=safety)
  stage(run,'BASIC_SAFETY_AND_EVIDENCE_REPORTS','completed',
     {'report':str(run/'research_reports/CANDIDATE_REPORT.md'),
      'tier_counts':rep_summary['evidence_tier_counts'],
      'public_label_lookups':public_labels_retrieved,
      'clinical_safety_certified':False})
  stage(run,'EXTENDED_PUBLIC_LABEL_REVIEW','running')
  extended_count=0
  if use_ai and safety_review_count>0:
   extra=expert_triage.safety_review_queue(pairs,set(safety),safety_review_count)
   for candidate in extra:
    if candidate in safety or candidate not in lookup:continue
    try:
     safety[candidate]=deep.safety_gate(candidate,lookup[candidate],private,
                                        prefmap,run/'safety_label_sources')
     extended_count+=1
    except Exception as exc:
     safety[candidate]={'state':'HOLD','flags':['label_lookup_unavailable'],
        'reason':type(exc).__name__,'label_query':candidate}
  save(run/'safety_label_sources/EXTENDED_REVIEW_SUMMARY.json',{
      'requested_additional_reviews':safety_review_count,
      'additional_successful_lookups':extended_count,
      'candidate_count_with_label_record':len(safety),
      'unreviewed_candidates_are_not_safe':True})
  stage(run,'EXTENDED_PUBLIC_LABEL_REVIEW','completed',{
      'additional_successful_lookups':extended_count,'clinical_clearance':False})
  stage(run,'EXPERT_RESEARCH_TRIAGE','running')
  expert_summary,expert_top=expert_triage.report(run/'expert_research',
     pairs,genomic_expert,safety=safety,top=30)
  stage(run,'EXPERT_RESEARCH_TRIAGE','completed',{
     'report':'expert_research/EXPERT_RESEARCH_REPORT.md',
     'bucket_counts':expert_summary['expert_bucket_counts'],
     'prescribing_advice_generated':False})
  stage(run,'PATIENT_FIRST_MECHANISM_PREDICTION','running')
  prediction_summary,top_predictions=predictive.build_reports(
     run/'patient_first_predictions',grouped,pairs,patient_mechanisms,
     registry_index,private,safety_by_name=safety,
     curated_pk_path=v286_pk_reference,top_n=30)
  pred_index={(x['candidate'],x['route_id']):x for x in top_predictions}
  for candidate in research_shortlist:
   found=pred_index.get((candidate['candidate'],candidate['route_id']))
   candidate['patient_first_prediction']=(predictive._safe_summary(found) if found else
      next((x.get('patient_first_prediction') for x in pairs
       if x['ingredient']==candidate['candidate'] and x['route_id']==candidate['route_id']),None))
   candidate['clinical_recommendation']=False
  save(run/'research_shortlist/RESEARCH_SHORTLIST.json',research_shortlist)
  stage(run,'PATIENT_FIRST_MECHANISM_PREDICTION','completed',{
    'approved_drug_groups':prediction_summary['approved_single_ingredient_drug_groups_considered'],
    'research_predictions_not_clinical_efficacy':True,
    'report':'patient_first_predictions/PATIENT_FIRST_REPORT.md'})
  promoted=[x for x in attack if x['verdict']=='PROMOTE_RESEARCH_LEAD'];holds=[x for x in attack if x['verdict']=='HOLD'];kills=[x for x in attack if x['verdict']=='KILL']
  selected_names={x['ingredient'] for x in finalists}|{x['candidate'] for x in research_shortlist}|{x['candidate'] for x in attack}
  selected_registry={n:registry_output.for_candidate(n,registry_index) for n in sorted(selected_names)}
  save(run/'REGULATORY_APPROVALS.json',{'snapshot':'2026-10-10','source':'pinned official registers','candidates':selected_registry})
  printable=['# Rare Disease MVA Pharma — V2.8.10','Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch','', '# Regulatory registers of AI-assessed research candidates','Snapshot: 2026-10-10. These are research hypotheses, not treatment recommendations.','']
  for name,approval in selected_registry.items():
   printable += ['## '+name,*registry_output.text_lines(approval),'']
  (run/'REGULATORY_APPROVALS.md').write_text('\n'.join(printable)+'\n',encoding='utf-8')
  result={'source_resource_bootstrap':'resource_bootstrap/gnomad/GNOMAD_AUTO_RESOURCE_AUDIT.json',
    'source_public_dependency_catalogue':'license_and_redistribution/FREE_LICENSE_AUDIT.json',
    'gnomad_download_status':bootstrap_report.get('completion'),
    'automated_public_resources_enabled':bool(auto_resources),
    'cost_free_default_ai':not use_ai,
    'restricted_illumina_spliceai_used':False,
    'open_source_splice_prediction_state':pangolin_status.get('prediction_state'),
    'open_source_splice_tool_license':'GPL-3.0 external program',
    'enrichr_not_used_by_default':True,
    'specialist_evidence_report':'specialist_reports/SPECIALIST_EVIDENCE_STATUS.md',
    'vep_annotation_state':vep_status['tool_state'],
    'spliceai_annotation_state':spliceai_status['state'],
    'gnomad_v41_coverage_state':gnomad_status['state'],
    'gnomad_v41_full_nuclear_coverage':gnomad_status['v4_1_full_coverage'],
    'population_pk_models_qualified':quant_pk_status.get('models_source_qualified',0),
    'specialist_genomic_variants_externalized':False,
    'source_qualification_is_not_clinical_efficacy':True,
    'uniprot_protein_context':'expert_genetics/UNIPROT_PROTEIN_CONTEXT.json',
    'additional_public_label_lookups':extended_count,
    'expert_genetics_dossier':'expert_genetics/GENETIC_EXPERT_SUMMARY.md',
    'expert_research_report':'expert_research/EXPERT_RESEARCH_REPORT.md',
    'expert_research_bucket_counts':expert_summary['expert_bucket_counts'],
    'patient_first_research_utility_not_clinical_prediction':True,
    'expert_computational_method':'V288_EXPERT_RESEARCH_TRIAGE_1',
    'population_af_source_status':genomic_expert['population_frequency_source']['state'],
    'clinical_genotype_confirmed':False,
    'patient_first_prediction_model':predictive.VERSION,
    'patient_first_report':'patient_first_predictions/PATIENT_FIRST_REPORT.md',
    'patient_first_summary':'patient_first_predictions/PREDICTOR_SUMMARY.json',
    'approved_drug_prediction_groups':prediction_summary['approved_single_ingredient_drug_groups_considered'],
    'patient_first_prediction_is_calibrated_clinical_probability':False,
    'patient_first_allele_mechanisms_remain_hypotheses':True,
    'research_evidence_report':'research_reports/CANDIDATE_REPORT.md',
    'all_research_candidates':'research_reports/CANDIDATE_REPORT_ALL.jsonl',
    'research_evidence_tier_counts':rep_summary['evidence_tier_counts'],
    'clinical_safety_clearance':False,
    'research_evidence_upgrades':len(upgrades),
    'patient_allele_effects_are_hypotheses':True,
    'exploratory_functional_search_count':len(exploration),
    'exploratory_functional_search_statuses':dict(Counter(x['evidence_state'] for x in exploration)),
    'science_gate_version':science.SCHEMA,'scientific_molecular_groups':len(grouped),
   'safety_gate_requires_real_functional_evidence':True,
   'research_signal_not_therapy_validation':True,
   'regulatory_snapshot':'2026-10-10','candidate_regulatory_approvals':selected_registry,'version':VERSION,'project_title':PROJECT_TITLE,'developer_credit':PROJECT_CREDIT,'completed_utc':now(),'case_id':private['case_id'],'primary_genes':sorted(private['genes']),
   'autonomous_from_raw_patient_inputs':case_path==run/'intake/case_model_v2.json','autonomous_from_case_input':True,
   'raw_intake_performed':case_path==run/'intake/case_model_v2.json','scientific_dependency_on_previous_version_outputs':False,
   'master_rebuilt_from_raw_public_sources':True,'master_entities':len(entities),'all_entities_screened':len(ledger),
   'screen_status_counts':ss['status_counts'],'routes_base':len(base),'routes_discovered':len(discovered),'routes_final':len(merged),
   'signal_entities':ss['signal_entities'],'finalist_count':len(finalists),'finalists':[x['ingredient'] for x in finalists],
   'deep_attack_results':attack,'promoted_research_leads':[x['candidate'] for x in promoted],
   'holds':[x['candidate'] for x in holds],'kills':[x['candidate'] for x in kills],
   'research_shortlist_count':len(research_shortlist),'research_shortlist':research_shortlist,
   'research_shortlist_summary':shortlist_summary,
   'enrichr_library_mode':ss.get('enrichr_library_mode'),'openai_budget':ai.budget_status() if use_ai else None,'patient_data_sent_external':False,'clinical_recommendation':False,
   'scalability_contract':{
    'case_specific_inputs':['raw VCF/VCF.GZ','optional clinical DOCX with HPO-coded features','optional prebuilt case JSON for controlled replay'],
    'case_agnostic_engines':['raw VCF QC','GENCODE/MANE local consequence annotation','ClinVar local evidence','HPO phenotype-gene ranking','UniProt context','SIGNOR route generation','FDA/PurpleBook master','ChEMBL mechanisms/activities','RDKit structure','LINCS/GEO','AI source-bound discovery','safety/deep attack'],
    'new_case_requires_code_changes':False}}
  save(run/'FINAL_RESULT.json',result)
  s=load(run/'STATUS.json');s['completed']=True;s['current_stage']='COMPLETED';s['final_result']='FINAL_RESULT.json';save(run/'STATUS.json',s)
  return result
 except Exception as exc:
  save(run/'FAILURE.json',{'utc':now(),'stage':load(run/'STATUS.json').get('current_stage') if (run/'STATUS.json').exists() else None,
   'error':type(exc).__name__,'message':str(exc)[:1000],'traceback':traceback.format_exc()[-5000:]})
  raise

def main():
 ap=argparse.ArgumentParser(description=PROJECT_TITLE+' | '+PROJECT_CREDIT);ap.add_argument('--vcf',default=DEFAULT_VCF);ap.add_argument('--clinical',default=DEFAULT_CLINICAL);ap.add_argument('--case',default=None,help='Optional prebuilt case JSON for controlled replay; raw intake is the default.');ap.add_argument('--out');ap.add_argument('--final-count',type=int,default=5);ap.add_argument('--no-ai',action='store_true');ap.add_argument('--functional-review-count',type=int,default=12,
  help='Number of diverse experimental hypotheses to audit with the public literature (default 12, max 50).')
 ap.add_argument('--safety-review-count',type=int,default=20,
   help='Additional public FDA label lookups for diverse research groups beyond the normal shortlist (default 20, max 100).')
 ap.add_argument('--gnomad-manifest',default=None,
   help='Local, indexed GRCh38 gnomAD v4.1 exome/genome shard manifest; never query a third party with patient alleles.')
 ap.add_argument('--spliceai-manifest',default=None,
   help='DISALLOWED in the free open-research distribution profile (Illumina restrictions).')
 ap.add_argument('--no-auto-resources',action='store_true',
   help='Skip all automatic dependency/data checks for reproducible offline-only tests.')
 ap.add_argument('--no-resource-download',action='store_true',
   help='Check and plan chromosomes without downloading; existing data remain usable.')
 ap.add_argument('--resource-max-gib',type=int,default=240,
   help='Maximum new gnomAD source data to fetch in this run (default 240 GiB; WSL free-space reserve enforced).')
 ap.add_argument('--enable-online-ai',action='store_true',
   help='Optional BYO paid API use only; free/redistributable baseline does not require it.')
 ap.add_argument('--pk-manifest',default=None,
   help='Independently evaluated PK-Sim/qualified IV model manifest. No defaults or child dosing.')
 ap.add_argument('--vep-on-replay',action='store_true',
   help='Run VEP offline even for explicitly prebuilt case JSON; default VEP always runs for raw VCF intake.')
 ap.add_argument('--population-af',default=None,
   help='Optional locally curated, allele-specific GRCh38 population-frequency JSON. Never treats VCF AF as population frequency.');ap.add_argument('--commercial-mode',action='store_true',help='Disable Enrichr LINCS/GEO signals, which require separate commercial permission; does not certify other licences or medical compliance.')
 a=ap.parse_args();
 if not a.vcf and not a.case: ap.error('Supply --vcf or --case with authorized local data')
 if not 0<=a.functional_review_count<=50: ap.error('--functional-review-count must be 0..50')
 if not 0<=a.safety_review_count<=100: ap.error('--safety-review-count must be 0..100')
 if not 1<=a.resource_max_gib<=1500: ap.error('--resource-max-gib must be 1..1500')
 if a.spliceai_manifest:ap.error('Restricted SpliceAI cannot be used in the free CC-BY redistribution profile')
 # Enrichr is restricted: always omit it from the redistributable research baseline.
 os.environ['RDMVA_COMMERCIAL_MODE']='1'
 os.environ['RDMVA_ENRICHR_MODE']='disabled'
 print(PROJECT_TITLE)
 print(PROJECT_CREDIT)
 print('This service/product is using the Human Phenotype Ontology (release v2026-09-02). See https://hpo.jax.org/.')
 result=run(a.vcf,a.clinical,a.case,a.out,a.final_count,(bool(a.enable_online_ai) and not a.no_ai),a.functional_review_count,a.population_af,a.safety_review_count,a.gnomad_manifest,None,a.pk_manifest,a.vep_on_replay,(not a.no_auto_resources),(not a.no_resource_download),a.resource_max_gib)
 # Console output must stay short and safely encodable on Windows cp1252.
 # All detailed research data already live in access-controlled local files.
 summary={'version':result.get('version'),'completed':True,
   'primary_gene_hypotheses':result.get('primary_genes',[]),
   'registered_drug_labels_screened':result.get('master_entities'),
   'research_tiers':result.get('research_evidence_tier_counts',{}),
   'expert_mechanistic_buckets':result.get('expert_research_bucket_counts',{}),
   'genetic_expert_report':result.get('expert_genetics_dossier'),
   'specialist_evidence_report':result.get('specialist_evidence_report'),
   'primary_research_report':result.get('expert_research_report'),
   'prediction_is_clinical_recommendation':False}
 print(json.dumps(summary,indent=2,ensure_ascii=True))
if __name__=='__main__':main()
