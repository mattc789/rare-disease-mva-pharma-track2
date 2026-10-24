"""Conservative causal, chemical-identity and evidence-provenance gates.

Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch

These gates prevent network proximity, gene expression, repeated salt forms,
and unverified model judgments from masquerading as disease functional rescue.
Neither drugs nor endpoints are hard-coded as positive discoveries.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from copy import deepcopy
import hashlib,re
from .core import norm_name

SCHEMA='SCIENCE_GATES_V284_1'
HIGH_RISK_CELL_CYCLE={'CDK1','CDK2','PLK1','AURKB','CDC20','AURKA','BUB1B'}
DISALLOWED_CONTEXT={'meiotic','oocyte_only','nonhuman_only'}
NEGATIVE_ACTION_WORDS=r'inhibit(?:s|ed|ing|ion)?|block(?:s|ed|ing)?|deplet(?:e|ed|ion)|knockdown|suppress(?:ion|es|ed)?|decreas(?:e|es|ed|ing)'

def _word_target(target):
 if not target:return None
 if not re.fullmatch(r'[A-Z0-9_-]{3,20}',target):return None
 return r'\b'+re.escape(target)+r'\b'

def adjudicate_routes(routes):
 """Use source gates without rewriting or fabricating counterfactual biology."""
 out=[];audits=[]
 for route in routes:
  r=deepcopy(route)
  reasons=set(r.get('case_gate_reasons') or [])
  state=r.get('state','conditional')
  target_gene=(r.get('targets') or [{}])[0].get('gene','').upper()
  desc=str(r.get('mechanism_concept') or '')
  direction=str(r.get('desired_action') or '')
  family=r.get('family')
  if direction not in ('increase','decrease'):reasons.add('unspecified_functional_direction')
  if family=='downstream_compensation':
   reasons.add('downstream_sign_is_not_validated_disease_functional_rescue')
  if family=='signed_upstream':
   reasons.add('signor_activity_sign_requires_context_specific_validation')
  if family=='literature_discovered':
   if not r.get('evidence'):reasons.add('literature_route_missing_quoted_evidence')
   verdict=r.get('route_direction_verification') or {}
   if verdict.get('verdict')=='REJECT':reasons.add('source_direction_rejected')
   if verdict.get('evidence_context')=='meiotic':reasons.add('meiotic_transfer_not_established')
  # Directional consistency is tested against actionable imperative in prose.
  pat=_word_target(target_gene)
  if pat and direction=='increase':
   harmful=rf'\b(?:{NEGATIVE_ACTION_WORDS})\s+(?:the\s+)?{pat}'
   if re.search(harmful,desc,re.I):
    reasons.add('mechanism_text_conflicts_with_increase_action')
  if pat and direction=='decrease':
   contradictory=rf'\b(?:increase|increase\s+the|overexpress|activate|stabilize|restore)\s+{pat}'
   if re.search(contradictory,desc,re.I):
    reasons.add('mechanism_text_conflicts_with_decrease_action')
  if target_gene=='CDC20' and family=='downstream_compensation':
   reasons.add('cdc20_activation_not_equivalent_to_checkpoint_rescue')
  if target_gene=='APC' and family=='downstream_compensation':
   reasons.add('APC_P25054_distinct_from_APC_C_complex')
  if r.get('evidence_context') in DISALLOWED_CONTEXT:
   reasons.add('nontransferrable_model_context')
  blockers={
   'source_direction_rejected','mechanism_text_conflicts_with_increase_action',
   'mechanism_text_conflicts_with_decrease_action',
   'downstream_sign_is_not_validated_disease_functional_rescue',
   'cdc20_activation_not_equivalent_to_checkpoint_rescue',
   'APC_P25054_distinct_from_APC_C_complex',
   'meiotic_transfer_not_established','nontransferrable_model_context',
   'unspecified_functional_direction','literature_route_missing_quoted_evidence'
  }
  if reasons & blockers or state=='hold':state='hold'
  # A route can remain conditional for exploratory hypothesis mining, but
  # only curated same-drug functional validation can raise tier A.
  r['state']=state
  primary_target_identity=(family=='direct_gene' and target_gene==str(r.get('target_concept') or '').upper())
  r['scientific_gate']={'schema':SCHEMA,
   'primary_gene_target_identity':primary_target_identity,
   'allele_specific_rescue_verified':False,
   'eligible_for_exploratory_screen':state!='hold',
   'functional_rescue_validated':False,
   'causal_direction_validated_for_disease':False,
   'requires_independent_endpoint':True,
   'hard_hold':state=='hold','issues':sorted(reasons)}
  out.append(r)
  audits.append({'route_id':r.get('route_id'),'target':r.get('target_concept'),
    'state':state,'issues':sorted(reasons)})
 return out,audits

def molecular_group_id(entry):
 parents=entry.get('chembl_parent_ids') or []
 if len(parents)==1 and str(parents[0]).startswith('CHEMBL'):
  return 'ChEMBL_PARENT:'+parents[0]
 # Do not guess salts, synonyms, stereochemistry or biosimilar identity.
 return 'NAME:'+norm_name(entry['canonical_label'])

def _representative(records):
 def quality(e):
  n=e['canonical_label']
  return (0 if len(e.get('chembl_parent_ids') or [])==1 else 1,
   0 if e.get('global_approved_as_single') else 1,
   len(n),n)
 return min(records,key=quality)

def _signal_key(s):
 prov=s.get('provider')
 rid=s.get('route_id')
 if prov=='ChEMBL_MoA_SIGNOR':
  # Multiple mechanisms of the same parent drug at the same signed path
  # do not constitute independent experimental proof.
  return (rid,prov,s.get('depth'))
 if prov=='perturbational_transcriptomics':
  return (rid,prov)
 if prov=='ChEMBL_measured_activity':
  return (rid,prov,s.get('assay_id'),s.get('value_nM'))
 if prov=='rdkit_morgan_similarity':
  return (rid,prov,s.get('reference_molecule'))
 return (rid,prov)

def canonicalize_ledger(ledger):
 groups=defaultdict(list)
 for e in ledger:groups[molecular_group_id(e)].append(e)
 output=[];duplicate_groups=[]
 for group_id,members in sorted(groups.items()):
  rep=_representative(members);r=deepcopy(rep)
  best={};rawsignals=0
  for m in members:
   for sig in m.get('signals',[]):
    rawsignals+=1;k=_signal_key(sig)
    current=best.get(k)
    if current is None or float(sig.get('tanimoto') or 0)>float(current.get('tanimoto') or 0):
     best[k]=deepcopy(sig)
  labels=sorted({m['canonical_label'] for m in members})
  signals=sorted(best.values(),key=lambda x:(str(x.get('route_id')),str(x.get('provider')),str(x.get('depth'))))
  r.update({'scientific_group_id':group_id,'member_labels':labels,
   'member_entity_ids':sorted({m['entity_id'] for m in members}),
   'formulation_sensitive':len(labels)>1,
   'scientific_identity_method':'ChEMBL exact parent hierarchy' if group_id.startswith('ChEMBL_PARENT:') else 'exact label only',
   'signals':signals,'signal_count':len(signals),'raw_formulation_signal_count':rawsignals,
   'approval_jurisdictions_by_form':{m['canonical_label']:m.get('approved_jurisdictions') or [] for m in members},
   'parent_id':group_id.split(':',1)[1] if group_id.startswith('ChEMBL_PARENT:') else None,
   'status':'SIGNAL' if signals else (
     'INCOMPATIBLE_MODALITY' if all(m.get('status')=='INCOMPATIBLE_MODALITY' for m in members) else
     'NO_SIGNAL' if any(m.get('status')=='NO_SIGNAL' for m in members) else 'INSUFFICIENT_DATA')})
  output.append(r)
  if len(labels)>1:duplicate_groups.append({'group_id':group_id,'name':r['canonical_label'],'labels':labels})
 summary={'schema':SCHEMA,'raw_labels':len(ledger),'grouped_molecular_families':len(output),
  'multi_form_groups':len(duplicate_groups),'duplicate_examples':duplicate_groups[:20],
  'signal_groups':sum(x['status']=='SIGNAL' for x in output),
  'grouping_note':'ChEMBL grouping is chemical-parent related, NOT clinical bioequivalence or dose interchangeability.'}
 return output,summary

def _study_groups(examples):
 """Conservatively count unresolved GEO GDS signals as ONE source family."""
 seen=set();programs=set();contexts=set();species=set()
 for x in examples or []:
  term=str(x.get('term') or '')
  lib=str(x.get('library') or '')
  lower=term.lower()
  if 'LINCS' in lib:
   m=re.match(r'^(\S+)\s+(\S+)\s+(\d+(?:\.\d+)?H)-',term,re.I)
   context=(m.group(2).upper() if m else 'unspecified')
   programs.add('LINCS_L1000')
   contexts.add('LINCS:'+context)
   # Distinct cell lines/perturbation conditions are not automatically
   # independent published validation studies.
   seen.add('LINCS:L1000:single_data_collection')
  elif 'GEO' in lib:
   m=re.search(r'\bgse\d+\b',term,re.I)
   if m:
    seen.add('GEO:'+m.group().upper())
   else:
    seen.add('GEO:UNKNOWN_PARENT_SERIES')
   contexts.add('GEO:'+re.search(r'\bgpl\d+\b',lower).group() if re.search(r'\bgpl\d+\b',lower) else 'GEO:unspecified')
   for organism in ('homo sapiens','mus musculus','rattus norvegicus'):
    if organism in lower:species.add(organism)
   programs.add('GEO')
 return {'source_families':sorted(seen),'source_programs':sorted(programs),
   'experiment_contexts':sorted(contexts),'species_observed':sorted(species),
   'independent_validation_studies':0,
   'reason':'Source series and replicate independence not fully established; no validated rescue study.'}

def enrich_perturbation_signal(signal,tx_record=None):
 s=dict(signal)
 if s.get('provider')!='perturbational_transcriptomics':return s
 tx_record=tx_record or {}
 evidence=_study_groups(tx_record.get('examples') or [])
 s['raw_signature_count']=tx_record.get('support_count',s.get('support_count',0))
 s['signature_evidence']=evidence
 s['independent_validation_studies']=0
 s['functional_rescue']=False
 s['experimental_gene_expression_only']=True
 s['study_independence_unresolved']=True
 return s

def attach_transcript_context(ledger,transcript_rows):
 refs={(x.get('ingredient'),x.get('route_id')):x for x in transcript_rows}
 for entry in ledger:
  forms=entry.get('member_labels') or [entry['canonical_label']]
  augmented=[]
  for s in entry.get('signals',[]):
   if s.get('provider')=='perturbational_transcriptomics':
    exemplars=[refs.get((f,s.get('route_id'))) for f in forms]
    exemplar=next((x for x in exemplars if x),None)
    augmented.append(enrich_perturbation_signal(s,exemplar))
   else:augmented.append(s)
  entry['signals']=augmented
 return ledger

def merge_registry_by_family(ledger,original_index):
 """Display source authorisations PER FORM, never claim form equivalence."""
 merged=dict(original_index)
 for e in ledger:
  labels=e.get('member_labels') or [e['canonical_label']]
  records=[];rnames=set();jurisdictions=set();solo=False;combo=False
  for label in labels:
   raw=original_index.get(norm_name(label)) or {}
   for row in raw.get('records',[]):
    item={**row,'specific_ingredient_form':label}
    records.append(item)
    rnames.add(row['register']);jurisdictions.add(row['jurisdiction'])
    solo |= row['product_context']=='single_ingredient_product'
    combo |= row['product_context']=='combination_component'
  # Dedupe only truly identical authorisation evidence tied to same form.
  unique={}
  for item in records:
   k=(item['specific_ingredient_form'],item['register'],item['jurisdiction'],
      item['status'],item['product_context'])
   unique[k]=item
  grouped_record={'query_name':e['canonical_label'],'match_method':'ChEMBL parent group with form-specific authorisation',
    'snapshot_date':'2026-10-10','records':sorted(unique.values(),key=lambda x:(x['register'],x['jurisdiction'],x['specific_ingredient_form'])),
    'registry_names':sorted(rnames),'jurisdictions':sorted(jurisdictions),
    'single_product_listed':solo,'combination_component_listed':combo,
    'overall_state':'single_ingredient_product_listed' if solo else 'combination_only' if combo else 'not_in_snapshot',
    'source_forms':labels,'scientific_group_id':e.get('scientific_group_id'),
    'warning':'Individual forms may have distinct authorisations, exposure, and routes; no blanket approval for the family.'}
  merged[norm_name(e['canonical_label'])]=grouped_record
 return merged

def diverse_exploratory_review_queue(pairs,count=5):
 """One plausible unvalidated query per route family before repeated analogues."""
 if count<=0:return []
 candidates=[p for p in pairs if p.get('tier') in ('C','D')]
 byroute=defaultdict(list)
 for p in candidates:byroute[p['route_id']].append(p)
 for group in byroute.values():
  group.sort(key=lambda p:(
    {'T1':0,'T2':1,'T3':2,'T4':3,'HOLD':4}.get(p.get('expert_investigation_bucket'),5),
    -p['score'],p['ingredient']))
 routes=sorted(byroute,key=lambda rid:(
    -max((float(p.get('patient_first_score',p.get('score',0))) for p in byroute[rid]),default=0.0),
    {'direct_gene':0,'literature_discovered':1,'proteostasis':2,'signed_upstream':3}.get(
      byroute[rid][0]['route'].get('family'),4),rid))
 queue=[];used=set()
 for rid in routes:
  for p in byroute[rid]:
   key=p.get('scientific_group_id') or p['ingredient']
   if key in used:continue
   used.add(key)
   queue.append(p)
   break
  if len(queue)>=count:return queue
 for p in candidates:
  key=p.get('scientific_group_id') or p['ingredient']
  if key not in used:
   queue.append(p);used.add(key)
   if len(queue)>=count:break
 return queue

def exploratory_functional_evidence_search(pairs,primary_gene,outdir,
                                            registry_index=None,count=5,use_ai=True):
 """Source-bound literature falsification, not clinical or automatic promotion.

 A missing qualifying paper never becomes positive evidence. Even an AI
 'SUPPORTED' response is only flagged for human source and assay audit.
 """
 from pathlib import Path
 from .core import save
 from . import ai,registry_output
 out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
 queue=diverse_exploratory_review_queue(pairs,count)
 results=[]
 for i,p in enumerate(queue,1):
  candidate=p['ingredient'];r=p['route']
  dirname=f"{i:02d}_{re.sub(r'[^A-Za-z0-9]+','_',candidate)[:55]}"
  place=out/dirname;place.mkdir(parents=True,exist_ok=True)
  regulatory=registry_output.for_candidate(candidate,registry_index or {})
  result={'candidate':candidate,'scientific_group_id':p.get('scientific_group_id'),
          'drug_forms':p.get('member_labels',[]),'source_route_id':p['route_id'],
          'priority_level':'EXPLORATORY','research_gene':r.get('case_gene') or primary_gene,
          'clinical_recommendation':False,
          'eligible_for_treatment_claim':False,
          'regulatory_approval':regulatory,
          'evidence_state':'NOT_YET_ASSESSED','assay_verification':'requires_external_validation',
          'original_network_hops':p.get('shortest_network_hops'),
          'sources_reviewed':0}
  if not use_ai:
   result['evidence_state']='REVIEW_DISABLED'
  else:
   try:
    target_names=[t.get('gene') for t in r.get('targets',[]) if t.get('gene')]
    primary_gene_for_route=str(r.get('case_gene') or primary_gene)
    public_docs=ai.literature(primary_gene_for_route,place/'literature',
       max_docs=5,candidate=p.get('member_labels') or candidate,targets=target_names)
    result['sources_reviewed']=len(public_docs)
    if not public_docs:
     result['evidence_state']='NO_RIGHTS_ADMITTED_SOURCE_FOUND'
    else:
     v=ai.deep_attack(candidate,r,public_docs,place,'verify',regulatory_approval=regulatory)
     result['ai_evaluation']=v
     claims=(v.get('adjudication')=='SUPPORTED' and
       v.get('same_drug_functional_rescue') is True and
       v.get('direct_route_modulation') is True and
       v.get('exposure_reported') is True and
       bool(v.get('support_quotes')))
     if claims:
      result['evidence_state']='POTENTIAL_SOURCE_TO_VALIDATE_INDEPENDENTLY'
     elif v.get('adjudication')=='OPPOSED':
      result['evidence_state']='COUNTEREVIDENCE_TO_EXAMINE'
     else:result['evidence_state']='DIRECT_FUNCTIONAL_EVIDENCE_NOT_VERIFIED'
   except Exception as exc:
    result['evidence_state']='RESEARCH_ACQUISITION_OR_AI_ERROR'
    result['error_type']=type(exc).__name__
    result['error']=str(exc)[:250]
  # AI source checks cannot substitute experimental validation, even when supported.
  save(place/'REVIEW.json',result)
  results.append(result)
 save(out/'RESEARCH_FALSIFICATION_QUEUE.json',results)
 save(out/'SUMMARY.json',{'count':len(results),
   'statuses':dict(Counter(r['evidence_state'] for r in results)),
   'promoted_validated_treatments':0,
   'no_patient_data_sent_to_external_models':True,
   'automated_source_search_is_not_independent_experimental_reproduction':True})
 return results
