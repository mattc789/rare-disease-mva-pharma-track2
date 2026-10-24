"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

"""Recover a completed V2.8 computation after a transient Windows status-write failure."""
from pathlib import Path
from .core import *
from . import ai
import argparse

def finalize(run):
 run=Path(run);s=load(run/'STATUS.json')
 for stage in ('RAW_INTAKE','CASE_MODEL','MASTER_UNIVERSE','ROUTE_GENERATION','UNIVERSAL_SCREEN','PRIORITY_AND_GATES'):
  if s['stages'][stage]['state']!='completed':raise RuntimeError('STAGE_INCOMPLETE:'+stage)
 private=load(run/'case_private.json')
 entities=list(read_jsonl(run/'master/entities.jsonl'))
 ledger=list(read_jsonl(run/'screen/universal_ledger.jsonl'))
 if len(entities)!=2414 or len(ledger)!=2414:raise RuntimeError('COVERAGE_MISMATCH')
 finalists=load(run/'finalists.json')
 attacks=load(run/'finalist_attack/finalist_attack_results.json')
 if len(finalists)!=len(attacks) or {x['ingredient'] for x in finalists}!={x['candidate'] for x in attacks}:
  raise RuntimeError('ATTACK_INCOMPLETE')
 screen=load(run/'screen/screen_summary.json')
 routes=load(run/'routes_final.json')
 base=load(run/'routes_base.json')
 discovered=load(run/'routes_discovered.json')
 result={'version':VERSION,'completed_utc':now(),'case_id':private['case_id'],
  'primary_genes':sorted(private['genes']),'autonomous_from_raw_patient_inputs':True,
  'raw_intake_performed':True,'scientific_dependency_on_previous_version_outputs':False,
  'master_rebuilt_from_raw_public_sources':True,'master_entities':len(entities),
  'all_entities_screened':len(ledger),'screen_status_counts':screen['status_counts'],
  'routes_base':len(base),'routes_discovered':len(discovered),'routes_final':len(routes),
  'signal_entities':screen['signal_entities'],'finalist_count':len(finalists),
  'finalists':[x['ingredient'] for x in finalists],'deep_attack_results':attacks,
  'promoted_research_leads':[x['candidate'] for x in attacks if x['verdict']=='PROMOTE_RESEARCH_LEAD'],
  'holds':[x['candidate'] for x in attacks if x['verdict']=='HOLD'],
  'kills':[x['candidate'] for x in attacks if x['verdict']=='KILL'],
  'openai_budget':ai.budget_status(),'patient_data_sent_external':False,
  'clinical_recommendation':False,'recovered_from_completed_artifacts':True}
 save(run/'FINAL_RESULT.json',result)
 s['stages']['FINALIST_DEEP_ATTACK']={'state':'completed','utc':now(),'detail':{'attacked':len(attacks)}}
 s['completed']=True;s['current_stage']='COMPLETED';s['final_result']='FINAL_RESULT.json'
 save(run/'STATUS_FINAL.json',s)
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',required=True);a=p.parse_args()
 x=finalize(a.run)
 print(dump({'completed':True,'screened':x['all_entities_screened'],'finalists':x['finalists'],'holds':x['holds'],'kills':x['kills']}))
