"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
from collections import defaultdict,Counter
import math,re,urllib.parse
from .core import *
from . import ai,registry_output,science

def prioritize(ledger,routes):
 # Transparent, conservative gates: ranked network hints are exploration,
 # not biologically validated functional restoration of a patient's allele.
 rmap={r['route_id']:r for r in routes}
 pairs=[]
 for entity in ledger:
  if entity.get('status')!='SIGNAL':continue
  groups=defaultdict(list)
  for signal in entity.get('signals',[]):
   rid=signal.get('route_id');route=rmap.get(rid)
   if not route or route.get('state')=='hold':continue
   groups[rid].append(signal)
  for rid,sigs in groups.items():
   route=rmap[rid];sci=route.get('scientific_gate') or {}
   refs={};indirect_hops=[];net=[];omics=[]
   for sig in sigs:
    p=sig.get('provider')
    if p=='ChEMBL_MoA_SIGNOR':
     hop=int(sig.get('depth',99))
     net.append(hop)
    if p=='perturbational_transcriptomics':
     omics.append(sig)
   shortest=min(net or [99])
   direct_moa=shortest==0
   assay_direct=any(s.get('provider')=='ChEMBL_measured_activity' and
      s.get('value_nM') is not None and float(s['value_nM'])<=1000 and
      s.get('direct_target_assay_verified') is True for s in sigs)
   real_rescue=any(s.get('provider')=='patient_disease_functional_rescue' and
       all(s.get(k) is True for k in ('same_drug','functional_endpoint_verified',
           'disease_relevant_model','exposure_supported','source_independently_verified'))
       for s in sigs)
   context_verified=(sci.get('causal_direction_validated_for_disease') is True or
      (route.get('family')=='direct_gene' and sci.get('primary_gene_target_identity') is True and assay_direct))
   if real_rescue and context_verified: tier='A'
   elif assay_direct and context_verified: tier='B'
   elif direct_moa and sci.get('eligible_for_exploratory_screen',True):tier='C'
   else:tier='D'
   # Only one independent contribution per evidence modality / chemical family.
   score=round(
       (8 if real_rescue and context_verified else 0)+
       (4 if assay_direct and context_verified else 0)+
       ({0:2.0,1:0.8,2:0.35,3:0.1}.get(shortest,0))+
       (0.15 if omics else 0)+
       (0.1 if any(s.get('provider')=='rdkit_morgan_similarity' for s in sigs) else 0),4)
   providers=sorted({s.get('provider') for s in sigs if s.get('provider')})
   reasons=[]
   if not real_rescue:reasons.append('no_same_drug_disease_functional_rescue')
   if not assay_direct:reasons.append('no_verified_direct_target_modulation')
   if not context_verified:reasons.append('target_to_disease_direction_not_independently_verified')
   if shortest>=2 and shortest<99:reasons.append('multi_hop_context_sensitive_inference')
   if omics:reasons.append('gene_expression_signature_is_not_functional_protein_rescue')
   pairs.append({'ingredient':entity['canonical_label'],'entity_id':entity['entity_id'],
     'scientific_group_id':entity.get('scientific_group_id'),'member_labels':entity.get('member_labels',[entity['canonical_label']]),
     'route_id':rid,'route':route,'tier':tier,'score':score,'providers':providers,
     'signals':sigs,'combo_only':entity['combo_only'],'modality':entity['modality'],
     'publication_strength':'experimental' if tier in ('A','B') else 'hypothesis_only',
     'evidence_blockers':reasons,'shortest_network_hops':shortest if shortest<99 else None,
     'scientific_gate':sci,'promotion_eligible':tier=='A'})
 pairs.sort(key=lambda x:({'A':0,'B':1,'C':2,'D':3}[x['tier']],
   -x['score'],x['ingredient'],x['route_id']))
 return pairs

def risk_domains(case_private):
 terms=' '.join(str(x.get('term') or '') for x in case_private.get('phenotype_features',[]) if x.get('presence')=='reported_present').lower()
 domains=set()
 if any(x in terms for x in ('renal','kidney','nephro')):domains.add('renal')
 if any(x in terms for x in ('cancer','tumor','tumour','sarcoma','carcinoma','leukemia','lymphoma')):domains.add('cancer_genome')
 if {'BUB1B','CEP57','TRIP13'} & {str(g).upper() for g in case_private.get('genes',{})}:
  domains.add('cancer_genome')
 if any(x in terms for x in ('short stature','failure to thrive','muscle','growth','weight')):domains.add('growth_muscle')
 if any(x in terms for x in ('cardiac','heart','arrhythm','cardiomy')):domains.add('cardiac')
 if any(x in terms for x in ('liver','hepatic')):domains.add('hepatic')
 if any(x in terms for x in ('seizure','neurolog','brain','developmental')):domains.add('neurologic')
 return sorted(domains)

def chembl_pref_names(entities):
 catalog=load(RAW/'chembl/catalog_phase4.json');byid={}
 for x in catalog:
  h=x.get('molecule_hierarchy') or {};pid=h.get('parent_chembl_id') or x['molecule_chembl_id']
  if x.get('pref_name'):byid.setdefault(pid,x['pref_name'])
 return byid

def label(candidate,entity,prefmap,cache_dir):
 cache=Path(cache_dir)/('label_'+re.sub(r'[^A-Za-z0-9]+','_',candidate)[:70]+'.json')
 if cache.exists():return load(cache)
 names=[candidate]+[prefmap[p] for p in entity.get('chembl_parent_ids',[]) if p in prefmap]
 for name in dict.fromkeys(names):
  for field in ['openfda.generic_name','openfda.substance_name']:
   q=f'{field}:"{name}"'
   u='https://api.fda.gov/drug/label.json?'+urllib.parse.urlencode({'search':q,'limit':1})
   try:j=fetch_json(u,15,5_000_000)
   except Exception:continue
   if j.get('results'):
    save(cache,{'query_name':name,'field':field,'result':j['results'][0]});return load(cache)
 save(cache,{'query_name':None,'result':None});return load(cache)

def _section_text(r):
 if not r:return ''
 keys=['boxed_warning','warnings','warnings_and_cautions','adverse_reactions','clinical_pharmacology','contraindications','drug_interactions','nonclinical_toxicology','carcinogenesis_and_mutagenesis_and_impairment_of_fertility']
 vals=[]
 for k in keys:
  v=r.get(k,[])
  if isinstance(v,list):vals+=v
  elif isinstance(v,str):vals.append(v)
 return ' '.join(vals).lower()

def _positive_risk_hit(txt,pattern):
 for m in re.finditer(pattern,txt,re.I):
  lo=max(0,m.start()-80);hi=min(len(txt),m.end()+80);ctx=txt[lo:hi]
  pre=ctx[:m.start()-lo]
  # Do not turn explicit negative toxicology into a positive risk flag.
  if re.search(r'\b(?:not|no|negative for|did not|does not|without)\b[^.;]{0,45}$',pre,re.I) or re.search(r'\bnon[- ]?$',pre,re.I):continue
  return True
 return False

def safety_gate(candidate,entity,case_private,prefmap,cache_dir):
 lab=label(candidate,entity,prefmap,cache_dir);r=lab.get('result');txt=_section_text(r);domains=risk_domains(case_private)
 hard=[];flags=[]
 if not r:return {'state':'HOLD','hard_flags':[],'flags':['public_label_not_resolved'],'risk_domains':domains,'label_query':lab.get('query_name')}
 # Phenotype-domain matches are risk flags, not automatic contraindications. Final KILL requires stronger evidence.
 if 'cancer_genome' in domains and _positive_risk_hit(txt,r'genotoxic|mutagenic|chromosome aberration|clastogenic|secondary malign|carcinogenic'):
  flags.append('genome_or_cancer_risk_requires_review')
 if 'cancer_genome' in domains and _positive_risk_hit(
       txt,r'myelosuppress|bone marrow suppress|neutropeni|cytotoxic|chemotherap|immunosuppress|pancytopeni'):
  flags.append('genome_stability_or_immune_suppression_risk_requires_review')
 if 'renal' in domains and _positive_risk_hit(txt,r'nephrotox|renal failure|acute kidney|kidney injury'):
  flags.append('renal_risk_requires_review')
 if 'growth_muscle' in domains and _positive_risk_hit(txt,r'myopathy|muscle wasting|growth retard|decreased growth|rhabdomyolysis'):
  flags.append('growth_or_muscle_risk_requires_review')
 if 'cardiac' in domains and _positive_risk_hit(txt,r'qt prolong|torsade|cardiomyopathy|heart failure'):
  flags.append('cardiac_risk_requires_review')
 if 'hepatic' in domains and _positive_risk_hit(txt,r'hepatic failure|hepatotoxic|liver failure'):
  flags.append('hepatic_risk_requires_review')
 if r.get('boxed_warning'):flags.append('boxed_warning_requires_review')
 if r.get('contraindications'):flags.append('contraindications_section_requires_case_review')
 state='HOLD' if flags else 'PASS'
 return {'state':state,'hard_flags':hard,'flags':sorted(set(flags)),'risk_domains':domains,'label_query':lab.get('query_name')}

def deterministic_gates(pairs,entities,case_private,cache_dir):
 emap={e['canonical_label']:e for e in entities};out=[]
 # Load the ChEMBL preferred-name catalogue only when an A/B candidate
 # actually requires immediate safety adjudication. Source-only/offline
 # regression runs with no A/B candidates must not depend on bootstrapped data.
 pref=chembl_pref_names(entities) if any(p.get('tier') in ('A','B') for p in pairs) else {}
 # only A/B need immediate deep safety, C/D remain discovery queue.
 seen=set()
 for p in pairs:
  if p['tier'] not in ('A','B') or p['ingredient'] in seen:continue
  seen.add(p['ingredient']);e=emap[p['ingredient']];sg=safety_gate(p['ingredient'],e,case_private,pref,cache_dir)
  route=p['route'];haz=[]
  genes={t.get('gene','').upper() for t in route.get('targets',[])}
  if genes & {'AURKB','PLK1','CDK1','CDC20'}:haz.append('global_cell_cycle_modulation_requires_functional_rescue')
  if route.get('state')=='hold' or route.get('case_gate_reasons'):haz.extend(route.get('case_gate_reasons',[]))
  evidence_state='HOLD' # until same-drug functional attack
  final='CLOSED' if sg['state']=='CLOSED' else 'HOLD'
  search_aliases=[p['ingredient']]+[pref[pid] for pid in e.get('chembl_parent_ids',[]) if pid in pref and pref[pid]!=p['ingredient']]
  out.append({**p,'candidate_search_aliases':list(dict.fromkeys(search_aliases)),'safety_gate':sg,'route_hazards':sorted(set(haz)),'evidence_gate':{'state':evidence_state,'reason':'same_drug_functional_rescue_not_yet_verified'},
   'final_status':final,'clinical_recommendation':False})
 return out

def choose_finalists(gated,count=5):
 eligible=[x for x in gated if x['final_status']=='HOLD']
 eligible.sort(key=lambda x:({'A':0,'B':1}.get(x['tier'],2),len(x['safety_gate']['flags'])+len(x['route_hazards']),-x['score'],x['ingredient']))
 return eligible[:count]

def attack_finalists(finalists,primary_genes,outdir,registry_index=None):
 out=Path(outdir);results=[]
 gene=primary_genes[0] if primary_genes else ''
 for i,c in enumerate(finalists,1):
  d=out/f'{i:02d}_{re.sub(r"[^A-Za-z0-9]+","_",c["ingredient"])[:50]}';d.mkdir(parents=True,exist_ok=True)
  approval=registry_output.for_candidate(c['ingredient'],registry_index or {})
  targets=[t.get('gene') for t in c['route'].get('targets',[]) if t.get('gene')]
  docs=ai.literature(gene,d/'literature',max_docs=5,candidate=c.get('candidate_search_aliases') or c['ingredient'],targets=targets)
  if not docs:
   v={'candidate':c['ingredient'],'adjudication':'INSUFFICIENT','same_drug_functional_rescue':False,'direct_route_modulation':False,
      'exposure_reported':False,'counterevidence_present':False,'support_quotes':[],'counter_quotes':[],
      'reason':'No rights-admitted candidate-specific public full text was acquired.','limitations':['literature_coverage_missing']}
   final=v;sol=None
  else:
   try:v=ai.deep_attack(c['ingredient'],c['route'],docs,d,'verify',regulatory_approval=approval)
   except Exception as exc:
    v={'candidate':c['ingredient'],'adjudication':'INSUFFICIENT','same_drug_functional_rescue':False,'direct_route_modulation':False,
       'exposure_reported':False,'counterevidence_present':False,'support_quotes':[],'counter_quotes':[],
       'reason':'Verifier output unavailable or invalid.','limitations':[type(exc).__name__]}
   sol=None;final=v
   if v['adjudication']=='SUPPORTED' and v['same_drug_functional_rescue']:
    try:sol=ai.deep_attack(c['ingredient'],c['route'],docs,d,'final',regulatory_approval=approval);final=sol
    except Exception:sol=None
  if c['safety_gate']['state']=='CLOSED':verdict='KILL'
  elif c.get('promotion_eligible') and final['adjudication']=='SUPPORTED' and final['same_drug_functional_rescue'] and not c['route_hazards']:verdict='PROMOTE_RESEARCH_LEAD'
  elif final['adjudication']=='OPPOSED':verdict='KILL'
  else:verdict='HOLD'
  rec={'candidate':c['ingredient'],'regulatory_approval':approval,'route_id':c['route_id'],'tier':c['tier'],'screen_score':c['score'],'safety_gate':c['safety_gate'],
    'route_hazards':c['route_hazards'],'terra':v,'sol':sol,'verdict':verdict,'clinical_recommendation':False}
  save(d/'final_attack.json',rec);results.append(rec)
 save(out/'finalist_attack_results.json',results);return results

def _public_blockers_from_attack(a):
 v=a.get('terra') or {}
 out=[]
 if v.get('adjudication') and v.get('adjudication')!='SUPPORTED':out.append('deep_adjudication_'+v['adjudication'].lower())
 if not v.get('same_drug_functional_rescue',False):out.append('same_drug_functional_rescue_absent')
 if not v.get('direct_route_modulation',False):out.append('direct_route_modulation_unverified')
 if not v.get('exposure_reported',False):out.append('exposure_relevance_not_established')
 if v.get('counterevidence_present',False):out.append('counterevidence_present')
 return out

def _fallback_shortlist_rationale(x):
 providers=', '.join(x.get('providers',[])) or 'computational evidence'
 if x['shortlist_class']=='DEEPEN-HIGH':
  why_not='The candidate reached deep review but did not demonstrate same-drug functional rescue of the nominated route.'
  why_yes=f"It remains worth investigating because it had one of the strongest surviving route-linked evidence packages ({providers}) and was not mechanistically disproven."
  priority='HIGH'
 elif x['shortlist_class']=='DEEPEN-MEDIUM':
  why_not='The candidate had A/B-level screening evidence but ranked below the finalists and has not yet received candidate-specific functional rescue verification.'
  why_yes=f"It remains worth investigating because its screening evidence is comparatively strong ({providers}) and survived deterministic prioritization."
  priority='MEDIUM'
 else:
  why_not='The evidence is currently indirect or below the A/B threshold and has not been candidate-specifically verified.'
  why_yes=f"It remains hypothesis-generating because at least one route-compatible signal was detected ({providers})."
  priority='EXPLORATORY'
 return {'candidate':x['candidate'],'route_id':x['route_id'],'why_not_first_choice':why_not,'why_worth_deepening':why_yes,
  'missing_evidence':['direct candidate-to-route modulation','same-drug functional rescue','exposure/concentration relevance'],
  'next_checks':['candidate-specific literature and assay review','functional rescue assay if mechanistic evidence strengthens'],
  'research_priority':priority,'confidence':'moderate' if priority!='EXPLORATORY' else 'low'}

def build_research_shortlist(pairs,gated,finalists,attacks,outdir,count=10,use_ai=True,entities=None,registry_index=None):
 out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
 entity_by_id={e['entity_id']:e for e in (entities or [])}
 pair_by={(p['ingredient'],p['route_id']):p for p in pairs}
 best_by={}
 for p in pairs:
  q=best_by.get(p['ingredient'])
  if q is None or ({'A':0,'B':1,'C':2,'D':3}[p['tier']],-p['score'])<({'A':0,'B':1,'C':2,'D':3}[q['tier']],-q['score']):best_by[p['ingredient']]=p
 gate_by={g['ingredient']:g for g in gated}
 attack_by={a['candidate']:a for a in attacks}
 finalist_names={x['ingredient'] for x in finalists}
 promoted={a['candidate'] for a in attacks if a['verdict']=='PROMOTE_RESEARCH_LEAD'}
 killed={a['candidate'] for a in attacks if a['verdict']=='KILL'}
 rows=[];seen=set();seen_identity=set()
 def add(candidate,klass,preferred_route_id=None):
  if candidate in seen or candidate in promoted or candidate in killed:return
  p0=pair_by.get((candidate,preferred_route_id)) if preferred_route_id else None
  p0=p0 or best_by.get(candidate)
  e=entity_by_id.get(p0.get('entity_id')) if p0 else None
  parents=tuple(sorted(e.get('chembl_parent_ids',[]))) if e else ()
  identity=('CHEMBL',parents) if len(parents)==1 else ('NAME',re.sub(r'-[A-Z]{4}$','',candidate.upper()))
  if identity in seen_identity:return
  a=attack_by.get(candidate);g=gate_by.get(candidate)
  rid=(a or g or p0 or {}).get('route_id')
  p=pair_by.get((candidate,rid)) or p0 or best_by.get(candidate)
  if not p:return
  blockers=_public_blockers_from_attack(a) if a else ['same_drug_functional_rescue_not_verified','candidate_specific_deep_review_not_completed']
  local_flags=['patient_specific_drug_safety_review_not_performed',
    'mechanistic_functional_rescue_not_verified'] if not g else []
  if g:
   local_flags=sorted(set((g.get('safety_gate') or {}).get('flags',[])+g.get('route_hazards',[])))
  deep_public=None
  if a:
   v=a.get('terra') or {}
   deep_public={k:v.get(k) for k in ['adjudication','same_drug_functional_rescue','direct_route_modulation','exposure_reported','counterevidence_present','reason','limitations']}
  rows.append({'candidate':candidate,'regulatory_approval':registry_output.for_candidate(candidate,registry_index or {}),'route_id':p['route_id'],'shortlist_class':klass,'tier':p['tier'],'screen_score':p['score'],
   'route':p['route'],'providers':p['providers'],'signals':p['signals'],'public_blockers':blockers,'deep_public':deep_public,
   'local_case_fit_caveats':local_flags,'case_specific_review_required':bool(local_flags),'clinical_recommendation':False})
  seen.add(candidate);seen_identity.add(identity)
 # Deep-attacked survivors are the strongest second choices.
 for a in attacks:
  if a['verdict']=='HOLD':add(a['candidate'],'DEEPEN-HIGH')
 # A/B candidates that were not among the finalists come next.
 for g in gated:
  if g['ingredient'] not in finalist_names and g.get('final_status')=='HOLD':add(g['ingredient'],'DEEPEN-MEDIUM')
 # Fill remaining slots with the best C-tier hypotheses, never D by default.
 for p in science.diverse_exploratory_review_queue(pairs,count=count):
  if len(rows)>=count:break
  if p['tier'] in ('C','D'):add(p['ingredient'],'DEEPEN-EXPLORATORY',p['route_id'])
 # If route diversity yields fewer than requested, fill from remaining groups.
 for p in pairs:
  if len(rows)>=count:break
  if p['tier'] in ('C','D'):add(p['ingredient'],'DEEPEN-EXPLORATORY',p['route_id'])
 rows=rows[:count]
 ai_rows=[]
 if use_ai and rows:
  try:ai_rows=ai.research_shortlist_rationales(rows,out)
  except Exception as exc:save(out/'AI_SHORTLIST_ERROR.json',{'error':type(exc).__name__,'message':str(exc)[:500]})
 amap={(x['candidate'],x['route_id']):x for x in ai_rows}
 final=[]
 for x in rows:
  rationale=amap.get((x['candidate'],x['route_id'])) or _fallback_shortlist_rationale(x)
  if x.get('tier')=='D':
   rationale={**rationale,'research_priority':'EXPLORATORY','confidence':'low',
     'disclaimer':'Indirect hypothesis only; no case-relevant functional rescue demonstrated.'}
  y={**x,'ai_rationale':{**rationale,'regulatory_approval':x['regulatory_approval']},
     'scientific_validation':'unverified' if x.get('tier') in ('C','D') else 'requires_external_replication',
     'ai_patient_context_used':False}
  final.append(y)
 save(out/'RESEARCH_SHORTLIST.json',final)
 summary={'count':len(final),'classes':dict(Counter(x['shortlist_class'] for x in final)),
  'candidates':[x['candidate'] for x in final],'clinical_recommendation':False,'patient_data_sent_to_ai':False}
 save(out/'RESEARCH_SHORTLIST_SUMMARY.json',summary)
 # Human-readable export for challenge/reporting.
 lines=['# Research Shortlist','', 'These are research hypotheses, not treatment recommendations.','']
 for i,x in enumerate(final,1):
  r=x['ai_rationale'];lines += [f"## {i}. {x['candidate']} — {x['shortlist_class']}",
   f"Route: {x['route_id']} / {x['route'].get('target_concept')}",
   *registry_output.text_lines(x["regulatory_approval"]),
   f"Why not first choice: {r['why_not_first_choice']}",
   f"Why deepen: {r['why_worth_deepening']}",
   'Missing evidence: '+'; '.join(r.get('missing_evidence',[])),
   'Next checks: '+'; '.join(r.get('next_checks',[])),
   'Local case-fit caveats: '+('; '.join(x['local_case_fit_caveats']) if x['local_case_fit_caveats'] else 'none flagged by current local gates'),'']
 (out/'RESEARCH_SHORTLIST.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 return final,summary
