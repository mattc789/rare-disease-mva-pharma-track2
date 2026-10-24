"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
from collections import defaultdict,Counter
import json,re,csv,urllib.parse,math,os,concurrent.futures
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from .core import *
from .routes import signor_edges

CHEMBL_CATALOG=RAW/'chembl/catalog_phase4.json'
CHEMBL_MECH=RAW/'chembl/mechanisms_records.json'
CHEMBL_TARGETS=RAW/'chembl/target_metadata_records.json'
FPGEN=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
SPECIAL={'vaccine','immune_globulin','blood_or_coagulation_product','gene_or_viral_therapy_review','cell_therapy_review'}
ACTION_SIGN={'INHIBITOR':-1,'ANTAGONIST':-1,'BLOCKER':-1,'DEGRADER':-1,'RNAI INHIBITOR':-1,'ANTISENSE INHIBITOR':-1,
             'AGONIST':1,'ACTIVATOR':1,'OPENER':1,'POSITIVE ALLOSTERIC MODULATOR':1,'NEGATIVE ALLOSTERIC MODULATOR':-1}
LIBS={'LINCS_up':('LINCS_L1000_Chem_Pert_up','up'),'LINCS_down':('LINCS_L1000_Chem_Pert_down','down'),'GEO':('Drug_Perturbations_from_GEO_2014','mixed')}

def raw_identity_maps(entities):
 parent_to_labels=defaultdict(set);alias=defaultdict(set)
 for e in entities:
  lab=e['canonical_label'];alias[norm_token(lab)].add(lab)
  for p in e.get('chembl_parent_ids',[]):parent_to_labels[p].add(lab)
 catalog=load(CHEMBL_CATALOG)
 for x in catalog:
  h=x.get('molecule_hierarchy') or {};pid=h.get('parent_chembl_id') or x['molecule_chembl_id']
  labs=parent_to_labels.get(pid,set())
  if not labs:continue
  vals=[x.get('pref_name')]+[z.get('molecule_synonym') for z in x.get('molecule_synonyms',[])]+[z.get('synonyms') for z in x.get('molecule_synonyms',[])]
  for a in vals:
   n=norm_token(a)
   if len(n)>=5:alias[n].update(labs)
 return parent_to_labels,alias

def target_index(cache_dir,routes):
 cache=Path(cache_dir)/'target_index.json'
 if cache.exists():return load(cache)
 idx={}
 for t in load(CHEMBL_TARGETS):
  if t.get('organism')!='Homo sapiens':continue
  comps=t.get('target_components') or []
  accs={x.get('accession') for x in comps if x.get('accession')}
  if len(accs)==1:idx.setdefault(next(iter(accs)),[]).append(t['target_chembl_id'])
 # fill missing route accessions via API
 needed={t.get('accession') for r in routes for t in r.get('targets',[]) if t.get('accession')}
 for acc in sorted(needed):
  if acc in idx:continue
  try:
   u='https://www.ebi.ac.uk/chembl/api/data/target.json?'+urllib.parse.urlencode({'target_components__accession':acc,'limit':100})
   j=fetch_json(u,20,8_000_000)
   ids=[]
   for t in j.get('targets',[]):
    comps={x.get('accession') for x in t.get('target_components',[]) if x.get('accession')}
    if t.get('tax_id')==9606 and t.get('target_type')=='SINGLE PROTEIN' and comps=={acc}:ids.append(t['target_chembl_id'])
   if ids:idx[acc]=sorted(set(ids))
  except Exception as exc:raise RuntimeError('V286_CHEMBL_TARGET_UNAVAILABLE:'+str(acc)) from exc
 save(cache,idx);return idx

def mechanisms():
 out=defaultdict(list);tmeta={x['target_chembl_id']:x for x in load(CHEMBL_TARGETS)}
 for m in load(CHEMBL_MECH):
  sign=ACTION_SIGN.get((m.get('action_type') or '').upper())
  if sign not in (-1,1):continue
  t=tmeta.get(m.get('target_chembl_id'));accs=[]
  if t:
   if t.get('organism')!='Homo sapiens':continue
   accs=[x.get('accession') for x in t.get('target_components',[]) if x.get('accession')]
   if len(set(accs))!=1:continue
  if not accs:continue
  p=m.get('parent_molecule_chembl_id') or m.get('molecule_chembl_id')
  out[p].append({'mechanism_id':m.get('mec_id'),'target_id':m.get('target_chembl_id'),'accession':accs[0],
    'action_sign':sign,'action_type':m.get('action_type'),'mechanism':m.get('mechanism_of_action'),'direct':bool(m.get('direct_interaction'))})
 return out

def signed_reverse(max_hops,routes):
 rev=defaultdict(list)
 for e in signor_edges():
  if e.get('direct'):rev[e['target_acc']].append(e)
 result={}
 for r in routes:
  m=defaultdict(list)
  for t in r.get('targets',[]):
   acc=t.get('accession')
   if not acc:continue
   m[acc].append({'depth':0,'sign':1})
   front=[(acc,1,{acc})]
   for d in range(1,max_hops+1):
    nxt=[]
    for node,sgn,seen in front:
     for e in rev.get(node,[]):
      if e['source_acc'] in seen:continue
      ns=sgn*e['sign'];m[e['source_acc']].append({'depth':d,'sign':ns});nxt.append((e['source_acc'],ns,seen|{e['source_acc']}))
    front=nxt
  result[r['route_id']]=m
 return result

def chembl_paged(endpoint,params,field,max_records=30000):
 out=[];url='https://www.ebi.ac.uk/chembl/api/data/'+endpoint+'.json?'+urllib.parse.urlencode(dict(params,limit=1000),safe=',')
 while url:
  j=fetch_json(url,30,20_000_000);out.extend(j.get(field,[]))
  if len(out)>max_records:raise RuntimeError('V28_CHEMBL_RECORD_CAP')
  nxt=j.get('page_meta',{}).get('next');url=urllib.parse.urljoin('https://www.ebi.ac.uk',nxt) if nxt else None
 return out

def route_activities(routes,target_idx,cache_dir):
 cache=Path(cache_dir)/'route_activities.json'
 if cache.exists():return load(cache)
 tid_routes=defaultdict(list)
 for r in routes:
  for t in r.get('targets',[]):
   for tid in target_idx.get(t.get('accession'),[]):tid_routes[tid].append(r['route_id'])
 acts=[];assays={}
 tids=sorted(tid_routes)
 def _target(tid):
  cp=Path(cache_dir)/f'act_{tid}.json'
  if cp.exists():return tid,load(cp)
  try:
   rr=chembl_paged('activity',{'target_chembl_id':tid,'standard_type__in':'IC50,EC50,AC50,Ki,Kd','standard_units':'nM','standard_value__lte':'10000'},'activities')
  except Exception as exc:
   raise RuntimeError('V286_CHEMBL_ACTIVITY_UNAVAILABLE:'+str(tid)) from exc
  save(cp,rr);return tid,rr
 # Bounded parallelism only changes I/O scheduling; processing below remains
 # deterministic in sorted target order and every raw target response is cached.
 with concurrent.futures.ThreadPoolExecutor(max_workers=min(6,max(1,len(tids)))) as ex:
  target_rows=list(ex.map(_target,tids))
 for tid,rr in target_rows:
  for a in rr:
   try:v=float(a.get('standard_value'))
   except:continue
   if a.get('standard_relation')!='=' or a.get('standard_units')!='nM' or not 0<v<=10000:continue
   acts.append({k:a.get(k) for k in ['activity_id','assay_chembl_id','molecule_chembl_id','parent_molecule_chembl_id','target_chembl_id','standard_type','standard_value','canonical_smiles','document_chembl_id','data_validity_comment']})
 aids=sorted({a['assay_chembl_id'] for a in acts if a.get('assay_chembl_id')})
 batches=[aids[off:off+100] for off in range(0,len(aids),100)]
 def _assay_batch(ids):
  try:return chembl_paged('assay',{'assay_chembl_id__in':','.join(ids)},'assays',5000)
  except Exception as exc:raise RuntimeError('V286_CHEMBL_ASSAY_UNAVAILABLE') from exc
 if batches:
  with concurrent.futures.ThreadPoolExecutor(max_workers=min(6,len(batches))) as ex:
   assay_batches=list(ex.map(_assay_batch,batches))
  for rr in assay_batches:
   for a in rr:assays[a['assay_chembl_id']]={'description':a.get('description'),'confidence_score':a.get('confidence_score'),'assay_type':a.get('assay_type')}
 route_by={r['route_id']:r for r in routes};rows=[]
 for a in acts:
  ass=assays.get(a.get('assay_chembl_id'),{})
  if ass.get('confidence_score') not in (9,'9') or ass.get('assay_type') not in ('B','F'):continue
  lo=(ass.get('description') or '').lower()
  # Binding Ki/Kd and generic EC50 are NOT signed activations.
  # A functional assay description must establish the pharmacological direction.
  typ=a.get('standard_type')
  if typ=='IC50' and ('agonist' not in lo and 'activate' not in lo):
   direction=-1
  elif typ in ('EC50','AC50') and any(k in lo for k in ('agonist','stimulat','activat','induc')):
   direction=1
  elif typ in ('EC50','AC50') and any(k in lo for k in ('inhibit','antagon','repress','suppress')):
   direction=-1
  else:direction=0
  direct_functional=(ass.get('assay_type')=='F' and direction in (-1,1) and
     a.get('standard_type') in ('IC50','EC50','AC50') and
     any(token in lo for token in ('activat','agonist','inhibit','antagon','stimulat','suppress')))
  for rid in tid_routes[a['target_chembl_id']]:
   desired=1 if route_by[rid]['desired_action']=='increase' else -1 if route_by[rid]['desired_action']=='decrease' else 0
   role='compatible' if direction and direction==desired else 'opposite' if direction and direction!=desired else 'unknown'
   rows.append({**a,'route_id':rid,'value_nM':float(a['standard_value']),'direction':direction,'role':role,
    'assay_description':ass.get('description'),'assay_type':ass.get('assay_type'),
    'direct_target_assay_verified':bool(direct_functional and role=='compatible'),
    'assay_target_confidence':ass.get('confidence_score')})
 save(cache,rows);return rows

def structures(parent_ids,cache_dir):
 cache=Path(cache_dir)/'structures.json'
 if cache.exists():return load(cache)
 out={}
 for off in range(0,len(parent_ids),50):
  ids=parent_ids[off:off+50]
  try:
   u='https://www.ebi.ac.uk/chembl/api/data/molecule.json?'+urllib.parse.urlencode({'molecule_chembl_id__in':','.join(ids),'limit':1000},safe=',')
   j=fetch_json(u,30,20_000_000)
  except Exception as exc:raise RuntimeError('V286_CHEMBL_STRUCTURE_UNAVAILABLE') from exc
  for x in j.get('molecules',[]):
   st=x.get('molecule_structures') or {};out[x['molecule_chembl_id']]={'smiles':st.get('canonical_smiles'),'molecule_type':x.get('molecule_type')}
 save(cache,out);return out
def fp(s):
 if not s:return None
 try:m=Chem.MolFromSmiles(s)
 except:return None
 return FPGEN.GetFingerprint(m) if m else None

def _lib(name,lib,cache_dir):
 p=Path(cache_dir)/f'{name}.txt'
 if not p.exists():
  u='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName='+urllib.parse.quote(lib)
  p.write_bytes(fetch_bytes(u,45,30_000_000))
 return p
def _compound(term,key):
 if key.startswith('LINCS'):
  m=re.search(r'\b(?:\d+(?:\.\d+)?H)-(.+?)-\d+(?:\.\d+)?$',term,re.I);return m.group(1) if m else None
 m=re.match(r'^(.*?)\s+(?:homo sapiens|mus musculus|rattus norvegicus)\s+',term,re.I);return m.group(1).strip() if m else None

def transcriptomics(routes,alias,cache_dir):
 # Maayan Lab Enrichr web service/libraries require a separate commercial licence.
 # Keep original research analysis unchanged unless the caller explicitly disables this channel.
 mode=os.environ.get('RDMVA_ENRICHR_MODE','enabled').strip().lower()
 if mode not in ('enabled','disabled'):
  raise RuntimeError('RDMVA_ENRICHR_MODE must be enabled or disabled')
 commercial=os.environ.get('RDMVA_COMMERCIAL_MODE','0').strip().lower() in ('1','true','yes')
 if commercial and mode!='disabled':
  raise RuntimeError('ENRICHR_COMMERCIAL_PERMISSION_REQUIRED: use --commercial-mode to disable the Enrichr channel, or obtain and document permission')
 if mode=='disabled':
  return {'rows':[],'coverage':{},'license_mode':'disabled','provider':'Enrichr_LINCS_GEO','limitations':['Enrichr signatures excluded from this analysis; screening results can differ']}
 cache=Path(cache_dir)/'transcriptomics.json'
 if cache.exists():return load(cache)
 needs=defaultdict(list)
 for r in routes:
  d='up' if r['desired_action']=='increase' else 'down' if r['desired_action']=='decrease' else None
  if not d:continue
  for t in r.get('targets',[]):
   if t.get('gene'):needs[(t['gene'].upper(),d)].append(r['route_id'])
 sig=defaultdict(list);coverage=Counter()
 for key,(lib,dd) in LIBS.items():
  with _lib(key,lib,cache_dir).open(encoding='utf-8',errors='replace') as f:
   for line in f:
    p=line.rstrip('\n').split('\t')
    if len(p)<3:continue
    term=p[0];comp=_compound(term,key);labs=sorted(alias.get(norm_token(comp),set())) if comp else []
    if not labs:continue
    for lab in labs:coverage[lab]+=1
    direction=dd
    if key=='GEO':
     lo=term.lower();direction='down' if ' chdir down' in lo else 'up' if ' chdir up' in lo else None
     if not direction:continue
    genes={x.upper() for x in p[2:] if x}
    for (g,d),rids in needs.items():
     if d!=direction or g not in genes:continue
     for lab in labs:
      for rid in rids:sig[(lab,rid)].append({'library':lib,'term':term,'gene':g,'direction':direction})
 rows=[{'ingredient':lab,'route_id':rid,'support_count':len(xs),'examples':xs[:6]} for (lab,rid),xs in sig.items()]
 obj={'rows':rows,'coverage':dict(coverage)};save(cache,obj);return obj

def universal_screen(entities,routes,outdir):
 # Held causal hypotheses never contribute positive signal evidence.
 held_routes=sum(r.get('state')=='hold' for r in routes)
 routes=[r for r in routes if r.get('state')!='hold']
 out=Path(outdir);cache=out/'sources';cache.mkdir(parents=True,exist_ok=True)
 parent_labels,alias=raw_identity_maps(entities);target_idx=target_index(cache,routes);mech=mechanisms();paths=signed_reverse(3,routes)
 acts=route_activities(routes,target_idx,cache);parents=sorted({p for e in entities for p in e.get('chembl_parent_ids',[])})
 structs=structures(parents,cache);tx=transcriptomics(routes,alias,cache)
 route_by={r['route_id']:r for r in routes}
 act_by=defaultdict(list)
 for a in acts:
  pid=a.get('parent_molecule_chembl_id') or a.get('molecule_chembl_id')
  for lab in parent_labels.get(pid,set()):act_by[lab].append(a)
 tx_by=defaultdict(list)
 for x in tx['rows']:tx_by[x['ingredient']].append(x)
 # structural refs
 refs=defaultdict(list)
 for a in acts:
  if a['role']!='compatible' or a['value_nM']>1000 or not a.get('canonical_smiles'):continue
  z=fp(a['canonical_smiles'])
  if z is not None:refs[a['route_id']].append((z,a))
 for rid in list(refs):refs[rid]=sorted(refs[rid],key=lambda z:z[1]['value_nM'])[:250]
 ledger=[];signalq=[]
 for e in entities:
  lab=e['canonical_label'];signals=[];audit=[];pids=e.get('chembl_parent_ids',[])
  mrecs=[m for p in pids for m in mech.get(p,[])]
  for m in mrecs:
   for r in routes:
    desired=1 if r['desired_action']=='increase' else -1 if r['desired_action']=='decrease' else 0
    for p in paths.get(r['route_id'],{}).get(m['accession'],[]):
     if desired and m['action_sign']*p['sign']==desired:
      signals.append({'route_id':r['route_id'],'provider':'ChEMBL_MoA_SIGNOR','depth':p['depth'],'mechanism_id':m['mechanism_id'],'functional_rescue':False});break
  audit.append({'channel':'mechanism_network','attempted':bool(pids),'data_records':len(mrecs)})
  good=[a for a in act_by.get(lab,[]) if a['role']=='compatible' and a['value_nM']<=1000]
  for a in good:signals.append({'route_id':a['route_id'],'provider':'ChEMBL_measured_activity','value_nM':a['value_nM'],
    'assay_id':a['assay_chembl_id'],'document_chembl_id':a.get('document_chembl_id'),
    'direct_target_assay_verified':a.get('direct_target_assay_verified',False),
    'assay_description':a.get('assay_description'),'assay_type':a.get('assay_type'),
    'functional_rescue':False})
  audit.append({'channel':'measured_activity','attempted':bool(pids),'data_records':len(act_by.get(lab,[]))})
  for x in tx_by.get(lab,[]):signals.append({'route_id':x['route_id'],'provider':'perturbational_transcriptomics','support_count':x['support_count'],'functional_rescue':False})
  txcov=int(tx['coverage'].get(lab,0));audit.append({'channel':'transcriptomics','attempted':txcov>0,'data_records':txcov})
  has_struct=False
  if len(pids)==1:
   sm=(structs.get(pids[0]) or {}).get('smiles');z=fp(sm);has_struct=z is not None
   if z is not None:
    for r in routes:
     rr=refs.get(r['route_id'],[])
     if not rr:continue
     sims=DataStructs.BulkTanimotoSimilarity(z,[q[0] for q in rr]);i=max(range(len(sims)),key=sims.__getitem__);sim=float(sims[i])
     if sim>=0.60:signals.append({'route_id':r['route_id'],'provider':'rdkit_morgan_similarity','tanimoto':round(sim,4),'reference_molecule':rr[i][1].get('parent_molecule_chembl_id'),'functional_rescue':False})
  audit.append({'channel':'rdkit_structure','attempted':has_struct,'data_records':1 if has_struct else 0})
  seen=set();uniq=[]
  for s in signals:
   k=(s['route_id'],s['provider'],s.get('mechanism_id'),s.get('assay_id'),s.get('reference_molecule'))
   if k not in seen:seen.add(k);uniq.append(s)
  signals=uniq
  if signals:status='SIGNAL'
  elif e['modality'] in SPECIAL:status='INCOMPATIBLE_MODALITY'
  elif any(a['attempted'] for a in audit):status='NO_SIGNAL'
  else:status='INSUFFICIENT_DATA'
  row={**e,'screened':True,'status':status,'signals':signals,'signal_count':len(signals),'channels':audit,'patient_data_external':False}
  ledger.append(row)
  if len(ledger)%500==0:
   print('V31_APPROVED_SCREEN_PROGRESS',len(ledger),'OF',len(entities),flush=True)
  if status=='SIGNAL':signalq.append(row)
 if len(ledger)!=len(entities):raise RuntimeError('V28_SCREEN_COVERAGE_FAILURE')
 jsonl(out/'universal_ledger.jsonl',ledger);save(out/'signal_queue.json',signalq)
 summary={'scientific_held_routes_excluded':held_routes,
  'enrichr_library_mode':os.environ.get('RDMVA_ENRICHR_MODE','enabled').strip().lower(),'screened':len(ledger),'status_counts':dict(Counter(x['status'] for x in ledger)),'signal_entities':len(signalq),
  'signal_providers':dict(Counter(s['provider'] for x in ledger for s in x['signals'])),'route_count':len(routes),'no_entity_skipped':True}
 save(out/'screen_summary.json',summary);return ledger,summary
