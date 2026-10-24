"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
from collections import defaultdict,deque
import csv,re,urllib.parse
from .core import *

SIGNOR=RAW/'signor/Oct2026_release.txt'
UNIPROT_DIR=RAW/'uniprot'

def _gene_symbols(u):
 out=[]
 for g in u.get('genes',[]):
  x=g.get('geneName',{}).get('value')
  if x:out.append(x.upper())
 return out
def uniprot_record(gene,cache_dir):
 cache=Path(cache_dir)/f'uniprot_{gene}.json'
 if cache.exists():return load(cache)
 # Prefer local public input when exact gene file exists.
 matches=list(UNIPROT_DIR.glob(gene+'_*.json'))
 if matches:
  u=load(matches[0]);save(cache,u);return u
 q=f'(gene_exact:{gene}) AND (organism_id:9606)'
 url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':q,'format':'json','size':5})
 try:j=fetch_json(url,20,5_000_000)
 except Exception:j={'results':[]}
 for u in j.get('results',[]):
  if gene.upper() in _gene_symbols(u):save(cache,u);return u
 # Fallback for literature protein names/synonyms to canonical human gene records.
 # Very short non-gene labels such as "PP1" are too ambiguous for automatic pharmacologic routing.
 needle=norm_token(gene)
 if len(needle)<4:raise RuntimeError('UNIPROT_NAME_TOO_AMBIGUOUS_'+str(gene))
 q2=f'({gene}) AND (organism_id:9606)'
 url2='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':q2,'format':'json','size':10})
 j2=fetch_json(url2,20,8_000_000);needle=norm_token(gene);scored=[]
 for u in j2.get('results',[]):
  vals=[]
  vals += _gene_symbols(u)
  for g in u.get('genes',[]):vals += [s.get('value') for s in g.get('synonyms',[]) if s.get('value')]
  pd=u.get('proteinDescription') or {}
  x=((pd.get('recommendedName') or {}).get('fullName') or {}).get('value')
  if x:vals.append(x)
  for a in pd.get('alternativeNames',[]):
   y=(a.get('fullName') or {}).get('value')
   if y:vals.append(y)
   vals += [s.get('value') for s in a.get('shortNames',[]) if s.get('value')]
  ntoks=re.findall(r'[a-z]+|\d+',str(gene).lower())
  scores=[]
  for v in vals:
   nv=norm_token(v);vtoks=re.findall(r'[a-z]+|\d+',str(v).lower())
   s=1.0 if nv==needle else 0.7 if needle and needle in nv else 0.0
   if ntoks and all(tok in vtoks for tok in ntoks):s=max(s,0.9)
   scores.append(s)
  score=max(scores or [0])
  if score:scored.append((score,u))
 if scored:
  scored.sort(key=lambda z:-z[0]);best=scored[0]
  if len(scored)==1 or best[0]>scored[1][0]:save(cache,best[1]);return best[1]
 raise RuntimeError('UNIPROT_NAME_NOT_UNIQUE_'+gene)

def accession_for_gene(gene,cache_dir):
 u=uniprot_record(gene,cache_dir)
 return u.get('primaryAccession')

def canonical_gene_for_name(name,cache_dir):
 u=uniprot_record(name,cache_dir)
 genes=_gene_symbols(u)
 return genes[0] if genes else str(name).upper()

def search_aliases(gene,cache_dir):
 u=uniprot_record(gene,cache_dir);vals=[gene]
 for g in u.get('genes',[]):
  x=g.get('geneName',{}).get('value')
  if x:vals.append(x)
  vals += [s.get('value') for s in g.get('synonyms',[]) if s.get('value')]
 pd=u.get('proteinDescription') or {}
 rn=((pd.get('recommendedName') or {}).get('fullName') or {}).get('value')
 if rn:vals.append(rn)
 for a in pd.get('alternativeNames',[]):
  x=(a.get('fullName') or {}).get('value')
  if x:vals.append(x)
  vals += [s.get('value') for s in a.get('shortNames',[]) if s.get('value')]
 out=[];seen=set()
 for v in vals:
  k=str(v).strip().lower()
  if k and k not in seen:seen.add(k);out.append(str(v).strip())
 return out[:12]

def uniprot_context(gene,cache_dir):
 u=uniprot_record(gene,cache_dir);texts=[]
 for c in u.get('comments',[]):
  if c.get('commentType') in ('FUNCTION','SUBUNIT','PTM','DISEASE'):
   for t in c.get('texts',[]):texts.append({'kind':c['commentType'],'text':t.get('value','')})
 return {'gene':gene,'accession':u.get('primaryAccession'),'texts':texts}

def signor_edges():
 out=[]
 with SIGNOR.open(encoding='utf-8',errors='replace') as f:
  rd=csv.DictReader(f,delimiter='\t')
  for r in rd:
   if r.get('TYPEA')!='protein' or r.get('TYPEB')!='protein' or str(r.get('TAX_ID','')).strip()!='9606':continue
   if (r.get('DATABASEA') or '').upper()!='UNIPROT' or (r.get('DATABASEB') or '').upper()!='UNIPROT':continue
   eff=(r.get('EFFECT') or '').strip().strip('"').lower()
   sign=1 if eff.startswith('up-regulates') else -1 if eff.startswith('down-regulates') else 0
   if not sign:continue
   out.append({'source_gene':r.get('ENTITYA'),'source_acc':r.get('IDA'),'target_gene':r.get('ENTITYB'),'target_acc':r.get('IDB'),
    'sign':sign,'effect':r.get('EFFECT'),'mechanism':(r.get('MECHANISM') or '').strip().strip('"'),'direct':str(r.get('DIRECT','')).lower() in ('t','true','1','yes'),
    'pmid':r.get('PMID'),'signor_id':r.get('SIGNOR_ID'),'sentence':r.get('SENTENCE')})
 return out

def _route(r):
 r=dict(r);r.setdefault('state','active');r.setdefault('score',0.0);r.setdefault('targets',[])
 return r

def base_routes(case_private,cache_dir):
 edges=signor_edges();routes=[]
 genes=sorted(case_private['genes'])
 for gene in genes:
  start_index=len(routes)
  acc=accession_for_gene(gene,cache_dir)
  # direct restoration routes
  routes.append(_route({'family':'direct_gene','target_concept':gene,'targets':[{'gene':gene,'accession':acc}],
    'desired_action':'increase','variable':'abundance','mechanism_concept':f'restore functional {gene} abundance','score':16.0,
    'requires_residual_product':True,'source':'case_gene'}))
  routes.append(_route({'family':'direct_gene','target_concept':gene,'targets':[{'gene':gene,'accession':acc}],
    'desired_action':'increase','variable':'activity','mechanism_concept':f'support residual {gene} function','score':16.0,
    'requires_residual_product':True,'source':'case_gene'}))
  classes={x['class'] for x in case_private['variant_classes'] if x['gene']==gene}
  if 'nonsense' in classes or 'frameshift' in classes:
   if 'nonsense' in classes:
    routes.append(_route({'family':'variant_rescue','target_concept':f'{gene} premature-stop readthrough','targets':[],
     'desired_action':'increase','variable':'full_length_translation','mechanism_concept':'promote premature-termination-codon readthrough','score':15.0,
     'state':'conditional','source':'variant_class'}))
   routes.append(_route({'family':'variant_rescue','target_concept':'nonsense-mediated decay','targets':[],
     'desired_action':'decrease','variable':f'NMD_on_{gene}_PTC','mechanism_concept':'reduce transcript loss only if allele is expressed','score':3.0,
     'state':'hold','source':'variant_class'}))
  if 'missense_or_other_coding' in classes:
   routes.append(_route({'family':'proteostasis','target_concept':f'{gene} folding/stability','targets':[{'gene':gene,'accession':acc}],
     'desired_action':'increase','variable':'folding_stability','mechanism_concept':'increase folding/stability of residual missense protein','score':13.0,
     'state':'conditional','source':'variant_class'}))
  if 'splice' in classes:
   routes.append(_route({'family':'splicing','target_concept':f'{gene} productive splicing','targets':[{'gene':gene,'accession':acc}],
     'desired_action':'increase','variable':'productive_transcript','mechanism_concept':'restore productive splicing','score':14.0,
     'state':'conditional','source':'variant_class'}))
  # Direct signed upstream regulators: manipulate regulator to increase primary gene state.
  for e in edges:
   if e['target_acc']!=acc or not e['direct']:continue
   desired='increase' if e['sign']>0 else 'decrease'
   routes.append(_route({'family':'signed_upstream','target_concept':e['source_gene'],'targets':[{'gene':e['source_gene'],'accession':e['source_acc']}],
     'desired_action':desired,'variable':'activity','mechanism_concept':f'{desired} upstream regulator to support {gene}','score':12.0,
     'state':'conditional','source':e['signor_id'],'pmid':e['pmid']}))
  # Direct downstream compensation: restore the normal direction normally imposed by primary gene.
  for e in edges:
   if e['source_acc']!=acc or not e['direct']:continue
   desired='increase' if e['sign']>0 else 'decrease'
   routes.append(_route({'family':'downstream_compensation','target_concept':e['target_gene'],'targets':[{'gene':e['target_gene'],'accession':e['target_acc']}],
     'desired_action':desired,'variable':'activity','mechanism_concept':f'compensate downstream state normally driven by {gene}','score':9.0,
     'state':'conditional','source':e['signor_id'],'pmid':e['pmid']}))
  # Preserve patient disease-gene provenance for EVERY route (including regulators).
  for added in routes[start_index:]:added.setdefault('case_gene',gene)
 # dedupe
 d={}
 for r in routes:
  k=(r.get('case_gene'),r['family'],r['target_concept'],r['desired_action'],r['variable'])
  if k not in d or r['score']>d[k]['score']:d[k]=r
 out=sorted(d.values(),key=lambda x:(-x['score'],x['family'],x['target_concept']))
 for i,r in enumerate(out,1):r['route_id']=f'V28R{i:03d}';r['rank']=i
 return out

def apply_case_gates(route_rows,case_private):
 out=[];residual=case_private.get('residual_protein_status','unknown')
 for r in route_rows:
  x=dict(r);reason=[]
  quote=' '.join(e.get('quote','') for e in x.get('evidence',[])).lower()
  context=x.get('evidence_context')
  if context=='meiotic' or any(k in quote for k in ('oocyte','meiotic','meiosis')):
   x['state']='hold';reason.append('meiotic_to_somatic_transfer_unresolved')
  elif x.get('state') not in ('active','conditional','hold'):x['state']='conditional'
  if residual=='unknown':
   if x.get('family') in ('direct_gene','proteostasis') or x.get('family')=='literature_discovered':
    if x['state']!='hold':x['state']='conditional'
    reason.append('residual_primary_protein_unknown')
  targets={t.get('gene','').upper() for t in x.get('targets',[])}
  if x.get('desired_action')=='decrease' and targets & {'AURKB','PLK1','CDK1','CDC20'}:
   if x['state']=='active':x['state']='hold'
   reason.append('global_cell_cycle_modulation_context_dependent')
  x['case_gate_reasons']=sorted(set(reason));out.append(x)
 return out

def merge_discovered(base,discovered):
 out=list(base);seen={(x.get('case_gene'),x['target_concept'].upper(),x['desired_action'],x['variable']) for x in out}
 for x in discovered:
  key=(x.get('case_gene'),str(x.get('target_concept','')).upper(),x.get('desired_action'),x.get('variable'))
  if not key[0] or key in seen:continue
  targets=x.get('targets') or []
  if not targets and x.get('target_gene') and x.get('target_accession'):
   targets=[{'gene':x['target_gene'],'accession':x['target_accession']}]
  out.append(_route({'family':'literature_discovered','target_concept':x['target_concept'],'targets':targets,
    'desired_action':x['desired_action'],'variable':x['variable'],'mechanism_concept':x.get('mechanism_concept') or x['target_concept'],
    'score':float(x.get('score',11.0)),'state':x.get('state','conditional'),'source':x.get('source_id'),'evidence':x.get('evidence',[]),
    'route_direction_verification':x.get('route_direction_verification'),
    'evidence_context':x.get('evidence_context'),'case_gene':x.get('case_gene')}));seen.add(key)
 out.sort(key=lambda x:(-x['score'],x['family'],x['target_concept']))
 for i,r in enumerate(out,1):r['route_id']=f'V28R{i:03d}';r['rank']=i
 return out
