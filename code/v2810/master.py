"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
from collections import defaultdict,Counter
import zipfile,csv,io,re,json,os
from .core import *
from . import regulatory

FDAZIP=RAW/'fda/drugsatfda.zip'
PURPLE=RAW/'fda/purplebook_2026-09.csv'
CHEMBL=RAW/'chembl/catalog_phase4.json'

INVERTED=[
 r'^ESTROGENS, CONJUGATED(?: SYNTHETIC [AB])?$',
 r'^ACETIC ACID, GLACIAL$',r'^NITROFURANTOIN, MACROCRYSTALLINE$',
 r'^GONADOTROPIN, CHORIONIC$',r'^KRYPTON, KR-81M$',r'^UREA, C-14$',
 r'^FERROUS CITRATE, FE-59$'
]
def fda_rows(z,name):
 raw=z.read(name)
 try:text=raw.decode('utf-8-sig')
 except UnicodeDecodeError:text=raw.decode('cp1252')
 return list(csv.DictReader(io.StringIO(text,newline=''),delimiter='\t'))
def fda_status(statuses,has_orig):
 if not has_orig:return 'approval_history_not_confirmed'
 if any(x in {'1','2'} for x in statuses):return 'active_rx_otc_with_orig_approval'
 if '4' in statuses:return 'tentative_not_sufficient'
 if '3' in statuses:return 'discontinued'
 return 'status_unknown'
def _inv(s):return any(re.match(p,norm_name(s)) for p in INVERTED)
def _split(seg,strength=''):
 s=norm_name(seg)
 if not s:return []
 if _inv(s):return [s]
 if ' AND ' in s:
  p=[norm_name(re.sub(r'^AND\s+','',x,flags=re.I)) for x in re.split(r'\s*,\s*|\s+AND\s+',s)]
  p=[x for x in p if x]
  if len(p)>1:return p
 cp=[norm_name(x) for x in re.split(r'\s*,\s*',s) if norm_name(x)]
 sp=[x.strip() for x in re.split(r'\s*,\s*',strength or '') if x.strip()]
 return cp if len(cp)>1 and len(cp)==len(sp) else [s]
def parse_components(active,strength=''):
 isegs=(active or '').split(';');ssegs=(strength or '').split(';');flat=[];groups=[]
 for i,seg in enumerate(isegs):
  comps=_split(seg,ssegs[i] if i<len(ssegs) else '');flat+=comps
  groups.append({'group_index':i,'literal':seg,'components':comps})
 seen=set();uniq=[]
 for x in flat:
  if x and x not in seen:seen.add(x);uniq.append(x)
 return {'literal':active,'components':uniq,'groups':groups,'combination':len(uniq)>1,'component_count':len(uniq)}
def purple_rows():
 lines=PURPLE.read_text(encoding='utf-8-sig',errors='replace').splitlines()
 heads=[i for i,l in enumerate(lines) if l.startswith('N/R/U,Applicant,BLA Number')]
 if len(heads)<2:raise RuntimeError('purple_header_missing')
 return list(csv.DictReader(io.StringIO('\n'.join(lines[heads[1]:]))))
def biologic_modality(name):
 n=norm_name(name)
 if 'VACCINE' in n:return 'vaccine'
 if 'IMMUNE GLOBULIN' in n or 'IMMUNOGLOBULIN' in n:return 'immune_globulin'
 if any(x in n for x in ('GENE THERAPY','VECTOR','VIRUS')) or re.search(r'GENE\b',n):return 'gene_or_viral_therapy_review'
 if any(x in n for x in ('CELL','T-CELL','T CELL')):return 'cell_therapy_review'
 if any(x in n for x in ('ANTIHEMOPHILIC FACTOR','COAGULATION FACTOR','ALBUMIN','PLASMA')):return 'blood_or_coagulation_product'
 return 'biologic_or_therapeutic_protein'

def chembl_alias():
 alias=defaultdict(set);meta={}
 for x in load(CHEMBL):
  h=x.get('molecule_hierarchy') or {};pid=h.get('parent_chembl_id') or x['molecule_chembl_id']
  meta.setdefault(pid,{'parent_id':pid,'molecule_types':set(),'names':set()})
  if x.get('molecule_type'):meta[pid]['molecule_types'].add(x['molecule_type'])
  vals=[x.get('pref_name')]
  vals += [z.get('molecule_synonym') for z in x.get('molecule_synonyms',[])]
  vals += [z.get('synonyms') for z in x.get('molecule_synonyms',[])]
  for a in vals:
   n=norm_name(a)
   if n:alias[n].add(pid);meta[pid]['names'].add(n)
 return alias,meta

def current_dataset_profile():
 # The environment flag is set only by the official current dataset launcher.
 # Fail closed on an incomplete/missing current receipt; never assume historic
 # scientific totals remain constant across updated official releases.
 return os.environ.get('RDMVA_CURRENT_PUBLIC_DATA')=='1'

def verify_universe_profile(us_count):
 if current_dataset_profile():
  receipt=ROOT/'public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json'
  manifest=ROOT/'public_inputs/current_v31/public_data_manifest_v28_current.json'
  if not receipt.is_file() or not manifest.is_file():
   raise RuntimeError('V31_CURRENT_DATASET_RECEIPT_MISSING')
  rec=load(receipt)
  if rec.get('state')!='READY' or rec.get('profile')!='current_official_sources':
   raise RuntimeError('V31_CURRENT_DATASET_NOT_READY')
  static={x['path']:x for x in load(manifest)['static_files']}
  for name,p in (('fda/drugsatfda.zip',FDAZIP),('fda/purplebook_2026-09.csv',PURPLE)):
   expected=static[name]
   if p.stat().st_size!=expected['bytes'] or sha(p)!=expected['sha256']:
    raise RuntimeError('V31_CURRENT_US_SOURCE_HASH_MISMATCH:'+name)
  # Defensive sanity threshold, NOT an asserted release-dependent count.
  # Actual count and all source SHA256s are persisted in master/summary.json.
  if us_count<2000:
   raise RuntimeError('V31_CURRENT_US_UNIVERSE_UNEXPECTEDLY_SMALL')
  return 'CURRENT_OFFICIAL_SOURCES'
 if us_count!=2414:
  raise RuntimeError('V283_BASELINE_US_DRIFT')
 return 'FROZEN_2026_10_09'

def build(outdir):
 out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
 alias,chemmeta=chembl_alias()
 with zipfile.ZipFile(FDAZIP) as z:
  prod=fda_rows(z,'Products.txt');states=fda_rows(z,'MarketingStatus.txt');subs=fda_rows(z,'Submissions.txt')
 ms=defaultdict(set)
 for x in states:ms[(x['ApplNo'],x['ProductNo'])].add(x['MarketingStatusID'])
 orig={x['ApplNo'] for x in subs if x['SubmissionType']=='ORIG' and x['SubmissionStatus']=='AP'}
 contexts=defaultdict(lambda:{'solo':set(),'combo':set(),'purple':set()});products=[]
 for x in sorted(prod,key=lambda z:(z['ApplNo'],z['ProductNo'])):
  key=(x['ApplNo'],x['ProductNo'])
  if fda_status(ms[key],x['ApplNo'] in orig)!='active_rx_otc_with_orig_approval':continue
  px=parse_components(x.get('ActiveIngredient'),x.get('Strength'))
  if not px['components']:continue
  pid='DAFDA:'+':'.join(key)
  products.append({'product_id':pid,'source':'Drugs@FDA','brand':x.get('DrugName'),'form':x.get('Form'),'strength':x.get('Strength'),
    'active_ingredient_literal':x.get('ActiveIngredient'),'component_parse':px,'approval_scope':'product_only'})
  for c in px['components']:(contexts[c]['combo'] if px['combination'] else contexts[c]['solo']).add(pid)
 for x in purple_rows():
  if (x.get('Licensure') or '').strip().lower()!='licensed':continue
  if (x.get('Marketing Status') or '').strip().lower() in ('disc','discontinued'):continue
  proper=norm_name(x.get('Proper Name'))
  if not proper:continue
  pid='PB:'+str(x.get('BLA Number','')).strip()+':'+str(x.get('Product Number','')).strip()
  products.append({'product_id':pid,'source':'PurpleBook','brand':x.get('Proprietary Name'),'proper_name':x.get('Proper Name'),
    'component_parse':{'components':[proper],'combination':False,'component_count':1},'approval_scope':'licensed_biologic_product_only'})
  contexts[proper]['purple'].add(pid)
 # The 2026-10-09 US universe stays strictly frozen, while an explicitly
 # verified current-source snapshot records its own observed ingredient count.
 baseline_labels=set(contexts)
 dataset_profile=verify_universe_profile(len(baseline_labels))
 nonus,regulatory_summary=regulatory.load_global_records()
 if (not nonus or
     not regulatory_summary.get('source_statistics') or
     any(v.get('products',0)==0 for v in regulatory_summary['source_statistics'].values())):
  raise RuntimeError('V31_GLOBAL_REGULATORY_SOURCES_INCOMPLETE')
 new_provenance=defaultdict(lambda:{'countries':set(),'sources':set(),'mono':set(),'combo':set()})
 for p in nonus:
  products.append(p)
  parsed=p['component_parse']
  for c in parsed['components']:
   kind='combo' if parsed['combination'] else 'mono'
   (contexts[c]['combo'] if parsed['combination'] else contexts[c]['solo']).add(p['product_id'])
   new_provenance[c][kind].add(p['product_id'])
   new_provenance[c]['countries'].add(p['country'])
   new_provenance[c]['sources'].add(p['source'])
 entities=[]
 for label,cx in sorted(contexts.items()):
  chem=sorted(alias.get(label,set()));purple=bool(cx['purple']);solo=bool(cx['solo']);co=bool(cx['combo'])
  types=sorted({t for p in chem for t in chemmeta.get(p,{}).get('molecule_types',set())})
  prov=new_provenance[label]
  us_solo=any(x.startswith('DAFDA:') for x in cx['solo'])
  us_combo=any(x.startswith('DAFDA:') for x in cx['combo'])
  juris=sorted(prov['countries']|({'United States'} if us_solo or us_combo or purple else set()))
  approved_sources=sorted(prov['sources']|({'Drugs@FDA'} if us_solo or us_combo else set())|({'FDA PurpleBook'} if purple else set()))
  entities.append({'entity_id':stable_id('V28USM:' if label in baseline_labels else 'V283GLOBAL:',label),
   'fda_solo_approved':us_solo,'fda_combination_approved':us_combo,
   'approved_jurisdictions':juris,'authorisation_sources':approved_sources,
   'non_us_mono_product_ids':sorted(prov['mono']),'non_us_combo_product_ids':sorted(prov['combo']),
   'global_approved_as_single':bool(cx['solo'] or cx['purple']),
   'global_combination_only':bool(cx['combo'] and not cx['solo'] and not purple),'canonical_label':label,
   'entity_kind':'biologic' if purple else 'drug_component',
   'modality':biologic_modality(label) if purple else 'small_or_other_drug_component',
   'solo_approved':solo,'combination_approved':co,'combo_only':co and not solo and not purple,'biologic_licensed':purple,
   'solo_product_ids':sorted(cx['solo']),'combination_product_ids':sorted(cx['combo']),'purple_product_ids':sorted(cx['purple']),
   'chembl_parent_ids':chem,'chembl_identity_state':'exact_normalized_alias' if len(chem)==1 else ('ambiguous_exact_alias' if len(chem)>1 else 'unlinked'),
   'chembl_molecule_types':types,'standalone_approval_inferred_from_combination':False})
 jsonl(out/'products.jsonl',products);jsonl(out/'entities.jsonl',entities)
 summary={'version':'V283_MASTER_1','dataset_profile':dataset_profile,'us_snapshot_entity_count':len(baseline_labels),'baseline_us_entities':len(baseline_labels),'new_entity_labels':len(entities)-len(baseline_labels),'regulatory_extensions':regulatory_summary,'products':len(products),'entities':len(entities),'solo':sum(x['solo_approved'] for x in entities),
  'combo_only':sum(x['combo_only'] for x in entities),'biologic':sum(x['biologic_licensed'] for x in entities),
  'chembl_exact':sum(x['chembl_identity_state']=='exact_normalized_alias' for x in entities),
  'source_hashes':{'drugsatfda':sha(FDAZIP),'purplebook':sha(PURPLE),'chembl_catalogue':sha(CHEMBL)}}
 save(out/'summary.json',summary);return entities,products,summary
