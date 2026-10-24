"""Deterministic provenance for drug authorisation; original code and data remain separate.

Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""
import json
from collections import defaultdict
from .core import norm_name

REGISTRIES={
 'Drugs@FDA':('FDA Drugs@FDA','United States','https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files'),
 'PurpleBook':('FDA Purple Book','United States','https://purplebooksearch.fda.gov/'),
 'EMA_ARTICLE57':('EMA Article 57','EEA-specific','https://www.ema.europa.eu/en/human-regulatory-overview/post-authorisation/data-medicines-iso-idmp-standards-post-authorisation/public-data-article-57-database'),
 'EMA_CENTRALISED':('EMA centralised medicines','European Union','https://www.ema.europa.eu/en/medicines/download-medicine-data'),
 'MHRA_CAT1':('MHRA Category 1','United Kingdom','https://www.gov.uk/government/publications/category-lists-following-implementation-of-the-windsor-framework/category-lists-of-products'),
 'MHRA_CAT2':('MHRA Category 2','United Kingdom','https://www.gov.uk/government/publications/category-lists-following-implementation-of-the-windsor-framework/category-lists-of-products'),
 'HEALTH_CANADA_DPD':('Health Canada DPD','Canada','https://health-products.canada.ca/dpd-bdpp/')
}
# Live profile is dated from the verified local regulatory snapshot manifest.
def snapshot_date():
  from . import regulatory
  try:
    return json.loads(regulatory.MANIFEST.read_text(encoding='utf-8')).get('snapshot_date','2026-10-10')
  except (OSError,ValueError):
    return '2026-10-10'
SNAPSHOT=snapshot_date()
def _key(p):
 return 'PurpleBook' if (p.get('product_id') or '').startswith('PB:') else p.get('source')
def _status(p,src):
 if src=='Drugs@FDA':return 'active FDA drug product in verified source'
 if src=='PurpleBook':return 'FDA-licensed biological product in verified source'
 return str(p.get('approval_status') or 'status_unknown')

def _ref(p,src):
 if src in ('Drugs@FDA','PurpleBook'):
  return p.get('source_record_identifier') or p.get('product_id')
 if src=='EMA_ARTICLE57':
  return 'Excel row '+str(p['source_row']) if p.get('source_row') else None
 return p.get('source_record_identifier')

def build_index(products):
 grouped=defaultdict(dict)
 for p in products:
  src=_key(p)
  if src not in REGISTRIES:continue
  parsed=p.get('component_parse') or {}
  if parsed.get('identity_state')=='identity_unresolved':continue
  ingredients=parsed.get('components') or []
  if not ingredients:continue
  rname,default_country,url=REGISTRIES[src]
  country=p.get('country') or default_country
  status=_status(p,src)
  kind='combination_component' if parsed.get('combination') else 'single_ingredient_product'
  for name in dict.fromkeys(norm_name(x) for x in ingredients if norm_name(x)):
   k=(src,country,status,kind)
   g=grouped[name].get(k)
   if g is None:
    g={'register':rname,'source_code':src,'jurisdiction':country,
       'status':status,'product_context':kind,'source_url':p.get('source_url') or url,
       'snapshot_date':SNAPSHOT,'source_product_count':0,'example_products':[]}
    grouped[name][k]=g
   g['source_product_count']+=1
   if len(g['example_products'])<2:
    sample={'brand':p.get('brand') or None,'original_reference':_ref(p,src)}
    if sample not in g['example_products']:g['example_products'].append(sample)
 out={}
 for name,buckets in grouped.items():
  rows=sorted(buckets.values(),key=lambda r:(r['register'],r['jurisdiction'],r['product_context'],r['status']))
  single=any(r['product_context']=='single_ingredient_product' for r in rows)
  combined=any(r['product_context']=='combination_component' for r in rows)
  out[name]={'query_name':name,'match_method':'exact_normalized_active_ingredient',
    'snapshot_date':SNAPSHOT,'records':rows,
    'registry_names':sorted({r['register'] for r in rows}),
    'jurisdictions':sorted({r['jurisdiction'] for r in rows}),
    'single_product_listed':single,'combination_component_listed':combined,
    'overall_state':'single_ingredient_product_listed' if single else 'combination_only' if combined else 'not_in_snapshot',
    'caveat':'Product/jurisdiction records, not evidence of treatment efficacy for rare disease.'}
 return out

def for_candidate(name,index):
 key=norm_name(name)
 if key in index:
  return index[key]
 return {'query_name':key,'records':[],'registry_names':[],'jurisdictions':[],'overall_state':'not_in_snapshot'}

def compact_for_ai(c):
 return {'state':c['overall_state'],'registers':c['registry_names'],
  'jurisdictions':c['jurisdictions'],
  'regulatory_facts':[{'register':r['register'],'jurisdiction':r['jurisdiction'],
   'status':r['status'],'context':r['product_context']} for r in c['records'][:48]]}

def text_lines(c):
 rows=c.get('records') or []
 if not rows:
  return ['Regulatory status: not identified in loaded snapshots; this does not prove absence of approval.']
 headline=('At least one single-ingredient registered product is recorded.'
   if c['overall_state']=='single_ingredient_product_listed'
   else 'Only a component of combination products is recorded; no standalone authorisation established.')
 lines=['Regulatory approvals and listings: '+headline]
 for r in rows:
  kind='standalone product' if r['product_context']=='single_ingredient_product' else 'combination component ONLY'
  sample=r['example_products'][0] if r['example_products'] else {}
  ref=sample.get('original_reference')
  evidence=(' | original reference: '+str(ref)) if ref else ''
  form=(' | documented form: '+r['specific_ingredient_form']) if r.get('specific_ingredient_form') else ''
  lines.append('  - '+r['register']+' ('+r['jurisdiction']+'): '+r['status']+'; '+kind+form+evidence+' | '+r['source_url'])
 lines.append('No approval for treatment of the rare disease is implied.')
 return lines
