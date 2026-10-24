"""Additional official marketing-authorisation registers, strictly source-bound.

Rare Disease MVA Pharma — V2.8.6 (separate research development line)
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch

Never infer approval in a different jurisdiction or individual-component
authorisation from a combination product. No raw patient data are accessed.
"""
from collections import defaultdict,Counter
from pathlib import Path
import datetime,hashlib,json,os,re,zipfile,xml.etree.ElementTree as ET
from .core import ROOT,RAW,norm_name

STORE=RAW/'regulatory_v283'
MANIFEST=(ROOT/'public_inputs/current_v31/regulatory_snapshot_v283_current.json'
  if os.environ.get('RDMVA_CURRENT_PUBLIC_DATA')=='1'
  else ROOT/'config/regulatory_snapshot_v283.json')
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
LOGIC_VERSION='V283_REGULATORY_SCHEMA_1'

def pinned_sources():
 if not MANIFEST.is_file():raise RuntimeError('V283_REGULATORY_MANIFEST_MISSING')
 spec=json.loads(MANIFEST.read_text(encoding='utf-8'))
 if spec.get('schema')!=LOGIC_VERSION:raise RuntimeError('V283_REGULATORY_MANIFEST_SCHEMA')
 files=spec['files']
 for name,item in files.items():
  p=STORE/name
  if not p.is_file():raise RuntimeError('V283_REGULATORY_FILE_MISSING:'+name)
  if p.stat().st_size!=item['bytes']:raise RuntimeError('V283_REGULATORY_SIZE_MISMATCH:'+name)
  h=hashlib.sha256()
  with p.open('rb') as f:
   for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
  if h.hexdigest()!=item['sha256']:raise RuntimeError('V283_REGULATORY_SHA256_MISMATCH:'+name)
 return files

def _letters_to_index(ref):
 col=re.match(r'[A-Z]+',ref).group();j=0
 for x in col:j=j*26+ord(x)-64
 return j-1

def xlsx_rows(path):
 """Yield (Excel row number, {zero-based col index: string}) from OOXML, stdlib only."""
 with zipfile.ZipFile(path) as z:
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   with z.open('xl/sharedStrings.xml') as f:
    for _,e in ET.iterparse(f,events=('end',)):
     if e.tag==NS+'si':
      strings.append(''.join(q.text or '' for q in e.iter(NS+'t')))
      e.clear()
  sheet='xl/worksheets/sheet1.xml'
  if sheet not in z.namelist():raise RuntimeError('V283_XLSX_SHEET_MISSING:'+path.name)
  with z.open(sheet) as f:
   for _,el in ET.iterparse(f,events=('end',)):
    if el.tag!=NS+'row':continue
    record={}
    for c in el.findall(NS+'c'):
     ref=c.get('r')
     if not ref:continue
     v=c.find(NS+'v');inline=c.find(NS+'is')
     if v is None and inline is None:continue
     value=v.text if v is not None else ''.join(x.text or '' for x in inline.iter(NS+'t'))
     if c.get('t')=='s' and value is not None:
      value=strings[int(value)]
     if value is not None and str(value).strip():
      record[_letters_to_index(ref)]=str(value).strip()
    idx=int(el.get('r','0'));el.clear()
    if record:yield idx,record

def _clean_one(s):
 """Strict mono-ingredient transcription. Do not guess at ambiguous mixtures."""
 s=' '.join(str(s or '').split()).strip()
 if not s or len(s)>115:return None
 if re.search(r'(?i)(?:\||\+|\s+AND\s+|\s+WITH\s+|\s*;\s*|\s*/\s*)',s):return None
 if re.search(r'(?i)\b(?:extract|tincture|homeopathic|anthroposoph|herbal mix|essential oil|mixture|mother tincture|technical grade)\b',s):return None
 if s.count(',')>=2:return None
 return norm_name(s)

def _record(source,country,approval,brand,parts,reg_id,raw=None,route=None):
 names=[];seen=set()
 for s in parts or []:
  n=_clean_one(s)
  if n is None:
   names=[]
   break  # Fail closed: never drop one ingredient of a combination.
  if n not in seen:seen.add(n);names.append(n)
 if not names:state='identity_unresolved'
 elif len(names)>1:state='components_explicit'
 else:state='single_substance'
 return {'product_id':source+':'+reg_id,'source':source,'country':country,
   'approval_status':approval,'brand':str(brand or ''),'source_active_ingredient_literal':str(raw or ''),
   'route':route,'component_parse':{'components':names,'combination':len(names)>1,
      'component_count':len(names),'identity_state':state},
   'approval_scope':'product_in_jurisdiction_only',
   'source_record_identifier':str(reg_id)}

def _rid(*parts):
 return hashlib.sha256('\x1f'.join(str(x) for x in parts).encode('utf-8')).hexdigest()[:20]

def _parse_ema57():
 fname='ema_article57_20261010.xlsx';counts=Counter()
 for row,cols in xlsx_rows(STORE/fname):
  if row<=20:continue
  brand=cols.get(0,'');raw=cols.get(1,'');country=cols.get(3,'')
  if not brand or not raw or not country:counts['missing_required']+=1;continue
  counts['all_valid_authorisation_rows']+=1
  # Article 57 supplies a free-text field without a structured component count.
  # Comma-separated formulations must never be interpreted as one medicine.
  parts=[] if ',' in raw else [raw]
  record=_record('EMA_ARTICLE57',country,'authorised_on_snapshot',brand,parts,
                 _rid(brand,raw,country,cols.get(2,''),cols.get(4,'')),raw,cols.get(2))
  record['source_row']=row
  record['source_record_identifier']=None
  if record['component_parse']['components']:counts['identity_accepted']+=1
  else:counts['identity_unresolved']+=1
  yield record

def _parse_ema_central():
 for row,cols in xlsx_rows(STORE/'ema_central_20261010.xlsx'):
  if row<=9 or cols.get(0,'').strip().lower()!='human':continue
  if cols.get(3,'').strip().lower()!='authorised':continue
  name=cols.get(1,'');pid=cols.get(2,'')
  raw=cols.get(7) or cols.get(6,'')
  if not name or not pid or not raw:continue
  parts=[] if ',' in raw else [raw]  # Cannot infer a single active from a free-text mixture
  entry=_record('EMA_CENTRALISED','European Union','authorised_on_snapshot',
    name,parts,_rid(pid,name),raw)
  entry['source_record_identifier']=pid
  entry['source_url']=cols.get(38)
  yield entry

def _parse_mhra():
 for n in ('1','2'):
  fname='uk_mhra_cat'+n+'_20261010.xlsx';group={}
  for row,c in xlsx_rows(STORE/fname):
   if row<=1:continue
   licence=c.get(0,'').strip();brand=c.get(2,'').strip();ingredient=c.get(3,'').strip()
   if not licence or not brand or not ingredient:continue
   key=(licence,brand)
   if key not in group:group[key]=[]
   group[key].append(ingredient)
  for (licence,brand),parts in sorted(group.items()):
   entry=_record('MHRA_CAT'+n,'United Kingdom','listed_in_published_uk_category_register',
      brand,parts,_rid(licence,brand),'; '.join(parts))
   entry['source_record_identifier']=licence
   yield entry

def _parse_canada():
 approved=json.loads((STORE/'ca_approved_20261010.json').read_text(encoding='utf-8'))
 marketed=json.loads((STORE/'ca_marketed_20261010.json').read_text(encoding='utf-8'))
 active={}
 for status,records in (('approved_not_marketed',approved),('marketed',marketed)):
  for r in records:
   if str(r.get('class_name','')).lower()!='human':continue
   code=int(r['drug_code'])
   if code not in active or status=='marketed':active[code]=(status,r)
 ingredient_map=defaultdict(list)
 for rec in json.loads((STORE/'ca_ingredients_20261010.json').read_text(encoding='utf-8')):
  code=int(rec['drug_code'])
  if code in active:ingredient_map[code].append(rec.get('ingredient_name',''))
 for code,(status,rec) in sorted(active.items()):
  ai=ingredient_map.get(code,[])
  listed=int(rec.get('number_of_ais') or 0)
  # Never present incomplete or discrepant component lists as single-agent approval.
  if listed<=0 or listed!=len(ai):parts=[]
  else:parts=ai
  entry=_record('HEALTH_CANADA_DPD','Canada',status,rec.get('brand_name'),parts,
    str(code),'; '.join(ai))
  entry['source_record_identifier']=rec.get('drug_identification_number')
  entry['source_drug_code']=code
  yield entry

def load_global_records():
 spec=pinned_sources()
 bysource={'EMA_ARTICLE57':_parse_ema57,'EMA_CENTRALISED':_parse_ema_central,
           'MHRA':_parse_mhra,'HEALTH_CANADA_DPD':_parse_canada}
 rows=[];statistics={}
 for name,fn in bysource.items():
  rc=list(fn())
  statuses=Counter(x['component_parse']['identity_state'] for x in rc)
  statistics[name]={'products':len(rc),'identity_states':dict(statuses),
   'active_component_mentions':sum(len(x['component_parse']['components']) for x in rc),
   'jurisdictions':len({x['country'] for x in rc})}
  rows.extend(rc)
 return rows,{'schema':LOGIC_VERSION,'source_statistics':statistics,
              'snapshot_sha256':{k:v['sha256'] for k,v in spec.items()},
              'total_additional_product_records':len(rows)}
