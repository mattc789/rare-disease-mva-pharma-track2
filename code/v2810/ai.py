"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
import os,sys,ctypes,sqlite3,decimal,json,re,urllib.request,urllib.error,urllib.parse,xml.etree.ElementTree as ET,datetime,concurrent.futures,math
from ctypes import wintypes
from .core import *
from . import registry_output

STATE=Path(os.environ.get('LOCALAPPDATA',Path.home()/'AppData/Local'))/'Examedix'/'RareDiseaseMVAPharma'
BUDGET_DB=STATE/'v2810_budget_final.sqlite3'
API_CACHE=STATE/'v2810_api_cache'
GRANT='v2810-research-ccby-autobootstrap-2026-10-11'
CAP_USD=decimal.Decimal('3.0')
MODELS={
 'scout':{'model':'gpt-5.4-mini-2026-03-17','input':decimal.Decimal('0.75'),'output':decimal.Decimal('4.50'),'reasoning':'none'},
 'verify':{'model':'gpt-5.6-terra','input':decimal.Decimal('2.00'),'output':decimal.Decimal('12.00'),'reasoning':'low'},
 'final':{'model':'gpt-5.6-sol','input':decimal.Decimal('4.00'),'output':decimal.Decimal('20.00'),'reasoning':'medium'},
}
ALLOWED_LICENSES={'CC-BY-4.0','CC-BY-3.0','CC-BY-2.0','CC0-1.0'}
MICRO=decimal.Decimal('1000000')

class Blob(ctypes.Structure):_fields_=[('cbData',wintypes.DWORD),('pbData',ctypes.POINTER(ctypes.c_ubyte))]
def _dpapi(data):
 buf=(ctypes.c_ubyte*len(data)).from_buffer_copy(data);src=Blob(len(data),buf);dst=Blob()
 crypt=ctypes.WinDLL('crypt32',use_last_error=True);kernel=ctypes.WinDLL('kernel32',use_last_error=True)
 fn=crypt.CryptUnprotectData;fn.argtypes=[ctypes.POINTER(Blob),ctypes.c_void_p,ctypes.POINTER(Blob),ctypes.c_void_p,ctypes.c_void_p,wintypes.DWORD,ctypes.POINTER(Blob)]
 if not fn(ctypes.byref(src),None,None,None,None,1,ctypes.byref(dst)):raise RuntimeError('V28_CREDENTIAL_UNAVAILABLE')
 kernel.LocalFree.argtypes=[ctypes.c_void_p]
 try:return ctypes.string_at(dst.pbData,dst.cbData)
 finally:kernel.LocalFree(dst.pbData)
def api_key():
 if 'OPENAI_API_KEY' in os.environ:
  k=os.environ['OPENAI_API_KEY']
 else:k=_dpapi((STATE/'openai_key.dpapi').read_bytes()).decode('ascii')
 if not k or any(c.isspace() for c in k):raise RuntimeError('V28_INVALID_API_KEY')
 return k

def _db():
 BUDGET_DB.parent.mkdir(parents=True,exist_ok=True);c=sqlite3.connect(BUDGET_DB)
 c.execute('CREATE TABLE IF NOT EXISTS grants(id TEXT PRIMARY KEY,cap INTEGER NOT NULL)')
 c.execute('CREATE TABLE IF NOT EXISTS req(hash TEXT PRIMARY KEY,reserve INTEGER,status TEXT,actual INTEGER,model TEXT)')
 cap=int(CAP_USD*MICRO);c.execute('INSERT OR IGNORE INTO grants VALUES(?,?)',(GRANT,cap));c.commit();return c
def _committed_microdollars(c):
 # Completed requests consume their measured token-cost estimate, not the much
 # larger conservative preflight reservation. Unresolved/unknown requests keep
 # the full reservation so transport uncertainty remains fail-safe.
 return c.execute("""SELECT COALESCE(SUM(
   CASE WHEN status='completed' AND actual IS NOT NULL THEN actual ELSE reserve END
 ),0) FROM req""").fetchone()[0]
def budget_status():
 with _db() as c:
  cap=c.execute('SELECT cap FROM grants WHERE id=?',(GRANT,)).fetchone()[0]
  r=c.execute('SELECT COALESCE(SUM(reserve),0),COALESCE(SUM(actual),0),COUNT(*) FROM req').fetchone()
  committed=_committed_microdollars(c)
 return {'cap_usd':cap/1e6,'reserved_usd':r[0]/1e6,'actual_estimate_usd':r[1]/1e6,
         'committed_usd':committed/1e6,'requests':r[2]}
def _reserve(h,usd,model):
 amt=int((usd*MICRO).to_integral_value(rounding=decimal.ROUND_CEILING))
 with _db() as c:
  if c.execute('SELECT 1 FROM req WHERE hash=?',(h,)).fetchone():raise RuntimeError('V28_DUPLICATE_REQUEST')
  cap=c.execute('SELECT cap FROM grants WHERE id=?',(GRANT,)).fetchone()[0];tot=_committed_microdollars(c)
  if tot+amt>cap:raise RuntimeError('V28_AI_BUDGET_HOLD')
  c.execute('INSERT INTO req VALUES(?,?,?,?,?)',(h,amt,'reserved',None,model))
 return amt
def _finish(h,status,actual):
 with _db() as c:c.execute('UPDATE req SET status=?,actual=? WHERE hash=?',(status,actual,h))
def _request(body,role,outfile):
 cfg=MODELS[role];outfile=Path(outfile)
 rawbody=dump(body).encode();h=hashlib.sha256(rawbody).hexdigest();cachefile=API_CACHE/(h+'.json')
 if outfile.exists():
  old=load(outfile)
  if old.get('request_sha256')==h:return old
  raise RuntimeError('V28_RESUME_HASH_MISMATCH')
 if cachefile.exists():
  rec=load(cachefile);replay=dict(rec);replay['replayed_from_global_cache']=True;replay['request_sha256']=h
  save(outfile,replay);return replay
 in_upper=2*len(rawbody)+16384;maxout=int(body['max_output_tokens'])
 reserve=(decimal.Decimal(in_upper)*cfg['input']+decimal.Decimal(maxout)*cfg['output'])/MICRO
 reserved=_reserve(h,reserve,cfg['model'])
 req=urllib.request.Request('https://api.openai.com/v1/responses',data=rawbody,method='POST',
  headers={'Authorization':'Bearer '+api_key(),'Content-Type':'application/json'})
 try:
  with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req,timeout=90) as r:raw=json.loads(r.read(2_000_000))
 except Exception:
  _finish(h,'transport_unknown',None);raise
 usage=raw.get('usage') or {};actual=int((decimal.Decimal(usage.get('input_tokens',0))*cfg['input']+decimal.Decimal(usage.get('output_tokens',0))*cfg['output']).to_integral_value(rounding=decimal.ROUND_CEILING))
 chunks=[]
 for o in raw.get('output',[]):
  for p in o.get('content',[]):
   if p.get('type')=='output_text':chunks.append(p.get('text',''))
 if raw.get('status')!='completed':
  _finish(h,'invalid_output',actual);raise RuntimeError('V28_AI_INCOMPLETE')
 try:value=json.loads(''.join(chunks))
 except:
  _finish(h,'invalid_output',actual);raise RuntimeError('V28_AI_INVALID_JSON')
 rec={'request_sha256':h,'model':raw.get('model'),'response_id':raw.get('id'),'usage':usage,'actual_microdollars':actual,'reserved_microdollars':reserved,
  'patient_context_sent':False,'store':False,'tools':[],'value':value}
 API_CACHE.mkdir(parents=True,exist_ok=True);save(cachefile,rec);save(outfile,rec);_finish(h,'completed',actual);return rec

def _txt(x):return ' '.join(''.join(x.itertext()).split()) if x is not None else ''
def _license(root):
 s=' '.join(_txt(x)+' '+' '.join(x.attrib.values()) for x in root.findall('.//permissions//*')).lower()
 if 'creativecommons.org/publicdomain/zero/1.0' in s:return 'CC0-1.0'
 for v in ('4.0','3.0','2.0'):
  if f'creativecommons.org/licenses/by/{v}' in s:return f'CC-BY-{v}'
 return None
def literature(gene,outdir,max_docs=8,candidate=None,targets=None):
 out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
 terms=[gene]
 if targets:
  short=[]
  for x in targets:
   s=str(x).strip()
   if s and len(s)<=24 and len(s.split())<=3:short.append(s)
  terms+=short[:6]
 # de-duplicate aliases case-insensitively
 seen_terms=set();terms=[x for x in terms if not (x.lower() in seen_terms or seen_terms.add(x.lower()))]
 tq=' OR '.join(f'"{x}"' if ' ' in x else x for x in terms)
 if candidate:
  cand_terms=[candidate] if isinstance(candidate,str) else [str(x) for x in candidate if str(x).strip()]
  cq=' OR '.join(f'"{x}"' if ' ' in x else x for x in dict.fromkeys(cand_terms))
  queries=[f'(({cq}) AND ({tq}) AND (rescue OR checkpoint OR chromosome OR stabiliz* OR degrad* OR ubiquit* OR functional OR toxicity)) AND OPEN_ACCESS:Y AND HAS_FT:Y']
 else:
  # Multiple orthogonal mechanism queries improve recall without any drug-name seed.
  queries=[
   f'(TITLE_ABS:({tq}) AND (ubiquitin* OR degrad* OR stability OR proteostasis)) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND ("E3 ligase" OR ubiquitin-ligase OR deubiquitin*)) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND (deacetylase OR deacetyl* OR acetylation OR acetyltransferase)) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND (chaperone OR folding OR proteostasis OR cochaperone)) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND (phosphatase OR localization OR trafficking OR recruitment OR "complex assembly")) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND (rescue OR suppressor OR compensation OR restore)) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND (pharmacolog* OR drug OR compound OR inhibitor OR activator)) AND OPEN_ACCESS:Y AND HAS_FT:Y',
   f'(TITLE_ABS:({tq}) AND (approved OR "FDA-approved" OR repurpos* OR reposition* OR "off-target") AND (drug OR compound OR inhibitor OR activator OR binding)) AND OPEN_ACCESS:Y AND HAS_FT:Y'
  ]
 # Deterministic public-literature cache. Reuse an exact prior bundle instead
 # of rewriting it; this also avoids transient Windows/AV locks on reruns.
 bundle=out/'literature_bundle.json'
 if bundle.is_file():
  try:
   cached=load(bundle)
   if (cached.get('gene')==gene and cached.get('candidate')==candidate and
       cached.get('queries')==queries and isinstance(cached.get('documents'),list)):
    return cached['documents'][:max_docs]
  except Exception:
   pass
 def _meta_score(r):
  title=(r.get('title') or '').lower();abstract=(r.get('abstractText') or '').lower();text=title+' '+abstract
  sc=0
  if any(t.lower() in title for t in terms):sc+=8
  if any(k in title for k in ('ubiquit','degrad','stabil','deacetyl','acetyl','chaperone','folding','phosphatase','localization','recruitment','rescue','suppress','compensation')):sc+=8
  if any(k in text for k in ('knockdown','overexpression','inhibition','inhibitor','deacetylation','polyubiquitination','stabilization','recruitment','localization','rescued','rescue','e3 ligase','phosphatase','chaperone')):sc+=5
  if any(k in title for k in ('bioinformatics','prognostic','signature','hub gene','machine learning')):sc-=8
  if 'review' in title or 'insights into' in title:sc-=3
  try:cites=int(r.get('citedByCount') or 0);sc+=min(8.0,1.5*math.log2(cites+1))
  except Exception:pass
  return sc
 def _search(q):
  u='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','pageSize':100,'resultType':'core'})
  try:res=fetch_json(u,20,10_000_000).get('resultList',{}).get('result',[])
  except Exception:return []
  return sorted(list(enumerate(res)),key=lambda z:(-_meta_score(z[1]),z[0]))
 with concurrent.futures.ThreadPoolExecutor(max_workers=min(4,len(queries))) as ex:
  scored=list(ex.map(_search,queries))
 query_lists=[(q,[r for original_rank,r in rows]) for q,rows in zip(queries,scored)]
 # Interleave result ranks across mechanism classes so one broad query cannot monopolize the bundle.
 pooled=[];seen_meta=set();mx=max([len(x[1]) for x in query_lists] or [0])
 for rank in range(mx):
  for q,res in query_lists:
   if rank>=len(res):continue
   r=res[rank];pmc=r.get('pmcid')
   if not pmc or pmc in seen_meta:continue
   seen_meta.add(pmc);pooled.append((rank,r,q))
 def _full(item):
  rank,r,q=item;pmc=r.get('pmcid')
  try:root=ET.fromstring(fetch_bytes(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',15,12_000_000))
  except Exception:return None
  lic=_license(root)
  if lic not in ALLOWED_LICENSES:return None
  passages=[];keywords=[x.lower() for x in terms]+['rescue','suppress','stabil','degrad','ubiquit','kinetochore','localiz','chaperone','fold','acetyl','phosphatase']
  allp=[('A',i,p) for i,p in enumerate(root.findall('.//abstract//p'),1)]+[('B',i,p) for i,p in enumerate(root.findall('.//body//p'),1)]
  for prefix,i,p in allp:
   t=_txt(p)
   if len(t)<80:continue
   lo=t.lower();score=sum(2 for k in keywords if k in lo)+sum(2 for k in ('knockdown','overexpression','inhibition','inhibitor','deacetylation','polyubiquitination','rescued','rescue') if k in lo)
   if any(term.lower() in lo for term in terms) and score>=6:passages.append((score,prefix,i,t))
  passages=sorted(sorted(passages,key=lambda z:(-z[0],z[1],z[2]))[:10],key=lambda z:(z[1],z[2]))
  if not passages:return None
  return {'source_id':pmc,'citation':_txt(root.find('.//article-title')) or pmc,'source_url':f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/',
   'license':lic,'search_query':q,'passages':[{'id':f'{prefix}{i}','text':t} for sc,prefix,i,t in passages]}
 # Take a diversified top slice from every mechanism query, then rank globally.
 candidate_map={}
 for q,res in query_lists:
  for rank,r in enumerate(res[:12]):
   pmc=r.get('pmcid')
   if not pmc:continue
   item=(-_meta_score(r),rank,r,q)
   if pmc not in candidate_map or item[0:2] < candidate_map[pmc][0:2]:candidate_map[pmc]=item
 candidates=[(rank,r,q) for neg,rank,r,q in sorted(candidate_map.values(),key=lambda z:(z[0],z[1],z[2].get('pmcid') or ''))[:60]]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  fetched=list(ex.map(_full,candidates))
 valid=[d for d in fetched if d];docs=[];seen_doc=set()
 if not candidate:
  # Reserve up to two admitted papers per mechanism class before global fill.
  for q in queries:
   for d in [x for x in valid if x.get('search_query')==q][:2]:
    if d['source_id'] not in seen_doc:docs.append(d);seen_doc.add(d['source_id'])
    if len(docs)>=max_docs:break
   if len(docs)>=max_docs:break
 for d in valid:
  if len(docs)>=max_docs:break
  if d['source_id'] in seen_doc:continue
  docs.append(d);seen_doc.add(d['source_id'])
 save(out/'literature_bundle.json',{'gene':gene,'candidate':candidate,'queries':queries,'documents':docs,'patient_data':False});return docs

def _route_schema():
 item={'type':'object','properties':{
  'target_concept':{'type':'string'},'target_gene':{'type':['string','null']},'desired_action':{'type':'string','enum':['increase','decrease']},
  'variable':{'type':'string','enum':['activity','abundance','expression','stability','degradation','folding','localization','complex_function','other']},'mechanism_concept':{'type':'string'},'state':{'type':'string','enum':['active','conditional','hold']},
  'score':{'type':'number'},'source_id':{'type':'string'},'passage_id':{'type':'string'},'quote':{'type':'string'}},
  'required':['target_concept','target_gene','desired_action','variable','mechanism_concept','state','score','source_id','passage_id','quote'],'additionalProperties':False}
 return {'type':'object','properties':{'routes':{'type':'array','items':item,'maxItems':12}},'required':['routes'],'additionalProperties':False}

ROUTE_DISCOVERY_INSTRUCTIONS="Use only the supplied public passages. Discover molecular rescue routes, not treatments.\nExtract a route ONLY for a causal molecular modifier explicitly manipulated or mechanistically shown to change the named gene/protein's abundance, stability, degradation, folding, localization, checkpoint/complex function, or a demonstrated functional rescue endpoint.\nThe target_gene field must be the MODIFIER gene/protein that one would pharmacologically increase or decrease, not merely the primary disease gene mentioned in the outcome. Do not return the primary gene itself as target_gene unless the passage describes a direct pharmacologic modifier of that primary protein.\nReject correlation-only, co-expression, prognostic association, generic cancer signatures, and mere pathway membership. Reject passages that only say the primary gene level was reduced without identifying a causal modifier.\nDo not use patient variants or phenotype. Do not invent targets. Preserve an exact contiguous source quote that contains the causal relationship. Prefer targets that can later be screened pharmacologically. Scores should reflect mechanistic strength on a 0-15 scale. Return only JSON."

def route_discovery_body(gene,docs,exclude_targets=None):
 cfg=MODELS['scout'];exclude_targets=list(exclude_targets or [])
 payload={'gene':gene,'objective':'discover pharmacologic rescue mechanisms for loss or insufficiency of this gene/protein','documents':docs,'already_extracted_modifiers':exclude_targets,'patient_context':False}
 instructions=ROUTE_DISCOVERY_INSTRUCTIONS
 if exclude_targets:instructions+='\nThis is an additional-recall pass. Do not return these already extracted modifiers: '+', '.join(exclude_targets)+'. Find only distinct additional causal mechanisms supported by the passages.'
 return {'model':cfg['model'],'input':dump(payload),'instructions':instructions,'store':False,'background':False,'tools':[],
  'max_output_tokens':2500,'reasoning':{'effort':cfg['reasoning']},
  'text':{'format':{'type':'json_schema','name':'v28_route_discovery','strict':True,'schema':_route_schema()}}}


def _gene_like_in_text(name,text):
 needle=norm_token(name)
 if not needle:return False
 if needle in norm_token(text):return True
 import difflib
 toks=re.findall(r'[A-Za-z0-9-]+',text)
 return any(difflib.SequenceMatcher(None,needle,norm_token(t)).ratio()>=0.72 for t in toks if len(norm_token(t))>=3)

def _causal_quote_ok(q,target_gene):
 lo=q.lower()
 causal=('inhibit','inhibition','knockdown','deplet','overexpress','exogenous expression','ameliorat','activat','deacetyl','ubiquit','degrad','stabil','recruit','rescue','suppres','phosphorylat','prevent','protect','binding')
 return len(q)>=30 and _gene_like_in_text(target_gene,q) and any(c in lo for c in causal)

def _repair_exact_quote(passage_text,model_quote,target_gene):
 if model_quote and model_quote in passage_text and _causal_quote_ok(model_quote,target_gene):return model_quote
 chunks=[x.strip(' .') for x in re.split(r'\s*(?:\.\.\.|…)+\s*',model_quote or '') if len(x.strip())>=30]
 for x in sorted(chunks,key=len,reverse=True):
  if x in passage_text and _causal_quote_ok(x,target_gene):return x
 # Last resort: exact sentence from source passage, ranked by causal content + target-name similarity.
 sents=[x.strip() for x in re.split(r'(?<=[.!?])\s+',passage_text) if len(x.strip())>=30]
 scored=[]
 for s in sents:
  if not _causal_quote_ok(s,target_gene):continue
  score=sum(1 for c in ('inhibit','knockdown','overexpress','deacetyl','ubiquit','degrad','stabil','recruit','rescue','prevent','protect') if c in s.lower())
  score+=min(4,len(set(re.findall(r'[a-z0-9]+',(model_quote or '').lower())) & set(re.findall(r'[a-z0-9]+',s.lower())))//3)
  scored.append((score,len(s),s))
 if scored:
  scored.sort(key=lambda z:(-z[0],z[1]));return scored[0][2]
 return None


def _repair_quote_in_document(doc,model_passage_id,model_quote,target_gene,mechanism_concept=''):
 passages=doc.get('passages',[])
 mech={w for w in re.findall(r'[a-z0-9]+',(mechanism_concept or '').lower()) if len(w)>=5}
 model_words={w for w in re.findall(r'[a-z0-9]+',(model_quote or '').lower()) if len(w)>=4}
 candidates={}
 def add(pid,q,bonus=0):
  if not q or not _causal_quote_ok(q,target_gene):return
  low=q.lower();words=set(re.findall(r'[a-z0-9]+',low))
  score=bonus+2*len(mech & words)+min(5,len(model_words & words)//3)
  if any(k in low for k in ('rescue','ameliorat','restore')):score+=8
  if any(k in low for k in ('prevent','suppress','protect')):score+=6
  if any(k in low for k in ('exogenous expression','overexpression','upregulat','increase')):score+=4
  if any(k in low for k in ('stabil','degrad','ubiquit','deacetyl','phosphorylat','recruit')):score+=3
  key=(pid,q)
  if key not in candidates or score>candidates[key]:candidates[key]=score
 # Exact model quote/chunks get provenance bonus, but do not automatically win.
 chunks=[x.strip(' .') for x in re.split(r'\s*(?:\.\.\.|…)+\s*',model_quote or '') if len(x.strip())>=30]
 for p in passages:
  if model_quote and model_quote in p['text']:add(p['id'],model_quote,5 if p['id']==model_passage_id else 3)
  for chunk in chunks:
   if chunk in p['text']:add(p['id'],chunk,3 if p['id']==model_passage_id else 2)
  for s in [x.strip() for x in re.split(r'(?<=[.!?])\s+',p['text']) if len(x.strip())>=30]:add(p['id'],s,1 if p['id']==model_passage_id else 0)
 if candidates:
  ranked=sorted(((score,-len(q),pid,q) for (pid,q),score in candidates.items()),reverse=True)
  return ranked[0][2],ranked[0][3]
 return None,None

def discover_routes(gene,docs,outdir,accession_resolver,exclude_targets=None,tag='route_discovery'):
 if not docs:return []
 body=route_discovery_body(gene,docs,exclude_targets)
 rec=_request(body,'scout',Path(outdir)/(tag+'_response.json'));value=rec['value'];bydoc={d['source_id']:d for d in docs};out=[]
 try:primary_acc=accession_resolver(gene)
 except:primary_acc=None
 for x in value.get('routes',[]):
  d=bydoc.get(x['source_id'])
  if not d:continue
  tg=(x.get('target_gene') or '').strip()
  if not tg or tg.upper()==gene.upper():continue
  repaired_pid,exact_quote=_repair_quote_in_document(d,x.get('passage_id'),x.get('quote',''),tg,x.get('mechanism_concept',''))
  if not exact_quote:continue
  qlo=exact_quote.lower()
  causal=('inhibit','inhibition','knockdown','deplet','overexpress','exogenous expression','ameliorat','activat','deacetyl','ubiquit','degrad','stabil','recruit','rescue','suppres','phosphorylat','prevent','protect','binding')
  if not any(c in qlo for c in causal):continue
  acc=None
  try:acc=accession_resolver(tg)
  except:acc=None
  if not acc or (primary_acc and acc==primary_acc):continue
  x['score']=min(15.0,max(8.0,float(x.get('score',10.0))))
  out.append({**x,'passage_id':repaired_pid,'quote':exact_quote,'model_passage_id_original':x.get('passage_id'),'model_quote_original':x.get('quote'),'target_accession':acc,'targets':[{'gene':tg,'accession':acc}],
    'evidence':[{'source_id':x['source_id'],'passage_id':repaired_pid,'quote':exact_quote,'license':d['license']}],
    'patient_context_used':False,'quote_repaired_deterministically':exact_quote!=x.get('quote') or repaired_pid!=x.get('passage_id')})
 save(Path(outdir)/(tag+'_routes.json'),out);return out

def _route_verify_schema():
 item={'type':'object','properties':{
  'source_id':{'type':'string'},'target_concept':{'type':'string'},
  'verdict':{'type':'string','enum':['ACCEPT','CORRECT','REJECT']},
  'corrected_action':{'type':['string','null'],'enum':['increase','decrease',None]},
  'evidence_context':{'type':'string','enum':['somatic','meiotic','mixed','unclear']},
  'reason':{'type':'string'}},
  'required':['source_id','target_concept','verdict','corrected_action','evidence_context','reason'],'additionalProperties':False}
 return {'type':'object','properties':{'decisions':{'type':'array','items':item,'maxItems':20}},'required':['decisions'],'additionalProperties':False}

ROUTE_VERIFY_INSTRUCTIONS="You are a strict causal-direction verifier for molecular rescue routes.\nUse only each supplied exact source quote and route description. The primary gene/protein is deficient or insufficient.\nFor each modifier, determine which manipulation of the modifier would move the primary gene/protein or its function toward rescue.\nIf the modifier promotes degradation/ubiquitination or suppresses the primary protein, rescue direction is usually decrease.\nIf the modifier protects/stabilizes/recruits/activates the primary protein, rescue direction is usually increase.\nExample: if X prevents degradation of the primary protein, increasing X is the rescue direction.\nReject association-only or logically unsupported routes. Identify meiotic/oocyte-only evidence separately from somatic evidence.\nDo not use patient data, outside knowledge, or treatment recommendations. Return JSON only."

def route_verify_body(gene,route_rows):
 cfg=MODELS['verify']
 payload={'primary_gene':gene,'routes':[{
   'source_id':x.get('source_id'),'target_concept':x.get('target_concept'),'target_gene':x.get('target_gene'),
   'proposed_action':x.get('desired_action'),'variable':x.get('variable'),'mechanism_concept':x.get('mechanism_concept'),
   'quote':x.get('quote')} for x in route_rows],
   'patient_context':False}
 return {'model':cfg['model'],'input':dump(payload),'instructions':ROUTE_VERIFY_INSTRUCTIONS,'store':False,'background':False,'tools':[],
  'max_output_tokens':1600,'reasoning':{'effort':cfg['reasoning']},
  'text':{'format':{'type':'json_schema','name':'v28_route_direction_verify','strict':True,'schema':_route_verify_schema()}}}

def verify_routes(gene,route_rows,outdir):
 if not route_rows:return []
 body=route_verify_body(gene,route_rows)
 rec=_request(body,'verify',Path(outdir)/'route_direction_verify_response.json')
 decisions=rec['value'].get('decisions',[]);keyed={(d['source_id'],d['target_concept']):d for d in decisions};out=[]
 for x in route_rows:
  d=keyed.get((x.get('source_id'),x.get('target_concept')))
  if not d or d['verdict']=='REJECT':continue
  y=dict(x)
  if d['verdict']=='CORRECT' and d.get('corrected_action') in ('increase','decrease'):
   y['desired_action']=d['corrected_action']
  y['route_direction_verification']=d
  y['evidence_context']=d.get('evidence_context')
  out.append(y)
 save(Path(outdir)/'verified_routes.json',out)
 return out

def _deep_schema():
 quote={'type':'object','properties':{'source_id':{'type':'string'},'passage_id':{'type':'string'},'quote':{'type':'string'}},
        'required':['source_id','passage_id','quote'],'additionalProperties':False}
 return {'type':'object','properties':{
  'candidate':{'type':'string'},'adjudication':{'type':'string','enum':['SUPPORTED','INSUFFICIENT','OPPOSED','CONFLICTED']},
  'same_drug_functional_rescue':{'type':'boolean'},'direct_route_modulation':{'type':'boolean'},'exposure_reported':{'type':'boolean'},
  'counterevidence_present':{'type':'boolean'},'support_quotes':{'type':'array','items':quote,'maxItems':8},
  'counter_quotes':{'type':'array','items':quote,'maxItems':8},'reason':{'type':'string'},'limitations':{'type':'array','items':{'type':'string'},'maxItems':10}},
  'required':['candidate','adjudication','same_drug_functional_rescue','direct_route_modulation','exposure_reported','counterevidence_present','support_quotes','counter_quotes','reason','limitations'],
  'additionalProperties':False}

def deep_attack(candidate,route,docs,outdir,role='verify',regulatory_approval=None):
 cfg=MODELS[role];payload={'candidate':candidate,'route':{k:route.get(k) for k in ['route_id','target_concept','desired_action','variable','mechanism_concept']},
  'documents':docs,'patient_context':False,'registry_facts':registry_output.compact_for_ai(regulatory_approval) if regulatory_approval else {}}
 instructions="""Strictly verify the named candidate against the named molecular route using only supplied public passages.
Do not infer rescue from binding, expression change, cytotoxicity, tumor growth inhibition, or cell-cycle arrest.
same_drug_functional_rescue=true only for a direct functional endpoint relevant to the route/gene dysfunction.
Copy every quote exactly and contiguously. Registry facts (if supplied) are read-only context and must never be extrapolated to other countries or therapeutic indications. If exact evidence is absent, return INSUFFICIENT. Return JSON only."""
 body={'model':cfg['model'],'input':dump(payload),'instructions':instructions,'store':False,'background':False,'tools':[],
  'max_output_tokens':1800,'reasoning':{'effort':cfg['reasoning']},
  'text':{'format':{'type':'json_schema','name':'v28_deep_attack','strict':True,'schema':_deep_schema()}}}
 rec=_request(body,role,Path(outdir)/f'{role}_response.json');v=rec['value'];bydoc={d['source_id']:d for d in docs}
 if v.get('candidate')!=candidate:raise RuntimeError('V28_DEEP_BINDING_MISMATCH')
 for q in v.get('support_quotes',[])+v.get('counter_quotes',[]):
  d=bydoc.get(q['source_id']);p=next((p for p in d.get('passages',[]) if p['id']==q['passage_id']),None) if d else None
  if not p or q['quote'] not in p['text']:raise RuntimeError('V28_DEEP_QUOTE_NOT_EXACT')
 if regulatory_approval is not None:
  v['regulatory_approval']=regulatory_approval
  save(Path(outdir)/f'{role}_annotated_research_view.json',v)
 return v

def _shortlist_schema():
 item={'type':'object','properties':{
  'candidate':{'type':'string'},
  'route_id':{'type':'string'},
  'why_not_first_choice':{'type':'string'},
  'why_worth_deepening':{'type':'string'},
  'missing_evidence':{'type':'array','items':{'type':'string'},'maxItems':6},
  'next_checks':{'type':'array','items':{'type':'string'},'maxItems':6},
  'research_priority':{'type':'string','enum':['HIGH','MEDIUM','EXPLORATORY']},
  'confidence':{'type':'string','enum':['low','moderate','high']}
 },'required':['candidate','route_id','why_not_first_choice','why_worth_deepening','missing_evidence','next_checks','research_priority','confidence'],
 'additionalProperties':False}
 return {'type':'object','properties':{'items':{'type':'array','items':item,'maxItems':12}},'required':['items'],'additionalProperties':False}

SHORTLIST_INSTRUCTIONS="""You are explaining a research shortlist produced by a rare-disease drug-discovery pipeline.
Use ONLY the supplied public/mechanistic evidence summary. Do not make treatment recommendations.
For every candidate explain two different things:
1) why it is NOT a first-choice validated research lead yet;
2) why it is still worth deeper investigation.
A candidate is not first choice when evidence is indirect, route modulation is unverified, same-drug functional rescue is absent, exposure relevance is missing, evidence is contradictory, or evidence strength/rank is below the finalist threshold.
Do not call NO_SIGNAL proof of ineffectiveness. Do not convert network proximity, transcriptomics, structural similarity, or target binding into functional rescue.
Prioritize missing evidence that could actually change the verdict: direct candidate-to-route modulation, same-drug rescue of the disease-relevant functional phenotype, concentration/exposure plausibility, and counterevidence.
The payload contains no patient-specific phenotype or variants; do not infer them. Regulatory authorisations must only be stated for exact provided official source facts and their own jurisdiction; no off-label efficacy inference. Return JSON only."""

def research_shortlist_rationales(rows,outdir):
 if not rows:return []
 payload={'patient_context':False,'clinical_recommendation':False,'candidates':[]}
 for x in rows[:12]:
  route=x.get('route') or {}
  deep=x.get('deep_public')
  payload['candidates'].append({
   'candidate':x['candidate'],'route_id':x['route_id'],'shortlist_class':x['shortlist_class'],
   'tier':x.get('tier'),'screen_score':x.get('screen_score'),
   'route':{k:route.get(k) for k in ['target_concept','desired_action','variable','mechanism_concept','family','state']},
   'providers':x.get('providers',[]),
   'signals':[{k:s.get(k) for k in ['provider','depth','value_nM','support_count','tanimoto','functional_rescue']} for s in x.get('signals',[])[:12]],
   'deep_verification':deep,
   'public_blockers':x.get('public_blockers',[]),
   'registry_facts':registry_output.compact_for_ai(x['regulatory_approval']) if x.get('regulatory_approval') else {}
  })
 cfg=MODELS['scout']
 body={'model':cfg['model'],'input':dump(payload),'instructions':SHORTLIST_INSTRUCTIONS,'store':False,'background':False,'tools':[],
  'max_output_tokens':3500,'reasoning':{'effort':cfg['reasoning']},
  'text':{'format':{'type':'json_schema','name':'v28_research_shortlist','strict':True,'schema':_shortlist_schema()}}}
 rec=_request(body,'scout',Path(outdir)/'research_shortlist_ai_response.json')
 vals=rec['value'].get('items',[]);expected={(x['candidate'],x['route_id']) for x in rows[:12]};out=[]
 for v in vals:
  if (v.get('candidate'),v.get('route_id')) in expected:
   original=next(x for x in rows[:12] if (x['candidate'],x['route_id'])==(v['candidate'],v['route_id']))
   out.append({**v,'regulatory_approval':original.get('regulatory_approval')})
 save(Path(outdir)/'research_shortlist_ai_rationales.json',out)
 return out
