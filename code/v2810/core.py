"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

from pathlib import Path
import json,hashlib,datetime,re,urllib.request,urllib.parse,os,tempfile
ROOT=Path(__file__).resolve().parents[2]
DATA_ROOT=ROOT.parent/'data'
PUBLIC=ROOT/'public_inputs'
RAW=PUBLIC/'third_party'
RUN_ROOT=ROOT/'experiments'
VERSION='2.8.10-open-research-auto-resources'
PROJECT_TITLE='Rare Disease MVA Pharma — V2.8.10'
PROJECT_CREDIT='Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch'

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def norm_name(s):
 s='' if s is None else str(s)
 s=s.replace('\u2010','-').replace('\u2011','-').replace('\u2013','-').replace('\u2014','-')
 return re.sub(r'\s+',' ',s.strip().upper())
def norm_token(s): return re.sub(r'[^a-z0-9]+','',str(s).lower())
def stable_id(prefix,label): return prefix+hashlib.sha256(label.encode('utf-8')).hexdigest()[:16]
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(x): return json.dumps(x,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+'\n'
def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def save(path,x):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 fd,tmp=tempfile.mkstemp(prefix='.v28-',dir=p.parent)
 try:
  with os.fdopen(fd,'w',encoding='utf-8',newline='\n') as f:
   f.write(dump(x));f.flush();os.fsync(f.fileno())
  for attempt in range(80):
   try:
    os.replace(tmp,p)
    break
   except PermissionError:
    if attempt==79:
     if p.name in ('STATUS.json','STATUS_V30_NATIVE.json'):
      os.replace(tmp,p.with_name(p.stem+'_FINAL'+p.suffix))
      return
     raise
    # Windows Defender/indexers and remote file readers can transiently hold
    # JSON outputs open. Back off without weakening atomic replacement.
    __import__('time').sleep(min(1.0,0.20+0.02*attempt))
 finally:
  if os.path.exists(tmp):os.unlink(tmp)
def jsonl(path,rows):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf-8',newline='\n') as f:
  for x in rows:f.write(json.dumps(x,ensure_ascii=False,sort_keys=True)+'\n')
def read_jsonl(path):
 with Path(path).open(encoding='utf-8') as f:
  for line in f:
   if line.strip():yield json.loads(line)

def fetch_json(url,timeout=30,cap=30_000_000,ua='RareDiseaseMVAPharma-V28/1.0'):
 if url.startswith('https://www.ebi.ac.uk/chembl/api/data/'):
  from .chembl_sqlite_v31 import configured,local_json,PROVENANCE
  if configured() or PROVENANCE.is_file():
   return local_json(url)
 req=urllib.request.Request(url,headers={'User-Agent':ua,'Accept':'application/json'})
 with urllib.request.urlopen(req,timeout=timeout) as r:
  b=r.read(cap+1)
  if len(b)>cap:raise RuntimeError('response_size_cap')
 return json.loads(b)

def fetch_bytes(url,timeout=30,cap=30_000_000,ua='RareDiseaseMVAPharma-V28/1.0'):
 req=urllib.request.Request(url,headers={'User-Agent':ua})
 with urllib.request.urlopen(req,timeout=timeout) as r:
  b=r.read(cap+1)
  if len(b)>cap:raise RuntimeError('response_size_cap')
 return b

def normalize_case(case_path):
 raw=load(case_path)
 primary={norm_name(x) for x in raw.get('primary_genes',[]) if norm_name(x)}
 genes={}
 for g in raw.get('genetic_candidates',[]):
  gene=norm_name(g.get('gene'))
  if not gene:continue
  if primary and gene not in primary:continue
  if g.get('selected_for_downstream') is False:continue
  genes.setdefault(gene,[]).append({
   'hgvsc':g.get('hgvsc'),'hgvsp':g.get('hgvsp'),'genotype_class':g.get('genotype_class'),
   'interpretation':g.get('interpretation'),'consequence':g.get('consequence'),'impact_score':g.get('impact_score'),
   'clinvar':g.get('clinvar'),'chrom':g.get('chrom'),'pos':g.get('pos'),'ref':g.get('ref'),'alt':g.get('alt'),
   'GT':g.get('GT'),'AD':g.get('AD'),'DP':g.get('DP'),'GQ':g.get('GQ'),'FILTER':g.get('FILTER'),
   'phase_tags':g.get('phase_tags'),'mane_tags':g.get('mane_tags'),'refseq_rna':g.get('refseq_rna')})
 phen=[{
  'id':x.get('id'),'hpo_id':x.get('HPO_id'),'term':x.get('HPO_term_verbatim') or x.get('clinical_feature_verbatim'),
  'presence':x.get('presence'),'subject':x.get('subject')} for x in raw.get('phenotype_features',[])]
 variant_classes=[]
 for gene,vs in genes.items():
  for v in vs:
   cons=(v.get('consequence') or '').lower();p=(v.get('hgvsp') or '').lower();c=(v.get('hgvsc') or '').lower()
   if 'frameshift' in cons or 'fs' in p:typ='frameshift'
   elif 'stop_gained' in cons or 'ter' in p or '*' in p:typ='nonsense'
   elif 'splice' in cons or 'splice' in p or re.search(r'[+-][12][a-z]?>',c):typ='splice'
   elif 'missense' in cons or (p and 'p.' in p):typ='missense_or_other_coding'
   elif cons:typ=cons
   else:typ='unknown'
   variant_classes.append({'gene':gene,'class':typ})
 private={
  'schema_version':'V28_CASE_PRIVATE_2','case_id':raw.get('case_id') or stable_id('CASE:',sha(case_path)),
  'source_path':str(Path(case_path)),'source_sha256':sha(case_path),'genes':genes,'variant_classes':variant_classes,
  'phenotype_features':phen,'clinical_narrative':raw.get('clinical_narrative',[]),
  'phase':raw.get('phase'),'phase_evidence':raw.get('phase_evidence'),'residual_protein_status':raw.get('residual_protein_status','unknown'),
  'expression_measured':bool(raw.get('expression_measured_in_this_case',False)),
  'compound_heterozygosity_confirmed':bool(raw.get('compound_heterozygosity_confirmed',False)),
  'phenotype_review_domains':raw.get('phenotype_review_domains',{}),
  'genome_wide_analysis':raw.get('genome_wide_analysis',{}),'genome_wide_candidate_genes':raw.get('genome_wide_candidate_genes',[]),
  'diagnostic_certainty':raw.get('diagnostic_certainty','unknown'),
  'preferred_drugs':raw.get('preferred_drugs',[]),'excluded_drug_names':raw.get('excluded_drug_names',[]),
  'patient_specific_api_payload_allowed':False}
 public={
  'schema_version':'V28_PUBLIC_SAFE_CASE_ENVELOPE_2','case_id_hash':hashlib.sha256(private['case_id'].encode()).hexdigest()[:12],
  'genes':sorted(genes),'variant_classes':sorted({x['class'] for x in variant_classes}),
  'generic_objectives':['restore or compensate primary gene/protein function','identify pharmacologic rescue mechanisms with functional endpoints'],
  'patient_variant_strings_included':False,'phenotype_combination_included':False,'genome_wide_variant_data_included':False}
 return private,public
