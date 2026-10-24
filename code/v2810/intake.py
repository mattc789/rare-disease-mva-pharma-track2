"""Rare Disease MVA Pharma — V2.8.6
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""

"""V2.8 Raw Intake Engine.

Purpose
-------
Turn raw patient inputs into a provenance-complete case model for the autonomous V2.8
pipeline. Patient variants never leave the workstation.

Inputs:
- germline/somatic VCF(.gz), single proband sample
- optional DOCX containing HPO-coded clinical phenotype
- public local resources (GENCODE MANE, HPO, ClinVar)

Outputs:
- case_model_v2.json: compact downstream case model
- genome_wide_variants.jsonl: every MANE coding/splice variant considered
- genome_wide_genes.json: ranked candidate genes
- intake_manifest.json: hashes/resource versions/QC/limitations
"""
from pathlib import Path
from collections import defaultdict,Counter
import gzip,csv,io,json,re,zipfile,xml.etree.ElementTree as ET,urllib.request,os,bisect,math,hashlib
from .core import *
from . import case_science

RES=PUBLIC/'third_party/intake_v28'
RESOURCE_URLS={
 'gencode.v50.annotation.gtf.gz':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz',
 'gencode.v50.pc_transcripts.fa.gz':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.pc_transcripts.fa.gz',
 'gencode.v50.metadata.RefSeq.gz':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.metadata.RefSeq.gz',
 'genes_to_phenotype.txt':'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-02/genes_to_phenotype.txt',
 'hp.obo':'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-02/hp.obo',
 'clinvar.vcf.gz':'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz',
}
PRIMARY_CHROMS={str(i) for i in range(1,23)}|{'X','Y','MT'}
BIN=100_000
SEVERITY={
 'stop_gained':12,'frameshift_variant':12,'splice_acceptor_variant':11,'splice_donor_variant':11,
 'start_lost':10,'stop_lost':9,'missense_variant':6,'inframe_indel':5,'splice_region_variant':4,
 'coding_unknown':4,'synonymous_variant':0,'utr_variant':0,'exonic_noncoding':0,
}
AA3={
'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly','H':'His','I':'Ile',
'L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','*':'Ter'
}
CODON={
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
 'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
 'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
COMP=str.maketrans('ACGTNacgtn','TGCANtgcan')

def chrom_norm(c):
 c=str(c).removeprefix('chr')
 return 'MT' if c in ('M','MT') else c

def _resource_valid(p):
 p=Path(p)
 if not p.exists() or p.stat().st_size<1000:return False
 if p.suffix=='.gz':
  try:
   with gzip.open(p,'rb') as f:
    for _ in iter(lambda:f.read(4*1024*1024),b''):pass
   return True
  except Exception:return False
 return True

def ensure_resources():
 RES.mkdir(parents=True,exist_ok=True);downloaded=[];stamp=RES/'resource_validation.json'
 if stamp.exists():
  try:
   s=load(stamp);files=s.get('files',{})
   if all((RES/n).exists() and (RES/n).stat().st_size==files.get(n,{}).get('size') for n in RESOURCE_URLS):
    return {'resource_dir':str(RES),'downloaded':[],'hashes':{n:files[n]['sha256'] for n in RESOURCE_URLS},'validated_from_stamp':True}
  except Exception:pass
 for name,url in RESOURCE_URLS.items():
  p=RES/name
  if _resource_valid(p):continue
  tmp=p.with_suffix(p.suffix+'.part')
  req=urllib.request.Request(url,headers={'User-Agent':'RareDiseaseMVAPharma-V28-Intake/1.0'})
  with urllib.request.urlopen(req,timeout=180) as r, tmp.open('wb') as f:
   while True:
    b=r.read(1024*1024)
    if not b:break
    f.write(b)
  os.replace(tmp,p);downloaded.append(name)
 files={n:{'size':(RES/n).stat().st_size,'sha256':sha(RES/n),'url':u} for n,u in RESOURCE_URLS.items()}
 save(stamp,{'version':'V28_INTAKE_RESOURCES_1','validated_utc':now(),'files':files})
 return {'resource_dir':str(RES),'downloaded':downloaded,'hashes':{n:files[n]['sha256'] for n in RESOURCE_URLS},'validated_from_stamp':False}

def _attrs(s):
 out=defaultdict(list)
 for piece in s.strip().strip(';').split(';'):
  piece=piece.strip()
  if not piece:continue
  m=re.match(r'(\S+)\s+"(.*)"$',piece)
  if m:out[m.group(1)].append(m.group(2))
 return out

def _open_text(p):
 p=Path(p)
 return gzip.open(p,'rt',encoding='utf-8',errors='replace') if p.suffix=='.gz' else p.open(encoding='utf-8',errors='replace')

def _split_txid(raw):
 raw=(raw or '').strip()
 m=re.match(r'^(ENST\d+)(?:\.(\d+))?$',raw)
 return (m.group(1),m.group(2)) if m else (raw,None)

def _transcript_key(a):
 raw=(a.get('transcript_id') or [None])[0]
 tid,ver0=_split_txid(raw)
 ver=(a.get('transcript_version') or [None])[0] or ver0
 return (tid+'.'+ver) if tid and ver else tid

def build_mane_cache():
 out=RES/'mane_models_v50.json.gz'
 if out.exists():return out
 gtf=RES/'gencode.v50.annotation.gtf.gz'
 selected=set();meta={}
 with _open_text(gtf) as f:
  for line in f:
   if line.startswith('#'):continue
   p=line.rstrip('\n').split('\t')
   if len(p)<9 or p[2]!='transcript':continue
   a=_attrs(p[8]);tags=set(a.get('tag',[]))
   if not tags.intersection({'MANE_Select','MANE_Plus_Clinical'}):continue
   raw_tid=(a.get('transcript_id') or [None])[0]
   tid,ver0=_split_txid(raw_tid)
   if not tid:continue
   ver=(a.get('transcript_version') or [None])[0] or ver0
   selected.add(tid);meta[tid]={
    'transcript_id':tid,'transcript_version':ver,
    'gene_id':(a.get('gene_id') or [None])[0],'gene_name':(a.get('gene_name') or [None])[0],
    'gene_type':(a.get('gene_type') or [None])[0],'transcript_type':(a.get('transcript_type') or [None])[0],
    'chrom':chrom_norm(p[0]),'strand':p[6],'tags':sorted(tags),'exons':[],'cds':[]}
 with _open_text(gtf) as f:
  for line in f:
   if line.startswith('#'):continue
   p=line.rstrip('\n').split('\t')
   if len(p)<9 or p[2] not in ('exon','CDS'):continue
   a=_attrs(p[8]);tid,_=_split_txid((a.get('transcript_id') or [None])[0])
   if tid not in selected:continue
   rec={'start':int(p[3]),'end':int(p[4]),'frame':p[7],
        'exon_number':int((a.get('exon_number') or [0])[0] or 0)}
   meta[tid]['exons' if p[2]=='exon' else 'cds'].append(rec)
 # RefSeq cross-reference
 refmap=defaultdict(list);rp=RES/'gencode.v50.metadata.RefSeq.gz'
 if rp.exists():
  with _open_text(rp) as f:
   for line in f:
    p=line.rstrip('\n').split('\t')
    if len(p)>=2:refmap[_split_txid(p[0])[0]].append(p[1])
 # sequences
 seqs={};fa=RES/'gencode.v50.pc_transcripts.fa.gz'
 wanted=set(selected)
 with _open_text(fa) as f:
  key=None;buf=[]
  for line in f:
   if line.startswith('>'):
    if key in wanted:seqs[key]=''.join(buf).upper()
    head=line[1:].split('|',1)[0];key=head.split('.')[0];buf=[]
   else:
    if key in wanted:buf.append(line.strip())
  if key in wanted:seqs[key]=''.join(buf).upper()
 models={}
 for tid,m in meta.items():
  ex=sorted(m['exons'],key=lambda x:x['start'],reverse=m['strand']=='-')
  offset=0
  for e in ex:e['cdna_start0']=offset;offset+=e['end']-e['start']+1
  m['exons']=ex;m['transcript_length']=offset;m['sequence']=seqs.get(tid)
  cpos=[]
  for c in m['cds']:
   vals=[]
   for pos in (c['start'],c['end']):
    z=_genomic_to_cdna(m,pos)
    if z is not None:vals.append(z)
   if vals:cpos+=vals
  m['cds_cdna_start0']=min(cpos) if cpos else None
  m['cds_cdna_end0']=max(cpos)+1 if cpos else None
  m['refseq_rna']=sorted(set(refmap.get(tid,[])))
  models[tid]=m
 payload={'version':'GENCODE50_MANE_LOCAL_1','source_hashes':{'gtf':sha(gtf),'pc_transcripts':sha(fa),'refseq':sha(rp) if rp.exists() else None},
          'transcripts':models}
 with gzip.open(out,'wt',encoding='utf-8') as f:json.dump(payload,f,ensure_ascii=False,separators=(',',':'))
 return out

def load_mane_models():
 p=build_mane_cache()
 with gzip.open(p,'rt',encoding='utf-8') as f:return json.load(f)

def _genomic_to_cdna(m,pos):
 for e in m.get('exons',[]):
  if e['start']<=pos<=e['end']:
   return e['cdna_start0']+(pos-e['start'] if m['strand']=='+' else e['end']-pos)
 return None

def _coding_overlap(m,pos,end):
 for c in m.get('cds',[]):
  if not (end<c['start'] or pos>c['end']):return True
 return False

def _splice_consequence(m,pos,end):
 # Only canonical +/-2 essential sites and +/-8 region; exon-overlapping variants are handled first.
 for e in m.get('exons',[]):
  for boundary,kind_plus,kind_minus in [(e['start'],'splice_acceptor_variant','splice_donor_variant'),(e['end'],'splice_donor_variant','splice_acceptor_variant')]:
   d=min(abs(pos-boundary),abs(end-boundary))
   if d<=2:return kind_plus if m['strand']=='+' else kind_minus
   if d<=8:return 'splice_region_variant'
 return None

def _translate(c):return CODON.get(c.upper(),'X')
def _revcomp(s):return s.translate(COMP)[::-1]

def annotate_variant_on_transcript(m,pos,ref,alt):
 end=pos+len(ref)-1;seq=m.get('sequence');cds0=m.get('cds_cdna_start0');cds1=m.get('cds_cdna_end0')
 cdna=_genomic_to_cdna(m,pos)
 in_exon=cdna is not None
 if not in_exon:
  sc=_splice_consequence(m,pos,end)
  if not sc:return None
  return {'consequence':sc,'severity':SEVERITY[sc],'transcript_id':m['transcript_id'],'transcript_version':m.get('transcript_version'),
   'refseq_rna':m.get('refseq_rna',[]),'gene':m['gene_name'],'gene_id':m['gene_id'],'mane_tags':m['tags'],'strand':m['strand']}
 if cds0 is None or cds1 is None or not _coding_overlap(m,pos,end):
  return {'consequence':'utr_variant' if m.get('transcript_type')=='protein_coding' else 'exonic_noncoding','severity':0,
   'transcript_id':m['transcript_id'],'transcript_version':m.get('transcript_version'),'refseq_rna':m.get('refseq_rna',[]),
   'gene':m['gene_name'],'gene_id':m['gene_id'],'mane_tags':m['tags'],'strand':m['strand']}
 rec={'transcript_id':m['transcript_id'],'transcript_version':m.get('transcript_version'),'refseq_rna':m.get('refseq_rna',[]),
      'gene':m['gene_name'],'gene_id':m['gene_id'],'mane_tags':m['tags'],'strand':m['strand']}
 if len(ref)!=1 or len(alt)!=1:
  delta=len(alt)-len(ref);cons='inframe_indel' if delta%3==0 else 'frameshift_variant'
  rec.update({'consequence':cons,'severity':SEVERITY[cons],'coding_position':cdna-cds0+1 if cdna>=cds0 else None,
              'protein_change':None,'annotation_note':'indel consequence class; exact HGVS protein not generated by lightweight local annotator'})
  return rec
 if not seq or not (cds0<=cdna<cds1):
  rec.update({'consequence':'coding_unknown','severity':SEVERITY['coding_unknown'],'annotation_note':'coding interval but transcript sequence unavailable'})
  return rec
 tref=ref.upper() if m['strand']=='+' else _revcomp(ref.upper());talt=alt.upper() if m['strand']=='+' else _revcomp(alt.upper())
 observed=seq[cdna] if cdna<len(seq) else None;coding0=cdna-cds0;codon0=cds0+(coding0//3)*3
 codon=seq[codon0:codon0+3]
 if len(codon)!=3:
  rec.update({'consequence':'coding_unknown','severity':SEVERITY['coding_unknown'],'reference_sequence_mismatch':observed!=tref});return rec
 mut=list(codon);mut[coding0%3]=talt;mut=''.join(mut);aa1=_translate(codon);aa2=_translate(mut);aa_pos=coding0//3+1
 if observed!=tref:cons='coding_unknown'
 elif aa1==aa2:cons='synonymous_variant'
 elif aa1=='M' and aa_pos==1 and aa2!='M':cons='start_lost'
 elif aa1!='*' and aa2=='*':cons='stop_gained'
 elif aa1=='*' and aa2!='*':cons='stop_lost'
 else:cons='missense_variant'
 hgvs_c=f"c.{coding0+1}{tref}>{talt}"
 hgvs_p='p.'+(AA3.get(aa1,aa1)+str(aa_pos)+AA3.get(aa2,aa2))
 rec.update({'consequence':cons,'severity':SEVERITY[cons],'coding_position':coding0+1,'protein_position':aa_pos,
   'reference_codon':codon,'alternate_codon':mut,'reference_aa':aa1,'alternate_aa':aa2,'hgvsc_local':hgvs_c,'hgvsp_local':hgvs_p,
   'reference_sequence_mismatch':observed!=tref})
 return rec

def build_interval_bins(models):
 # Precompute only transcript IDs per 100 kb bucket. Variant lookup is O(k) for the
 # small local bucket rather than O(all transcripts); exact exon +/-8 overlap is checked later.
 bins=defaultdict(lambda:defaultdict(set))
 for tid,m in models.items():
  if m['chrom'] not in PRIMARY_CHROMS:continue
  starts=[];ends=[]
  for e in m.get('exons',[]):starts.append(e['start']);ends.append(e['end'])
  if not starts:continue
  a=max(1,min(starts)-8);b=max(ends)+8
  # Do NOT fill every bucket across long introns: add only buckets touched by an exon.
  touched=set()
  for e in m.get('exons',[]):
   x=max(1,e['start']-8);y=e['end']+8
   touched.update(range(x//BIN,y//BIN+1))
  for bucket in touched:bins[m['chrom']][bucket].add(tid)
 return {c:{b:sorted(v) for b,v in d.items()} for c,d in bins.items()}

def parse_clinical_docx(path):
 if not path:return {'features':[],'narrative':[],'source_present':False}
 p=Path(path);features=[];narr=[]
 with zipfile.ZipFile(p) as z:
  xml=z.read('word/document.xml')
 root=ET.fromstring(xml);ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 for pi,par in enumerate(root.findall('.//w:p',ns),1):
  t=''.join(x.text or '' for x in par.findall('.//w:t',ns)).strip()
  if t:narr.append({'locator':f'word/document.xml p{pi}','text':t})
 for ti,table in enumerate(root.findall('.//w:tbl',ns),1):
  for ri,row in enumerate(table.findall('w:tr',ns),1):
   cells=[' '.join(''.join(t.itertext()) for t in c.findall('.//w:t',ns)).strip() for c in row.findall('w:tc',ns)]
   ids=[m for x in cells for m in re.findall(r'HP:\d{7}',x)]
   if len(ids)!=1:continue
   subject='family_history' if re.search(r'parent|famil',cells[0] if cells else '',re.I) else 'proband'
   features.append({'id':'F'+str(len(features)+1).zfill(2),'clinical_feature_verbatim':cells[0] if cells else None,
    'HPO_term_verbatim':cells[1] if len(cells)>1 else None,'HPO_id':ids[0],'notes_verbatim':cells[3] if len(cells)>3 else None,
    'subject':subject,'presence':case_science.presence_from_row(cells),'locator':f'word/document.xml table{ti} row{ri}','source_sha256':sha(p)})
 patient={x['HPO_id'] for x in features if x['subject']=='proband'}
 review={'oncologic_context_requires_review':'HP:0002859' in patient,'renal_history_requires_review':'HP:0000121' in patient,
  'growth_and_nutrition_requires_review':bool(patient&{'HP:0004322','HP:0001508','HP:0003202','HP:0001518'}),
  'current_medications_known':False,'renal_failure_inferred':False,'current_cancer_activity_inferred':False,
  'current_age_weight_known':False,'pharmacological_vetoes':[]}
 return {'features':features,'narrative':narr,'source_present':True,'source_sha256':sha(p),'phenotype_review_domains':review}

def parse_hpo_resources():
 gp=RES/'genes_to_phenotype.txt';obo=RES/'hp.obo'
 gene_terms=defaultdict(set)
 with gp.open(encoding='utf-8',errors='replace') as f:
  rd=csv.DictReader(f,delimiter='\t')
  for r in rd:
   g=(r.get('gene_symbol') or r.get('Gene Symbol') or '').upper();h=r.get('hpo_id') or r.get('HPO_ID')
   if g and h:gene_terms[g].add(h)
 parents=defaultdict(set);cur=None
 with obo.open(encoding='utf-8',errors='replace') as f:
  for line in f:
   line=line.rstrip()
   if line=='[Term]':cur=None
   elif line.startswith('id: HP:'):cur=line.split('id: ',1)[1]
   elif cur and line.startswith('is_a: HP:'):parents[cur].add(line.split()[1])
 return gene_terms,parents

def _anc(term,parents,maxdepth=8):
 seen={term};front={term}
 for _ in range(maxdepth):
  nxt=set()
  for x in front:nxt|=parents.get(x,set())
  nxt-=seen
  if not nxt:break
  seen|=nxt;front=nxt
 return seen

def phenotype_gene_score(gene,features,gene_terms,parents):
 terms=gene_terms.get(gene.upper(),set())
 if not terms:return {'score':0.0,'matched':[],'annotated_gene':False}
 ganc={g:_anc(g,parents) for g in terms};matches=[];total=0.0;den=0.0
 for f in features:
  if f.get('presence')!='reported_present':continue
  h=f['HPO_id'];w=1.0 if f.get('subject')=='proband' else 0.4;den+=w
  pa=_anc(h,parents);best=0.0;hit=None
  for g in terms:
   if h==g:s=1.0
   elif g in pa:s=0.75
   elif h in ganc[g]:s=0.6
   else:
    common=(pa & ganc[g])-{'HP:0000001'}
    s=0.25 if common else 0.0
   if s>best:best=s;hit=g
  total+=w*best
  if best>0:matches.append({'patient_hpo':h,'gene_hpo':hit,'similarity':best,'subject':f.get('subject')})
 return {'score':round(10*total/max(den,1),3),'matched':matches,'annotated_gene':True}

def _parse_info(s):
 d={}
 for x in s.split(';'):
  if '=' in x:
   k,v=x.split('=',1);d[k]=v
  elif x:d[x]=True
 return d

def _gt_alt_indices(gt):
 if not gt or gt in ('.','./.','.|.'):return []
 out=[]
 for x in re.split(r'[/|]',gt):
  try:i=int(x)
  except:continue
  if i>0:out.append(i)
 return sorted(set(out))

def detect_build(header):
 ref=' '.join(x for x in header if x.startswith('##reference='))
 lengths={}
 for x in header:
  m=re.match(r'##contig=<ID=([^,>]+),length=(\d+)',x)
  if m:lengths[chrom_norm(m.group(1))]=int(m.group(2))
 if lengths.get('1')==248956422 and lengths.get('15')==101991189:return 'GRCh38'
 if 'GRCh38' in ref or 'GCA_000001405.15' in ref:return 'GRCh38'
 return 'unknown'

def scan_vcf(vcf_path,models,bins,out_jsonl):
 qc={'records':0,'PASS':0,'nonref_records':0,'symbolic_sv_records':0,'candidate_overlap_records':0,'filters':Counter(),'GT':Counter()}
 rows=[];header=[];sample_names=[]
 op=gzip.open if str(vcf_path).endswith('.gz') else open
 with op(vcf_path,'rt',encoding='utf-8',errors='replace') as f:
  for line in f:
   if line.startswith('##'):header.append(line.rstrip());continue
   if line.startswith('#CHROM'):
    cols=line.rstrip().split('\t');sample_names=cols[9:]
    if len(sample_names)!=1:raise RuntimeError('V28_INTAKE_REQUIRES_ONE_PROBAND_SAMPLE')
    continue
   p=line.rstrip('\n').split('\t')
   if len(p)<10:continue
   qc['records']+=1
   if qc['records']%500000==0:print('V28_INTAKE_VCF_RECORDS',qc['records'],flush=True)
   chrom=chrom_norm(p[0]);pos=int(p[1]);ref=p[3];alts=p[4].split(',');flt=p[6]
   qc['filters'][flt]+=1
   if flt=='PASS':qc['PASS']+=1
   info=_parse_info(p[7]);fmt=p[8].split(':');vals=p[9].split(':');sd=dict(zip(fmt,vals));gt=sd.get('GT');qc['GT'][gt]+=1
   idxs=_gt_alt_indices(gt)
   if not idxs:continue
   qc['nonref_records']+=1
   if any(a.startswith('<') or '[' in a or ']' in a for a in alts) or 'SVTYPE' in info:qc['symbolic_sv_records']+=1
   tids=bins.get(chrom,{}).get(pos//BIN,[])
   if not tids:continue
   per_alt=[]
   for ai in idxs:
    if ai>len(alts):continue
    alt=alts[ai-1];end=pos+len(ref)-1;anns=[]
    for tid in tids:
     m=models[tid]
     # cheap overlap test
     if not any(not (end<e['start']-8 or pos>e['end']+8) for e in m.get('exons',[])):continue
     a=annotate_variant_on_transcript(m,pos,ref,alt)
     if a:anns.append(a)
    if not anns:continue
    qc['candidate_overlap_records']+=1
    anns.sort(key=lambda x:(-x['severity'],0 if 'MANE_Select' in x.get('mane_tags',[]) else 1,x['transcript_id']))
    best=anns[0]
    ad=sd.get('AD');ab=None
    if ad:
     try:
      counts=[int(x) for x in ad.split(',')];den=sum(counts)
      if den and ai<len(counts):ab=counts[ai]/den
     except:pass
    row={'chrom':chrom,'pos':pos,'ref':ref,'alt':alt,'id':None if p[2]=='.' else p[2],'qual':None if p[5]=='.' else float(p[5]),
      'filter':flt,'GT':gt,'AD':ad,'DP':int(sd['DP']) if sd.get('DP','').isdigit() else None,'GQ':int(sd['GQ']) if sd.get('GQ','').isdigit() else None,
      'PGT':sd.get('PGT'),'PID':sd.get('PID'),'PS':sd.get('PS'),'allele_balance':round(ab,4) if ab is not None else None,
      'gene':best['gene'],'gene_id':best['gene_id'],'best_annotation':best,'all_annotations':anns,'clinvar':None}
    rows.append(row)
 jsonl(out_jsonl,rows)
 qc['filters']=dict(qc['filters']);qc['GT']=dict(qc['GT']);qc['build']=detect_build(header);qc['samples']=sample_names
 qc['sv_callset_interpretation']='symbolic_SV_records_present' if qc['symbolic_sv_records'] else 'no_symbolic_SV_records_in_this_VCF; CNV/SV absence is NOT established'
 return rows,qc

def add_clinvar(rows):
 wanted=defaultdict(list)
 for i,r in enumerate(rows):wanted[(r['chrom'],r['pos'],r['ref'],r['alt'])].append(i)
 if not wanted:return
 p=RES/'clinvar.vcf.gz'
 # Cache only exact alleles requested by previous intake runs; full ClinVar remains the source of truth.
 cpath=RES/'clinvar_exact_cache.json';cache=load(cpath) if cpath.exists() else {}
 unresolved=set(wanted)
 for key,idxs in wanted.items():
  sk=':'.join(map(str,key))
  if sk in cache:
   for i in idxs:rows[i]['clinvar']=cache[sk]
   unresolved.discard(key)
 if not unresolved:return
 with gzip.open(p,'rt',encoding='utf-8',errors='replace') as f:
  for line in f:
   if line.startswith('#'):continue
   x=line.rstrip('\n').split('\t')
   if len(x)<8:continue
   chrom=chrom_norm(x[0]);pos=int(x[1]);ref=x[3];alts=x[4].split(',');info=_parse_info(x[7])
   for alt in alts:
    key=(chrom,pos,ref,alt)
    if key not in wanted:continue
    cv={'variation_id':x[2] if x[2]!='.' else None,'CLNSIG':info.get('CLNSIG'),'CLNREVSTAT':info.get('CLNREVSTAT'),
        'CLNDN':info.get('CLNDN'),'GENEINFO':info.get('GENEINFO'),'CLNVC':info.get('CLNVC'),'source':'ClinVar_GRCh38_local_vcf'}
    for i in wanted[key]:rows[i]['clinvar']=cv

def _clinvar_weight(cv):
 # Clinical-significance consensus and review status must be accounted for.
 # No automatic ACMG diagnosis; conflicts are never positive evidence.
 return case_science.clinvar_interpretation(cv)['weight']

def rank_genes(rows,features):
 gene_terms,parents=parse_hpo_resources();by=defaultdict(list)
 for r in rows:
  # Exact benign/likely-benign ClinVar evidence is a veto for that allele in rare-disease ranking,
  # but the allele remains in the genome-wide ledger.
  if r['best_annotation']['severity']>=4 and _clinvar_weight(r.get('clinvar'))>-10:
   by[r['gene'].upper()].append(r)
 ranked=[]
 for gene,vs in by.items():
  ph=phenotype_gene_score(gene,features,gene_terms,parents)
  best=max(v['best_annotation']['severity'] for v in vs)
  q=0
  for v in vs:
   if v['filter']=='PASS':q=max(q,2)
   if (v.get('GQ') or 0)>=60:q=max(q,3)
   if v.get('allele_balance') is not None and 0.2<=v['allele_balance']<=0.8:q=max(q,4)
  cv=max([_clinvar_weight(v.get('clinvar')) for v in vs] or [0])
  def usable(v):
   if v['best_annotation']['severity']<6 or v['filter']!='PASS' or _clinvar_weight(v.get('clinvar'))<=-10:return False
   if v['best_annotation'].get('reference_sequence_mismatch'):return False
   if v.get('GQ') is not None and v['GQ']<20:return False
   if v['GT'] in ('0/1','1/0','0|1','1|0') and v.get('allele_balance') is not None and not (0.15<=v['allele_balance']<=0.85):return False
   return True
  damaging=[v for v in vs if usable(v)]
  het=sum(v['GT'] in ('0/1','1/0','0|1','1|0') for v in damaging);hom=sum(v['GT'] in ('1/1','1|1') for v in damaging)
  artifact_review=len(damaging)>4 or len(vs)>8
  phase_eval=case_science.phase_assessment(damaging)
  inheritance_boost=0 if artifact_review else (7 if hom else 6 if phase_eval['trans_alleles_established'] else 0)
  artifact_penalty=12 if artifact_review else 0
  score=best+q+cv+ph['score']+inheritance_boost-artifact_penalty
  ranked.append({'gene':gene,'score':round(score,3),'variant_count':len(vs),'damaging_pass_count':len(damaging),
    'mapping_artifact_review':artifact_review,'artifact_penalty':artifact_penalty,
    'heterozygous_damaging_count':het,'homozygous_damaging_count':hom,'inheritance_boost':inheritance_boost,
    'allelic_phase_evidence':phase_eval,
    'phenotype_score':ph['score'],'phenotype_matches':ph['matched'],'clinvar_weight':cv,
    'best_severity':best,'variant_keys':[f"{v['chrom']}:{v['pos']}:{v['ref']}:{v['alt']}" for v in vs]})
 ranked.sort(key=lambda x:(-x['score'],-x['damaging_pass_count'],x['gene']))
 return ranked

def build_case_model(vcf_path,clinical_path,outdir,max_primary_genes=3):
 out=Path(outdir);out.mkdir(parents=True,exist_ok=True);res=ensure_resources();mane=load_mane_models();models=mane['transcripts'];bins=build_interval_bins(models)
 clinical=parse_clinical_docx(clinical_path)
 variants,qc=scan_vcf(Path(vcf_path),models,bins,out/'genome_wide_variants.preclinvar.jsonl')
 if qc['build']!='GRCh38':raise RuntimeError('V28_INTAKE_REQUIRES_GRCH38_OR_EXPLICIT_LIFTOVER')
 add_clinvar(variants);jsonl(out/'genome_wide_variants.jsonl',variants)
 ranked=rank_genes(variants,clinical['features']);save(out/'genome_wide_genes.json',ranked)
 # Select genes conservatively: strongest phenotype+variant evidence; preserve full differential separately.
 selected=[]
 if ranked:
  top=ranked[0]['score']
  for g in ranked:
   if len(selected)>=max_primary_genes:break
   # retain close alternatives rather than pretending diagnostic certainty
   if not selected or g['score']>=max(top-6,18):selected.append(g['gene'])
 bykey={(r['chrom'],r['pos'],r['ref'],r['alt']):r for r in variants}
 genetic=[]
 for rank,gene in enumerate(selected,1):
  gr=next(x for x in ranked if x['gene']==gene)
  for key in gr['variant_keys']:
   c,pos,ref,alt=key.split(':');v=bykey[(c,int(pos),ref,alt)]
   if v['best_annotation']['severity']<4:continue
   a=v['best_annotation'];refseq=(a.get('refseq_rna') or [None])[0]
   hgvsc=(refseq+':'+a['hgvsc_local']) if refseq and a.get('hgvsc_local') else ((a['transcript_id']+'.'+str(a.get('transcript_version'))+':'+a['hgvsc_local']) if a.get('hgvsc_local') else None)
   hgvsp=a.get('hgvsp_local')
   genetic.append({'gene':gene,'selected_gene_rank':rank,'selected_for_downstream':True,'chrom':v['chrom'],'pos':v['pos'],'ref':v['ref'],'alt':v['alt'],
    'GT':v['GT'],'AD':v['AD'],'DP':v['DP'],'GQ':v['GQ'],'FILTER':v['filter'],'allele_balance':v['allele_balance'],
    'phase_tags':{'PS':v['PS'],'PID':v['PID'],'PGT':v['PGT']},'genotype_class':'heterozygous_unphased' if v['GT'] in ('0/1','1/0') and not any((v['PS'],v['PID'],v['PGT'])) else v['GT'],
    'transcript_id':a['transcript_id'],'transcript_version':a.get('transcript_version'),'refseq_rna':a.get('refseq_rna',[]),'mane_tags':a.get('mane_tags',[]),
    'hgvsc':hgvsc,'hgvsp':hgvsp,'consequence':a['consequence'],'impact_score':a['severity'],'clinvar':v.get('clinvar'),
    'interpretation':'genome_wide_candidate_not_functionally_validated','source':'raw_vcf_local_annotation'})
 gene_allelic_phase={}
 for g in selected:
  vs=[x for x in genetic if x['gene']==g]
  gene_allelic_phase[g]=case_science.phase_assessment(vs)
 phase_resolved=bool(gene_allelic_phase) and all(v['state'] in ('cis_supported','trans_supported_by_vcf','homozygous_observed') for v in gene_allelic_phase.values())
 case={
  'schema_version':2,'case_model_version':'V28_RAW_INTAKE_2','case_id':'CASE_'+sha(Path(vcf_path))[:12],'origin':'raw_patient_inputs',
  'build':qc['build'],'genetic_candidates':genetic,'primary_genes':selected,'genome_wide_candidate_genes':ranked[:50],
  'genome_wide_analysis':{'performed':True,'mane_resource':'GENCODE 50 MANE Select + MANE Plus Clinical','candidate_variant_ledger':'genome_wide_variants.jsonl',
    'ranked_gene_ledger':'genome_wide_genes.json','vcf_qc':qc,'limitations':[
      'VCF callset analysis does not prove absence of variants in poorly covered regions','No BAM/CRAM read-level validation was performed',
      'CNV/SV absence is not inferred unless represented by the supplied VCF','Population-frequency filtering requires a separate local population database and is not currently used',
      'Lightweight local consequence engine generates coding/splice classes and SNV HGVS; complex indel HGVS is not generated']},
  'phase':'vcf_phase_supported' if phase_resolved else 'unresolved',
  'phase_evidence':gene_allelic_phase,'phase_assessment_version':case_science.GATE_VERSION,
  'compound_heterozygosity_confirmed':False,'residual_protein_status':'unknown','expression_measured_in_this_case':False,
  'clinical_phenotype_source_present':clinical['source_present'],'phenotype_features':clinical['features'],'clinical_narrative':clinical['narrative'],
  'phenotype_review_domains':clinical.get('phenotype_review_domains',{}),'current_medication_and_renal_function':'not_provided_unless_explicitly_present',
  'preferred_drugs':[],'excluded_drug_names':[],'public_export_allowed':False,'api_payload_allowed':False,
  'official_source_hashes':[{'file':Path(vcf_path).name,'sha256':sha(Path(vcf_path))}]+([{'file':Path(clinical_path).name,'sha256':sha(Path(clinical_path))}] if clinical_path else []),
  'resource_hashes':res['hashes'],'patient_data_external':False,'diagnostic_certainty':'hypothesis_ranked_not_clinically_confirmed'}
 save(out/'case_model_v2.json',case)
 manifest={'version':'V28_INTAKE_MANIFEST_2','created_utc':now(),'inputs':case['official_source_hashes'],'resources':res,
   'outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='intake_manifest.json'},
   'privacy':{'raw_vcf_external':False,'variant_strings_external':False,'clinical_text_external':False}}
 save(out/'intake_manifest.json',manifest)
 return case,manifest

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--vcf',required=True);ap.add_argument('--clinical');ap.add_argument('--out',required=True)
 a=ap.parse_args();case,_=build_case_model(a.vcf,a.clinical,a.out)
 print(json.dumps({'primary_genes':case['primary_genes'],'top10':[(x['gene'],x['score']) for x in case['genome_wide_candidate_genes'][:10]],'genetic_candidates':len(case['genetic_candidates']),'hpo':len(case['phenotype_features']),'build':case['build']},indent=2))
