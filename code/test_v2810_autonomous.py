"""Independent P0/P1 tests of autonomous V2.8.5 drug-report process."""
import json,unittest,sys,zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
from collections import defaultdict
from unittest.mock import patch
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from v2810 import case_science,intake,science,deep,screen,evidence,reporting,pipeline
from v2810.core import load

def variant(pos,gt,ps=None,consequence='missense_variant',hgvsp='p.Ala1Thr'):
 return {'chrom':'1','pos':pos,'ref':'A','alt':'G','gene':'DEMO1','GT':gt,
  'phase_tags':{'PS':ps} if ps else {},'filter':'PASS','GQ':99,
  'allele_balance':0.48,'clinvar':None,
  'hgvsp':hgvsp,'consequence':consequence,
  'best_annotation':{'severity':12,'reference_sequence_mismatch':False}}
def route(rid='R1',family='direct_gene',gene='DEMO1'):
 return {'route_id':rid,'family':family,'target_concept':gene,
  'targets':[{'gene':gene,'accession':'P00001'}],
  'desired_action':'increase','variable':'activity','mechanism_concept':'restore DEMO1 abundance',
  'score':12,'state':'conditional'}
def entity(name,signals=None):
 signals=signals or []
 return {'canonical_label':name,'entity_id':'SYN:'+name,
  'chembl_parent_ids':[],'global_approved_as_single':True,
  'combo_only':False,'modality':'small_molecule','approved_jurisdictions':['Fictionland'],
  'signals':signals,'signal_count':len(signals),
  'status':'SIGNAL' if signals else 'NO_SIGNAL',
  'screened':True,'patient_data_external':False}

class V286GeneticAudits(unittest.TestCase):
 def test_conflicting_clinvar_does_not_score_pathogenic(self):
  self.assertEqual(case_science.clinvar_interpretation({'CLNSIG':'Pathogenic|Benign'})['weight'],0)
  self.assertEqual(case_science.clinvar_interpretation({'CLNSIG':'conflicting_classifications_of_pathogenicity'})['weight'],0)
  self.assertEqual(intake._clinvar_weight({'CLNSIG':'Pathogenic|Benign'}),0)

 def test_clinvar_without_review_not_automatically_pathogenic(self):
  self.assertEqual(case_science.clinvar_interpretation({'CLNSIG':'Pathogenic','CLNREVSTAT':'no_assertion_criteria_provided'})['weight'],0)

 def test_absent_hpo_not_counted(self):
  self.assertEqual(case_science.presence_from_row(['Seizures absent','Seizures','HP:0001250','explicitly denied']),'reported_absent')
  self.assertEqual(case_science.presence_from_row(['Short stature','Short stature','HP:0004322','present']),'reported_present')
  with patch.object(intake,'parse_hpo_resources',return_value=({'DEMO1':{'HP:0001250'}},defaultdict(set))):
   out=intake.phenotype_gene_score('DEMO1',[{'HPO_id':'HP:0001250','presence':'reported_absent','subject':'proband'}],
     {'DEMO1':{'HP:0001250'}},defaultdict(set))
  self.assertEqual(out['score'],0)

 def test_clinical_docx_parser_explicit_absence(self):
  NS="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
  w=lambda x:'{'+NS+'}'+x
  doc=ET.Element(w('document')); body=ET.SubElement(doc,w('body')); tab=ET.SubElement(body,w('tbl'))
  for parts in [['Seizures absent','Seizures','HP:0001250','denied'],
                ['Short stature','Short stature','HP:0004322','present']]:
   tr=ET.SubElement(tab,w('tr'))
   for val in parts:
    c=ET.SubElement(tr,w('tc'));p=ET.SubElement(c,w('p'));r=ET.SubElement(p,w('r'));ET.SubElement(r,w('t')).text=val
  with TemporaryDirectory() as d:
   p=Path(d)/'fic.docx'
   with zipfile.ZipFile(p,'w') as z:z.writestr('word/document.xml',ET.tostring(doc,encoding='utf-8'))
   features=intake.parse_clinical_docx(p)['features']
  self.assertEqual([x['presence'] for x in features],['reported_absent','reported_present'])

 def test_cis_variants_no_recessive_bonus(self):
  v=[variant(100,'0|1','phase1'),variant(200,'0|1','phase1')]
  with patch.object(intake,'parse_hpo_resources',return_value=(defaultdict(set),defaultdict(set))):
   ranked=intake.rank_genes(v,[])
  self.assertEqual(ranked[0]['inheritance_boost'],0)
  self.assertEqual(ranked[0]['allelic_phase_evidence']['state'],'cis_supported')

 def test_trans_variants_bonus_without_confirming_clinical_diagnosis(self):
  v=[variant(100,'0|1','phase1'),variant(200,'1|0','phase1')]
  with patch.object(intake,'parse_hpo_resources',return_value=(defaultdict(set),defaultdict(set))):
   ranked=intake.rank_genes(v,[])
  self.assertEqual(ranked[0]['inheritance_boost'],6)
  self.assertEqual(ranked[0]['allelic_phase_evidence']['state'],'trans_supported_by_vcf')

 def test_allele_specific_objectives_are_separate(self):
  private={'genes':{'DEMO1':[{'gene':'DEMO1','hgvsp':'p.Ala20Thr',
     'GT':'0/1','consequence':'missense_variant'},
    {'gene':'DEMO1','hgvsp':'p.Arg70Ter','GT':'0/1','consequence':'stop_gained'}]},
    'phenotype_features':[],'diagnostic_certainty':'unconfirmed'}
  view=case_science.build_case_science(private)
  profiles=view['primary_gene_hypotheses'][0]['variant_effect_models']
  self.assertNotEqual(profiles[0]['variant_ref'],profiles[1]['variant_ref'])
  self.assertNotEqual(profiles[0]['effect_class'],profiles[1]['effect_class'])
  self.assertFalse(any(x['mutation_specific_therapy_direction_known'] for x in profiles))
  self.assertFalse(view['genomic_allele_strings_exported'])

class V286PharmaAudits(unittest.TestCase):
 def test_real_functional_chEMBL_assay_can_enter_b(self):
  r=route();g=science.adjudicate_routes([r])[0][0]
  activity={'assay_chembl_id':'AA','molecule_chembl_id':'M1',
    'parent_molecule_chembl_id':'M1','target_chembl_id':'T1',
    'standard_relation':'=','standard_units':'nM','standard_type':'EC50',
    'standard_value':'50','canonical_smiles':None}
  assay={'assay_chembl_id':'AA','confidence_score':9,'assay_type':'F',
    'description':'Direct activation of DEMO1 protein functional activity'}
  def fetch(endpoint,*a,**kw):return [activity] if endpoint=='activity' else [assay]
  with TemporaryDirectory() as td,patch.object(screen,'chembl_paged',side_effect=fetch):
   records=screen.route_activities([g],{'P00001':['T1']},Path(td))
  self.assertEqual(len(records),1)
  self.assertTrue(records[0]['direct_target_assay_verified'])
  self.assertEqual(records[0]['role'],'compatible')
  sig={'route_id':'R1','provider':'ChEMBL_measured_activity',
    'value_nM':50,'direct_target_assay_verified':True}
  grade=deep.prioritize([entity('MEDICINE',[sig])],[g])
  self.assertEqual(grade[0]['tier'],'B')
  self.assertFalse(grade[0]['promotion_eligible'])

 def test_binding_ki_never_activation(self):
  self.assertEqual(case_science.clinvar_interpretation({})['weight'],0)
  g=science.adjudicate_routes([route()])[0][0]
  sig={'route_id':'R1','provider':'ChEMBL_measured_activity','value_nM':5,
    'direct_target_assay_verified':False}
  self.assertEqual(deep.prioritize([entity('MEDICINE',[sig])],[g])[0]['tier'],'D')

 def test_chEMBL_unavailable_not_implicit_zero_effect(self):
  g=route()
  with TemporaryDirectory() as td,patch.object(screen,'chembl_paged',side_effect=TimeoutError("synthetic outage")):
   with self.assertRaisesRegex(RuntimeError,'V286_CHEMBL_ACTIVITY_UNAVAILABLE'):
    screen.route_activities([g],{'P00001':['T1']},Path(td))
   self.assertFalse((Path(td)/'route_activities.json').exists())

 def test_source_quote_two_publications_can_be_provisional_a(self):
  base=deep.prioritize([entity('EXAMPLE COMPOUND',[{'route_id':'R1','provider':'ChEMBL_MoA_SIGNOR','depth':2}])],
      science.adjudicate_routes([route()])[0])[0]
  base['scientific_group_id']='G1'
  quote="Example compound rescues the mitotic checkpoint and improves DEMO1 functional activity after exposure."
  research={'scientific_group_id':'G1','source_route_id':'R1','ai_evaluation':{
    'adjudication':'SUPPORTED','exposure_reported':True,'direct_route_modulation':True,
    'same_drug_functional_rescue':True,'counterevidence_present':False,
    'support_quotes':[{'source_id':'PMC1234567','quote':quote},
                      {'source_id':'PMC2234567','quote':quote}]}}
  rows,events=evidence.enrich_pairs_from_literature([base],[research],{'DEMO1':[]})
  self.assertEqual(rows[0]['tier'],'A')
  self.assertFalse(rows[0]['promotion_eligible'])
  self.assertEqual(len(events),1)
  research['ai_evaluation']['support_quotes']=research['ai_evaluation']['support_quotes'][:1]
  one,_=evidence.enrich_pairs_from_literature([base],[research],{'DEMO1':[]})
  self.assertEqual(one[0]['tier'],'B')
  research['ai_evaluation']['counterevidence_present']=True
  contradicted,_=evidence.enrich_pairs_from_literature([base],[research],{'DEMO1':[]})
  self.assertEqual(contradicted[0]['tier'],'D')

 def test_report_describes_all_tiers_and_unknowns(self):
  private={'genes':{},'phenotype_features':[]}
  grouped=[entity(name) for name in ['EXAMPLE A','EXAMPLE B','EXAMPLE C','EXAMPLE D','EXAMPLE U']]
  for e in grouped:e['scientific_group_id']=e['canonical_label']
  pairs=[]
  for grade,name in [('A','EXAMPLE A'),('B','EXAMPLE B'),('C','EXAMPLE C'),('D','EXAMPLE D')]:
   p={'ingredient':name,'scientific_group_id':name,'entity_id':'SYN:'+name,
     'route_id':name,'route':route(name),'tier':grade,'score':5,
     'providers':['ChEMBL_MoA_SIGNOR'],'signals':[{'provider':'ChEMBL_MoA_SIGNOR','depth':2}],
     'scientific_gate':{}}
   pairs.append(p)
  with TemporaryDirectory() as td:
   summary,top=reporting.create(td,grouped,pairs,private,
       {'primary_gene_hypotheses':[],'phenotype_counts':{'present':0}},registry_index={})
   self.assertEqual(summary['evidence_tier_counts'],{x:1 for x in 'ABCDU'})
   markdown=(Path(td)/'CANDIDATE_REPORT.md').read_text(encoding='utf-8')
   for grade in 'ABCDU':self.assertIn('**'+grade+'**',markdown)
   self.assertIn('Limitations',markdown)
   self.assertIn('not prescriptions',markdown)
   self.assertEqual(len((Path(td)/'CANDIDATES_ALL.csv').read_text(encoding='utf-8').splitlines()),6)

 def test_full_pipeline_synthetic_case_to_report(self):
  with TemporaryDirectory() as td:
   root=Path(td)
   case=root/'fictional_case.json'
   case.write_text(json.dumps({'case_id':'DEMO_ONLY','primary_genes':['DEMO1'],
     'genetic_candidates':[{'gene':'DEMO1','hgvsp':'p.Gly4Ser',
       'GT':'0/1','phase_tags':{},'consequence':'missense_variant',
       'selected_for_downstream':True}],
     'phenotype_features':[],'phase':'unresolved','residual_protein_status':'unknown'}),encoding='utf-8')
   entities=[entity('DEMO DRUG '+str(n)) for n in range(2414)]
   activated={**entities[0],'signals':[{'provider':'ChEMBL_MoA_SIGNOR','route_id':'V28R001','depth':2}],
      'status':'SIGNAL','signal_count':1}
   ledger=[activated]+entities[1:]
   info={'screened':2414,'status_counts':{'SIGNAL':1,'NO_SIGNAL':2413},
      'signal_entities':1,'no_entity_skipped':True,'enrichr_library_mode':'disabled'}
   with patch.object(pipeline.master,'build',return_value=(entities,[],{'entities':2414})),\
       patch.object(pipeline.routes,'base_routes',return_value=[route()]),\
       patch.object(pipeline.screen,'universal_screen',return_value=(ledger,info)):
    result=pipeline.run(case_path=case,outdir=root/'run',use_ai=False,functional_review_count=1)
   self.assertEqual(result['research_evidence_tier_counts'].get('D',0),0)
   self.assertEqual(result['research_evidence_tier_counts']['U'],2414)
   self.assertTrue((root/'run/research_reports/CANDIDATE_REPORT.md').exists())
   self.assertTrue((root/'run/research_reports/CANDIDATES_ALL.csv').exists())
   self.assertTrue((root/'run/PATIENT_VARIANT_MECHANISM_MODEL.json').exists())
   self.assertFalse(result['clinical_safety_clearance'])
   self.assertTrue(load(root/'run/STATUS.json')['completed'])


class V286RawVCFIntegration(unittest.TestCase):
 def test_true_vcf_input_to_research_report_entire_pipeline_without_network(self):
  import gzip
  with TemporaryDirectory() as td:
   root=Path(td)
   vcf=root/'fictional_variants.vcf.gz'
   header=("##fileformat=VCFv4.2\n##reference=GRCh38\n"
     "##contig=<ID=1,length=248956422>\n##contig=<ID=15,length=101991189>\n"
     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
   rows=("1\t96\t.\tA\tG\t99\tPASS\t.\tGT:PS:AD:DP:GQ\t0|1:77:22,19:41:90\n"
     "1\t99\t.\tA\tG\t99\tPASS\t.\tGT:PS:AD:DP:GQ\t1|0:77:22,19:41:90\n")
   with gzip.open(vcf,'wt') as f:f.write(header+rows)
   model={'TX1':{'transcript_id':'TX1','transcript_version':'1',
     'gene_name':'DEMO1','gene_id':'ENSG000TEST',
     'chrom':'1','strand':'+','tags':['MANE_Select'],'refseq_rna':[],
     'exons':[{'start':90,'end':120,'cdna_start0':0}],
     'cds':[{'start':90,'end':120}],
     'cds_cdna_start0':0,'cds_cdna_end0':31,'sequence':'A'*31}}
   entities=[entity('DEMO DRUG '+str(n)) for n in range(2414)]
   ledger=[{**entities[0],'signals':[{'provider':'ChEMBL_MoA_SIGNOR',
     'route_id':'V28R001','depth':2}],'status':'SIGNAL','signal_count':1}]+entities[1:]
   screening={'screened':2414,'status_counts':{'SIGNAL':1,'NO_SIGNAL':2413},
       'signal_entities':1,'no_entity_skipped':True,'enrichr_library_mode':'disabled'}
   with patch.object(pipeline.intake,'ensure_resources',return_value={'hashes':{},'resource_dir':td}),\
        patch.object(pipeline.intake,'load_mane_models',return_value={'transcripts':model}),\
        patch.object(pipeline.intake,'add_clinvar',return_value=None),\
        patch.object(pipeline.intake,'parse_hpo_resources',return_value=(defaultdict(set),defaultdict(set))),\
        patch.object(pipeline.master,'build',return_value=(entities,[],{'entities':2414})),\
        patch.object(pipeline.routes,'base_routes',return_value=[route()]),\
        patch.object(pipeline.screen,'universal_screen',return_value=(ledger,screening)):
    result=pipeline.run(vcf_path=vcf,outdir=root/'end_to_end',
      use_ai=False,functional_review_count=1)
   self.assertTrue(result['raw_intake_performed'])
   self.assertEqual(result['master_entities'],2414)
   self.assertFalse(result['clinical_recommendation'])
   case=load(root/'end_to_end/intake/case_model_v2.json')
   self.assertEqual(case['primary_genes'],['DEMO1'])
   self.assertEqual(case['phase_evidence']['DEMO1']['state'],'trans_supported_by_vcf')
   self.assertFalse(case['compound_heterozygosity_confirmed'])
   self.assertEqual(len(case['genetic_candidates']),2)
   self.assertTrue((root/'end_to_end/research_reports/CANDIDATE_REPORT.md').exists())
   self.assertTrue(load(root/'end_to_end/STATUS.json')['completed'])

if __name__=='__main__': unittest.main()
