"""Independent V2.8.8 adversarial and synthetic multi-gene end-to-end validation.

Zero network, zero patient records, zero drug recommendations.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from collections import Counter
import json,unittest
from v2810 import expert_genomics,expert_protein,expert_triage,routes,science,evidence

def vc(gt='0/1',clinvar=None,gene='DEMO1'):
    return {'gene':gene,'chrom':'1','pos':123,'ref':'A','alt':'G',
     'GT':gt,'DP':38,'GQ':90,'FILTER':'PASS','allele_balance':1.0 if gt=='1/1' else 0.46,
     'phase_tags':{},'consequence':'missense_variant',
     'hgvsp':'p.Arg123Cys','clinvar':clinvar}
def demo_dossier(gene="DEMO1"):
    return {"genes":{gene:{"gene_mechanism_confidence":"SUPPORTED_VARIANT_BUT_PHASE_OR_SECOND_ALLELE_UNRESOLVED",
        "phase_state":"unresolved","protein_rescue_endpoint_measured":False,
        "variants":[{"consequence":"missense_variant"}]}}}
def pair(name="CANDIDATE",gene="DEMO1",sig=None,family="direct_gene",routeid="R1",desc=None,tier="D",score=0.15):
    return {"ingredient":name,"scientific_group_id":"G:"+name,"route_id":routeid,
        "route":{"case_gene":gene,"route_id":routeid,"family":family,"target_concept":gene,
        "targets":[{"gene":gene,"accession":"P00001"}],
        "desired_action":"increase","variable":"activity","state":"conditional",
        "mechanism_concept":desc or "increase "+gene+" protein function",
        "scientific_gate":{"hard_hold":False,"issues":[]}},
        "member_labels":[name],"tier":tier,"score":score,
        "signals":sig or [],
        "patient_first_prediction":{"primary_gene_hypothesis":gene,"prediction_class":"HYPOTHESIS_ONLY"}}

class LocalClinicalEvidence(unittest.TestCase):
 def test_conflicting_or_benign_never_causally_established(self):
    v=vc('1/1',{'CLNSIG':'Benign/Likely_benign','CLNREVSTAT':'criteria_provided,_multiple_submitters,_no_conflicts'})
    item=expert_genomics.inspect_allele(v)
    self.assertEqual(item['clinvar_assertion'],'reported_benign')
    self.assertFalse(item['functional_effect_observed_in_patient'])
    self.assertIn('clinical_assertion_not_supportive_of_pathogenicity',item['gaps'])
 def test_quality_failure_and_not_pop_AF_from_generic_AF(self):
    v={**vc(),'DP':3,'GQ':12,'allele_balance':0.05,'AF':0.01}
    e=expert_genomics.inspect_allele(v)
    self.assertFalse(e['qc_pass_for_exploration'])
    self.assertIsNone(e['population_af'])
    self.assertIn('source_verified_population_af_unavailable',e['gaps'])
 def test_population_annotation_source_and_build_controls(self):
    with TemporaryDirectory() as td:
       p=Path(td)/'pop.json'
       data={'schema':'V288_CURATED_POPULATION_AF_1','assembly':'GRCh38','entries':[
         {'variant_key':'1:123:A:G','source_reference':'SYNTHETIC_GNOMAD_TEST','population':'test',
          'allele_frequency':0.03,'allele_number':500,'source_checked':True}]}
       p.write_text(json.dumps(data),encoding='utf-8')
       pop,summary=expert_genomics.load_local_population(p)
       self.assertEqual(summary['entries'],1)
       self.assertIn('allele_frequency_above_1pct_requires_disease_model_review',
        expert_genomics.inspect_allele(vc(),pop)['gaps'])
       data['assembly']='GRCh37';p.write_text(json.dumps(data))
       with self.assertRaisesRegex(ValueError,'BUILD_MISMATCH'):expert_genomics.load_local_population(p)
 def test_prepared_dossier_contains_no_allelic_strings(self):
    private={'genes':{'DEMO1':[vc()]},'phenotype_features':[],'genome_wide_candidate_genes':[]}
    doc=expert_genomics.build(private)
    s=json.dumps(doc)
    self.assertNotIn('"pos"',s)
    self.assertNotIn('"ref"',s)
    self.assertFalse(doc['clinical_diagnosis_confirmed'])
    self.assertEqual(doc['genes']['DEMO1']['variants'][0]['acmg_amp_classification'],
     'NOT_ESTABLISHED_BY_AUTOMATION')
 def test_uniprot_features_map_residue_but_not_function(self):
    with TemporaryDirectory() as td:
      p=Path(td)/'uniprot';p.mkdir()
      payload={'primaryAccession':'P12345','organism':{'taxonId':9606},
       'genes':[{'geneName':{'value':'DEMO1'}}],
       'features':[{'type':'Domain','description':'Catalytic region',
                   'location':{'start':{'value':101},'end':{'value':160}}}]}
      (p/'DEMO1_P12345.json').write_text(json.dumps(payload),encoding='utf-8')
      private={'genes':{'DEMO1':[vc()]}}
      x=expert_protein.build_from_private(private,{},td)
      features=x['genes']['DEMO1']['annotated_protein_positions'][0]
      self.assertEqual(features['feature_status'],'overlap_observed')
      self.assertEqual(features['feature_hits'][0]['feature_type'],'Domain')
      self.assertFalse(features['domain_overlap_is_functional_proof'])
      self.assertFalse(features['transcript_isoform_to_uniprot_alignment_validated'])

class RouteAndEvidenceProvenance(unittest.TestCase):
 def test_every_base_route_has_disease_gene(self):
    private={'genes':{'DEMO1':[], 'DEMO2':[]},'variant_classes':[]}
    with patch.object(routes,'signor_edges',return_value=[]),\
         patch.object(routes,'accession_for_gene',side_effect=lambda g,d:'P'+g):
      rows=routes.base_routes(private,Path('.'))
    self.assertEqual({x['case_gene'] for x in rows},{'DEMO1','DEMO2'})
    for x in rows:self.assertEqual(x['case_gene'],x['target_concept'])
 def test_same_target_alternative_genes_not_merged(self):
    candidates=[]
    for g in ('DEMO1','DEMO2'):
      candidates.append({'case_gene':g,'target_concept':'BRD4',
       'desired_action':'decrease','variable':'activity','target_gene':'BRD4',
       'target_accession':'P00004','source_id':'PMC1234567',
       'evidence':[{'quote':'synthetic'}]})
    merged=routes.merge_discovered([],candidates)
    self.assertEqual(len(merged),2)
    self.assertEqual({x['case_gene'] for x in merged},{'DEMO1','DEMO2'})
 def test_literature_lookup_uses_corresponding_gene(self):
    rows=[pair('C1','DEMO1',sig=[{'provider':'perturbational_transcriptomics'}],routeid='G1'),
          pair('C2','DEMO2',sig=[{'provider':'perturbational_transcriptomics'}],routeid='G2')]
    visited=[]
    def fake(gene,*args,**kw):visited.append(gene);return []
    with TemporaryDirectory() as td,\
      patch.object(science.ai if hasattr(science,'ai') else evidence,'_nonexistent',create=True),\
      patch('v2810.ai.literature',side_effect=fake):
       results=science.exploratory_functional_evidence_search(rows,'WRONG',td,
           count=2,use_ai=True,registry_index={})
    self.assertEqual(set(visited),{'DEMO1','DEMO2'})
    self.assertEqual(len(results),2)
 def test_source_check_uses_route_gene_not_alphabetical_first(self):
    doc=pair("EXAMPLE COMPOUND","DEMO2",sig=[{'provider':'perturbational_transcriptomics'}])
    doc['scientific_group_id']='G222'
    quote="Example compound restores DEMO2 mitotic checkpoint function and improves chromosome segregation."
    review={'scientific_group_id':'G222','source_route_id':'R1',
     'ai_evaluation':{'adjudication':'SUPPORTED','exposure_reported':False,
        'direct_route_modulation':True,'same_drug_functional_rescue':True,
        'counterevidence_present':False,
        'support_quotes':[{'source_id':'PMC1234567','quote':quote}]}}
    r,changed=evidence.enrich_pairs_from_literature([doc],[review],{'DEMO1':[],'DEMO2':[]})
    self.assertEqual(r[0]['tier'],'B')
    self.assertEqual(r[0]['functional_evidence']['validated_source_ids'],['PMC1234567'])

class MechanismSpecificity(unittest.TestCase):
 def test_expression_only_is_not_a_patient_drug_efficacy_prediction(self):
    sig={'provider':'perturbational_transcriptomics','route_id':'R1','support_count':14}
    x=expert_triage.assessment(pair(sig=[sig],score=22.5),demo_dossier())
    self.assertEqual(x['investigation_bucket'],'T3')
    self.assertIsNone(x['calibrated_efficacy_probability'])
    self.assertTrue(x['mechanistic_compatibility_not_computable'])
    self.assertIn('RNA_perturbation_or_similarity_does_not_establish_direction_of_protein_effect',x['hypothesis_blockers'])
 def test_direct_measured_target_assay_is_separate_from_disease_rescue(self):
    sig={'provider':'ChEMBL_measured_activity','route_id':'R1','value_nM':40,
      'direct_target_assay_verified':True,'assay_id':'CHEMBL-TEST'}
    x=expert_triage.assessment(pair(sig=[sig],score=2),demo_dossier())
    self.assertEqual(x['investigation_bucket'],'T1')
    self.assertTrue(x['positive_experimental_disease_rescue_not_established'])
 def test_direct_moa_not_mistaken_for_same_drug_functional_rescue(self):
    x=expert_triage.assessment(pair(sig=[{'provider':'ChEMBL_MoA_SIGNOR','depth':0}]),demo_dossier())
    self.assertEqual(x['investigation_bucket'],'T2')
    self.assertTrue(x['positive_experimental_disease_rescue_not_established'])
 def test_cancer_suppression_of_primary_gene_holds_prediction(self):
    x=expert_triage.assessment(pair('CANCER DRUG','ERCC2',
        sig=[{'provider':'ChEMBL_MoA_SIGNOR','depth':1}],routeid='R1',
        family='literature_discovered',
        desc="BRD4 inhibition downregulates ERCC2 expression in carcinoma cells"),
        demo_dossier('ERCC2'))
    self.assertEqual(x['investigation_bucket'],'HOLD')
    self.assertIn('source_describes_suppression_of_the_gene_we_aim_to_restore',x['hypothesis_blockers'])
 def test_held_route_and_wrong_patient_gene_are_separately_identified(self):
    p=pair('A','DEMO1',sig=[{'provider':'ChEMBL_measured_activity','direct_target_assay_verified':True,'value_nM':80}])
    p['route']['state']='hold'
    self.assertEqual(expert_triage.assessment(p,demo_dossier())['investigation_bucket'],'HOLD')
    p['route']['state']='conditional'
    self.assertEqual(expert_triage.assessment(p,demo_dossier('OTHER'))['investigation_bucket'],'T4')
 def test_low_evidence_numerical_score_does_not_determine_research_bucket(self):
    x=pair("HIGH NUMBER",sig=[{'provider':'perturbational_transcriptomics'}],score=60)
    x['patient_first_score']=45
    y=pair("LOW NUMBER",sig=[{'provider':'ChEMBL_measured_activity','value_nM':90,
                             'direct_target_assay_verified':True}],score=0.1)
    y['patient_first_score']=1
    records,summary=expert_triage.order_pairs([x,y],demo_dossier())
    self.assertEqual(records[0]['ingredient'],"LOW NUMBER")
    self.assertEqual(records[0]['expert_investigation_bucket'],"T1")
    self.assertEqual(records[1]['expert_investigation_bucket'],"T3")
 def test_safety_review_diversifies_routes_and_disallows_worse_links(self):
    rows=[]
    for i in range(9):
      p=pair("C"+str(i),routeid="R"+str(i%3),sig=[{'provider':'perturbational_transcriptomics'}])
      p['scientific_group_id']="G"+str(i)
      p['expert_investigation_bucket']="T3"
      p['expert_triage']={'investigation_bucket':'T3'}
      rows.append(p)
    chosen=expert_triage.safety_review_queue(rows,{"C0"},4)
    self.assertEqual(len(chosen),4)
    self.assertNotIn("C0",chosen)
 def test_gene_suppression_language_is_source_direction_specific(self):
    self.assertTrue(expert_triage.is_gene_suppression_mechanism(
      "HDAC3 inhibition downregulates ERCC2 expression in cancer cells","ERCC2"))
    self.assertTrue(expert_triage.is_gene_suppression_mechanism(
      "ERCC2 downregulation in metastatic colon cancer","ERCC2"))
    self.assertFalse(expert_triage.is_gene_suppression_mechanism(
      "NOP2 increases ERCC2 mRNA stability and suppresses malignant progression","ERCC2"))
    self.assertFalse(expert_triage.is_gene_suppression_mechanism(
      "BRD4 inhibition restores ERCC2 protein and suppresses tumor growth","ERCC2"))
 def test_expert_reports_export_tiers_and_falsification_not_doses(self):
    p=pair('EXAMPLE',sig=[{'provider':'perturbational_transcriptomics','depth':1}])
    with TemporaryDirectory() as td:
      summary,top=expert_triage.report(td,[p],demo_dossier())
      self.assertEqual(summary['expert_bucket_counts']['T3'],1)
      content=(Path(td)/'EXPERT_RESEARCH_REPORT.md').read_text()
      self.assertIn('NOT a validated efficacy estimate',content)
      self.assertIn('Critical limitations',content)
      self.assertTrue((Path(td)/'EXPERT_DRUG_CATALOGUE.jsonl').exists())
      self.assertNotIn('prescription',content.lower())


class ConsoleContract(unittest.TestCase):
 def test_windows_cp1252_cli_does_not_dump_full_patient_or_drug_result(self):
    import sys,io
    from unittest.mock import patch
    from v2810 import pipeline
    fake={'version':'V2.8.8','primary_genes':['DEMO1'],
      'master_entities':12,'research_evidence_tier_counts':{'D':2},
      'expert_research_bucket_counts':{'T3':2},
      'expert_genetics_dossier':'expert_genetics/GENETIC_EXPERT_SUMMARY.md',
      'expert_research_report':'expert_research/EXPERT_RESEARCH_REPORT.md',
      'research_shortlist':[{'candidate':'TEST_UNICODE_Δ»'}]*10000,
      'patient_specific_variants':'MUST_NOT_PRINT'}
    sink=io.BytesIO()
    wrapper=io.TextIOWrapper(sink,encoding='cp1252',write_through=True)
    with patch.object(sys,'argv',['pipeline.py','--case','synthetic_case.json','--no-ai']),\
         patch.object(pipeline,'run',return_value=fake),\
         patch.object(sys,'stdout',wrapper):
        pipeline.main()
    wrapper.flush()
    printed=sink.getvalue().decode('cp1252')
    self.assertIn('"completed": true',printed.lower())
    self.assertIn('EXPERT_RESEARCH_REPORT.md',printed)
    self.assertNotIn('MUST_NOT_PRINT',printed)
    self.assertNotIn('TEST_UNICODE',printed)

if __name__=='__main__':unittest.main()
