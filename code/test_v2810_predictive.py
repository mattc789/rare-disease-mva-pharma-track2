"""Synthetic, no-patient, no-network scientific tests for V2.8.6 predictor.

Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch
"""
from __future__ import annotations
import importlib.util
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

try:
    from v2810 import predictive
except ImportError:
    # Local bundle audit, before installation into a specific development workstation.
    path=Path(__file__).with_name('predictive.py')
    if not path.is_file():
        raise
    spec=importlib.util.spec_from_file_location('predictive_bundle_test',path)
    predictive=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(predictive)


def model(effects=('missense',),phase='unresolved', gene='DEMO1'):
    return {'primary_gene_hypotheses':[
        {'gene':gene,'phase':{'state':phase},'actual_disease_mechanism_confirmed':False,
         'functional_objectives':[],
         'variant_effect_models':[{'variant_ref':'SHA256:synthetic'+str(i),
                                   'effect_class':effect,'hgvsp_local':'PRIVATE_SHOULD_NOT_EXPORT'}
                                  for i,effect in enumerate(effects)]}],
        'phenotype_counts':{'present':1,'absent':1,'unknown':0}}


def route(family='direct_gene',action='increase',gene='DEMO1',variable='activity',
          target=None, state='conditional'):
    return {'route_id':'R1','family':family,'desired_action':action,
            'target_concept':target or gene,'mechanism_concept':'restore '+gene+' function',
            'variable':variable,'state':state,
            'targets':[{'gene': target or gene}],
            'scientific_gate':{'hard_hold':state=='hold','issues':[]}}


def pair(family='direct_gene',action='increase',variable='activity',target=None,
         tier='D',signals=None,name='APPROVED DEMODRUG',group='ChEMBL_PARENT:CHEMBL123'):
    return {'ingredient':name,'member_labels':[name],
            'scientific_group_id':group,'route_id':'R1',
            'route':route(family,action,variable=variable,target=target),
            'tier':tier,'score':1.0,'signals':signals or [
                {'route_id':'R1','provider':'ChEMBL_MoA_SIGNOR','depth':2,'mechanism_id':901}]}


def approved():
    return {'overall_state':'single_ingredient_product_listed','single_product_listed':True,
            'registry_names':['FDA Drugs@FDA'], 'jurisdictions':['United States'],
            'records':[{'register':'FDA Drugs@FDA','jurisdiction':'United States',
                        'status':'approved_on_snapshot','product_context':'single_ingredient_product',
                        'specific_ingredient_form':'APPROVED DEMODRUG'}]}


def combos():
    return {'overall_state':'combination_only','single_product_listed':False,
            'records':[{'product_context':'combination_component','status':'approved_on_snapshot'}]}


class V286PredictionTest(unittest.TestCase):
    def test_registry_requires_unwithdrawn_single_ingredient(self):
        self.assertTrue(predictive.registry_eligibility(approved())['eligible_for_research_search'])
        self.assertFalse(predictive.registry_eligibility(combos())['eligible_for_research_search'])
        self.assertFalse(predictive.registry_eligibility({})['eligible_for_research_search'])
        old=approved()
        old['records'][0]['status']='withdrawn'
        self.assertFalse(predictive.registry_eligibility(old)['eligible_for_research_search'])

    def test_patient_variant_changes_rank_for_same_drug_and_same_network(self):
        drug=pair(variable='activity')
        missense=predictive.predict_pair(drug,model(('missense',)),approved())
        stop=predictive.predict_pair(drug,model(('premature_stop',)),approved())
        self.assertGreater(missense['mechanistic_priority_score'],stop['mechanistic_priority_score'])
        self.assertIn('missense',missense['affected_variant_effect_classes'])
        self.assertIn('premature_stop',stop['affected_variant_effect_classes'])

    def test_readthrough_applies_only_to_nonsense_not_frameshift(self):
        drug=pair(family='variant_rescue',variable='full_length_translation',
                  target='DEMO1 premature-stop readthrough')
        nons=predictive.predict_pair(drug,model(('premature_stop',)),approved())
        frame=predictive.predict_pair(drug,model(('frameshift',)),approved())
        self.assertGreater(nons['mechanistic_priority_score'],frame['mechanistic_priority_score'])
        self.assertEqual(frame['prediction_class'],'CONFLICTING_DIRECTION')
        self.assertEqual(frame['mechanistic_priority_score'],0)

    def test_decreasing_defective_gene_fails_directionality(self):
        down=pair(action='decrease')
        x=predictive.predict_pair(down,model(('missense',)),approved())
        self.assertEqual(x['prediction_class'],'CONFLICTING_DIRECTION')
        self.assertEqual(x['mechanistic_priority_score'],0)

    def test_uncertain_phase_does_not_claim_confirmed_recessive_diagnosis(self):
        x=predictive.predict_pair(pair(),model(('missense','premature_stop'),phase='unresolved'),approved())
        self.assertEqual(x['phase_assessment'],'unresolved')
        self.assertIn('allelic_phase_not_confirmed_as_recessive_causal_evidence',
                      x['potential_risks_and_gaps'])
        self.assertFalse(x['clinical_suitability_verified'])

    def test_no_approved_drugs_in_prediction_or_AI_priority_queue(self):
        candidates=[pair(name='APPROVED DEMODRUG'),
                    pair(name='UNAPPROVED DEMODRUG',group='GOTHER')]
        catalog={'APPROVED DEMODRUG':approved(),'UNAPPROVED DEMODRUG':{}}
        ranked,summary=predictive.rank_pairs(candidates,model(),catalog,
            resolver=lambda name,registry:registry.get(name,{}))
        self.assertEqual(len(ranked),1)
        self.assertEqual(ranked[0]['ingredient'],'APPROVED DEMODRUG')
        self.assertEqual(summary['excluded_unregistered_or_combo_pairs'],1)

    def test_no_patient_exposure_yields_no_simulated_target_engagement(self):
        x=predictive.predict_pair(pair(signals=[{'route_id':'R1','provider':'ChEMBL_measured_activity',
                            'assay_id':'A1','value_nM':45,'direct_target_assay_verified':True}]),
                            model(),approved())
        self.assertEqual(x['exposure_target_assay_context']['state'], 'HUMAN_UNBOUND_EXPOSURE_UNKNOWN')
        self.assertIsNone(x['predicted_probability_of_clinical_benefit'])
        self.assertIsNone(x['predicted_degree_of_patient_rescue'])

    def test_curated_exposure_only_illustrates_assay_vs_plasma_context(self):
        drug=pair(signals=[{'route_id':'R1','provider':'ChEMBL_measured_activity',
                            'assay_id':'A1','value_nM':50,'direct_target_assay_verified':True}])
        with TemporaryDirectory() as td:
            path=Path(td)/'pk.json'
            path.write_text(json.dumps({'schema':'V286_CURATED_PK_1','entries':[
                {'compound':'APPROVED DEMODRUG','unbound_cmax_nM':100,
                 'data_type':'measured_unbound_plasma_cmax','units':'nM',
                 'source_reference':'synthetic-paper','population':'synthetic-adult',
                 'source_checked':True}]}),encoding='utf-8')
            index=predictive._load_exposure(path)
            x=predictive.predict_pair(drug,model(),approved(),pk_index=index)
        info=x['exposure_target_assay_context']
        self.assertAlmostEqual(info['functional_target_engagement'],0.6667,places=3)
        self.assertEqual(info['clinical_or_paediatric_exposure'],'not_estimated')
        self.assertFalse(x['clinical_suitability_verified'])
        self.assertIsNone(x['predicted_probability_of_clinical_benefit'])

    def test_unverified_or_mislabeled_pk_is_discarded(self):
        with TemporaryDirectory() as td:
            path=Path(td)/'pk.json'
            path.write_text(json.dumps({'schema':'V286_CURATED_PK_1','entries':[
                {'compound':'APPROVED DEMODRUG','unbound_cmax_nM':100,
                 'data_type':'measured_unbound_plasma_cmax','units':'nM',
                 'source_reference':'synthetic','population':'adult',
                 'source_checked':False}]}))
            self.assertEqual(predictive._load_exposure(path),{})

    def test_evidence_A_does_not_become_a_clinical_efficacy_percentage(self):
        p=pair(tier='A')
        x=predictive.predict_pair(p,model(),approved())
        self.assertEqual(x['research_evidence_tier'],'A')
        self.assertIsNone(x['predicted_probability_of_clinical_benefit'])
        self.assertFalse(x['score_is_calibrated_efficiency_probability'])

    def test_held_route_not_promoted_even_if_direct_assay_exists(self):
        p=pair(signals=[{'route_id':'R1','provider':'ChEMBL_measured_activity',
                            'value_nM':5,'direct_target_assay_verified':True}])
        p['route']['state']='hold'
        p['route']['scientific_gate']['hard_hold']=True
        x=predictive.predict_pair(p,model(),approved())
        self.assertEqual(x['mechanistic_priority_score'],0)
        self.assertEqual(x['prediction_class'],'CONFLICTING_DIRECTION')

    def test_multi_gene_case_has_no_guess_when_route_does_not_name_gene(self):
        p=pair(family='signed_upstream',target='CHEK1')
        p['route']['mechanism_concept']='increase upstream regulator activity'
        g=model(gene='DEMO1')
        g['primary_gene_hypotheses']+=model(gene='DEMO2')['primary_gene_hypotheses']
        x=predictive.predict_pair(p,g,approved())
        self.assertIsNone(x['primary_gene_hypothesis'])
        self.assertEqual(x['prediction_class'],'CASE_LINK_UNRESOLVED')

    def test_network_indirect_cannot_be_predicted_strongly(self):
        p=pair(signals=[{'route_id':'R1','provider':'ChEMBL_MoA_SIGNOR','depth':3}])
        x=predictive.predict_pair(p,model(),approved())
        self.assertLessEqual(x['mechanistic_priority_score'],45)
        self.assertEqual(x['prediction_class'],'HYPOTHESIS_ONLY')

    def test_variants_not_exported_as_sequences(self):
        x=predictive.predict_pair(pair(),model(),approved())
        text=json.dumps(x)
        self.assertNotIn('PRIVATE_SHOULD_NOT_EXPORT',text)
        self.assertTrue(x['patient_sequence_included'] is False)
        self.assertTrue(all(v.startswith('SHA256:') for v in x['affected_variant_refs']))

    def test_model_predicts_ranking_not_calibrated_clinical_effect(self):
        eligible,summary=predictive.rank_pairs([pair()],model(),{'APPROVED DEMODRUG':approved()},
                                               resolver=lambda name,catalog:catalog.get(name,{}))
        self.assertEqual(len(eligible),1)
        self.assertFalse(summary['clinical_prediction_validated'])
        self.assertEqual(eligible[0]['tier'],'D')
        self.assertGreater(eligible[0]['patient_first_score'],0)
        self.assertEqual(summary['top_candidate_route_pairs'][0]['research_evidence_tier'],'D')

    def test_single_drug_multiroute_conflicts_flagged_in_reports(self):
        a=pair(action='increase')
        b=pair(action='decrease')
        b['route_id']='R2'
        b['route']['route_id']='R2'
        with TemporaryDirectory() as td:
            path=Path(td)/'out'
            summ,top=predictive.build_reports(path,[],[a,b],model(),
                {'APPROVED DEMODRUG':approved()},resolver=lambda name,catalog:catalog.get(name,{}))
            self.assertEqual(summ['approved_single_ingredient_drug_groups_considered'],1)
            self.assertEqual(top[0]['cross_route_direction_conflicts'],['DEMO1'])
            self.assertTrue((path/'PATIENT_FIRST_REPORT.md').exists())
            self.assertTrue((path/'PATIENT_FIRST_RANKING.csv').exists())
            self.assertTrue((path/'PATIENT_FIRST_CANDIDATES.jsonl').exists())

    def test_score_uses_phase_as_uncertainty_not_mendelian_certainty(self):
        p=pair()
        phase_unknown=predictive.predict_pair(p,model(phase='unresolved'),approved())
        phase_trans=predictive.predict_pair(p,model(phase='trans_supported_by_vcf'),approved())
        self.assertGreaterEqual(phase_trans['mechanistic_priority_score'],phase_unknown['mechanistic_priority_score'])
        self.assertFalse(phase_trans['clinical_suitability_verified'])


if __name__=='__main__':
    unittest.main()
