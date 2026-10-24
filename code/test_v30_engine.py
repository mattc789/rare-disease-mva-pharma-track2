import json,unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from v30 import license_gate,patient_state,mechanism_engine,rescue_objectives
from v30 import candidate_graph,falsification,dossier,openai_reasoner,literature_attack,api_budget

def private_case():
    return {
      "case_id":"SYNTHETIC","genes":{"GENE1":[{},{}]},
      "phenotype_features":[],"diagnostic_certainty":"unknown",
      "residual_protein_status":"unknown","compound_heterozygosity_confirmed":False,
      "expression_measured":False,
    }

def v28_mechanisms():
    return {"primary_gene_hypotheses":[{
      "gene":"GENE1","phase":{"state":"unresolved","trans_alleles_established":False},
      "variant_effect_models":[
        {"variant_ref":"SHA256:abc","effect_class":"premature_stop",
         "clinical_assertion":"reported_pathogenic","clinvar_review":"multiple submitters","clinvar_heuristic_weight":6},
        {"variant_ref":"SHA256:def","effect_class":"missense",
         "clinical_assertion":"not_available","clinvar_review":"unknown","clinvar_heuristic_weight":0}],
      "actual_disease_mechanism_confirmed":False}]}

class V30StateMechanismTests(unittest.TestCase):
    def test_state_to_mechanisms_to_objectives(self):
        state=patient_state.build(private_case(),v28_mechanisms(),{}, {}, {})
        self.assertEqual(state["genes"][0]["gene"],"GENE1")
        self.assertFalse(state["clinical_diagnosis_asserted"])
        mech=mechanism_engine.build(state)
        kinds={x["mechanism_type"] for x in mech["hypotheses"]}
        self.assertIn("loss_of_function_or_truncated_product",kinds)
        self.assertIn("protein_stability_or_folding_defect",kinds)
        gof=[x for x in mech["hypotheses"] if x["mechanism_type"]=="gain_of_function"][0]
        self.assertEqual(gof["status"],"NOT_INFERRED_WITHOUT_DIRECT_EVIDENCE")
        obj=rescue_objectives.build(mech)
        types={x["objective_type"] for x in obj["objectives"]}
        self.assertIn("restore_functional_gene_product",types)
        self.assertIn("restore_full_length_translation",types)
        self.assertNotIn("decrease_abnormal_activity",types)

class V30GraphFalsificationTests(unittest.TestCase):
    def _objects(self):
        state=patient_state.build(private_case(),v28_mechanisms())
        mech=mechanism_engine.build(state);obj=rescue_objectives.build(mech)
        route={"route_id":"R1","case_gene":"GENE1","family":"direct_gene",
          "target_concept":"GENE1","desired_action":"increase","variable":"activity",
          "scientific_gate":{"hard_hold":False,"issues":[]}}
        pair={"ingredient":"DEMODRUG","scientific_group_id":"G1","route_id":"R1","route":route,
          "tier":"B","signals":[{"provider":"ChEMBL_measured_activity","direct_target_assay_verified":True}],
          "functional_evidence":{"same_drug_functional_rescue":True},
          "expert_triage":{"investigation_bucket":"T1","label_state":"PASS"},
          "evidence_blockers":[]}
        pred={"candidate":"DEMODRUG","route_id":"R1","prediction_class":"DIRECTIONALLY_CONSISTENT",
          "mechanistic_priority_score":50,"public_label_safety_state":"PASS",
          "registry":{"eligible_for_research_search":True}}
        return state,mech,obj,route,pair,pred
    def test_conditional_lead_requires_real_support(self):
        state,mech,obj,route,pair,pred=self._objects()
        graph=candidate_graph.build(obj,[route],[pair],[pred])
        f=falsification.evaluate(graph)
        self.assertEqual(f["candidates"][0]["decision"],"CONDITIONAL_RESEARCH_LEAD")
        self.assertFalse(f["candidates"][0]["research_priority_score_is_efficacy_probability"])
    def test_direction_conflict_is_killed(self):
        state,mech,obj,route,pair,pred=self._objects()
        pred["prediction_class"]="CONFLICTING_DIRECTION"
        graph=candidate_graph.build(obj,[route],[pair],[pred])
        f=falsification.evaluate(graph)
        self.assertEqual(f["candidates"][0]["decision"],"KILLED")

class V30LiteratureAttackTests(unittest.TestCase):
    def test_independent_direct_and_rescue_can_promote_research_only(self):
        base={"schema":"X","candidates":[{
          "candidate":"DEMODRUG","gene":"GENE1","route_id":"R1","objective_id":"O1",
          "decision":"HOLD","research_priority_score":20,
          "supporting_factors":[],"hold_reasons":[
            "published_or_assay_evidence_not_yet_A_or_B",
            "direct_directional_target_assay_missing",
            "same_drug_disease_relevant_functional_rescue_missing"],
          "kill_reasons":[],"falsification_tests":[]}]}
        attack={"results":[{
          "candidate":"DEMODRUG","evidence_level":"P2_INDEPENDENT_DIRECT_AND_RESCUE",
          "direct_sources":["PMC1"],"rescue_sources":["PMC2"],
          "counterevidence_sources":[]}]}
        out=literature_attack.apply(base,attack)
        self.assertEqual(out["candidates"][0]["decision"],"CONDITIONAL_RESEARCH_LEAD")
        self.assertFalse(out["candidates"][0].get("clinical_recommendation",False))
    def test_single_source_does_not_promote(self):
        base={"candidates":[{
          "candidate":"DEMODRUG","decision":"HOLD","research_priority_score":20,
          "supporting_factors":[],"hold_reasons":[],"kill_reasons":[]}]}
        attack={"results":[{
          "candidate":"DEMODRUG","evidence_level":"P1_SINGLE_SOURCE_DIRECT_AND_RESCUE",
          "direct_sources":["PMC1"],"rescue_sources":["PMC1"],"counterevidence_sources":[]}]}
        out=literature_attack.apply(base,attack)
        self.assertEqual(out["candidates"][0]["decision"],"HOLD")

class V30PrivacyAndDossierTests(unittest.TestCase):
    def test_api_envelope_contains_no_allele_or_phenotype(self):
        state=patient_state.build(private_case(),v28_mechanisms())
        mech=mechanism_engine.build(state);obj=rescue_objectives.build(mech)
        fals={"candidates":[{"candidate":"DEMODRUG","gene":"GENE1","route_target":"TARGET1",
          "objective_type":"restore_functional_gene_product","decision":"HOLD",
          "supporting_factors":[],"hold_reasons":["direct_assay_missing"],"kill_reasons":[]}]}
        payload=openai_reasoner.safe_envelope(state,mech,obj,fals)
        raw=json.dumps(payload)
        self.assertNotIn("SHA256:",raw)
        self.assertNotIn("variant_ref",raw)
        self.assertNotIn("phenotype",raw.lower())
        self.assertNotIn("case_id",raw.lower())
    def test_privacy_guard_rejects_raw_variant_fields(self):
        with self.assertRaises(ValueError):
            openai_reasoner._assert_safe({"chrom":"1","pos":123,"ref":"A","alt":"G"})
    def test_no_lead_is_explicit(self):
        state=patient_state.build(private_case(),v28_mechanisms())
        mech=mechanism_engine.build(state);obj=rescue_objectives.build(mech)
        f={"candidates":[{"candidate":"X","decision":"HOLD","hold_reasons":["missing"],
                         "falsification_tests":[],"gene":"GENE1"}]}
        d=dossier.build(state,mech,obj,f)
        self.assertEqual(d["research_conclusion"],"NO_VALIDATED_LEAD")
        self.assertIsNone(d["proposed_research_lead"])

class V30BudgetTests(unittest.TestCase):
    def test_budget_guard_blocks_conservative_overrun(self):
        import os
        from unittest.mock import patch
        with TemporaryDirectory() as td, patch.dict(os.environ,{"RDMVA_OPENAI_BUDGET_USD":"0.01"}):
            with self.assertRaisesRegex(RuntimeError,"V30_OPENAI_BUDGET_HOLD"):
                api_budget.reserve(td,b"x"*100000,2200,2.0,10.0,"gpt-6.1-sol")
    def test_budget_actual_cost_is_recorded(self):
        import os
        from unittest.mock import patch
        with TemporaryDirectory() as td, patch.dict(os.environ,{"RDMVA_OPENAI_BUDGET_USD":"3.0"}):
            h,_,_=api_budget.reserve(td,b"{}",100,2.0,10.0,"gpt-6.1-sol")
            api_budget.finish(td,h,{"input_tokens":1000,"output_tokens":100},2.0,10.0)
            s=api_budget.status(td)
            self.assertEqual(s["requests"],1)
            self.assertGreater(s["actual_estimate_usd"],0)
            self.assertLessEqual(s["actual_estimate_usd"],s["reserved_usd"])

class V30LicenseTests(unittest.TestCase):
    def test_default_source_policy(self):
        result=license_gate.assert_default_profile()
        self.assertTrue(result["default_profile_valid"])
        self.assertFalse(result["red_sources_enabled"])

if __name__=="__main__":unittest.main()


class V30NativeRouteTests(unittest.TestCase):
    def test_native_routes_are_derived_before_drug_screen(self):
        from unittest.mock import patch
        from v30 import native_routes
        private={"genes":{"GENE1":[{}]}}
        objectives={"objectives":[
          {"objective_id":"OBJ1","gene":"GENE1","objective_type":"increase_residual_activity",
           "desired_action":"increase","variable":"protein_specific_activity",
           "priority_tier":"O1","drug_search_eligible":True,
           "requirements_before_promotion":[]}
        ]}
        edges=[{"source_gene":"REG1","source_acc":"PREG1","target_gene":"GENE1",
                "target_acc":"PGENE1","sign":1,"direct":True,
                "signor_id":"S1","pmid":"1"}]
        with patch("v30.native_routes.v28_routes.signor_edges",return_value=edges), \
             patch("v30.native_routes.v28_routes.accession_for_gene",return_value="PGENE1"):
            result=native_routes.build(private,objectives,"unused")
        self.assertTrue(result["derived_before_drug_screen"])
        self.assertFalse(result["route_generation_uses_drug_names"])
        self.assertTrue(any(x["family"]=="v30_direct_gene" for x in result["routes"]))
        self.assertTrue(any(x["family"]=="v30_signed_upstream" for x in result["routes"]))

class V30CombinationTests(unittest.TestCase):
    def _graph(self):
        def p(drug,oid,otype,gene,tier="B",score=8):
            return {"candidate":drug,"objective_id":oid,"objective_type":otype,
                    "objective_priority":"O1","gene":gene,"route_id":drug+"_"+oid,
                    "route_target":gene,"route_desired_action":"increase",
                    "evidence_tier":tier,"objective_match_score":score,
                    "direct_target_assay_count":1,"same_drug_functional_rescue":False,
                    "functional_rescue_validated":False,"scientific_hard_hold":False,
                    "registry_eligible":True,"market_approved_single_agent":True,
                    "public_label_safety_state":"PASS","source_channels":["test"],
                    "evidence_blockers":[],"prediction_class":"DIRECTIONALLY_CONSISTENT",
                    "shortest_network_hops":0}
        return {"candidate_paths":[
          p("DRUG_A","O1","restore_functional_gene_product","GENE1"),
          p("DRUG_B","O2","restore_productive_splicing","GENE2"),
          p("DRUG_C","O3","increase_residual_activity","GENE3","C",6),
        ]}
    def _objectives(self):
        return {"objectives":[
          {"objective_id":"O1","gene":"GENE1","objective_type":"restore_functional_gene_product","priority_tier":"O1"},
          {"objective_id":"O2","gene":"GENE2","objective_type":"restore_productive_splicing","priority_tier":"O1"},
          {"objective_id":"O3","gene":"GENE3","objective_type":"increase_residual_activity","priority_tier":"O2"},
        ]}
    def test_combinations_require_unique_complementary_coverage(self):
        from v30 import combination_engine
        result=combination_engine.discover(self._graph(),self._objectives(),max_components=3,shortlist=10,max_results=20)
        pairs=[x for x in result["combinations"] if set(x["components"])=={"DRUG_A","DRUG_B"}]
        self.assertTrue(pairs)
        self.assertFalse(pairs[0]["synergy_claimed"])
        self.assertEqual(pairs[0]["component_objective_contributions"]["DRUG_A"],1)
        self.assertEqual(pairs[0]["component_objective_contributions"]["DRUG_B"],1)

    def test_challenge_hypothesis_does_not_require_proven_rescue(self):
        from v30 import combination_engine,challenge_hypotheses
        graph=self._graph();objs=self._objectives()
        combos=combination_engine.discover(graph,objs,max_components=3,shortlist=10,max_results=20)
        result=challenge_hypotheses.build(graph,objs,combos)
        singles={x["candidate"]:x for x in result["submission_level_single_agents"]}
        self.assertIn("DRUG_A",singles)
        self.assertFalse(singles["DRUG_A"]["same_drug_functional_rescue"])
        self.assertTrue(singles["DRUG_A"]["hypothesis_for_follow_up_not_evidence_it_works"])
        self.assertTrue(result["validated_disease_rescue_required_for_hypothesis"] is False)


class V30Track2EvidenceTests(unittest.TestCase):
    def test_source_bound_evidence_reorders_without_becoming_clinical_claim(self):
        from v30 import challenge_hypotheses
        track={
          "single_agent_hypotheses":[
            {"candidate":"A","challenge_priority_score":10.0,"clinical_recommendation":False},
            {"candidate":"B","challenge_priority_score":12.0,"clinical_recommendation":False}],
          "submission_level_single_agents":[
            {"candidate":"A","challenge_priority_score":10.0,"clinical_recommendation":False},
            {"candidate":"B","challenge_priority_score":12.0,"clinical_recommendation":False}],
          "combination_hypotheses":[],
          "has_submission_level_hypothesis":True}
        attack={"results":[
          {"candidate":"A","evidence_level":"P2_INDEPENDENT_DIRECT_AND_RESCUE",
           "overall_assessment":"DIRECT_RESCUE_EVIDENCE_FOUND",
           "direct_sources":["S1"],"rescue_sources":["S2"],
           "independent_support_sources":["S1","S2"],"counterevidence_sources":[],
           "quotes_verified_exactly":True},
          {"candidate":"B","evidence_level":"P0_NO_SOURCE_BOUND_SUPPORT",
           "overall_assessment":"INSUFFICIENT","direct_sources":[],"rescue_sources":[],
           "independent_support_sources":[],"counterevidence_sources":[],
           "quotes_verified_exactly":True}]}
        out=challenge_hypotheses.apply_source_bound_evidence(track,attack)
        self.assertEqual(out["submission_level_single_agents"][0]["candidate"],"A")
        self.assertFalse(out["submission_level_single_agents"][0]["clinical_recommendation"])
        self.assertEqual(out["submission_level_single_agents"][0]["source_bound_evidence_state"],
                         "P2_INDEPENDENT_DIRECT_AND_RESCUE")

    def test_combination_safety_never_treats_no_pair_mention_as_clearance(self):
        from v30 import safety_review
        track={
          "single_agent_hypotheses":[{"candidate":"A"},{"candidate":"B"}],
          "submission_level_single_agents":[{"candidate":"A"},{"candidate":"B"}],
          "combination_hypotheses":[{"components":["A","B"]}]}
        safety={"candidates":{
          "A":{"state":"PASS","flags":[]},
          "B":{"state":"PASS","flags":[]}}}
        out=safety_review.annotate_hypotheses(track,safety,{"A":{},"B":{}})
        c=out["combination_hypotheses"][0]["combination_safety_review"]
        self.assertEqual(c["state"],"NO_EXPLICIT_PAIR_MENTION_NOT_SAFETY_CLEARANCE")
        self.assertFalse(c["clinical_clearance"])
        self.assertTrue(c["no_interaction_found_does_not_mean_safe"])


class V30EvidenceSpecificityTests(unittest.TestCase):
    def test_generic_route_rescue_without_primary_gene_is_rejected(self):
        from v30 import literature_attack
        path={"gene":"LZTR1","route_target":"KRAS","member_labels":["SOTORASIB"]}
        passages=[{"source_id":"S1","passage_id":"P1",
                   "text":"Sotorasib inhibited KRAS and rescued CCDC6 protein levels."}]
        model={"value":{"evidence_items":[{
          "source_id":"S1","passage_id":"P1",
          "exact_quote":"Sotorasib inhibited KRAS and rescued CCDC6 protein levels.",
          "evidence_type":"SAME_DRUG_FUNCTIONAL_RESCUE","direction":"restore",
          "mechanistically_relevant":True,"reason":"generic downstream rescue"}]}}
        accepted,rejected=literature_attack._validate_items("SOTORASIB",path,passages,model)
        self.assertEqual(accepted,[])
        self.assertEqual(rejected[0]["rejection_reason"],
                         "primary_gene_absent_from_functional_rescue_quote")

    def test_class_counterevidence_can_be_accepted_without_candidate_name(self):
        from v30 import literature_attack
        path={"gene":"LZTR1","route_target":"KRAS","member_labels":["SOTORASIB"]}
        quote="Loss of LZTR1 promotes resistance to G12C inhibitors through reactivation of RAS signaling."
        passages=[{"source_id":"S2","passage_id":"P2","text":quote}]
        model={"value":{"evidence_items":[{
          "source_id":"S2","passage_id":"P2","exact_quote":quote,
          "evidence_type":"MECHANISTIC_COUNTEREVIDENCE","direction":"increase",
          "mechanistically_relevant":True,"reason":"class-level resistance"}]}}
        accepted,rejected=literature_attack._validate_items("SOTORASIB",path,passages,model)
        self.assertEqual(len(accepted),1)
        self.assertEqual(rejected,[])


class V30DeterministicCounterScanTests(unittest.TestCase):
    def test_primary_gene_route_resistance_is_caught_if_model_omits_it(self):
        from v30 import literature_attack
        path={"gene":"LZTR1","route_target":"KRAS"}
        text=("Loss of LZTR1 would promote resistance to G12C inhibitors, "
              "through feedback reactivation of KRAS and RAS signaling.")
        items=literature_attack._deterministic_counter_items(
          [{"source_id":"S","passage_id":"P","text":text}],path)
        self.assertEqual(len(items),1)
        self.assertEqual(items[0]["evidence_type"],"MECHANISTIC_COUNTEREVIDENCE")
        self.assertTrue(items[0]["deterministic_counter_scan"])


class V31DiseaseStateGraphTests(unittest.TestCase):
    def test_secondary_genes_expand_network_without_directional_claim(self):
        from v30 import disease_state_graph
        private={"genes":{"GENE1":[{}]},
                 "genome_wide_candidate_genes":[
                   {"gene":"GENE1","score":30,"phenotype_score":5,"damaging_pass_count":2},
                   {"gene":"GENE2","score":27,"phenotype_score":4,"damaging_pass_count":2,
                    "mapping_artifact_review":False},
                   {"gene":"WEAK","score":5,"phenotype_score":0.2,"damaging_pass_count":1}]}
        obj={"objectives":[{"objective_id":"O1","gene":"GENE1","objective_type":"x",
                            "priority_tier":"O1","desired_action":"increase","variable":"activity"}]}
        g=disease_state_graph.build(private,{"genes":[]},{"hypotheses":[]},obj,signor_edges=[])
        rows={x["gene"]:x for x in g["nodes"] if x.get("type")=="gene"}
        self.assertIn("GENE2",rows)
        self.assertFalse(rows["GENE2"]["directional_pharmacology_allowed"])
        self.assertFalse(g["secondary_genes_can_create_directional_drug_claims"])

class V31NetworkExpansionTests(unittest.TestCase):
    def test_unsigned_network_never_claims_direction(self):
        from unittest.mock import patch
        from v30 import network_expansion
        dsg={"network_seed_genes":["GENE1"]}
        fake={"genes":{"GENE1":{"open_targets":{
          "interactions":[{"partner_gene":"MOD1","score":0.9,"source_database":"string",
                           "pmids":["1"],"evidences":[{"evidenceScore":0.8,
                           "interactionDetectionMethodShortName":"experimental"}]}],
          "reactome_pathways_from_opentargets":[]}}}}
        with TemporaryDirectory() as td, patch(
            "v30.network_expansion.public_enrichment.enrich_genes",return_value=fake):
            x=network_expansion.build(dsg,{},td,online=False)
        edges=[e for e in x["edges"] if e["type"]=="open_targets_interaction_unsigned"]
        self.assertEqual(len(edges),1)
        self.assertFalse(edges[0]["direction_known"])
        self.assertFalse(x["direction_inferred_from_unsigned_network"])

class V31VariantContextTests(unittest.TestCase):
    def test_variant_predictions_do_not_upgrade_evidence(self):
        from v30 import variant_structure
        obj={"objectives":[{"objective_id":"O1","gene":"GENE1",
                            "objective_type":"stabilize_residual_protein",
                            "priority_tier":"O1"}]}
        structure={"genes":{"GENE1":{"variants":[{"protein_position":10,
                    "alphafold_plddt":95.0,"alphafold_confidence_class":"VERY_HIGH"}]}}}
        am={"variants":[{"gene":"GENE1","prediction":{"protein_variant":"A10V",
              "am_pathogenicity":0.95,"am_class":"likely_pathogenic"}}]}
        out=variant_structure.annotate_objectives(obj,structure,am)
        ctx=out["objectives"][0]["variant_specific_computational_context"]
        self.assertFalse(ctx["context_is_functional_proof"])
        self.assertFalse(ctx["context_changes_evidence_tier"])

class V31DGIdbTests(unittest.TestCase):
    def test_dgidb_pair_stays_discovery_only(self):
        from v30 import dgidb_discovery
        routes=[{"route_id":"R1","case_gene":"GENE1","target_concept":"TARGET1",
                 "targets":[{"gene":"TARGET1","accession":"P1"}],
                 "desired_action":"decrease","variable":"activity",
                 "scientific_gate":{}}]
        entities=[{"canonical_label":"DRUGA","entity_id":"E1",
                   "scientific_group_id":"G1","member_labels":["DRUGA"],
                   "global_approved_as_single":True,"combo_only":False,
                   "modality":"small_or_other_drug_component"}]
        hints={"hints":[{"candidate":"DRUGA","route_id":"R1","direction_match":True,
                         "interaction_score":1.0,"pmids":["1"],"upstream_sources":["SRC"]}]}
        pairs,a=dgidb_discovery.merge_pairs([],hints,routes,entities)
        self.assertEqual(a["added_pairs"],1)
        self.assertEqual(pairs[0]["tier"],"D")
        self.assertIn("dgidb_aggregator_requires_primary_source_verification",
                      pairs[0]["evidence_blockers"])
        self.assertFalse(pairs[0]["promotion_eligible_validated_lead"])

class V31MultiAxisTests(unittest.TestCase):
    def test_single_agent_records_multiple_rescue_axes(self):
        from v30 import challenge_hypotheses
        base={"candidate":"DRUGA","scientific_hard_hold":False,
              "market_approved_single_agent":True,"evidence_tier":"B",
              "objective_match_score":8,"direct_target_assay_count":1,
              "shortest_network_hops":0,"route_family":"v30_direct_gene",
              "same_drug_functional_rescue":False,"functional_rescue_validated":False,
              "evidence_blockers":[],"prediction_class":"DIRECTIONALLY_CONSISTENT",
              "route_desired_action":"increase","source_channels":["test"],
              "polypharmacology":{"multi_target":True,"multi_rescue_axis":True}}
        p1={**base,"gene":"GENE1","objective_axis":"AX1","objective_type":"OBJ1",
            "objective_id":"O1","objective_priority":"O1","route_id":"R1","route_target":"GENE1"}
        p2={**base,"gene":"GENE2","objective_axis":"AX2","objective_type":"OBJ2",
            "objective_id":"O2","objective_priority":"O2","route_id":"R2","route_target":"GENE2"}
        rows=challenge_hypotheses.classify_singles({"candidate_paths":[p1,p2]},max_results=10)
        self.assertEqual(rows[0]["candidate"],"DRUGA")
        self.assertEqual(rows[0]["covered_rescue_axis_count"],2)
        self.assertIn("multiple_rescue_axes_covered",rows[0]["why_plausible"])


class V31AlphaMissenseMappingTests(unittest.TestCase):
    def test_isoform_discordance_is_not_exposed_as_prediction(self):
        from unittest.mock import patch
        import gzip, json
        from v30 import variant_structure
        private={"genes":{"GENE1":[{
          "chrom":"1","pos":100,"ref":"A","alt":"G",
          "consequence":"missense_variant","hgvsp":"p.Asn10Lys"}]}}
        with TemporaryDirectory() as td:
            p=Path(td)/"am.tsv.gz"
            with gzip.open(p,"wt",encoding="utf-8") as f:
                f.write("#CHROM\tPOS\tREF\tALT\tgenome\tuniprot\ttranscript\tprotein_variant\tam_pathogenicity\tam_class\n")
                f.write("chr1\t100\tA\tG\thg38\tP1\tT1\tN9K\t0.9\tlikely_pathogenic\n")
            resource={"state":"READY","path":str(p)}
            out=variant_structure.lookup_alphamissense(private,resource,Path(td)/"out")
            row=out["variants"][0]
            self.assertEqual(row["state"],"MATCHED_GENOMIC_ISOFORM_DISCORDANT")
            self.assertFalse(row["protein_mapping_concordant"])
            self.assertIsNone(row["prediction"])
            self.assertEqual(out["matched"],0)
            self.assertEqual(out["mapping_discordant"],1)


class V31ObjectiveBindingRegressionTests(unittest.TestCase):
    def test_degradation_route_never_binds_nmd_or_generic_activity(self):
        from v30 import native_routes
        objectives={"objectives":[
          {"objective_id":"O_ACT","gene":"BUB1B","objective_type":"increase_residual_activity","drug_search_eligible":True},
          {"objective_id":"O_NMD","gene":"BUB1B","objective_type":"reduce_nmd_conditionally","drug_search_eligible":True},
          {"objective_id":"O_STAB","gene":"BUB1B","objective_type":"stabilize_residual_protein","drug_search_eligible":True},
          {"objective_id":"O_DEG","gene":"BUB1B","objective_type":"reduce_causal_degradation","drug_search_eligible":True},
          {"objective_id":"O_PROD","gene":"BUB1B","objective_type":"restore_functional_gene_product","drug_search_eligible":True},
        ]}
        route={"case_gene":"BUB1B","target_concept":"Pellino 1","variable":"activity",
          "mechanism_concept":"Pellino 1 drives ubiquitination and degradation of BubR1; reducing Pellino 1 preserves BubR1 stability."}
        types={x["objective_type"] for x in native_routes._bind_literature_objectives(route,objectives)}
        self.assertIn("stabilize_residual_protein",types)
        self.assertIn("reduce_causal_degradation",types)
        self.assertNotIn("reduce_nmd_conditionally",types)
        self.assertNotIn("increase_residual_activity",types)

    def test_unbound_native_route_cannot_attach_heuristically(self):
        from v30 import candidate_graph
        objectives={"objectives":[{"objective_id":"O_NMD","gene":"GENE1",
          "objective_type":"reduce_nmd_conditionally","objective_axis":"rna_surveillance",
          "priority_tier":"O3","desired_action":"decrease","variable":"nmd","drug_search_eligible":True}]}
        route={"route_id":"R1","case_gene":"GENE1","family":"v30_literature_discovered",
          "target_concept":"MOD1","desired_action":"decrease","variable":"activity",
          "objective_binding_state":"UNBOUND","rescue_objective_ids":[],"targets":[]}
        pair={"ingredient":"DRUGA","route_id":"R1","route":route,"tier":"C","signals":[],
              "combo_only":False,"global_approved_as_single":True}
        graph=candidate_graph.build(objectives,[route],[pair],[])
        self.assertEqual(graph["candidate_paths"],[])

class V31SignedNetworkRegressionTests(unittest.TestCase):
    def test_signed_network_objective_types_are_conservative(self):
        from unittest.mock import patch
        from v30 import signed_network_expansion
        objectives={"objectives":[
          {"objective_id":"A","gene":"GENE1","objective_type":"increase_residual_activity"},
          {"objective_id":"S","gene":"GENE1","objective_type":"stabilize_residual_protein"},
          {"objective_id":"D","gene":"GENE1","objective_type":"compensate_downstream_deficit"},
        ]}
        fake_edges=[
          {"source_acc":"P_UP","target_acc":"P_G","source_gene":"UP","target_gene":"GENE1",
           "sign":1,"direct":True,"pmid":"1","signor_id":"S1"},
          {"source_acc":"P_G","target_acc":"P_DOWN","source_gene":"GENE1","target_gene":"DOWN",
           "sign":-1,"direct":True,"pmid":"2","signor_id":"S2"},
        ]
        with patch("v30.signed_network_expansion.v28_routes.signor_edges",return_value=fake_edges), \
             patch("v30.signed_network_expansion.v28_routes.accession_for_gene",return_value="P_G"):
            out=signed_network_expansion.build({"genes":{"GENE1":[{}]}},objectives,Path("."),max_hops=2)
        ups={r["rescue_objective_type"] for r in out["routes"] if r["family"]=="v31_signed_upstream_multihop"}
        downs={r["rescue_objective_type"] for r in out["routes"] if r["family"]=="v31_signed_downstream_multihop"}
        self.assertEqual(ups,{"increase_residual_activity"})
        self.assertEqual(downs,{"compensate_downstream_deficit"})
        self.assertNotIn("stabilize_residual_protein",ups)

class V31Track2DistanceTests(unittest.TestCase):
    def test_one_total_hop_direct_drug_route_can_be_plausible_but_two_hops_cannot(self):
        from v30 import challenge_hypotheses
        base={"candidate":"DRUGA","scientific_hard_hold":False,
          "market_approved_single_agent":True,"evidence_tier":"C",
          "objective_match_score":8,"direct_target_assay_count":0,
          "route_family":"v30_downstream_compensation",
          "same_drug_functional_rescue":False,"functional_rescue_validated":False,
          "evidence_blockers":[],"prediction_class":"DIRECTIONALLY_CONSISTENT",
          "route_desired_action":"decrease","source_channels":["ChEMBL_MoA_SIGNOR"],
          "polypharmacology":{},"gene":"GENE1","objective_axis":"downstream_compensation",
          "objective_type":"compensate_downstream_deficit","objective_id":"O1",
          "objective_priority":"O2","route_id":"R1","route_target":"TARGET1",
          "drug_to_route_hops":0}
        one={**base,"shortest_network_hops":1,"route_to_case_hops":1}
        rows=challenge_hypotheses.classify_singles({"candidate_paths":[one]},max_results=10)
        self.assertEqual(rows[0]["status"],"PLAUSIBLE_CHALLENGE_HYPOTHESIS")
        two={**base,"shortest_network_hops":2,"route_to_case_hops":2}
        rows=challenge_hypotheses.classify_singles({"candidate_paths":[two]},max_results=10)
        self.assertEqual(rows[0]["status"],"EXPLORATORY_ONLY")


class V31ModifierContextTests(unittest.TestCase):
    def test_opposite_primary_overexpression_context_is_rejected_for_restore_objective(self):
        from v30 import modifier_discovery
        quote=("Targeting CEP170 by siRNA in CAG and H929 BUB1B-OE cells significantly "
               "decreased chromosomal plate width and increased mitotic bipolar spindle length.")
        passages=[{"source_id":"S1","passage_id":"P1","text":quote}]
        rec={"value":{"modifiers":[{
          "modifier_gene":"CEP170","intervention_action":"decrease",
          "effect_on_primary_defect":"improves","mechanism_class":"GENETIC_SUPPRESSOR",
          "functional_rescue_observed":True,"experiment_context":"BUB1B overexpression",
          "source_id":"S1","passage_id":"P1","exact_quote":quote,"reason":"test"}]}}
        accepted,rejected=modifier_discovery._validate(
            "BUB1B",passages,rec,{"desired_actions":["increase","restore"]})
        self.assertEqual(accepted,[])
        self.assertEqual(rejected[0]["rejection_reason"],
                         "opposite_primary_gene_overexpression_context")


class V31ModifierTranslationalContextTests(unittest.TestCase):
    def test_nonhuman_rescue_remains_exploratory_not_validated(self):
        from unittest.mock import patch
        from v30 import modifier_discovery
        findings={"accepted_modifiers":[{
          "primary_gene":"LZTR1","modifier_gene":"NF1",
          "intervention_action":"increase","effect_on_primary_defect":"improves",
          "mechanism_class":"GENETIC_SUPPRESSOR","functional_rescue_observed":True,
          "translational_context":"NONHUMAN_ONLY",
          "source_id":"S1","passage_id":"P1","exact_quote":"Nf1 overexpression rescued nowl loss in Drosophila.",
          "experiment_context":"Drosophila nowl loss model"}]}
        objectives={"objectives":[{"objective_id":"O1","gene":"LZTR1",
          "objective_type":"compensate_downstream_deficit","priority_tier":"O2",
          "drug_search_eligible":True}]}
        with patch("v30.modifier_discovery.v28_routes.accession_for_gene",return_value="P1"):
            routes=modifier_discovery.to_routes(findings,"unused",objectives)
        self.assertEqual(len(routes),1)
        self.assertFalse(routes[0]["scientific_gate"]["functional_rescue_validated"])
        self.assertFalse(routes[0]["scientific_gate"]["causal_direction_validated_for_disease"])
        self.assertIn("nonhuman_model_requires_human_validation",
                      routes[0]["scientific_gate"]["issues"])


class V31ExactGeneSymbolModifierTests(unittest.TestCase):
    def test_pou6f2_does_not_match_pou6f2_as1(self):
        from v30 import modifier_discovery
        self.assertTrue(modifier_discovery._gene_exact("POU6F2 expression was reduced.", "POU6F2"))
        self.assertFalse(modifier_discovery._gene_exact("POU6F2-AS1 binds YBX1.", "POU6F2"))

    def test_modifier_requires_primary_gene_in_same_exact_quote(self):
        from v30 import modifier_discovery
        quote="YBX1 knockdown rescued the growth phenotype caused by POU6F2-AS1 overexpression."
        passages=[{"source_id":"S1","passage_id":"P1","text":quote}]
        rec={"value":{"modifiers":[{
          "modifier_gene":"YBX1","intervention_action":"decrease",
          "effect_on_primary_defect":"improves","mechanism_class":"GENETIC_SUPPRESSOR",
          "functional_rescue_observed":True,"experiment_context":"human cells",
          "source_id":"S1","passage_id":"P1","exact_quote":quote,"reason":"rescue"}]}}
        accepted,rejected=modifier_discovery._validate("POU6F2",passages,rec,{})
        self.assertEqual(accepted,[])
        self.assertEqual(rejected[0]["rejection_reason"],"primary_gene_absent_from_exact_quote_and_passage_linkage_insufficient")


class V31ModifierSourceRepairTests(unittest.TestCase):
    def test_exact_source_sentence_repairs_ai_omitted_figure_reference(self):
        from v30 import modifier_discovery
        text=("Moreover, the half-life of BUB1B was markedly extended in cells exhibiting "
              "elevated levels of KIFC1 expression (Figure 7C,D). Our evidence indicates "
              "that KIFC1 influences the protein expression level of BUB1B by stabilising "
              "the process of BUB1B degradation.")
        ai_quote=("Moreover, the half-life of BUB1B was markedly extended in cells exhibiting "
                  "elevated levels of KIFC1 expression. Our evidence indicates that KIFC1 "
                  "influences the protein expression level of BUB1B by stabilising the process "
                  "of BUB1B degradation.")
        rec={"value":{"modifiers":[{
          "modifier_gene":"KIFC1","intervention_action":"increase",
          "effect_on_primary_defect":"improves","mechanism_class":"DEGRADATION",
          "functional_rescue_observed":False,"experiment_context":"human pancreatic cells",
          "source_id":"S","passage_id":"P","exact_quote":ai_quote,"reason":"stability"}]}}
        acc,rej=modifier_discovery._validate(
            "BUB1B",[{"source_id":"S","passage_id":"P","text":text}],rec,
            {"desired_actions":["increase","restore"]})
        self.assertEqual(rej,[])
        self.assertEqual(len(acc),1)
        self.assertTrue(acc[0]["quote_repaired_from_source"])
        self.assertIn("(Figure 7C,D)",acc[0]["exact_quote"])

    def test_same_verified_passage_can_link_ortholog_rescue_but_stays_nonhuman(self):
        from v30 import modifier_discovery
        text=("Human LZTR1 regulates Ras signaling. Nf1 is another regulator of Ras signaling, "
              "and mutations in Nf1 and LZTR1 both give rise to neurofibromatosis. "
              "Furthermore, the reduction in total sleep caused by loss of nowl is partially "
              "rescued by Nf1 overexpression in Drosophila.")
        quote=("Furthermore, the reduction in total sleep caused by loss of nowl is partially "
               "rescued by Nf1 overexpression in Drosophila.")
        rec={"value":{"modifiers":[{
          "modifier_gene":"NF1","intervention_action":"increase",
          "effect_on_primary_defect":"improves","mechanism_class":"GENETIC_SUPPRESSOR",
          "functional_rescue_observed":True,"experiment_context":"Drosophila nowl loss model",
          "source_id":"S","passage_id":"P","exact_quote":quote,"reason":"genetic rescue"}]}}
        acc,rej=modifier_discovery._validate(
            "LZTR1",[{"source_id":"S","passage_id":"P","text":text}],rec,
            {"desired_actions":["increase","restore"]})
        self.assertEqual(rej,[])
        self.assertEqual(acc[0]["primary_gene_linkage"],"same_verified_passage")
        self.assertEqual(acc[0]["translational_context"],"NONHUMAN_ONLY")
        self.assertFalse(acc[0]["nonhuman_rescue_is_human_validation"])


class V31DirectDrugRescueTests(unittest.TestCase):
    def test_exact_name_boundaries_reject_substrings(self):
        from v30 import direct_drug_rescue
        self.assertFalse(direct_drug_rescue._term_in("N6-methyladenosine modification","adenosine"))
        self.assertFalse(direct_drug_rescue._gene_in("POU6F2-AS1 drives cancer","POU6F2"))
        self.assertTrue(direct_drug_rescue._gene_in("Lztr1R409C/+ mice","LZTR1"))

    def test_adjudicated_mutant_rescue_pair_is_research_only(self):
        from v30 import direct_drug_rescue
        payload={"adjudication_completed":True,"accepted_findings":[{
          "candidate":"DRUGX","entity_id":"E1","scientific_group_id":"G1",
          "member_labels":["DRUGX"],"gene":"GENE1","objective_id":"O1",
          "objective_type":"bypass_or_compensate_downstream","objective_priority":"O2",
          "source_id":"PMC1","passage_id":"B1","license":"CC-BY-4.0",
          "exact_quote":"DrugX treatment improved the phenotype in Gene1 mutant mice.",
          "supporting_passages":[{"passage_id":"B1","exact_quote":"DrugX treatment improved the phenotype in Gene1 mutant mice.","license":"CC-BY-4.0"}],
          "evidence_class":"MUTANT_MODEL_FUNCTIONAL_RESCUE",
          "source_bound_primary_gene_rescue":False,
          "source_bound_mutant_model_functional_rescue":True,
          "translational_context":"NONHUMAN_ONLY"}]}
        pairs,a=direct_drug_rescue.merge_pairs([],payload,{"objectives":[]})
        self.assertEqual(a["added_pairs"],1)
        self.assertEqual(pairs[0]["tier"],"C")
        self.assertTrue(pairs[0]["source_bound_mutant_model_functional_rescue"])
        self.assertFalse(pairs[0]["promotion_eligible_validated_lead"])
        self.assertFalse(pairs[0]["clinical_recommendation"])

    def test_source_bound_mutant_rescue_can_be_track2_plausible_not_validated(self):
        from v30 import candidate_graph,challenge_hypotheses,direct_drug_rescue
        objectives={"objectives":[{
          "objective_id":"O1","gene":"GENE1","objective_type":"bypass_or_compensate_downstream",
          "objective_axis":"downstream_compensation","priority_tier":"O2",
          "desired_action":"restore","variable":"disease_relevant_pathway_output",
          "drug_search_eligible":True}]}
        payload={"adjudication_completed":True,"accepted_findings":[{
          "candidate":"DRUGX","entity_id":"E1","scientific_group_id":"G1",
          "member_labels":["DRUGX"],"gene":"GENE1","objective_id":"O1",
          "objective_type":"bypass_or_compensate_downstream","objective_priority":"O2",
          "source_id":"PMC1","passage_id":"B1","license":"CC-BY-4.0",
          "exact_quote":"DrugX treatment improved the phenotype in Gene1 mutant mice.",
          "evidence_class":"MUTANT_MODEL_FUNCTIONAL_RESCUE",
          "source_bound_primary_gene_rescue":False,
          "source_bound_mutant_model_functional_rescue":True,
          "translational_context":"NONHUMAN_ONLY"}]}
        pairs,a=direct_drug_rescue.merge_pairs([],payload,objectives)
        graph=candidate_graph.build(objectives,a["routes"],pairs,[])
        t=challenge_hypotheses.build(graph,objectives,{"combinations":[]})
        rows={x["candidate"]:x for x in t["submission_level_single_agents"]}
        self.assertIn("DRUGX",rows)
        self.assertEqual(rows["DRUGX"]["status"],"PLAUSIBLE_CHALLENGE_HYPOTHESIS")
        self.assertTrue(rows["DRUGX"]["source_bound_mutant_model_functional_rescue"])
        self.assertEqual(rows["DRUGX"]["translational_context"],"NONHUMAN_ONLY")
        self.assertFalse(rows["DRUGX"]["clinical_recommendation"])


class V31ContextCounterEvidenceTests(unittest.TestCase):
    def test_generic_resistance_is_not_deterministic_counter_to_direct_drug_rescue(self):
        from v30 import literature_attack
        path={"gene":"LZTR1","route_target":"LZTR1 protein abundance",
              "route_family":"v31_direct_source_drug_rescue",
              "member_labels":["BORTEZOMIB"]}
        passages=[{"source_id":"S1","passage_id":"P1",
                   "text":"Loss of LZTR1 is associated with resistance to EGFR inhibitors."}]
        items=literature_attack._deterministic_counter_items(
            passages,path,candidate="BORTEZOMIB")
        self.assertEqual(items,[])

    def test_direct_rescue_with_context_counterevidence_remains_track2_with_penalty(self):
        from v30 import challenge_hypotheses
        track={
          "single_agent_hypotheses":[{
            "candidate":"DRUGX","challenge_priority_score":30.0,
            "clinical_recommendation":False}],
          "submission_level_single_agents":[{
            "candidate":"DRUGX","challenge_priority_score":30.0,
            "clinical_recommendation":False}],
          "combination_hypotheses":[],"has_submission_level_hypothesis":True}
        attack={"results":[{
          "candidate":"DRUGX","evidence_level":"P0_CONFLICTED_SOURCE_BOUND_EVIDENCE",
          "overall_assessment":"DIRECT_RESCUE_EVIDENCE_FOUND",
          "direct_sources":[],"rescue_sources":["S1"],
          "independent_support_sources":["S1"],"counterevidence_sources":["S2"],
          "missing_evidence":[],"quotes_verified_exactly":True}]}
        out=challenge_hypotheses.apply_source_bound_evidence(track,attack)
        self.assertEqual(len(out["submission_level_single_agents"]),1)
        row=out["submission_level_single_agents"][0]
        self.assertTrue(row["context_conflict_requires_resolution"])
        self.assertFalse(row["submission_hold_from_source_evidence"])
        self.assertEqual(row["challenge_priority_score_source_adjusted"],25.0)


class V31CrossSourceBridgeTests(unittest.TestCase):
    def test_cross_source_bridge_is_track2_plausible_not_validated(self):
        from v30 import cross_source_bridge,candidate_graph,challenge_hypotheses
        objectives={"objectives":[{
          "objective_id":"O1","gene":"GENE1",
          "objective_type":"bypass_or_compensate_downstream",
          "objective_axis":"downstream_compensation","priority_tier":"O2",
          "desired_action":"restore","variable":"disease_relevant_pathway_output",
          "drug_search_eligible":True}]}
        result={"bridges":[{
          "candidate":"DRUGX","entity_id":"E1","scientific_group_id":"G1",
          "member_labels":["DRUGX"],"primary_gene":"GENE1",
          "target_symbol":"TARGET1","target_aliases":["TARGET1"],
          "desired_action":"decrease","objective_id":"O1",
          "objective_type":"bypass_or_compensate_downstream",
          "objective_priority":"O2",
          "rescue_sources":[{"pmid":"1","title":"rescue","text_redistributed":False}],
          "drug_target_source":{"pmid":"2","title":"assay","text_redistributed":False},
          "bridge_class":"CROSS_SOURCE_MECHANISTIC_RESCUE",
          "same_drug_disease_rescue_demonstrated":False,
          "clinical_recommendation":False}]}
        pairs,a=cross_source_bridge.merge_pairs([],result)
        graph=candidate_graph.build(objectives,a["routes"],pairs,[])
        t=challenge_hypotheses.build(graph,objectives,{"combinations":[]})
        rows={x["candidate"]:x for x in t["submission_level_single_agents"]}
        self.assertIn("DRUGX",rows)
        self.assertEqual(rows["DRUGX"]["status"],"PLAUSIBLE_CHALLENGE_HYPOTHESIS")
        self.assertTrue(rows["DRUGX"]["source_bound_cross_source_bridge"])
        self.assertFalse(rows["DRUGX"]["same_drug_functional_rescue"])
        self.assertFalse(rows["DRUGX"]["clinical_recommendation"])


class V31CrossSourceBridgeTests(unittest.TestCase):
    def test_merge_requires_structured_bridge_and_stays_research_only(self):
        from v30 import cross_source_bridge,candidate_graph,challenge_hypotheses
        objectives={"objectives":[{
          "objective_id":"O1","gene":"GENE1",
          "objective_type":"bypass_or_compensate_downstream",
          "objective_axis":"downstream_compensation","priority_tier":"O2",
          "desired_action":"restore","variable":"pathway_output",
          "drug_search_eligible":True}]}
        result={"bridges":[{
          "candidate":"DRUGX","entity_id":"E1","scientific_group_id":"G1",
          "member_labels":["DRUGX"],"primary_gene":"GENE1",
          "target_symbol":"TARGET1","desired_action":"decrease",
          "objective_id":"O1","objective_type":"bypass_or_compensate_downstream",
          "objective_priority":"O2",
          "rescue_sources":[{"pmid":"1","title":"rescue source","text_redistributed":False}],
          "drug_target_source":{"pmid":"2","title":"drug-target source","text_redistributed":False},
          "bridge_class":"CROSS_SOURCE_MECHANISTIC_RESCUE",
          "same_drug_disease_rescue_demonstrated":False,
          "clinical_recommendation":False}]}
        pairs,a=cross_source_bridge.merge_pairs([],result)
        self.assertEqual(a["added_pairs"],1)
        self.assertEqual(pairs[0]["tier"],"C")
        self.assertFalse(pairs[0]["promotion_eligible_validated_lead"])
        graph=candidate_graph.build(objectives,a["routes"],pairs,[])
        t=challenge_hypotheses.build(graph,objectives,{"combinations":[]})
        rows={x["candidate"]:x for x in t["submission_level_single_agents"]}
        self.assertIn("DRUGX",rows)
        self.assertTrue(rows["DRUGX"]["source_bound_cross_source_bridge"])
        self.assertFalse(rows["DRUGX"]["clinical_recommendation"])

    def test_batch_summary_has_no_abstract_text_flag(self):
        from v30 import cross_source_bridge
        rows=cross_source_bridge._route_target_candidates([{
          "route_id":"R1","case_gene":"GENE1","family":"v31_signed_upstream_multihop",
          "state":"conditional","score":11.0,"network_hops":2,
          "rescue_priority_tier":"O1","desired_action":"decrease",
          "targets":[{"gene":"TARGET1","accession":"P1"}]}],max_targets=10)
        self.assertEqual(rows[0]["target_symbol"],"TARGET1")
        self.assertEqual(rows[0]["desired_action"],"decrease")


class V31ExposureFeasibilityTests(unittest.TestCase):
    def test_borderline_exposure_penalizes_but_does_not_kill(self):
        from v30 import exposure_feasibility
        track={
          "single_agent_hypotheses":[{"candidate":"DRUGX","challenge_priority_score":40.0}],
          "submission_level_single_agents":[{"candidate":"DRUGX","challenge_priority_score":40.0}],
          "combination_hypotheses":[],"has_submission_level_hypothesis":True}
        assessment={"candidate":"DRUGX","state":"EXPOSURE_BORDERLINE_PROTEIN_BINDING",
                    "clinical_recommendation":False,"dose_recommendation":False}
        out=exposure_feasibility.annotate_track2(track,[assessment])
        self.assertEqual(len(out["submission_level_single_agents"]),1)
        row=out["submission_level_single_agents"][0]
        self.assertEqual(row["challenge_priority_score_exposure_adjusted"],35.0)
        self.assertTrue(row["exposure_uncertainty_requires_resolution"])
        self.assertFalse(row["exposure_submission_hold"])

    def test_implausible_exposure_holds_submission(self):
        from v30 import exposure_feasibility
        track={
          "single_agent_hypotheses":[{"candidate":"DRUGX","challenge_priority_score":40.0}],
          "submission_level_single_agents":[{"candidate":"DRUGX","challenge_priority_score":40.0}],
          "combination_hypotheses":[],"has_submission_level_hypothesis":True}
        assessment={"candidate":"DRUGX","state":"EXPOSURE_IMPLAUSIBLE_FROM_LABEL_CMAX",
                    "clinical_recommendation":False,"dose_recommendation":False}
        out=exposure_feasibility.annotate_track2(track,[assessment])
        self.assertEqual(out["submission_level_single_agents"],[])
        self.assertFalse(out["has_submission_level_hypothesis"])


class V31RightsAuditTests(unittest.TestCase):
    def test_default_rights_profile_passes_and_conflicted_source_is_off(self):
        from v30 import rights_audit
        a=rights_audit.build()
        self.assertTrue(a["challenge_public_profile_pass"])
        self.assertEqual(a["failures"],[])
        self.assertEqual(
            a["special_findings"]["alphafold_db"]["status"],
            "ALLOWED_WITH_ATTRIBUTION")
        self.assertTrue(
            a["special_findings"]["alphafold_db"]["commercial_use_permitted"])
        self.assertEqual(
            a["special_findings"]["alphamissense"]["status"],
            "DISABLED_LICENSE_CONFLICT")
        self.assertFalse(
            a["special_findings"]["alphamissense"]["default_enabled"])

class V31ChallengeDossierTests(unittest.TestCase):
    def test_dossier_remains_research_only(self):
        from v30 import challenge_dossier
        track={"submission_level_single_agents":[{
          "candidate":"DRUGX","status":"PLAUSIBLE_CHALLENGE_HYPOTHESIS",
          "gene":"GENE1","objective_type":"bypass_or_compensate_downstream",
          "route_target":"TARGET1","route_desired_action":"decrease",
          "source_bound_evidence_state":"P2_CROSS_SOURCE_MECHANISTIC_SUPPORT",
          "source_bound_cross_source_bridge":True,
          "cross_source_bridge_rescue_sources":[{"pmid":"1","title":"rescue","license":"CC BY 4.0"}],
          "cross_source_bridge_drug_target_source":{"pmid":"2","title":"target","license":"reference-only"},
          "challenge_priority_score":10.0,"challenge_priority_score_source_adjusted":10.0,
          "public_label_safety_review":{"state":"HOLD"},
          "exposure_feasibility":{"state":"EXPOSURE_BORDERLINE_PROTEIN_BINDING"},
          "main_blockers":["same_drug_disease_model_rescue_not_demonstrated"],
          "clinical_recommendation":False}],
          "combination_hypotheses":[]}
        ps={"genes":[{"gene":"GENE1","variant_count":1,"effect_classes":["missense"],
                      "evidence_flags":{},"phase_state":"unresolved",
                      "mechanism_confirmed":False}]}
        dsg={"nodes":[{"id":"GENE:GENE1","genome_rank":1,"rank_score_not_diagnostic":10,
                       "phenotype_score":3}]}
        d=challenge_dossier.build(track,ps,dsg,{"version":"test"})
        self.assertEqual(d["primary_hypothesis"]["candidate"],"DRUGX")
        self.assertFalse(d["clinical_recommendation"])
        self.assertFalse(d["dose_recommendation"])
        self.assertFalse(d["clinical_safety_clearance"])
