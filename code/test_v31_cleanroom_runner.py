"""Deterministic independent verifier tests with synthetic-only fixtures."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("v31_verify", ROOT/"verify_challenge_run_v31.py")
module=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def put(path, obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj),encoding="utf-8")

class V31CleanroomVerifierTests(unittest.TestCase):
    def setup_fixture(self,tmp,with_receipt):
        run=Path(tmp)
        put(run/"FINAL_RESULT_V31_NATIVE.json",{
            "all_approved_entities_screened":7523,
            "track2_submission_hypotheses":"TRACK2_HYPOTHESES_V31_FINAL.json",
            "patient_raw_data_sent_to_ai":False,
            "clinical_recommendation":False,
            "openai_enabled":True,
        })
        put(run/"TRACK2_HYPOTHESES_V31_FINAL.json",{
            "submission_level_single_agents":[{"candidate":"OTHERDRUG","gene":"OTHERGENE",
                                                "clinical_recommendation":False}],
            "combination_hypotheses":[],
        })
        put(run/"STATUS_V30_NATIVE.json",{
            "completed":True,"final_result":"FINAL_RESULT_V31_NATIVE.json"})
        put(run/"V31_THIRD_PARTY_RIGHTS_AUDIT.json",{"schema":"TEST"})
        put(run/"GENOMIC_REQUIRED_EVIDENCE_V31.json",{
            "schema":"RDMVA_CASE_GENOMIC_EVIDENCE_AUDIT_V31_1",
            "required_science_pass":True})
        put(run/"test_ai_response.json",{
            "request_sha256":"a"*64,"usage":{"input_tokens":1,"output_tokens":1},
            "response_id":"resp_test_QA_example","model":"gpt-6-luna","store":False
        } if with_receipt else {})
        put(run/"expected.json",{
            "expected_minimum_approved_entities_screened":7000})
        return run

    def test_technical_pass_does_not_require_known_winning_drug(self):
        with tempfile.TemporaryDirectory() as temp:
            run=self.setup_fixture(temp,True)
            ok,_,count=module.verify(run,run/"expected.json",True)
            self.assertTrue(ok)
            self.assertEqual(count,1)
            context=module.load(run/"RUN_REPRODUCIBILITY_CONTEXT_V31.json")
            self.assertEqual(context["primary_hypothesis"]["candidate"],"OTHERDRUG")
            self.assertIn("No particular drug is required",context["note"])

    def test_online_ai_requires_real_response_receipt(self):
        with tempfile.TemporaryDirectory() as temp:
            run=self.setup_fixture(temp,False)
            ok,_,count=module.verify(run,run/"expected.json",True)
            self.assertFalse(ok)
            self.assertEqual(count,0)

    def test_critic_summary_without_store_must_reference_verified_receipt(self):
        with tempfile.TemporaryDirectory() as temp:
            run=self.setup_fixture(temp,True)
            put(run/"AI_RESEARCH_CRITIC_V31.json",{
                "schema":"RDMVA_OPENAI_CRITIC_V30_1",
                "response_id":"resp_test_QA_example",
                "model":"gpt-6-luna","value":{"hypothesis":False}})
            ok,_,count=module.verify(run,run/"expected.json",True)
            self.assertTrue(ok)
            self.assertEqual(count,1)
            audit=module.load(run/"CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json")
            self.assertEqual(audit["ai_provenance_metadata"]["summary_references_verified"],1)
            self.assertEqual(audit["ai_provenance_metadata"]["response_records_without_store_false"],0)

    def test_unbacked_summary_response_id_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            run=self.setup_fixture(temp,True)
            put(run/"AI_RESEARCH_CRITIC_V31.json",{
                "response_id":"resp_QA_NOT_IN_ACTUAL_RECEIPT",
                "model":"gpt-6-luna"})
            ok,_,count=module.verify(run,run/"expected.json",True)
            self.assertFalse(ok)
            self.assertEqual(count,1)
            a=module.load(run/"CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json")
            self.assertFalse(a["checks"]["ai_response_store_false"])
            self.assertEqual(a["ai_provenance_metadata"]["orphan_summary_response_references"],1)

    def test_direct_response_store_true_or_absent_is_rejected(self):
        for value in (True,None):
            with self.subTest(store=value),tempfile.TemporaryDirectory() as temp:
                run=self.setup_fixture(temp,True)
                receipt=module.load(run/"test_ai_response.json")
                if value is None:receipt.pop("store")
                else:receipt["store"]=value
                put(run/"test_ai_response.json",receipt)
                ok,_,count=module.verify(run,run/"expected.json",True)
                self.assertFalse(ok)
                self.assertEqual(count,1)

    def test_summary_cannot_claim_store_true_for_a_verified_request(self):
        with tempfile.TemporaryDirectory() as temp:
            run=self.setup_fixture(temp,True)
            put(run/"AI_RESEARCH_CRITIC_V31.json",{
                "response_id":"resp_test_QA_example","store":True})
            ok,_,count=module.verify(run,run/"expected.json",True)
            self.assertFalse(ok)
            self.assertEqual(count,1)

if __name__=="__main__":unittest.main()
