"""Synthetic-only tests for optional live paid API model preflight."""
from pathlib import Path
from unittest.mock import patch
import importlib.util,json,tempfile,unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("jury_model_check_v31",
                                             ROOT/"JURY_CHECK_AI_V31.py")
modelcheck=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(modelcheck)

class JuryModelValidationTests(unittest.TestCase):
    def test_v30_and_v28_all_configured_model_roles_are_enumerated(self):
        m=modelcheck.required_models()
        for role in ("v30_default","v30_critic","v30_literature",
                     "v2810_scout","v2810_verify","v2810_final"):
            self.assertIn(role,m)
            self.assertTrue(isinstance(m[role],str) and m[role])
        self.assertGreaterEqual(len(set(m.values())),2)

    def test_mocked_public_synthetic_responses_pass_without_secrets(self):
        ids=sorted(set(modelcheck.required_models().values()))
        payloads=[]
        def fake_request(method,url,key,body=None,timeout=60):
            self.assertEqual(key,"TEST_ONLY_FAKE_NOT_A_REAL_SECRET")
            if method=="GET":
                self.assertTrue(url.endswith("/models"))
                return {"data":[{"id":m} for m in ids]}
            self.assertEqual(method,"POST")
            self.assertTrue(url.endswith("/responses"))
            self.assertEqual(body["store"],False)
            self.assertNotIn("variant",json.dumps(body).lower())
            self.assertNotIn("hpo",json.dumps(body).lower())
            self.assertIn("Synthetic public connectivity test",body["input"])
            payloads.append(body)
            return {"id":"resp_synthetic_ABC","status":"completed",
                    "model":body["model"],
                    "usage":{"input_tokens":12,"output_tokens":2}}
        with tempfile.TemporaryDirectory() as tmp:
            report=Path(tmp)/"check.json"
            with patch.object(modelcheck,"api_key",
                        return_value="TEST_ONLY_FAKE_NOT_A_REAL_SECRET"),\
                 patch.object(modelcheck,"request_json",side_effect=fake_request):
                good,receipt=modelcheck.validate(test_responses=True,report_path=report)
            self.assertTrue(good)
            self.assertTrue(report.is_file())
            self.assertEqual(len(payloads),len(ids))
            text=report.read_text()
            self.assertNotIn("TEST_ONLY_FAKE_NOT_A_REAL_SECRET",text)
            self.assertTrue(receipt["passed"])
            self.assertFalse(receipt["tested_using_patient_data"])

    def test_failed_or_unlisted_live_model_never_gets_pass(self):
        ids=sorted(set(modelcheck.required_models().values()))
        def fake_request(method,url,key,body=None,timeout=60):
            if method=="GET": return {"data":[{"id":m} for m in ids[:-1]]}
            if body["model"]==ids[-1]:
                return {"id":"resp_incomplete","status":"incomplete","model":ids[-1]}
            return {"id":"resp_other","status":"completed","model":body["model"]}
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(modelcheck,"api_key",return_value="fake"),\
                 patch.object(modelcheck,"request_json",side_effect=fake_request):
                ok,receipt=modelcheck.validate(test_responses=True,
                                    report_path=Path(tmp)/"bad.json")
                self.assertFalse(ok)
                self.assertFalse(receipt["passed"])
                ok_listing,_=modelcheck.validate(test_responses=False)
                self.assertFalse(ok_listing)

if __name__=="__main__":
    unittest.main()
