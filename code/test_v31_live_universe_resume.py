"""Regression: current official data is independently source-verified.

Synthetic test inputs only, never Challenge patient files.
"""
import importlib.util,hashlib,json,os,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from test_v31_cleanroom_runner import V31CleanroomVerifierTests,module as verifier,put

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("rdmva_e2e_live_test",ROOT/"reproduce_v31_pipeline.py")
runner=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)
from v2810 import master

def sha(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class LiveUSCountTests(unittest.TestCase):
 def test_historic_release_must_still_be_exact(self):
  with patch.dict(os.environ,{"RDMVA_CURRENT_PUBLIC_DATA":"0"}):
   self.assertEqual(master.verify_universe_profile(2414),"FROZEN_2026_10_09")
   with self.assertRaisesRegex(RuntimeError,"V283_BASELINE_US_DRIFT"):
    master.verify_universe_profile(2415)

 def test_verified_live_release_accepts_observed_count_but_detects_tampering(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)
   pub=root/"public_inputs/current_v31"
   pub.mkdir(parents=True)
   fda=root/"new_fda.zip";fda.write_bytes(b"synthetic FDA")
   purple=root/"new_purple.csv";purple.write_bytes(b"synthetic Purple")
   put(pub/"CURRENT_PUBLIC_DATA_V31.json",
       {"state":"READY","profile":"current_official_sources"})
   put(pub/"public_data_manifest_v28_current.json",{
      "static_files":[
       {"path":"fda/drugsatfda.zip","bytes":fda.stat().st_size,"sha256":sha(fda)},
       {"path":"fda/purplebook_2026-09.csv","bytes":purple.stat().st_size,"sha256":sha(purple)}]})
   with patch.dict(os.environ,{"RDMVA_CURRENT_PUBLIC_DATA":"1"}),\
        patch.object(master,"ROOT",root),patch.object(master,"FDAZIP",fda),\
        patch.object(master,"PURPLE",purple):
    self.assertEqual(master.verify_universe_profile(2415),"CURRENT_OFFICIAL_SOURCES")
    with self.assertRaisesRegex(RuntimeError,"UNEXPECTEDLY_SMALL"):
     master.verify_universe_profile(1999)
    purple.write_bytes(b"tampered Purple")
    with self.assertRaisesRegex(RuntimeError,"SOURCE_HASH_MISMATCH"):
     master.verify_universe_profile(2415)

class LiveV31VerificationTests(unittest.TestCase):
 def test_new_source_count_is_measured_not_frozen(self):
  with tempfile.TemporaryDirectory() as tmp:
   run=V31CleanroomVerifierTests().setup_fixture(tmp,True)
   original=json.loads((run/"FINAL_RESULT_V31_NATIVE.json").read_text())
   original["all_approved_entities_screened"]=7533
   put(run/"FINAL_RESULT_V31_NATIVE.json",original)
   put(run/"master/summary.json",{
      "entities":7533,"us_snapshot_entity_count":2415,"new_entity_labels":5118,
      "dataset_profile":"CURRENT_OFFICIAL_SOURCES"})
   put(run/"E2E_EXECUTION_V31.json",{
      "public_data_profile":"CURRENT_OFFICIAL_SOURCES",
      "current_public_data_receipt_sha256":"0"*64})
   put(run/"expected.json",{
      "expected_minimum_approved_entities_screened":9000,
      "reference_snapshot_approved_entities":7523})
   ok,_,num=verifier.verify(run,run/"expected.json",True)
   self.assertTrue(ok)
   self.assertEqual(num,1)
   context=json.loads((run/"RUN_REPRODUCIBILITY_CONTEXT_V31.json").read_text())
   self.assertEqual(context["approved_entities"]["observed"],7533)
   self.assertTrue(context["current_public_data_profile"])
   corrupted=json.loads((run/"master/summary.json").read_text())
   corrupted["entities"]=7600
   put(run/"master/summary.json",corrupted)
   ok,_,_=verifier.verify(run,run/"expected.json",True)
   self.assertFalse(ok)

class VerifiedVCFIntakeResumeTests(unittest.TestCase):
 def _create(self,tmp):
  t=Path(tmp)
  vcf=t/"test.vcf.gz";vcf.write_bytes(b"synthetic-private-genome")
  doc=t/"clinical.docx";doc.write_bytes(b"synthetic-private-docx")
  run=t/"run";intake=run/"intake";intake.mkdir(parents=True)
  source=[{"file":p.name,"sha256":sha(p)} for p in (vcf,doc)]
  case=intake/"case_model_v2.json"
  put(case,{"official_source_hashes":source})
  put(intake/"intake_manifest.json",{
       "inputs":source,"outputs":{"case_model_v2.json":sha(case)}})
  receipt={"status":"FAILED","stages":{"stage1_native_core":{"exit_code":1}},
       "input_file_type":"authorized_vcf","input_file_sha256":sha(vcf),
       "public_data_profile":"FROZEN_CHALLENGE_BASELINE"}
  return run,vcf,doc,case,receipt
 def test_original_vcf_and_phenotype_hashes_are_required(self):
  with tempfile.TemporaryDirectory() as tmp:
   run,vcf,doc,case,receipt=self._create(tmp)
   with patch.dict(os.environ,{"RDMVA_CURRENT_PUBLIC_DATA":"0"}):
    self.assertEqual(runner.verify_failed_intake_for_resume(run,vcf,doc,receipt),case)
    doc.write_bytes(b"changed")
    with self.assertRaisesRegex(RuntimeError,"INTAKE_INPUT_PROVENANCE_MISMATCH"):
     runner.verify_failed_intake_for_resume(run,vcf,doc,receipt)
 def test_broken_intake_artifact_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   run,vcf,doc,case,receipt=self._create(tmp)
   case.write_text('{"tampered":true}',encoding="utf-8")
   with patch.dict(os.environ,{"RDMVA_CURRENT_PUBLIC_DATA":"0"}):
    with self.assertRaisesRegex(RuntimeError,"INTAKE_INPUT_PROVENANCE_MISMATCH"):
     runner.verify_failed_intake_for_resume(run,vcf,doc,receipt)

if __name__=="__main__":unittest.main()
