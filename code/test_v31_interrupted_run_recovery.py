"""Verify interrupted-run provenance and public ChEMBL cache recovery gates."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import json,os,tempfile,unittest
from test_v31_live_universe_resume import VerifiedVCFIntakeResumeTests,runner,put,sha
from v30 import fast_screen

class InterruptedRunRecoveryTests(unittest.TestCase):
 def test_new_run_verified_intake_and_untouched_previous(self):
  with tempfile.TemporaryDirectory() as td:
   t=Path(td)
   prior,vcf,doc,case,initial=VerifiedVCFIntakeResumeTests()._create(td)
   put(prior/"E2E_EXECUTION_V31_PREVIOUS_FAILED.json",initial)
   current={**initial,
     "status":"RUNNING","resumed_with_verified_original_intake":True}
   put(prior/"E2E_EXECUTION_V31.json",current)
   put(prior/"STATUS_V30_NATIVE.json",{"current_stage":"V30_UNIVERSAL_SCREEN"})
   original=sha(prior/"E2E_EXECUTION_V31.json")
   with patch.dict(os.environ,{"RDMVA_CURRENT_PUBLIC_DATA":"0"}),\
        patch.object(runner.subprocess,"run",return_value=SimpleNamespace(returncode=0,stdout="[]")):
    recovered=t/"recovered_new_run"
    receipt=runner.prepare_recovery_from_interrupted_run(prior,recovered,vcf,doc)
    self.assertTrue((recovered/"intake/case_model_v2.json").is_file())
    self.assertEqual(receipt["source_status"],"RUNNING")
    self.assertFalse(receipt["fresh_end_to_end"])
    self.assertEqual(sha(prior/"E2E_EXECUTION_V31.json"),original)
    self.assertEqual(sha(case),sha(recovered/"intake/case_model_v2.json"))
    with self.assertRaisesRegex(RuntimeError,"REFUSES_EXISTING_OUTPUT"):
     runner.prepare_recovery_from_interrupted_run(prior,recovered,vcf,doc)
 def test_recovery_rejects_active_previous_run(self):
  with tempfile.TemporaryDirectory() as td:
   prior,vcf,doc,case,initial=VerifiedVCFIntakeResumeTests()._create(td)
   put(prior/"E2E_EXECUTION_V31_PREVIOUS_FAILED.json",initial)
   put(prior/"E2E_EXECUTION_V31.json",{
     **initial,"status":"RUNNING","resumed_with_verified_original_intake":True})
   put(prior/"STATUS_V30_NATIVE.json",{"current_stage":"V30_UNIVERSAL_SCREEN"})
   process=[{"ProcessId":123,"CommandLine":
         "python -m v30.pipeline --out "+str(prior)}]
   with patch.dict(os.environ,{"RDMVA_CURRENT_PUBLIC_DATA":"0"}),\
        patch.object(runner.subprocess,"run",return_value=SimpleNamespace(returncode=0,stdout=json.dumps(process))):
    with self.assertRaisesRegex(RuntimeError,"ORIGINAL_PROCESS_STILL_ACTIVE"):
     runner.prepare_recovery_from_interrupted_run(prior,Path(td)/"new",vcf,doc)
   self.assertFalse((Path(td)/"new").exists())

class PublicChemblCacheTests(unittest.TestCase):
 def _fixture(self,td):
  old=Path(td)/"old";new=Path(td)/"new"
  summary={"entities":7533,"products":233433,
           "baseline_us_entities":2415,"dataset_profile":"CURRENT_OFFICIAL_SOURCES",
           "source_hashes":{"drugsatfda":"a"*64,"chembl_catalogue":"b"*64}}
  for run in (old,new):
   put(run/"master/summary.json",summary)
   (run/"master/entities.jsonl").write_text(
     '{"canonical_label":"SYNTHETIC"}\n',encoding="utf-8")
   put(run/"NATIVE_ROUTES_V31.json",{"routes":[{"route_id":"r1","desired_action":"increase"}]})
  oldsrc=old/"screen/sources"
  put(oldsrc/"act_CHEMBL123.json",[{"target_chembl_id":"CHEMBL123","activity_id":"a1"}])
  put(oldsrc/"structures.json",{"CHEMBL123":{"smiles":"CC"}})
  put(oldsrc/"target_index.json",{"ACC":["CHEMBL123"]})
  put(oldsrc/"route_activities.json",[{"route_id":"r1","target_chembl_id":"CHEMBL123"}])
  return old,new
 def test_exact_routes_reuse_and_sha_audit(self):
  with tempfile.TemporaryDirectory() as t:
   old,new=self._fixture(t)
   a=fast_screen.restore_public_chembl_cache(old,new)
   self.assertEqual(a["reused_public_target_files"],1)
   self.assertTrue(a["routes_identical"])
   self.assertTrue((new/"screen/sources/route_activities.json").exists())
   self.assertTrue((new/"V31_CHEMBL_CACHE_REUSE_AUDIT.json").exists())
 def test_new_routes_do_not_reuse_route_dependent_caches(self):
  with tempfile.TemporaryDirectory() as t:
   old,new=self._fixture(t)
   put(new/"NATIVE_ROUTES_V31.json",{"routes":[{"route_id":"r2","desired_action":"decrease"}]})
   a=fast_screen.restore_public_chembl_cache(old,new)
   self.assertEqual(a["reused_public_target_files"],1)
   self.assertFalse(a["routes_identical"])
   self.assertFalse((new/"screen/sources/route_activities.json").exists())
 def test_different_master_refuses_every_cache(self):
  with tempfile.TemporaryDirectory() as t:
   old,new=self._fixture(t)
   s=json.loads((new/"master/summary.json").read_text())
   s["products"]+=1;put(new/"master/summary.json",s)
   with self.assertRaisesRegex(RuntimeError,"UNIVERSE_SNAPSHOT_MISMATCH"):
    fast_screen.restore_public_chembl_cache(old,new)
   self.assertFalse((new/"screen/sources").exists())
 def test_private_or_wrong_target_cache_is_rejected(self):
  with tempfile.TemporaryDirectory() as t:
   old,new=self._fixture(t)
   put(old/"screen/sources/act_CHEMBL123.json",
       [{"target_chembl_id":"CHEMBL999","patient_genome":"illegal"}])
   with self.assertRaisesRegex(RuntimeError,"PUBLIC_TARGET_PAYLOAD_INVALID"):
    fast_screen.restore_public_chembl_cache(old,new)

if __name__=="__main__":unittest.main()
