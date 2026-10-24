"""Offline regression tests for transient EBI ChEMBL errors and verified resume."""
from __future__ import annotations
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase,main
from unittest.mock import patch
import io,json,os,urllib.error,urllib.parse,hashlib,importlib.util

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("v28_chembl_resilience",ROOT/"bootstrap_v28.py")
v28=importlib.util.module_from_spec(spec)
spec.loader.exec_module(v28)

URL="https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4&limit=1000&offset=0"
SCOPE="1234567890abcdef1234"

def payload(items=None,next_url=None):
    return {"molecules":[{"molecule_chembl_id":x} for x in (items or ["CHEMBL1"])],
            "page_meta":{"next":next_url}}

def as_response(value):
    return io.BytesIO(json.dumps(value).encode("utf-8"))

def error(code=500,headers=None):
    return urllib.error.HTTPError(URL,code,"test upstream unavailable",headers or {},None)

class RetryAndPageCache(TestCase):
    def setUp(self):
        rest=patch.dict(os.environ,{"RDMVA_CHEMBL_BACKEND":"rest"})
        rest.start()
        self.addCleanup(rest.stop)

    def test_public_source_release_includes_network_failure_regressions(self):
        manifest=(ROOT/"build_source_release_v31.py").read_text(encoding="utf-8")
        self.assertIn("code/test_v31_chembl_resilience.py",manifest)

    def test_single_upstream_http_500_retries_then_caches(self):
        with TemporaryDirectory() as folder,\
             patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),\
             patch.object(v28,"PUBLIC",Path(folder)),\
             patch.object(v28.time,"sleep") as sleep,\
             patch.object(v28.urllib.request,"urlopen",
                   side_effect=[error(500),as_response(payload())]) as network:
            self.assertEqual(v28.fetch_json(URL),payload())
            self.assertEqual(network.call_count,2)
            sleep.assert_called_once_with(2)
            cached=v28._chembl_cache_path(URL)
            self.assertTrue(cached.is_file())
            self.assertEqual(json.loads(cached.read_text())["schema"],
                             v28.CHEMBL_CACHE_SCHEMA)
            with patch.object(v28.urllib.request,"urlopen",
                               side_effect=AssertionError("no network on resume")):
                self.assertEqual(v28.fetch_json(URL),payload())

    def test_retry_after_obeyed_with_upper_bound(self):
        with patch.object(v28.time,"sleep") as pause,\
             patch.object(v28.urllib.request,"urlopen",
                          side_effect=[error(429,{"Retry-After":"19"}),
                                       as_response(payload())]):
            self.assertEqual(v28.fetch_json(URL),payload())
            pause.assert_called_once_with(19)

    def test_invalid_cache_sha_is_not_trusted_and_is_replaced(self):
        with TemporaryDirectory() as folder,\
             patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),\
             patch.object(v28,"PUBLIC",Path(folder)):
            path=v28._chembl_cache_path(URL)
            path.parent.mkdir(parents=True)
            path.write_text(json.dumps({
                "schema":v28.CHEMBL_CACHE_SCHEMA,
                "url":URL,"response_sha256":"f"*64,
                "response":payload(["TAMPERED"])}))
            with patch.object(v28.urllib.request,"urlopen",
                               return_value=as_response(payload(["REAL"]))):
                self.assertEqual(v28.fetch_json(URL),payload(["REAL"]))
            self.assertEqual(json.loads(path.read_text())["response"],
                             payload(["REAL"]))

    def test_permanent_404_does_not_retry_or_cache(self):
        with TemporaryDirectory() as folder,\
             patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),\
             patch.object(v28,"PUBLIC",Path(folder)),\
             patch.object(v28.time,"sleep") as pause,\
             patch.object(v28.urllib.request,"urlopen",
                          side_effect=error(404)) as request:
            with self.assertRaises(urllib.error.HTTPError):
                v28.fetch_json(URL)
            request.assert_called_once()
            pause.assert_not_called()
            self.assertFalse(v28._chembl_cache_path(URL).is_file())

    def test_exhausted_500_fails_closed_without_phantom_data(self):
        with TemporaryDirectory() as folder,\
             patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),\
             patch.object(v28,"PUBLIC",Path(folder)),\
             patch.object(v28.time,"sleep") as pause,\
             patch.object(v28.urllib.request,"urlopen",
                          side_effect=[error(500) for _ in range(5)]) as request:
            with self.assertRaisesRegex(RuntimeError,
                          "V31_CHEMBL_UPSTREAM_TEMPORARY_HTTP_EXHAUSTED"):
                v28.fetch_json(URL)
            self.assertEqual(request.call_count,5)
            self.assertEqual(pause.call_count,4)
            self.assertFalse(v28._chembl_cache_path(URL).is_file())

    def test_transient_timeout_recovered_with_bounded_retry(self):
        with patch.object(v28.time,"sleep") as pause,\
             patch.object(v28.urllib.request,"urlopen",side_effect=[
                 urllib.error.URLError("temporary network timeout"),
                 as_response(payload())]):
            self.assertEqual(v28.fetch_json(URL),payload())
            pause.assert_called_once_with(2)

    def test_untrusted_url_is_rejected_without_network(self):
        for url in (
            "http://www.ebi.ac.uk/chembl/api/data/molecule.json",
            "https://example.org/chembl/api/data/molecule.json",
            "https://www.ebi.ac.uk/admin/secrets"):
            with self.assertRaisesRegex(RuntimeError,"V31_CHEMBL_UNTRUSTED"):
                v28.fetch_json(url)

    def test_corrupted_json_or_invalid_cache_scope_cannot_pass(self):
        with patch.object(v28.time,"sleep") as pause,             patch.object(v28.urllib.request,"urlopen",
                          side_effect=[io.BytesIO(b"not json") for _ in range(5)]) as calls:
            with self.assertRaisesRegex(RuntimeError,"BAD_RESPONSE_EXHAUSTED"):
                v28.fetch_json(URL)
            self.assertEqual(calls.call_count,5)
            self.assertEqual(pause.call_count,4)
        with patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":"../bad"}):
            with self.assertRaisesRegex(RuntimeError,"UNSAFE_CACHE_SCOPE"):
                v28.fetch_json(URL)

    def test_http_200_json_error_envelope_does_not_create_fake_no_signal(self):
        fake={"detail":"Internal Server Error"}
        with TemporaryDirectory() as folder,             patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),             patch.object(v28,"PUBLIC",Path(folder)),             patch.object(v28.time,"sleep") as pause,             patch.object(v28.urllib.request,"urlopen",
                   side_effect=[as_response(fake) for _ in range(5)]) as requests:
            with self.assertRaisesRegex(RuntimeError,"BAD_RESPONSE_EXHAUSTED"):
                v28.fetch_json(URL)
            self.assertEqual(requests.call_count,5)
            self.assertEqual(pause.call_count,4)
            self.assertFalse(v28._chembl_cache_path(URL).is_file())

    def test_native_windows_installer_reports_upstream_outage_as_resumable(self):
        installer=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        launcher=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("Record 'UPSTREAM_CHEMBL_RETRY_LATER'",installer)
        self.assertIn("exit 11",installer)
        self.assertIn("if ($status -eq 11)",launcher)
        self.assertIn("RDMVA_JURY_PUBLIC_UPSTREAM_SERVICE_TEMPORARILY_UNAVAILABLE",launcher)
        self.assertIn("if ($Action -eq 'install') {exit 11}",launcher)

    def test_independent_capture_scope_does_not_reuse_previous_snapshot(self):
        with TemporaryDirectory() as folder,\
             patch.object(v28,"PUBLIC",Path(folder)):
            with patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),\
                 patch.object(v28.urllib.request,"urlopen",
                              return_value=as_response(payload(["CAPTURE_A"]))):
                self.assertEqual(v28.fetch_json(URL),payload(["CAPTURE_A"]))
            with patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":"abcdef0987654321abcd"}),\
                 patch.object(v28.urllib.request,"urlopen",
                              return_value=as_response(payload(["CAPTURE_B"]))):
                self.assertEqual(v28.fetch_json(URL),payload(["CAPTURE_B"]))

    def test_pagination_saves_and_restores_exact_pages_with_progress(self):
        a=URL
        b="https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit=1000&offset=1000"
        mockresp=[
          as_response(payload(["C1"],"/chembl/api/data/molecule.json?limit=1000&offset=1000")),
          as_response(payload(["C2"],None))]
        with TemporaryDirectory() as folder,\
             patch.dict(os.environ,{"RDMVA_CHEMBL_CACHE_SCOPE":SCOPE}),\
             patch.object(v28,"PUBLIC",Path(folder)),\
             patch.object(v28.urllib.request,"urlopen",side_effect=mockresp) as net:
            self.assertEqual([x["molecule_chembl_id"] for x in
                v28.paged("molecule",{"max_phase":4},"molecules")],["C1","C2"])
            self.assertEqual(net.call_count,2)
            with patch.object(v28.urllib.request,"urlopen",
                              side_effect=AssertionError("must use cached pages")):
                self.assertEqual(len(v28.paged("molecule",{"max_phase":4},"molecules")),2)

    def test_pagination_loop_is_rejected(self):
        first="https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4&limit=1000&offset=0"
        with patch.object(v28,"fetch_json",return_value=payload(["C1"],first)):
            with self.assertRaisesRegex(RuntimeError,"V31_CHEMBL_PAGINATION_CYCLE"):
                v28.paged("molecule",{"max_phase":4},"molecules")

class SnapshotInterruptedResumeTests(TestCase):
    def test_verified_prior_public_downloads_and_chembl_capture_scope_persist(self):
        # Test one failed acquisition, then resume, without a single real
        # patient input or any network call.
        import importlib.util,sys
        spec=importlib.util.spec_from_file_location(
            "rdmva_current_chembl_resume_fixture",
            ROOT/"bootstrap_current_public_v31.py")
        live=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(live)
        with TemporaryDirectory() as tmp:
            folder=Path(tmp)
            manifest=folder/"original_static.json"
            regulatory=folder/"original_regulatory.json"
            raw=folder/"current_references"
            audit=folder/"audit"
            manifest.write_text(json.dumps({
                "strict_snapshot":True,"version":"SYNTHETIC_ONLY",
                "static_files":[{"path":"fda/synthetic_public.json",
                                 "url":"https://www.fda.gov/synthetic",
                                 "bytes":1,"sha256":"original_hash"}],
                "generated_files":[{"path":"chembl/catalog_phase4.json",
                                    "method":"synthetic_catalog",
                                    "bytes":1,"sha256":"original_hash"}]
            }),encoding="utf-8")
            regulatory.write_text(json.dumps({
                "schema":"V283_REGULATORY_SCHEMA_1","files":{}
            }),encoding="utf-8")
            total_downloads=[]
            captured_scopes=[]
            calls=[]
            def fake_public_download(url,dst,old_bytes,name):
                total_downloads.append(name)
                dst.parent.mkdir(parents=True,exist_ok=True)
                dst.write_bytes(b'{"synthetic":true}')
                return {"source_url":url,
                        "sha256":live.sha(dst),
                        "bytes":dst.stat().st_size,
                        "retrieved_utc":live.utc()}
            def fake_chembl_build(model):
                captured_scopes.append(os.environ.get("RDMVA_CHEMBL_CACHE_SCOPE"))
                calls.append(1)
                if len(calls)==1:
                    raise urllib.error.HTTPError(URL,500,"synthetic backend 500",{},None)
                for entry in model["generated_files"]:
                    target=raw/entry["path"]
                    target.parent.mkdir(parents=True,exist_ok=True)
                    target.write_bytes(b'[{"chembl":"synthetic"}]')
                    entry["bytes"]=target.stat().st_size
                    entry["sha256"]=live.sha(target)
            with patch.object(live,"STATIC_ORIGINAL",manifest),                 patch.object(live,"REG_ORIGINAL",regulatory),                 patch.object(live,"RAW",raw),                 patch.object(live,"AUDIT",audit),                 patch.object(live,"RECEIPT",audit/"CURRENT_PUBLIC_DATA_V31.json"),                 patch.object(live,"STATIC_LIVE",audit/"public_data_manifest_v28_current.json"),                 patch.object(live,"REG_LIVE",audit/"regulatory_snapshot_v283_current.json"),                 patch.object(live,"curl_download",side_effect=fake_public_download),                 patch.object(v28,"build_chembl",side_effect=fake_chembl_build),                 patch.dict(sys.modules,{"bootstrap_v28":v28}), patch.dict(os.environ,{"RDMVA_CHEMBL_BACKEND":"rest"}):
                with self.assertRaises(urllib.error.HTTPError):
                    live.bootstrap()
                failed=live.read(live.RECEIPT)
                self.assertEqual(failed["state"],"FAILED")
                self.assertEqual(len(total_downloads),1)
                self.assertFalse(live.verify())
                live.bootstrap()
                resumed=live.read(live.RECEIPT)
                self.assertEqual(resumed["state"],"READY")
                self.assertEqual(len(total_downloads),1,
                    "all prior publicly downloaded resources must be reused")
                self.assertEqual(len(captured_scopes),2)
                self.assertEqual(captured_scopes[0],captured_scopes[1])
                self.assertTrue(live.verify())
                self.assertEqual(
                    resumed["capture_started_utc"],failed["capture_started_utc"])
                self.assertEqual(
                    resumed["chembl_capture_scope"],failed["chembl_capture_scope"])

if __name__=="__main__":main()
