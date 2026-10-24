"""Live public-source acquisition checks: no network, no patient data."""
import importlib.util,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def module_at(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

live=module_at(ROOT/"bootstrap_current_public_v31.py","rdmva_current_sources_test")
v28=module_at(ROOT/"bootstrap_v28.py","rdmva_v28_source_test")

class CurrentSnapshotTests(unittest.TestCase):
    def test_only_https_allowlisted_sources(self):
        self.assertEqual(live.validate_source_url(
            "https://www.fda.gov/media/89850/download?attachment").hostname,
            "www.fda.gov")
        for url in ("http://www.fda.gov/file.zip",
                    "https://evil.invalid/file.zip",
                    "file:///C:/secret"):
            with self.assertRaises(RuntimeError):live.validate_source_url(url)

    def test_content_type_and_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"sample.json"
            p.write_text('{"safe":true,"padding":"'+('x'*80)+'"}',encoding="utf-8")
            live.validate_download(p,"sample.json")
            receipt={"bytes":p.stat().st_size,"sha256":live.sha(p)}
            self.assertTrue(live.validate_record(p,receipt))
            p.write_text('{"safe":false,"padding":"'+('x'*80)+'"}',encoding="utf-8")
            self.assertFalse(live.validate_record(p,receipt))

    def test_current_generated_hash_captured(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"generated.json"
            spec={"path":"chembl/generated.json","bytes":1,"sha256":"incorrect",
                  "capture_current":True}
            v28.write_exact(p,'[{"id":"x"}]',spec)
            self.assertEqual(spec["bytes"],p.stat().st_size)
            self.assertEqual(spec["sha256"],v28.sha256(p))

    def test_original_frozen_manifest_not_updated_by_current_download(self):
        st=json.loads((ROOT/"config/public_data_manifest_v28.json").read_text(encoding="utf-8"))
        r=json.loads((ROOT/"config/regulatory_snapshot_v283.json").read_text(encoding="utf-8"))
        self.assertEqual(st["version"],"V28_PUBLIC_SNAPSHOT_2026-10-09")
        self.assertEqual(r["snapshot_date"],"2026-10-10")

    def test_regulatory_bootstrap_manifest_follows_live_profile(self):
        command=[sys.executable,"-c",
                 "import bootstrap_regulatory_v283 as r; print(r.SPEC.relative_to(r.ROOT).as_posix())"]
        for profile,expected in (
            (None,"config/regulatory_snapshot_v283.json"),
            ("1","public_inputs/current_v31/regulatory_snapshot_v283_current.json"),
        ):
            env=os.environ.copy()
            env.pop("RDMVA_CURRENT_PUBLIC_DATA",None)
            if profile is not None:env["RDMVA_CURRENT_PUBLIC_DATA"]=profile
            value=subprocess.check_output(command,cwd=ROOT,env=env,text=True,timeout=20).strip()
            self.assertEqual(value,expected)

    def test_regulatory_bootstrap_validates_current_manifest_contents(self):
        bootstrap=module_at(ROOT/"bootstrap_regulatory_v283.py","rdmva_regulatory_bootstrap_test")
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            data=root/"approved.json"
            data.write_text('{"drug":"synthetic"}',encoding="utf-8")
            manifest=root/"current-regulatory.json"
            manifest.write_text(json.dumps({
                "schema":"V283_REGULATORY_SCHEMA_1",
                "snapshot_date":"2026-10-19",
                "files":{"approved.json":{
                    "bytes":data.stat().st_size,
                    "sha256":live.sha(data),
                    "url":"https://health-products.canada.ca/synthetic"
                }}}),encoding="utf-8")
            with patch.object(bootstrap,"SPEC",manifest),patch.object(bootstrap,"DIR",root),patch.object(sys,"argv",["bootstrap_regulatory_v283.py","--verify-only"]):
                bootstrap.main()
            data.write_text('{"drug":"changed"}',encoding="utf-8")
            with patch.object(bootstrap,"SPEC",manifest),patch.object(bootstrap,"DIR",root),patch.object(sys,"argv",["bootstrap_regulatory_v283.py","--verify-only"]):
                with self.assertRaises(RuntimeError):bootstrap.main()

if __name__=="__main__":unittest.main()
