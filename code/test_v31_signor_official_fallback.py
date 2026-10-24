"""Official Oct2026 SIGNOR fallback; offline, no AI, no restricted inputs."""
from __future__ import annotations
import contextlib,hashlib,importlib.util,io,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
def load():
    spec=importlib.util.spec_from_file_location("v312_signor_bootstrap",
                                                ROOT/"bootstrap_current_public_v31.py")
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

class OfficialSignorFallbackTests(unittest.TestCase):
    def test_exact_seed_release_and_attribution(self):
        x=load()
        self.assertEqual(x.SIGNOR_SEED.stat().st_size,6722506)
        self.assertEqual(x.sha(x.SIGNOR_SEED),x.SIGNOR_SEED_SHA256)
        self.assertEqual(x.SIGNOR_OFFICIAL_URL,
                         "https://signor.uniroma2.it/releases/Oct2026_release.txt")
        notice=(ROOT/"reference_seed/signor/ATTRIBUTION.md").read_text(encoding="utf-8")
        self.assertIn("CC BY 4.0",notice)
        self.assertIn("Lo Surdo",notice)
    def test_official_pinned_fallback_restores_exact_original(self):
        x=load()
        with tempfile.TemporaryDirectory() as folder,contextlib.redirect_stdout(io.StringIO()) as log:
            dst=Path(folder)/"signor"/"Oct2026_release.txt"
            meta=x.recover_official_signor_snapshot(
                x.SIGNOR_RELEASE,x.SIGNOR_OFFICIAL_URL,dst,x.SIGNOR_RELEASE_BYTES)
            self.assertEqual(dst.stat().st_size,x.SIGNOR_RELEASE_BYTES)
            self.assertEqual(x.sha(dst),x.SIGNOR_RELEASE_SHA256)
            self.assertEqual(meta["sha256"],x.SIGNOR_RELEASE_SHA256)
            self.assertTrue(meta["official_pinned_fallback"])
            self.assertTrue(meta["upstream_capture_attempted_and_failed"])
            self.assertEqual(meta["source_url"],x.SIGNOR_OFFICIAL_URL)
            self.assertIn("V312_SIGNOR_PINNED_RELEASE_RESTORED",log.getvalue())
    def test_wrong_source_or_release_rejected_without_creating_file(self):
        x=load()
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/"file.txt"
            for name,url,length in (
                    ("other/File.tsv",x.SIGNOR_OFFICIAL_URL,x.SIGNOR_RELEASE_BYTES),
                    (x.SIGNOR_RELEASE,"https://evil.invalid/mirror",x.SIGNOR_RELEASE_BYTES),
                    (x.SIGNOR_RELEASE,x.SIGNOR_OFFICIAL_URL,x.SIGNOR_RELEASE_BYTES+1)):
                with self.assertRaisesRegex(RuntimeError,"REFUSES_UNPINNED"):
                    x.recover_official_signor_snapshot(name,url,p,length)
                self.assertFalse(p.exists())
    def test_corrupted_seed_fail_closed(self):
        x=load()
        with tempfile.TemporaryDirectory() as folder:
            seed=Path(folder)/"fake.gz"
            seed.write_bytes(b"broken")
            dest=Path(folder)/"target.tsv"
            with patch.object(x,"SIGNOR_SEED",seed):
                with self.assertRaisesRegex(RuntimeError,"MISSING_OR_CORRUPT"):
                    x.recover_official_signor_snapshot(
                        x.SIGNOR_RELEASE,x.SIGNOR_OFFICIAL_URL,
                        dest,x.SIGNOR_RELEASE_BYTES)
            self.assertFalse(dest.exists())
    def test_installer_records_real_upstream_failure_as_pinned_exception(self):
        x=load()
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            raw=root/"third_party"
            audit=root/"receipt"
            audit.mkdir()
            mock_files=[]
            def curl_stub(url,dst,old_bytes,name):
                mock_files.append(name)
                if name==x.SIGNOR_RELEASE:
                    raise RuntimeError("V31_CURRENT_DOWNLOAD_FAILED:"+name+
                                       ":curl_exit=7:http=000")
                if len(mock_files)>3:
                    raise RuntimeError("SYNTHETIC_STOP_AFTER_SIGNOR")
                dst.parent.mkdir(parents=True,exist_ok=True)
                data=(b"synthetic only "+name.encode())*3
                dst.write_bytes(data)
                return {"bytes":len(data),
                        "sha256":hashlib.sha256(data).hexdigest(),
                        "source_url":url,"retrieved_utc":x.utc(),
                        "download_transport":"synthetic"}
            with patch.object(x,"RAW",raw),\
                 patch.object(x,"AUDIT",audit),\
                 patch.object(x,"RECEIPT",audit/"CURRENT_PUBLIC_DATA_V31.json"),\
                 patch.object(x,"STATIC_LIVE",audit/"static.json"),\
                 patch.object(x,"REG_LIVE",audit/"reg.json"),\
                 patch.object(x,"curl_download",side_effect=curl_stub),\
                 contextlib.redirect_stdout(io.StringIO()) as output:
                with self.assertRaisesRegex(RuntimeError,"SYNTHETIC_STOP_AFTER_SIGNOR"):
                    x.bootstrap()
                receipt=json.loads(x.RECEIPT.read_text(encoding="utf-8"))
                self.assertEqual(receipt["state"],"FAILED")
                self.assertEqual(receipt["files"][x.SIGNOR_RELEASE]["sha256"],
                                 x.SIGNOR_RELEASE_SHA256)
                self.assertEqual(receipt["pinned_public_source_exception"]["reason"],
                                 "official_server_unavailable")
                self.assertIn("NOT a fresh",output.getvalue())
                self.assertEqual(x.sha(raw/x.SIGNOR_RELEASE),x.SIGNOR_RELEASE_SHA256)
    def test_source_builder_has_only_one_binary_data_allowlist_exception(self):
        builder=(ROOT/"build_source_release_v31.py").read_text(encoding="utf-8")
        self.assertIn('SIGNED_SIGNOR_SEED="reference_seed/signor/Oct2026_release.txt.gz"',builder)
        self.assertIn("V312_PUBLIC_SIGNOR_SEED_PIN_MISMATCH",builder)
        self.assertIn('"version":"3.1.4"',builder)
        self.assertIn('"contains_patient_raw_genome":False',builder)
        self.assertIn('"contains_api_keys_or_secrets":False',builder)

if __name__=="__main__":unittest.main()
