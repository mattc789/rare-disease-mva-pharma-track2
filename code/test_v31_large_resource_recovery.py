"""Regression tests for repaired offline VEP and large gnomAD integrity checks."""
from pathlib import Path
from unittest.mock import patch
import hashlib,os,re,tempfile,unittest
from v2810 import resource_bootstrap as res

ROOT=Path(__file__).resolve().parents[1]

class LargeDownloadIntegrityTests(unittest.TestCase):
    def test_windows_mounted_volume_uses_streamed_sha256(self):
        if os.name!="nt":self.skipTest("Windows drive test")
        with tempfile.TemporaryDirectory(dir=str(ROOT)) as folder:
            p=Path(folder)/"synthetic-chromosome.bgz"
            p.write_bytes(b"synthetic test data"*100000)
            as_wsl="/mnt/"+p.drive[0].lower()+p.as_posix()[2:]
            self.assertEqual(res._sha256(as_wsl),hashlib.sha256(p.read_bytes()).hexdigest())
            # Never call WSL sha256sum for mounted Windows storage.
            with patch.object(res,"_run_wsl",side_effect=RuntimeError("WSL should not be invoked")):
                self.assertEqual(res._sha256(as_wsl),hashlib.sha256(p.read_bytes()).hexdigest())

    def test_unix_native_sha_has_long_timeout(self):
        with patch.object(res,"_run_wsl",return_value=("c"*64+"  /home/test/file")) as call:
            self.assertEqual(res._sha256("/home/test/file"),"c"*64)
            self.assertGreaterEqual(call.call_args.kwargs["timeout"],1800)

    def test_atomic_download_audit_updates(self):
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder)
            res._write_gnomad_audit(out,{"active_job":{"phase":"HASHING_SHA256"},
                                         "completion":"NOT_STARTED","records":[]})
            self.assertIn('"HASHING_SHA256"',(out/"GNOMAD_AUTO_RESOURCE_AUDIT.json").read_text())
            self.assertFalse((out/"GNOMAD_AUTO_RESOURCE_AUDIT.json.partial").exists())
            res._write_gnomad_audit(out,{"active_job":{"phase":"VERIFIED"},
                                         "completion":"READY","records":[]})
            self.assertIn('"READY"',(out/"GNOMAD_AUTO_RESOURCE_AUDIT.json").read_text())

class VepDependencyTests(unittest.TestCase):
    def test_complete_required_zip_extractor_dependencies(self):
        text=(ROOT/"code/v2810/resource_bootstrap.py").read_text(encoding="utf-8")
        for pkg in ("unzip","tar","gzip","libarchive-extract-perl"):
            self.assertIn('"'+pkg+'"',text)
        self.assertIn('"-MArchive::Extract"',text)
        self.assertIn('"-MBio::EnsEMBL::Registry"',text)
        self.assertIn('toolroot+"/ensembl-vep-116"',text)
        vep_source=(ROOT/"code/v2810/vep_local.py").read_text(encoding="utf-8")
        self.assertIn('str(PurePosixPath(vep).parent)',vep_source)
        self.assertIn('args+=["--offline","--cache","--dir_cache",cache_dir]',vep_source)
        self.assertIn('args+=["--gtf",gtf]',vep_source)
        self.assertIn('.rdmva_api116_ready',text)
        self.assertIn('"--DESTDIR"',text)
        self.assertIn('V31_REQUIRED_VEP_SETUP_INCOMPLETE',
                      (ROOT/"code/v30/native_pipeline.py").read_text(encoding="utf-8"))

    def test_incomplete_gnomad_does_not_prove_full_profile(self):
        s=(ROOT/"code/v30/native_pipeline.py").read_text(encoding="utf-8")
        self.assertIn('V31_REQUIRED_GNOMAD_INCOMPLETE',s)
        self.assertIn('gnomad_audit.get("completion")!="READY"',s)

if __name__=="__main__":unittest.main()
