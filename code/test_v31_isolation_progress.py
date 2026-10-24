"""V3.1.1 clean-PC isolation and truthful long-task progress regression tests.

Synthetic only: does not install Ubuntu, download data, inspect patients or call AI.
"""
from __future__ import annotations
import contextlib,hashlib,importlib.util,io,subprocess,sys,tempfile,time,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from v2810 import resource_bootstrap as rb
from v2810 import chembl_sqlite_download_v31 as ch37

def load_setup():
    spec=importlib.util.spec_from_file_location("rdmva_v311_isolated_setup",ROOT/"setup_v2810.py")
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

class IsolationAndProgressTests(unittest.TestCase):
    def test_setup_probe_uses_only_path_derived_distro(self):
        setup=load_setup()
        fake=SimpleNamespace(returncode=0,
            stdout="RDMVA:ubuntu:24.04:rdmva:/home/rdmva",stderr="")
        with patch.object(setup.subprocess,"run",return_value=fake) as run:
            result=setup.probe_isolated_wsl()
        argv=run.call_args.args[0]
        self.assertEqual(argv[:3],["wsl.exe","-d",rb.project_wsl_distribution(ROOT)])
        self.assertEqual(argv[3:5],["--exec","bash"])
        self.assertEqual(result["user"],"rdmva")
        self.assertEqual(result["distribution"],rb.project_wsl_distribution(ROOT))
        self.assertNotEqual(result["distribution"],"Ubuntu")
    def test_setup_rejects_unexpected_wsl_release_or_default_user(self):
        setup=load_setup()
        for reported in ("RDMVA:ubuntu:26.04:rdmva:/home/rdmva",
                         "RDMVA:ubuntu:24.04:root:/root",
                         "RDMVA:debian:24.04:rdmva:/home/rdmva"):
            fake=SimpleNamespace(returncode=0,stdout=reported,stderr="")
            with patch.object(setup.subprocess,"run",return_value=fake):
                with self.assertRaisesRegex(RuntimeError,"ISOLATED_UBUNTU"):
                    setup.probe_isolated_wsl()
    def test_installer_checks_registry_ownership_and_image(self):
        src=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("if(!(UbuntuRegistered) -or !(UbuntuOwnedByProject))",src)
        self.assertIn("$release -ne 'ubuntu:24.04'",src)
        self.assertIn("$arch -ne 'x86_64'",src)
        self.assertIn("$name -ne 'rdmva'",src)
        self.assertNotIn("wsl.exe --unregister",src)
    def test_no_implicit_default_wsl_probe(self):
        src=(ROOT/"setup_v2810.py").read_text(encoding="utf-8")
        self.assertNotIn('["wsl.exe","--exec"',src)
        self.assertIn("probe_isolated_wsl()",src)
        self.assertIn("wsl_command(",src)
    def test_quick_check_heartbeat_does_not_invent_percent(self):
        class FakeDb:
            def __init__(self):self.handler=None;self.end_handlers=0
            def set_progress_handler(self,handler,steps):
                self.handler=handler
                if handler is None:self.end_handlers+=1
            def execute(self,sql):
                if sql!="PRAGMA quick_check(1)":raise AssertionError(sql)
                time.sleep(.055)
                if self.handler:self.handler()
                return self
            def fetchone(self):return ("ok",)
        db=FakeDb()
        out=io.StringIO()
        with contextlib.redirect_stdout(out):
            answer=ch37._sqlite_quick_check_with_activity(db,interval_s=.01)
        self.assertEqual(answer,("ok",))
        self.assertEqual(db.end_handlers,1)
        self.assertIn("V31_CH37_SQLITE_QUICK_CHECK_ACTIVITY",out.getvalue())
        self.assertIn("percent=unavailable",out.getvalue())
    def test_true_hash_still_computes_real_digest(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/"synthetic.txt"
            p.write_bytes(b"synthetic-only integrity test")
            self.assertEqual(rb._sha256_native_file(p,"TEST"),
                             hashlib.sha256(p.read_bytes()).hexdigest())
            # Scientific bootstrap receives WSL /mnt/<drive>/ paths,
            # not native C:\\... paths. Reproduce that real call contract.
            linux_path="/mnt/"+p.drive[0].lower()+"/"+p.as_posix()[3:]
            self.assertEqual(rb._sha256(linux_path),
                             hashlib.sha256(p.read_bytes()).hexdigest())
    def test_long_stages_emit_status_without_false_success(self):
        source=(ROOT/"reproduce_v31_pipeline.py").read_text(encoding="utf-8")
        self.assertIn('long_stage_activity("stage1_native_core")',source)
        self.assertIn('long_stage_activity("stage2_v31_research")',source)
        self.assertIn("percent=unavailable",source)
        self.assertIn("E2E_V31_REPRODUCTION_PASS",source)
    def test_downloads_and_reference_hashes_report_real_progress(self):
        sources=[ROOT/"code/v2810/chembl_sqlite_download_v31.py",
                 ROOT/"code/v2810/resource_bootstrap.py",
                 ROOT/"bootstrap_current_public_v31.py"]
        all_sources="\n".join(f.read_text(encoding="utf-8") for f in sources)
        for marker in ("V31_CH37_SHA256_PROGRESS",
                       "V31_CH37_DB_EXTRACTION_PROGRESS",
                       "V31_PUBLIC_SHA256_PROGRESS",
                       "V31_PUBLIC_DOWNLOAD_PROGRESS",
                       "V31_CURRENT_SHA256_PROGRESS",
                       "V31_CURRENT_DOWNLOAD_ACTIVITY"):
            self.assertIn(marker,all_sources)
    def test_isolated_python_venv_and_hash_pins(self):
        source=(ROOT/"setup_v2810.py").read_text(encoding="utf-8")
        self.assertIn('ENV=ROOT/"env/v28_runtime"',source)
        self.assertIn('"--require-hashes"',source)
        self.assertIn('"--only-binary=:all:"',source)
        self.assertIn('version_info[:2]!=(3,12)',source)
        self.assertIn('V31_SETUP_STAGE_ACTIVITY',source)

if __name__=="__main__":unittest.main()
