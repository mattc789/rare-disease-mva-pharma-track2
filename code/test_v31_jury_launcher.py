"""Jury Windows interface and release-readiness tests (synthetic only)."""
from pathlib import Path
import json,os,shutil,subprocess,sys,tempfile,unittest
import importlib.util

ROOT=Path(__file__).resolve().parents[1]

class WindowsJuryLaunchTests(unittest.TestCase):
    def test_launchers_and_docs_are_in_source_only_release_manifest(self):
        spec=importlib.util.spec_from_file_location(
            "build_v31_jury_source_release",ROOT/"build_source_release_v31.py")
        builder=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(builder)
        required=[
          "JURY_START_V31.cmd","JURY_START_V31.ps1",
          "JURY_MONITOR_V31.py","README.md",
          "code/test_v31_jury_launcher.py"]
        for file in required:
            self.assertTrue((ROOT/file).is_file(),file)
            self.assertIn(file,builder.FILES)
        self.assertNotIn("JURY_MONITOR_V31.ps1",builder.FILES)
        self.assertIn(".cmd",builder.TEXT_SUFFIXES)
        self.assertIn("V31_API_CREDENTIAL_ALREADY_PRESENT=YES",
          (ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8"))
        self.assertIn(b"\r\n",builder.canonical_bytes(
            "JURY_START_V31.cmd",ROOT/"JURY_START_V31.cmd"))

    @unittest.skipUnless(os.name=="nt","PowerShell parser on Windows")
    def test_windows_launcher_parses_without_any_side_effect(self):
        ps=shutil.which("powershell.exe")
        for name in ("JURY_START_V31.ps1",):
            script=ROOT/name
            command=(
             "$tokens=$null;$errors=$null;"
             "[System.Management.Automation.Language.Parser]::ParseFile("
             "'"+str(script).replace("'","''")+"',"
             "[ref]$tokens,[ref]$errors)|Out-Null;"
             "if(@($errors).Count){$errors|ForEach-Object{$_.Message};exit 7};"
             "Write-Output 'PARSE_OK'")
            p=subprocess.run([ps,"-NoProfile","-Command",command],
                capture_output=True,text=True,timeout=15)
            self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn("PARSE_OK",p.stdout)

    def test_monitor_observes_synthetic_run_once_without_patient_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            run=Path(tmp)
            files={
              "E2E_EXECUTION_V31.json":{"status":"RUNNING",
                  "online_ai_requested":True,"full_scientific_modules":True},
              "STATUS_V30_NATIVE.json":{"current_stage":"GENOMIC_EVIDENCE",
                   "stages":{"GENOMIC_EVIDENCE":{
                       "utc":"2026-10-20T10:00:00Z","state":"running"}}},
              "resource_bootstrap/gnomad/GNOMAD_AUTO_RESOURCE_AUDIT.json":{
                  "completion":"NOT_STARTED","active_job":{
                      "cohort":"exomes","chrom":"7","phase":"DOWNLOADING"}},
            }
            for name,obj in files.items():
                p=run/name
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_text(json.dumps(obj),encoding="utf-8")
            p=subprocess.run([sys.executable,str(ROOT/"JURY_MONITOR_V31.py"),
                   "--out",str(run),"--interval","0"],
                   capture_output=True,text=True,timeout=20)
            self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn("NATIVE_STAGE=GENOMIC_EVIDENCE",p.stdout)
            self.assertIn("GNOMAD_ACTIVE_JOB=exomes/chr7/DOWNLOADING",p.stdout)
            self.assertIn("RUN_STATUS=RUNNING",p.stdout)

    def test_monitor_never_inspects_processes_or_runs_shell_commands(self):
        viewer=(ROOT/"JURY_MONITOR_V31.py").read_text(encoding="utf-8")
        launcher=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        builder=(ROOT/"build_source_release_v31.py").read_text(encoding="utf-8")
        self.assertIn("MAX_STATUS_BYTES",viewer)
        self.assertIn("STATUS_FILES",viewer)
        self.assertNotIn("subprocess",viewer)
        self.assertNotIn("Get-CimInstance",viewer)
        self.assertNotIn(".rglob(",viewer)
        self.assertNotIn("os.system(",viewer)
        self.assertIn("Start-Process -FilePath $Runtime",launcher)
        self.assertNotIn("Start-Process -FilePath 'powershell.exe'",launcher)
        self.assertNotIn("JURY_MONITOR_V31.ps1",launcher)
        self.assertNotIn("JURY_MONITOR_V31.ps1",builder)

    def test_monitor_never_displays_patient_payload(self):
        with tempfile.TemporaryDirectory() as temp:
            run=Path(temp)
            private_marker="TEST_ONLY_PRIVATE_GENOME_NEVER_SHOW_THIS"
            (run/"case_private.json").write_text(
                json.dumps({"private_variant":private_marker}),encoding="utf-8")
            (run/"E2E_EXECUTION_V31.json").write_text(
                json.dumps({"status":"RUNNING","online_ai_requested":False}),
                encoding="utf-8")
            p=subprocess.run([sys.executable,str(ROOT/"JURY_MONITOR_V31.py"),
                "--out",str(run),"--interval","0"],
                capture_output=True,text=True,timeout=20)
            self.assertEqual(p.returncode,0,p.stdout+p.stderr)
            self.assertIn("RUN_STATUS=RUNNING",p.stdout)
            self.assertNotIn(private_marker,p.stdout)
            self.assertNotIn("case_private",p.stdout)

    def test_jury_workflow_order_and_case_insensitive_consent(self):
        source=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("Write-Host ' 3  Test live OpenAI models",source)
        self.assertIn("Write-Host ' 4  Run full Challenge analysis'",source)
        self.assertIn("'3' {Test-LiveModels; break}",source)
        # Option 4 must launch the run AND clear stale selections even after errors.
        run_action=next((line for line in source.splitlines() if "'4' {" in line),"")
        self.assertIn("try { Start-Research } finally {",run_action)
        for variable in ("Out","VcfPath","ClinicalPath"):
            self.assertIn(f"$script:{variable}=''",run_action)
        self.assertNotIn("-cne 'YES'",source)
        self.assertEqual(source.count(".Trim() -ine 'YES'"),2)
        self.assertIn(".Trim() -ieq 'YES'",source)

    def test_authorized_inputs_discovered_by_filename_only_with_confirmation(self):
        source=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("function Discover-AuthorizedInputs",source)
        self.assertIn("RDMVA_CHALLENGE_INPUT_DIR",source)
        self.assertIn("RDMVA_AUTHORIZED_INPUT_PAIR_CONFIRMED",source)
        self.assertIn("Multiple",source.replace("MULTIPLE","Multiple"))
        self.assertIn("if (!$VcfPath -and !$ClinicalPath)",source)
        self.assertIn("$null=Discover-AuthorizedInputs",source)
        self.assertIn("Get-ChildItem -LiteralPath $directory -Filter '*.vcf.gz'",source)

    def test_long_public_download_progress_is_source_only(self):
        source=(ROOT/"code/v2810/resource_bootstrap.py").read_text(encoding="utf-8")
        self.assertIn("V31_LARGE_PUBLIC_DOWNLOAD_NOTICE",source)
        self.assertIn("V31_PUBLIC_DOWNLOAD_PROGRESS",source)
        self.assertIn("current/total",source)
        self.assertIn("--continue-at",source)
        chembl=(ROOT/"code/v2810/chembl_sqlite_download_v31.py").read_text(encoding="utf-8")
        self.assertIn("V31_CH37_DOWNLOAD_NOTICE",chembl)

    def test_jury_launcher_never_contains_authorized_case_file_identifiers(self):
        for name in ("JURY_START_V31.ps1","JURY_MONITOR_V31.py",
                     "JURY_START_V31.cmd","README.md"):
            text=(ROOT/name).read_text(encoding="utf-8")
            self.assertNotRegex(text,r"(?i)WGS[_-][A-Z0-9]{8,}")
            self.assertNotRegex(text,r"(?i)clinical_phenotype_\d+")
            self.assertNotRegex(text,r"\b\d{1,3}(?:\.\d{1,3}){3}\b")
            self.assertNotIn("BEGIN PRIVATE KEY",text)

if __name__=="__main__":
    unittest.main()
