"""Windows jury input auto-discovery: only unambiguous authorized local files.

Fixtures are synthetic and deliberately contain no patient information.
"""
from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.name=="nt","Windows PowerShell launcher integration")
class JuryAutomaticInputTests(unittest.TestCase):
    def test_unique_separate_folders_identical_copies_and_ambiguity(self):
        source=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        start=source.index("function Discover-AuthorizedInputs {")
        end=source.index("function Choose-File(",start)
        function=source[start:end]
        with tempfile.TemporaryDirectory(prefix="RDMVA-AUTO-") as temp:
            base=Path(temp)
            ps=base/"test-discovery.ps1"
            for mode in ("adjacent","separate","duplicates","different"):
                for existing in base.iterdir():
                    if existing.is_dir():
                        import shutil
                        shutil.rmtree(existing)
                pair1=base/"CaseA"
                pair1.mkdir()
                if mode=="separate":
                    doc_dir=base/"CaseB"
                    doc_dir.mkdir()
                else:
                    doc_dir=pair1
                (pair1/"Synthetic.vcf.gz").write_bytes(b"simulated-vcf-1")
                (doc_dir/"Challenge_Clinical_Phenotype.docx").write_bytes(b"simulated-clinical-1")
                if mode in ("duplicates","different"):
                    pair2=base/"CaseC"
                    pair2.mkdir()
                    (pair2/"Synthetic.vcf.gz").write_bytes(
                        b"simulated-vcf-2" if mode=="different" else b"simulated-vcf-1")
                    (pair2/"Challenge_Clinical_Phenotype.docx").write_bytes(
                        b"simulated-clinical-2" if mode=="different" else b"simulated-clinical-1")
                location=str(base).replace("'","''")
                script=(
                    "$ErrorActionPreference='Stop'\n"
                    "$Root='"+location+"'\n"
                    "$env:USERPROFILE='"+location+"'\n"
                    "$env:RDMVA_CHALLENGE_INPUT_DIR='"+location+"'\n"
                    "function Get-PSDrive { @() }\n"
                    "function Read-Host { param([string]$Prompt); return 'YES' }\n"
                    "$script:VcfPath=''; $script:ClinicalPath=''\n"
                    +function+
                    "\n$r=Discover-AuthorizedInputs\n"
                    "Write-Output ('AUTO_INPUT_RESULT='+[string]$r)\n"
                    "Write-Output ('AUTO_VCF_PRESENT='+[string][bool]$script:VcfPath)\n"
                    "Write-Output ('AUTO_DOC_PRESENT='+[string][bool]$script:ClinicalPath)\n"
                )
                ps.write_text(script,encoding="utf-8")
                process=subprocess.run(["powershell.exe","-NoProfile","-NonInteractive",
                    "-ExecutionPolicy","Bypass","-File",str(ps)],
                    capture_output=True,text=True,timeout=35,errors="replace")
                with self.subTest(case=mode):
                    self.assertEqual(process.returncode,0,process.stdout+process.stderr)
                    expected=mode!="different"
                    self.assertIn("AUTO_INPUT_RESULT="+str(expected),process.stdout)
                    self.assertIn("AUTO_VCF_PRESENT="+str(expected),process.stdout)
                    self.assertIn("AUTO_DOC_PRESENT="+str(expected),process.stdout)
                    if mode=="duplicates":
                        self.assertIn("RDMVA_DUPLICATE_AUTHORIZED_PAIRS_DEDUPLICATED",process.stdout)
                    if mode=="separate":
                        self.assertIn("RDMVA_UNIQUE_AUTHORIZED_SEPARATE_FOLDER_PAIR_FOUND",process.stdout)


if __name__=="__main__":
    unittest.main()
