"""No privileged PowerShell automatic variables may be assigned by the
reviewer installer. Verifies the existing Ubuntu distro recovery branch with
synthetic commands; never calls real WSL or downloads any patient data.
"""
from __future__ import annotations
import os,re,subprocess,tempfile,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
INSTALLER=ROOT/"INSTALL_WINDOWS_V31.ps1"

class PowerShellReservedVariablesTests(unittest.TestCase):
    def test_readonly_powershell_automatic_variables_never_assigned(self):
        source=INSTALLER.read_text(encoding="utf-8")
        # Matching is deliberately case-insensitive: $home == $HOME in
        # Windows PowerShell, and automatic variable HOME is ReadOnly.
        forbidden=("home","host","pid","pshome","false","true","null")
        for name in forbidden:
            assignment=r"(?im)^\s*\$"+name+r"\s*(?:=|\+=|-=|\+\+|--)"
            self.assertNotRegex(source,assignment,
                                "Do not overwrite the PowerShell automatic variable: "+name)
        self.assertIn("$linuxHomePath=(",source)
        self.assertIn("$linuxHomePath -ne '/home/rdmva'",source)
        self.assertNotRegex(source,r"(?im)^\s*\$profile\s*=")

    @unittest.skipUnless(os.name=="nt","The Windows WSL readiness branch is PowerShell-specific")
    def test_readonly_wsl_probe_accepts_existing_distro_without_reinstall(self):
        source=INSTALLER.read_text(encoding="utf-8")
        start=source.index("function UbuntuReady {")
        end=source.index("function InitializeResearchUser {",start)
        function=source[start:end]
        # A fake per-installation WSL-owned folder forces the branch that used
        # to assign to $home and throw. Every WSL response is mocked; no call
        # reaches real Windows Linux networking, installers or credentials.
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp)/"env/wsl_ubuntu_24_04").mkdir(parents=True)
            escaped=tmp.replace("'","''")
            powershell="\n".join([
                "$ErrorActionPreference='Stop'",
                "$Root='"+escaped+"'",
                "$Ubuntu='RDMVA-Ubuntu-24.04-SYNTHETIC'",
                "function FeatureOn([string]$feature){return $true}",
                "function UbuntuRegistered {return $true}",
                "function UbuntuOwnedByProject {return $true}",
                "function wsl.exe {",
                "  $called=($args -join ' ')",
                "  if($called -match 'os-release') { $global:LASTEXITCODE=0;Write-Output 'ubuntu:24.04';return }",
                "  if($called -match 'uname -m') { $global:LASTEXITCODE=0;Write-Output 'x86_64';return }",
                "  if($called -match 'id -un') { $global:LASTEXITCODE=0;Write-Output 'rdmva';return }",
                "  if($called -match 'printf %s') { $global:LASTEXITCODE=0;Write-Output '/home/rdmva';return }",
                "  throw ('FAKE_WSL_UNEXPECTED_COMMAND '+$called)",
                "}",
                function,
                "if(!(UbuntuReady)){throw 'SYNTHETIC_EXISTING_WSL_SHOULD_BE_READY'}",
                "Write-Output 'JURY_EXISTING_ISOLATED_UBUNTU_READY_WITHOUT_READONLY_VARIABLE_ASSIGNMENT'",
            ])
            proc=subprocess.run(["powershell.exe","-NoProfile","-NonInteractive",
                "-ExecutionPolicy","Bypass","-Command",powershell],
                text=True,errors="replace",capture_output=True,timeout=25)
            self.assertEqual(proc.returncode,0,proc.stderr+"\n"+proc.stdout)
            self.assertIn("JURY_EXISTING_ISOLATED_UBUNTU_READY",proc.stdout)

if __name__=="__main__":
    unittest.main()
