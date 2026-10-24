"""Non-destructive storage auto-configuration and installer integration tests."""
from pathlib import Path
import re,subprocess,unittest,os

ROOT=Path(__file__).resolve().parents[1]
SELECTOR=ROOT/"CONFIGURE_STORAGE_V31.ps1"
INSTALLER=ROOT/"INSTALL_WINDOWS_V31.ps1"
RELEASE=ROOT/"build_source_release_v31.py"

class StorageAutomationTests(unittest.TestCase):
    def test_installer_invokes_real_storage_preflight(self):
        self.assertTrue(SELECTOR.is_file())
        install=INSTALLER.read_text(encoding="utf-8")
        self.assertIn("CONFIGURE_STORAGE_V31.ps1",install)
        self.assertIn("RDMVA_STORAGE_AUTO_SELECTION_START",install)
        self.assertIn("RDMVA_SCIENTIFIC_STORAGE_WSL_READY",install)
        self.assertIn("'-DataDrive'",install)
        self.assertLess(install.index("RDMVA_STORAGE_AUTO_SELECTION_START"),
                        install.index("$version=InstalledVersion 'py.exe'"))
        self.assertIn("SINGLE_PROJECT_ROOT",SELECTOR.read_text(encoding="utf-8"))
        self.assertIn("Join-Path $Root 'public_inputs\\third_party'",
                      SELECTOR.read_text(encoding="utf-8"))
        release=RELEASE.read_text(encoding="utf-8")
        self.assertIn('"CONFIGURE_STORAGE_V31.ps1"',release)
        self.assertIn('"SETUP_OPENAI_KEY_V31.ps1"',release)

    def test_no_storage_initialization_or_formatting(self):
        s=SELECTOR.read_text(encoding="utf-8")
        for unsafe in ("Initialize-Disk","Clear-Disk","New-Partition",
                       "Format-Volume","New-VHD","Remove-Partition",
                       "Set-Disk","Mount-VHD","Remove-Item"):
            self.assertNotRegex(s,r"(?im)^\s*(?:&\s*)?"+re.escape(unsafe)+r"\b")
        self.assertIn("340",s)
        self.assertIn("Get-Volume",s)
        self.assertIn("Get-Partition",s)
        self.assertIn("expected_disk_name",s)
        self.assertIn("V31_STORAGE_EXISTING_CONFIG_VERIFIED",s)

    @unittest.skipUnless(os.name=="nt","PowerShell parser only available on Windows")
    def test_windows_powershell_parser(self):
        for path in (SELECTOR,INSTALLER,
                     ROOT/"REPRODUCE_CHALLENGE_V31.ps1",
                     ROOT/"SETUP_OPENAI_KEY_V31.ps1"):
            literal=str(path).replace("'","''")
            script=("$e=$null;$t=$null;[Management.Automation.Language.Parser]::ParseFile('"
                    +literal+"',[ref]$t,[ref]$e)|Out-Null;"
                    +"if(@($e).Count){$e|ForEach-Object{$_.Message};exit 3}")
            proc=subprocess.run(["powershell.exe","-NoProfile","-NonInteractive",
                "-Command",script],capture_output=True,text=True,timeout=30)
            self.assertEqual(proc.returncode,0,(path.name,proc.stdout,proc.stderr))

class CredentialContractTests(unittest.TestCase):
    def test_credential_is_automatically_derived_and_no_plaintext_env_file_required(self):
        launcher=(ROOT/"REPRODUCE_CHALLENGE_V31.ps1").read_text(encoding="utf-8")
        secret_helper=(ROOT/"SETUP_OPENAI_KEY_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("openai_key.dpapi",launcher)
        self.assertIn("LOCALAPPDATA",launcher)
        self.assertIn("ProtectedData",secret_helper)
        self.assertIn("Read-Host",secret_helper)
        self.assertIn("-AsSecureString",secret_helper)
        self.assertIn("CurrentUser",secret_helper)
        self.assertNotIn("Set-Content",secret_helper)
        self.assertNotIn("Add-Content",secret_helper)
        doc=(ROOT/"README.md").read_text(encoding="utf-8")
        self.assertIn("does not automatically read a .env file",doc)
        self.assertIn("SETUP_OPENAI_KEY_V31.ps1",doc)

if __name__=="__main__":unittest.main()
