"""Windows jury preflight contracts. Synthetic tests only; no credentials or installs."""
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch
import os,tempfile
from v2810 import resource_bootstrap as rb
from v2810 import ai
ROOT=Path(__file__).resolve().parents[1]

class WindowsInstallCompatibilityTests(TestCase):
    def test_installer_handles_old_git_and_python(self):
        text=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        for token in (
            "WingetEnsure 'Git.Git' 'upgrade'",
            "WingetEnsure 'Python.Python.3.12' 'upgrade'",
            "WingetEnsure 'Git.Git' 'install'",
            "GitPatch","PythonPatch",
            "$GitMinVersion=[version]'2.55.0'",
            "$PythonMinVersion=[version]'3.12.10'",
        ):
            self.assertIn(token,text)
        self.assertNotIn("if($version -notmatch '^git version 2",text)

    def test_windows_required_reboot_is_not_a_generic_install_failure(self):
        launcher=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        installer=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("if ($status -eq 10)",launcher)
        self.assertIn("RDMVA_JURY_INSTALL_ACTION_REQUIRED=",launcher)
        self.assertIn("if ($Action -eq 'install') {exit 10}",launcher)
        self.assertIn("restart Windows when safe",installer)
        self.assertIn("Record 'REBOOT_REQUIRED'",installer)

    def test_new_ubuntu_distribution_stored_beside_project_not_on_C(self):
        text=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("Join-Path $Root 'env\\wsl_ubuntu_24_04'",text)
        self.assertIn("--location $distroLocation",text)
        self.assertIn("$UbuntuSource='Ubuntu-24.04'",text)
        self.assertIn("'RDMVA-Ubuntu-24.04-'+$distroSuffix",text)
        self.assertIn("--name $Ubuntu",text)
        self.assertIn("RDMVA_ISOLATED_WSL_DISTRIBUTION=",text)

    def test_two_independent_checkouts_have_separate_research_distros(self):
        a=rb.project_wsl_distribution("T:/first_challenge_checkout")
        b=rb.project_wsl_distribution("T:/second_challenge_checkout")
        self.assertNotEqual(a,b)
        self.assertTrue(a.startswith("RDMVA-Ubuntu-24.04-"))
        self.assertTrue(b.startswith("RDMVA-Ubuntu-24.04-"))
        self.assertNotEqual(a,"Ubuntu-24.04")

    @__import__("unittest").skipUnless(os.name=="nt","Windows root hash agreement")
    def test_windows_powershell_and_python_derive_same_distro_name(self):
        import subprocess
        src=ROOT.resolve()
        formula=(
          "$p=[IO.Path]::GetFullPath('"+str(src).replace("'","''")+"').TrimEnd('\\').ToLowerInvariant();"
          "$h=[Security.Cryptography.SHA256]::Create();"
          "try{$b=$h.ComputeHash([Text.Encoding]::UTF8.GetBytes($p));"
          "$s=([BitConverter]::ToString($b)).Replace('-','').ToLowerInvariant().Substring(0,8);"
          "Write-Output ('RDMVA-Ubuntu-24.04-'+$s)}finally{$h.Dispose()}"
        )
        observed=subprocess.check_output(["powershell.exe","-NoProfile","-Command",formula],
                                         text=True,timeout=20).strip()
        self.assertEqual(observed,rb.project_wsl_distribution(src))

    def test_linux_user_provisioning_noninteractive_and_nonprivileged(self):
        source=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        for required in (
            "--no-launch -d $UbuntuSource --name $Ubuntu --location $distroLocation",
            "$researchUser='rdmva'",
            "wsl.exe -d $Ubuntu --user root --exec bash -c",
            "useradd -m -s /bin/bash rdmva",
            "wsl.exe --manage $Ubuntu --set-default-user $researchUser",
            "wsl.exe --terminate $Ubuntu",
            "$actualUser -ne $researchUser",
            "$actualHome -ne ('/home/'+$researchUser)",
            "RDMVA_WSL_UNPRIVILEGED_ACCOUNT_READY=",
        ):
            self.assertIn(required,source)
        self.assertNotIn("Ubuntu may request its standard local Linux user setup once",source)
        quickstart=(ROOT/"README.md").read_text(encoding="utf-8")
        self.assertIn("The reviewer is not asked for a Linux",quickstart)
        self.assertIn("no Linux username/password prompt or shell exit is required",
                      (ROOT/"README.md").read_text(encoding="utf-8"))

    def test_interrupted_first_run_resumes_without_reinstalling_ubuntu(self):
        source=(ROOT/"INSTALL_WINDOWS_V31.ps1").read_text(encoding="utf-8")
        for required in (
            "function UbuntuRegistration",
            "function UbuntuOwnedByProject",
            "function InitializeResearchUser",
            "RDMVA_WSL_RESUMING_INTERRUPTED_NONINTERACTIVE_SETUP",
            "if(UbuntuRegistered)",
            "if(!(UbuntuOwnedByProject))",
            "if(!(InitializeResearchUser)){exit 10}",
            "if($actual.StartsWith('",
        ):
            self.assertIn(required,source)
        self.assertEqual(source.count("wsl.exe --install --web-download --no-launch"),1)
        self.assertNotIn("wsl.exe --unregister",source)

    def test_wsl_scientific_execution_never_uses_windows_default(self):
        name=rb.project_wsl_distribution()
        self.assertTrue(name.startswith("RDMVA-Ubuntu-24.04-"))
        self.assertEqual(len(name.split("-")[-1]),8)
        self.assertEqual(rb.wsl_command("tabix","file","chr1:1-1"),
            ["wsl.exe","-d",name,"--exec","tabix","file","chr1:1-1"])
        self.assertEqual(rb.wsl_command("apt-get","update",root=True),
            ["wsl.exe","-d",name,"--user","root","--exec","apt-get","update"])
        self.assertEqual(rb._default_wsl_distribution(),name)
        for name in ("gnomad_local","open_splicing","spliceai_local","vep_local"):
            s=(ROOT/"code/v2810"/(name+".py")).read_text(encoding="utf-8")
            self.assertNotIn('["wsl.exe","--exec",',s)
            self.assertIn("wsl_command",s)

class NativeKeyDiscoveryTests(TestCase):
    def test_new_reviewer_without_key_has_no_inferred_secret(self):
        # Simulates a new Windows user with NO native key (without accessing
        # the real user's credential or network).
        with tempfile.TemporaryDirectory() as td,\
             patch.dict(os.environ,{"OPENAI_API_KEY":""}),\
             patch.object(ai,"STATE",Path(td)):
            with self.assertRaises(Exception):
                ai.api_key()

    def test_reviewers_may_use_their_own_environment_key(self):
        fake="TEST_ONLY_NOT_A_REAL_CREDENTIAL"
        with tempfile.TemporaryDirectory() as td,\
             patch.dict(os.environ,{"OPENAI_API_KEY":fake}),\
             patch.object(ai,"STATE",Path(td)):
            self.assertEqual(ai.api_key(),fake)
        launcher=(ROOT/"JURY_START_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("$env:LOCALAPPDATA",launcher)
        self.assertIn("Configure-Key",launcher)
        self.assertNotIn("C:\\Users\\matte",launcher)
        secure=(ROOT/"SETUP_OPENAI_KEY_V31.ps1").read_text(encoding="utf-8")
        self.assertIn("Read-Host",secure)
        self.assertIn("-AsSecureString",secure)
        self.assertIn("DataProtectionScope]::CurrentUser",secure)
