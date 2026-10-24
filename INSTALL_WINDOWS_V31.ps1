# Rare Disease MVA Pharma V3.1: clean Windows 11 setup from GitHub source.
# Installs pinned Windows prerequisites, then the hash-locked scientific runtime.
param([switch]$PrerequisitesOnly,[switch]$SkipReferenceDownload,[switch]$CurrentPublicData,[ValidatePattern('^[A-Za-z]$')][string]$DataDrive)
$ErrorActionPreference='Stop'
$Root=Split-Path -Parent $MyInvocation.MyCommand.Path
$PythonVersion='3.12.10'
$GitVersion='2.55.0'
$PythonMinVersion=[version]'3.12.10'
$GitMinVersion=[version]'2.55.0'
# The Ubuntu 24.04 distribution is isolated PER source installation.
# Its name is derived deterministically from the absolute project path;
# other tools' Ubuntu distros and older RDMVA installations are untouched.
$UbuntuSource='Ubuntu-24.04'
$CanonicalProjectPath=[IO.Path]::GetFullPath($Root).TrimEnd('\').ToLowerInvariant()
$hashAlgorithm=[Security.Cryptography.SHA256]::Create()
try {
  $hashBytes=$hashAlgorithm.ComputeHash([Text.Encoding]::UTF8.GetBytes($CanonicalProjectPath))
  $distroSuffix=([BitConverter]::ToString($hashBytes)).Replace('-','').ToLowerInvariant().Substring(0,8)
} finally {$hashAlgorithm.Dispose()}
$Ubuntu='RDMVA-Ubuntu-24.04-'+$distroSuffix
$Status=Join-Path $Root 'INSTALL_WINDOWS_V31_STATUS.json'

function Record([string]$state,[string]$message) {
  $obj=[ordered]@{schema='RDMVA_WINDOWS_INSTALL_V31_1';state=$state;message=$message;
    python_pin=$PythonVersion;git_pin=$GitVersion;ubuntu_pin=$Ubuntu;
    utc=[datetime]::UtcNow.ToString('o');contains_secrets=$false}
  $obj|ConvertTo-Json|Set-Content -LiteralPath $Status -Encoding UTF8
  Write-Host ('RDMVA_INSTALL_'+$state+': '+$message)
}
function Abort([string]$msg) { throw $msg }
function InstalledVersion([string]$cmd,[string[]]$arguments) {
  if(!(Get-Command $cmd -ErrorAction SilentlyContinue)){return ''}
  $value=(& $cmd @arguments 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0){return ''}
  return $value
}
function WingetEnsure([string]$pkg,[ValidateSet('install','upgrade')][string]$mode) {
  # WinGet's package patch version can change (e.g., Git.Git 2.55.0.5).
  # Select its current official release, then verify compatible minimums.
  Write-Host ('WINDOWS_DEPENDENCY_ACTION='+$mode+' PACKAGE='+$pkg)
  $args=@($mode,'--id',$pkg,'--exact','--source','winget',
    '--accept-source-agreements','--accept-package-agreements',
    '--disable-interactivity','--silent')
  & winget.exe @args
  if($LASTEXITCODE -ne 0){Abort ('V31_WINGET_'+$mode+'_FAILED:'+ $pkg+':'+$LASTEXITCODE)}
  $env:Path=$env:Path+';'+[Environment]::GetEnvironmentVariable('Path','Machine')+';'+[Environment]::GetEnvironmentVariable('Path','User')
}
function PythonPatch([string]$actual) {
  if($actual -notmatch '^Python (3\.12\.\d+)$') {
    Abort ('V31_PYTHON_3_12_REQUIRED: '+$actual)
  }
  return [version]$Matches[1]
}
function GitPatch([string]$actual) {
  if($actual -notmatch '^git version (\d+\.\d+\.\d+)(?:\.|$)') {
    Abort ('V31_GIT_VERSION_UNRECOGNIZED: '+$actual)
  }
  return [version]$Matches[1]
}
function FeatureOn([string]$name) {
  $f=Get-WindowsOptionalFeature -Online -FeatureName $name -ErrorAction Stop
  return ($f.State -eq 'Enabled')
}
function UbuntuRegistration {
  @(Get-ChildItem 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Lxss' -ErrorAction SilentlyContinue |
      Get-ItemProperty -ErrorAction SilentlyContinue |
      Where-Object {$_.DistributionName -eq $Ubuntu -and $_.Version -eq 2})
}
function UbuntuRegistered {
  return (@(UbuntuRegistration).Count -eq 1)
}
function UbuntuOwnedByProject {
  $keys=@(UbuntuRegistration)
  if($keys.Count -ne 1){return $false}
  $actual=[string]$keys[0].BasePath
  if([string]::IsNullOrWhiteSpace($actual)){return $false}
  if($actual.StartsWith('\\?\')){$actual=$actual.Substring(4)}
  $expected=Join-Path $Root 'env\wsl_ubuntu_24_04'
  try {
    return ([IO.Path]::GetFullPath($actual).TrimEnd('\') -ieq
            [IO.Path]::GetFullPath($expected).TrimEnd('\'))
  } catch {return $false}
}
function UbuntuReady {
  if(!(FeatureOn 'Microsoft-Windows-Subsystem-Linux')){return $false}
  if(!(FeatureOn 'VirtualMachinePlatform')){return $false}
  if(!(UbuntuRegistered) -or !(UbuntuOwnedByProject)){return $false}
  # Verify the actual OS image as well as its registry version, Windows
  # default distro, and any existing Linux on the host are irrelevant.
  $release=(& wsl.exe -d $Ubuntu --exec bash -lc '. /etc/os-release; printf "%s:%s" "$ID" "$VERSION_ID"' 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0 -or $release -ne 'ubuntu:24.04'){return $false}
  $arch=(& wsl.exe -d $Ubuntu --exec uname -m 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0 -or $arch -ne 'x86_64'){return $false}
  $name=(& wsl.exe -d $Ubuntu --exec id -un 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0 -or $name -ne 'rdmva'){return $false}
  # PowerShell's $HOME is a read-only automatic variable.
  $linuxHomePath=(& wsl.exe -d $Ubuntu --exec bash -lc 'printf %s "$HOME"' 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0 -or $linuxHomePath -ne '/home/rdmva'){return $false}
  return $true
}
function InitializeResearchUser {
  $researchUser='rdmva'
  Write-Host ('RDMVA_WSL_CREATE_UNPRIVILEGED_USER='+$researchUser)
  & wsl.exe -d $Ubuntu --user root --exec bash -c 'id -u rdmva >/dev/null 2>&1 || useradd -m -s /bin/bash rdmva'
  if($LASTEXITCODE -ne 0){
    Record 'WSL_SETUP_REQUIRED' 'Research account creation failed; the installed distro is preserved for a safe retry.'
    return $false
  }
  & wsl.exe --manage $Ubuntu --set-default-user $researchUser
  if($LASTEXITCODE -ne 0){
    Record 'WSL_SETUP_REQUIRED' 'Setting the default research account failed; safe retry is supported.'
    return $false
  }
  # Restart only the isolated research distro to apply its default account.
  & wsl.exe --terminate $Ubuntu
  if($LASTEXITCODE -ne 0){
    Record 'WSL_SETUP_REQUIRED' ('Restarting isolated '+$Ubuntu+' failed; safe retry is supported.')
    return $false
  }
  $actualUser=(& wsl.exe -d $Ubuntu --exec id -un 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0 -or $actualUser -ne $researchUser){
    Record 'WSL_SETUP_REQUIRED' 'Default research account verification failed.'
    return $false
  }
  $actualHome=(& wsl.exe -d $Ubuntu --exec bash -lc 'printf %s "$HOME"' 2>&1 | Out-String).Trim()
  if($LASTEXITCODE -ne 0 -or $actualHome -ne ('/home/'+$researchUser)){
    Record 'WSL_SETUP_REQUIRED' 'Research account home verification failed.'
    return $false
  }
  Write-Host ('RDMVA_WSL_UNPRIVILEGED_ACCOUNT_READY='+$actualUser)
  return $true
}

try {
  Record 'START' 'Check clean Windows 11 x64 prerequisites.'
  if(![Environment]::Is64BitOperatingSystem){Abort 'Windows x64 required'}
  $os=Get-CimInstance Win32_OperatingSystem
  if([int]$os.BuildNumber -lt 22000){Abort 'Tested profile requires Windows 11'}
  $principal=New-Object Security.Principal.WindowsPrincipal([Security.Principal.WindowsIdentity]::GetCurrent())
  if(!$principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)){Abort 'Open PowerShell as Administrator, then rerun installer'}
  if(!(Get-Command winget.exe -ErrorAction SilentlyContinue)){Abort 'WinGet missing: install Microsoft App Installer'}
  Write-Host ('WIN_BUILD='+$os.BuildNumber+' WINGET='+((& winget.exe --version)|Out-String).Trim())

  # Fail early on insufficient scientific resource storage, before expensive
  # downloads. We never create/format disks or relocate an existing checkout.
  if(!$PrerequisitesOnly){
    $storageScript=Join-Path $Root 'CONFIGURE_STORAGE_V31.ps1'
    if(!(Test-Path -LiteralPath $storageScript -PathType Leaf)){Abort 'Missing storage auto-configuration script'}
    Write-Host 'RDMVA_STORAGE_AUTO_SELECTION_START'
    $storageArgs=@('-NoProfile','-ExecutionPolicy','Bypass','-File',$storageScript)
    if($DataDrive){$storageArgs+=@('-DataDrive',$DataDrive)}
    & powershell.exe @storageArgs
    if($LASTEXITCODE -ne 0){Abort 'Automatic public research storage selection failed: need an existing formatted drive with >=340 GiB available'}
  }
  $version=InstalledVersion 'py.exe' @('-3.12','--version')
  if(!$version){
    WingetEnsure 'Python.Python.3.12' 'install'
    $version=InstalledVersion 'py.exe' @('-3.12','--version')
    if(!$version){Record 'REOPEN_SHELL' 'Python installed; reopen administrator PowerShell and rerun installer.';exit 10}
  }
  $v=PythonPatch $version
  if($v -lt $PythonMinVersion){
    WingetEnsure 'Python.Python.3.12' 'upgrade'
    $version=InstalledVersion 'py.exe' @('-3.12','--version')
    if(!$version){Record 'REOPEN_SHELL' 'Python updated; reopen administrator PowerShell and rerun installer.';exit 10}
    $v=PythonPatch $version
  }
  if($v -lt $PythonMinVersion){Abort ('V31_PYTHON_UPGRADE_NOT_EFFECTIVE: '+$version)}
  Write-Host ('PYTHON_OK='+$version+' MIN_REQUIRED='+$PythonMinVersion)

  $version=InstalledVersion 'git.exe' @('--version')
  if(!$version){
    WingetEnsure 'Git.Git' 'install'
    $version=InstalledVersion 'git.exe' @('--version')
    if(!$version){Record 'REOPEN_SHELL' 'Git installed; reopen administrator PowerShell and rerun installer.';exit 10}
  }
  $v=GitPatch $version
  if($v -lt $GitMinVersion){
    # A 2.52 installation is NOT "absent"; upgrading is necessary, and
    # no older global Git installation is silently accepted or deleted.
    WingetEnsure 'Git.Git' 'upgrade'
    $version=InstalledVersion 'git.exe' @('--version')
    if(!$version){Record 'REOPEN_SHELL' 'Git updated; reopen administrator PowerShell and rerun installer.';exit 10}
    $v=GitPatch $version
  }
  if($v -lt $GitMinVersion){
    Abort ('V31_GIT_UPGRADE_NOT_EFFECTIVE: '+$version+' REQUIRED_MIN='+$GitVersion)
  }
  Write-Host ('GIT_OK='+$version+' MIN_REQUIRED='+$GitMinVersion)

  $reboot=$false
  foreach($feature in @('Microsoft-Windows-Subsystem-Linux','VirtualMachinePlatform')){
    if(!(FeatureOn $feature)){
      Write-Host ('ENABLE_FEATURE='+$feature)
      $result=Enable-WindowsOptionalFeature -Online -FeatureName $feature -All -NoRestart -ErrorAction Stop
      if($result.RestartNeeded){$reboot=$true}
    }
  }
  if($reboot){Record 'REBOOT_REQUIRED' 'Windows optional features changed. Save your work and restart Windows when safe, then rerun the same installer. Existing downloaded files are preserved.';exit 10}
  if(!(UbuntuReady)){
    $distroLocation=Join-Path $Root 'env\wsl_ubuntu_24_04'
    Write-Host ('RDMVA_ISOLATED_WSL_DISTRIBUTION='+$Ubuntu)
    if(UbuntuRegistered){
      # A registered distro may have been installed before interruption.
      # Resume only if the WSL registry binds it to THIS source directory.
      if(!(UbuntuOwnedByProject)){
        Record 'WSL_SETUP_REQUIRED' ('Isolated research distro '+$Ubuntu+' is registered but does not belong to this source folder. Existing distro preserved.')
        exit 10
      }
      Write-Host 'RDMVA_WSL_RESUMING_INTERRUPTED_NONINTERACTIVE_SETUP'
    } else {
      # On fresh Windows, optional features do not install the WSL engine.
      Write-Host 'RDMVA_WSL_ENGINE_UPDATE_START'
      & wsl.exe --update --web-download
      if($LASTEXITCODE -ne 0){
        Record 'WSL_UPDATE_FAILED' 'Official WSL engine update failed. Check connectivity and retry.'
        exit 10
      }
      Write-Host 'RDMVA_WSL_ENGINE_UPDATE_OK'
      & wsl.exe --set-default-version 2
      if($LASTEXITCODE -ne 0){
        Record 'WSL_SETUP_REQUIRED' 'WSL2 version selection failed. Retry after Windows reboot if prompted.'
        exit 10
      }
      if(Test-Path -LiteralPath $distroLocation) {
        Abort 'V31_NEW_WSL_DISTRIBUTION_LOCATION_ALREADY_EXISTS: refusing to overwrite'
      }
      Write-Host ('INSTALLING_WSL2_UBUNTU_24_04_LOCATION='+$distroLocation)
      Write-Host 'RDMVA_WSL_NONINTERACTIVE_FIRST_LAUNCH=YES'
      & wsl.exe --install --web-download --no-launch -d $UbuntuSource --name $Ubuntu --location $distroLocation
      if($LASTEXITCODE -ne 0){
        Record 'WSL_SETUP_REQUIRED' ('Research distro '+$Ubuntu+' installation incomplete. No existing distribution deleted; retry after checking Windows requirements.')
        exit 10
      }
      if(!(UbuntuOwnedByProject)){
        Record 'WSL_SETUP_REQUIRED' 'New Ubuntu registration does not match expected project storage path. No distro changes performed.'
        exit 10
      }
    }
    if(!(InitializeResearchUser)){exit 10}
    if(!(UbuntuReady)){
      Record 'WSL_SETUP_REQUIRED' ('Isolated research distro '+$Ubuntu+' failed the post-provisioning functional check.')
      exit 10
    }
  }
  Write-Host 'WSL2_UBUNTU_24_04_READY'
  if($PrerequisitesOnly){Record 'PREREQUISITES_PASS' 'Python, Git and Ubuntu WSL2 installed and verified.';exit 0}
  $setup=Join-Path $Root 'setup_v2810.py'
  if(!(Test-Path -LiteralPath $setup -PathType Leaf)){Abort 'Missing setup_v2810.py in GitHub source'}
  $params=@('-3.12',$setup)
  if($SkipReferenceDownload){$params+='--skip-reference-download'}
  if($CurrentPublicData){$params+='--current-public-data'}
  $publicDataProfile=if($CurrentPublicData){'CURRENT_PUBLIC_DATA_WITH_NEW_SHA256_PROVENANCE'}else{'FROZEN_CHALLENGE_INPUTS'}
  Record 'SCIENTIFIC_INSTALL_START' ('Downloading/verifying '+$publicDataProfile+' and installing hash-pinned Python packages.')
  & py.exe @params
  if($LASTEXITCODE -ne 0){
    $setupExit=$LASTEXITCODE
    if($CurrentPublicData){
      $sourceReceipt=Join-Path $Root 'public_inputs\current_v31\CURRENT_PUBLIC_DATA_V31.json'
      if(Test-Path -LiteralPath $sourceReceipt -PathType Leaf){
        try{
          $sourceState=Get-Content -LiteralPath $sourceReceipt -Raw | ConvertFrom-Json
          if($sourceState.state -eq 'FAILED' -and
            [string]$sourceState.failure_summary -match
              'V31_CHEMBL_(UPSTREAM_TEMPORARY_HTTP_EXHAUSTED|UPSTREAM_TIMEOUT_EXHAUSTED|BAD_RESPONSE_EXHAUSTED)'){
            Record 'UPSTREAM_CHEMBL_RETRY_LATER' ('Official ChEMBL API unavailable after bounded retries. Verified public snapshots and ChEMBL responses remain cached, SHA-256 checked, and can resume safely by selecting Install again. No scientific run started.')
            exit 11
          }
        }catch{}
      }
    }
    Abort ('Scientific runtime setup failed; exit='+$setupExit+'; inspect INSTALL_V2810_STATUS.json and public source receipt')
  }
  $scientificPy=Join-Path $Root 'env\v28_runtime\Scripts\python.exe'
  $probe="import sys;from pathlib import Path;sys.path.insert(0,str(Path.cwd()/'code'));from v2810.resource_bootstrap import wsl_context,_public_storage_location;s=_public_storage_location(wsl_context(),Path.cwd());assert s['windows_root'] and s['free_bytes']>s['reserve_bytes'];print('RDMVA_SCIENTIFIC_STORAGE_WSL_READY='+str(s['windows_root']),flush=True)"
  Push-Location -LiteralPath $Root
  try{
    & $scientificPy -c $probe
    if($LASTEXITCODE -ne 0){Abort 'WSL2 cannot access the automatically selected storage volume'}
  } finally { Pop-Location }
  $receiptFile=Join-Path $Root 'public_inputs\current_v31\CURRENT_PUBLIC_DATA_V31.json'
  $sourceInfo=if(Test-Path $receiptFile){Get-Content $receiptFile -Raw | ConvertFrom-Json}else{$null}
  if($sourceInfo -and $sourceInfo.state -eq 'DEGRADED'){
    Write-Host 'PUBLIC_DATABASE_WARNING: Setup completed WITH MISSING DATABASES; see INSTALL_WARNINGS.json and menu.' -ForegroundColor Red
    Record 'PASS_WITH_SKIPPED_PUBLIC_SOURCES' 'Runtime available. Downloaded sources verified; missing scientific databases explicitly recorded.'
  } else {
    Record 'PASS' 'Pinned scientific environment and automatically verified external-capacity storage are ready.'
  }
  exit 0
} catch {
  Record 'FAILED' $_.Exception.Message
  Write-Host ('ERROR='+$_.Exception.Message) -ForegroundColor Red
  exit 1
}
