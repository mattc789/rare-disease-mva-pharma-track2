# Rare Disease MVA Pharma V3.1 — autonomous Challenge reproduction
# Windows PowerShell 5.1+, Windows x64, CPython 3.12, WSL2 Ubuntu.
param(
  [string]$CasePath,
  [string]$VcfPath,
  [string]$ClinicalPath,
  [string]$Out = "",
  [switch]$SkipSetup,
  [switch]$NoOnlineAI,
  [switch]$NoResourceDownload,
  [switch]$SkipFullScientificModules,
  [switch]$ResumeFailedRun,
  [string]$RecoverFrom = ""
)

$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Py = Join-Path $Root 'env\v28_runtime\Scripts\python.exe'
$runner = Join-Path $Root 'reproduce_v31_pipeline.py'
$currentReceipt = Join-Path $Root 'public_inputs\current_v31\CURRENT_PUBLIC_DATA_V31.json'
$currentMode = $false
if (Test-Path -LiteralPath $currentReceipt -PathType Leaf) {
  $meta = Get-Content -LiteralPath $currentReceipt -Raw | ConvertFrom-Json
  if ($meta.state -eq 'DEGRADED') {
    Write-Host 'DOWNLOAD WARNING: Public sources missing. Attempting best-effort research, with limited evidence coverage.' -ForegroundColor Red
    $currentMode = $true
    $script:IncompleteSources = $true
    $env:RDMVA_CURRENT_PUBLIC_DATA = '1'
  }
  if ($meta.state -eq 'READY') {
    $currentMode = $true
    $env:RDMVA_CURRENT_PUBLIC_DATA = '1'
    Write-Host 'V31_CURRENT_PUBLIC_DATA_PROFILE_ACTIVE'
  }
}

if (![bool]$CasePath -and ![bool]$VcfPath) {
  throw 'Supply -CasePath or -VcfPath with an authorized local input file.'
}
if ($CasePath -and $VcfPath) {
  throw 'Supply only one of -CasePath and -VcfPath.'
}
if (!(Test-Path -LiteralPath $runner -PathType Leaf)) {
  throw 'Incomplete source release: reproduce_v31_pipeline.py missing.'
}
if (!$NoOnlineAI -and [string]::IsNullOrWhiteSpace($env:OPENAI_API_KEY)) {
  $localEncryptedKey = Join-Path $env:LOCALAPPDATA 'Examedix\RareDiseaseMVAPharma\openai_key.dpapi'
  if (!(Test-Path -LiteralPath $localEncryptedKey -PathType Leaf)) {
    throw 'OPENAI_API_KEY is not set, and this Windows account has no local encrypted QA credential.'
  }
}

if (!$SkipSetup) {
  Write-Host 'V31_SETUP_START'
  $setup = Join-Path $Root 'setup_v2810.py'
  $setupArgs = @($setup)
  if ($currentMode) { $setupArgs += '--current-public-data' }
  if (Test-Path -LiteralPath $Py -PathType Leaf) {
    $setupArgs += '--verify-only'
  } elseif ($NoResourceDownload) {
    $setupArgs += '--skip-reference-download'
  }
  $old = $ErrorActionPreference
  try {
    $ErrorActionPreference = 'Continue'
    & py -3.12 @setupArgs
    $exit = $LASTEXITCODE
  } finally {
    $ErrorActionPreference = $old
  }
  if ($exit -ne 0) { throw "V31 setup failed with exit code $exit" }
}

if (!(Test-Path -LiteralPath $Py -PathType Leaf)) {
  throw 'Missing Python runtime. Run: py -3.12 .\setup_v2810.py'
}

$inputPath = if ($CasePath) { (Resolve-Path -LiteralPath $CasePath).Path } else { (Resolve-Path -LiteralPath $VcfPath).Path }
# Resolve default output against the source root, not the caller's working directory.
$outPath = if ([string]::IsNullOrWhiteSpace($Out)) {
  Join-Path $Root 'experiments\CLEANROOM_V31'
} else {
  [IO.Path]::GetFullPath($Out)
}
$arguments = @($runner,'--out',$outPath)
if ($CasePath) { $arguments += @('--case',$inputPath) } else { $arguments += @('--vcf',$inputPath) }
if ($ClinicalPath) { $arguments += @('--clinical',(Resolve-Path -LiteralPath $ClinicalPath).Path) }
if ($NoOnlineAI) { $arguments += '--no-online-ai' }
if ($NoResourceDownload) { $arguments += '--no-resource-download' }
if ($SkipFullScientificModules) { $arguments += '--skip-full-scientific-modules' }
if ($ResumeFailedRun) { $arguments += '--resume-failed-run' }
if (![string]::IsNullOrWhiteSpace($RecoverFrom)) {
  $arguments += @('--recover-from', (Resolve-Path -LiteralPath $RecoverFrom).Path)
}

Write-Host 'V31_COMPLETE_CHALLENGE_START'
$old = $ErrorActionPreference
try {
  $ErrorActionPreference = 'Continue'
  & $Py @arguments
  $exit = $LASTEXITCODE
} finally {
  $ErrorActionPreference = $old
}
if ($exit -ne 0) { throw "V31 Challenge reproduction failed, exit code $exit. Inspect E2E_EXECUTION_V31.json in output directory." }
if ($script:IncompleteSources) {
  Write-Host 'DOWNLOAD WARNING: Analysis used incomplete public sources; do not interpret it as a fully validated science run.' -ForegroundColor Red
  Write-Host 'V31_COMPLETE_CHALLENGE_WITH_MISSING_PUBLIC_DATA' -ForegroundColor Red
} else {
  Write-Host 'V31_COMPLETE_CHALLENGE_PASS'
}
