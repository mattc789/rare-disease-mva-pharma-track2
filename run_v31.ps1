# Rare Disease MVA Pharma — V3.1 Expanded Patient-First Rescue Engine
$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Python = Join-Path $Root 'env\v28_runtime\Scripts\python.exe'
if (!(Test-Path $Python)) { throw 'V3.1 runtime missing. Run: py -3.12 .\setup_v2810.py' }
if (!(Test-Path (Join-Path $Root 'config\v30_source_policy.json'))) { throw 'V3.1 source-rights policy missing' }
Set-Location (Join-Path $Root 'code')
$NativeExitCode = 1
$PreviousErrorPreference = $ErrorActionPreference
try {
  $ErrorActionPreference = 'Continue'
  & $Python -m v30.pipeline @args
  $NativeExitCode = $LASTEXITCODE
} finally {
  $ErrorActionPreference = $PreviousErrorPreference
}
exit $NativeExitCode
