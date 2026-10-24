# Rare Disease MVA Pharma V3.1: save an OpenAI API key encrypted under
# the currently signed-in Windows user (DPAPI CurrentUser).
# Do not run this as a different administrator account from the account
# that will execute REPRODUCE_CHALLENGE_V31.ps1.
# No .env file and no plaintext credential is saved in the source directory.
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Security -ErrorAction Stop
$destination = Join-Path $env:LOCALAPPDATA 'Examedix\RareDiseaseMVAPharma\openai_key.dpapi'
if (Test-Path -LiteralPath $destination -PathType Leaf) {
    throw 'V31_LOCAL_API_CREDENTIAL_ALREADY_EXISTS: not overwriting encrypted key.'
}
$secret = Read-Host 'Enter your billable OpenAI API key (input hidden)' -AsSecureString
$pointer = [IntPtr]::Zero
$raw = $null
try {
  $pointer=[Runtime.InteropServices.Marshal]::SecureStringToBSTR($secret)
  $text=[Runtime.InteropServices.Marshal]::PtrToStringBSTR($pointer)
  if ([string]::IsNullOrWhiteSpace($text) -or $text.Length -lt 20 -or $text.Length -gt 1024) {
    throw 'V31_API_KEY_FORMAT_INVALID'
  }
  $raw=[System.Text.Encoding]::UTF8.GetBytes($text)
  $sealed=[Security.Cryptography.ProtectedData]::Protect(
      $raw,$null,[Security.Cryptography.DataProtectionScope]::CurrentUser)
  $parent=Split-Path -Parent $destination
  New-Item -Path $parent -ItemType Directory -Force -ErrorAction Stop | Out-Null
  if (Test-Path -LiteralPath $destination -PathType Leaf) {
    throw 'V31_API_CREDENTIAL_CONCURRENTLY_CREATED'
  }
  [System.IO.File]::WriteAllBytes($destination,$sealed)
  Write-Host 'V31_OPENAI_USER_DPAPI_KEY_STORED=YES'
  Write-Host 'V31_API_KEY_PLAINTEXT_FILE_CREATED=NO'
} finally {
  if ($raw) {[Array]::Clear($raw,0,$raw.Length)}
  if ($pointer -ne [IntPtr]::Zero) {
    [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($pointer)
  }
  $secret.Dispose()
}
