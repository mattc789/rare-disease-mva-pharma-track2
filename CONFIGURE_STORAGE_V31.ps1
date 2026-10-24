# Rare Disease MVA Pharma V3.1: non-destructive automatic public-data storage selection.
# Never initializes, repartitions, formats, resizes or assigns disk letters.
# On a new machine, selects an existing mounted NTFS/exFAT volume with enough
# measured free space for the full 240 GiB gnomAD budget plus 100 GiB reserve.
# Machine-specific configuration is never bundled in the GitHub source release.
param([ValidatePattern('^[A-Za-z]$')][string]$DataDrive)
$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Config = Join-Path $Root 'config\local_public_storage_v2810.json'
$MaxGnomadGiB = 240
$ReservedGiB = 100
$RequiredFreeBytes = ($MaxGnomadGiB + $ReservedGiB) * 1GB
$ProjectRootLetter = ([System.IO.Path]::GetPathRoot($Root)).Substring(0,1).ToUpperInvariant()

function Get-Candidate([string]$letter) {
  $volume = Get-Volume -DriveLetter $letter -ErrorAction SilentlyContinue
  if (!$volume -or $volume.FileSystem -notin @('NTFS','exFAT')) { return $null }
  $partition = Get-Partition -DriveLetter $letter -ErrorAction SilentlyContinue
  if (!$partition -or @($partition).Count -ne 1) { return $null }
  $disk = $partition | Get-Disk -ErrorAction SilentlyContinue
  if (!$disk -or $disk.IsOffline -or $disk.IsReadOnly) { return $null }
  if ($volume.SizeRemaining -lt $RequiredFreeBytes) { return $null }
  if ($volume.Size -lt 100GB) { return $null }
  return [pscustomobject]@{
    Letter = $letter.ToUpperInvariant()
    DiskName = [string]$disk.FriendlyName
    DiskNumber = [int]$disk.Number
    FreeBytes = [Int64]$volume.SizeRemaining
    SystemDisk = [bool]($disk.IsBoot -or $disk.IsSystem)
    ProjectDrive = [bool]($letter.ToUpperInvariant() -eq $ProjectRootLetter)
    FileSystem = [string]$volume.FileSystem
  }
}

if (Test-Path -LiteralPath $Config -PathType Leaf) {
  $existing = Get-Content -LiteralPath $Config -Raw | ConvertFrom-Json
  if ($existing.schema -ne 'V2810_LOCAL_PUBLIC_STORAGE_1') {
    throw 'V31_STORAGE_CONFIG_UNKNOWN_SCHEMA: manual review required; not overwriting'
  }
  $path = [string]$existing.windows_root
  if ($path -notmatch '^[A-Za-z]:\\[^\\/?*<>|:]+') {
    throw 'V31_STORAGE_CONFIG_PATH_INVALID'
  }
  if ($DataDrive -and $path.Substring(0,1).ToUpperInvariant() -ne $DataDrive.ToUpperInvariant()) {
    throw 'V31_STORAGE_EXISTING_SELECTION_CONFLICT: will not overwrite'
  }
  $letter = $path.Substring(0,1).ToUpperInvariant()
  $volume = Get-Volume -DriveLetter $letter -ErrorAction Stop
  $disk = Get-Partition -DriveLetter $letter -ErrorAction Stop | Get-Disk -ErrorAction Stop
  if ($volume.FileSystem -notin @('NTFS','exFAT') -or $disk.IsOffline -or $disk.IsReadOnly -or
      ($existing.expected_disk_name -and $disk.FriendlyName -ne $existing.expected_disk_name)) {
    throw 'V31_STORAGE_EXISTING_DRIVE_CHANGED_OR_INVALID'
  }
  if ([Int64]$volume.SizeRemaining -lt [int]$existing.reserve_free_gib * 1GB) {
    throw 'V31_STORAGE_BELOW_RESERVED_FREE_SPACE'
  }
  Write-Host ('V31_STORAGE_EXISTING_CONFIG_VERIFIED='+$path)
  exit 0
}

$candidates = @()
$volumes = @(Get-Volume -ErrorAction Stop | Where-Object {
    $_.DriveLetter -and $_.FileSystem -in @('NTFS','exFAT')
  })
foreach ($v in $volumes) {
  $letter = ([string]$v.DriveLetter).ToUpperInvariant()
  if ($DataDrive -and $letter -ne $DataDrive.ToUpperInvariant()) { continue }
  $candidate = Get-Candidate $letter
  if ($candidate) { $candidates += $candidate }
}
if ($candidates.Count -eq 0) {
  Write-Host 'V31_STORAGE_NO_ELIGIBLE_VOLUME'
  Write-Host ('V31_STORAGE_REQUIRES_FREE_GIB='+($MaxGnomadGiB+$ReservedGiB))
  throw 'V31_STORAGE_REQUIRED: Attach or select an already formatted NTFS/exFAT drive with at least 340 GiB free, then rerun the installer. No existing disk will be formatted.'
}
# Prefer the project checkout volume when it already has sufficient capacity:
# this makes all project assets and large public references a single local
# folder tree. Otherwise safely choose a separate formatted data volume.
$chosen = $candidates | Sort-Object -Property @{
    Expression='ProjectDrive';Descending=$true
  },@{
    Expression='SystemDisk';Descending=$false
  },@{
    Expression='FreeBytes';Descending=$true
  },@{
    Expression='Letter';Descending=$false
  } | Select-Object -First 1
$StorageRoot = if ($chosen.ProjectDrive) {
    Join-Path $Root 'public_inputs\third_party'
  } else {
    $chosen.Letter + ':\RDMVA_Public_References'
  }
$settings = [ordered]@{
  schema='V2810_LOCAL_PUBLIC_STORAGE_1'
  windows_root=$StorageRoot
  expected_disk_name=$chosen.DiskName
  reserve_free_gib=$ReservedGiB
  require_external_usb=$false
  source_zip_include=$false
  configured_by='INSTALL_WINDOWS_V31_AUTO_STORAGE'
}
if (Test-Path -LiteralPath $Config) { throw 'V31_STORAGE_CONFIG_CREATED_CONCURRENTLY' }
New-Item -Path $StorageRoot -ItemType Directory -Force -ErrorAction Stop | Out-Null
$json = $settings | ConvertTo-Json -Depth 4
$utf8 = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText($Config,$json + [Environment]::NewLine,$utf8)
Write-Host ('V31_STORAGE_AUTOMATIC_SELECTED='+$StorageRoot)
Write-Host ('V31_STORAGE_LAYOUT='+$(if($chosen.ProjectDrive){'SINGLE_PROJECT_ROOT'}else{'SEPARATE_HIGH_CAPACITY_DATA_VOLUME'}))
Write-Host ('V31_STORAGE_FREE_GIB='+[math]::Round($chosen.FreeBytes/1GB,1))
Write-Host ('V31_STORAGE_RESERVED_FREE_GIB='+$ReservedGiB)
Write-Host 'V31_STORAGE_NO_DISK_FORMATTING_PERFORMED=YES'
