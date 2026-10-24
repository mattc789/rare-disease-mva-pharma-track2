# Single Windows entry point for the Rare Disease MVA Pharma V3.1 judge.
# Local patient data never leave the computer through this launcher.
[CmdletBinding()]
param(
    [ValidateSet('menu','check','install','key','run','monitor','test','models')]
    [string]$Action='menu',
    [string]$VcfPath='',
    [string]$ClinicalPath='',
    [string]$Out='',
    [string]$DataDrive=''
)
$ErrorActionPreference='Stop'
$Root=Split-Path -Parent $MyInvocation.MyCommand.Path
$Installer=Join-Path $Root 'INSTALL_WINDOWS_V31.ps1'
$Runner=Join-Path $Root 'REPRODUCE_CHALLENGE_V31.ps1'
$KeySetup=Join-Path $Root 'SETUP_OPENAI_KEY_V31.ps1'
$Monitor=Join-Path $Root 'JURY_MONITOR_V31.py'
$Runtime=Join-Path $Root 'env\v28_runtime\Scripts\python.exe'
function Show-DownloadWarnings {
    # Persistent red warning inside the existing 0-7 menu (no new selection).
    $statusFile=Join-Path $Root 'INSTALL_WARNINGS.json'
    if (!(Test-Path -LiteralPath $statusFile -PathType Leaf)) {return}
    try {
        $data=Get-Content -LiteralPath $statusFile -Raw | ConvertFrom-Json
        $missing=@($data.skipped_sources)
        if ($missing.Count -gt 0) {
            Write-Host ''
            Write-Host ('DOWNLOAD WARNING: '+$missing.Count+' public dataset(s) could not be downloaded.') -ForegroundColor Red
            Write-Host 'Setup CONTINUED WITHOUT these datasets. Scientific coverage is incomplete:' -ForegroundColor Red
            foreach ($entry in $missing) {
                Write-Host ('  MISSING: '+[string]$entry.name) -ForegroundColor Red
                Write-Host ('  OFFICIAL LINK: '+[string]$entry.official_url) -ForegroundColor Red
            }
            Write-Host 'Unavailable reference data are never silently treated as successfully downloaded.' -ForegroundColor Red
        }
    } catch {
        Write-Host 'DOWNLOAD WARNING: INSTALL_WARNINGS.json could not be read.' -ForegroundColor Red
    }
}
function Check-Prerequisites {
    Write-Host ""
    Write-Host "=== Rare Disease MVA Pharma V3.1 - Judge preflight ===" -ForegroundColor Cyan
    $ok=$true
    if (-not [Environment]::Is64BitOperatingSystem) {Write-Host "WINDOWS_X64_REQUIRED" -ForegroundColor Red;$ok=$false}
    $os=Get-CimInstance Win32_OperatingSystem
    Write-Host ("WINDOWS_BUILD="+$os.BuildNumber)
    if ([int]$os.BuildNumber -lt 22000) {Write-Host "WINDOWS_11_REQUIRED" -ForegroundColor Red;$ok=$false}
    foreach ($p in @($Installer,$Runner,$KeySetup,$Monitor)) {
        if (!(Test-Path -LiteralPath $p -PathType Leaf)) {
            Write-Host ("MISSING_BUNDLED_SCRIPT="+[IO.Path]::GetFileName($p)) -ForegroundColor Red
            $ok=$false
        }
    }
    Write-Host ("PINNED_PYTHON_RUNTIME="+(Test-Path -LiteralPath $Runtime -PathType Leaf))
    if (!(Test-Path -LiteralPath $Runtime -PathType Leaf)) {$ok=$false}
    $receipt=Join-Path $Root 'public_inputs\current_v31\CURRENT_PUBLIC_DATA_V31.json'
    if (Test-Path -LiteralPath $receipt -PathType Leaf) {
        try {
            $j=Get-Content -LiteralPath $receipt -Raw | ConvertFrom-Json
            Write-Host ("PUBLIC_REFERENCE_SNAPSHOT="+$j.state)
            if ($j.state -ne 'READY') {$ok=$false}
        } catch { Write-Host 'PUBLIC_REFERENCE_SNAPSHOT=INVALID' -ForegroundColor Red;$ok=$false }
    } else { Write-Host 'PUBLIC_REFERENCE_SNAPSHOT=NOT_INSTALLED' -ForegroundColor Yellow;$ok=$false }
    $key=Join-Path $env:LOCALAPPDATA 'Examedix\RareDiseaseMVAPharma\openai_key.dpapi'
    if (($env:OPENAI_API_KEY -and $env:OPENAI_API_KEY.Trim()) -or
        (Test-Path -LiteralPath $key -PathType Leaf)) {
        Write-Host 'API_CREDENTIAL_PRESENT=YES (usable status checked at runtime)'
    } else {Write-Host 'API_CREDENTIAL_PRESENT=NO' -ForegroundColor Yellow;$ok=$false}
    $storage=Join-Path $Root 'config\local_public_storage_v2810.json'
    if (Test-Path -LiteralPath $storage -PathType Leaf) {
        try {
            $settings=Get-Content -LiteralPath $storage -Raw | ConvertFrom-Json
            Write-Host ("SCIENTIFIC_STORAGE="+$settings.windows_root)
        } catch {Write-Host 'SCIENTIFIC_STORAGE_CONFIG=INVALID' -ForegroundColor Red;$ok=$false}
    } else {Write-Host 'SCIENTIFIC_STORAGE=NOT_CONFIGURED' -ForegroundColor Yellow;$ok=$false}
    if ($ok) {
        Write-Host 'RDMVA_JURY_PREFLIGHT_RUNTIME_PREREQUISITES_PRESENT' -ForegroundColor Green
    } else {Write-Host 'RDMVA_JURY_PREFLIGHT_INCOMPLETE' -ForegroundColor Yellow}
    Write-Host 'NOTE: full VEP/gnomAD evidence is verified during the scientific analysis, not here.'
    return $ok
}
function Install {
    $principal=New-Object Security.Principal.WindowsPrincipal([Security.Principal.WindowsIdentity]::GetCurrent())
    if (!$principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
        throw 'Reopen this launcher as Administrator in the SAME Windows account.'
    }
    $argsList=@('-NoProfile','-ExecutionPolicy','Bypass','-File',$Installer,'-CurrentPublicData')
    if ($DataDrive) {$argsList+=@('-DataDrive',$DataDrive)}
    & powershell.exe @argsList
    $status=$LASTEXITCODE
    if ($status -eq 11) {
        Write-Host 'RDMVA_JURY_PUBLIC_UPSTREAM_SERVICE_TEMPORARILY_UNAVAILABLE' -ForegroundColor Yellow
        Write-Host 'Official ChEMBL did not recover within the automated retry limit.'
        Write-Host 'Verified public inputs and ChEMBL pages are preserved. When the service is available, select Install again; do NOT delete anything.'
        Write-Host 'The scientific run has not started and no incomplete source data are approved.'
        if ($Action -eq 'install') {exit 11}
        return
    }
    if ($status -eq 10) {
        $installReceipt=Join-Path $Root 'INSTALL_WINDOWS_V31_STATUS.json'
        $reason='SETUP_REQUIRED'
        if (Test-Path -LiteralPath $installReceipt -PathType Leaf) {
            try {
                $reason=(Get-Content -LiteralPath $installReceipt -Raw |
                    ConvertFrom-Json).state
            } catch {
                $reason='SETUP_REQUIRED'
            }
        }
        Write-Host ("RDMVA_JURY_INSTALL_ACTION_REQUIRED="+$reason) -ForegroundColor Yellow
        if ($reason -eq 'REBOOT_REQUIRED') {
            Write-Host 'Windows must restart once to activate WSL. Save your work and restart when safe.' -ForegroundColor Yellow
            Write-Host 'After restarting, select Install again. Already installed files and downloads are preserved.'
        } else {
            Write-Host 'Complete the Windows or Ubuntu setup step indicated in INSTALL_WINDOWS_V31_STATUS.json, then select Install again.'
        }
        if ($Action -eq 'install') {exit 10}
        return
    }
    if ($status -ne 0) {throw "RDMVA_JURY_INSTALL_FAILED_EXIT_"+$status}
    Show-DownloadWarnings
    $srcReceipt=Join-Path $Root 'public_inputs\current_v31\CURRENT_PUBLIC_DATA_V31.json'
    if (Test-Path $srcReceipt -PathType Leaf) {
        $srcStatus=Get-Content $srcReceipt -Raw | ConvertFrom-Json
        if ($srcStatus.state -eq 'DEGRADED') {
            Write-Host 'RDMVA_JURY_INSTALL_COMPLETED_WITH_MISSING_PUBLIC_SOURCES' -ForegroundColor Red
            return
        }
    }
    Write-Host 'RDMVA_JURY_INSTALL_PASS' -ForegroundColor Green
}
function Configure-Key {
    $key=Join-Path $env:LOCALAPPDATA 'Examedix\RareDiseaseMVAPharma\openai_key.dpapi'
    if (Test-Path -LiteralPath $key -PathType Leaf) {
        Write-Host 'V31_API_CREDENTIAL_ALREADY_PRESENT=YES; use model preflight (option 3) to check it.'
        return
    }
    Write-Host 'API credential stays private to your Windows user via DPAPI.'
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $KeySetup
    if ($LASTEXITCODE -ne 0) {throw "RDMVA_KEY_SETUP_FAILED_EXIT_"+$LASTEXITCODE}
}
function Run-SelfTests {
    if (!(Test-Path -LiteralPath $Runtime -PathType Leaf)) {
        throw 'RDMVA_SELFTEST_INSTALL_SCIENTIFIC_RUNTIME_FIRST'
    }
    $codeRoot=Join-Path $Root 'code'
    Push-Location $codeRoot
    try {
        Write-Host 'Starting autonomous regression tests (synthetic fixtures; no Challenge patient files).' -ForegroundColor Cyan
        & $Runtime -m unittest discover -s . -p 'test_*.py' -q
        if ($LASTEXITCODE -ne 0) { throw 'RDMVA_JURY_SELF_TEST_FAILED' }
        Write-Host 'RDMVA_JURY_SELF_TEST_PASS' -ForegroundColor Green
    } finally {
        Pop-Location
    }
}
function Test-LiveModels {
    if (!(Test-Path -LiteralPath $Runtime -PathType Leaf)) {
        throw 'Install the pinned scientific runtime first.'
    }
    Write-Host 'This test makes a small, PAID OpenAI Responses call for each configured model.'
    Write-Host 'It uses synthetic public text only, never the Challenge patient genome.'
    if ((Read-Host 'Type YES (any capitalization) to authorize the small paid model tests').Trim() -ine 'YES') {
        Write-Host 'No paid API request has been made.'
        return
    }
    & $Runtime (Join-Path $Root 'JURY_CHECK_AI_V31.py') --verify-calls
    if ($LASTEXITCODE -ne 0) { throw 'One or more required AI models is not usable; do not start the full run.' }
    Write-Host 'ALL_REQUIRED_OPENAI_MODELS_VERIFIED' -ForegroundColor Green
}
# Search organizer-provided input filenames locally; only SHA-256 duplicate
# comparisons read file bytes before consent, never parse or upload patient
# data. Distinct pairs still require manual selection. Inputs stay local.
function Discover-AuthorizedInputs {
    Write-Host 'RDMVA_LOCAL_INPUT_DISCOVERY_START (local names and duplicate SHA-256 only; no patient content uploaded)'
    $roots=@()
    foreach ($location in @(
        $env:RDMVA_CHALLENGE_INPUT_DIR,
        (Join-Path $Root 'private_inputs'),
        (Join-Path $Root 'inputs'),
        (Join-Path $env:USERPROFILE 'Downloads'),
        (Join-Path $env:USERPROFILE 'Documents'),
        (Join-Path $env:USERPROFILE 'Desktop')
    )) {
        if ($location -and (Test-Path -LiteralPath $location -PathType Container)) {
            $roots+=([IO.Path]::GetFullPath($location))
        }
    }
    # Recognize other local research checkouts, regardless of drive letter.
    # Do not scan all files on every drive or include system/source databases.
    foreach ($disk in @(Get-PSDrive -PSProvider FileSystem -ErrorAction SilentlyContinue)) {
        if ($disk.Name -eq 'C') {continue}
        foreach ($folder in @(Get-ChildItem -LiteralPath $disk.Root -Directory -ErrorAction SilentlyContinue)) {
            if ($folder.Name -match '(?i)(rare.*disease|rdmva|challenge)') {
                $roots+=$folder.FullName
            }
        }
    }
    $seen=@{}
    $pairs=@()
    foreach ($directory in @($roots | Select-Object -Unique)) {
        foreach ($vcf in @(Get-ChildItem -LiteralPath $directory -Filter '*.vcf.gz' -Recurse -Depth 6 -File -ErrorAction SilentlyContinue)) {
            if ($vcf.FullName -match '(?i)\\(public_inputs|env|experiments|release|runs|\.git|gnomad|vep116|site-packages)\\') {continue}
            $id=$vcf.FullName.ToLowerInvariant()
            if ($seen.ContainsKey($id)) {continue}
            $seen[$id]=$true
            $docs=@(Get-ChildItem -LiteralPath $vcf.DirectoryName -Filter '*.docx' -File -ErrorAction SilentlyContinue)
            if ($docs.Count -gt 1) {
                $docs=@($docs | Where-Object {$_.Name -match '(?i)(clinical|phenotype|challenge|patient|case)'})
            }
            if ($docs.Count -eq 1) {
                $pairs+=([pscustomobject]@{VCF=$vcf.FullName;Clinical=$docs[0].FullName})
            }
        }
    }
    # Multiple paths may contain byte-identical copies of the same authorized
    # pair. SHA-256 comparison is local only; no clinical content is parsed or
    # uploaded. Never silently choose between genuinely different pairs.
    if ($pairs.Count -gt 1) {
        $originalPairs=$pairs.Count
        $unique=@{}
        foreach ($pair in $pairs) {
            $vcfPin=(Get-FileHash -LiteralPath $pair.VCF -Algorithm SHA256).Hash
            $docPin=(Get-FileHash -LiteralPath $pair.Clinical -Algorithm SHA256).Hash
            $key=$vcfPin+'|'+$docPin
            if (-not $unique.ContainsKey($key)) {$unique[$key]=$pair}
        }
        $pairs=@($unique.Values)
        if ($pairs.Count -eq 1) {
            Write-Host ('RDMVA_DUPLICATE_AUTHORIZED_PAIRS_DEDUPLICATED='+$originalPairs)
        }
    }
    # The VCF and clinical DOCX may be in DIFFERENT folders. When no adjacent
    # pair exists, search only authorized candidate locations; select a unique
    # one-of-each by content, not merely by whichever file was found first.
    if ($pairs.Count -eq 0) {
        $vcfCandidates=@{}
        $docCandidates=@{}
        foreach ($directory in @($roots | Select-Object -Unique)) {
            foreach ($f in @(Get-ChildItem -LiteralPath $directory -Filter '*.vcf.gz' -Recurse -Depth 6 -File -ErrorAction SilentlyContinue)) {
                if ($f.FullName -match '(?i)\\(public_inputs|env|experiments|release|runs|\.git|gnomad|vep116|site-packages)\\') {continue}
                $id=$f.FullName.ToLowerInvariant()
                if (-not $vcfCandidates.ContainsKey($id)) {$vcfCandidates[$id]=$f.FullName}
            }
            foreach ($f in @(Get-ChildItem -LiteralPath $directory -Filter '*.docx' -Recurse -Depth 6 -File -ErrorAction SilentlyContinue)) {
                if ($f.FullName -match '(?i)\\(public_inputs|env|experiments|release|runs|\.git|site-packages)\\') {continue}
                if ($f.Name -notmatch '(?i)(clinical|phenotype|challenge|patient|case)') {continue}
                $id=$f.FullName.ToLowerInvariant()
                if (-not $docCandidates.ContainsKey($id)) {$docCandidates[$id]=$f.FullName}
            }
        }
        if ($vcfCandidates.Count -gt 0 -and $docCandidates.Count -gt 0) {
            $vcfUnique=@{}
            foreach ($f in $vcfCandidates.Values) {
                $pin=(Get-FileHash -LiteralPath $f -Algorithm SHA256).Hash
                if (-not $vcfUnique.ContainsKey($pin)) {$vcfUnique[$pin]=$f}
            }
            $docUnique=@{}
            foreach ($f in $docCandidates.Values) {
                $pin=(Get-FileHash -LiteralPath $f -Algorithm SHA256).Hash
                if (-not $docUnique.ContainsKey($pin)) {$docUnique[$pin]=$f}
            }
            if ($vcfUnique.Count -eq 1 -and $docUnique.Count -eq 1) {
                $pairs=@([pscustomobject]@{
                    VCF=@($vcfUnique.Values)[0]
                    Clinical=@($docUnique.Values)[0]
                })
                Write-Host 'RDMVA_UNIQUE_AUTHORIZED_SEPARATE_FOLDER_PAIR_FOUND'
            } else {
                Write-Host ('RDMVA_AMBIGUOUS_SOURCE_INPUTS_VCF='+$vcfUnique.Count+' DOCX='+$docUnique.Count)
            }
        }
    }
    if ($pairs.Count -eq 1) {
        Write-Host ('RDMVA_LOCAL_INPUT_PAIR_FOUND_VCF='+$pairs[0].VCF)
        Write-Host ('RDMVA_LOCAL_INPUT_PAIR_FOUND_DOCX='+$pairs[0].Clinical)
        Write-Host 'Confirm that you are authorized to use these organizer-provided files.'
        if ((Read-Host 'Use this VCF and clinical DOCX? Type YES (any capitalization)').Trim() -ieq 'YES') {
            $script:VcfPath=$pairs[0].VCF
            $script:ClinicalPath=$pairs[0].Clinical
            Write-Host 'RDMVA_AUTHORIZED_INPUT_PAIR_CONFIRMED'
            return $true
        }
    } elseif ($pairs.Count -gt 1) {
        Write-Host ('RDMVA_MULTIPLE_INPUT_PAIRS_FOUND='+$pairs.Count)
        Write-Host 'Automatic selection skipped to avoid mixing patients; select the authorized files manually.'
    } else {
        Write-Host 'RDMVA_NO_LOCAL_ORGANIZER_PAIR_FOUND; use the file picker.'
    }
    return $false
}
function Choose-File([string]$Title,[string]$Filter,[string]$InitialDirectory='') {
    try {
        Add-Type -AssemblyName System.Windows.Forms -ErrorAction Stop
        $dialog=New-Object Windows.Forms.OpenFileDialog
        $dialog.Title=$Title
        $dialog.Filter=$Filter
        if ($InitialDirectory -and (Test-Path -LiteralPath $InitialDirectory -PathType Container)) {
            $dialog.InitialDirectory=$InitialDirectory
        }
        $dialog.CheckFileExists=$true
        $dialog.Multiselect=$false
        $result=$dialog.ShowDialog()
        if ($result -eq [Windows.Forms.DialogResult]::OK) {return $dialog.FileName}
    } catch {
        Write-Host 'Windows file picker unavailable, using manual path selection.'
    }
    return [string](Read-Host $Title)
}
function Start-Research {
    if (!(Test-Path -LiteralPath $Runtime -PathType Leaf)) {
        throw 'Install the scientific environment first (menu option 1).'
    }
    $modelProof=Join-Path $env:LOCALAPPDATA 'Examedix\RareDiseaseMVAPharma\JURY_AI_MODEL_CHECK_V31.json'
    if (!(Test-Path -LiteralPath $modelProof -PathType Leaf)) {
        throw 'Run menu option 3 (live API model check) before long scientific analysis.'
    }
    try {
        $modelStatus=Get-Content -LiteralPath $modelProof -Raw | ConvertFrom-Json
        $age=([DateTime]::UtcNow-[DateTime]::Parse([string]$modelStatus.utc).ToUniversalTime()).TotalDays
        if ($modelStatus.passed -ne $true -or $modelStatus.billable_synthetic_calls -ne $true -or
            $age -lt 0 -or $age -gt 7) {
            throw 'stale or failed model check'
        }
    } catch {
        throw 'Model readiness not verified in last 7 days; run menu option 3.'
    }
    if (!$VcfPath -and !$ClinicalPath) {
        $null=Discover-AuthorizedInputs
    }
    if (!$VcfPath) {
        $script:VcfPath=Choose-File 'Choose the authorized compressed VCF' 'Compressed VCF (*.vcf.gz)|*.vcf.gz|All files (*.*)|*.*'
    }
    if (!$ClinicalPath) {
        $preferredDirectory=''
        if ($VcfPath -and (Test-Path -LiteralPath $VcfPath -PathType Leaf)) {
            $preferredDirectory=Split-Path -Parent $VcfPath
        }
        $script:ClinicalPath=Choose-File 'Choose the authorized clinical DOCX' 'Word files (*.docx)|*.docx|All files (*.*)|*.*' $preferredDirectory
    }
    if (!(Test-Path -LiteralPath $VcfPath -PathType Leaf) -or !(Test-Path -LiteralPath $ClinicalPath -PathType Leaf)) {
        throw 'RDMVA_AUTHORIZED_INPUTS_NOT_FOUND'
    }
    if (!$VcfPath.EndsWith('.vcf.gz',[StringComparison]::OrdinalIgnoreCase)) {
        throw 'RDMVA_REQUIRES_COMPRESSED_VCF_GRCH38'
    }
    if (!$ClinicalPath.EndsWith('.docx',[StringComparison]::OrdinalIgnoreCase)) {
        throw 'RDMVA_REQUIRES_CLINICAL_DOCX'
    }
    if (!$Out) {
        $script:Out=Join-Path $Root ('experiments\JURY_V31_'+(Get-Date -Format 'yyyyMMdd_HHmmss'))
        Write-Host ("DEFAULT_OUTPUT="+$Out)
        $override=Read-Host 'Press Enter to use this NEW folder, or enter another NEW folder'
        if (![string]::IsNullOrWhiteSpace($override)) {$script:Out=$override.Trim('"')}
    }
    $script:Out=[IO.Path]::GetFullPath($Out)
    if (Test-Path -LiteralPath $Out) {throw 'Output directory already exists; choose a NEW directory.'}
    Write-Host 'Starting full genomic screening and paid OpenAI API research.'
    Write-Host 'Experimental research hypotheses ONLY: not medical recommendations.'
    if ((Read-Host 'Type YES (any capitalization) to start the scientific analysis').Trim() -ine 'YES') {
        Write-Host 'No analysis started.'
        return
    }
    Write-Host ("RUN_OUTPUT="+$Out)
    try {
        # The optional viewer reads only known status JSON through the
        # pinned Python runtime. Its failure never interrupts the science.
        # Do not launch a separate unsigned PowerShell monitor process.
        try {
            $monitorArgs='"{0}" --out "{1}" --interval 60 --wait-for-start' -f $Monitor,$Out
            Start-Process -FilePath $Runtime -ArgumentList $monitorArgs -ErrorAction Stop | Out-Null
        } catch {
            Write-Host 'MONITOR_UNAVAILABLE_SCIENTIFIC_RUN_CONTINUES' -ForegroundColor Yellow
            Write-Host 'Progress will remain visible in this main console.'
        }
        & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $Runner -VcfPath $VcfPath -ClinicalPath $ClinicalPath -Out $Out
        if ($LASTEXITCODE -ne 0) {throw "RDMVA_JURY_SCIENTIFIC_EXECUTION_INCOMPLETE_EXIT_"+$LASTEXITCODE}
        $srcReceipt=Join-Path $Root 'public_inputs\current_v31\CURRENT_PUBLIC_DATA_V31.json'
        if (Test-Path $srcReceipt -PathType Leaf) {
            $s=Get-Content $srcReceipt -Raw | ConvertFrom-Json
            if ($s.state -eq 'DEGRADED') {
                Write-Host 'DOWNLOAD WARNING: Research completed with MISSING public datasets; results are limited.' -ForegroundColor Red
                return
            }
        }
        Write-Host 'RDMVA_JURY_SCIENTIFIC_AND_INDEPENDENT_VERIFICATION_PASS' -ForegroundColor Green
    } finally {
        Write-Host ("LOCAL_RUN_OUTPUT="+$Out)
    }
}
function Show-Monitor {
    if (!(Test-Path -LiteralPath $Runtime -PathType Leaf)) {
        throw 'Install the scientific runtime before opening the progress monitor.'
    }
    if (!$Out) {$script:Out=(Read-Host 'Enter the existing research output directory').Trim('"')}
    & $Runtime $Monitor --out $Out --interval 60
    if ($LASTEXITCODE -ne 0) {throw 'RDMVA_MONITOR_CANNOT_READ_RUN_STATUS'}
}
function Menu {
    while ($true) {
        Show-DownloadWarnings
        Write-Host ''
        Write-Host '=====================================' -ForegroundColor Cyan
        Write-Host ' RARE DISEASE MVA PHARMA V3.1 - JURY '
        Write-Host '=====================================' -ForegroundColor Cyan
        Write-Host ' 1  Install and verify research environment'
        Write-Host ' 2  Set private OpenAI API key'
        Write-Host ' 3  Test live OpenAI models (small billable synthetic requests)'
        Write-Host ' 4  Run full Challenge analysis'
        Write-Host ' 5  Monitor existing analysis'
        Write-Host ' 6  Check installation readiness'
        Write-Host ' 7  Run automatic regression tests'
        Write-Host ' 0  Exit'
        $choice=Read-Host 'Choose 0-7'
        try {
            switch($choice) {
                '1' {Install; break}
                '2' {Configure-Key; break}
                '3' {Test-LiveModels; break}
                '4' { try { Start-Research } finally { $script:Out=''; $script:VcfPath=''; $script:ClinicalPath='' }; break }
                '5' {Show-Monitor; break}
                '6' {$null=Check-Prerequisites; break}
                '7' {Run-SelfTests; break}
                '0' {return}
                default {Write-Host 'Choose an option 0-7'}
            }
        } catch {
            Write-Host ('OPERATION_FAILED='+$_.Exception.Message) -ForegroundColor Red
        }
    }
}
switch($Action) {
    'check' {$ok=Check-Prerequisites;if (!$ok){exit 2}}
    'install' {Install}
    'key' {Configure-Key}
    'run' {Start-Research}
    'monitor' {Show-Monitor}
    'test' {Run-SelfTests}
    'models' {Test-LiveModels}
    'menu' {Menu}
}
