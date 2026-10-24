# Rare Disease MVA Pharma — V3.1 Challenge Edition

**ExaMedix Technologies** | Developed by Matteo Caprani | [examedix.ch](https://www.examedix.ch)

**Rare Disease, Real Kid: MVA Hackathon 2026 — Track 2, Drug Repurposing**

## Track 2 submission: pitch video and scientific report

- **[Watch the 3-minute scientific pitch on YouTube](https://youtu.be/_r-avURNnsY)** — English narration, unlisted and viewable by anyone with the link.
- **[Video pitch information for the judging panel](submission/VIDEO_PITCH.md)** — title, proposed medicines, and review links.
- **[Scientific report (PDF)](submission/examedix-technologies_track2_report.pdf)** and **[report in Markdown](submission/examedix-technologies_track2_report.md)** — detailed mechanism-based proposals, supporting literature, limitations, and validation plan.

*The official Track 2 submission form requires the YouTube video URL and the repository URL in separate fields.*


Rare Disease MVA Pharma is a patient-first research program that interprets authorized genomic and clinical inputs, builds biological rescue hypotheses, screens approved medicines and checks supporting and contradictory evidence. It runs locally through a guided Windows menu; an external ChatGPT session is not needed to operate it.

> **Research only.** The output is not a diagnosis, prescription, dose recommendation, clinical efficacy estimate or safety clearance. A successful software run can legitimately conclude **NO_VALIDATED_LEAD**.

## Contents

1. [System requirements and what you need](#1-system-requirements-and-what-you-need)
2. [Download the repository and extract the ZIP](#2-download-the-repository-and-extract-the-zip)
3. [Install the research environment step by step](#3-install-the-research-environment-step-by-step)
4. [Configure and verify OpenAI access](#4-configure-and-verify-openai-access)
5. [Check readiness and start the full analysis](#5-check-readiness-and-start-the-full-analysis)
6. [Verify completion and read the results](#6-verify-completion-and-read-the-results)
7. [Troubleshooting and safe recovery](#7-troubleshooting-and-safe-recovery)

**Guided menu order:** `1 Install → 2 API key → 3 Model checks → 6 Readiness → 7 Self-tests → 4 Full analysis → 0 Exit after completion`.

Menu numbers are choices to enter **one at a time**, not a single command. Each operation must finish before you select the next one. Option `5` is a read-only progress monitor, not another analysis.

---

## 1. System requirements and what you need

### 1.1 Computer and access

| Requirement | What to prepare |
|---|---|
| Operating system | **Windows 11 x64**. This guide describes the Windows/WSL workflow. |
| Windows account | Administrator permission in the **same Windows user account** that will configure the API key and run the research. |
| Package installation | Microsoft **App Installer / WinGet** must be available. The installer checks it before proceeding. |
| Virtualization | Hardware virtualization and **WSL2** support. A Windows virtual machine also needs nested virtualization exposed by its host. |
| Scientific storage | An existing, writable, formatted **NTFS or exFAT data volume with at least 340 GiB free at initial selection**. Leave additional headroom for environment installation and outputs. |
| Internet | Access to official software repositories, biomedical data providers and the configured OpenAI API models. |
| Authorized case | The organizer-authorized **GRCh38 compressed VCF (`.vcf.gz`)** and its corresponding clinical phenotype **Word document (`.docx`)**. |
| API account | Your own **billable OpenAI API key**, with access to the models configured in the release. |

The 340 GiB storage check comprises a configured **240 GiB gnomAD resource budget plus a 100 GiB protected free-space reserve**. It is an admission requirement, not a promise that every run will download 240 GiB or occupy exactly 340 GiB. Actual case-specific requirements vary, and additional space checks apply during execution.

The release does not specify a validated numeric minimum for RAM, CPU cores or GPU hardware. Do not infer such a requirement from a particular development workstation. Storage speed, available memory, public servers, network throughput and API access all affect execution time; no fixed completion time is guaranteed.

### 1.2 What the installer supplies

You do **not** need to install the project's development tools manually before downloading the ZIP. Option `1` installs or verifies:

- Python **3.12**, patch **3.12.10 or later**, and Git **2.55.0 or later**, as checked by the bundled installer;
- the Windows Subsystem for Linux engine and a **project-specific Ubuntu 24.04 WSL2** distribution;
- the isolated Python research environment and hash-pinned dependencies;
- the initial public scientific references, local ChEMBL database and scientific-storage configuration.

WSL is the Linux environment used by the genomic tools. Its research account is created automatically with limited privileges. The reviewer is not asked for a Linux username or password: **no Linux username/password prompt or shell exit is required** for the guided setup.

### 1.3 What you must provide separately

The repository does **not** contain the participant's genome, clinical document or an API credential. Obtain case access through the organizer's authorized process. Keep the VCF and matching DOCX in a private directory **outside the extracted software folder**.

An illustrative layout is:

```text
E:\ResearchSoftware\rare-disease-mva-pharma-main\   Software and generated environment
E:\AuthorizedCase\case.vcf.gz                      Authorized genome
E:\AuthorizedCase\clinical.docx                    Matching clinical document
```

These are examples only. Your drive letter and filenames may differ. Keep the `.vcf.gz` compressed, keep the clinical document in `.docx` format, and use the matching genome reference build. Renaming a different file format does not convert it.

API calls in options `3` and `4` can incur charges. A ChatGPT subscription is not the API credential used by this application. No exact cost estimate is specified by the release; usage depends on the configured models and analysis.

## 2. Download the repository and extract the ZIP

1. Open [the GitHub repository](https://github.com/mattc789/rare-disease-mva-pharma-track2). If access is private, sign in with an authorized GitHub account.
2. Select **Code → Download ZIP**. Downloading the ZIP through the browser does not require Git to be installed.
3. Choose the final installation location on your intended research drive, checking its available capacity first.
4. Extract **all** ZIP contents there. Do not launch the program from inside the compressed ZIP window.
5. Open the extracted directory containing `JURY_START_V31.cmd`, `README.md`, `INSTALL_WINDOWS_V31.ps1`, `code` and `config`. This directory is the **project root**. Depending on extraction settings, it may be nested one level below the directory you selected.
6. Leave the authorized patient files in their separate private directory.

**Do not rename or move the project after installation starts.** The isolated WSL registration and generated local storage settings depend on the project location. Do not copy another computer's `env` directory, credentials, runtime status files or `config/local_*.json` into a new installation.

The ZIP is a **source package**, not a preinstalled environment. Large scientific databases are acquired separately. The attributed, hash-verified SIGNOR official-release seed is the documented exception, explained below.

## 3. Install the research environment step by step

### 3.1 Open the guided launcher

In the extracted project root, right-click **`JURY_START_V31.cmd` → Run as administrator**. On Windows 11 the command may appear under **Show more options**. Approve the Windows User Account Control prompt only for the launcher you deliberately opened.

Use your normal research Windows account with elevated permissions, **not a different administrator identity**. WSL registration and the encrypted API credential are user-specific.

The console menu is:

```text
1  Install and verify research environment
2  Set private OpenAI API key
3  Test live OpenAI models (small billable synthetic requests)
4  Run full Challenge analysis
5  Monitor existing analysis
6  Check installation readiness
7  Run automatic regression tests
0  Exit
Choose 0-7:
```

Type `1` and press **Enter**. Do not select the analysis option yet.

### 3.2 Verify the chosen scientific-data drive

Early in setup, look for these messages:

```text
RDMVA_STORAGE_AUTO_SELECTION_START
V31_STORAGE_AUTOMATIC_SELECTED=<selected scientific-data folder>
V31_STORAGE_LAYOUT=<selected layout>
V31_STORAGE_FREE_GIB=<measured free space>
V31_STORAGE_RESERVED_FREE_GIB=100
V31_STORAGE_NO_DISK_FORMATTING_PERFORMED=YES
```

The installer **prefers the volume containing the project** when that volume qualifies. With `SINGLE_PROJECT_ROOT`, public references use `<project root>\public_inputs\third_party` and the project keeps its other working folders under the same root.

If that volume does not qualify, the installer can select another eligible volume and report `SEPARATE_HIGH_CAPACITY_DATA_VOLUME`. That layout redirects the configured large genomic-reference storage; it does not relocate the entire project or its existing files. **Read the selected path rather than assuming a particular drive was used.**

The installer never formats, repartitions or erases disks. No eligible volume produces `V31_STORAGE_REQUIRED` instead of an incomplete installation.

<details>
<summary>Optional: explicitly choose a data drive before the first installation</summary>

Open **Windows PowerShell as administrator**, using the same Windows user. Change into the actual extracted project root and start the menu with `-DataDrive`:

```powershell
Set-Location -LiteralPath 'E:\ResearchSoftware\rare-disease-mva-pharma-main'
powershell.exe -NoProfile -STA -ExecutionPolicy Bypass -File .\JURY_START_V31.ps1 -Action menu -DataDrive E
```

Replace the example path and drive letter. Then select `1`. This selection does **not** bypass the space requirement. A conflicting existing storage configuration is rejected rather than silently moved or overwritten.

</details>

### 3.3 Let the installer complete its stages

Setup checks Windows and installs or reuses compatible Python/Git versions, initializes WSL2 Ubuntu, installs the research packages and acquires scientific references. A project-specific WSL name has the form `RDMVA-Ubuntu-24.04-<project-path hash>`; there is no fixed name to enter manually.

If it reports `REBOOT_REQUIRED`, save your work and restart Windows. Afterward reopen **the same launcher in the same folder and Windows account** and select `1` again. If it reports `REOPEN_SHELL`, close and reopen the administrator launcher, then select `1` again. Existing verified installation work is preserved; these messages are not a completed install.

### 3.4 Understand which datasets arrive when

| Phase | Scientific work |
|---|---|
| **Option 1: initial installation** | Public regulatory and annotation sources, including FDA/international records, ClinVar, GENCODE, HPO, SIGNOR and the local ChEMBL 37 SQLite database. |
| **Option 4: after genomic intake** | Preparation and verification of required genomic resources, including offline VEP/GRCh38 and the **gnomAD exome/genome chromosome files selected from the initial case analysis**, with their indexes. |

**gnomAD is deliberately not downloaded in full before interpreting the case.** The later selection can require tens of gigabytes or more. Successful installation alone does not certify the case's VEP or gnomAD evidence; the analysis checks those separately.

ChEMBL's official compressed archive is approximately **5.76 GB**, expanding into a database of approximately **30.48 GB**. Download, SHA-256 verification, extraction and database checks are different stages. A download or extraction message reaching 100% does **not** mean the whole installation is finished.

Large scientific files use parallel HTTPS byte-range transfers where the provider supports trustworthy ranges and source identity checks. Other files use a single HTTPS stream. `segmented_https` and `single_https` are both normal modes; parallelism does not guarantee a particular speed. Existing public files are reused only under the applicable integrity checks.

SIGNOR can be reconstructed from the included attributed October 2026 official-release seed when its official upstream cannot be reached. The software verifies the release hashes and explicitly labels that fallback; it is not a claim that a new download occurred.

### Automatic recovery when an official public-data source is unavailable

For each public-science download, the installer attempts the official HTTPS source, then **three additional attempts at 120-second intervals**. If the provider remains unavailable, it tries an eligible, verified alternative. For example, Health Canada's DPD REST API has an automatic non-API fallback to Health Canada's own official drug/ingredient ZIP extracts; SIGNOR has a pinned, attributed official-release seed. No unverified third-party source is silently substituted.

If all available methods fail, the installer **records the missing dataset and continues with the next source**. The warning appears in red during the attempts, after acquisition and whenever the existing menu returns. The machine-readable record is **INSTALL_WARNINGS.json** in the installation folder, with official URLs and explicit SKIPPED_CONTINUED states. The published dataset manifest includes **only successfully downloaded and SHA-256-verified material**.

A successful *software installation* with missing databases is marked PASS_WITH_SKIPPED_PUBLIC_SOURCES, not full-data PASS. Research can attempt the available evidence on a best-effort basis, but the missing coverage remains explicit; a scientifically required reference can still prevent full validation. This policy applies to public reference data, **not** mandatory Windows/WSL components or data-integrity checks. No menu option was added.



<details>
<summary>Optional: adjust scientific download concurrency</summary>

The default is **four connections per eligible file**; the supported range is **1–8**. To set a value, launch from the same administrator PowerShell session:

```powershell
Set-Location -LiteralPath 'E:\ResearchSoftware\rare-disease-mva-pharma-main'
$env:RDMVA_DATA_CONNECTIONS = '4'
powershell.exe -NoProfile -STA -ExecutionPolicy Bypass -File .\JURY_START_V31.ps1 -Action menu
```

The setting affects scientific-data transfers, not installation of Windows components or numerical computation. More connections need not be faster. Do not disable integrity or TLS checks to obtain a speedup.

</details>

### 3.5 Confirm installation success

Proceed only when the console prints **both**:

```text
RDMVA_INSTALL_PASS: ...
RDMVA_JURY_INSTALL_PASS
```

The menu then returns. A returned menu by itself is **not** proof of success: read any `OPERATION_FAILED=...` or action-required message first. At this point the environment is installed, but no complete patient analysis has been certified. If an upstream scientific dataset remains unavailable, the existing menu prominently displays red warnings and the installer records PASS_WITH_SKIPPED_PUBLIC_SOURCES instead of a complete-data success.

## 4. Configure and verify OpenAI access

### 4.1 Option 2 — store or reuse your private credential

Type `2` and press **Enter**. For a new account, enter the API key at the **hidden input prompt** and confirm it. A successful save prints:

```text
V31_OPENAI_USER_DPAPI_KEY_STORED=YES
V31_API_KEY_PLAINTEXT_FILE_CREATED=NO
```

The encrypted credential is stored for the current Windows user at:

```text
%LOCALAPPDATA%\Examedix\RareDiseaseMVAPharma\openai_key.dpapi
```

If it already exists, the menu reports `V31_API_CREDENTIAL_ALREADY_PRESENT=YES` and preserves it. Presence does not prove the key works: option `3` performs that check. The setup helper does not silently overwrite an existing encrypted key.

**The program does not automatically read a .env file.** Its credential inputs are the process environment's `OPENAI_API_KEY` or the user-scoped encrypted credential; use the guided encrypted-key workflow for this guide. Do not put a real key into README examples, screenshots, source files, GitHub, support messages or videos. Do not switch Windows user accounts between configuration and execution.

The bundled secure helper is `SETUP_OPENAI_KEY_V31.ps1`; it is called by option `2`. No plaintext credential file is needed.

### 4.2 Option 3 — verify the configured models

Type `3`, press **Enter**, and type `YES` at the explicit consent prompt to authorize the **small billable synthetic API tests**. These checks do not use the patient's genome. Declining consent makes no test request and does not establish model readiness.

Wait for the model checks to finish and confirm:

```text
ALL_REQUIRED_OPENAI_MODELS_VERIFIED
```

The required models come from the released code/configuration; do not substitute a model merely to suppress a failure. The full-run launcher requires a **successful billable model check from the last seven days**. Repeat option `3` when requested because the saved check is stale, absent or failed.

If a check fails, resolve credential, billing, model-availability or network errors before starting the long research run. API keys and raw responses can be sensitive; share only redacted diagnostic messages.

## 5. Check readiness and start the full analysis

### 5.1 Option 6 — installation preflight

Type `6` and press **Enter**. Confirm these conditions:

```text
PINNED_PYTHON_RUNTIME=True
PUBLIC_REFERENCE_SNAPSHOT=READY
API_CREDENTIAL_PRESENT=YES ...
SCIENTIFIC_STORAGE=<configured location>
RDMVA_JURY_PREFLIGHT_RUNTIME_PREREQUISITES_PRESENT
```

`RDMVA_JURY_PREFLIGHT_INCOMPLETE` means at least one prerequisite needs attention. Read the individual lines instead of proceeding. This preflight is separate from the genomic evidence gate later in the analysis.

### 5.2 Option 7 — automated software tests

Type `7` and press **Enter**. The suite uses **synthetic fixtures**, not the authorized patient case. Let it finish and look for:

```text
OK
RDMVA_JURY_SELF_TEST_PASS
```

Some tests deliberately exercise errors or unavailable services. Judge the **final suite result**, not an isolated message from a synthetic scenario. A nonzero exit or `RDMVA_JURY_SELF_TEST_FAILED` requires investigation. These tests establish software behavior, not a valid clinical lead and not completion of the actual case analysis.

### 5.3 Option 4 — confirm the authorized genome and clinical document

Type `4` and press **Enter**. The launcher searches configured candidate locations for local inputs. Its discovery uses filenames and local duplicate SHA-256 checks; it does not upload the source files.

- **One unambiguous pair:** the console displays the VCF/DOCX paths and asks whether you are authorized to use them. Check both paths and type `YES` to accept.
- **Byte-identical copies:** local hashing can recognize duplicates as the same pair. Hashing a large VCF can take time.
- **Different possible pairs, missing files or declined suggestion:** select the intended `.vcf.gz` and `.docx` in the Windows file dialogs. These dialogs may also be used for files in separate directories. If a dialog is unavailable, enter the complete file path at the console prompt.

Automatic discovery is a convenience, not authorization or proof that any two files belong together. **You must verify that the selected genome and phenotype document describe the same authorized case.**

<details>
<summary>Optional: point discovery to your private input directory</summary>

Before opening the menu, set this in the same PowerShell session:

```powershell
Set-Location -LiteralPath 'E:\ResearchSoftware\rare-disease-mva-pharma-main'
$env:RDMVA_CHALLENGE_INPUT_DIR = 'E:\AuthorizedCase'
powershell.exe -NoProfile -STA -ExecutionPolicy Bypass -File .\JURY_START_V31.ps1 -Action menu
```

Other recognized candidate locations may still be checked. This setting does not authorize automatic guessing between different cases. The file picker remains available when discovery is ambiguous.

</details>

### 5.4 Accept a new results path and authorize the run

The program proposes a timestamped path:

```text
DEFAULT_OUTPUT=<project root>\experiments\JURY_V31_<timestamp>
Press Enter to use this NEW folder, or enter another NEW folder:
```

Press **Enter** to accept it. Alternatively enter a new full path. **The output directory must not already exist, even if it is empty; do not create it manually.** The program creates it. An existing path is rejected to protect earlier results.

Next, type `YES` at the scientific-analysis consent prompt. This starts the complete genomic and pharmacological workflow, including its billable online research. Save the displayed `RUN_OUTPUT=...` path for monitoring and final verification.

### 5.5 Keep the main research console open

The workflow proceeds through local patient intake, required genomic resources and evidence, patient-state mechanisms, rescue objectives, biological routes, approved-medicine screening, source review, safety/exposure and combination checks, and independent verification.

The main window may show:

```text
E2E_STAGE_1_START
V31_NATIVE_STAGE_START <stage>
V31_RESEARCH_STAGE_ACTIVITY <stage> elapsed_s=... percent=unavailable state=running
V31_NATIVE_STAGE_DONE <stage>
E2E_STAGE_1_PASS
E2E_STAGE_2_START
```

`percent=unavailable` means the program has no honest total-work percentage for that stage; it is not itself an error. A heartbeat or a saved `RUNNING` status alone is not proof that every underlying operation is advancing. Do not close the main console, change the program files, relocate its folders or switch source configurations during the run.

A separate read-only monitor may open automatically. Its failure does not stop the scientific process. Option `5` can display an existing output directory from another launcher instance; **do not start a second option-4 analysis just to inspect progress**.

<details>
<summary>Advanced: direct command-line execution and monitoring</summary>

The guided menu above is the normal path. All commands below run in Windows PowerShell from the extracted project root. Replace the example paths with your own. The output directory must be new.

```powershell
Set-Location -LiteralPath 'E:\ResearchSoftware\rare-disease-mva-pharma-main'

# Alternative to option 1; use an administrator shell.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\INSTALL_WINDOWS_V31.ps1 -CurrentPublicData

# Alternative to option 2, only when no encrypted key is already stored.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\SETUP_OPENAI_KEY_V31.ps1

# Perform model checks, readiness and self-tests through the guided menu.
powershell.exe -NoProfile -STA -ExecutionPolicy Bypass -File .\JURY_START_V31.ps1 -Action menu

# Alternative to option 4 after the required checks have passed.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\REPRODUCE_CHALLENGE_V31.ps1 -VcfPath 'E:\AuthorizedCase\case.vcf.gz' -ClinicalPath 'E:\AuthorizedCase\clinical.docx' -Out 'E:\AuthorizedResults\NEW_RUN'
```

To monitor that run from a **separate** console, without starting research again:

```powershell
Set-Location -LiteralPath 'E:\ResearchSoftware\rare-disease-mva-pharma-main'
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\JURY_START_V31.ps1 -Action monitor -Out 'E:\AuthorizedResults\NEW_RUN'
```

The direct research wrapper does not replace the guided readiness sequence. Do not use reduced-science, no-download or offline-AI switches as a substitute for the complete Challenge workflow described here.

</details>

## 6. Verify completion and read the results

### 6.1 Confirm the whole run, not just an individual stage

For a successful guided full run, the console reports:

```text
E2E_STAGE_1_PASS
E2E_STAGE_2_PASS
CLEANROOM_REPRODUCTION_PASS True
AI_RESPONSE_RECEIPTS <verified number>
ATTESTATION <output directory>\CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json
E2E_V31_REPRODUCTION_PASS
V31_COMPLETE_CHALLENGE_PASS
RDMVA_JURY_SCIENTIFIC_AND_INDEPENDENT_VERIFICATION_PASS
LOCAL_RUN_OUTPUT=<output directory>
```

The response count is **measured**, not a fixed target. The final overall markers matter: stage 1 alone, installation success or a passed test suite does not certify the full research run.

Open the directory printed after `LOCAL_RUN_OUTPUT=`. Then inspect the following files; JSON is readable text and can be opened in a text editor **without editing or saving changes**.

| File in the output directory | What to check |
|---|---|
| `E2E_EXECUTION_V31.json` | Overall `"status": "PASS"` and successful native, V3.1 research and independent-verification stages. Fresh/recovery provenance is recorded separately. |
| `GENOMIC_REQUIRED_EVIDENCE_V31.json` | `"required_science_pass": true`, actual offline VEP evidence and successful required gnomAD cohort/chromosome lookups. Case-focused coverage is not global whole-genome coverage. |
| `CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json` | `"passed": true` and successful individual checks. For the full AI-enabled profile, genuine API receipts and their `store=false` consistency must be verified. |
| `FINAL_RESULT_V31_NATIVE.json` | Actual screened medicine count, route statistics, strict research conclusion and pointers to the result artifacts. |
| `TRACK2_HYPOTHESES_V31_FINAL.json` | The broader `single_agent_hypotheses` list and the narrower `submission_level_single_agents` list, including evidence states, blockers and combination results. |
| `RUN_REPRODUCIBILITY_CONTEXT_V31.json` | Observations and source-coverage context from this run, without an imposed expected medicine. |
| `challenge_submission_v31/TRACK2_RESEARCH_DOSSIER_V31.md` | The generated, human-readable hypothesis dossier. Read its limitations as well as its candidate names. |
| `STATUS_V30_NATIVE.json` | Detailed native-stage completion/progress information, useful when diagnosing an incomplete run. |

A receipt's name is not proof of an independent virgin operating system. Use the recorded source/input/stage provenance to describe what was actually executed, and do not represent a recovered analysis as a new installation.

### 6.2 Interpret the scientific outcome correctly

**Technical PASS means completed execution of the implemented checks, not clinical validation.** A study in another disease, species or cellular context is not proof of MVA rescue or of correction of the participant's genotype.

An exploratory candidate is not an approved treatment for this case. `HOLD` means a relevant limitation or concern remains. `NOT_ASSESSED`, `NOT_REVIEWED` and `NOT_ATTACKED` indicate that the corresponding assessment was not performed, not that it was passed. `NO_VALIDATED_LEAD` is a legitimate scientific conclusion, not a program crash.

The numerical research-priority scores are **not percentages of efficacy**. Missing evidence is not proof of ineffectiveness, and a publication or a chemical/network association is not proof of clinical benefit. A fresh run is not required to reproduce a particular medicine, ranking or count: official source snapshots, source availability and model outputs can differ. Preserve the actual run's provenance.

For the biological method, disclosure of AI tools, candidate evidence and limitations, see the [final Track 2 scientific report](submission/examedix-technologies_track2_report.md), which includes the official methods-template answers. The report describes the completed research execution; it is not a preselected answer fed to the program.

### 6.3 Exit and keep the private results

When the analysis has finished and the menu has returned with the final success markers, type `0` and press **Enter**. This exits the launcher; it does not delete the saved results. Selecting `4` again starts a new analysis rather than opening the completed report.

Keep the authorized inputs, full output directory and source/data provenance available under the organizer's privacy rules. **Do not upload the entire run folder to GitHub.** It can contain genotype/clinical information, local paths, source excerpts and API metadata that are not suitable for public redistribution. Use only reviewed, appropriately redacted submission materials.

The repository contains the program, installation/reproducibility tools, software tests, licensing information, and the [complete scientific report with integrated methods disclosure in Markdown](submission/examedix-technologies_track2_report.md) and [PDF](submission/examedix-technologies_track2_report.pdf). The editable Word version is kept privately, not published to GitHub. The separate three-minute video URL belongs in the Challenge submission process; video scripts and production assets are not part of this repository.

## 7. Troubleshooting and safe recovery

### 7.1 Find the first relevant error

Returning to the menu does not necessarily mean an operation passed. Read the message immediately preceding it, especially `OPERATION_FAILED=...`, and preserve the existing files. Do not repeatedly reinstall, delete databases or rerun paid research without first identifying the failed stage.

| Symptom or message | Meaning and next action |
|---|---|
| Launcher says to reopen as administrator | Close it, then use **Run as administrator in the same Windows account**. Do not use a different user's API/WSL environment. |
| WinGet is missing | Install or restore Microsoft App Installer/WinGet through the appropriate Windows administration process, then reopen the launcher and retry option `1`. |
| `REOPEN_SHELL` after Python/Git installation | Reopen the administrator launcher in the same project directory; select `1` again to verify the installed version. |
| `REBOOT_REQUIRED` | Save work, restart Windows, reopen the same launcher and select `1` again. Do not proceed to analysis before installation PASS. |
| `WSL_UPDATE_FAILED` or `WSL_SETUP_REQUIRED` | Read `INSTALL_WINDOWS_V31_STATUS.json`; check the reported Windows/WSL prerequisite, connectivity and ownership/location issue. A VM host may need to enable nested virtualization. Do not unregister unrelated WSL distributions. |
| `V31_STORAGE_REQUIRED` / `V31_STORAGE_NO_ELIGIBLE_VOLUME` | Initial storage selection did not find enough usable capacity. Provide an eligible formatted volume with at least 340 GiB free, or select one explicitly before setup. No disk is automatically formatted. |
| `V31_STORAGE_EXISTING_SELECTION_CONFLICT`, changed drive or reserve error | The saved storage configuration no longer matches the requested drive/device or protected capacity. Restore the intended available storage or investigate the configuration; do not alter the JSON simply to bypass the check. |
| `PUBLIC_REFERENCE_SNAPSHOT` is not `READY` | Installation has not established a complete verified snapshot. Inspect the public-data receipt and rerun option `1` after resolving the reported source issue. |
| Upstream provider temporarily unavailable | Keep verified files and retry the affected installation after source connectivity is restored. The program does not approve an incomplete snapshot. |
| Download/extraction reaches 100%, but no installation PASS | Hashing, extraction, database validation or another setup stage may still be running. Use the current stage message; a per-file percentage is not global completion. |
| SHA-256, official-size, archive or index check fails | The file is incomplete, damaged, changed or inconsistent with its source. Do not replace expected hashes or force a PASS. Retain the diagnostic and resolve the source/download problem. |
| API credential exists, but option `3` fails | Credential presence alone does not validate billing, key access or model availability. Check the actual error in the same Windows account. The key helper intentionally will not overwrite an existing encrypted key. |
| Model readiness not verified in the last seven days | Run option `3` again and obtain `ALL_REQUIRED_OPENAI_MODELS_VERIFIED` before option `4`. |
| No input pair found or multiple input pairs found | Use the file selectors to choose the authorized matching `.vcf.gz` and `.docx`. This is a safety control, not permission to use another case. |
| `RDMVA_REQUIRES_COMPRESSED_VCF_GRCH38` / `RDMVA_REQUIRES_CLINICAL_DOCX` | Supply the original required format/build. Do not solve the problem by renaming a file extension. |
| `Output directory already exists; choose a NEW directory.` | Choose a nonexistent output path or accept the next timestamped default. Do not pre-create the directory or delete completed results to reuse its name. |
| `MONITOR_UNAVAILABLE_SCIENTIFIC_RUN_CONTINUES` | The optional viewer failed to open; use the main scientific console or option `5` from another instance. Do not restart the science merely to restore the viewer. |
| Heartbeat continues with no percentage | The stage lacks a reliable total-work estimate. Check recent stage/output activity; a heartbeat alone cannot rule out a slow or blocked dependency. |
| `RDMVA_JURY_SELF_TEST_FAILED` or full analysis exit error | Inspect the final test summary or the scientific execution/failure receipts. Do not label the run complete until the appropriate independent checks pass. |

Do not disable Windows security, TLS validation, source hashes or mandatory genomic gates to make setup appear successful.

### 7.2 Inspect the appropriate local status files

For installation problems, use `INSTALL_WINDOWS_V31_STATUS.json`, `INSTALL_V2810_STATUS.json` and `public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json` under the project root. For research problems, use `E2E_EXECUTION_V31.json`, `STATUS_V30_NATIVE.json` and any failure/genomic-evidence files in the **specific output directory** printed by the launcher.

Report the software commit/version, the failing menu option, the exact redacted error and relevant status fields. Never publish raw VCF/DOCX contents, genomic coordinates, complete private logs or API secrets. The project contact/presentation page is [examedix.ch](https://www.examedix.ch).

### 7.3 Interrupted installation versus interrupted research

**Installation:** rerun option `1` in the same location after resolving the cause. The setup revalidates existing public inputs and supports applicable download checkpoints. Checkpoint reuse depends on a consistent source identity and integrity; it is not a blanket promise that every partial file is reusable.

**Research:** retain the incomplete output and original inputs. The program includes explicit recovery paths guarded by matching input hashes and source provenance; a recovered run must remain labeled as such. For an independent fresh run, use the original authorized VCF/DOCX and a **new output directory**. This repeats scientific work and may repeat billable API requests. Do not change receipt status fields or use scientific-skip flags to convert an incomplete run into a PASS.

Do not replace or move the source tree, WSL environment or data-storage configuration during execution. Do not combine outputs from different source profiles or patients.

### 7.4 Privacy, licenses and source attribution

The runtime uses the OpenAI Responses API with `store=false` for configured online tasks. Its data-handling design excludes raw VCF, exact genomic coordinates/HGVS, genotypes and identifying clinical records from model payloads. These are source-level/program declarations, **not an independent network-traffic audit**. ChatGPT Pro was used externally during engineering/research review; it is not required to orchestrate execution. Full disclosure is in the [scientific report and integrated methods section](submission/examedix-technologies_track2_report.md).

Original code and documentation are distributed under the project's [LICENSE](LICENSE). External datasets, packages and pretrained resources retain **their own licenses**. ChEMBL remains external CC BY-SA data; the included SIGNOR official-release seed retains its own attribution. See [THIRD_PARTY_NOTICES_V31.md](THIRD_PARTY_NOTICES_V31.md), [SIGNOR attribution](reference_seed/signor/ATTRIBUTION.md) and [source policy](config/v30_source_policy.json). Gated patient data remain subject to the organizer's access agreement and are not redistributed.

**No result from this software authorizes a medication, dosage or treatment for the participant.**
