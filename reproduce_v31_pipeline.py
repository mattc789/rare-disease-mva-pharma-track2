#!/usr/bin/env python3
"""Autonomous two-stage V3.1 runner: released source + authorized local input.

No ChatGPT orchestration at runtime. The API key is read from the environment
by the original V3 modules and is never written to this runner's logs.
"""
from __future__ import annotations
import argparse, datetime as dt, hashlib, json, os, shutil, subprocess, sys, traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE = ROOT / "code"

from contextlib import contextmanager
import threading,time

@contextmanager
def long_stage_activity(label,interval_seconds=15):
    """Console heartbeat for stages without a truthful total-work percent."""
    started=time.monotonic()
    completed=threading.Event()
    def emit():
        while not completed.wait(interval_seconds):
            print(f"V31_RESEARCH_STAGE_ACTIVITY {label} "
                  f"elapsed_s={int(time.monotonic()-started)} "
                  "percent=unavailable state=running",flush=True)
    watcher=threading.Thread(target=emit,daemon=True)
    watcher.start()
    try:
        yield
    finally:
        completed.set()
        watcher.join(timeout=2)

def write_json(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temp.replace(path)

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def verify_failed_intake_for_resume(out: Path, src: Path, clinical: Path | None,
                                    previous: dict) -> Path:
    # This is a scientific *resume*, not a new clean-room run. Reuse only the
    # software-generated intake of the SAME original, SHA-256-confirmed VCF,
    # phenotype and frozen/updated public snapshot. Never trust a supplied
    # patient case from a different input or a half-written intake directory.
    if previous.get("status")!="FAILED" or not previous.get("stages",{}).get("stage1_native_core"):
        raise RuntimeError("V31_RESUME_REQUIRES_FAILED_STAGE1")
    if previous.get("input_file_type")!="authorized_vcf" or previous.get("input_file_sha256")!=sha256(src):
        raise RuntimeError("V31_RESUME_VCF_SOURCE_MISMATCH")
    if not clinical:
        raise RuntimeError("V31_RESUME_REQUIRES_CLINICAL_DOCX")
    if previous.get("public_data_profile") != (
            "CURRENT_OFFICIAL_SOURCES" if os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1"
            else "FROZEN_CHALLENGE_BASELINE"):
        raise RuntimeError("V31_RESUME_PUBLIC_PROFILE_MISMATCH")
    live=ROOT/"public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json"
    if os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1" and (
            not live.is_file() or
            previous.get("current_public_data_receipt_sha256")!=sha256(live)):
        raise RuntimeError("V31_RESUME_PUBLIC_SNAPSHOT_CHANGED")
    intake=out/"intake"
    case=intake/"case_model_v2.json"
    manifest=intake/"intake_manifest.json"
    if not case.is_file() or not manifest.is_file():
        raise RuntimeError("V31_RESUME_INTAKE_FILES_MISSING")
    original=json.loads(case.read_text(encoding="utf-8"))
    metadata=json.loads(manifest.read_text(encoding="utf-8"))
    expected_inputs={src.name:sha256(src),clinical.name:sha256(clinical)}
    original_hashes={x["file"]:x["sha256"] for x in original.get("official_source_hashes",[])}
    manifest_hashes={x["file"]:x["sha256"] for x in metadata.get("inputs",[])}
    if expected_inputs != original_hashes or expected_inputs != manifest_hashes:
        raise RuntimeError("V31_RESUME_INTAKE_INPUT_PROVENANCE_MISMATCH")
    outputs=metadata.get("outputs") or {}
    if "case_model_v2.json" not in outputs or not isinstance(outputs,dict):
        raise RuntimeError("V31_RESUME_INTAKE_OUTPUT_MANIFEST_INVALID")
    for filename,checksum in outputs.items():
        p=intake/filename
        if (not p.is_file() or p.name!=filename or
                p.parent!=intake or sha256(p)!=checksum):
            raise RuntimeError("V31_RESUME_INTAKE_OUTPUT_SHA256_MISMATCH:"+filename)
    return case


def prepare_recovery_from_interrupted_run(
        previous: Path, destination: Path, src: Path, clinical: Path | None) -> dict:
    # Preserve the original attempt without changing its incomplete receipt.
    # Recover ONLY a verified software-generated intake, not a supplied case
    # whose scientific origin has not been independently established.
    if not clinical:raise RuntimeError("V31_RECOVER_REQUIRES_CLINICAL_DOCX")
    if previous.resolve()==destination.resolve():
        raise RuntimeError("V31_RECOVER_REQUIRES_NEW_OUTPUT_DIRECTORY")
    if destination.exists():
        raise RuntimeError("V31_RECOVER_REFUSES_EXISTING_OUTPUT")
    old_file=previous/"E2E_EXECUTION_V31.json"
    archive=previous/"E2E_EXECUTION_V31_PREVIOUS_FAILED.json"
    if not old_file.is_file() or not archive.is_file():
        raise RuntimeError("V31_RECOVER_PREVIOUS_RUN_RECEIPTS_MISSING")
    current=json.loads(old_file.read_text(encoding="utf-8"))
    initial=json.loads(archive.read_text(encoding="utf-8"))
    if (current.get("status") not in ("RUNNING","FAILED") or
        current.get("resumed_with_verified_original_intake") is not True):
        raise RuntimeError("V31_RECOVER_EXPECTS_INTERRUPTED_RESUMED_RUN")
    if (current.get("input_file_sha256")!=sha256(src) or
        current.get("current_public_data_receipt_sha256")!=
            initial.get("current_public_data_receipt_sha256")):
        raise RuntimeError("V31_RECOVER_INPUT_OR_SNAPSHOT_CHANGED")
    prior_status_file=previous/"STATUS_V30_NATIVE.json"
    if not prior_status_file.is_file():
        raise RuntimeError("V31_RECOVER_PREVIOUS_STAGE_LOG_MISSING")
    previous_stage=json.loads(prior_status_file.read_text(encoding="utf-8"))
    if previous_stage.get("current_stage") not in (
            "GENOMIC_EVIDENCE","V30_UNIVERSAL_SCREEN","APPROVED_DRUG_UNIVERSE",
            "CAUSAL_GRAPH_AND_COMBINATIONS","V31_NATIVE_ROUTE_GENERATION",
            "V31_SOURCE_BOUND_ATTACK","SOURCE_BOUND_ATTACK"):
        raise RuntimeError("V31_RECOVER_UNEXPECTED_PREVIOUS_STAGE")
    case=verify_failed_intake_for_resume(previous,src,clinical,initial)
    # No new scientific process may be running against the original output.
    if os.name=="nt":
        ps=(
          "Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
          "Select-Object ProcessId,CommandLine | ConvertTo-Json -Compress")
        check=subprocess.run(["powershell.exe","-NoProfile","-NonInteractive",
                              "-Command",ps],capture_output=True,text=True,
                              errors="replace",timeout=45)
        if check.returncode:
            raise RuntimeError("V31_RECOVER_CANNOT_VERIFY_OLD_RUN_STOPPED")
        processes=json.loads(check.stdout) if check.stdout.strip() else []
        if isinstance(processes,dict):processes=[processes]
        for proc in processes:
            cmd=str(proc.get("CommandLine") or "").lower()
            if (str(previous).lower() in cmd and
                ("v30.pipeline" in cmd or
                 ("reproduce_v31_pipeline.py" in cmd and "--recover-from" not in cmd))):
                raise RuntimeError("V31_RECOVER_ORIGINAL_PROCESS_STILL_ACTIVE")
    old_intake=case.parent
    metadata=json.loads((old_intake/"intake_manifest.json").read_text(encoding="utf-8"))
    source_outputs=metadata["outputs"]
    dest_intake=destination/"intake"
    dest_intake.mkdir(parents=True,exist_ok=False)
    paths={"intake_manifest.json":old_intake/"intake_manifest.json"}
    for name,checksum in source_outputs.items():
        f=old_intake/name
        if f.parent!=old_intake or not f.is_file() or sha256(f)!=checksum:
            raise RuntimeError("V31_RECOVER_UNVERIFIED_INTAKE_FILE")
        paths[name]=f
    for name,f in paths.items():
        dst=dest_intake/name
        shutil.copy2(f,dst)
        if sha256(dst)!=sha256(f):
            raise RuntimeError("V31_RECOVER_COPY_SHA256_MISMATCH:"+name)
    output={
      "schema":"RDMVA_RECOVERY_PROVENANCE_V31_1",
      "source_run":str(previous),"source_status":current["status"],
      "source_stage":previous_stage["current_stage"],
      "source_initial_failed_receipt_sha256":sha256(archive),
      "source_interrupted_receipt_sha256":sha256(old_file),
      "authorized_input_vcf_sha256":sha256(src),
      "authorized_clinical_docx_sha256":sha256(clinical),
      "intake_files_verified_and_copied":len(paths),
      "reused_genomic_annotations":False,
      "reuse_of_chembl_target_cache":"only_after_source_master_sha256_check",
      "fresh_end_to_end":False
    }
    write_json(destination/"RECOVERY_PROVENANCE_V31.json",output)
    return output

def main() -> int:
    p=argparse.ArgumentParser(description="Rare Disease MVA Pharma V3.1 full Challenge runner")
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--case")
    g.add_argument("--vcf")
    p.add_argument("--clinical")
    p.add_argument("--out",required=True)
    p.add_argument("--resume-failed-run",action="store_true",
                   help="Resume an existing failed stage-1 run using its SHA256-verified, software-generated VCF intake.")
    p.add_argument("--recover-from",help="Recovery-only: reuse verified intake from prior interrupted run into a NEW output.")
    p.add_argument("--no-online-ai",action="store_true",help="Diagnostics only: not the Challenge reproduction")
    p.add_argument("--no-resource-download",action="store_true",help="Use verified pinned public snapshots already installed")
    p.add_argument("--literature-attack-count",type=int,default=20)
    p.add_argument("--skip-full-scientific-modules",action="store_true",help="Diagnostics only")
    a=p.parse_args()
    if a.clinical and not a.vcf:p.error("--clinical requires --vcf")
    if not 0 <= a.literature_attack_count <= 100:p.error("invalid literature attack count")
    if not a.no_online_ai and not os.environ.get("OPENAI_API_KEY"):
        # Optional Windows-user DPAPI credential for this local QA machine.
        # Public reviewers should normally set OPENAI_API_KEY in their environment.
        sys.path.insert(0,str(CODE))
        try:
            from v2810.ai import api_key
            api_key()   # Check usability without logging or copying the secret
        except (OSError,ValueError,RuntimeError) as exc:
            p.error("Missing/unavailable OPENAI_API_KEY or local Windows DPAPI credential.")
    src=Path(a.case or a.vcf).resolve(strict=True)
    clinical=Path(a.clinical).resolve(strict=True) if a.clinical else None
    out=Path(a.out).resolve()
    reused_intake=None
    recovery_provenance=None
    if a.recover_from and a.resume_failed_run:
        p.error("Choose only --recover-from or --resume-failed-run")
    if a.recover_from:
        if a.case:p.error("--recover-from requires original VCF, not --case")
        previous=Path(a.recover_from).resolve(strict=True)
        recovery_provenance=prepare_recovery_from_interrupted_run(
            previous,out,src,clinical)
        reused_intake=out/"intake/case_model_v2.json"
        # Transfer only public ChEMBL target-cache objects; route-dependent
        # reuse is individually gated when the new master ledger is built.
        os.environ["RDMVA_RECOVERY_CHEMBL_FROM"]=str(previous)
    elif out.exists():
        if not a.resume_failed_run:
            p.error("Refusing to overwrite an existing run directory")
        receipt_file=out/"E2E_EXECUTION_V31.json"
        if not receipt_file.is_file():
            p.error("Cannot resume without an existing original execution receipt")
        prev=json.loads(receipt_file.read_text(encoding="utf-8"))
        if (prev.get("online_ai_requested") is not (not a.no_online_ai) or
                prev.get("no_resource_download") is not a.no_resource_download or
                prev.get("full_scientific_modules") is not (not a.skip_full_scientific_modules) or
                prev.get("literature_attack_count")!=a.literature_attack_count):
            p.error("V31_RESUME_ORIGINAL_FLAGS_DIFFER")
        reused_intake=verify_failed_intake_for_resume(out,src,clinical,prev)
        archive=out/"E2E_EXECUTION_V31_PREVIOUS_FAILED.json"
        if archive.exists():p.error("V31_RESUME_ALREADY_ATTEMPTED")
        archive.write_text(json.dumps(prev,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    elif a.resume_failed_run:
        p.error("Cannot resume a nonexistent run directory")
    else:
        out.mkdir(parents=True)
    # One common V3 AI cost ledger for all agent/submodule calls in this run.
    # Legacy V2.8 calls additionally retain their own built-in global cap.
    os.environ['RDMVA_E2E_GLOBAL_AI_BUDGET_ROOT']=str(out/'api_budget_global')
    started=dt.datetime.now(dt.timezone.utc)
    receipt={
        "schema":"RDMVA_E2E_EXECUTION_V31_1",
        "started_utc":started.isoformat(),
        "version":"3.1",
        "source_repository":str(ROOT),
        "input_file_sha256":sha256(src),
        "input_file_type":"prebuilt_case" if a.case else "authorized_vcf",
        "resumed_with_verified_original_intake":bool(reused_intake),
        "fresh_clean_room_run":not bool(reused_intake),
        "recovered_from_interrupted_run":bool(recovery_provenance),
        "recovery_provenance_file":"RECOVERY_PROVENANCE_V31.json" if recovery_provenance else None,
        "online_ai_requested":not a.no_online_ai,
        "no_resource_download":a.no_resource_download,
        "full_scientific_modules":not a.skip_full_scientific_modules,
        "literature_attack_count":a.literature_attack_count,
        "status":"RUNNING",
        "public_data_profile":("CURRENT_OFFICIAL_SOURCES"
             if os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1"
             else "FROZEN_CHALLENGE_BASELINE"),
        "current_public_data_manifest":"public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json"
             if os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1" else None,
        "current_public_data_receipt_sha256":(
            sha256(ROOT/"public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json")
            if os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1"
            and (ROOT/"public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json").is_file()
            else None),
        "stages":{},
        "api_key_included":False
    }
    report=out/"E2E_EXECUTION_V31.json"
    write_json(report,receipt)
    try:
        cmd=[sys.executable,"-m","v30.pipeline"]
        # Optional user-supplied *public* gnomAD manifest, never bundled.
        # Fresh jury installations download the official case-required shards.
        # The scientific stage independently verifies source+index SHA-256
        # and runs real local index queries before accepting such a manifest.
        manifest_setting=os.environ.get("RDMVA_LOCAL_GNOMAD_MANIFEST","").strip()
        local_manifest=(Path(manifest_setting).expanduser() if manifest_setting
                        else ROOT/"config/local_public_gnomad_manifest.json")
        if manifest_setting or local_manifest.exists():
            if not local_manifest.is_file():
                raise RuntimeError("V31_SUPPLIED_GNOMAD_MANIFEST_NOT_FOUND")
            sys.path.insert(0,str(CODE))
            from v2810.gnomad_local import validate_manifest
            supplied=json.loads(local_manifest.read_text(encoding="utf-8"))
            shards,coverage=validate_manifest(local_manifest)
            if not shards or coverage.get("invalid_or_missing_shards"):
                raise RuntimeError("V31_SUPPLIED_GNOMAD_MANIFEST_INVALID_OR_INCOMPLETE")
            for entry in supplied.get("shards",[]):
                for pin in ("source_sha256","index_sha256"):
                    value=str(entry.get(pin) or "")
                    if len(value)!=64 or any(c not in "0123456789abcdef" for c in value.lower()):
                        raise RuntimeError("V31_SUPPLIED_GNOMAD_MISSING_SHA256_PIN:"+pin)
            cmd+=["--gnomad-manifest",str(local_manifest)]
            receipt["external_public_gnomad_manifest"]={
                "source":"EXPLICIT_LOCAL_PUBLIC_SHA256_PINNED_MANIFEST",
                "manifest_sha256":sha256(local_manifest),
                "available_shards":len(shards),
                "full_coverage":coverage.get("v4_1_full_coverage",False),
                "required_case_source_and_index_SHA256":"VERIFIED_IN_NATIVE_STAGE",
                "large_gnomad_downloads_skipped":True}
            write_json(report,receipt)
            print("V31_PUBLIC_GNOMAD_LOCAL_MANIFEST_REGISTERED",
                  "shards="+str(len(shards)),"CASE_SHA256_VERIFICATION_FOLLOWS",
                  flush=True)
        if reused_intake:
            # The program itself produced this model from the original VCF.
            # Its input and output hashes have just been independently verified.
            cmd+=["--case",str(reused_intake)]
            print("E2E_STAGE_1_VERIFIED_INTAKE_RECOVERY" if recovery_provenance
                  else "E2E_STAGE_1_VERIFIED_INTAKE_RESUME",flush=True)
        else:
            cmd+=["--case",str(src)] if a.case else ["--vcf",str(src)]
            if clinical:cmd+=["--clinical",str(clinical)]
        cmd+=["--out",str(out),"--literature-attack-count",str(a.literature_attack_count)]
        if not a.no_online_ai:cmd.append("--enable-online-ai")
        if a.no_resource_download:cmd.append("--no-resource-download")
        if a.skip_full_scientific_modules:
            cmd+=["--no-auto-resources","--no-public-gene-enrichment","--no-network-expansion",
                  "--no-variant-structure","--no-modifier-discovery","--no-dgidb","--no-public-safety-review"]
        print("E2E_STAGE_1_START",flush=True)
        with long_stage_activity("stage1_native_core"):
            result=subprocess.run(cmd,cwd=CODE,check=False)
        receipt["stages"]["stage1_native_core"]={"exit_code":result.returncode,
                                                     "result_file":"FINAL_RESULT_V30_NATIVE.json"}
        write_json(report,receipt)
        if result.returncode:
            raise RuntimeError("STAGE1_NATIVE_PIPELINE_EXIT_"+str(result.returncode))
        if not (out/"FINAL_RESULT_V30_NATIVE.json").is_file():
            raise RuntimeError("STAGE1_RESULT_MISSING")
        print("E2E_STAGE_1_PASS",flush=True)
        print("E2E_STAGE_2_START",flush=True)
        sys.path.insert(0,str(CODE))
        from v30.native_pipeline import resume_from_screen_stage
        with long_stage_activity("stage2_v31_research"):
            final=resume_from_screen_stage(
                out,
                literature_attack_count=a.literature_attack_count,
                enable_online_ai=not a.no_online_ai,
                enable_public_safety_review=not a.skip_full_scientific_modules,
                enable_dgidb=not a.skip_full_scientific_modules,
                enable_cross_source_bridge=not a.skip_full_scientific_modules,
                revalidate_modifiers=not a.skip_full_scientific_modules,
            )
        if not (out/"FINAL_RESULT_V31_NATIVE.json").is_file():
            raise RuntimeError("STAGE2_RESULT_MISSING")
        receipt["stages"]["stage2_v31_research"]={"exit_code":0,"result_file":"FINAL_RESULT_V31_NATIVE.json",
                                                    "approved_entities_screened":final.get("all_approved_entities_screened"),
                                                    "track2_has_hypothesis":bool(final.get("track2_has_hypothesis"))}
        print("E2E_STAGE_2_PASS",flush=True)
        write_json(report,receipt)
        with long_stage_activity("independent_verification"):
            result=subprocess.run([sys.executable,str(ROOT/"verify_challenge_run_v31.py"),
                                   "--run-dir",str(out),
                                   "--require-ai" if not a.no_online_ai else "--no-require-ai"],cwd=ROOT)
        receipt["stages"]["independent_verification"]={"exit_code":result.returncode}
        if result.returncode:raise RuntimeError("E2E_VERIFIER_EXIT_"+str(result.returncode))
        receipt["status"]="PASS"
        receipt["ended_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
        write_json(report,receipt)
        print("E2E_V31_REPRODUCTION_PASS",flush=True)
        return 0
    except BaseException as exc:
        receipt["status"]="FAILED"
        receipt["failure_type"]=type(exc).__name__
        receipt["failure_stage_message"]=str(exc)[:350]
        receipt["ended_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
        write_json(report,receipt)
        print("E2E_V31_REPRODUCTION_FAIL",type(exc).__name__,str(exc)[:200],file=sys.stderr,flush=True)
        traceback.print_exc(file=sys.stderr)
        return 1

if __name__=="__main__":
    raise SystemExit(main())
