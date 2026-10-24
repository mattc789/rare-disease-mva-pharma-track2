#!/usr/bin/env python3
"""Portable installer/verifier for Rare Disease MVA Pharma V2.8.10.

Supported baseline: Windows x64 + CPython 3.12 + WSL2 Ubuntu/Debian family.
The checkout may live on any drive/folder. No workstation username, drive
letter, external disk identity or pre-existing virtual environment is assumed.

Patient/gated data are never downloaded by this installer.
"""
from __future__ import annotations
from pathlib import Path
import argparse,subprocess,sys,json,os,venv,re,threading,time

ROOT=Path(__file__).resolve().parent
REQ=ROOT/"requirements-v2810.txt"
ENV=ROOT/"env/v28_runtime"

def checked(args,cwd=ROOT,label="step",timeout=None):
    print("V2810_SETUP_START",label,flush=True)
    started=time.monotonic()
    finished=threading.Event()
    def activity():
        while not finished.wait(10):
            # Package installations and independent verifier subprocesses do
            # not expose a trustworthy denominator. Do NOT invent a percent.
            print("V31_SETUP_STAGE_ACTIVITY",label,
                  "elapsed_s="+str(round(time.monotonic()-started)),
                  "percent=unavailable state=working",flush=True)
    observer=threading.Thread(target=activity,daemon=True)
    observer.start()
    try:
        p=subprocess.run(args,cwd=str(cwd),capture_output=True,text=True,
                         errors="replace",timeout=timeout)
    finally:
        finished.set()
        observer.join(timeout=2)
    if p.returncode:
        print("V2810_SETUP_ERROR",label,p.returncode,p.stderr[-1100:],flush=True)
        raise RuntimeError("V2810_SETUP_FAILED_"+label)
    print("V2810_SETUP_OK",label,"elapsed_s="+str(round(time.monotonic()-started)),flush=True)
    return p.stdout[-1200:]

def checked_stream(args,cwd=ROOT,label="step",timeout=None):
    # Large current-dataset acquisitions may run for many minutes. Keep
    # per-file integrity/provenance progress visible in the PowerShell console.
    print("V2810_SETUP_START",label,flush=True)
    p=subprocess.run(args,cwd=str(cwd),timeout=timeout)
    if p.returncode:
        raise RuntimeError("V2810_SETUP_FAILED_"+label+":exit="+str(p.returncode))
    print("V2810_SETUP_OK",label,flush=True)
    return "verified_live_source_data"

def runtime_python():
    return ENV/("Scripts/python.exe" if os.name=="nt" else "bin/python")

def locked_versions():
    out={}
    for line in REQ.read_text(encoding="utf-8").splitlines():
        line=line.strip()
        if not line or line.startswith("#"):continue
        m=re.match(r"^([A-Za-z0-9_.-]+)==([^\s]+)",line)
        if not m:raise RuntimeError("V2810_REQUIREMENT_NOT_EXACTLY_PINNED:"+line)
        out[m.group(1).lower()]=m.group(2)
    return out

def ensure_runtime(verify_only=False):
    py=runtime_python()
    if not py.exists():
        if verify_only:raise RuntimeError("V2810_CORE_PYTHON_ENV_MISSING")
        print("V2810_CREATING_ENV",ENV,flush=True)
        venv.EnvBuilder(with_pip=True,clear=False,symlinks=False).create(ENV)
    if not verify_only:
        checked([str(py),"-m","pip","install","--disable-pip-version-check",
                 "--require-hashes","--only-binary=:all:","-r",str(REQ)],
                label="locked_python_dependencies",timeout=1800)
    probe=(
      "import json,sys,importlib.metadata as im,numpy,PIL,rdkit;"
      "print(json.dumps({'python':sys.version.split()[0],"
      "'numpy':im.version('numpy'),'pillow':im.version('pillow'),'rdkit':im.version('rdkit')}))"
    )
    p=checked([str(py),"-c",probe],label="runtime_probe",timeout=60)
    line=[x for x in p.splitlines() if x.strip().startswith("{")][-1]
    actual=json.loads(line)
    expected=locked_versions()
    for pkg,key in (("numpy","numpy"),("pillow","pillow"),("rdkit","rdkit")):
        if actual[key]!=expected[pkg]:
            raise RuntimeError("V2810_RUNTIME_VERSION_DRIFT:"+pkg+":"+actual[key])
    if not actual["python"].startswith("3.12."):
        raise RuntimeError("V2810_RUNTIME_PYTHON_DRIFT:"+actual["python"])
    return actual

def probe_isolated_wsl():
    """Check the *project-owned* Ubuntu 24.04, never Windows' default distro."""
    code_root=str(ROOT/"code")
    if code_root not in sys.path:sys.path.insert(0,code_root)
    from v2810.resource_bootstrap import project_wsl_distribution, wsl_command
    distro=project_wsl_distribution(ROOT)
    command=wsl_command(
        "bash","-lc",
        'source /etc/os-release && printf "RDMVA:%s:%s:%s:%s" "$ID" "$VERSION_ID" "$(id -un)" "$HOME"'
    )
    # An altered module or mixed checkout must not point to another distro.
    if command[:3]!=["wsl.exe","-d",distro]:
        raise RuntimeError("V31_WSL_PROJECT_DISTRIBUTION_MISMATCH")
    proc=subprocess.run(command,capture_output=True,text=True,
                        errors="replace",timeout=60)
    expected="RDMVA:ubuntu:24.04:rdmva:/home/rdmva"
    if proc.returncode or proc.stdout.strip()!=expected:
        raise RuntimeError("V31_ISOLATED_UBUNTU_24_04_REQUIRED:"+
                           distro+":status="+str(proc.returncode)+
                           ":stdout="+proc.stdout[-120:])
    print("V31_ISOLATED_WSL_READY "+distro+" user=rdmva release=24.04",flush=True)
    return {"distribution":distro,"user":"rdmva","release":"Ubuntu 24.04",
            "other_distributions_untouched":True}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--current-public-data",action="store_true",
      help="Use newly downloaded official public data, recording a new local SHA-256 snapshot.")
    ap.add_argument("--verify-only",action="store_true")
    ap.add_argument("--skip-reference-download",action="store_true",
      help="Do not download missing public reference snapshots; verify only.")
    args=ap.parse_args()
    if sys.version_info[:2]!=(3,12):
        raise RuntimeError("V2810_REQUIRES_CPYTHON_3_12")
    if os.name!="nt":raise RuntimeError("V2810_DEFAULT_WORKSTATION_WINDOWS_WSL2_REQUIRED")
    v28=ROOT/"bootstrap_v28.py";regs=ROOT/"bootstrap_regulatory_v283.py"
    if not all(p.is_file() for p in (REQ,v28,regs)):
        raise RuntimeError("V2810_SOURCE_RELEASE_INCOMPLETE")
    result={
      "schema":"V2810_PUBLIC_SOURCE_SETUP_2",
      "research_installation_cost":"free_source_and_public_data_only",
      "source_code_redistribution_licence":"CC-BY-4.0",
      "private_patient_data_in_setup":False,
      "paid_API_or_licence_required":False,
      "project_root_is_runtime_derived":True,
      "success":False,
      "source_inputs_checked":{},
      "public_data_profile":"CURRENT_OFFICIAL_SOURCES" if args.current_public_data else "FROZEN_CHALLENGE_BASELINE",
    }
    state=ROOT/"INSTALL_V2810_STATUS.json"
    try:
        if args.current_public_data:
            current=ROOT/"bootstrap_current_public_v31.py"
            if not current.is_file():raise RuntimeError("V31_CURRENT_BOOTSTRAP_MISSING")
            current_args=[sys.executable,str(current)]
            if args.verify_only or args.skip_reference_download:
                current_args.append("--verify-only")
            result["current_public_data"]=checked_stream(
                current_args,label="current_official_public_data",
                timeout=9000)
            current_receipt=json.loads((ROOT/"public_inputs/current_v31/CURRENT_PUBLIC_DATA_V31.json").read_text(encoding="utf-8"))
            result["public_data_state"]=current_receipt.get("state")
            result["skipped_source_count"]=len(current_receipt.get("skipped_sources",[]))
            os.environ["RDMVA_CURRENT_PUBLIC_DATA"]="1"
        verify_sources=(args.verify_only or args.skip_reference_download or
                        args.current_public_data)
        v28_args=[sys.executable,str(v28),"--skip-env"]
        regs_args=[sys.executable,str(regs)]
        if verify_sources:
            v28_args.append("--verify-only")
            regs_args.append("--verify-only")
        result["source_inputs_checked"]["public_v28"]=checked(
            v28_args,
            label="public_v28_current_verified" if args.current_public_data else "public_v28_pinned_snapshot",
            timeout=1800)
        result["source_inputs_checked"]["public_v283"]=checked(
            regs_args,
            label="public_v283_current_verified" if args.current_public_data else "public_v283_pinned_snapshot",
            timeout=1200)
        result["runtime"]=ensure_runtime(verify_only=args.verify_only)
        result["wsl"]=probe_isolated_wsl()
        py=runtime_python()
        checked([str(py),"-m","v2810.pipeline","--help"],
                cwd=ROOT/"code",label="cli_smoke",timeout=60)
        result["success"]=True
        print("V2810_SOURCE_INSTALL_SUCCESS",flush=True)
    finally:
        state.write_text(json.dumps(result,indent=2,ensure_ascii=True)+"\n",encoding="utf-8")

if __name__=="__main__":main()
