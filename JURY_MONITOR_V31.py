#!/usr/bin/env python3
"""Unprivileged, read-only console monitor for an RDMVA V3.1 scientific run.

Uses Python's standard library and reads ONLY known, bounded status summaries.
Does not enumerate processes, execute commands, scan folders, access patient
variant/clinical files, open network connections, or alter the scientific run.
The main runner writes its own progress directly to its console independently.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import time
from pathlib import Path

MAX_STATUS_BYTES = 2_000_000
STATUS_FILES = {
    "execution": "E2E_EXECUTION_V31.json",
    "native": "STATUS_V30_NATIVE.json",
    "gnomad": "resource_bootstrap/gnomad/GNOMAD_AUTO_RESOURCE_AUDIT.json",
    "vep": "resource_bootstrap/vep/VEP_BOOTSTRAP_STATUS.json",
    "genomics": "GENOMIC_REQUIRED_EVIDENCE_V31.json",
    "screen": "screen/screen_summary.json",
    "revalidation": "E2E_VERIFIER_REVALIDATION_V31.json",
}


def read_status(root: Path, name: str) -> dict | None:
    """Only a predefined, small status JSON, never patient payloads."""
    if name not in STATUS_FILES:
        raise ValueError("Unrecognized status file")
    path = root / STATUS_FILES[name]
    try:
        # Do not follow symlinks that lead outside the selected run directory.
        if path.resolve() != root and root not in path.resolve().parents:
            return None
        if not path.is_file() or path.stat().st_size > MAX_STATUS_BYTES:
            return None
        with path.open("r", encoding="utf-8-sig") as stream:
            value = json.load(stream)
        return value if isinstance(value, dict) else None
    except (OSError, ValueError, UnicodeError, TypeError):
        return None


def show_progress(root: Path) -> bool:
    utc = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n=== RDMVA V3.1 SCIENTIFIC RUN MONITOR ({utc} UTC) ===")
    execution = read_status(root, "execution")
    native = read_status(root, "native")
    if execution:
        print("RUN_STATUS=" + str(execution.get("status", "UNKNOWN")))
        print("SCIENCE_FULL_PROFILE=" + str(execution.get("full_scientific_modules")))
        print("ONLINE_AI_REQUESTED=" + str(execution.get("online_ai_requested")))
        if execution.get("failure_stage_message"):
            print("FAILURE=" + str(execution["failure_stage_message"])[:240])
    else:
        print("RUN_STATUS=NOT_STARTED_OR_RECEIPT_UNREADABLE")

    if native:
        stage = str(native.get("current_stage", "NOT_YET_REPORTED"))
        print("NATIVE_STAGE=" + stage)
        entry = native.get("stages") or {}
        stage_info = entry.get(stage, {}) if isinstance(entry, dict) else {}
        if isinstance(stage_info, dict):
            if "state" in stage_info:
                print("NATIVE_STAGE_STATE=" + str(stage_info["state"]))
            if "utc" in stage_info:
                print("NATIVE_STAGE_STARTED_UTC=" + str(stage_info["utc"]))
    else:
        print("NATIVE_STAGE=NOT_YET_REPORTED")

    gnomad = read_status(root, "gnomad")
    if gnomad:
        print("GNOMAD_STATE=" + str(gnomad.get("completion", "UNKNOWN")))
        if gnomad.get("installed_jobs") is not None:
            print("GNOMAD_SHARDS_VERIFIED=" + str(gnomad["installed_jobs"])
                  + "/" + str(gnomad.get("required_jobs", "?")))
        active = gnomad.get("active_job")
        if isinstance(active, dict):
            print("GNOMAD_ACTIVE_JOB=" + str(active.get("cohort", "?"))
                  + "/chr" + str(active.get("chrom", "?")) + "/"
                  + str(active.get("phase", "?")))
        if gnomad.get("error_summary"):
            print("GNOMAD_ERROR=" + str(gnomad["error_summary"])[:240])

    vep = read_status(root, "vep")
    if vep:
        print("VEP_BOOTSTRAP=" + str(vep.get("state", "UNKNOWN")))
    genomics = read_status(root, "genomics")
    if genomics:
        print("GENOMICS_EVIDENCE_GATE=" + str(genomics.get("required_science_pass")))
    screen = read_status(root, "screen")
    if screen:
        print("APPROVED_ENTITIES_SCREENED=" + str(screen.get("screened", "UNKNOWN")))
    revalidated = read_status(root, "revalidation")
    if revalidated:
        print("POST_RUN_REVALIDATION=" + str(revalidated.get("technical_verification", "UNKNOWN")))

    # A receipt written before an unexpected interruption may be stale.
    # Do not infer process liveness from status files or read process lists.
    if execution and execution.get("status") == "RUNNING":
        print("NOTE: RUNNING is a saved receipt, not proof the process is still live.")
    completed = bool(execution and execution.get("status") in ("PASS", "FAILED", "INTERRUPTED"))
    print("TERMINAL_EXECUTION_STATUS=" + str(completed))
    return completed


def main() -> int:
    parser = argparse.ArgumentParser(description="Read-only V3.1 run monitor")
    parser.add_argument("--out", required=True)
    parser.add_argument("--interval", type=int, default=60)
    parser.add_argument("--wait-for-start", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.interval <= 3600:
        parser.error("--interval must be between 0 and 3600 seconds")
    root = Path(args.out).expanduser().resolve()
    if not root.is_dir():
        if not args.wait_for_start:
            parser.error("Output directory not found")
        print("Waiting up to 120 seconds for the research runner to create its directory.", flush=True)
        for _ in range(120):
            if root.is_dir():
                break
            time.sleep(1)
        if not root.is_dir():
            parser.error("Research run did not start within 120 seconds")
    try:
        while True:
            finished = show_progress(root)
            if finished or args.interval == 0:
                return 0
            print(f"Next status update in {args.interval} seconds. "
                  "Ctrl+C stops ONLY this monitor.", flush=True)
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\nOnly the monitoring window was stopped.", flush=True)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
