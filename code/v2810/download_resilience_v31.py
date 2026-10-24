#!/usr/bin/env python3
"""Resilient *public-data* download policy for the V3.1 jury installer.

One official attempt + three retries, 120 seconds apart; then a verified
alternative, if one exists. Missing databases are recorded and shown in the
existing menu. Never convert a missing source into a fabricated scientific hit.
No patient data and no API credentials are read here.
"""
from __future__ import annotations
from pathlib import Path
import datetime as dt
import json
import os
import sys
import tempfile
import time

SCHEMA = "RDMVA_PUBLIC_DOWNLOAD_WARNINGS_V31_1"
RETRIES = 3
RETRY_DELAY_SECONDS = 120


def _now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def red_warning(message):
    """Show an actual red Windows-console warning, with a plain log fallback."""
    prefix = "[RDMVA DOWNLOAD WARNING] "
    rendered = prefix + str(message)
    if os.name == "nt":
        try:
            import ctypes
            k32 = ctypes.windll.kernel32
            handle = k32.GetStdHandle(-11)
            mode = ctypes.c_uint()
            if k32.GetConsoleMode(handle, ctypes.byref(mode)):
                k32.SetConsoleMode(handle, mode.value | 0x0004)
                print("\x1b[91m" + rendered + "\x1b[0m", flush=True)
                return
        except (OSError, AttributeError, ValueError):
            pass
    if hasattr(sys.stdout, "isatty") and sys.stdout.isatty():
        print("\x1b[91m" + rendered + "\x1b[0m", flush=True)
    else:
        print(rendered, flush=True)


def _error(exc):
    return (type(exc).__name__ + ": " + str(exc).replace("\r", " ").replace("\n", " "))[:320]


class DownloadWarnings:
    """Persist missing public resources without ever storing patient material."""
    def __init__(self, project_root):
        self.path = Path(project_root) / "INSTALL_WARNINGS.json"
        if self.path.is_file():
            try:
                self.data = json.loads(self.path.read_text(encoding="utf-8"))
            except (OSError, ValueError, UnicodeError):
                self.data = {}
        else:
            self.data = {}
        if self.data.get("schema") != SCHEMA:
            self.data = {"schema": SCHEMA, "skipped_sources": [], "recovered_sources": []}

    def _save(self):
        self.data["updated_utc"] = _now()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=self.path.parent, prefix=".rdmva-download-warnings-")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(self.data, f, indent=2, ensure_ascii=False)
                f.write("\n")
            os.replace(tmp, self.path)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)

    def resolved(self, label):
        before = len(self.data["skipped_sources"])
        self.data["skipped_sources"] = [
            entry for entry in self.data["skipped_sources"] if entry["name"] != label
        ]
        if before != len(self.data["skipped_sources"]):
            self._save()

    def recovered(self, label, official_url, alternative_label, error):
        self.resolved(label)
        self.data["recovered_sources"] = [
            x for x in self.data["recovered_sources"] if x["name"] != label
        ]
        self.data["recovered_sources"].append({
            "name": label, "official_url": official_url,
            "alternative": alternative_label,
            "official_failure": error, "recorded_utc": _now(),
        })
        self._save()

    def skipped(self, label, official_url, error, alternative_label="not_available"):
        self.data["skipped_sources"] = [
            x for x in self.data["skipped_sources"] if x["name"] != label
        ]
        self.data["skipped_sources"].append({
            "name": label, "official_url": official_url,
            "status": "SKIPPED_CONTINUED", "reason": error,
            "alternative": alternative_label, "recorded_utc": _now(),
        })
        self._save()
        red_warning("SKIPPED " + label + " after all attempts. "
                    "Installation continues WITHOUT this database.")

    def print_summary(self):
        entries = self.data["skipped_sources"]
        if entries:
            red_warning("DOWNLOAD SUMMARY: " + str(len(entries)) +
                        " database(s) NOT acquired. Scientific coverage is incomplete.")
            for e in entries:
                red_warning("MISSING: " + e["name"] + " | " + e["reason"])
        return len(entries)


def resilient_download(label, official_url, primary, *, fallback=None,
                       fallback_name=None, warnings=None, sleeper=time.sleep,
                       retry_count=RETRIES, retry_delay=RETRY_DELAY_SECONDS):
    """Return verified primary/fallback data or None. Does not mask errors as success.

    Source validation must happen inside primary / fallback. The first attempt
    is immediate, followed by precisely three retries separated by 120 s.
    Test harnesses inject sleeper; production uses real time.sleep.
    """
    if retry_count < 0 or retry_delay < 0:
        raise ValueError("negative retry policy not allowed")
    last_error = "source unavailable"
    for attempt in range(retry_count + 1):
        if attempt:
            red_warning(label + ": retry " + str(attempt) + "/" +
                        str(retry_count) + " in " + str(retry_delay) +
                        " seconds after upstream failure.")
            sleeper(retry_delay)
        try:
            result = primary()
            if result is None:
                raise RuntimeError("source returned no verified download result")
            if warnings:
                warnings.resolved(label)
            if attempt:
                print("RDMVA_DOWNLOAD_RECOVERED_OFFICIAL", label,
                      "attempt=" + str(attempt + 1), flush=True)
            return result
        except Exception as exc:
            last_error = _error(exc)
            red_warning(label + ": official attempt " + str(attempt + 1) +
                        "/" + str(retry_count + 1) + " failed: " + last_error)
    if fallback is not None:
        print("RDMVA_DOWNLOAD_ALTERNATIVE_START", label,
              fallback_name or "verified_alternative", flush=True)
        try:
            result = fallback()
            if result is None:
                raise RuntimeError("alternative returned no verified data")
            if warnings:
                warnings.recovered(label, official_url,
                                   fallback_name or "verified_alternative", last_error)
            print("RDMVA_DOWNLOAD_ALTERNATIVE_PASS", label, flush=True)
            return result
        except Exception as exc:
            last_error += " | alternative failed: " + _error(exc)
            red_warning(label + ": verified alternative failed: " + _error(exc))
    if warnings:
        warnings.skipped(label, official_url, last_error,
                         fallback_name or "not_available")
    else:
        red_warning(label + ": unavailable after all attempts; proceeding")
    return None
