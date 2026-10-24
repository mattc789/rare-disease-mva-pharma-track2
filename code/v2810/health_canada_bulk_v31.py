#!/usr/bin/env python3
"""Health Canada official DPD bulk-export fallback for public-source snapshots.

Activated ONLY when the Health Canada DPD REST API responds with HTTP 503.
Downloads fresh official ZIP extracts, converts the two relevant product-status
categories to the DPD API-compatible JSON subset used by this application,
and records origin URLs and cryptographic hashes. No patient data or keys.
The official extracts are published separately from the DPD REST service:
https://www.canada.ca/en/health-canada/services/drugs-health-products/drug-products/drug-product-database/what-data-extract-drug-product-database.html
"""
from __future__ import annotations
import csv
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile

BASE = "https://www.canada.ca/content/dam/hc-sc/documents/services/drug-product-database/"
BULK_FILES = ("drug_ap.zip", "drug.zip", "ingred_ap.zip", "ingred.zip")
STATUS_FILES = {
    "ca_approved_20261010.json": ("drug_ap.zip",),
    "ca_marketed_20261010.json": ("drug.zip",),
    "ca_ingredients_20261010.json": ("ingred_ap.zip", "ingred.zip"),
}
FIELDS_DRUG = 14
FIELDS_INGREDIENT = 15


def _fetch_small_official_zip(url: str, path: Path) -> dict:
    # Health Canada's static ZIP host is reachable with Windows Schannel even
    # when the Python urllib opener encounters upstream TCP resets.
    # Use strict TLS and byte/HTTP limits; never use --insecure.
    if not url.startswith(BASE) or not url.endswith(".zip"):
        raise RuntimeError("HEALTH_CANADA_DPD_UNAPPROVED_BULK_URL")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(prefix=".official-dpd-transfer-", suffix=".zip", dir=path.parent)
    os.close(fd)
    try:
        proc = subprocess.run(
            ["curl.exe", "--fail", "--location", "--retry", "2",
             "--retry-all-errors", "--silent", "--show-error",
             "--connect-timeout", "20", "--max-time", "120",
             "--max-filesize", "15000000",
             "--output", temp_path, "--write-out", "%{http_code}", url],
            capture_output=True, text=True, timeout=165,
        )
        if proc.returncode or proc.stdout.strip() != "200":
            raise RuntimeError("HEALTH_CANADA_DPD_OFFICIAL_BULK_HTTPS_FAILED:"
                               + str(proc.returncode) + ":" + proc.stdout[-3:]
                               + ":" + proc.stderr[-140:])
        size = os.path.getsize(temp_path)
        if size < 1000 or size > 15_000_000:
            raise RuntimeError("HEALTH_CANADA_DPD_ARCHIVE_SIZE_INVALID")
        h = hashlib.sha256()
        with open(temp_path, "rb") as fp:
            for block in iter(lambda: fp.read(1024 * 1024), b""):
                h.update(block)
        os.replace(temp_path, path)
        return {"bytes": size, "sha256": h.hexdigest(),
                "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def _parse_zip(path: Path, csv_name: str, width: int) -> list[list[str]]:
    with zipfile.ZipFile(path, "r") as z:
        entries = z.namelist()
        if entries != [csv_name]:
            raise RuntimeError("HEALTH_CANADA_DPD_UNEXPECTED_ZIP_ENTRIES:" + path.name)
        info = z.getinfo(csv_name)
        if info.file_size > 80_000_000 or info.file_size < 1000:
            raise RuntimeError("HEALTH_CANADA_DPD_BULK_FILE_SIZE_INVALID:" + path.name)
        with z.open(info) as raw, io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as stream:
            records = list(csv.reader(stream))
    if not records or any(len(row) != width for row in records):
        raise RuntimeError("HEALTH_CANADA_DPD_BULK_SCHEMA_MISMATCH:" + path.name)
    return records


def _write_json_atomic(path: Path, records: list[dict]) -> tuple[int, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_path = tempfile.mkstemp(prefix=".official-dpd-", suffix=".json", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fp:
            json.dump(records, fp, ensure_ascii=False, separators=(",", ":"))
            fp.write("\n")
        h = hashlib.sha256()
        size = 0
        with open(temp_path, "rb") as fp:
            for block in iter(lambda: fp.read(1024 * 1024), b""):
                size += len(block)
                h.update(block)
        os.replace(temp_path, path)
        return size, h.hexdigest()
    finally:
        if os.path.exists(temp_path):
            os.unlink(temp_path)


def recover_official_dpd_bulk(raw_root: Path, download_fn, api_urls: dict) -> dict:
    """Fetch all four official extracts, create three API-compatible JSON snapshots.

    download_fn is bootstrap_current_public_v31.curl_download, maintaining HTTPS
    certificate verification and receipt logging. Do not enable on parse errors.
    """
    archives: dict[str, dict] = {}
    bulk_paths: dict[str, Path] = {}
    for archive in BULK_FILES:
        url = BASE + archive
        path = raw_root / "canada_official_dpd_bulk" / archive
        provenance = _fetch_small_official_zip(url, path)
        archives[archive] = {
            "url": url,
            "bytes": provenance["bytes"],
            "sha256": provenance["sha256"],
            "retrieved_utc": provenance["retrieved_utc"],
        }
        bulk_paths[archive] = path

    drug_ap_rows = _parse_zip(bulk_paths["drug_ap.zip"], "drug_ap.txt", FIELDS_DRUG)
    drug_rows = _parse_zip(bulk_paths["drug.zip"], "drug.txt", FIELDS_DRUG)
    ingredient_ap_rows = _parse_zip(bulk_paths["ingred_ap.zip"], "ingred_ap.txt", FIELDS_INGREDIENT)
    ingredient_rows = _parse_zip(bulk_paths["ingred.zip"], "ingred.txt", FIELDS_INGREDIENT)
    codes = set()
    products: dict[str, list[dict]] = {}
    for name, rows in (
        ("ca_approved_20261010.json", drug_ap_rows),
        ("ca_marketed_20261010.json", drug_rows),
    ):
        out = []
        for row in rows:
            code = int(row[0])
            ais = int(row[8]) if row[8] else 0
            if code <= 0 or ais < 0:
                raise RuntimeError("HEALTH_CANADA_DPD_INVALID_DRUG_VALUE")
            codes.add(code)
            out.append({
                "drug_code": code,
                "class_name": row[2],
                "drug_identification_number": row[3],
                "brand_name": row[4],
                "number_of_ais": str(ais),
            })
        products[name] = out

    # Retain every official ingredient row. Deduplicating identical
    # ingredient/strength fields would make number_of_ais disagree with
    # official product records and hide otherwise-complete combinations.
    ingredients = []
    for row in ingredient_ap_rows + ingredient_rows:
        code = int(row[0])
        if code not in codes:
            continue
        ingredients.append({
            "drug_code": code,
            "ingredient_name": row[2],
        })
    products["ca_ingredients_20261010.json"] = ingredients

    response = {}
    for name, rows in products.items():
        if len(rows) < 200:
            raise RuntimeError("HEALTH_CANADA_DPD_TOO_FEW_RECORDS:" + name)
        path = raw_root / "regulatory_v283" / name
        size, digest = _write_json_atomic(path, rows)
        origins = STATUS_FILES[name]
        response[name] = {
            "bytes": size,
            "sha256": digest,
            "source_url": archives[origins[0]]["url"],
            "requested_source_url": api_urls[name],
            "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "download_transport": "health_canada_official_bulk_export_https_converted",
            "fallback_reason": "HEALTH_CANADA_DPD_API_HTTP_503",
            "fallback_source": "Health Canada Drug Product Database (DPD) official UTF-8 data extract",
            "fallback_conversion": "approved_marketed_api_compatible_field_subset_v1",
            "fallback_ingredient_coverage": "approved_and_marketed_only",
            "bulk_archive_sha256_and_origin": {fn: archives[fn] for fn in origins},
            "record_count": len(rows),
            "not_a_direct_api_response": True,
        }
        print("V31_HEALTH_CANADA_OFFICIAL_BULK_FALLBACK",
              name, len(rows), "sha256=" + digest, flush=True)
    return response
