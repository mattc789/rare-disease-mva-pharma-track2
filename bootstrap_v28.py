#!/usr/bin/env python3
"""Standalone bootstrap for Rare Disease MVA Pharma V2.8.x.

Creates the clean V2.8 runtime and reconstructs the pinned public-data snapshot.
No challenge/patient data are downloaded or copied by this script.

Usage:
  python bootstrap_v28.py
  python bootstrap_v28.py --verify-only
  python bootstrap_v28.py --skip-env
"""
from __future__ import annotations
from pathlib import Path
import argparse, hashlib, json, os, ssl, subprocess, sys, tempfile, urllib.error, urllib.parse, urllib.request, venv, time, datetime as dt

ROOT = Path(__file__).resolve().parent
PUBLIC = ROOT / "public_inputs" / "third_party"
MANIFEST_PATH = (
    ROOT/"public_inputs/current_v31/public_data_manifest_v28_current.json"
    if os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1"
    else ROOT/"config/public_data_manifest_v28.json"
)
REQ = ROOT / "requirements-v28.txt"
ENV = ROOT / "env" / "v28_runtime"
UA = "RareDiseaseMVAPharma-V28-Bootstrap/1.0"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def load_manifest():
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

def check_file(path: Path, spec: dict) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing"
    size = path.stat().st_size
    if size != int(spec["bytes"]):
        return False, f"size {size} != {spec['bytes']}"
    got = sha256(path)
    if got != spec["sha256"]:
        return False, f"sha256 {got} != {spec['sha256']}"
    return True, "ok"

def download_exact(spec: dict):
    dst = PUBLIC / spec["path"]
    ok, why = check_file(dst, spec)
    if ok:
        print("OK", spec["path"], "(cached)")
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    # Stable stage path enables verified segment resume after an interruption.
    tmp = dst.with_name(dst.name+".rdmva-download")
    try:
        # Pinned historical scientific archives: accelerated if upstream
        # supports validated byte ranges, without relaxing exact SHA-256 pins.
        code=str(ROOT/"code")
        if code not in sys.path:sys.path.insert(0,code)
        from v2810.scientific_download_v31 import download_public_data
        download_public_data(spec["url"], tmp,
                             expected_bytes=int(spec["bytes"]),
                             max_bytes=int(spec["bytes"])+16*1024*1024,
                             expected_sha256=spec["sha256"],
                             label=spec["path"])
        ok, why = check_file(tmp, spec)
        if not ok:
            raise RuntimeError(
                f"Pinned public snapshot drift for {spec['path']}: {why}. "
                "Refusing to silently change the scientific input snapshot."
            )
        os.replace(tmp, dst)
        print("FETCHED", spec["path"])
    finally:
        if tmp.exists():
            tmp.unlink()

# Temporary network errors are not scientific evidence of absence.
# Checkpoint every successful public ChEMBL API response inside one exact
# acquisition attempt. No patient data, keys or API response bodies go to Git.
RETRYABLE_CHEMBL_HTTP={408,425,429,500,502,503,504}
CHEMBL_MAX_ATTEMPTS=5
CHEMBL_CACHE_SCHEMA="RDMVA_V31_CHEMBL_VERIFIED_PAGE_1"
CHEMBL_BASE="https://www.ebi.ac.uk/chembl/api/data/"

class ChEMBLUpstreamBadResponse(RuntimeError):
    """HTTP success but not a validated ChEMBL API response."""


def _chembl_verified_url(url: str):
    parsed=urllib.parse.urlsplit(url)
    if (parsed.scheme!="https" or parsed.hostname!="www.ebi.ac.uk" or
        not parsed.path.startswith("/chembl/api/data/") or
        parsed.username or parsed.password):
        raise RuntimeError("V31_CHEMBL_UNTRUSTED_API_URL")
    return url

def _validate_chembl_api_response(url,data):
    # Prevent a transient HTML/JSON "error" envelope with HTTP 200 from
    # turning into an empty drug/target universe. All current calls use
    # ChEMBL list resources with mandatory page_meta.
    resource=urllib.parse.urlsplit(url).path.rsplit("/",1)[-1].split(".",1)[0]
    field={"molecule":"molecules","mechanism":"mechanisms",
           "target":"targets"}.get(resource)
    if field is None:
        raise RuntimeError("V31_CHEMBL_UNSUPPORTED_RESOURCE:"+resource)
    if (not isinstance(data,dict) or
        not isinstance(data.get("page_meta"),dict) or
        not isinstance(data.get(field),list)):
        raise ChEMBLUpstreamBadResponse(
            "V31_CHEMBL_HTTP_200_BUT_NOT_A_VALID_"+field.upper()+"_PAGE")
    return data

def _chembl_cache_path(url: str):
    scope=os.environ.get("RDMVA_CHEMBL_CACHE_SCOPE") or ""
    if not scope:
        return None
    if len(scope)!=20 or any(c not in "0123456789abcdef" for c in scope):
        raise RuntimeError("V31_CHEMBL_UNSAFE_CACHE_SCOPE")
    _chembl_verified_url(url)
    digest=hashlib.sha256(url.encode("utf-8")).hexdigest()
    return PUBLIC/"chembl"/".v31-page-cache"/scope/(digest+".json")

def _chembl_payload_sha(obj):
    raw=json.dumps(obj,ensure_ascii=False,sort_keys=True,separators=(",",":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()

def _chembl_read_verified_cache(path,url,cap):
    if path is None or not path.is_file():
        return None
    try:
        if path.stat().st_size>cap*2+1_000_000:
            raise ValueError("cached response unusually large")
        record=json.loads(path.read_text(encoding="utf-8"))
        payload=record.get("response")
        if (record.get("schema")!=CHEMBL_CACHE_SCHEMA or
            record.get("url")!=url or not isinstance(payload,dict) or
            record.get("response_sha256")!=_chembl_payload_sha(payload)):
            raise ValueError("cached response provenance mismatch")
        _validate_chembl_api_response(url,payload)
        print("V31_CHEMBL_CACHED_VERIFIED",
              hashlib.sha256(url.encode("utf-8")).hexdigest()[:12],flush=True)
        return payload
    except (OSError,ValueError,KeyError,TypeError) as exc:
        # Never use invalid pages; safely fetch again. The corrupt cache is
        # retained for forensic diagnosis until an exact verified replacement
        # succeeds. No broad directory deletion is performed.
        print("V31_CHEMBL_PAGE_CACHE_INVALID_REFETCH",
              type(exc).__name__,flush=True)
        return None

def _chembl_checkpoint(path,url,obj):
    if path is None:return
    path.parent.mkdir(parents=True,exist_ok=True)
    payload={
        "schema":CHEMBL_CACHE_SCHEMA,
        "url":url,
        "captured_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "response_sha256":_chembl_payload_sha(obj),
        "response":obj}
    fd,tmpname=tempfile.mkstemp(prefix=".v31-chembl-api-",dir=path.parent)
    os.close(fd)
    temp=Path(tmpname)
    try:
        temp.write_text(json.dumps(payload,ensure_ascii=False)+"\n",encoding="utf-8")
        os.replace(temp,path)
    finally:
        if temp.exists():temp.unlink()

def _retry_after_seconds(exc,backoff):
    header=exc.headers.get("Retry-After") if exc.headers is not None else None
    try:
        seconds=int(header)
        return max(1,min(90,seconds))
    except (ValueError,TypeError):
        return backoff

def _fetch_json_once(url: str, cap: int):
    req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
    try:
        with urllib.request.urlopen(req,timeout=90) as r:
            raw=r.read(cap+1)
    except urllib.error.HTTPError:
        raise
    except urllib.error.URLError as exc:
        # Local Windows root certificates may be absent from Python OpenSSL.
        # For ONLY verified certificate-chain errors use Schannel curl, which
        # performs proper TLS validation. Never use -k or verify=False.
        if os.name!="nt" or "CERTIFICATE_VERIFY_FAILED" not in str(exc):
            raise
        with tempfile.TemporaryDirectory(prefix="rdmva-verified-chembl-") as folder:
            path=Path(folder)/"chembl.json"
            cmd=["curl.exe","--fail","--location","--retry","2",
                 "--retry-all-errors","--connect-timeout","25",
                 "--max-time","120","--max-filesize",str(cap),
                 "--user-agent",UA,"--silent","--show-error",
                 "--output",str(path),"--write-out","%{http_code}",url]
            done=subprocess.run(cmd,capture_output=True,text=True,timeout=140)
            if done.returncode!=0 or done.stdout.strip()!="200":
                raise RuntimeError("V31_VERIFIED_CURL_FALLBACK_FAILED:"+str(done.returncode)) from exc
            raw=path.read_bytes()
        print("V31_CHEMBL_WINDOWS_TLS_FALLBACK_VERIFIED",flush=True)
    if len(raw)>cap:
        raise RuntimeError("ChEMBL response exceeded cap")
    data=json.loads(raw)
    return _validate_chembl_api_response(url,data)

def fetch_json(url: str, cap: int = 50_000_000) -> dict:
    _chembl_verified_url(url)
    # Local licensed database first whenever the official SQLite profile is
    # installed. The REST backend is retained only as an explicit alternative.
    backend=os.environ.get("RDMVA_CHEMBL_BACKEND","sqlite").strip().lower()
    if backend=="sqlite" and (PUBLIC/"chembl37/CHEMBL_SQLITE_PROVENANCE_V31.json").is_file():
        import sys
        code_root=str(ROOT/"code")
        if code_root not in sys.path:sys.path.insert(0,code_root)
        from v2810.chembl_sqlite_v31 import local_json
        return local_json(url)
    if backend=="sqlite" and os.environ.get("RDMVA_CURRENT_PUBLIC_DATA")=="1":
        raise RuntimeError("V31_CHEMBL_SQLITE_REQUIRED_NOT_READY")
    cached=_chembl_cache_path(url)
    prior=_chembl_read_verified_cache(cached,url,cap)
    if prior is not None:return prior
    for attempt in range(1,CHEMBL_MAX_ATTEMPTS+1):
        try:
            result=_fetch_json_once(url,cap)
            _chembl_checkpoint(cached,url,result)
            return result
        except urllib.error.HTTPError as exc:
            # A 4xx such as 401/403/404 is NOT a transient 500.
            # Do not retry permission/configuration/endpoint errors.
            if exc.code not in RETRYABLE_CHEMBL_HTTP:
                print("V31_CHEMBL_NONRETRYABLE_HTTP",exc.code,flush=True)
                raise
            if attempt==CHEMBL_MAX_ATTEMPTS:
                raise RuntimeError(
                    "V31_CHEMBL_UPSTREAM_TEMPORARY_HTTP_EXHAUSTED:"
                    f"status={exc.code}:attempts={attempt}:rerun_installer_to_resume_verified_pages"
                ) from exc
            delay=_retry_after_seconds(exc,min(30,2**attempt))
            print("V31_CHEMBL_RETRY_HTTP",exc.code,
                  "attempt",attempt,"of",CHEMBL_MAX_ATTEMPTS,
                  "wait_seconds",delay,flush=True)
        except (urllib.error.URLError,TimeoutError,ConnectionError) as exc:
            if attempt==CHEMBL_MAX_ATTEMPTS:
                raise RuntimeError(
                    "V31_CHEMBL_UPSTREAM_TIMEOUT_EXHAUSTED:"
                    f"attempts={attempt}:rerun_installer_to_resume_verified_pages"
                ) from exc
            delay=min(30,2**attempt)
            print("V31_CHEMBL_RETRY_NETWORK",type(exc).__name__,
                  "attempt",attempt,"of",CHEMBL_MAX_ATTEMPTS,
                  "wait_seconds",delay,flush=True)
        except (json.JSONDecodeError,ChEMBLUpstreamBadResponse) as exc:
            if attempt==CHEMBL_MAX_ATTEMPTS:
                raise RuntimeError(
                    "V31_CHEMBL_BAD_RESPONSE_EXHAUSTED:"
                    f"attempts={attempt}:refusing_empty_scientific_evidence"
                ) from exc
            delay=min(30,2**attempt)
            print("V31_CHEMBL_RETRY_BAD_RESPONSE",type(exc).__name__,
                  "attempt",attempt,"of",CHEMBL_MAX_ATTEMPTS,
                  "wait_seconds",delay,flush=True)
        time.sleep(delay)
    raise RuntimeError("V31_CHEMBL_RETRY_INTERNAL_FAILURE")

def write_exact(path: Path, text: str, spec: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".v28-generated-", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        tmp.write_text(text, encoding="utf-8", newline="\n")
        if spec.get("capture_current") is True:
            spec["bytes"]=tmp.stat().st_size
            spec["sha256"]=sha256(tmp)
        else:
            ok, why = check_file(tmp, spec)
            if not ok:
                raise RuntimeError(
                    f"Pinned ChEMBL snapshot drift for {spec['path']}: {why}. "
                    "The upstream API no longer reproduces the 2026-10-09 snapshot."
                )
        os.replace(tmp, path)
        print("GENERATED", spec["path"])
    finally:
        if tmp.exists():
            tmp.unlink()

def paged(endpoint: str, params: dict, field: str) -> list:
    rows = []
    params = dict(params)
    params.setdefault("limit", 1000)
    params.setdefault("offset", 0)
    url = "https://www.ebi.ac.uk/chembl/api/data/" + endpoint + ".json?" + urllib.parse.urlencode(params, safe=",")
    seen=set()
    page=0
    while url:
        _chembl_verified_url(url)
        if url in seen:
            raise RuntimeError("V31_CHEMBL_PAGINATION_CYCLE")
        seen.add(url)
        j=fetch_json(url)
        data=j.get(field)
        if not isinstance(data,list):
            raise RuntimeError("V31_CHEMBL_PAGINATION_MISSING_ROWS:"+field)
        metadata=j.get("page_meta")
        if not isinstance(metadata,dict):
            raise RuntimeError("V31_CHEMBL_PAGINATION_METADATA_MISSING:"+field)
        rows.extend(data)
        page+=1
        print("V31_CHEMBL_PAGE_VERIFIED",endpoint,page,
              "rows_total",len(rows),flush=True)
        nxt=metadata.get("next")
        if nxt is not None and not isinstance(nxt,str):
            raise RuntimeError("V31_CHEMBL_NEXT_PAGE_INVALID")
        url=urllib.parse.urljoin("https://www.ebi.ac.uk",nxt) if nxt else None
        if len(rows)>200_000 or page>500:
            raise RuntimeError("unexpected ChEMBL pagination expansion")
    return rows

def build_chembl(manifest: dict):
    specs = {x["path"]: x for x in manifest["generated_files"]}
    cdir = PUBLIC / "chembl"
    catalog_path = cdir / "catalog_phase4.json"
    mech_path = cdir / "mechanisms_records.json"
    target_path = cdir / "target_metadata_records.json"

    checks = [check_file(cdir / p.split("/",1)[1], s)[0] for p,s in specs.items()]
    if all(checks):
        for p in specs:
            print("OK", p, "(cached)")
        return

    fields = (
        "molecule_chembl_id,pref_name,molecule_hierarchy,molecule_type,"
        "molecule_synonyms,max_phase,withdrawn_flag,availability_type,veterinary,"
        "first_approval,oral,parenteral,topical,prodrug,atc_classifications"
    )
    catalog = paged("molecule", {"max_phase": 4, "only": fields}, "molecules")
    write_exact(
        catalog_path,
        json.dumps(catalog, ensure_ascii=False),
        specs["chembl/catalog_phase4.json"],
    )

    mechanisms = paged("mechanism", {}, "mechanisms")
    write_exact(
        mech_path,
        json.dumps(mechanisms, indent=2, ensure_ascii=False) + "\n",
        specs["chembl/mechanisms_records.json"],
    )

    catalogue_ids = {x["molecule_chembl_id"] for x in catalog}
    parent_ids = {
        ((x.get("molecule_hierarchy") or {}).get("parent_chembl_id") or x["molecule_chembl_id"])
        for x in catalog
    }
    target_ids = set()
    for m in mechanisms:
        did = m.get("molecule_chembl_id")
        pid = m.get("parent_molecule_chembl_id") or did
        if pid not in parent_ids and did not in catalogue_ids:
            continue
        if m.get("target_chembl_id"):
            target_ids.add(m["target_chembl_id"])

    records = {}
    ids = sorted(target_ids)
    for i in range(0, len(ids), 40):
        batch = ids[i:i+40]
        print("V31_CHEMBL_TARGET_METADATA_PROGRESS",
              min(i+40,len(ids)),"of",len(ids),flush=True)
        j = fetch_json(
            "https://www.ebi.ac.uk/chembl/api/data/target.json?"
            + urllib.parse.urlencode({"limit": 1000, "target_chembl_id__in": ",".join(batch)}, safe=",")
        )
        if (j.get("page_meta") or {}).get("next"):
            raise RuntimeError("unexpected ChEMBL target pagination")
        for row in j.get("targets", []):
            if row.get("target_chembl_id") in batch:
                records[row["target_chembl_id"]] = row
        time.sleep(0.05)
    ordered = [records[k] for k in sorted(records)]
    write_exact(
        target_path,
        json.dumps(ordered, indent=2),
        specs["chembl/target_metadata_records.json"],
    )

def env_python() -> Path:
    return ENV / ("Scripts/python.exe" if os.name == "nt" else "bin/python")

def ensure_env():
    py = env_python()
    if not py.exists():
        print("CREATING", ENV)
        venv.EnvBuilder(with_pip=True, clear=False, symlinks=False).create(ENV)
    subprocess.check_call([str(py), "-m", "pip", "install", "--disable-pip-version-check", "-r", str(REQ)])
    probe = (
        "import sys,numpy,PIL,rdkit;"
        "print(sys.version.split()[0],numpy.__version__,PIL.__version__,rdkit.__version__)"
    )
    subprocess.check_call([str(py), "-c", probe])
    print("ENV_OK", py)

def verify_all(manifest: dict):
    bad = []
    for spec in manifest["static_files"] + manifest["generated_files"]:
        p = PUBLIC / spec["path"]
        ok, why = check_file(p, spec)
        print("VERIFY", spec["path"], "OK" if ok else why)
        if not ok:
            bad.append((spec["path"], why))
    if bad:
        raise RuntimeError("Public snapshot verification failed: " + repr(bad))
    print("PUBLIC_SNAPSHOT_OK", manifest["version"])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--skip-env", action="store_true")
    args = ap.parse_args()
    manifest = load_manifest()

    if not args.verify_only:
        if not args.skip_env:
            ensure_env()
        for spec in manifest["static_files"]:
            download_exact(spec)
        build_chembl(manifest)

    verify_all(manifest)
    if not args.skip_env:
        if not env_python().exists():
            raise RuntimeError("V2.8 runtime environment missing; run without --verify-only")
        print("READY", env_python())

if __name__ == "__main__":
    main()
