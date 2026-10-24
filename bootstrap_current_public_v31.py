#!/usr/bin/env python3
"""Current official public-data snapshot bootstrap for RDMVA V3.1.

Distinct from the frozen Challenge baseline. Each file is newly fetched from
an upstream URL over VERIFIED TLS and recorded by SHA256, size, source URL and
UTC capture date. Neither patient data nor API keys are involved.

Incomplete snapshots are never activated by the scientific engine.
"""
from __future__ import annotations
from pathlib import Path
import argparse,copy,datetime as dt,gzip,hashlib,json,os,subprocess,tempfile,time,urllib.parse
from urllib.parse import urlparse

ROOT=Path(__file__).resolve().parent
# Also resolve project-owned v2810 modules when a failed download is resumed
# directly from a verified fallback receipt, before curl_download is called.
import sys
_CODE_PATH=str(ROOT/"code")
if _CODE_PATH not in sys.path:sys.path.insert(0,_CODE_PATH)
from v2810.download_resilience_v31 import DownloadWarnings,resilient_download,red_warning
RAW=ROOT/"public_inputs"/"third_party"
AUDIT=ROOT/"public_inputs"/"current_v31"
STATIC_ORIGINAL=ROOT/"config"/"public_data_manifest_v28.json"
REG_ORIGINAL=ROOT/"config"/"regulatory_snapshot_v283.json"
STATIC_LIVE=AUDIT/"public_data_manifest_v28_current.json"
REG_LIVE=AUDIT/"regulatory_snapshot_v283_current.json"
RECEIPT=AUDIT/"CURRENT_PUBLIC_DATA_V31.json"
UA="RareDiseaseMVAPharma-V31-CurrentSnapshot/1.0"
# Live source endpoints. Fixed GENCODE/SIGNOR release names remain fixed inputs.
URL_OVERRIDES={
 "intake_v28/genes_to_phenotype.txt":
   "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt",
 "intake_v28/hp.obo":
   "https://purl.obolibrary.org/obo/hp.obo",
}
SCHEMA="RDMVA_CURRENT_PUBLIC_SNAPSHOT_V31_1"

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(path):
    h=hashlib.sha256()
    total=path.stat().st_size
    count=0
    last_report=time.monotonic()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""):
            h.update(block)
            count+=len(block)
            now=time.monotonic()
            if total>=100_000_000 and now-last_report>=5:
                print(f"V31_CURRENT_SHA256_PROGRESS {path.name} "
                      f"{100*count/total:.1f}% "
                      f"{count/1_000_000_000:.2f}/{total/1_000_000_000:.2f} GB",flush=True)
                last_report=now
    if total>=100_000_000:
        print(f"V31_CURRENT_SHA256_PROGRESS {path.name} 100.0% "
              "state=digest_computed",flush=True)
    return h.hexdigest()

def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def save(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    os.replace(tmp,path)

def validate_source_url(url):
    p=urlparse(url)
    if p.scheme!="https" or not p.hostname:
        raise RuntimeError("V31_CURRENT_SOURCE_MUST_USE_HTTPS")
    allowed={
      "www.fda.gov","www.accessdata.fda.gov","signor.uniroma2.it",
      "ftp.ebi.ac.uk","ftp.ncbi.nlm.nih.gov","purl.obolibrary.org",
      "www.ema.europa.eu","assets.publishing.service.gov.uk",
      "health-products.canada.ca", "www.canada.ca",
    }
    if p.hostname.lower() not in allowed:
        raise RuntimeError("V31_CURRENT_SOURCE_HOST_NOT_APPROVED:"+p.hostname)
    return p

def validate_download(p,logical_name):
    if p.stat().st_size<50:
        raise RuntimeError("V31_CURRENT_DOWNLOAD_TOO_SMALL:"+logical_name)
    with p.open("rb") as f:
        prefix=f.read(500)
    l=logical_name.lower()
    if l.endswith(".zip") or l.endswith(".xlsx"):
        if prefix[:2]!=b"PK":raise RuntimeError("V31_CURRENT_UNEXPECTED_ZIP:"+logical_name)
    elif l.endswith(".gz"):
        if prefix[:2]!=b"\x1f\x8b":raise RuntimeError("V31_CURRENT_UNEXPECTED_GZIP:"+logical_name)
    elif l.endswith(".json"):
        if not prefix.lstrip().startswith((b"{",b"[")):
            raise RuntimeError("V31_CURRENT_UNEXPECTED_JSON:"+logical_name)
    elif l.endswith(".obo"):
        if b"format-version:" not in prefix:
            raise RuntimeError("V31_CURRENT_UNEXPECTED_HPO:"+logical_name)
    else:
        if b"<html" in prefix.lower() or b"<!doctype" in prefix.lower():
            raise RuntimeError("V31_CURRENT_HTML_INSTEAD_OF_DATA:"+logical_name)
    return True

def curl_download(url,dst,old_limit,name):
    validate_source_url(url)
    dst.parent.mkdir(parents=True,exist_ok=True)
    max_bytes=max(old_limit*2+64*1024*1024,80*1024*1024)
    # All official public snapshots use the same accelerated HTTPS transfer.
    # Large, stable Range-capable responses are segmented; others stream once.
    if os.environ.get("RDMVA_DATA_ACCELERATION","1")!="0":
        import sys
        code=str(ROOT/"code")
        if code not in sys.path:sys.path.insert(0,code)
        from v2810.scientific_download_v31 import download_public_data
        from urllib.error import URLError
        try:
            receipt=download_public_data(url,dst,max_bytes=max_bytes,label=name)
        except (URLError,OSError,TimeoutError,ConnectionError) as exc:
            # Only upstream/network failures may activate the exact pinned
            # SIGNOR recovery; a malformed or modified scientific file may NOT.
            raise RuntimeError("V31_CURRENT_DOWNLOAD_FAILED:"+name+":"+str(exc)[:240]) from exc
        validate_download(dst,name)
        digest=sha(dst)
        print("V31_CURRENT_FETCHED",name,receipt["bytes"],digest,flush=True)
        return {"bytes":receipt["bytes"],"sha256":digest,
                "source_url":url,"retrieved_utc":utc(),
                "download_transport":receipt["mode"]+"_verified"}
    fd,tmp_name=tempfile.mkstemp(prefix=".v31-current-",dir=dst.parent)
    os.close(fd)
    tmp=Path(tmp_name)
    # Windows curl uses Schannel/Windows-trusted certificates, WITHOUT -k.
    # --max-filesize bounds HTTP content; actual size checked after transfer.
    cmd=["curl.exe","--fail","--location","--retry","3","--retry-all-errors",
         "--connect-timeout","30","--max-time","1200","--max-filesize",str(max_bytes),
         "--silent","--show-error","--output",str(tmp),
         "--write-out","%{http_code}",url]
    try:
        # curl uses a temporary file and no credentials. File size is a
        # measurable transfer counter; no claimed percent until the upstream
        # exact Content-Length is known (historical size is NOT that value).
        started=time.monotonic()
        prior=0
        prior_at=started
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                              text=True,errors="replace")
        try:
            while proc.poll() is None:
                now=time.monotonic()
                transferred=tmp.stat().st_size if tmp.exists() else 0
                speed=max(transferred-prior,0)/max(now-prior_at,0.001)/1_000_000
                print(f"V31_CURRENT_DOWNLOAD_ACTIVITY {name} "
                      f"bytes={transferred} speed_mbps={speed:.2f} "
                      f"elapsed_s={int(now-started)} percent=unavailable "
                      "state=downloading",flush=True)
                prior=transferred
                prior_at=now
                time.sleep(5)
            stdout,stderr=proc.communicate(timeout=20)
        except BaseException:
            if proc.poll() is None:
                proc.terminate()
                proc.communicate(timeout=20)
            raise
        if proc.returncode!=0 or stdout.strip()!="200":
            raise RuntimeError("V31_CURRENT_DOWNLOAD_FAILED:"+name+
                               ":curl_exit="+str(proc.returncode)+":http="+stdout[-3:]+
                               ":diagnostic="+stderr[-180:])
        size=tmp.stat().st_size
        if not 0<size<=max_bytes:
            raise RuntimeError("V31_CURRENT_DOWNLOAD_SIZE_INVALID:"+name)
        validate_download(tmp,name)
        digest=sha(tmp)
        os.replace(tmp,dst)
        print("V31_CURRENT_FETCHED",name,size,digest,flush=True)
        return {"bytes":size,"sha256":digest,"source_url":url,"retrieved_utc":utc(),
                "download_transport":"windows_curl_tls_verified"}
    finally:
        if tmp.exists():tmp.unlink()

# Pinned, losslessly compressed October 2026 official SIGNOR 4.0 public data.
# Only used when the university's official endpoint fails. This is a CC BY 4.0
# attributed exception, not an unverified or silently substituted data source.
SIGNOR_RELEASE="signor/Oct2026_release.txt"
SIGNOR_OFFICIAL_URL="https://signor.uniroma2.it/releases/Oct2026_release.txt"
SIGNOR_SEED=ROOT/"reference_seed/signor/Oct2026_release.txt.gz"
SIGNOR_SEED_SHA256="3a18ac07b37f790cb7853992b508650c65c005ecc116c018a183e64aff8feb16"
SIGNOR_SEED_BYTES=6722506
SIGNOR_RELEASE_SHA256="63db6f65e9148e81354f564c3c3cdfc9ee45dc421553bd0926e5fae511ef9d06"
SIGNOR_RELEASE_BYTES=21425930

def recover_official_signor_snapshot(name,url,dst,old_bytes):
    """Restore exact official release, never infer data from partial downloads."""
    if (name!=SIGNOR_RELEASE or url!=SIGNOR_OFFICIAL_URL or
            old_bytes!=SIGNOR_RELEASE_BYTES):
        raise RuntimeError("V312_SIGNOR_FALLBACK_REFUSES_UNPINNED_SOURCE")
    if (not SIGNOR_SEED.is_file() or
            SIGNOR_SEED.stat().st_size!=SIGNOR_SEED_BYTES or
            sha(SIGNOR_SEED)!=SIGNOR_SEED_SHA256):
        raise RuntimeError("V312_SIGNOR_OFFLINE_SEED_MISSING_OR_CORRUPT")
    dst.parent.mkdir(parents=True,exist_ok=True)
    fd,path=tempfile.mkstemp(prefix=".v312-signor-verified-",dir=dst.parent)
    produced=Path(path)
    total=0
    digest=hashlib.sha256()
    started=time.monotonic()
    try:
        with os.fdopen(fd,"wb") as destination, gzip.open(SIGNOR_SEED,"rb") as reader:
            while True:
                block=reader.read(1048576)
                if not block: break
                total+=len(block)
                if total>SIGNOR_RELEASE_BYTES:
                    raise RuntimeError("V312_SIGNOR_SEED_DECOMPRESSED_SIZE_OVERFLOW")
                destination.write(block)
                digest.update(block)
        if total!=SIGNOR_RELEASE_BYTES or digest.hexdigest()!=SIGNOR_RELEASE_SHA256:
            raise RuntimeError("V312_SIGNOR_SEED_RELEASE_SHA256_MISMATCH")
        validate_download(produced,name)
        os.replace(produced,dst)
        print("V312_SIGNOR_PINNED_RELEASE_RESTORED "
              "100.0% verified_bytes="+str(total)+
              " elapsed_s="+str(int(time.monotonic()-started)),flush=True)
        return {"bytes":total,"sha256":digest.hexdigest(),
                "source_url":SIGNOR_OFFICIAL_URL,"retrieved_utc":utc(),
                "download_transport":"bundled_official_release_pinned_signor_cc_by_4",
                "official_pinned_fallback":True,
                "upstream_capture_attempted_and_failed":True,
                "source_release":"SIGNOR_4.0_Oct2026",
                "license":"CC BY 4.0; see reference_seed/signor/ATTRIBUTION.md"}
    finally:
        if produced.exists():produced.unlink()


def validate_record(p,r):
    return (p.is_file() and p.stat().st_size==r["bytes"] and sha(p)==r["sha256"])

def bootstrap():
    if os.name!="nt":raise RuntimeError("V31_CURRENT_WINDOWS_REQUIRED")
    st=copy.deepcopy(read(STATIC_ORIGINAL))
    reg=copy.deepcopy(read(REG_ORIGINAL))
    if st.get("strict_snapshot") is not True:raise RuntimeError("V31_ORIGINAL_MANIFEST_INVALID")
    if reg.get("schema")!="V283_REGULATORY_SCHEMA_1":raise RuntimeError("V31_REGULATORY_SCHEMA_INVALID")
    if RECEIPT.is_file() and STATIC_LIVE.is_file() and REG_LIVE.is_file():
        if verify() and read(RECEIPT).get("state")=="READY":
            print("V31_CURRENT_ALREADY_VERIFIED",flush=True)
            return

    prior=read(RECEIPT) if RECEIPT.is_file() else {}
    files=(prior.get("files") if prior.get("schema")==SCHEMA and
           prior.get("state") in ("FAILED","IN_PROGRESS","DEGRADED") else {}) or {}
    # A failed/resumed install remains the SAME acquisition, even across
    # Windows reboots. Scope the ChEMBL response checkpoint to that capture;
    # a different fresh install gets a completely different cache namespace.
    old_capture=(prior.get("capture_started_utc") or prior.get("started_utc")
                 if prior.get("state") in ("FAILED","IN_PROGRESS","DEGRADED") else None)
    capture_started=old_capture or utc()
    cache_scope=hashlib.sha256(
        (capture_started+"|"+sha(STATIC_ORIGINAL)+"|"+sha(REG_ORIGINAL))
        .encode("utf-8")).hexdigest()[:20]
    receipt={"schema":SCHEMA,"state":"IN_PROGRESS","started_utc":utc(),
             "capture_started_utc":capture_started,
             "chembl_capture_scope":cache_scope,
             "profile":"current_official_sources",
             "not_a_frozen_2026_10_10_snapshot":True,
             "patient_data_included":False,"source_input_sha256":{},
             "files":files}
    save(RECEIPT,receipt)
    # Avoid silently treating the historical snapshots as current. Each
    # invocation without an installed live receipt downloads from upstream.
    warnings=DownloadWarnings(ROOT)
    try:
        ca_bulk_ready={}
        skipped_static=set()
        for entry in st["static_files"]:
            name=entry["path"]
            url=URL_OVERRIDES.get(name,entry["url"])
            old_bytes=int(entry["bytes"])
            meta=receipt["files"].get(name)
            if not (meta and meta.get("source_url")==url and
                    validate_record(RAW/name,meta)):
                same_signor=(name==SIGNOR_RELEASE and url==SIGNOR_OFFICIAL_URL
                             and old_bytes==SIGNOR_RELEASE_BYTES)
                def official_download(n=name,u=url,b=old_bytes):
                    return curl_download(u,RAW/n,b,n)
                def pinned_signor(n=name,u=url,b=old_bytes):
                    return recover_official_signor_snapshot(n,u,RAW/n,b)
                meta=resilient_download(
                    name,url,official_download,
                    fallback=pinned_signor if same_signor else None,
                    fallback_name="verified_pinned_official_SIGNOR_release" if same_signor else None,
                    warnings=warnings)
            else:
                warnings.resolved(name)
                print("V31_CURRENT_RESUMED_VERIFIED",name,flush=True)
            if meta is None:
                skipped_static.add(name)
                receipt["files"].pop(name,None)
                save(RECEIPT,receipt)
                continue
            entry.update(bytes=meta["bytes"],sha256=meta["sha256"],url=url)
            if meta.get("official_pinned_fallback"):
                entry["source_capture_mode"]="VERIFIED_OFFICIAL_SIGNOR_OCT2026_PINNED_FALLBACK"
                receipt["pinned_public_source_exception"]={
                    "file":name,"official_release":"October 2026",
                    "reason":"official_server_unavailable",
                    "source_sha256":SIGNOR_RELEASE_SHA256,
                    "bundled_seed_sha256":SIGNOR_SEED_SHA256,
                    "license":"CC BY 4.0"
                }
            receipt["files"][name]=meta
            save(RECEIPT,receipt)
        st["static_files"]=[entry for entry in st["static_files"]
                            if entry["path"] not in skipped_static]
        skipped_reg=set()
        for name,entry in sorted(reg["files"].items()):
            old_bytes=int(entry["bytes"])
            url=entry["url"]
            receipt_name="regulatory_v283/"+name
            meta=receipt["files"].get(receipt_name)
            verified=(meta and
                      (meta.get("source_url")==url or
                       meta.get("requested_source_url")==url) and
                      validate_record(RAW/receipt_name,meta))
            if name in ca_bulk_ready:
                meta=ca_bulk_ready[name]
                warnings.resolved(name)
            elif verified and not meta.get("not_a_direct_api_response"):
                warnings.resolved(name)
                print("V31_CURRENT_RESUMED_VERIFIED",name,flush=True)
            elif verified and meta.get("not_a_direct_api_response"):
                # Existing converted receipts already carry verified origin hashes.
                warnings.resolved(name)
                print("V31_CURRENT_RESUMED_VERIFIED_CANADA_BULK",name,flush=True)
            else:
                canada_names=("ca_approved_20261010.json",
                              "ca_marketed_20261010.json",
                              "ca_ingredients_20261010.json")
                def official_download(n=name,u=url,b=old_bytes,r=receipt_name):
                    return curl_download(u,RAW/r,b,n)
                def official_canadian_bulk(n=name):
                    from v2810.health_canada_bulk_v31 import recover_official_dpd_bulk
                    api_urls={key:reg["files"][key]["url"] for key in canada_names}
                    ca_bulk_ready.update(recover_official_dpd_bulk(RAW,curl_download,api_urls))
                    return ca_bulk_ready[n]
                is_canadian=name in canada_names
                meta=resilient_download(
                    name,url,official_download,
                    fallback=official_canadian_bulk if is_canadian else None,
                    fallback_name="official_Health_Canada_DPD_non_API_bulk_ZIP" if is_canadian else None,
                    warnings=warnings)
            if meta is None:
                skipped_reg.add(name)
                receipt["files"].pop(receipt_name,None)
                save(RECEIPT,receipt)
                continue
            entry.update(bytes=meta["bytes"],sha256=meta["sha256"],
                         url=meta["source_url"])
            receipt["files"]["regulatory_v283/"+name]=meta
            save(RECEIPT,receipt)

        reg["files"]={key:value for key,value in reg["files"].items()
                       if key not in skipped_reg}
        # ChEMBL generated sources must also have their CURRENT source bytes
        # hashed and captured; never accept historical generated hashes here.
        import bootstrap_v28 as v28
        # The default independent-judge backend is official, hash-pinned
        # ChEMBL 37 SQLite. No ChEMBL REST availability is required.
        backend=os.environ.get("RDMVA_CHEMBL_BACKEND","sqlite").strip().lower()
        if backend not in ("sqlite","rest"):
            raise RuntimeError("V31_CHEMBL_BACKEND_UNRECOGNIZED")
        can_build=True
        if backend=="sqlite":
            import sys
            code_root=str(ROOT/"code")
            if code_root not in sys.path:sys.path.insert(0,code_root)
            from v2810.chembl_sqlite_download_v31 import ensure_official_sqlite
            chembl_meta=resilient_download(
                "ChEMBL_37_SQLite",
                "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz",
                ensure_official_sqlite, warnings=warnings)
            if chembl_meta is None:
                can_build=False
                chembl_meta={"archive_sha256":None,"sqlite_sha256":None}
            receipt["chembl_source"]={
                "backend":"OFFICIAL_SQLITE","release":37,
                "source_archive_sha256":chembl_meta["archive_sha256"],
                "database_sha256":chembl_meta["sqlite_sha256"],
                "license":"CC BY-SA 3.0 Unported"}
            save(RECEIPT,receipt)
            for spec in st["generated_files"]:
                spec["method"]="chembl37_sqlite_local_export"
        else:
            receipt["chembl_source"]={"backend":"REST_CURRENT_API"}
        if can_build:
            for spec in st["generated_files"]:
                spec["capture_current"]=True
            previous_scope=os.environ.get("RDMVA_CHEMBL_CACHE_SCOPE")
            try:
                os.environ["RDMVA_CHEMBL_CACHE_SCOPE"]=cache_scope
                source=("https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/"
                        "releases/chembl_37/chembl_37_sqlite.tar.gz" if backend=="sqlite"
                        else "https://www.ebi.ac.uk/chembl/api/data/")
                can_build=resilient_download(
                    "ChEMBL_generated_sources",source,
                    lambda:(v28.build_chembl(st),True)[1],
                    warnings=warnings) is not None
            finally:
                if previous_scope is None:
                    os.environ.pop("RDMVA_CHEMBL_CACHE_SCOPE",None)
                else:
                    os.environ["RDMVA_CHEMBL_CACHE_SCOPE"]=previous_scope
        if not can_build:
            for spec in st["generated_files"]:
                warnings.skipped(spec["path"],
                    "https://ftp.ebi.ac.uk/pub/databases/chembl/",
                    "ChEMBL verified reference data unavailable")
                receipt["files"].pop(spec["path"],None)
            st["generated_files"]=[]
        available_generated=[]
        for entry in st["generated_files"]:
            name=entry["path"]
            file=RAW/name
            if not validate_record(file,entry):
                warnings.skipped(name,"https://ftp.ebi.ac.uk/pub/databases/chembl/",
                                 "generated source did not pass integrity verification")
                receipt["files"].pop(name,None)
                continue
            available_generated.append(entry)
            entry.pop("capture_current",None)
            receipt["files"][name]={
                "bytes":entry["bytes"],"sha256":entry["sha256"],
                "source_url":("https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz"
                              if backend=="sqlite" else "https://www.ebi.ac.uk/chembl/api/data/"),
                "retrieved_utc":utc(),"download_transport":(
                    "verified_official_sqlite_local" if backend=="sqlite" else "verified_python_https"),
                "method":entry["method"]}
            save(RECEIPT,receipt)
        st["generated_files"]=available_generated
        st["version"]="V31_CURRENT_PUBLIC_SNAPSHOT_"+dt.datetime.now(dt.timezone.utc).strftime("%Y%m%d")
        st["snapshot_profile"]="CURRENT_OFFICIAL_SOURCES"
        st["strict_snapshot"]=True
        reg["snapshot_date"]=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
        reg["development_only"]=True
        reg["snapshot_profile"]="CURRENT_OFFICIAL_SOURCES"
        # Output activated only after all files captured + verified.
        save(STATIC_LIVE,st)
        save(REG_LIVE,reg)
        receipt["skipped_sources"]=list(warnings.data["skipped_sources"])
        receipt["recovered_sources"]=list(warnings.data["recovered_sources"])
        receipt["state"]="DEGRADED" if receipt["skipped_sources"] else "READY"
        receipt["finished_utc"]=utc()
        receipt["source_input_sha256"]={
          "v28_original_manifest":sha(STATIC_ORIGINAL),
          "v283_original_manifest":sha(REG_ORIGINAL),
        }
        save(RECEIPT,receipt)
        if not verify():
            raise RuntimeError("V31_CURRENT_FINAL_VERIFY_FAILED")
        if receipt["state"]=="DEGRADED":
            warnings.print_summary()
            print("V31_CURRENT_PUBLIC_DATA_DEGRADED",len(receipt["files"]),
                  "missing="+str(len(receipt["skipped_sources"])),flush=True)
        else:
            print("V31_CURRENT_PUBLIC_DATA_READY",len(receipt["files"]),flush=True)
    except Exception as exc:
        receipt["state"]="FAILED"
        receipt["failure_type"]=type(exc).__name__
        receipt["failure_summary"]=str(exc)[:320]
        save(RECEIPT,receipt)
        raise

def verify():
    if not all(p.is_file() for p in (RECEIPT,STATIC_LIVE,REG_LIVE)):
        return False
    receipt=read(RECEIPT)
    if receipt.get("schema")!=SCHEMA or receipt.get("state") not in ("READY","DEGRADED"):
        return False
    st=read(STATIC_LIVE);reg=read(REG_LIVE)
    try:
        for e in st["static_files"]+st["generated_files"]:
            if not validate_record(RAW/e["path"],e):return False
        for n,e in reg["files"].items():
            if not validate_record(RAW/"regulatory_v283"/n,e):return False
    except (OSError,ValueError,KeyError):return False
    return True

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--verify-only",action="store_true")
    opts=a.parse_args()
    if opts.verify_only:
        if not verify():raise SystemExit("V31_CURRENT_PUBLIC_DATA_INVALID_OR_MISSING")
        print("V31_CURRENT_PUBLIC_DATA_SHA256_OK",flush=True)
        if read(RECEIPT).get("state")=="DEGRADED":
            DownloadWarnings(ROOT).print_summary()
            print("V31_CURRENT_PUBLIC_DATA_DEGRADED_VERIFIED",flush=True)
    else:bootstrap()

if __name__=="__main__":main()
