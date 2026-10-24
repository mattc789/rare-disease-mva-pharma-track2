"""Download and verify the exact official ChEMBL 37 SQLite database.

Source: EMBL-EBI ChEMBL FTP. Data are CC BY-SA 3.0; never bundled in
the ExaMedix source distribution. Install on the customer's own T:/etc drive.
The 5.76 GB archive is resumable, hash pinned and not silently replaced.
"""
from __future__ import annotations
import datetime,hashlib,json,os,shutil,sqlite3,subprocess,tarfile,urllib.parse,contextlib,time
from pathlib import Path

BASE="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/"
ARCHIVE_NAME="chembl_37_sqlite.tar.gz"
ARCHIVE_SHA256="33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281"
ARCHIVE_BYTES=5764252857
SCHEMA="RDMVA_CHEMBL_SQLITE_OFFICIAL_V31_1"
ROOT=Path(__file__).resolve().parents[2]/"public_inputs/third_party/chembl37"
DB=ROOT/"chembl_37.db"
AUDIT=ROOT/"CHEMBL_SQLITE_PROVENANCE_V31.json"
REQUIRED={
  "molecule_dictionary":{"molregno","chembl_id","max_phase","pref_name","molecule_type"},
  "molecule_hierarchy":{"molregno","parent_molregno","active_molregno"},
  "molecule_synonyms":{"molregno","synonyms","syn_type"},
  "molecule_atc_classification":{"molregno","level5"},
  "drug_mechanism":{"mec_id","molregno","record_id","tid","action_type"},
  "compound_records":{"record_id","molregno"},
  "target_dictionary":{"tid","chembl_id","organism","tax_id","target_type"},
  "target_components":{"tid","component_id"},
  "component_sequences":{"component_id","accession"},
  "activities":{"activity_id","assay_id","molregno","standard_value","standard_type"},
  "assays":{"assay_id","chembl_id","tid","confidence_score","assay_type"},
  "compound_structures":{"molregno","canonical_smiles"},
  "compound_properties":{"molregno","full_mwt","mw_freebase"},
  "docs":{"doc_id","chembl_id"},
}
def sha(path):
    path=Path(path)
    h=hashlib.sha256()
    size=path.stat().st_size
    show_progress=size>=500_000_000
    seen=0
    last_report=time.monotonic()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""):
            h.update(block)
            if show_progress:
                seen+=len(block)
                now=time.monotonic()
                if now-last_report>=5:
                    print(f"V31_CH37_SHA256_PROGRESS {path.name} "
                          f"{(100*seen/size):.1f}% "
                          f"{seen/1_000_000_000:.2f}/{size/1_000_000_000:.2f} GB",
                          flush=True)
                    last_report=now
    if show_progress:
        print(f"V31_CH37_SHA256_PROGRESS {path.name} 100.0% "
              "state=digest_computed",flush=True)
    return h.hexdigest()

def _curl(url,destination,*,resume=False,timeout=180):
    p=urllib.parse.urlsplit(url)
    if p.scheme!="https" or p.hostname!="ftp.ebi.ac.uk" or not p.path.startswith("/pub/databases/chembl/ChEMBLdb/releases/chembl_37/"):
        raise RuntimeError("V31_CH37_UNAPPROVED_OFFICIAL_SOURCE")
    args=["curl.exe","--fail","--location","--retry","5","--retry-delay","5",
          "--retry-all-errors","--connect-timeout","30","--max-time",str(timeout),
          "--silent","--show-error","--output",str(destination)]
    if resume:args+=["--continue-at","-"]
    args.append(url)
    print("V31_CH37_OFFICIAL_DOWNLOAD",Path(p.path).name,
          "resuming" if resume else "official",flush=True)
    # ChEMBL is a 5.76 GB download: emit real, monotonic-file-size progress.
    # Do not infer completion from the printed percentage: checksum validation
    # and SQLite schema validation remain mandatory after curl finishes.
    is_archive=(Path(p.path).name==ARCHIVE_NAME and
                Path(destination).name==ARCHIVE_NAME+".part")
    if not is_archive:
        proc=subprocess.run(args,check=False)
        if proc.returncode:
            raise RuntimeError("V31_CH37_VERIFIED_CURL_FAILED:"+str(proc.returncode))
        return
    print("V31_CH37_DOWNLOAD_NOTICE: The 5.76 GB ChEMBL 37 download and its "
          "verification/extraction may take over 60 minutes (sometimes hours), "
          "depending on connection, disk and workstation speed. "
          "Do not interrupt; resume is supported.",flush=True)
    import time
    dest=Path(destination)
    previous_bytes=dest.stat().st_size if dest.exists() else 0
    previous_at=time.monotonic()
    last_change=previous_at
    proc=subprocess.Popen(args)
    try:
        while proc.poll() is None:
            now=time.monotonic()
            size=dest.stat().st_size if dest.exists() else 0
            delta_seconds=max(now-previous_at,0.001)
            mbps=max(size-previous_bytes,0)/delta_seconds/1_000_000
            if size!=previous_bytes:
                last_change=now
            pct=min(100.0,max(0.0,100.0*size/ARCHIVE_BYTES))
            state="downloading" if now-last_change<30 else "waiting_for_data"
            print(f"V31_CH37_DOWNLOAD_PROGRESS {pct:.1f}% "
                  f"{size/1_000_000_000:.2f}/{ARCHIVE_BYTES/1_000_000_000:.2f} GB "
                  f"{mbps:.2f} MB/s state={state}",flush=True)
            previous_bytes=size
            previous_at=now
            time.sleep(5)
        if proc.returncode:
            raise RuntimeError("V31_CH37_VERIFIED_CURL_FAILED:"+str(proc.returncode))
        final_bytes=dest.stat().st_size if dest.exists() else 0
        print(f"V31_CH37_DOWNLOAD_PROGRESS "
              f"{min(100.0,100.0*final_bytes/ARCHIVE_BYTES):.1f}% "
              f"{final_bytes/1_000_000_000:.2f}/{ARCHIVE_BYTES/1_000_000_000:.2f} GB "
              "state=download_finished_verifying_SHA256",flush=True)
    except BaseException:
        if proc.poll() is None:
            proc.terminate()
            proc.wait()
        raise

def _sqlite_quick_check_with_activity(db, *, interval_s=10):
    """Keep the console responsive while SQLite performs its *real* check.

    SQLite does not expose a reliable total-work or page-completion percent for
    PRAGMA quick_check. Report observed VM instructions and elapsed time, NEVER
    a guessed percent or a 'pass' until SQLite actually returns 'ok'.
    """
    from threading import Event, Thread
    started=time.monotonic()
    finished=Event()
    instructions=[0]
    stride=10000
    def counted():
        instructions[0]+=stride
        return 0
    def heartbeat():
        while not finished.wait(interval_s):
            print("V31_CH37_SQLITE_QUICK_CHECK_ACTIVITY "
                  f"elapsed_s={int(time.monotonic()-started)} "
                  f"sqlite_vm_steps={instructions[0]} "
                  "percent=unavailable state=checking",flush=True)
    db.set_progress_handler(counted,stride)
    watcher=Thread(target=heartbeat,daemon=True)
    watcher.start()
    try:
        return db.execute("PRAGMA quick_check(1)").fetchone()
    finally:
        finished.set()
        db.set_progress_handler(None,0)
        watcher.join(timeout=2)


def _check_schema(path):
    if not path.is_file() or path.stat().st_size<1_000_000_000:
        raise RuntimeError("V31_CH37_SQLITE_TOO_SMALL")
    uri=path.as_uri()+"?mode=ro&immutable=1"
    with contextlib.closing(sqlite3.connect(uri,uri=True,timeout=90)) as db:
        db.execute("PRAGMA query_only=ON")
        print("V31_CH37_SQLITE_QUICK_CHECK_START "
              "progress=indeterminate heartbeat_every_s=10",flush=True)
        status=_sqlite_quick_check_with_activity(db)
        if not status or status[0]!="ok":
            raise RuntimeError("V31_CH37_SQLITE_QUICK_CHECK_FAILED")
        print("V31_CH37_SQLITE_QUICK_CHECK_PASS",flush=True)
        for table,columns in REQUIRED.items():
            available={r[1].lower() for r in db.execute("PRAGMA table_info('"+table+"')")}
            if not columns<=available:
                raise RuntimeError("V31_CH37_SQLITE_SCHEMA_INCOMPLETE:"+table+
                                   ":"+",".join(sorted(columns-available)))
        if not db.execute("SELECT 1 FROM molecule_dictionary WHERE max_phase=4 LIMIT 1").fetchone():
            raise RuntimeError("V31_CH37_NO_APPROVED_DRUGS")
        if not db.execute("SELECT 1 FROM drug_mechanism LIMIT 1").fetchone():
            raise RuntimeError("V31_CH37_NO_DRUG_MECHANISMS")
        if not db.execute("SELECT 1 FROM activities LIMIT 1").fetchone():
            raise RuntimeError("V31_CH37_NO_BIOACTIVITY")

def _read_archive():
    archive=ROOT/ARCHIVE_NAME
    part=ROOT/(ARCHIVE_NAME+".part")
    if archive.is_file() and archive.stat().st_size==ARCHIVE_BYTES and sha(archive)==ARCHIVE_SHA256:
        print("V31_CH37_OFFICIAL_ARCHIVE_VERIFIED_CACHED",flush=True)
        return archive
    if part.exists() and part.stat().st_size>ARCHIVE_BYTES:
        raise RuntimeError("V31_CH37_PARTIAL_ARCHIVE_EXCEEDS_OFFICIAL_LENGTH")
    if part.exists() and part.stat().st_size==ARCHIVE_BYTES:
        if sha(part)!=ARCHIVE_SHA256:
            raise RuntimeError("V31_CH37_FULL_PARTIAL_ARCHIVE_SHA_MISMATCH")
    else:
        if part.exists() and part.stat().st_size>0:
            # Legacy sequential partial has no source-validator checkpoint.
            # Keep the existing resumable curl path; SHA-256 remains mandatory.
            if not shutil.which("curl.exe"):
                raise RuntimeError("V31_CH37_WINDOWS_CURL_NOT_AVAILABLE")
            _curl(BASE+ARCHIVE_NAME,part,resume=True,timeout=0)
        else:
            from .scientific_download_v31 import download_public_data
            download_public_data(BASE+ARCHIVE_NAME,part,
                                 expected_bytes=ARCHIVE_BYTES,
                                 expected_sha256=ARCHIVE_SHA256,
                                 label=ARCHIVE_NAME)
    if part.stat().st_size!=ARCHIVE_BYTES or sha(part)!=ARCHIVE_SHA256:
        raise RuntimeError("V31_CH37_OFFICIAL_ARCHIVE_SHA256_OR_LENGTH_MISMATCH")
    os.replace(part,archive)
    print("V31_CH37_OFFICIAL_ARCHIVE_SHA256_PASS",archive.stat().st_size,flush=True)
    return archive

def _extract_single_verified_db(archive):
    stage=ROOT/("chembl_37.db"+".partial")
    if stage.exists():
        print("V31_CH37_DISCARD_INCOMPLETE_LOCAL_DB_EXTRACTION",flush=True)
        stage.unlink()
    digest=hashlib.sha256()
    found=0
    try:
        with tarfile.open(archive,"r|gz") as z:
            for member in z:
                if not member.isfile() or not member.name.lower().endswith((".db",".sqlite",".sqlite3")):
                    continue
                if found:raise RuntimeError("V31_CH37_MULTIPLE_DATABASE_MEMBERS")
                if not 1_000_000_000<=member.size<=100_000_000_000:
                    raise RuntimeError("V31_CH37_UNEXPECTED_DB_MEMBER_SIZE")
                source=z.extractfile(member)
                if source is None:raise RuntimeError("V31_CH37_DB_MEMBER_UNREADABLE")
                found+=1
                total=0
                last_report=time.monotonic()
                print(f"V31_CH37_DB_EXTRACTION_START {member.size/1_000_000_000:.2f} GB",flush=True)
                with source,stage.open("wb") as dest:
                    while True:
                        chunk=source.read(8*1024*1024)
                        if not chunk:break
                        total+=len(chunk)
                        if total>member.size:raise RuntimeError("V31_CH37_EXTRACTED_DB_OVERSIZED")
                        dest.write(chunk);digest.update(chunk)
                        now=time.monotonic()
                        if now-last_report>=5:
                            print(f"V31_CH37_DB_EXTRACTION_PROGRESS "
                                  f"{100*total/member.size:.1f}% "
                                  f"{total/1_000_000_000:.2f}/{member.size/1_000_000_000:.2f} GB",
                                  flush=True)
                            last_report=now
                if total!=member.size:raise RuntimeError("V31_CH37_DB_MEMBER_TRUNCATED")
                print("V31_CH37_DB_EXTRACTION_PROGRESS 100.0% state=checking_database",flush=True)
        if found!=1:raise RuntimeError("V31_CH37_DB_MEMBER_NOT_FOUND")
        _check_schema(stage)
        os.replace(stage,DB)
        return digest.hexdigest()
    except Exception:
        if stage.exists():stage.unlink()
        raise

def ensure_official_sqlite():
    ROOT.mkdir(parents=True,exist_ok=True)
    if AUDIT.is_file() and DB.is_file():
        try:
            metadata=json.loads(AUDIT.read_text(encoding="utf-8"))
            if (metadata.get("schema")==SCHEMA and
                metadata.get("chembl_release")==37 and
                metadata.get("archive_sha256")==ARCHIVE_SHA256 and
                metadata.get("database_bytes")==DB.stat().st_size and
                metadata.get("sqlite_sha256")==sha(DB)):
                print("V31_CH37_OFFICIAL_SQLITE_ALREADY_VERIFIED",flush=True)
                return metadata
        except (ValueError,KeyError,OSError):
            pass
    archive=_read_archive()
    if DB.is_file() and AUDIT.is_file():
        # Invalid existing database is never silently overwritten.
        raise RuntimeError("V31_CH37_EXISTING_DB_AUDIT_CONFLICT")
    digest=_extract_single_verified_db(archive)
    for name in ("README","LICENSE","checksums.txt"):
        dst=ROOT/name
        _curl(BASE+name,dst,timeout=75)
        if dst.stat().st_size<100:raise RuntimeError("V31_CH37_LEGAL_METADATA_INCOMPLETE:"+name)
    checks=(ROOT/"checksums.txt").read_text(encoding="utf-8")
    if ARCHIVE_NAME not in checks or ARCHIVE_SHA256 not in checks:
        raise RuntimeError("V31_CH37_PUBLISHED_CHECKSUMS_MISMATCH")
    result={
      "schema":SCHEMA,"chembl_release":37,
      "archive_source_url":BASE+ARCHIVE_NAME,"archive_bytes":ARCHIVE_BYTES,
      "archive_sha256":ARCHIVE_SHA256,
      "sqlite_sha256":digest,"database_bytes":DB.stat().st_size,
      "retrieved_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "license":"CC BY-SA 3.0 Unported","data_provider":"ChEMBL EMBL-EBI",
      "license_source_url":BASE+"LICENSE","references_downloaded":True,
      "database_is_external_not_relicensed_under_app_license":True,
      "source_run_patient_data_included":False
    }
    temp=AUDIT.with_suffix(".json.tmp")
    temp.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    os.replace(temp,AUDIT)
    print("V31_CH37_LOCAL_SQLITE_READY",DB.stat().st_size,flush=True)
    return result

if __name__=="__main__":ensure_official_sqlite()
