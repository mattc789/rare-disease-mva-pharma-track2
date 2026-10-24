"""Build deterministic Rare Disease MVA Pharma V3.1 source-only release.

Strict whitelist: original source/tests/docs plus public metadata/policy only.
Patient/gated data, generated runs, bulk third-party data, secrets, API caches,
model weights, external binaries and local machine configuration are excluded.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,json,re,zipfile

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"release/RARE_DISEASE_MVA_PHARMA_V31_SOURCE_ONLY.zip"
REPORT=ROOT/"release/V31_SOURCE_RELEASE_MANIFEST.json"

FILES=[
 ".gitignore",".gitattributes","README.md","LICENSE",
 "THIRD_PARTY_NOTICES_V31.md","requirements-v2810.txt",
 "reference_seed/signor/Oct2026_release.txt.gz","reference_seed/signor/ATTRIBUTION.md",
 "bootstrap_v28.py","bootstrap_regulatory_v283.py","bootstrap_current_public_v31.py","setup_v2810.py",
 "run_v31.ps1","REPRODUCE_CHALLENGE_V31.ps1","INSTALL_WINDOWS_V31.ps1","CONFIGURE_STORAGE_V31.ps1","SETUP_OPENAI_KEY_V31.ps1",
 "JURY_START_V31.cmd","JURY_START_V31.ps1","JURY_MONITOR_V31.py","JURY_CHECK_AI_V31.py",
 "reproduce_v31_pipeline.py","verify_challenge_run_v31.py",
 "build_source_release_v31.py",
 "config/public_data_manifest_v28.json",
 "config/regulatory_snapshot_v283.json",
 "config/curated_human_pk_v286.example.json",
 "config/curated_population_af_v288.example.json",
 "config/local_public_storage_v2810.example.json",
 "config/source_registry.json",
 "config/openai_v30.json",
 "config/v30_source_policy.json",
 "config/v31_curated_quantitative_evidence.json",
 "submission/examedix-technologies_track2_report.md",
 "submission/examedix-technologies_track2_report.pdf",
 "submission/VIDEO_PITCH.md",
]
for base in (ROOT/"code/v2810",ROOT/"code/v30"):
    FILES += [p.relative_to(ROOT).as_posix() for p in sorted(base.glob("*.py"))]
FILES += [p.relative_to(ROOT).as_posix()
          for p in sorted((ROOT/"code/v2810/externals").glob("*.py"))]
FILES += [p.relative_to(ROOT).as_posix()
          for p in sorted((ROOT/"code/v2810/externals").glob("*.sh"))]
FILES += [p.relative_to(ROOT).as_posix()
          for p in sorted((ROOT/"code").glob("test_v2810_*.py"))]
FILES += ["code/test_v31_signor_official_fallback.py","code/test_v30_engine.py","code/test_v31_cleanroom_runner.py","code/test_v31_live_public_data.py","code/test_v31_auto_storage.py","code/test_v31_live_universe_resume.py","code/test_v31_large_resource_recovery.py","code/test_v31_interrupted_run_recovery.py","code/test_v31_required_genomic_evidence.py","code/test_v31_fasta_reuse.py","code/test_v31_jury_launcher.py","code/test_v31_jury_ai.py","code/test_v31_fresh_jury_install.py","code/test_v31_chembl_resilience.py","code/test_v31_powershell_reserved_vars.py","code/test_v31_chembl_sqlite.py",
         "code/test_v31_isolation_progress.py","code/test_v31_scientific_download.py",
         "code/test_v31_jury_input_auto.py","code/test_v31_local_gnomad_manifest.py",
          "code/test_v31_download_resilience.py"]

DENY_DIRS={
 "env","experiments","public_inputs","runs","frozen","state","tools","audits",
 "private","dataset","genomes","model_weights","rdmva_resources",".git",
 "__pycache__","release","work"}
DENY_SUFFIX={
 ".vcf",".bcf",".bam",".cram",".fasta",".fa",".fq",".fastq",".gz",".bgz",
 ".zip",".db",".sqlite",".pickle",".pkl",".pth",".pt",".onnx",".h5",".hdf5",
 ".ckpt",".parquet",".xlsx",".xls",".exe",".dll",".bin",".dpapi"}

FORBIDDEN_TEXT=[
 (re.compile(re.escape("E:"+chr(92)+"Rare Disease MVA Pharma"),re.I),
  "development_project_absolute_path"),
 (re.compile(re.escape("/home/"+"mc789"),re.I),
  "development_wsl_username"),
 (re.compile(re.escape("TOSHIBA"+" EXTERNAL_USB"),re.I),
  "development_external_disk_identity"),
 (re.compile(re.escape("Big"+"Main"),re.I),
  "development_workstation_name"),
 (re.compile(r"(?i)\bWGS[_-]EX\d{7,}\b"),"patient_sample_identifier"),
 (re.compile(r"(?i)\bChallenge_Clinical_Phenotype_\d+\b"),"specific_clinical_file_identifier"),
 (re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),"possible_api_key"),
 (re.compile(r"(?i)(api[_-]?key|secret|token)\s*[:=]\s*['\"][A-Za-z0-9_./+\-=]{16,}['\"]"),
  "possible_embedded_secret"),
]
TEXT_SUFFIXES=(".py",".ps1",".cmd",".json",".txt",".md",".sh",".gitignore",".gitattributes")
FIXED_ZIP_TIME=(1980,1,1,0,0,0)
SIGNED_SIGNOR_SEED="reference_seed/signor/Oct2026_release.txt.gz"
SIGNED_SIGNOR_SEED_SHA256="3a18ac07b37f790cb7853992b508650c65c005ecc116c018a183e64aff8feb16"
SIGNED_SIGNOR_SEED_BYTES=6722506
PUBLIC_REPORT_PDF="submission/examedix-technologies_track2_report.pdf"
PUBLIC_REPORT_PDF_SHA256="7178f6403ffac36af511a4e216f0ef0ee51712753691ce62f1d17fb47e12c501"

def canonical_bytes(rel,p):
    b=p.read_bytes()
    if rel.endswith(TEXT_SUFFIXES):
        text=b.decode("utf-8-sig").replace("\r\n","\n").replace("\r","\n")
        if rel.lower().endswith(".cmd"):
            return text.replace("\n","\r\n").encode("utf-8")
        return text.encode("utf-8")
    return b

def sha(blob):return hashlib.sha256(blob).hexdigest()

def _zip_write(z,name,data):
    info=zipfile.ZipInfo(name,FIXED_ZIP_TIME)
    info.compress_type=zipfile.ZIP_DEFLATED
    info.create_system=3
    info.external_attr=(0o100644 & 0xFFFF)<<16
    z.writestr(info,data,compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)

def _validate(rel,p,b):
    # Sole permitted binary reference exception: the exact CC BY 4.0
    # compressed official SIGNOR Oct2026 release, not arbitrary .gz files.
    if rel==SIGNED_SIGNOR_SEED:
        if len(b)!=SIGNED_SIGNOR_SEED_BYTES or sha(b)!=SIGNED_SIGNOR_SEED_SHA256:
            raise RuntimeError("V312_PUBLIC_SIGNOR_SEED_PIN_MISMATCH")
        return
    # One specifically reviewed public PDF report; never allow arbitrary binary files.
    if rel==PUBLIC_REPORT_PDF:
        if (len(b)>1_000_000 or not b.startswith(b"%PDF-") or
            not b.rstrip().endswith(b"%%EOF") or sha(b)!=PUBLIC_REPORT_PDF_SHA256):
            raise RuntimeError("V31_PUBLIC_REPORT_PDF_INTEGRITY_REVIEW_REQUIRED")
        return
    if any(part.lower() in DENY_DIRS for part in Path(rel).parts):
        raise RuntimeError("V31_DENIED_PATH_IN_ARCHIVE:"+rel)
    if p.suffix.lower() in DENY_SUFFIX:
        raise RuntimeError("V31_DENIED_FILE_TYPE_IN_ARCHIVE:"+rel)
    if p.stat().st_size>2_000_000:
        raise RuntimeError("V31_UNEXPECTED_LARGE_SOURCE_FILE:"+rel)
    if ".." in Path(rel).parts or p.is_symlink():
        raise RuntimeError("V31_ARCHIVE_PATH_ESCAPE:"+rel)
    if b"\x00" in b:
        raise RuntimeError("V31_BINARY_CONTENT_DENIED:"+rel)
    if rel.endswith(TEXT_SUFFIXES):
        text=b.decode("utf-8-sig")
        for rx,label in FORBIDDEN_TEXT:
            if rx.search(text):
                raise RuntimeError("V31_NONPORTABLE_OR_SECRET_SOURCE:"+label+":"+rel)

def build():
    OUT.parent.mkdir(parents=True,exist_ok=True)
    items=[];seen=set()
    for rel in sorted(FILES):
        rel=rel.replace("\\","/")
        if rel in seen:continue
        seen.add(rel);p=ROOT/rel
        if not p.is_file():
            raise FileNotFoundError("Required public source file missing: "+rel)
        b=canonical_bytes(rel,p);_validate(rel,p,b)
        items.append({"path":rel,"bytes":len(b),"sha256":sha(b),"data":b})
    public_manifest={
      "schema":"RDMVA_V31_SOURCE_ONLY_MANIFEST_1",
      "version":"3.1.4",
      "copyright":"Matteo Caprani | ExaMedix Technologies",
      "original_application_license":"CC-BY-4.0",
      "third_party_terms_preserved":True,
      "contains_patient_raw_genome":False,
      "contains_private_variant_records":False,
      "contains_gated_challenge_data":False,
      "contains_public_bulk_reference_database":True,
      "only_public_bulk_dataset":"SIGNOR_4.0_Oct2026_6.72MB_gzip_CC_BY_4.0",
      "bundled_public_signor_official_gzip_sha256":SIGNED_SIGNOR_SEED_SHA256,
      "contains_external_binaries_or_model_weights":False,
      "contains_local_machine_configuration":False,
      "contains_api_keys_or_secrets":False,
      "contains_api_cache":False,
      "contains_case_specific_research_outputs":False,
      "contains_source_code_tests_docs_public_metadata_only":True,
      "alphamissense_default_profile_enabled":False,
      "deterministic_archive":True,
      "files":[{k:x[k] for k in ("path","bytes","sha256")} for x in items],
    }
    content=(json.dumps(public_manifest,indent=2,ensure_ascii=False,
                        sort_keys=True)+"\n").encode("utf-8")
    tmp=OUT.with_suffix(OUT.suffix+".tmp")
    if tmp.exists():tmp.unlink()
    with zipfile.ZipFile(tmp,"w",allowZip64=False) as z:
        _zip_write(z,"SOURCE_MANIFEST.json",content)
        for item in items:_zip_write(z,item["path"],item["data"])
    with zipfile.ZipFile(tmp) as z:
        bad=z.testzip()
        if bad:raise RuntimeError("V31_ZIP_CRC_FAILURE:"+bad)
        expected={x["path"] for x in items}|{"SOURCE_MANIFEST.json"}
        if set(z.namelist())!=expected:
            raise RuntimeError("V31_ZIP_WHITELIST_MISMATCH")
        for item in items:
            if sha(z.read(item["path"]))!=item["sha256"]:
                raise RuntimeError("V31_ZIP_HASH_MISMATCH:"+item["path"])
    tmp.replace(OUT)
    REPORT.write_bytes(content)
    print("V31_PUBLIC_SOURCE_RELEASE_CREATED",OUT)
    print("FILES",len(items),"ZIP_BYTES",OUT.stat().st_size)
    print("PATIENT_DATA_BUNDLED",False)
    print("LICENSE_ATTRIBUTED_SIGNOR_PUBLIC_DATA_BUNDLED",True)
    print("API_KEYS_OR_CACHES_BUNDLED",False)
    print("DETERMINISTIC_SHA256",sha(OUT.read_bytes()))

if __name__=="__main__":
    build()
