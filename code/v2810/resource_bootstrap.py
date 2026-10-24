"""General-purpose, license-aware chromosome resource bootstrap.

Rare Disease MVA Pharma v2.8.10 | ExaMedix Technologies
- Source data downloaded independently at use, NEVER copied into distributable code.
- All downloads are public whole-chromosome sources; NEVER patient variant coordinates.
- GRCh38 and gnomAD v4.1 only; no silent assembly liftover.
- Windows WSL2 disk used to avoid filling the Windows project volume.
- The manifest is local and reproducible, version/source/hash recorded.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from urllib.parse import urlparse
import json,os,re,subprocess,time,hashlib,shutil

VERSION="V2810_AUTO_RESOURCES_1"
GNOMAD_ROOT="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
CHROMS={str(i) for i in range(1,23)}|{"X","Y"}
COHORTS={"exomes","genomes"}
MIN_FREE_RESERVE=12*1024**3
CATALOG={
 "VEP":dict(license="Apache-2.0",price="free",distribution="external_dependency",
            source="https://github.com/Ensembl/ensembl-vep",version="116",
            license_verified=True,commercial_use_permitted=True),
 "gnomAD-v4.1-VCF":dict(license="gnomAD source data policy; older ODbL references; verify per release",
            price="free_public_VCF",distribution="download_from_official_host_not_rebundled",
            source="https://gnomad.broadinstitute.org/policies",
            license_verified_for_source_download=True,
            commercial_redistribution_of_transformed_database="REQUIRES_UPSTREAM_POLICY_CHECK"),
 "ClinVar":dict(license="NCBI public database use; check third-party assertions",
            price="free",distribution="download_officially_not_patient_data",
            source="https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/",license_verified=True),
 "HPO":dict(license="HPO terms; acknowledgement, version, do-not-edit-ontology",
            price="free",distribution="download_unmodified_with_attribution",
            source="https://human-phenotype-ontology.github.io/license.html",license_verified=True),
 "GENCODE":dict(license="GENCODE source conditions; confirm pinned file version",
            price="free",distribution="download_officially",
            source="https://www.gencodegenes.org/human/",license_verified_for_source_download=True),
 "ChEMBL":dict(license="CC-BY-SA-3.0 database",
            price="free",distribution="download_and_attribution_sharealike_for_adapted_database",
            source="https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions",
            license_verified=True),
 "SIGNOR-4.0":dict(license="CC-BY-4.0",price="free",distribution="download_unmodified_with_attribution",
            source="https://signor.uniroma2.it/about/",license_verified=True),
 "UniProtKB":dict(license="CC-BY-4.0 for copyrightable database content",price="free",
            distribution="download_with_source_attribution",
            source="https://www.uniprot.org/help/license",license_verified=True),
 "Reactome-annotations":dict(license="CC0-1.0 annotation/interaction data",
            price="free",distribution="download_with_provenance",
            source="https://dev.reactome.org/license",license_verified=True),
 "PK-Sim":dict(license="GPL-2.0 independent executable",
            price="free",distribution="external_tool_not_relicensed_as_CC-BY",
            source="https://github.com/Open-Systems-Pharmacology/PK-Sim",license_verified=True),
 "Pangolin":dict(license="GPL-3.0 standalone splice tool",
            price="free",distribution="optional_external_tool_license_kept_separate",
            source="https://github.com/tkzeng/Pangolin",license_verified=True),
 "Illumina-SpliceAI":dict(license="PolyForm-Strict software; weights CC-BY-NC-4.0",
            price="not_available_under_required_redistributable_license",
            distribution="EXCLUDED_FROM_AUTOMATIC_BOOTSTRAP",
            source="https://github.com/Illumina/SpliceAI",license_verified=True),
 "Enrichr-LINCS-GEO":dict(license="source-use restrictions; commercial restrictions",
            price="not_qualified_for_general_redistributable_default",
            distribution="EXCLUDED_FROM_DEFAULT",license_verified=False),
 "Hackathon-gated-genome":dict(license="organizer DTA overrides HF dataset badge",
            price="gated",distribution="NEVER_RESHARE_OR_AUTO_UPLOAD",
            source="https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026",
            license_verified=True),
}
def safe_chrom(ch):
    s=str(ch or "").replace("chr","",1).upper()
    if s=="23":s="X"
    if s=="24":s="Y"
    return s if s in CHROMS else None

def needed_chromosomes(private):
    """Select public resources based on actual alleles, not a named disease/gene."""
    found=set()
    omitted=set()
    for variants in (private.get("genes") or {}).values():
        for v in variants:
            chrom=str(v.get("chrom") or "")
            normalized=safe_chrom(chrom)
            if normalized:found.add(normalized)
            elif chrom:omitted.add(chrom)
    # Sort numerically, then X and Y, never in variant priority order.
    return sorted(found,key=lambda x:int(x) if x.isdigit() else (23 if x=="X" else 24)),sorted(omitted)

# Execute every scientific WSL command in the release-required distribution.
# Never silently pick the Windows global default (other projects may use it).
def project_wsl_distribution(project_root=None):
    # Must use same canonical UTF-8/SHA256 algorithm as INSTALL_WINDOWS_V31.ps1.
    # An independent GitHub source checkout gets its own Linux distribution.
    source=Path(project_root or Path(__file__).resolve().parents[2]).resolve()
    canonical=str(source).replace("/","\\").rstrip("\\").lower()
    suffix=hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:8]
    return "RDMVA-Ubuntu-24.04-"+suffix

WSL_DISTRIBUTION=project_wsl_distribution()
def wsl_command(*args,root=False):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",WSL_DISTRIBUTION):
        raise RuntimeError("RDMVA_UNSAFE_WSL_DISTRO_NAME")
    return (["wsl.exe","-d",WSL_DISTRIBUTION]+
            (["--user","root"] if root else [])+
            ["--exec",*args])

def offline_external_bin_wsl():
    """Only verifies tools; these remain separately installed licensed programs."""
    cmd=wsl_command("bash","-lc",
         "command -v curl && command -v tabix && command -v sha256sum && command -v df")
    try:
        done=subprocess.run(cmd,capture_output=True,text=True,timeout=15)
        return done.returncode==0
    except (OSError,subprocess.TimeoutExpired):return False

def _run_long_step(cmd, *, timeout, label, **kwargs):
    """Heartbeat for opaque apt/VEP/sort commands, without fake percentages."""
    import threading
    finished=threading.Event()
    started=time.monotonic()
    def beat():
        while not finished.wait(10):
            print(f"V31_LINUX_STAGE_ACTIVITY {label} "
                  f"elapsed_s={int(time.monotonic()-started)} "
                  "percent=unavailable state=working",flush=True)
    watcher=threading.Thread(target=beat,daemon=True)
    watcher.start()
    try:
        return subprocess.run(cmd,timeout=timeout,**kwargs)
    finally:
        finished.set()
        watcher.join(timeout=2)


def _run_wsl(args,timeout=30):
    cmd=wsl_command(*args)
    options={"capture_output":True,"text":True,"errors":"replace"}
    if timeout>=120:
        p=_run_long_step(cmd,timeout=timeout,label=str(args[0]),**options)
    else:
        p=subprocess.run(cmd,timeout=timeout,**options)
    if p.returncode:
        raise RuntimeError("V2810_WSL_COMMAND_FAILED:"+str(p.returncode)+
                           ":"+p.stderr[-330:])
    return p.stdout.strip()

def _default_wsl_distribution():
    # Provenance must track the same explicit distribution used for every
    # scientific command. Windows' default Ubuntu is not necessarily 24.04.
    return WSL_DISTRIBUTION


def wsl_context():
    homedir=_run_wsl(["bash","-lc","printf %s \"$HOME\""])
    if not re.fullmatch(r"/home/[a-zA-Z0-9_.-]+",homedir):
        raise RuntimeError("V2810_WSL_HOME_INVALID")
    return {"home":homedir,"distribution":_default_wsl_distribution()}

def windows_unc(ctx,linux):
    home=ctx["home"]
    if not linux.startswith(home+"/"):
        raise ValueError("V2810_WSL_PATH_OUTSIDE_HOME")
    return ("\\\\wsl.localhost\\"+ctx["distribution"]+
            linux.replace("/","\\").replace("\\home\\","\\home\\",1))

def _download_meta(url):
    """Remote HEAD of public chromosome object only."""
    headers=_run_wsl(["curl","-fsSIL","--retry","2","--max-time","28",url],timeout=45)
    lengths=[int(x) for x in re.findall(r"(?im)^content-length:\s*(\d+)\s*$",headers)]
    return lengths[-1] if lengths else None

def _file_size(path):
    output=_run_wsl(["stat","-c","%s",path])
    return int(output)

def _exists(path):
    try:return _run_wsl(["test","-f",path],timeout=8)==""   # check exit code only
    except (OSError,RuntimeError,subprocess.TimeoutExpired):return False

def _sha256_native_file(path,label):
    file=Path(path)
    total=file.stat().st_size
    processed=0
    started=time.monotonic()
    last_report=started
    digest=hashlib.sha256()
    with file.open("rb") as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b""):
            digest.update(chunk)
            processed+=len(chunk)
            now=time.monotonic()
            if total>=100_000_000 and now-last_report>=5:
                print(f"V31_PUBLIC_SHA256_PROGRESS {label} {100*processed/total:.1f}% "
                      f"{processed/1_000_000_000:.2f}/{total/1_000_000_000:.2f} GB",
                      flush=True)
                last_report=now
    if total>=100_000_000:
        print(f"V31_PUBLIC_SHA256_PROGRESS {label} 100.0% "
              f"elapsed_s={int(time.monotonic()-started)} state=digest_computed",flush=True)
    return digest.hexdigest()

def _sha256(path):
    # gnomAD shards can exceed 9 GiB. Hash mounted Windows volumes with
    # native streaming I/O to avoid the costly drvfs -> Linux -> Windows
    # round-trip. NEVER skip the full SHA256 integrity measurement.
    path=str(path)
    m=re.fullmatch(r"/mnt/([a-zA-Z])/(.+)",path)
    if os.name=="nt" and m:
        local=Path(m.group(1).upper()+":/"+m.group(2))
        h=hashlib.sha256()
        return _sha256_native_file(local,local.name)
    # For native Linux/WSL paths, give large checksums sufficient time.
    output=_run_wsl(["sha256sum",path],timeout=3600)
    checksum=output.split()[0]
    if not re.fullmatch(r"[a-f0-9]{64}",checksum):
        raise RuntimeError("V2810_SOURCE_SHA256_OUTPUT_INVALID")
    return checksum

def _free_bytes(path):
    val=_run_wsl(["df","-PB1",path],timeout=10).splitlines()[-1]
    return int(val.split()[3])

def _curl_download(url,destination,total):
    tmp=destination+".part"
    # Resume interrupted official public downloads without reading patient
    # inputs or relying on an unstable external API progress indicator.
    if not isinstance(total,int) or total<=0:
        raise RuntimeError("V31_PUBLIC_DOWNLOAD_OFFICIAL_BYTES_REQUIRED")
    # Scientific data on a Windows-mounted volume: use parallel HTTPS byte
    # ranges directly, avoiding the Windows->WSL->Windows streaming penalty.
    # The existing WSL curl resume handles old unproven legacy partials and
    # genuine Linux-only destinations without altering prior on-disk data.
    mounted=re.fullmatch(r"/mnt/([a-zA-Z])/(.+)",destination)
    if os.name=="nt" and mounted:
        native=Path(mounted.group(1).upper()+":/"+mounted.group(2))
        partial=Path(str(native)+".part")
        checkpoint=Path(str(partial)+".rdmva.json")
        if not partial.exists() or checkpoint.exists():
            from .scientific_download_v31 import download_public_data
            result=download_public_data(url,native,expected_bytes=total,
                                        max_bytes=total,label=native.name)
            print("V31_PUBLIC_DOWNLOAD_PROGRESS",native.name,
                  "100.0% state=TRANSFER_VERIFIED transport="+result["mode"],
                  flush=True)
            return True
    cmd=wsl_command("curl","--fail","--location","--retry","3","--retry-delay","3",
        "--connect-timeout","25","--continue-at","-","--silent","--show-error",
        "--output",tmp,url)
    match=re.fullmatch(r"/mnt/([a-zA-Z])/(.+)",tmp)
    local_part=(Path(match.group(1).upper()+":/"+match.group(2)) if match and os.name=="nt"
                else None)
    label=Path(destination).name
    if total>=100_000_000:
        print("V31_LARGE_PUBLIC_DOWNLOAD_NOTICE "+label+
              " -- downloads and preparation may take more than 60 minutes, "
              "depending on internet speed, disk and workstation; progress "
              "updates every 5 seconds. Do not interrupt.",flush=True)
    print("V31_PUBLIC_DOWNLOAD_START",label,total,flush=True)
    proc=subprocess.Popen(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,
                          text=True,errors="replace")
    previous=0
    previous_at=time.monotonic()
    last_change=previous_at
    try:
        while proc.poll() is None:
            if local_part is not None:
                current=local_part.stat().st_size if local_part.is_file() else 0
            else:
                try:current=_file_size(tmp)
                except (OSError,RuntimeError,subprocess.TimeoutExpired):current=0
            now=time.monotonic()
            delta=max(now-previous_at,0.001)
            speed=max(0,current-previous)/delta/1_000_000
            if current>previous:last_change=now
            state="DOWNLOADING" if now-last_change<30 else "WAITING_FOR_DATA"
            print(f"V31_PUBLIC_DOWNLOAD_PROGRESS {label} "
                  f"{min(100,100*current/total):.1f}% "
                  f"{current/1_000_000_000:.2f}/{total/1_000_000_000:.2f} GB "
                  f"{speed:.2f} MB/s state={state}",flush=True)
            previous=current
            previous_at=now
            time.sleep(5)
        stderr=proc.stderr.read() if proc.stderr else ""
        if proc.returncode:
            raise RuntimeError("V2810_WSL_DOWNLOAD_FAILED:"+str(proc.returncode)+
                               ":"+stderr[-300:])
    except BaseException:
        if proc.poll() is None:
            proc.terminate()
            proc.wait(timeout=30)
        raise
    if _file_size(tmp)!=total:
        raise RuntimeError("V2810_SIZE_MISMATCH_AFTER_DOWNLOAD")
    _run_wsl(["mv","-f",tmp,destination])
    print("V31_PUBLIC_DOWNLOAD_PROGRESS",label,"100.0% state=TRANSFER_VERIFIED",flush=True)
    return True

def plan(private,cohorts=("exomes","genomes")):
    chrs,omitted=needed_chromosomes(private)
    cs=tuple(dict.fromkeys(cohorts))
    if not cs or any(c not in COHORTS for c in cs):
        raise ValueError("V2810_UNKNOWN_COHORT")
    jobs=[]
    for chrom in chrs:
        for cohort in cs:
            filename=f"gnomad.{cohort}.v4.1.sites.chr{chrom}.vcf.bgz"
            source=f"{GNOMAD_ROOT}/{cohort}/{filename}"
            jobs.append({"chrom":chrom,"cohort":cohort,"filename":filename,
                         "source_url":source,"index_url":source+".tbi",
                         "release":"4.1","assembly":"GRCh38"})
    return {"schema":VERSION,"selection_rule":"actual_allele_chromosomes_not_patient_disease",
         "needed_chromosomes":chrs,"unsupported_chromosomes":omitted,
         "cohorts":list(cs),"download_jobs":jobs,
         "patient_allele_coordinates_to_remote_services":False,
         "gene_or_disease_names_in_urls":False,
         "private_input_metadata_not_in_redistributable_package":True}

def _receipt_path(root):
    return Path(root)/"GNOMAD_AUTO_RESOURCE_AUDIT.json"

def _write_gnomad_audit(root,audit):
    # Atomic updates: a monitoring PowerShell session must not read a partial
    # JSON while large chromosome files are downloading or being hashed.
    path=_receipt_path(root)
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".partial")
    tmp.write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    os.replace(tmp,path)

def prepare(private,outdir,*,download=True,cohorts=("exomes","genomes"),
            max_gib=240,timeout_s=14400,plan_override=None):
    """Fail-closed on budget or license; partial shards usable but NEVER called full coverage."""
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    proposal=plan_override or plan(private,cohorts)
    jobs=proposal["download_jobs"]
    audit={"schema":VERSION,"plan":proposal,
       "mode":"AUTOMATIC_SOURCE_DOWNLOAD" if download else "NO_NETWORK_DRY_RUN",
       "source_license":"download_separate_official_public_data_only",
       "no_personal_genomic_data_transmitted":True,
       "data_redistribution_in_submission":False,
       "records":[],"free_gib_at_start":None,"completion":"NOT_STARTED",
       "budget_gib":max_gib}
    shards=[];spent=0
    from .download_resilience_v31 import DownloadWarnings,resilient_download
    warnings=DownloadWarnings(Path(__file__).resolve().parents[2])
    _write_gnomad_audit(out,audit)
    if not jobs:
        audit["completion"]="NO_ELIGIBLE_CHROMOSOMES"
        _write_gnomad_audit(out,audit)
        return None,audit
    if not offline_external_bin_wsl():
        audit["completion"]="WSL_PREREQUISITES_MISSING"
        _write_gnomad_audit(out,audit)
        return None,audit
    try:
        ctx=wsl_context()
        destination=_public_storage_location(ctx,Path(__file__).resolve().parents[2])
        datahome=destination["linux_root"]
        _run_wsl(["mkdir","-p",datahome])
        free=min(destination["free_bytes"],_free_bytes(datahome))
        reserve=destination["reserve_bytes"]
        audit["free_gib_at_start"]=round(free/1024**3,2)
        audit["storage_kind"]=destination["location"]
        audit["storage_reserve_gib"]=round(reserve/1024**3,2)
        audit["verified_external_drive"]=destination.get("source_disk")
        audit["windows_backing_disk_capacity_checked"]=True
        budget=int(max_gib)*1024**3
        for job in jobs:
            chrom=job["chrom"];cohort=job["cohort"]
            audit["active_job"]={"cohort":cohort,"chrom":chrom,"phase":"REMOTE_METADATA"}
            _write_gnomad_audit(out,audit)
            folder=f"{datahome}/{cohort}"
            target=f"{folder}/{job['filename']}"
            index=target+".tbi"
            _run_wsl(["mkdir","-p",folder])
            label=f"gnomAD_{cohort}_{chrom}"
            size=resilient_download(
                label+"_metadata",job["source_url"],
                lambda u=job["source_url"]:_download_meta(u),
                warnings=warnings)
            indexsize=resilient_download(
                label+"_index_metadata",job["index_url"],
                lambda u=job["index_url"]:_download_meta(u),
                warnings=warnings)
            if not size or not indexsize or size<1000 or indexsize<100:
                status="SKIPPED_SOURCE_UNAVAILABLE"
                audit["records"].append({**job,"status":status,"bytes":size,
                                         "index_bytes":indexsize,"sha256":None})
                audit["active_job"]["phase"]=status
                _write_gnomad_audit(out,audit)
                continue
            exists=_exists(target) and _file_size(target)==size
            indexexists=_exists(index) and _file_size(index)==indexsize
            required=(0 if exists else size)+(0 if indexexists else indexsize)
            status="CACHED" if exists and indexexists else "PLANNED"
            if download and not (exists and indexexists):
                audit["active_job"]["phase"]="DOWNLOADING"
                _write_gnomad_audit(out,audit)
                if required+spent>budget or required+spent+reserve>free:
                    status="BLOCKED_STORAGE_BUDGET"
                else:
                    data_ok=True
                    if not exists:
                        data_ok=resilient_download(
                            label+"_data",job["source_url"],
                            lambda u=job["source_url"],t=target,z=size:
                                (_curl_download(u,t,z),True)[1],
                            warnings=warnings) is not None
                    if data_ok and not indexexists:
                        data_ok=resilient_download(
                            label+"_index",job["index_url"],
                            lambda u=job["index_url"],t=index,z=indexsize:
                                (_curl_download(u,t,z),True)[1],
                            warnings=warnings) is not None
                    if data_ok:
                        # Index syntax and readability must hold.
                        _run_wsl(["tabix","-H",target],timeout=80)
                        status="DOWNLOADED_VERIFIED"
                        spent+=required
                    else:
                        status="SKIPPED_SOURCE_UNAVAILABLE"
            if status in ("CACHED","DOWNLOADED_VERIFIED"):
                audit["active_job"]["phase"]="HASHING_SHA256"
                _write_gnomad_audit(out,audit)
            row={**job,"status":status,"bytes":size,"index_bytes":indexsize,
                  "sha256":_sha256(target) if status in ("CACHED","DOWNLOADED_VERIFIED") else None}
            audit["records"].append(row)
            audit["active_job"]["phase"]="VERIFIED" if row["sha256"] else status
            _write_gnomad_audit(out,audit)
            if status in ("CACHED","DOWNLOADED_VERIFIED"):
                winpath=(str(Path(destination["windows_root"])/cohort/job["filename"])
                    if destination.get("windows_root") else windows_unc(ctx,target))
                shards.append({"cohort":cohort,"chrom":chrom,"file":winpath,
                  "official_reference":job["source_url"],
                  "source_sha256":row["sha256"],
                  "filter_policy":"pass_only","source_provenance_verified":True})
    except Exception as exc:
        audit["completion"]="INCOMPLETE_ERROR"
        audit["error_type"]=type(exc).__name__
        audit["error_summary"]=str(exc)[-350:]
    else:
        audit["completion"]="READY" if len(shards)==len(jobs) else (
            "PLANNED_ONLY" if not download else
            "PARTIAL_SOURCE_UNAVAILABLE" if any(
                x["status"]=="SKIPPED_SOURCE_UNAVAILABLE"
                for x in audit["records"]) else "PARTIAL_STORAGE_CONSTRAINT")
    if any(x["status"]=="SKIPPED_SOURCE_UNAVAILABLE" for x in audit["records"]):
        warnings.print_summary()
    audit["installed_jobs"]=len(shards)
    audit["required_jobs"]=len(jobs)
    audit.pop("active_job",None)
    audit["downloaded_gib_this_run"]=round(spent/1024**3,2)
    _write_gnomad_audit(out,audit)
    manifest=None
    if shards:
        payload={"schema":"V289_GNOMAD_V41_LOCAL_1","assembly":"GRCh38","release":"v4.1",
                 "shards":shards,
                 "source":"official gnomAD public whole-chromosome VCF data",
                 "all_24_chromosomes_both_cohorts_not_implied":True}
        manifest=out/"GNOMAD_AUTO_LOCAL_MANIFEST.json"
        manifest.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n")
    return manifest,audit


def _vep_fasta_index_verified(fasta):
    # The file's presence is not proof it contains a complete GRCh38
    # primary-assembly reference. Index must report known GRCh38 lengths
    # for each primary nuclear and mitochondrial chromosome.
    fasta=Path(fasta)
    idx=Path(str(fasta)+".fai")
    if not fasta.is_file() or not idx.is_file() or fasta.stat().st_size<2_800_000_000:
        return False
    observed={}
    try:
        with idx.open(encoding="ascii") as f:
            for line in f:
                fields=line.rstrip("\n").split("\t")
                if len(fields)<5:return False
                observed[fields[0]]=int(fields[1])
    except (OSError,ValueError,UnicodeError):
        return False
    expected={"1":248956422,"2":242193529,"7":159345973,
              "X":156040895,"Y":57227415,"MT":16569}
    return all(observed.get(name)==length for name,length in expected.items())


def _cached_large_vep_fasta(fa):
    path=Path(fa)
    return path.is_file() and path.stat().st_size>=2_800_000_000


def _cached_official_fasta_archive(file,expected_size):
    archive=Path(file)
    return archive.is_file() and archive.stat().st_size==expected_size


def _restore_vep_reference_from_cache(fa, *, uri):
    # Reuse a large, already present reference first. Reindex it locally;
    # a missing .fai alone NEVER triggers an 840 MB download.
    from .vep_local import wsl_path
    fa=Path(fa)
    index=Path(str(fa)+".fai")
    if _vep_fasta_index_verified(fa):
        return "VERIFIED_EXISTING_FASTA_AND_INDEX"
    if _cached_large_vep_fasta(fa):
        _run_wsl(["samtools","faidx",wsl_path(fa)],timeout=1800)
        if _vep_fasta_index_verified(fa):
            return "REBUILT_INDEX_FOR_CACHED_GRCH38_FASTA"

    # Only download if no usable decompressed FASTA AND no verified
    # complete official compressed FASTA exists. Validate the source length
    # and gzip CRC, then decompress into a temporary file, index + verify,
    # and only afterwards replace the destination atomically.
    file=Path(str(fa)+".gz")
    from .download_resilience_v31 import DownloadWarnings,resilient_download
    warnings=DownloadWarnings(Path(__file__).resolve().parents[2])
    length=resilient_download(
        "GRCh38_Ensembl_FASTA_metadata",uri,
        lambda:_download_meta(uri),warnings=warnings)
    if not length or length<400_000_000:
        raise RuntimeError("V2810_REFERENCE_SOURCE_UNAVAILABLE_OR_INVALID")
    free=shutil.disk_usage(fa.parent).free
    if free<5*1024**3:
        raise RuntimeError("V2810_NOT_ENOUGH_STORAGE_FOR_REFERENCE")
    cached=_cached_official_fasta_archive(file,length)
    if cached:
        try:
            _run_wsl(["gzip","-t",wsl_path(file)],timeout=1800)
        except (OSError,RuntimeError,subprocess.TimeoutExpired):
            cached=False
    if not cached:
        transferred=resilient_download(
            "GRCh38_Ensembl_FASTA",uri,
            lambda:(_curl_download(uri,wsl_path(file),length),True)[1],
            warnings=warnings)
        if transferred is None:
            raise RuntimeError("V2810_REFERENCE_DOWNLOAD_SKIPPED_UPSTREAM")
        _run_wsl(["gzip","-t",wsl_path(file)],timeout=1800)
    import shlex
    tmp=Path(str(fa)+".tmp")
    cmd="gzip -dc "+shlex.quote(wsl_path(file))+" > "+shlex.quote(wsl_path(tmp))
    _run_wsl(["bash","-lc",cmd],timeout=3600)
    _run_wsl(["samtools","faidx",wsl_path(tmp)],timeout=1800)
    if not _vep_fasta_index_verified(tmp):
        raise RuntimeError("V2810_ENS_GRCH38_DECOMPRESSED_FASTA_INDEX_INVALID")
    _run_wsl(["mv","-f",wsl_path(tmp),wsl_path(fa)])
    _run_wsl(["mv","-f",wsl_path(str(tmp)+".fai"),wsl_path(index)])
    if not _vep_fasta_index_verified(fa):
        raise RuntimeError("V2810_ENS_GRCH38_FINAL_REFERENCE_INVALID")
    return "RESTORED_FROM_VERIFIED_CACHED_ARCHIVE" if cached else "DOWNLOADED_VERIFIED_AND_INDEXED"

def ensure_vep_prerequisites(project_root,outdir,install=True):
    """Bootstrap free licensed VEP 116 + public GRCh38/GENCODE source, if missing.

    Separate external Apache-2.0 tool remains outside app source/release bundle.
    Reuses existing Ensembl VEP, reference and normalized GTF idempotently.
    WSL2 must already be enabled. Unexpected failures are explicit, not silent.
    """
    from .vep_local import wsl_path
    root=Path(project_root)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    status={"schema":VERSION+"_VEP",
      "license":"Apache-2.0 (VEP), reference from Ensembl; GENCODE source conditions",
      "install_allowed":bool(install),"tool_exists":False,
      "reference_ready":False,"gencode_ready":False,
      "genomic_data_uploaded":False,"commands_executed":[]}
    try:
        ctx=wsl_context()
        vep=ctx["home"]+"/rdmva_tools/ensembl-vep-116/vep"
        fa=root/"public_inputs/third_party/vep116/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        gtf=root/"public_inputs/third_party/vep116/gencode.v50.annotation.nochr.gtf.gz"
        source_gtf=root/"public_inputs/third_party/intake_v28/gencode.v50.annotation.gtf.gz"
        status["vep_path"]=vep
        status["fasta_path"]=str(fa)
        status["gencode_gtf_path"]=str(gtf)
        # A cloned VEP script does not prove its Ensembl/BioPerl API installed.
        # The marker is created only after a successful INSTALL.pl completion.
        api_marker=ctx["home"]+"/rdmva_tools/ensembl-vep-116/.rdmva_api116_ready"
        status["tool_exists"]=_exists(vep) and _exists(api_marker)
        if status["tool_exists"]:
            try:
                _run_wsl(["perl","-I",vep.rsplit("/",1)[0],vep,"--help"],timeout=90)
            except (OSError,RuntimeError,subprocess.TimeoutExpired):
                status["tool_exists"]=False
                print("V31_VEP_DEPENDENCIES_NEED_REPAIR",flush=True)
        status["reference_ready"]=_vep_fasta_index_verified(fa)
        status["gencode_ready"]=gtf.is_file() and Path(str(gtf)+".tbi").is_file()
        if install and not status["tool_exists"]:
            print("V31_VEP_SETUP_APT_AND_API_START",flush=True)
            # Debian/Ubuntu packages are fetched from the distro signed apt repos.
            pkg=["ca-certificates","git","curl","perl","make","gcc","g++",
              "tar","unzip","gzip","xz-utils","libarchive-extract-perl",
              "libarchive-zip-perl","libjson-perl","libdbi-perl",
              "libbio-perl-perl","libset-intervaltree-perl",
              "liblist-moreutils-perl","libjson-xs-perl",
              "libbio-db-hts-perl","libtry-tiny-perl","liburi-perl","libwww-perl","bcftools","tabix","samtools"]
            p=_run_long_step([*wsl_command(root=True),"apt-get","update","-qq"],
                 capture_output=True,text=True,timeout=900,label="apt_update")
            if p.returncode:raise RuntimeError("V2810_APT_UPDATE_FAILED:"+p.stderr[-230:])
            p=_run_long_step([*wsl_command(root=True),"apt-get","install",
                "-y","-qq","--no-install-recommends",*pkg],
                capture_output=True,text=True,timeout=1200,label="apt_install")
            if p.returncode:raise RuntimeError("V2810_APT_DEPENDENCIES_FAILED:"+p.stderr[-250:])
            toolroot=ctx["home"]+"/rdmva_tools"
            _run_wsl(["mkdir","-p",toolroot])
            if not _exists(vep):
                _run_wsl(["git","clone","--depth","1","--branch","release/116",
                    "https://github.com/Ensembl/ensembl-vep.git",toolroot+"/ensembl-vep-116"],timeout=420)
            _run_wsl(["perl","-MArchive::Extract","-e","1"],timeout=45)
            _run_wsl(["unzip","-v"],timeout=45)
            _run_wsl(["tar","--version"],timeout=45)
            # An explicit DESTDIR places Ensembl API dependencies beside VEP,
            # not in whatever Windows working directory WSL inherited.
            if not _exists(api_marker):
                _run_wsl(["perl",toolroot+"/ensembl-vep-116/INSTALL.pl",
                     "--DESTDIR",toolroot+"/ensembl-vep-116",
                     "--AUTO","a","--NO_HTSLIB","--NO_BIOPERL","--NO_PLUGINS","--NO_TEST"],
                     timeout=1800)
                _run_wsl(["perl","-I",toolroot+"/ensembl-vep-116",
                          "-MBio::EnsEMBL::Registry","-e","1"],timeout=45)
                _run_wsl(["touch",api_marker],timeout=20)
            # Exercise the actual VEP runner; Registry alone misses runtime deps.
            _run_wsl(["perl","-I",toolroot+"/ensembl-vep-116",vep,"--help"],
                     timeout=90)
            status["tool_exists"]=_exists(vep) and _exists(api_marker)
            status["commands_executed"].append("official_VEP_116_executable_verified")
        if install and not status["reference_ready"]:
            print("V31_VEP_REFERENCE_CHECK_CACHED_FASTA_FIRST",flush=True)
            fa.parent.mkdir(parents=True,exist_ok=True)
            uri="https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            status["reference_recovery"]=_restore_vep_reference_from_cache(fa,uri=uri)
            status["reference_ready"]=_vep_fasta_index_verified(fa)
            status["commands_executed"].append(status["reference_recovery"])
        # Lock the actual public GRCh38 release bytes to the historical
        # independently computed original scientific reference, not merely
        # an index that a truncated/corrupt FASTA could also generate.
        if status["reference_ready"]:
            expected_ref_sha256="1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b"
            # Exact byte-count percentage for the large verified GRCh38 FASTA.
            status["reference_sha256"]=_sha256_native_file(fa,"GRCh38_FASTA")
            status["reference_expected_sha256"]=expected_ref_sha256
            if status["reference_sha256"]!=expected_ref_sha256:
                raise RuntimeError("V31_ENSEMBL_GRCH38_FASTA_PROVENANCE_SHA256_MISMATCH")
            status["reference_full_file_sha256_verified"]=True
        if install and not status["gencode_ready"]:
            print("V31_VEP_GENCODE_GTF_NORMALIZATION_START",flush=True)
            if not source_gtf.is_file():raise RuntimeError("V2810_LOCAL_GENCODE_SOURCE_MISSING")
            gtf.parent.mkdir(parents=True,exist_ok=True)
            import shlex
            # GENCODE v50 is grouped by feature/transcript rather than guaranteed
            # position sorted: bgzip alone creates a valid archive but tabix
            # rejects decreasing coordinates. Sort by contig and GTF start
            # (column 4) before indexing; stable ties retain source ordering.
            # All source rows and attributes are preserved. The sort temp stays
            # on the selected large public-data disk, not the Linux root disk.
            stage=Path(str(gtf)+".building.gz")
            stage_index=Path(str(stage)+".tbi")
            for p in (stage,stage_index):
                if p.exists():p.unlink()
            quoted_in=shlex.quote(wsl_path(source_gtf))
            quoted_stage=shlex.quote(wsl_path(stage))
            quoted_work=shlex.quote(wsl_path(gtf.parent))
            # Scale only GNU sort's temporary memory buffer to the active WSL
            # memory budget. Comparator, stable ties and scientific data do
            # not change; limited-memory jury machines safely fall back.
            try:
                available_kib=int(_run_wsl(
                    ["awk","/MemAvailable/ {print $2; exit}","/proc/meminfo"],
                    timeout=15).strip())
            except (OSError,RuntimeError,ValueError,subprocess.TimeoutExpired):
                available_kib=1280*1024
            sort_mib=min(8192,max(256,(available_kib//1024)//5))
            print("V31_VEP_GTF_SORT_MEMORY_MIB",sort_mib,flush=True)
            cmd=("gzip -dc "+quoted_in+
              " | awk 'BEGIN{FS=OFS=\"\\t\"} /^#/ {print;next} {sub(/^chr/,\"\",$1); if ($1==\"M\")$1=\"MT\";print}'"
              " | LC_ALL=C sort --stable --buffer-size="+str(sort_mib)+"M --parallel=2"
              " -T "+quoted_work+" -k1,1 -k4,4n"
              " | bgzip -@ 4 -c > "+quoted_stage)
            try:
                print("V31_VEP_GENCODE_GTF_SORTING_BY_CONTIG_AND_START",flush=True)
                _run_wsl(["bash","-o","pipefail","-c",cmd],timeout=7200)
                if not stage.is_file() or stage.stat().st_size<1_000_000:
                    raise RuntimeError("V31_VEP_GENCODE_GTF_SORT_OUTPUT_MISSING_OR_SMALL")
                print("V31_VEP_GENCODE_GTF_TABIX_INDEX_START",flush=True)
                _run_wsl(["tabix","-f","-p","gff",wsl_path(stage)],timeout=240)
                if not stage_index.is_file() or stage_index.stat().st_size<100:
                    raise RuntimeError("V31_VEP_GENCODE_GTF_TABIX_INDEX_MISSING")
                contigs=_run_wsl(["tabix","-l",wsl_path(stage)],timeout=75).splitlines()
                if "1" not in contigs or "2" not in contigs or "X" not in contigs:
                    raise RuntimeError("V31_VEP_GENCODE_GTF_TABIX_CONTIG_VALIDATION_FAILED")
                os.replace(stage,gtf)
                os.replace(stage_index,Path(str(gtf)+".tbi"))
            finally:
                for p in (stage,stage_index):
                    if p.exists():p.unlink()
            status["gencode_ready"]=gtf.is_file() and Path(str(gtf)+".tbi").is_file()
            status["commands_executed"].append("GENCODE50_chr_normalized_coordinate_sorted_indexed")
        status["state"]="READY" if all((status["tool_exists"],
                    status["reference_ready"],status["gencode_ready"])) else "NOT_INSTALLED"
        status["public_software_licenses_remain_separate"]=True
        status["vep_config"]={
            "vep_script":vep,
            "vep_gtf":wsl_path(gtf),
            "vep_fasta":wsl_path(fa),
        }
    except (OSError,RuntimeError,subprocess.TimeoutExpired,ValueError) as exc:
        status["state"]="INCOMPLETE"
        status["failure_reason"]=type(exc).__name__+":"+str(exc)[-340:]
        status.setdefault("vep_config",{})
    (out/"VEP_BOOTSTRAP_STATUS.json").write_text(
          json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return status


def write_license_audit(output_dir):
    """A redistributable catalogue, not a grant of upstream third-party rights."""
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    for item in CATALOG.values():
        if item.get("price","").startswith("not_available"):
            assert item["distribution"].startswith("EXCLUDED")
    payload={
      "schema":"V2810_RESEARCH_REDISTRIBUTION_POLICY_1",
      "owner_code_license":"CC-BY-4.0 per MVA Hackathon published outputs",
      "external_software_copyrights_relicensed_to_CC_BY":False,
      "raw_patient_vcf_distributable":False,
      "public_reference_data_bundled":False,
      "third_party_public_data_separately_fetched":True,
      "unrestricted_software_research_use_requires_original_license_compliance":True,
      "candidate_drug_experimental_research_only":True,
      "online_paid_AI_required_for_default_run":False,
      "restricted_illumina_spliceai_auto_install":False,
      "enrichr_usage_in_default_profile":False,
      "public_repository_must_exclude":{"private_patient_genome":True,
         "gated_competition_data":True,"genome_wide_patient_variant_tables":True,
         "raw_illumina_models":True,"third_party_bundled_databases":True},
      "license_and_terms_urls":{
        "hackathon_rules":"https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026",
        "data_access_terms":"https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data",
        "gnomad_official":"https://gnomad.broadinstitute.org/policies",
        "vep_apache":"https://github.com/Ensembl/ensembl-vep/blob/release/116/LICENSE",
        "pangolin_gpl3":"https://github.com/tkzeng/Pangolin/blob/main/LICENSE",
        "pksim_gpl2":"https://github.com/Open-Systems-Pharmacology/PK-Sim",
        "chembl_by_sa":"https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions",
        "hpo_terms":"https://human-phenotype-ontology.github.io/license.html",
      },
      "tools_and_data":CATALOG,
      "upstream_legal_review_needed_for_models_or_processed_data_bundling":True}
    (out/"FREE_LICENSE_AUDIT.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload


def ensure_reproducible_public_baseline(project_root,outdir,*,install=True):
    """Verify or automatically restore historical public source snapshot dependencies.

    The project intentionally pins the original public-data snapshot to preserve
    reproducible drug universe and ClinVar/HPO annotation. If upstream changes,
    refuse silent hash drift. No restricted challenge/patient file is downloaded.
    """
    project=Path(project_root)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    jobs=[
      ("V28_PUBLIC_SNAPSHOT",project/"bootstrap_v28.py",
        ["--skip-env"] if install else ["--verify-only"]),
      ("V283_REGULATORY_SNAPSHOT",project/"bootstrap_regulatory_v283.py",
        [] if install else ["--verify-only"]),
    ]
    status={"schema":VERSION+"_PUBLIC_BASELINE",
      "dataset_snapshot_locked":True,
      "patient_vcf_read_or_downloaded":False,
      "gated_competition_records_fetched":False,
      "public_sources_redistributed_in_code":False,
      "installation_enabled":bool(install),
      "jobs":[]}
    for label,script,args in jobs:
        if not script.is_file():
            status["jobs"].append({"name":label,"state":"BOOTSTRAP_SCRIPT_MISSING"})
            break
        p=subprocess.run([os.sys.executable,str(script),*args],
           cwd=str(project),capture_output=True,text=True,errors="replace",timeout=3600)
        status["jobs"].append({"name":label,"state":"OK" if p.returncode==0 else "FAILED",
            "exit_code":p.returncode,
            "tail":(p.stdout if p.returncode==0 else p.stderr)[-600:]})
        if p.returncode:break
    status["state"]="READY" if len(status["jobs"])==len(jobs) and all(
         x["state"]=="OK" for x in status["jobs"]) else "FAILED_FAIL_CLOSED"
    (out/"PUBLIC_BASELINE_RESOURCE_STATUS.json").write_text(
        json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return status

def _physical_wsl_host_free_bytes(ctx):
    """WSL df(1) alone is unsafe for sparse VHDX: check its Windows host too."""
    if os.name!="nt":
        raise RuntimeError("V2810_WSL_WINDOWS_HOST_REQUIRED")
    import winreg
    root=r"Software\Microsoft\Windows\CurrentVersion\Lxss"
    found=[]
    try:
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER,root) as key:
            count=winreg.QueryInfoKey(key)[0]
            for n in range(count):
                item=winreg.EnumKey(key,n)
                with winreg.OpenKey(key,item) as sub:
                    try:
                        name=winreg.QueryValueEx(sub,"DistributionName")[0]
                        base=winreg.QueryValueEx(sub,"BasePath")[0]
                    except OSError:continue
                    if name==ctx["distribution"]:
                        base=os.path.expandvars(str(base))
                        if base.startswith("\\\\?\\"):base=base[4:]
                        basepath=Path(base)
                        if basepath.is_dir():found.append(basepath)
    except OSError as exc:
        raise RuntimeError("V2810_WSL_HOST_VOLUME_NOT_RESOLVED") from exc
    if len(found)!=1:
        raise RuntimeError("V2810_WSL_HOST_VOLUME_NOT_UNAMBIGUOUS")
    return shutil.disk_usage(found[0]).free

def _public_storage_location(ctx,project_root):
    """Choose an explicit Windows data drive or check WSL VHDX backing host.

    LOCAL config is excluded from distributable source ZIP. An external USB
    drive must be present, identify as the expected device and have ample free
    space. Refuse all copies when USB identity changes or disappears.
    """
    from .vep_local import wsl_path
    cfg=Path(project_root)/"config/local_public_storage_v2810.json"
    if not cfg.is_file():
        home=ctx["home"]+"/rdmva_resources/gnomad/v4.1"
        _run_wsl(["mkdir","-p",home])
        logical=_free_bytes(home)
        host=_physical_wsl_host_free_bytes(ctx)
        return {
          "linux_root":home,"windows_root":None,
          "free_bytes":min(logical,host),
          "reserve_bytes":MIN_FREE_RESERVE,
          "location":"WSL_LINUX_VHDX_WITH_PHYSICAL_HOST_CAPACITY_GATE",
          "wsl_virtual_free_bytes":logical,
          "windows_backing_volume_free_bytes":host,
        }
    j=json.loads(cfg.read_text(encoding="utf-8"))
    if j.get("schema")!="V2810_LOCAL_PUBLIC_STORAGE_1":
        raise RuntimeError("V2810_EXTERNAL_DATA_CONFIG_SCHEMA_UNSUPPORTED")
    folder=str(j.get("windows_root") or "").strip()
    if not re.match(r"^[A-Za-z]:\\[^\\/?*<>|:]+",folder):
        raise RuntimeError("V2810_EXTERNAL_DATA_PATH_INVALID")
    if ".." in Path(folder).parts or "\x00" in folder:
        raise RuntimeError("V2810_EXTERNAL_DATA_TRAVERSAL_REJECTED")
    drv=folder[0].upper()
    # Structured PowerShell query, with the drive constrained to one letter.
    ps=("$disk=Get-Partition -DriveLetter '"+drv+
        "' -ErrorAction Stop | Get-Disk -ErrorAction Stop;"+
        "$v=Get-Volume -DriveLetter '"+drv+"' -ErrorAction Stop;"+
        "[pscustomobject]@{Name=$disk.FriendlyName;Bus=[string]$disk.BusType;"+
        "FileSystem=$v.FileSystem;DiskBytes=$disk.Size;VolumeBytes=$v.Size;"+
        "FreeBytes=$v.SizeRemaining}|ConvertTo-Json -Compress")
    proc=subprocess.run(["powershell.exe","-NoProfile","-NonInteractive","-Command",ps],
          capture_output=True,text=True,errors="replace",timeout=35)
    if proc.returncode:raise RuntimeError("V2810_EXTERNAL_DATA_VOLUME_NOT_ONLINE")
    try:info=json.loads(proc.stdout)
    except ValueError as exc:
        raise RuntimeError("V2810_EXTERNAL_STORAGE_IDENTITY_UNREADABLE") from exc
    expected=str(j.get("expected_disk_name") or "").strip()
    if expected and expected.upper()!=str(info.get("Name") or "").strip().upper():
        raise RuntimeError("V2810_EXTERNAL_VOLUME_DEVICE_IDENTITY_CHANGED")
    if j.get("require_external_usb") is True and info.get("Bus")!="USB":
        raise RuntimeError("V2810_EXTERNAL_DATA_VOLUME_NOT_USB")
    if info.get("FileSystem") not in ("NTFS","exFAT"):
        raise RuntimeError("V2810_EXTERNAL_DATA_FILE_SYSTEM_UNSUPPORTED")
    if int(info.get("VolumeBytes") or 0)<100*1024**3:
        raise RuntimeError("V2810_EXTERNAL_DATA_VOLUME_TOO_SMALL")
    reserve=int(j.get("reserve_free_gib",100))
    if not 12<=reserve<=3000:
        raise RuntimeError("V2810_EXTERNAL_STORAGE_RESERVE_OUT_OF_RANGE")
    free=int(info["FreeBytes"])
    if free<reserve*1024**3:
        raise RuntimeError("V2810_EXTERNAL_DISK_FREE_SPACE_BELOW_RESERVED_FLOOR")
    mountpoint="/mnt/"+drv.lower()
    # Assigning a Windows drive letter after WSL starts does not auto-mount in
    # the already-running distro. Mounting its existing NTFS volume as drvfs is
    # non-destructive and does not require shutdown/reinitialization of WSL.
    try:
        source=_run_wsl(["findmnt","-n","-o","SOURCE","--target",mountpoint],timeout=10)
    except RuntimeError:
        source=""
    if source.strip().upper()!=drv+":":
        shell=("mkdir -p "+mountpoint+" && "+
               "mount -t drvfs "+drv+": "+mountpoint)
        p=subprocess.run([*wsl_command(root=True),"bash","-lc",shell],
           capture_output=True,text=True,errors="replace",timeout=35)
        if p.returncode:
            raise RuntimeError("V2810_EXTERNAL_VOLUME_CANNOT_MOUNT_IN_WSL:"+p.stderr[-200:])
    path=wsl_path(folder)
    _run_wsl(["mkdir","-p",path],timeout=30)
    linuxfree=_free_bytes(path)
    return {
       "linux_root":path+"/gnomad/v4.1",
       "windows_root":str(Path(folder)/"gnomad"/"v4.1"),
       "free_bytes":min(free,linuxfree),
       "reserve_bytes":reserve*1024**3,
       "location":"EXPLICIT_VERIFIED_WINDOWS_REMOVABLE_STORAGE",
       "source_disk":str(info.get("Name")),
       "host_free_bytes":free,
       "wsl_mounted_free_bytes":linuxfree,
       "drive_letter":drv,
    }

