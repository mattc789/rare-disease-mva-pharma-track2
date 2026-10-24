"""Private local Ensembl VEP 116 / GRCh38 integration.

No allele coordinates leave the workstation. VEP is run via WSL with --offline.
VEP GENCODE GTF+FASTA annotation is distinct from full 27.6 GB VEP cache.
"""
from __future__ import annotations
import hashlib,json,os,re,subprocess,tempfile
from .resource_bootstrap import wsl_command
from pathlib import Path,PurePosixPath
from collections import Counter

SCHEMA="V289_VEP_OFFLINE_1"
PROJECT_ROOT=Path(__file__).resolve().parents[2]
MAX_VARIANTS=3000
def wsl_path(path):
    s=str(path)
    m=re.match(r"^([A-Za-z]):[\\/](.*)$",s)
    if m:return "/mnt/"+m.group(1).lower()+"/"+m.group(2).replace("\\","/")
    if s.startswith("\\\\wsl.localhost\\") or s.startswith("\\\\wsl$\\"):
        parts=s.lstrip("\\").split("\\")
        if len(parts)>=4 and parts[1] and all(p not in ("", ".", "..") for p in parts[2:]):
            return "/"+"/".join(parts[2:])
    if s.startswith("/") and "\x00" not in s:return s
    raise ValueError("V2810_LOCAL_FILE_PATH_EXPECTED")

def _wsl_home():
    """Resolve the active WSL user's home without assuming a username."""
    # The isolated research distro lives on the selected research drive.
    # Its first start on a slow USB disk can exceed 20 seconds.
    p=subprocess.run([*wsl_command(),"bash","-lc",'printf %s "$HOME"'],
        capture_output=True,text=True,errors="replace",timeout=90)
    home=p.stdout.strip()
    if p.returncode or not re.fullmatch(r"/home/[A-Za-z0-9_.-]+",home):
        raise RuntimeError("V2810_WSL_HOME_UNAVAILABLE")
    return home

def default_paths():
    """Machine-independent VEP paths derived from this checkout + active WSL."""
    project_wsl=wsl_path(PROJECT_ROOT)
    home=_wsl_home()
    return {
      "vep":home+"/rdmva_tools/ensembl-vep-116/vep",
      "gtf":project_wsl+"/public_inputs/third_party/vep116/gencode.v50.annotation.nochr.gtf.gz",
      "fasta":project_wsl+"/public_inputs/third_party/vep116/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    }

def _project_local_from_wsl(path):
    """Map only this checkout's /mnt/<drive>/... path back to Windows."""
    root=wsl_path(PROJECT_ROOT).rstrip("/")
    p=str(path)
    if p==root:return PROJECT_ROOT
    if p.startswith(root+"/"):
        return PROJECT_ROOT/Path(p[len(root)+1:].replace("/",os.sep))
    return None

def _resolved(vep=None,gtf=None,fasta=None):
    d=default_paths()
    return vep or d["vep"],gtf or d["gtf"],fasta or d["fasta"]

def validate_input(v):
    chrom=str(v.get("chrom") or "").removeprefix("chr")
    if not chrom or not re.fullmatch(r"(?:[1-9]|1\d|2[0-2]|X|Y|MT)",chrom):return False
    if not isinstance(v.get("pos"),int) or v["pos"]<1:return False
    ref=str(v.get("ref") or "").upper()
    alt=str(v.get("alt") or "").upper()
    return bool(ref and alt and len(ref)<=500 and len(alt)<=500
       and re.fullmatch("[ACGTN]+",ref) and re.fullmatch("[ACGTN]+",alt))

def variant_id(v):
    from .expert_genomics import safe_hash
    return "SHA256:"+safe_hash({**v,"gene":str(v.get("gene") or "").upper()})

def generate_vcf(variants,dest):
    """Only selected candidate alleles; no personal identifiers or genotype in VCF."""
    unique={}
    for v in variants:
        if validate_input(v):
            key=(str(v["chrom"]).removeprefix("chr"),v["pos"],v["ref"].upper(),v["alt"].upper())
            unique.setdefault(key,v)
    if len(unique)>MAX_VARIANTS:raise ValueError("V289_SELECTED_VARIANT_CAP_EXCEEDED")
    target=Path(dest)
    target.parent.mkdir(parents=True,exist_ok=True)
    lines=["##fileformat=VCFv4.2","##reference=GRCh38",
     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for key in sorted(unique,key=lambda x:(x[0],x[1],x[2],x[3])):
        lines.append("%s\t%d\t.\t%s\t%s\t.\tPASS\t."%key)
    target.write_text("\n".join(lines)+"\n",encoding="ascii")
    return unique

def prerequisite_status(vep=None,gtf=None,fasta=None,cache_dir=None):
    vep,gtf,fasta=_resolved(vep,gtf,fasta)
    locations=[vep,gtf,fasta]
    if cache_dir:locations.append(cache_dir)
    status={}
    for p in locations:
        local=_project_local_from_wsl(p)
        if local is not None:
            status[p]={"exists":local.exists(),"bytes":local.stat().st_size if local.is_file() else None}
        else:
            status[p]={"exists":None,"validation":"wsl_executable_or_external_path_verified_by_launch"}
    return status

def _run_vep(input_file,output_file,*,vep=None,gtf=None,
             fasta=None,cache_dir=None,timeout_s=1200):
    vep,gtf,fasta=_resolved(vep,gtf,fasta)
    # The official INSTALL.pl --DESTDIR writes the API to <VEP ROOT>/Bio.
    # Perl -I must reference its parent, NOT the Bio subdirectory.
    # Always use offline mode for patient-derived alleles: never ask an
    # external Ensembl server to annotate patient genomic coordinates.
    args=[*wsl_command(),"perl","-I",str(PurePosixPath(vep).parent),vep,
        "--species","homo_sapiens","--assembly","GRCh38",
        "--format","vcf","--json","--force_overwrite","--no_stats",
        "--input_file",wsl_path(input_file),"--output_file",wsl_path(output_file),
        "--fasta",fasta,"--symbol"]
    if cache_dir:
        args+=["--offline","--cache","--dir_cache",cache_dir]
    else:
        # Ensembl VEP 116 requires a populated ~/.vep cache when --offline
        # is specified, even if --gtf is also supplied. Official Ensembl
        # documentation confirms --gtf + --fasta alone is a local-only
        # annotation mode, without --cache or --database. Do NOT add either
        # remote database mode or a patient-coordinate network lookup here.
        args+=["--gtf",gtf]
    done=subprocess.run(args,capture_output=True,text=True,errors="replace",timeout=timeout_s)
    if done.returncode:
        raise RuntimeError("V289_VEP_LOCAL_FAILED:"+str(done.returncode)+
            " STDOUT="+done.stdout[-500:]+" STDERR="+done.stderr[-650:])
    if not Path(output_file).is_file():
        raise RuntimeError("V289_VEP_NO_OUTPUT")
    return {"run_command_kind":"wsl_local_perl_offline",
       "exit_code":done.returncode,
       "stderr_tail":done.stderr[-350:]}

def parse_json_lines(path):
    out=[]
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if not line.strip():continue
            j=json.loads(line)
            text=str(j.get("input") or "")
            aa=re.split(r"\s+",text)
            # VEP VCF JSON input generally CHROM POS ID REF ALT ...
            ref,alt=(aa[3].upper(),aa[4].upper()) if len(aa)>=5 else ("","")
            # VEP normalizes internal start/end for indels; selected input VCF
            # uses the original POS/REF/ALT. Join on ORIGINAL input allele
            # representation, otherwise frameshift indels silently disappear.
            orig_chrom=aa[0] if aa else str(j.get("seq_region_name") or "")
            orig_pos=int(aa[1]) if len(aa)>1 and aa[1].isdigit() else j.get("start")
            key=(orig_chrom,orig_pos,ref,alt)
            consequences=[]
            for x in (j.get("transcript_consequences") or []):
                consequences.append({
                  "gene_symbol":x.get("gene_symbol"),
                  "gene_id":x.get("gene_id"),"transcript_id":x.get("transcript_id"),
                  "consequence_terms":x.get("consequence_terms") or [],
                  "impact":x.get("impact"),"biotype":x.get("biotype"),
                  "canonical":x.get("canonical"),
                  "protein_start":x.get("protein_start"),
                  "hgvsc":x.get("hgvsc"),"hgvsp":x.get("hgvsp"),
                })
            out.append({"chrom":key[0],"pos":key[1],"ref":key[2],"alt":key[3],
              "most_severe_consequence":j.get("most_severe_consequence"),
              "transcript_annotations":consequences,
              "annotation_source":"Ensembl VEP offline; GTF mode does not imply gnomAD/ClinVar cache coverage",
              "clinical_pathogenicity_classification":"not_inferred",
            })
    return out

def annotation_stage(private,outdir,config=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    c=config or {}
    src=[{**v,"gene":g} for g,vs in (private.get("genes") or {}).items() for v in vs]
    valid=[v for v in src if validate_input(v)]
    source_mode="full_vep_cache" if c.get("vep_cache") else "offline_gencode_gtf"
    manifest={"schema":SCHEMA,"mode":source_mode,
      "total_patient_case_variants":len(src),"suitable_variant_records":len(valid),
      "raw_patient_vcf_reopened":False,"external_http_requests_for_patient_variants":0,
      "clinically_qualified_interpretation":False,
      "genome_coordinates_in_public_summary":False,
      "vep_version_requested":"116","fasta":"GRCh38 primary assembly",
      "annotation_source":"GENCODE v50" if not c.get("vep_cache") else "VEP cache manifest supplied",
      "tool_state":"not_checked"}
    if not valid:
        manifest["tool_state"]="NO_SUITABLE_VARIANTS"
        (out/"VEP_STATUS.json").write_text(json.dumps(manifest,indent=2)+"\n")
        return manifest,[]
    input_file=out/"vep_selected_PRIVATE.vcf"
    output_file=out/"vep_annotated_PRIVATE.json"
    generate_vcf(valid,input_file)
    try:
        provenance=_run_vep(input_file,output_file,vep=c.get("vep_script"),
                 gtf=c.get("vep_gtf"),
                 fasta=c.get("vep_fasta"),
                 cache_dir=c.get("vep_cache"),timeout_s=int(c.get("timeout_s",1200)))
        annotated=parse_json_lines(output_file)
        manifest["tool_state"]="COMPLETED"
        manifest["annotation_records_returned"]=len(annotated)
        manifest["tool_provenance"]=provenance
    except (RuntimeError,FileNotFoundError,subprocess.TimeoutExpired,ValueError,OSError) as exc:
        annotated=[];manifest["tool_state"]="UNAVAILABLE_OR_FAILED"
        manifest["warning"]=type(exc).__name__+":"+str(exc)[-350:]
    # Private VEP output stays local; public report gets aggregate and source flags only.
    (out/"VEP_STATUS.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return manifest,annotated
