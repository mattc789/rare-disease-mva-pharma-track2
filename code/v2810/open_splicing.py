"""GPLv3 Pangolin independent executable bridge (free splice prediction).

NO Illumina source/models are used. Pangolin itself remains under GPLv3.
GPL source and model weights are downloaded as a separate independent tool.
Predictions never substitute for patient RNA splice assay confirmation.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,json,subprocess,re
from .resource_bootstrap import wsl_context,_run_wsl,_exists,wsl_command
from .vep_local import generate_vcf,wsl_path

SCHEMA="V2810_FREE_PANGOLIN_GPL3_1"
UPSTREAM="https://github.com/tkzeng/Pangolin.git"

def _license_check(repo):
    s=_run_wsl(["head","-n","15",repo+"/LICENSE"])
    return ("GNU GENERAL PUBLIC LICENSE" in s and "Version 3" in s)

def ensure_source(output_dir,*,install=False):
    """Use pinned provenance; downloading a GPL3 source is not bundling it."""
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    status={"schema":SCHEMA,
      "license":"GPL-3.0 independent third-party program",
      "copyright_reassigned_to_CC_BY":False,
      "source_url":UPSTREAM,
      "public_source":True,"source_installed":False,
      "trained_models_license_needs_source_file_notices":True,
      "private_genotypes_submitted_to_external_hosts":False,
      "read_only_when_not_installed":True}
    try:
        ctx=wsl_context()
        repo=ctx["home"]+"/rdmva_tools/Pangolin"
        runtime=ctx["home"]+"/rdmva_tools/pangolin_py312"
        python=runtime+"/bin/python"
        cli=runtime+"/bin/pangolin"
        status["repo_wsl_path"]=repo
        status["source_installed"]=_exists(repo+"/LICENSE")
        if install and not status["source_installed"]:
            _run_wsl(["mkdir","-p",ctx["home"]+"/rdmva_tools"])
            _run_wsl(["git","clone","--depth","1",UPSTREAM,repo],timeout=900)
            status["source_installed"]=_exists(repo+"/LICENSE")
        if status["source_installed"]:
            if not _license_check(repo):
                raise RuntimeError("PANGOLIN_SOURCE_NOT_GPL3")
            status["git_commit"]=_run_wsl(["git","-C",repo,"rev-parse","HEAD"])
            model_check=_run_wsl(["bash","-lc",
               "find '"+repo+"/pangolin/models' -maxdepth 1 -type f | wc -l"])
            status["trained_model_files"]=int(model_check)
            if status["trained_model_files"]<12:
                raise RuntimeError("PANGOLIN_MODEL_FILES_INCOMPLETE")
        status["runtime_python_312_present"]=_exists(python)
        status["runtime_cli_present"]=_exists(cli)
        if install and status["source_installed"] and not status["runtime_cli_present"]:
            installer=Path(__file__).resolve().parent/"externals/install_pangolin_gpl_wsl.sh"
            if not installer.is_file():
                raise RuntimeError("V2810_FREE_PANGOLIN_INSTALLER_NOT_IN_SOURCE_BUNDLE")
            p=subprocess.run([*wsl_command(),"bash",wsl_path(installer)],
               capture_output=True,text=True,errors="replace",timeout=2400)
            if p.returncode:
                raise RuntimeError("V2810_FREE_PANGOLIN_RUNTIME_SETUP_FAILED:"+p.stderr[-260:])
            status["runtime_python_312_present"]=_exists(python)
            status["runtime_cli_present"]=_exists(cli)
        if status["source_installed"] and status["runtime_cli_present"]:
            # Independently installed Pangolin is GPL-3.0. Current PyVCF3
            # requires a seventh INFO type_code argument; use a pinned,
            # idempotent compatibility patch inside the GPL source tree.
            patch_file=Path(__file__).resolve().parent/"externals/patch_pangolin_pyvcf3.py"
            if not patch_file.is_file():
                raise RuntimeError("V2810_PANGOLIN_GPL_PATCH_NOT_INCLUDED")
            shim=subprocess.run([*wsl_command(),"python3",
                wsl_path(patch_file),repo],capture_output=True,text=True,
                timeout=30,errors="replace")
            if shim.returncode:
                raise RuntimeError("V2810_PANGOLIN_PYVCF_COMPAT_FAILED:"+shim.stderr[-180:])
            status["pyvcf3_compatibility"]=shim.stdout.strip()
        status["state"]=(
           "SOURCE_AND_RUNTIME_READY" if status["runtime_cli_present"] and status["source_installed"]
            else "GPL3_SOURCE_VERIFIED_RUNTIME_NOT_READY" if status["source_installed"]
            else "NOT_INSTALLED")
        status["runtime_python"]=python
        status["runtime_cli"]=cli
        status["gene_or_allele_data_transmitted_to_upstream_git"]=False
    except (RuntimeError,OSError,subprocess.TimeoutExpired,ValueError) as exc:
        status["state"]="UNAVAILABLE_OR_UNVERIFIED"
        status["reason"]=type(exc).__name__+":"+str(exc)[-240:]
    (out/"PANGOLIN_GPL3_SOURCE_STATUS.json").write_text(
       json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return status


def _write_selected_gene_gtf(source_gtf,genes,destination):
    """Minimize Pangolin annotation DB to each case's expressed gene models.

    The input is the independently downloaded, public GENCODE 50 GTF.
    Only transcript records tagged with the case-selected gene symbols are
    retained; no genotype coordinates are transmitted remotely.
    """
    import gzip,re
    wanted={str(g).upper() for g in genes if g}
    if not wanted:raise ValueError("PANGOLIN_NO_GENES_TO_ANNOTATE")
    pat=re.compile(r'(?:^|;\s*)gene_name\s+"([^"]+)"')
    selected=0
    dst=Path(destination)
    with gzip.open(source_gtf,"rt",encoding="utf-8",errors="replace") as inp, \
         gzip.open(dst,"wt",encoding="utf-8") as output:
        for line in inp:
            if line.startswith("#"):
                output.write(line)
                continue
            m=pat.search(line)
            if m and m.group(1).upper() in wanted:
                output.write(line)
                selected+=1
    if selected<3:
        dst.unlink(missing_ok=True)
        raise RuntimeError("PANGOLIN_SELECTED_GTF_TRANSCRIPTS_NOT_FOUND")
    return {"gene_count":len(wanted),"public_annotation_records_selected":selected,
      "private_variant_coordinates_in_gtf":False}


def parse_private_pangolin_scores(private,model_vcf,outdir):
    """Extract local Pangolin VCF INFO into hashed, variant-bound hypotheses.

    SCORE means an in-silico RNA splicing signal, not an observed phenotype.
    Genomic coordinates and alleles remain inside private VCF; they are not
    written to any source-code or public report.
    """
    from .expert_genomics import safe_hash
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    bykey={}
    for gene,vs in (private.get("genes") or {}).items():
        for v in vs:
            key=(str(v.get("chrom") or "").removeprefix("chr").upper(),
                int(v.get("pos") or 0),str(v.get("ref") or "").upper(),
                str(v.get("alt") or "").upper())
            bykey.setdefault(key,[]).append({"gene":gene,"variant_ref":"SHA256:"+safe_hash({**v,"gene":gene})})
    measured={}
    if Path(model_vcf).is_file():
        with Path(model_vcf).open(encoding="utf-8",errors="replace") as f:
            for line in f:
                if not line or line.startswith("#"):continue
                c=line.rstrip("\n").split("\t")
                if len(c)<8:continue
                try:pos=int(c[1])
                except ValueError:continue
                chrom=c[0].removeprefix("chr").upper()
                k=(chrom,pos,c[3].upper(),c[4].upper())
                rawinfo=c[7]
                pangolin=next((x.split("=",1)[1] for x in rawinfo.split(";")
                               if x.startswith("Pangolin=")),None)
                scores=[]
                if pangolin and pangolin!=".":
                    for v in re.findall(r":([-+]?(?:\d+(?:\.\d+)?|\.\d+))(?=$|[|,])",
                                        pangolin):
                        try:s=float(v)
                        except ValueError:continue
                        if -1.00001<=s<=1.00001:scores.append(s)
                for entry in bykey.get(k,[]):
                    measured[entry["variant_ref"]]={
                      "variant_ref":entry["variant_ref"],"gene":entry["gene"],
                      "state":"COMPUTATIONAL_SPLICE_SCORE" if scores else "NO_APPLICABLE_SPLICE_SCORE",
                      "number_of_scored_positions":len(scores),
                      "max_absolute_splice_score_change":round(max(map(abs,scores)),6) if scores else None,
                      "signed_max_position_change":round(max(scores,key=abs),6) if scores else None,
                      "source_model":"Pangolin GPL-3.0 original pretrained model weights",
                      "RNA_splicing_validated":False,
                      "allele_specific_functional_rescue_validated":False,
                      "clinical_pathogenicity_probability":None,
                    }
    records=[]
    for entries in bykey.values():
        for entry in entries:
            records.append(measured.get(entry["variant_ref"]) or {
                "variant_ref":entry["variant_ref"],"gene":entry["gene"],
                "state":"NO_APPLICABLE_SPLICE_SCORE",
                "number_of_scored_positions":0,
                "max_absolute_splice_score_change":None,
                "source_model":"Pangolin GPL-3.0 external executable",
                "RNA_splicing_validated":False,
                "allele_specific_functional_rescue_validated":False})
    target=out/"PANGOLIN_PRIVATE_HASHED_SCORES.jsonl"
    with target.open("w",encoding="utf-8") as f:
        for r in records:f.write(json.dumps(r,ensure_ascii=False,sort_keys=True)+"\n")
    return records

def score_case(private,outdir,*,install=False,run_if_ready=True):
    """Run on private selected variants, never via public API.

    Program returns SOURCE/DEPENDENCIES_MISSING if no runnable Pangolin executable.
    ML statistical predictions are marked *in silico*, not RNA functional rescue.
    """
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    status=ensure_source(out,install=install)
    status["case_variants_with_splice_prediction"]=0
    if not run_if_ready or status.get("state")!="SOURCE_AND_RUNTIME_READY":
        status["prediction_state"]="SKIPPED_RUNTIME_UNAVAILABLE" if run_if_ready else "SKIPPED_REQUESTED"
        (out/"PANGOLIN_PREDICTION_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        return status
    root=Path(__file__).resolve().parents[2]
    fasta=root/"public_inputs/third_party/vep116/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gtf=root/"public_inputs/third_party/vep116/gencode.v50.annotation.nochr.gtf.gz"
    if not fasta.is_file() or not gtf.is_file():
        status["prediction_state"]="GRCH38_REFERENCE_NOT_READY"
        (out/"PANGOLIN_PREDICTION_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
        return status
    ctx=wsl_context()
    genes=sorted(str(g).upper() for g in (private.get("genes") or {}))
    key=hashlib.sha256(("GENCODE50|"+",".join(genes)).encode("utf-8")).hexdigest()[:16]
    destroot=ctx["home"]+"/rdmva_resources/pangolin_annot/gene_groups/"+key
    modeldb=destroot+"/selected_gene_models.db"
    try:
        if not _exists(modeldb):
            # An all-genome GENCODE gffutils build can take many minutes and
            # hundreds of MB. Build a precise per-gene selection automatically.
            public_subset=out/"selected_gene_models_PRIVATE.gtf.gz"
            sampling=_write_selected_gene_gtf(gtf,genes,public_subset)
            status["annotation_subset"]=sampling
            _run_wsl(["mkdir","-p",destroot])
            db_gtf=destroot+"/selected_gene_models.gtf.gz"
            _run_wsl(["cp","-f",wsl_path(public_subset),db_gtf],timeout=90)
            _run_wsl([status["runtime_python"],status["repo_wsl_path"]+"/scripts/create_db.py",
                      db_gtf,"--filter","None"],timeout=480)
            if not _exists(modeldb):
                raise RuntimeError("PANGOLIN_ANNOTATION_DB_NOT_CREATED")
        status["gene_group_db_cached"]=True
        variants=[v for gene,vs in (private.get("genes") or {}).items() for v in vs]
        selected_vcf=out/"pangolin_selected_PRIVATE.vcf"
        selected=generate_vcf(variants,selected_vcf)
        if not selected:
            status["prediction_state"]="NO_SUPPORTED_SELECTED_ALLELES"
        else:
            prefix=out/"pangolin_model_PRIVATE"
            p=subprocess.run([*wsl_command(),status["runtime_cli"],
               wsl_path(selected_vcf),wsl_path(fasta),modeldb,wsl_path(prefix)],
               timeout=3600,capture_output=True,text=True,errors="replace")
            status["executable_exit_code"]=p.returncode
            if p.returncode:
                status["prediction_state"]="PREDICTION_FAILED"
                status["failure_reason"]=p.stderr[-500:]
            else:
                status["prediction_state"]="COMPUTATIONAL_SPLICING_SCORES_WRITTEN"
                status["source_output_private_prefix"]=str(prefix)
                normalized=parse_private_pangolin_scores(private,str(prefix)+".vcf",out)
                status["case_variants_evaluated"]=len(normalized)
                status["case_variants_with_splice_prediction"]=sum(
                   x.get("state")=="COMPUTATIONAL_SPLICE_SCORE" for x in normalized)
                status["missing_splice_score_is_not_normal_splicing"]=True
                status["hashed_scores_file"]="PANGOLIN_PRIVATE_HASHED_SCORES.jsonl"
    except (RuntimeError,OSError,ValueError,subprocess.TimeoutExpired) as exc:
        status["prediction_state"]="PREDICTION_FAILED_OR_SKIPPED"
        status["failure_reason"]=type(exc).__name__+":"+str(exc)[-350:]
    status["rna_functional_rescue_confirmed"]=False
    status["trained_model_is_not_illumina"]=True
    status["external_genome_transmission"]=False
    (out/"PANGOLIN_PREDICTION_STATUS.json").write_text(
        json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return status
