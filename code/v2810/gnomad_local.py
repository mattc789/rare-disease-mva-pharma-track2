"""gnomAD v4.1 local-only indexed exome+genome allele-frequency lookup.

Never send variant coordinates to the network. No local VCF means UNKNOWN,
not AF=0. Completeness is true only if all 24 nuclear chromosomes are indexed
for BOTH exomes and genomes in the verified local manifest.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import json,re,subprocess
from .resource_bootstrap import wsl_command

SCHEMA="V289_GNOMAD_V41_LOCAL_1"
RELEASE="v4.1"
CHROMS={str(i) for i in range(1,23)}|{"X","Y"}
COHORTS={"exomes","genomes"}

def _wsl_file(path):
    from .vep_local import wsl_path
    return wsl_path(path)

def validate_manifest(path):
    if path is None or not Path(path).is_file():
        return {},{"state":"NOT_INSTALLED","v4_1_full_coverage":False,
                   "reason":"Complete GRCh38 gnomAD v4.1 exome/genome VCF shards not configured"}
    obj=json.loads(Path(path).read_text(encoding="utf-8"))
    if obj.get("schema")!=SCHEMA or obj.get("assembly")!="GRCh38" or obj.get("release")!=RELEASE:
        raise ValueError("V289_GNOMAD_SCHEMA_ASSEMBLY_OR_RELEASE_MISMATCH")
    manifest={}
    bad=[]
    for item in obj.get("shards",[]):
        cohort=str(item.get("cohort") or "")
        chrom=str(item.get("chrom") or "").removeprefix("chr").upper()
        file=str(item.get("file") or "")
        reference=str(item.get("official_reference") or "")
        if cohort not in COHORTS or chrom not in CHROMS or not file or not reference:
            bad.append("invalid_cohort_chrom_file_or_reference");continue
        if (cohort,chrom) in manifest:
            raise ValueError("V289_GNOMAD_DUPLICATE_SHARD")
        # The local query engine uses tabix on BGZF-compressed VCF, not BCF.
        # BCF inputs would require a separate bcftools query implementation.
        # Official gnomAD v4.1 shards are named *.vcf.bgz; tabix supports BGZF
        # independently of the .gz/.bgz extension. Accept both, not plain .vcf.
        if not file.endswith((".vcf.gz",".vcf.bgz")):
            bad.append(cohort+":"+chrom+":unsupported_format_expected_indexed_bgzf_vcf");continue
        ix=file+".tbi" if Path(file+".tbi").exists() else file+".csi"
        if not Path(file).is_file() or not Path(ix).is_file():
            bad.append(cohort+":"+chrom+":missing_index_or_data");continue
        manifest[(cohort,chrom)]={"file":file,"reference":reference,
            "filter_policy":item.get("filter_policy","pass_only"),
            "assembly":"GRCh38","release":RELEASE}
    expected={(c,ch) for c in COHORTS for ch in CHROMS}
    present=set(manifest)
    summary={"state":"LOCAL_COMPLETE" if present==expected else "LOCAL_PARTIAL",
      "release":RELEASE,"assembly":"GRCh38","shards_found":len(present),
      "shards_expected":len(expected),"v4_1_full_coverage":present==expected,
      "cohort_chromosome_completeness":{c:sorted(ch for (co,ch) in present if co==c)
        for c in sorted(COHORTS)},
      "missing_shard_count":len(expected-present),
      "invalid_or_missing_shards":bad[:30],
      "no_remote_variant_lookups":True}
    return manifest,summary

def parse_info(s):
    data={}
    for part in str(s or "").split(";"):
        if "=" in part:
            a,b=part.split("=",1);data[a]=b
    return data

def _number_for_allele(v,index):
    values=str(v or "").split(",")
    if index>=len(values) or values[index] in (".",""):
        return None
    try:return float(values[index])
    except (TypeError,ValueError):return None

def local_tabix(path,chrom,pos):
    # GRCh38 public datasets can use either chr21 or 21 sequence naming.
    # Looking up both styles never exposes private positions to the network:
    # tabix operates on the previously downloaded local whole-chromosome VCF.
    issues=[]
    for contig in (str(chrom),"chr"+str(chrom)):
        cmd=[*wsl_command(),"tabix",_wsl_file(path),f"{contig}:{pos}-{pos}"]
        p=subprocess.run(cmd,timeout=40,capture_output=True,text=True,errors="replace")
        if p.returncode:
            issues.append(p.stderr[-160:])
            continue
        if p.stdout.strip():return p.stdout.splitlines()
    if len(issues)==2:
        raise RuntimeError("V2810_GNOMAD_LOCAL_INDEX_QUERY_FAILED:"+issues[-1])
    return []

def lookup_allele(variant,shards,*,query=local_tabix):
    from .expert_genomics import safe_hash
    chrom=str(variant.get("chrom") or "").removeprefix("chr").upper()
    pos=variant.get("pos")
    ref=str(variant.get("ref") or "").upper()
    alt=str(variant.get("alt") or "").upper()
    if chrom not in CHROMS or not isinstance(pos,int) or not ref or not alt:
        return {"state":"VARIANT_KEY_UNSUPPORTED","measurements":[]}
    rows=[]
    for cohort in sorted(COHORTS):
        shard=shards.get((cohort,chrom))
        if not shard:
            rows.append({"cohort":cohort,"status":"DATASET_OR_CHROMOSOME_MISSING",
                         "allele_frequency":None})
            continue
        try:lines=query(shard["file"],chrom,pos)
        except (RuntimeError,OSError,subprocess.TimeoutExpired) as exc:
            rows.append({"cohort":cohort,"status":"INDEX_QUERY_FAILED",
                         "allele_frequency":None,"error":type(exc).__name__})
            continue
        found=False
        for line in lines:
            vals=line.split("\t")
            if len(vals)<8:continue
            if str(vals[0]).removeprefix("chr").upper()!=chrom:continue
            try:observed_pos=int(vals[1])
            except ValueError:continue
            if observed_pos!=pos or vals[3].upper()!=ref:continue
            alleles=vals[4].upper().split(",")
            if alt not in alleles:continue
            index=alleles.index(alt)
            info=parse_info(vals[7])
            af=_number_for_allele(info.get("AF"),index)
            ac=_number_for_allele(info.get("AC"),index)
            an=_number_for_allele(info.get("AN"),0)
            # Recompute AF if AF absent, provided AC/AN are trustworthy.
            if af is None and ac is not None and an is not None and an>0:
                af=ac/an
            flags=[]
            if vals[6] not in ("PASS","."):flags.append("source_filter_not_pass")
            if an is None or an<=0:flags.append("allele_number_missing")
            if af is None or not 0<=af<=1:flags.append("AF_unavailable_or_invalid")
            if flags:af=None
            rows.append({"cohort":cohort,"status":"LOCALLY_MATCHED" if af is not None else "SOURCE_RECORD_NONQUALIFYING",
              "allele_frequency":af,"allele_count":ac,"allele_number":an,
              "source_reference":shard["reference"],
              "source_filter":vals[6],"record_flags":flags,
              "allele_counts_are_specific_to_this_alt":True,
              "population_subcohort":"global_aggregate_only",
              "individual_population_frequencies_not_inferred":True})
            found=True
        if not found:rows.append({"cohort":cohort,"status":"NOT_OBSERVED_IN_THIS_INDEXED_LOCAL_SHARD",
           "allele_frequency":None,
           "absence_is_not_proof_of_zero_frequency":True})
    valid=[r["allele_frequency"] for r in rows if r.get("status")=="LOCALLY_MATCHED"]
    return {"variant_ref":"SHA256:"+safe_hash(variant),
       "state":"INDEXED_POPULATION_EVIDENCE" if valid else "FREQUENCY_UNDETERMINED",
       "measurements":rows,
       "max_valid_global_af":max(valid) if valid else None,
       "population_frequency_not_disease_causality":True}

def run(private,outdir,manifest_path=None,query=local_tabix):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    shards,status=validate_manifest(manifest_path)
    rows=[]
    for gene,vs in sorted((private.get("genes") or {}).items()):
        for v in vs:
            r=lookup_allele({**v,"gene":gene},shards,query=query)
            r["gene"]=gene
            rows.append(r)
    result={**status,"schema":SCHEMA,
      "variant_records_checked":len(rows),
      "variants_with_at_least_one_qualified_population_frequency":sum(
       r.get("state")=="INDEXED_POPULATION_EVIDENCE" for r in rows),
      "variants_without_population_frequency":sum(
       r.get("state")!="INDEXED_POPULATION_EVIDENCE" for r in rows),
      "client_genetic_variants_sent_to_remote_apis":False,
      "absence_of_record_interpreted_as_af_zero":False,
      "source_coverage_limitation":"24 nuclear chromosomes per cohort; mitochondrial genome not supplied by this gnomAD nuclear dataset",
      "clinical_pathogenicity_from_af_alone":False}
    (out/"GNOMAD_V41_STATUS.json").write_text(json.dumps(result,indent=2)+"\n")
    # Private local allele-hash evidence (does not include raw genomic variants).
    with (out/"GNOMAD_V41_PRIVATE_LOOKUPS.jsonl").open("w",encoding="utf-8") as f:
        for r in rows:f.write(json.dumps(r,sort_keys=True)+"\n")
    return result,rows
