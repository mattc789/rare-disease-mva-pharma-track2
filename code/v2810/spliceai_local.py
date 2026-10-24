"""License-gated LOCAL SpliceAI annotations; no unauthorized model weights.

Illumina SpliceAI trained models and precomputed predictions may require a
commercial license. This project never silently downloads or runs them.
Selected VCF alleles are compared locally to supplied indexed score VCFs.
"""
from __future__ import annotations
import json,re,subprocess
from .resource_bootstrap import wsl_command
from pathlib import Path
from collections import Counter

SCHEMA="V289_SPLICEAI_LOCAL_LICENSED_1"
SOURCE_REF="https://github.com/Illumina/SpliceAI"
FIELDS=("DS_AG","DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")
CATEGORIES=((0.8,"high_precision_prediction"),(0.5,"recommended_threshold"),(0.2,"high_recall_threshold"))
def validate_manifest(path):
    if not path or not Path(path).is_file():
        return None,{"state":"LICENSE_AND_LOCAL_PREDICTIONS_NOT_CONFIGURED",
          "commercial_use_automatically_permitted":False,
          "reason":"SpliceAI trained models and precomputed scores require a suitable license for commercial use"}
    data=json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema")!=SCHEMA or data.get("assembly")!="GRCh38":
        raise ValueError("V289_SPLICEAI_MANIFEST_SCHEMA_OR_ASSEMBLY_MISMATCH")
    if (data.get("licensee")!="ExaMedix Technologies"
        or data.get("commercial_use_permitted") is not True
        or data.get("license_scope") not in ("commercial", "commercial_research")
        or not str(data.get("license_reference") or "").strip()):
        return None,{"state":"LICENSE_NOT_ATTESTED",
             "commercial_use_automatically_permitted":False}
    shards=[]
    for c in data.get("score_files",[]):
        path=str(c.get("file") or "")
        if not path.endswith(".vcf.gz") or not Path(path).is_file() or not (
           Path(path+".tbi").is_file() or Path(path+".csi").is_file()):
            continue
        if c.get("assembly")!="GRCh38" or c.get("source") not in ("Illumina_SpliceAI","Ensembl_SpliceAI_scores"):
            continue
        shards.append(c)
    state="ACTIVE_WITH_LICENSED_LOCAL_DATA" if shards else "LICENSE_ATTESTED_BUT_DATA_NOT_INSTALLED"
    return shards,{"state":state,"commercial_use_automatically_permitted":bool(shards),
       "source":SOURCE_REF,"score_files_loaded":len(shards),
       "deidentified_alleles_transmitted_to_external_services":False,
       "prediction_is_clinical_diagnosis":False}

def parse_score_text(encoded,alt,gene=None):
    """SpliceAI=ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL."""
    vals=[]
    for entry in str(encoded or "").split(","):
        fields=entry.split("|")
        if len(fields)!=10:continue
        if fields[0].upper()!=str(alt).upper():continue
        if gene and fields[1].upper()!=str(gene).upper():continue
        scores={}
        for k,v in zip(FIELDS,fields[2:]):
            try:scores[k]=float(v)
            except (TypeError,ValueError):scores[k]=None
        if any(scores[k] is None or not 0<=scores[k]<=1 for k in FIELDS[:4]):continue
        maximum=max(scores[k] for k in FIELDS[:4])
        category="below_exploratory_threshold"
        for minimum,label in CATEGORIES:
            if maximum>=minimum:category=label;break
        vals.append({"gene":fields[1],"maximum_delta_score":maximum,
            "scores":scores,"prediction_band":category,
            "source_model":"SpliceAI","splicing_effect_verified_by_RNA":False})
    return vals

def _query(file,chrom,pos):
    from .vep_local import wsl_path
    proc=subprocess.run([*wsl_command(),"tabix",wsl_path(file),f"{chrom}:{pos}-{pos}"],
       capture_output=True,text=True,errors="replace",timeout=50)
    if proc.returncode:raise RuntimeError("V289_SPLICEAI_LOCAL_INDEX_QUERY_FAILED:"+proc.stderr[-120:])
    return proc.stdout.splitlines()

def run(private,outdir,manifest_path=None,query=_query):
    from .expert_genomics import safe_hash
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    shards,status=validate_manifest(manifest_path)
    results=[];checked=0
    for gene,vs in sorted((private.get("genes") or {}).items()):
        for v in vs:
            ref="SHA256:"+safe_hash({**v,"gene":gene})
            scores=[]
            if shards:
                for shard in shards:
                    chrom=str(v.get("chrom") or "").removeprefix("chr")
                    pos=v.get("pos")
                    if not chrom or not isinstance(pos,int):continue
                    try:lines=query(shard["file"],chrom,pos)
                    except (RuntimeError,subprocess.TimeoutExpired,OSError):continue
                    for line in lines:
                        fields=line.split("\t")
                        if len(fields)<8:continue
                        if fields[1]!=str(pos) or fields[3].upper()!=str(v.get("ref") or "").upper():continue
                        if str(v.get("alt") or "").upper() not in fields[4].upper().split(","):continue
                        info={}
                        for item in fields[7].split(";"):
                            if "=" in item:
                                a,b=item.split("=",1);info[a]=b
                        scores.extend(parse_score_text(info.get("SpliceAI") or info.get("SPLICEAI"),
                            str(v.get("alt") or ""),gene))
                    checked+=1
            result={"variant_ref":ref,"gene":gene,
                "state":"MATCHED_SOURCE_SCORE" if scores else (
                   "NO_QUALIFIED_SCORE_IN_LOCAL_REFERENCE" if shards else status["state"]),
                "scores":scores,
                "maximum_delta_score":max((x["maximum_delta_score"] for x in scores),default=None),
                "nucleotide_coordinates_in_public_summary":False,
                "RNA_level_validation_required":True,
                "patient_specific_splicing_rescue_demonstrated":False,
                "multiple_transcripts_can_have_distinct_outcomes":True}
            results.append(result)
    summary={**status,"schema":SCHEMA,
       "variant_records":len(results),"successful_source_queries":checked,
       "variant_scores_resolved":sum(x["maximum_delta_score"] is not None for x in results),
       "score_absence_is_not_prediction_of_normal_splicing":True,
       "patient_alleles_sent_to_external_apis":False}
    (out/"SPLICEAI_STATUS.json").write_text(json.dumps(summary,indent=2)+"\n")
    with (out/"SPLICEAI_PRIVATE_MATCHES.jsonl").open("w",encoding="utf-8") as f:
        for x in results:f.write(json.dumps(x,sort_keys=True)+"\n")
    return summary,results
