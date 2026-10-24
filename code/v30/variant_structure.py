"""Variant-specific structural context using AlphaFold DB / AlphaMissense.

AFDB predictions are used as theoretical structural context only. AlphaMissense
is optional computational evidence. Neither source establishes diagnosis,
causality, clinical efficacy, or treatment suitability.
"""
from __future__ import annotations
from pathlib import Path
import gzip,hashlib,json,os,re,urllib.request

SCHEMA="RDMVA_VARIANT_STRUCTURE_V31_1"
AFDB_API="https://alphafold.ebi.ac.uk/api/prediction/"
AM_URL="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
AM_MD5="9fd167735f16a1b87da6eb3e4c25fcb5"
HGVSP_RX=re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)")

AA3={"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
     "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
     "Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}

def _json(url,timeout=30,cap=4_000_000):
    req=urllib.request.Request(url,headers={"Accept":"application/json",
                                           "User-Agent":"RareDiseaseMVAPharma-V31/1.0"})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        b=r.read(cap+1)
        if len(b)>cap:raise RuntimeError("V31_AFDB_RESPONSE_TOO_LARGE")
        return json.loads(b)

def _bytes(url,timeout=60,cap=30_000_000):
    req=urllib.request.Request(url,headers={"User-Agent":"RareDiseaseMVAPharma-V31/1.0"})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        b=r.read(cap+1)
        if len(b)>cap:raise RuntimeError("V31_AFDB_STRUCTURE_TOO_LARGE")
        return b

def afdb_structure(accession,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    meta=_json(AFDB_API+accession)
    rows=meta if isinstance(meta,list) else [meta]
    if not rows:return {"state":"NOT_FOUND","accession":accession}
    m=rows[0] or {}
    pdb_url=m.get("pdbUrl")
    if not pdb_url:return {"state":"NO_PDB_URL","accession":accession,"metadata":m}
    pdb_path=out/f"{accession}.pdb"
    if not pdb_path.is_file():pdb_path.write_bytes(_bytes(pdb_url))
    residues={}
    for line in pdb_path.read_text(encoding="utf-8",errors="replace").splitlines():
        if not line.startswith("ATOM"):continue
        if line[12:16].strip()!="CA":continue
        try:
            pos=int(line[22:26].strip())
            b=float(line[60:66].strip())
            x=float(line[30:38]);y=float(line[38:46]);z=float(line[46:54])
        except Exception:continue
        residues[pos]={"position":pos,"residue3":line[17:20].strip(),
                       "plddt":b,"ca":[x,y,z]}
    return {"state":"RESOLVED","accession":accession,
            "model_version":m.get("latestVersion"),"pdb_url":pdb_url,
            "pae_url":m.get("paeDocUrl"),"residue_count":len(residues),
            "residues":residues,
            "license":"CC BY 4.0","commercial_use_permitted":True,
            "clinical_use_validated":False}

def _parse_hgvsp(h):
    if not h:return None
    m=HGVSP_RX.search(str(h))
    if not m:return None
    return {"ref":AA3.get(m.group(1)), "position":int(m.group(2)),
            "alt":"*" if m.group(3) in ("Ter","*") else AA3.get(m.group(3))}

def structural_variant_context(private,protein_context,accession_by_gene,outdir,online=True):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    result={}
    for gene,variants in sorted((private.get("genes") or {}).items()):
        acc=accession_by_gene.get(gene)
        af={"state":"NO_ACCESSION"}
        if acc and online:
            try:af=afdb_structure(acc,out/"afdb")
            except Exception as exc:af={"state":"UNAVAILABLE","error":type(exc).__name__+":"+str(exc)[:300]}
        vrows=[]
        residues=af.get("residues") or {}
        for v in variants:
            pv=_parse_hgvsp(v.get("hgvsp"))
            if not pv or pv.get("alt")=="*":continue
            rr=residues.get(pv["position"]) or residues.get(str(pv["position"])) or {}
            plddt=rr.get("plddt")
            if plddt is None:conf="UNKNOWN"
            elif plddt>=90:conf="VERY_HIGH"
            elif plddt>=70:conf="CONFIDENT"
            elif plddt>=50:conf="LOW"
            else:conf="VERY_LOW_DISORDER_POSSIBLE"
            vrows.append({
              "variant_ref":hashlib.sha256(
                f"{v.get('chrom')}:{v.get('pos')}:{v.get('ref')}:{v.get('alt')}".encode()
              ).hexdigest()[:16],
              "protein_variant":f"{pv['ref'] or '?'}{pv['position']}{pv['alt'] or '?'}",
              "protein_position":pv["position"],"alphafold_plddt":plddt,
              "alphafold_confidence_class":conf,
              "structure_prediction_is_functional_proof":False,
            })
        result[gene]={"gene":gene,"uniprot":acc,
                      "afdb_state":af.get("state"),
                      "afdb_model_version":af.get("model_version"),
                      "variants":vrows}
    payload={"schema":SCHEMA,"genes":result,
             "alphafold_db_license":"CC BY 4.0",
             "alphafold_db_commercial_use_permitted":True,
             "alphafold_clinical_use_validated":False,
             "raw_variant_coordinates_in_public_output":False,
             "clinical_recommendation":False}
    (out/"VARIANT_STRUCTURE_CONTEXT_V31.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload

def _resource_root(project_root):
    cfg=Path(project_root)/"config/local_public_storage_v2810.json"
    if cfg.is_file():
        try:
            j=json.loads(cfg.read_text(encoding="utf-8"))
            root=Path(j.get("windows_root") or "")
            if root.anchor and root.parent.exists():return root/"alphamissense"
        except Exception:pass
    return Path.home()/".rdmva_resources"/"alphamissense"

def ensure_alphamissense(project_root,download=False):
    root=_resource_root(project_root);root.mkdir(parents=True,exist_ok=True)
    p=root/"AlphaMissense_hg38.tsv.gz"
    if not p.is_file() and download:
        tmp=p.with_suffix(".download")
        with urllib.request.urlopen(AM_URL,timeout=60) as src,tmp.open("wb") as dst:
            while True:
                chunk=src.read(4*1024*1024)
                if not chunk:break
                dst.write(chunk)
        md5=hashlib.md5()
        with tmp.open("rb") as f:
            for chunk in iter(lambda:f.read(8*1024*1024),b""):md5.update(chunk)
        if md5.hexdigest()!=AM_MD5:
            tmp.unlink(missing_ok=True);raise RuntimeError("V31_ALPHAMISSENSE_MD5_MISMATCH")
        tmp.replace(p)
    return {"state":"READY" if p.is_file() else "NOT_INSTALLED",
            "path":str(p) if p.is_file() else None,
            "source_url":AM_URL,"expected_md5":AM_MD5,
            "license_status":"CONFLICTING_OFFICIAL_NOTICES",
            "distributed_file_header_license":"CC BY-NC-SA 4.0",
            "current_web_documentation_license":"CC BY 4.0",
            "commercial_use_permitted":False,
            "challenge_baseline_enabled":False,
            "clinical_use_validated":False}

def lookup_alphamissense(private,resource,outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    wanted={}
    for gene,vs in (private.get("genes") or {}).items():
        for v in vs:
            cons=str(v.get("consequence") or "").lower()
            if "missense" not in cons:continue
            key=("chr"+str(v.get("chrom")).removeprefix("chr"),int(v.get("pos")),
                 str(v.get("ref")),str(v.get("alt")))
            parsed=_parse_hgvsp(v.get("hgvsp"))
            expected=None
            if parsed and parsed.get("ref") and parsed.get("alt") and parsed.get("alt")!="*":
                expected=f"{parsed['ref']}{parsed['position']}{parsed['alt']}"
            wanted[key]={"gene":gene,"expected_protein_variant":expected,
                         "case_hgvsp":v.get("hgvsp")}
    hits={}
    if resource.get("state")=="READY" and wanted:
        path=Path(resource["path"])
        def chr_rank(c):
            s=str(c).removeprefix("chr")
            if s.isdigit():return int(s)
            return {"X":23,"Y":24,"M":25,"MT":25}.get(s,99)
        max_rank=max(chr_rank(k[0]) for k in wanted)
        max_pos_on_last=max(k[1] for k in wanted if chr_rank(k[0])==max_rank)
        with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
            header=None
            for line in f:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):header=line.lstrip("#").rstrip("\n").split("\t")
                    continue
                p=line.rstrip("\n").split("\t")
                if len(p)<10:continue
                try:
                    pos=int(p[1]);rank=chr_rank(p[0])
                    key=(p[0],pos,p[2],p[3])
                except:continue
                # Official GRCh38 file is coordinate sorted. Stop once all
                # possible requested loci are behind the current row.
                if rank>max_rank or (rank==max_rank and pos>max_pos_on_last):
                    break
                if key in wanted:
                    meta=wanted[key]
                    hits[key]={"gene":meta["gene"],"uniprot_id":p[5],
                               "transcript_id":p[6],"protein_variant":p[7],
                               "am_pathogenicity":float(p[8]),"am_class":p[9]}
                    if len(hits)==len(wanted):break
    rows=[]
    for key,meta in wanted.items():
        gene=meta["gene"];h=hits.get(key)
        expected=meta.get("expected_protein_variant")
        observed=(h or {}).get("protein_variant")
        concordant=bool(h and expected and observed and expected==observed)
        if h and expected and not concordant:state="MATCHED_GENOMIC_ISOFORM_DISCORDANT"
        elif h and concordant:state="MATCHED_CONSISTENT"
        elif h:state="MATCHED_GENOMIC_NO_CASE_PROTEIN_MAPPING"
        else:state="NOT_FOUND_OR_RESOURCE_MISSING"
        rows.append({"gene":gene,"variant_ref":hashlib.sha256(
            f"{key[0]}:{key[1]}:{key[2]}:{key[3]}".encode()).hexdigest()[:16],
            "state":state,"case_hgvsp":meta.get("case_hgvsp"),
            "expected_protein_variant":expected,
            "protein_mapping_concordant":concordant,
            "prediction":h if concordant else None,
            "raw_prediction_for_mapping_audit":h if h and not concordant else None,
            "prediction_is_causality_proof":False})
    payload={"schema":"RDMVA_ALPHAMISSENSE_V31_2","resource":resource,
             "variants":rows,
             "matched":sum(x["state"]=="MATCHED_CONSISTENT" for x in rows),
             "mapping_discordant":sum(x["state"]=="MATCHED_GENOMIC_ISOFORM_DISCORDANT" for x in rows),
             "license_status":"CONFLICTING_OFFICIAL_NOTICES",
             "distributed_file_header_license":"CC BY-NC-SA 4.0",
             "current_web_documentation_license":"CC BY 4.0",
             "challenge_baseline_enabled":False,
             "clinical_use_validated":False,
             "clinical_recommendation":False}
    (out/"ALPHAMISSENSE_V31.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload


def annotate_objectives(objectives,structure_context,alphamissense):
    """Attach variant-specific computational context without changing evidence tier."""
    am_by={}
    for x in alphamissense.get("variants",[]) if alphamissense else []:
        am_by.setdefault(x.get("gene"),[]).append(x)
    st_by=(structure_context or {}).get("genes") or {}
    rows=[]
    for o0 in objectives.get("objectives",[]):
        o=dict(o0);gene=o.get("gene")
        sv=[]
        for v in (st_by.get(gene) or {}).get("variants",[]):
            if v.get("alphafold_plddt") is not None:
                sv.append({"protein_position":v.get("protein_position"),
                           "plddt":v.get("alphafold_plddt"),
                           "confidence":v.get("alphafold_confidence_class")})
        av=[]
        for x in am_by.get(gene,[]):
            p=(x.get("prediction") or {})
            if p:
                av.append({"protein_variant":p.get("protein_variant"),
                           "am_pathogenicity":p.get("am_pathogenicity"),
                           "am_class":p.get("am_class")})
        o["variant_specific_computational_context"]={
          "alphafold":sv,"alphamissense":av,
          "context_is_functional_proof":False,
          "context_changes_evidence_tier":False}
        rows.append(o)
    return {**objectives,"objectives":rows,
            "variant_specific_context_attached":True,
            "variant_specific_context_is_clinical_validation":False}
