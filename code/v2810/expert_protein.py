"""Source-grounded, offline UniProt residue-domain mapping for local genotypes.

Evidence is only region overlap, not proof of variant pathogenicity or drug rescue.
"""
from __future__ import annotations
from pathlib import Path
import json,re

SCHEMA="V288_PROTEIN_CONTEXT_1"
FEATURE_TYPES={
 "Domain","Region","Repeat","Binding site","Active site","Site",
 "Metal binding","Motif","DNA binding","Disulfide bond",
 "Modified residue","Transmembrane","Signal peptide",
 "Coiled coil","Zinc finger",
}
AA_POS=re.compile(r"\bp\.\(?[A-Za-z]{3}(\d+)(?:[A-Za-z]{3}|Ter|\*|=|fs)")
def parse_protein_position(hgvsp):
    if not hgvsp:return None
    m=AA_POS.search(str(hgvsp))
    return int(m.group(1)) if m else None

def position(x):
    if x is None:return None
    if isinstance(x,(int,float)) and int(x)==x:return int(x)
    if isinstance(x,dict):
        try:return int(x.get("value"))
        except (ValueError,TypeError):return None
    return None

def source_record(gene,public_root):
    p=Path(public_root)/"uniprot"
    found=sorted(p.glob(str(gene).upper()+"_*.json"))
    if not found:return None,{"state":"local_source_missing"}
    for file in found:
        try:j=json.loads(file.read_text(encoding="utf-8"))
        except (OSError,ValueError):continue
        name={str(x.get("geneName",{}).get("value","")).upper()
              for x in j.get("genes",[])}
        tax=(j.get("organism") or {}).get("taxonId")
        if gene.upper() not in name or tax not in (9606,"9606"):
            continue
        return j,{"state":"local_uniprot_human_gene_match",
                  "uniprot_acc":j.get("primaryAccession"),
                  "source_file_sha256_not_exposed":True}
    return None,{"state":"local_uniprot_gene_or_species_mismatch"}

def annotate_gene(variants,source):
    out=[]
    if not source:
        for v in variants:out.append({"variant_ref":v.get("variant_ref"),
                "position":parse_protein_position(v.get("hgvsp")),
                "feature_status":"unavailable", "domain_overlap_is_functional_proof":False})
        return out
    features=source.get("features",[])
    for v in variants:
        p=parse_protein_position(v.get("hgvsp"))
        hits=[]
        if p is not None:
            for feature in features:
                typ=str(feature.get("type") or "")
                if typ not in FEATURE_TYPES:continue
                loc=feature.get("location") or {}
                st=position(loc.get("start"));ed=position(loc.get("end"))
                if st is None or ed is None:continue
                if st<=p<=ed:
                    hits.append({"feature_type":typ,
                        "description":str(feature.get("description") or "")[:180],
                        "start":st,"end":ed,
                        "sequence_coordinate_reference":"UniProt canonical sequence",
                        "interpretation":"position_overlap_only"})
        out.append({"variant_ref":v.get("variant_ref"),
          "position":p,"feature_status":"overlap_observed" if hits else "no_known_overlap_or_unmapped",
          "feature_hits":hits[:20],
          "transcript_isoform_to_uniprot_alignment_validated":False,
          "domain_overlap_is_functional_proof":False})
    return out

def build_from_private(private,model,public_root):
    from .expert_genomics import safe_hash
    result={"schema":SCHEMA,"sources_only_human_gene_verified":True,"genes":{},
      "external_patient_variant_export":False,
      "protein_structural_effect_predicted":False}
    for gene,variants in sorted((private.get("genes") or {}).items()):
        j,meta=source_record(gene,public_root)
        synthetic=[{"variant_ref":"SHA256:"+safe_hash({**v,"gene":gene}),
                  "hgvsp":v.get("hgvsp")} for v in variants]
        result["genes"][gene]={"source":meta,
          "annotated_protein_positions":annotate_gene(synthetic,j),
          "limitation":"UniProt canonical sequence position may differ from variant transcript isoform; requires alignment. Region overlap does not establish loss of function."}
    return result

def save_report(result,destination):
    root=Path(destination);root.mkdir(parents=True,exist_ok=True)
    (root/"UNIPROT_PROTEIN_CONTEXT.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return result
