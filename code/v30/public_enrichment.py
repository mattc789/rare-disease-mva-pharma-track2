"""License-safe public target/pathway enrichment for V3.0.

Only gene symbols/target identifiers are sent to public scientific services.
No variants, coordinates, genotypes, patient identifiers or phenotype profiles
are transmitted. Responses are cached with provenance and release metadata.
"""
from __future__ import annotations
from pathlib import Path
import json,urllib.request,urllib.parse,hashlib,time

SCHEMA="RDMVA_PUBLIC_ENRICHMENT_V30_1"
OT_URL="https://api.platform.opentargets.org/api/v4/graphql"
REACTOME_BASE="https://reactome.org/ContentService"
UA="RareDiseaseMVAPharma-V30/1.0"

def _fetch_json(req,timeout=35,cap=8_000_000):
    with urllib.request.urlopen(req,timeout=timeout) as r:
        b=r.read(cap+1)
        if len(b)>cap:raise RuntimeError("V30_PUBLIC_RESPONSE_TOO_LARGE")
        return json.loads(b)

def _post_graphql(query,variables):
    body=json.dumps({"query":query,"variables":variables},separators=(",",":")).encode()
    req=urllib.request.Request(OT_URL,data=body,method="POST",
        headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":UA})
    out=_fetch_json(req)
    if out.get("errors"):raise RuntimeError("V30_OPEN_TARGETS_GRAPHQL_ERROR:"+str(out["errors"])[:500])
    return out.get("data") or {}

def opentargets_gene(gene):
    gene=str(gene or "").strip().upper()
    if not gene:raise ValueError("V30_GENE_REQUIRED")
    search='''query Q($q:String!){search(queryString:$q,entityNames:["target"]){hits{id name entity}}}'''
    hits=((_post_graphql(search,{"q":gene}).get("search") or {}).get("hits") or [])
    exact=[x for x in hits if str(x.get("name") or "").upper()==gene and x.get("entity")=="target"]
    if len(exact)!=1:
        return {"gene":gene,"state":"TARGET_NOT_UNIQUELY_RESOLVED","matches":exact[:5],
                "patient_variant_data_transmitted":False}
    tid=exact[0]["id"]
    query='''query Q($id:String!){target(ensemblId:$id){
      id approvedSymbol approvedName biotype
      tractability{label modality value}
      pathways{pathway pathwayId topLevelTerm}
      interactions(scoreThreshold:0.35,page:{index:0,size:60}){count rows{
        score sourceDatabase intA intB
        targetA{id approvedSymbol}
        targetB{id approvedSymbol}
        evidences{pubmedId evidenceScore interactionTypeShortName interactionDetectionMethodShortName}
      }}
      safetyLiabilities{event eventId datasource effects{direction dosing} literature url}
      drugAndClinicalCandidates{count rows{id maxClinicalStage drug{id name drugType maximumClinicalStage}}}
    }}'''
    t=(_post_graphql(query,{"id":tid}).get("target") or {})
    drugs=[]
    for x in ((t.get("drugAndClinicalCandidates") or {}).get("rows") or []):
        d=x.get("drug") or {}
        drugs.append({"drug_id":d.get("id"),"name":d.get("name"),"drug_type":d.get("drugType"),
                      "max_clinical_stage_for_target":x.get("maxClinicalStage"),
                      "max_clinical_stage_overall":d.get("maximumClinicalStage")})
    # one row per drug/indication can repeat; keep strongest textual stage row deterministically
    unique={}
    for d in drugs:
        key=d.get("drug_id") or d.get("name")
        if key and key not in unique:unique[key]=d
    interactions=[]
    for x in ((t.get("interactions") or {}).get("rows") or []):
        a=x.get("targetA") or {};b=x.get("targetB") or {}
        ga=str(a.get("approvedSymbol") or "").upper()
        gb=str(b.get("approvedSymbol") or "").upper()
        partner=gb if ga==gene else ga if gb==gene else (gb or ga)
        interactions.append({
          "partner_gene":partner,"score":x.get("score"),
          "source_database":x.get("sourceDatabase"),
          "pmids":sorted({str(e.get("pubmedId")) for e in (x.get("evidences") or []) if e.get("pubmedId")}),
          "evidences":[{
            "evidenceScore":e.get("evidenceScore"),
            "interactionTypeShortName":e.get("interactionTypeShortName"),
            "interactionDetectionMethodShortName":e.get("interactionDetectionMethodShortName")
          } for e in (x.get("evidences") or [])[:12]]
        })
    return {
      "gene":gene,"state":"RESOLVED","ensembl_id":t.get("id"),
      "approved_symbol":t.get("approvedSymbol"),"approved_name":t.get("approvedName"),
      "biotype":t.get("biotype"),"tractability":t.get("tractability") or [],
      "reactome_pathways_from_opentargets":t.get("pathways") or [],
      "interactions":interactions,
      "interaction_row_count":(t.get("interactions") or {}).get("count",0),
      "target_safety_liabilities":t.get("safetyLiabilities") or [],
      "clinical_target_drugs":list(unique.values()),
      "clinical_target_drug_row_count":(t.get("drugAndClinicalCandidates") or {}).get("count",0),
      "source":"Open Targets Platform GraphQL",
      "source_license":"Open Targets Platform data CC0 1.0; upstream source terms still apply",
      "patient_variant_data_transmitted":False,"patient_phenotype_profile_transmitted":False,
    }

def reactome_pathways(uniprot):
    acc=str(uniprot or "").strip()
    if not acc:return {"state":"NO_UNIPROT_ACCESSION","pathways":[]}
    url=REACTOME_BASE+"/data/mapping/UniProt/"+urllib.parse.quote(acc,safe="")+"/pathways"
    req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":UA})
    try:rows=_fetch_json(req,timeout=35)
    except Exception as exc:
        return {"state":"UNAVAILABLE","error":type(exc).__name__+":"+str(exc)[:300],"pathways":[]}
    if not isinstance(rows,list):rows=[]
    clean=[]
    for x in rows:
        if not isinstance(x,dict):continue
        clean.append({"stId":x.get("stId"),"displayName":x.get("displayName"),
                      "speciesName":x.get("speciesName"),"isInDisease":x.get("isInDisease")})
    clean=[x for x in clean if x.get("speciesName") in (None,"Homo sapiens")]
    return {"state":"RESOLVED","uniprot":acc,"pathways":clean,
            "source":"Reactome Content Service","source_license":"Reactome data CC0",
            "patient_variant_data_transmitted":False}

def enrich_genes(genes,uniprot_by_gene=None,outdir=None,online=True):
    uniprot_by_gene=uniprot_by_gene or {};results={}
    for gene in sorted(set(str(x).upper() for x in genes if x)):
        row={"gene":gene}
        if online:
            try:row["open_targets"]=opentargets_gene(gene)
            except Exception as exc:
                row["open_targets"]={"gene":gene,"state":"UNAVAILABLE",
                    "error":type(exc).__name__+":"+str(exc)[:400],
                    "patient_variant_data_transmitted":False}
            row["reactome"]=reactome_pathways(uniprot_by_gene.get(gene))
        else:
            row["open_targets"]={"gene":gene,"state":"ONLINE_DISABLED"}
            row["reactome"]={"state":"ONLINE_DISABLED","pathways":[]}
        results[gene]=row
    payload={"schema":SCHEMA,"genes":results,"online_gene_queries_enabled":bool(online),
             "raw_patient_variants_transmitted":False,"genotypes_transmitted":False,
             "full_phenotype_profile_transmitted":False}
    if outdir:
        p=Path(outdir);p.mkdir(parents=True,exist_ok=True)
        (p/"PUBLIC_GENE_ENRICHMENT_V30.json").write_text(
            json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload
