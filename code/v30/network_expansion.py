"""Open network/pathway expansion for V3.1.

Uses only gene/protein identifiers with Open Targets and Reactome. The output
is an unsigned discovery graph. It does NOT create pharmacologic direction by
itself; direction must come from signed causal evidence or source-bound
modifier literature.
"""
from __future__ import annotations
from pathlib import Path
import json,re,urllib.parse,urllib.request
from . import public_enrichment

SCHEMA="RDMVA_NETWORK_EXPANSION_V31_1"
UA="RareDiseaseMVAPharma-V31/1.0"
GENE_RX=re.compile(r"^[A-Z][A-Z0-9-]{1,19}$")

def _fetch(url,cap=8_000_000,timeout=30):
    req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":UA})
    with urllib.request.urlopen(req,timeout=timeout) as r:
        b=r.read(cap+1)
        if len(b)>cap:raise RuntimeError("V31_NETWORK_RESPONSE_TOO_LARGE")
        return json.loads(b)

def reactome_participants(pathway_id,max_rows=50):
    if not pathway_id:return []
    url=("https://reactome.org/ContentService/data/event/"+
         urllib.parse.quote(str(pathway_id),safe="")+
         "/participatingPhysicalEntities")
    try:rows=_fetch(url)
    except Exception:return []
    out=[];seen=set()
    for x in rows if isinstance(rows,list) else []:
        if not isinstance(x,dict):continue
        if x.get("speciesName") not in (None,"Homo sapiens"):continue
        names=list(x.get("name") or [])
        disp=str(x.get("displayName") or "").split(" [",1)[0].strip()
        if disp:names.append(disp)
        for name in names:
            g=str(name or "").upper().strip()
            if GENE_RX.fullmatch(g) and g not in seen:
                seen.add(g);out.append(g)
                if len(out)>=max_rows:return out
    return out

def _experimental_weight(interaction):
    ev=interaction.get("evidences") or []
    ex=[float(x.get("evidenceScore") or 0) for x in ev
        if str(x.get("interactionDetectionMethodShortName") or "").lower()=="experimental"]
    return max(ex or [0.0])

def build(disease_graph,uniprot_by_gene,outdir,online=True,
          max_interactions_per_gene=30,max_pathways_per_gene=6,
          max_participants_per_pathway=35):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    seeds=list(disease_graph.get("network_seed_genes") or [])
    enriched=public_enrichment.enrich_genes(
        seeds,uniprot_by_gene,out/"gene_context",online=online)
    nodes={};edges=[];seedset=set(seeds)
    for gene in seeds:
        nodes["GENE:"+gene]={"id":"GENE:"+gene,"type":"gene","seed":True}
        row=(enriched.get("genes") or {}).get(gene) or {}
        ot=row.get("open_targets") or {}
        interactions=list(ot.get("interactions") or [])
        interactions.sort(key=lambda x:(-float(x.get("score") or 0),
                                        -_experimental_weight(x),
                                        str(x.get("partner_gene") or "")))
        for x in interactions[:max_interactions_per_gene]:
            partner=str(x.get("partner_gene") or "").upper()
            if not partner or partner==gene:continue
            nid="GENE:"+partner
            nodes.setdefault(nid,{"id":nid,"type":"gene","seed":partner in seedset})
            edges.append({
              "source":"GENE:"+gene,"target":nid,
              "type":"open_targets_interaction_unsigned",
              "score":x.get("score"),
              "experimental_evidence_score":_experimental_weight(x),
              "source_database":x.get("source_database"),
              "pmids":x.get("pmids") or [],
              "direction_known":False,
            })
        pathways=(ot.get("reactome_pathways_from_opentargets") or [])[:max_pathways_per_gene]
        for p in pathways:
            pid=p.get("pathwayId") or p.get("pathway_id")
            if not pid:continue
            pn="PATHWAY:"+pid
            nodes[pn]={"id":pn,"type":"reactome_pathway",
                       "label":p.get("pathway"),"top_level":p.get("topLevelTerm"),
                       "license":"CC0"}
            edges.append({"source":"GENE:"+gene,"target":pn,
                          "type":"participates_in_pathway"})
            if online:
                for participant in reactome_participants(
                    pid,max_rows=max_participants_per_pathway):
                    if participant==gene:continue
                    gid="GENE:"+participant
                    nodes.setdefault(gid,{"id":gid,"type":"gene",
                                          "seed":participant in seedset})
                    edges.append({"source":pn,"target":gid,
                                  "type":"pathway_participant_unsigned",
                                  "direction_known":False})
    modifier_candidates={}
    for e in edges:
        if e["type"]=="open_targets_interaction_unsigned":
            gene=e["target"].split(":",1)[1]
            rec=modifier_candidates.setdefault(gene,{
              "gene":gene,"open_targets_links":0,"pathway_links":0,
              "max_interaction_score":0.0,"max_experimental_score":0.0,
              "seed_gene":gene in seedset})
            rec["open_targets_links"]+=1
            rec["max_interaction_score"]=max(rec["max_interaction_score"],
                                              float(e.get("score") or 0))
            rec["max_experimental_score"]=max(rec["max_experimental_score"],
                                               float(e.get("experimental_evidence_score") or 0))
        elif e["type"]=="pathway_participant_unsigned":
            gene=e["target"].split(":",1)[1]
            rec=modifier_candidates.setdefault(gene,{
              "gene":gene,"open_targets_links":0,"pathway_links":0,
              "max_interaction_score":0.0,"max_experimental_score":0.0,
              "seed_gene":gene in seedset})
            rec["pathway_links"]+=1
    modifiers=sorted(modifier_candidates.values(),
      key=lambda x:(x["seed_gene"],-x["max_experimental_score"],
                    -x["max_interaction_score"],-x["pathway_links"],x["gene"]))
    payload={
      "schema":SCHEMA,"seed_genes":seeds,
      "nodes":list(nodes.values()),"edges":edges,
      "modifier_candidates":modifiers,
      "modifier_candidate_count":len(modifiers),
      "direction_inferred_from_unsigned_network":False,
      "raw_variants_transmitted":False,"genotypes_transmitted":False,
      "full_phenotype_transmitted":False,
      "open_targets_data_license":"CC0 1.0",
      "reactome_data_license":"CC0",
      "clinical_recommendation":False,
    }
    (out/"NETWORK_EXPANSION_V31.json").write_text(
        json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload
