"""DGIdb v5 discovery-only adapter for V3.1.

DGIdb software is MIT, but imported interaction claims can inherit restrictions
from upstream sources. Therefore results are discovery pointers only.
"""
from __future__ import annotations
from pathlib import Path
import json,urllib.request

SCHEMA="RDMVA_DGIDB_DISCOVERY_V31_1"
URL="https://dgidb.org/api/graphql"
UA="RareDiseaseMVAPharma-V31/1.0"

QUERY='''query Q($genes:[String!],$n:Int!){interactions(geneNames:$genes,first:$n){totalCount nodes{
 id interactionScore evidenceScore
 drug{name approved conceptId}
 gene{name conceptId}
 interactionTypes{type directionality}
 publications{pmid citation}
 sources{sourceDbName fullName license licenseLink pmid}
}}}'''

def query_genes(genes,limit=200):
    body=json.dumps({"query":QUERY,"variables":{"genes":list(genes),"n":int(limit)}}).encode()
    req=urllib.request.Request(URL,data=body,method="POST",
      headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":UA})
    with urllib.request.urlopen(req,timeout=35) as r:
        raw=json.loads(r.read(5_000_000))
    if raw.get("errors"):raise RuntimeError("V31_DGIDB_GRAPHQL:"+str(raw["errors"])[:500])
    c=((raw.get("data") or {}).get("interactions") or {})
    rows=[]
    for x in c.get("nodes") or []:
        sources=[]
        all_source_terms_compatible=True
        for s in x.get("sources") or []:
            lic=s.get("license")
            sources.append({
              "source_db":s.get("sourceDbName"),"source_name":s.get("fullName"),
              "license":lic,"license_link":s.get("licenseLink"),"pmid":s.get("pmid")})
            if not lic:all_source_terms_compatible=False
        rows.append({
          "drug":(x.get("drug") or {}).get("name"),
          "drug_approved_dgidb":(x.get("drug") or {}).get("approved"),
          "drug_concept_id":(x.get("drug") or {}).get("conceptId"),
          "gene":(x.get("gene") or {}).get("name"),
          "gene_concept_id":(x.get("gene") or {}).get("conceptId"),
          "interaction_score":x.get("interactionScore"),
          "evidence_score":x.get("evidenceScore"),
          "interaction_types":x.get("interactionTypes") or [],
          "pmids":sorted({str(p.get("pmid")) for p in x.get("publications") or [] if p.get("pmid")}),
          "sources":sources,
          "all_source_licenses_explicit":all_source_terms_compatible,
          "submission_evidence_eligible":False,
        })
    return {"total_count":c.get("totalCount",0),"rows":rows}

def build(genes,outdir,online=True):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    try:data=query_genes(sorted(set(genes))) if online and genes else {"total_count":0,"rows":[]}
    except Exception as exc:data={"total_count":0,"rows":[],"error":type(exc).__name__+":"+str(exc)[:400]}
    payload={
      "schema":SCHEMA,"queried_genes":sorted(set(genes)),
      "interaction_count":len(data.get("rows") or []),
      "interactions":data.get("rows") or [],
      "upstream_total_count":data.get("total_count"),
      "error":data.get("error"),
      "use_class":"DISCOVERY_INDEX_ONLY",
      "dgidb_software_license":"MIT",
      "aggregated_claim_license_status":"SOURCE_SPECIFIC",
      "dgidb_only_claim_can_support_submission":False,
      "primary_source_verification_required":True,
      "patient_variant_data_transmitted":False,
      "clinical_recommendation":False,
    }
    (out/"DGIDB_DISCOVERY_V31.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload


def _interaction_sign(types):
    vals=[]
    for x in types or []:
        typ=str(x.get("type") or "").lower()
        direction=str(x.get("directionality") or "").lower()
        text=typ+" "+direction
        if any(k in text for k in ("negative","inhibit","antagon","block","suppress","decreas")):
            vals.append(-1)
        elif any(k in text for k in ("positive","activat","agon","stimulat","induc","increas")):
            vals.append(1)
    return vals[0] if vals and all(v==vals[0] for v in vals) else 0

def discover(routes,entities,outdir):
    from v2810.core import norm_name
    gene_routes={}
    for r in routes or []:
        for t in r.get("targets") or []:
            g=str(t.get("gene") or "").upper()
            if g:gene_routes.setdefault(g,[]).append(r)
    payload=build(sorted(gene_routes),outdir,online=True)
    emap={}
    for e in entities or []:
        names=[e.get("canonical_label"),*(e.get("member_labels") or [])]
        for n in names:
            if n:emap.setdefault(norm_name(n),e)
    hints=[]
    for x in payload.get("interactions",[]):
        drug=x.get("drug");gene=str(x.get("gene") or "").upper()
        if not drug or gene not in gene_routes:continue
        entity=emap.get(norm_name(drug))
        if not entity:continue
        sign=_interaction_sign(x.get("interaction_types"))
        for r in gene_routes[gene]:
            desired=1 if r.get("desired_action")=="increase" else -1 if r.get("desired_action")=="decrease" else 0
            compatible=bool(sign and desired and sign==desired)
            hints.append({
              "candidate":entity.get("canonical_label"),
              "entity_id":entity.get("entity_id"),
              "route_id":r.get("route_id"),"gene":gene,
              "interaction_sign":sign,"desired_sign":desired,
              "direction_compatible":compatible,
              "interaction_score":x.get("interaction_score"),
              "evidence_score":x.get("evidence_score"),
              "interaction_types":x.get("interaction_types") or [],
              "pmids":x.get("pmids") or [],
              "sources":x.get("sources") or [],
              "source_licenses_explicit":x.get("all_source_licenses_explicit"),
              "submission_evidence_eligible":False})
    payload["hints"]=hints
    payload["hint_count"]=len(hints)
    payload["direction_matched_hint_count"]=sum(x["direction_compatible"] for x in hints)
    Path(outdir,"DGIDB_DISCOVERY_V31.json").write_text(
      json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return payload

def merge_pairs(pairs,payload,routes,entities):
    rmap={r.get("route_id"):r for r in routes or []}
    emap={e.get("entity_id"):e for e in entities or []}
    existing={(x.get("ingredient"),x.get("route_id")) for x in pairs or []}
    out=list(pairs or []);added=0
    for h in payload.get("hints",[]):
        if not (h.get("direction_compatible") or h.get("direction_match")):continue
        key=(h.get("candidate"),h.get("route_id"))
        if key in existing:continue
        route=rmap.get(h.get("route_id"))
        if not route:continue
        entity=emap.get(h.get("entity_id")) or {}
        out.append({
          "schema":"RDMVA_DGIDB_POINTER_PAIR_V31_1",
          "ingredient":h.get("candidate"),"entity_id":h.get("entity_id"),
          "scientific_group_id":entity.get("scientific_group_id"),
          "member_labels":entity.get("member_labels") or [h.get("candidate")],
          "route_id":h.get("route_id"),"route":route,
          "tier":"D","score":0.25,
          "providers":["DGIdb_discovery_pointer"],
          "signals":[{
            "route_id":h.get("route_id"),"provider":"DGIdb_discovery_pointer",
            "interaction_score":h.get("interaction_score"),
            "evidence_score":h.get("evidence_score"),
            "pmids":h.get("pmids"),"sources":h.get("sources"),
            "direct_target_assay_verified":False,"functional_rescue":False,
            "submission_evidence_eligible":False}],
          "combo_only":bool(entity.get("combo_only")),
          "global_approved_as_single":bool(entity.get("global_approved_as_single") or
                                            entity.get("solo_approved") or
                                            entity.get("biologic_licensed") or
                                            not entity.get("combo_only",False)),
          "modality":entity.get("modality"),
          "publication_strength":"discovery_pointer_only",
          "evidence_blockers":[
            "dgidb_aggregator_requires_primary_source_verification",
            "no_verified_direct_target_modulation",
            "no_same_drug_disease_functional_rescue"],
          "shortest_network_hops":route.get("network_hops"),
          "drug_to_route_hops":0,
          "route_to_case_hops":int(route.get("network_hops") or 0),
          "scientific_gate":route.get("scientific_gate") or {},
          "promotion_eligible_validated_lead":False,
          "clinical_recommendation":False})
        existing.add(key);added+=1
    out.sort(key=lambda x:({"A":0,"B":1,"C":2,"D":3}.get(x.get("tier"),4),
                           -float(x.get("score") or 0),str(x.get("ingredient")),str(x.get("route_id"))))
    return out,{"added_pairs":added,"aggregator_only_pairs_are_submission_eligible":False}
