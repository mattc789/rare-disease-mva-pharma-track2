"""V3-native rescue-route generator.

The route universe is derived from explicit patient rescue objectives BEFORE
any drug screen. No drug names are used to construct the objectives/routes.
"""
from __future__ import annotations
from collections import defaultdict
import re
from v2810 import routes as v28_routes, science as v28_science

SCHEMA="RDMVA_NATIVE_ROUTES_V30_1"

def _route(**kw):
    r=dict(kw)
    r.setdefault("state","conditional")
    r.setdefault("score",10.0)
    r.setdefault("targets",[])
    return r

def build(private,objectives,cache_dir):
    edges=v28_routes.signor_edges()
    by_gene=defaultdict(list)
    for e in edges:
        by_gene[e.get("target_gene","").upper()].append(e)
    out=[]
    accessions={}
    genes=sorted(private.get("genes") or {})
    for gene in genes:
        try:accessions[gene]=v28_routes.accession_for_gene(gene,cache_dir)
        except Exception:accessions[gene]=None

    for obj in objectives.get("objectives",[]):
        if not obj.get("drug_search_eligible"):continue
        gene=obj["gene"];acc=accessions.get(gene)
        typ=obj["objective_type"];action=obj["desired_action"];var=obj["variable"]
        base={
          "case_gene":gene,"rescue_objective_id":obj["objective_id"],
          "rescue_objective_type":typ,"rescue_priority_tier":obj["priority_tier"],
          "objective_requirements":obj.get("requirements_before_promotion",[]),
        }

        # Direct primary-gene routes.
        if typ in ("restore_functional_gene_product","increase_residual_activity",
                   "stabilize_residual_protein","restore_complex_or_interaction_function"):
            variable={
              "restore_functional_gene_product":"abundance",
              "increase_residual_activity":"activity",
              "stabilize_residual_protein":"folding_stability",
              "restore_complex_or_interaction_function":"complex_function",
            }[typ]
            out.append(_route(**base,family="v30_direct_gene",target_concept=gene,
              targets=[{"gene":gene,"accession":acc}] if acc else [],
              desired_action="increase",variable=variable,
              mechanism_concept=f"{typ} for {gene}",score=18.0,
              requires_residual_product=True,source="v30_rescue_objective"))

        # RNA-level objectives remain explicit even when no protein target exists.
        if typ=="restore_full_length_translation":
            out.append(_route(**base,family="v30_translation_rescue",
              target_concept=f"{gene} premature-stop readthrough",targets=[],
              desired_action="increase",variable="full_length_translation",
              mechanism_concept=f"restore full-length {gene} translation",score=17.0,
              source="v30_rescue_objective"))
        elif typ=="reduce_nmd_conditionally":
            out.append(_route(**base,family="v30_nmd_rescue",
              target_concept=f"NMD of {gene} transcript",targets=[],
              desired_action="decrease",variable="nmd_of_implicated_transcript",
              mechanism_concept=f"reduce NMD only if it increases functional {gene} product",
              score=7.0,state="hold",source="v30_rescue_objective"))
        elif typ=="restore_productive_splicing":
            out.append(_route(**base,family="v30_splicing_rescue",
              target_concept=f"{gene} productive splicing",
              targets=[{"gene":gene,"accession":acc}] if acc else [],
              desired_action="increase",variable="productive_transcript",
              mechanism_concept=f"restore productive {gene} splicing",score=17.0,
              source="v30_rescue_objective"))

        # Signed upstream regulators are objective-specific. A generic SIGNOR
        # edge is relevant to residual activity, but must NOT be reused to
        # claim splicing, folding, complex rescue, etc.
        if acc:
            for e in edges:
                if e.get("target_acc")!=acc or not e.get("direct"):continue
                mech=(" ".join(str(e.get(k) or "") for k in ("effect","mechanism","sentence"))).lower()
                gene_lo=gene.lower()
                gene_stability=bool(
                    re.search(rf"\b{re.escape(gene_lo)}\b[^.;]{{0,45}}\b(?:protein )?stability\b",mech) or
                    re.search(rf"\bstabili(?:ze|zes|zed|zing|sation|zation)\b[^.;]{{0,35}}\b{re.escape(gene_lo)}\b",mech))
                degradation=bool(re.search(r"\b(?:degrad|ubiquit|proteas|proteol)",mech))
                expression=bool(re.search(r"\b(?:express|transcription|promoter|quantity by expression)",mech))
                signed_ok=(typ=="increase_residual_activity")
                if typ=="restore_functional_gene_product":
                    signed_ok=expression or degradation or gene_stability
                elif typ=="stabilize_residual_protein":
                    signed_ok=degradation or gene_stability
                if not signed_ok:continue
                desired="increase" if e["sign"]>0 else "decrease"
                out.append(_route(**base,family="v30_signed_upstream",
                  target_concept=e["source_gene"],
                  targets=[{"gene":e["source_gene"],"accession":e["source_acc"]}],
                  desired_action=desired,variable="activity",
                  mechanism_concept=f"{desired} {e['source_gene']} to support objective {typ} for {gene}",
                  score=13.0,source=e["signor_id"],pmid=e.get("pmid"),network_hops=1))

        # Downstream compensation is only relevant for restore/bypass objectives.
        if typ in ("compensate_downstream_deficit","bypass_or_compensate_downstream"):
            if acc:
                for e in edges:
                    if e.get("source_acc")!=acc or not e.get("direct"):continue
                    desired="increase" if e["sign"]>0 else "decrease"
                    out.append(_route(**base,family="v30_downstream_compensation",
                      target_concept=e["target_gene"],
                      targets=[{"gene":e["target_gene"],"accession":e["target_acc"]}],
                      desired_action=desired,variable="activity",
                      mechanism_concept=f"restore downstream state normally controlled by {gene}",
                      score=11.0,source=e["signor_id"],pmid=e.get("pmid"),network_hops=1))

        # Proteostasis objective can also be satisfied by inhibiting a direct
        # negative regulator of the primary protein, but it stays conditional.
        if typ=="reduce_causal_degradation" and acc:
            for e in edges:
                if e.get("target_acc")!=acc or not e.get("direct") or e.get("sign")>=0:continue
                out.append(_route(**base,family="v30_degradation_relief",
                  target_concept=e["source_gene"],
                  targets=[{"gene":e["source_gene"],"accession":e["source_acc"]}],
                  desired_action="decrease",variable="activity",
                  mechanism_concept=f"decrease negative regulator {e['source_gene']} to preserve {gene}",
                  score=14.0,source=e["signor_id"],pmid=e.get("pmid")))

    # Dedupe without losing objective provenance: same pharmacologic action can
    # cover multiple rescue objectives.
    merged={}
    for r in out:
        k=(r["case_gene"],r["family"],str(r["target_concept"]).upper(),
           r["desired_action"],r["variable"])
        if k not in merged:
            merged[k]=r;merged[k]["rescue_objective_ids"]=[r["rescue_objective_id"]]
            merged[k]["rescue_objective_types"]=[r["rescue_objective_type"]]
        else:
            x=merged[k]
            x["rescue_objective_ids"]=sorted(set(x["rescue_objective_ids"]+[r["rescue_objective_id"]]))
            x["rescue_objective_types"]=sorted(set(x["rescue_objective_types"]+[r["rescue_objective_type"]]))
            x["score"]=max(x["score"],r["score"])
    rows=sorted(merged.values(),key=lambda x:(-x["score"],x["case_gene"],x["family"],str(x["target_concept"])))
    for i,r in enumerate(rows,1):
        r["route_id"]=f"V30R{i:04d}";r["rank"]=i
    # Reuse conservative V2 scientific gates, then attach V3 provenance.
    gated,audit=v28_science.adjudicate_routes(rows)
    return {
      "schema":SCHEMA,"routes":gated,"audit":audit,
      "route_count":len(gated),
      "route_generation_uses_drug_names":False,
      "derived_before_drug_screen":True,
    }


def _bind_literature_objectives(route,objectives):
    """Bind a source-bound literature route only to semantically supported objectives."""
    gene=str(route.get("case_gene") or "").upper()
    text=" ".join(str(route.get(k) or "") for k in
                  ("target_concept","mechanism_concept")).lower()
    compatible=set()
    if re.search(r"\b(?:nmd|nonsense[- ]mediated decay)\b",text):
        compatible.add("reduce_nmd_conditionally")
    if re.search(r"\b(?:splic|spliceosome|exon|intron)\w*",text):
        compatible.add("restore_productive_splicing")
    if re.search(r"\b(?:readthrough|premature[- ]stop|stop codon|full[- ]length translation)\b",text):
        compatible.add("restore_full_length_translation")
    if re.search(r"\b(?:degrad|ubiquit|proteas|proteol|stabili|fold|chaperon|proteostasis)\w*",text):
        compatible.update({"stabilize_residual_protein","reduce_causal_degradation",
                           "restore_functional_gene_product"})
    if re.search(r"\b(?:abundance|expression|protein level|protein levels|half[- ]life)\b",text):
        compatible.add("restore_functional_gene_product")
    if re.search(r"\b(?:activity|activate|activation|catalytic|kinase function|enzymatic)\b",text):
        compatible.add("increase_residual_activity")
    if re.search(r"\b(?:complex|assembly|interaction|localization|localisation|kinetochore)\b",text):
        compatible.add("restore_complex_or_interaction_function")
    if re.search(r"\b(?:downstream|bypass|compensat|pathway|signaling|signalling|phenotype|functional rescue)\b",text):
        compatible.update({"compensate_downstream_deficit","bypass_or_compensate_downstream"})
    rows=[o for o in objectives.get("objectives",[])
          if str(o.get("gene") or "").upper()==gene
          and o.get("drug_search_eligible")
          and o.get("objective_type") in compatible]
    return rows


def merge_literature_discovered(native_bundle,discovered,objectives=None):
    """Add source-bound public-literature modifiers before the drug screen."""
    rows=list(native_bundle.get("routes") or [])
    seen={(r.get("case_gene"),str(r.get("target_concept") or "").upper(),
           r.get("desired_action"),r.get("variable")) for r in rows}
    for x in discovered or []:
        gene=x.get("case_gene")
        key=(gene,str(x.get("target_concept") or "").upper(),
             x.get("desired_action"),x.get("variable"))
        if not gene or key in seen:continue
        targets=x.get("targets") or []
        if not targets and x.get("target_gene") and x.get("target_accession"):
            targets=[{"gene":x["target_gene"],"accession":x["target_accession"]}]
        bound=_bind_literature_objectives({**x,"case_gene":gene},objectives or {"objectives":[]})
        rows.append(_route(
          family="v30_literature_discovered",case_gene=gene,
          target_concept=x.get("target_concept"),targets=targets,
          desired_action=x.get("desired_action"),variable=x.get("variable"),
          mechanism_concept=x.get("mechanism_concept") or x.get("target_concept"),
          score=float(x.get("score",12.0)),state=x.get("state","conditional"),
          source=x.get("source_id"),evidence=x.get("evidence",[]),
          route_direction_verification=x.get("route_direction_verification"),
          evidence_context=x.get("evidence_context"),
          rescue_objective_ids=[o["objective_id"] for o in bound],
          rescue_objective_types=[o["objective_type"] for o in bound],
          objective_binding_state="SEMANTIC_SOURCE_BOUND" if bound else "UNBOUND"))
        seen.add(key)
    rows.sort(key=lambda r:(-float(r.get("score") or 0),str(r.get("case_gene")),
                            str(r.get("family")),str(r.get("target_concept"))))
    for i,r in enumerate(rows,1):
        r["route_id"]=f"V30R{i:04d}";r["rank"]=i
    gated,audit=v28_science.adjudicate_routes(rows)
    return {"schema":SCHEMA,"routes":gated,"audit":audit,
            "route_count":len(gated),"route_generation_uses_drug_names":False,
            "derived_before_drug_screen":True,
            "literature_discovered_route_count":sum(r.get("family")=="v30_literature_discovered" for r in gated)}


def merge_objective_routes(native_bundle,extra_routes,label="v31_extra"):
    """Merge V3.1 routes that already carry explicit rescue-objective provenance."""
    rows=[dict(x) for x in (native_bundle.get("routes") or [])]
    index={}
    for r in rows:
        key=(r.get("case_gene"),str(r.get("target_concept") or "").upper(),
             r.get("desired_action"),r.get("variable"))
        index[key]=r
    for x0 in extra_routes or []:
        x=dict(x0)
        key=(x.get("case_gene"),str(x.get("target_concept") or "").upper(),
             x.get("desired_action"),x.get("variable"))
        if not key[0] or not key[1]:continue
        if key in index:
            r=index[key]
            r["rescue_objective_ids"]=sorted(set((r.get("rescue_objective_ids") or [])+
                                                  (x.get("rescue_objective_ids") or [])))
            r["rescue_objective_types"]=sorted(set((r.get("rescue_objective_types") or [])+
                                                    (x.get("rescue_objective_types") or [])))
            r["evidence"]=(r.get("evidence") or [])+(x.get("evidence") or [])
            r["score"]=max(float(r.get("score") or 0),float(x.get("score") or 0))
            r.setdefault("merged_route_sources",[]).append(label)
            continue
        x.setdefault("state","conditional")
        x.setdefault("score",12.0)
        x.setdefault("targets",[])
        x.setdefault("merged_route_sources",[]).append(label)
        rows.append(x);index[key]=x
    rows.sort(key=lambda r:(-float(r.get("score") or 0),str(r.get("case_gene")),
                            str(r.get("family")),str(r.get("target_concept"))))
    for i,r in enumerate(rows,1):
        r["route_id"]=f"V31R{i:04d}";r["rank"]=i
    custom_gates={
      r["route_id"]:dict(r.get("scientific_gate") or {})
      for r in rows
      if str((r.get("scientific_gate") or {}).get("schema") or "").startswith("V31_")
    }
    gated,audit=v28_science.adjudicate_routes(rows)
    # V2 adjudication remains useful for legacy/new generic routes, but an
    # explicit V3.1 gate carries additional translational constraints that
    # must survive the merge (e.g. non-human rescue is not human validation).
    for r in gated:
        if r.get("route_id") in custom_gates:
            r["scientific_gate"]=custom_gates[r["route_id"]]
    return {"schema":SCHEMA,"routes":gated,"audit":audit,
            "route_count":len(gated),"route_generation_uses_drug_names":False,
            "derived_before_drug_screen":True,
            "v31_extra_route_count":sum(label in (r.get("merged_route_sources") or [])
                                        for r in gated)}
