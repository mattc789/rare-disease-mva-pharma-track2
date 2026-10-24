"""Genome-wide causal disease-state graph for Rare Disease MVA Pharma V3.1.

Primary hypotheses retain the only patient-specific directional rescue claims.
Genome-wide alternatives expand network/pathway discovery but are explicitly
non-diagnostic and cannot by themselves create a treatment hypothesis.
"""
from __future__ import annotations
from collections import defaultdict
import hashlib,json

SCHEMA="RDMVA_DISEASE_STATE_GRAPH_V31_1"

def _hid(x):
    return hashlib.sha256(json.dumps(x,sort_keys=True,ensure_ascii=True).encode()).hexdigest()[:16]

def select_genes(private,secondary_max=12,relative_margin=8.0,min_score=18.0,
                 min_phenotype_score=1.5,min_damaging=1):
    primary={str(x).upper() for x in (private.get("genes") or {})}
    ranked=list(private.get("genome_wide_candidate_genes") or [])
    top=float(ranked[0].get("score") or 0) if ranked else 0.0
    rows=[]
    for rank,x in enumerate(ranked,1):
        gene=str(x.get("gene") or "").upper()
        if not gene:continue
        role="PRIMARY" if gene in primary else "SECONDARY"
        if role=="SECONDARY":
            if len([r for r in rows if r["role"]=="SECONDARY"])>=secondary_max:continue
            if float(x.get("score") or 0)<max(min_score,top-relative_margin):continue
            if float(x.get("phenotype_score") or 0)<min_phenotype_score:continue
            if int(x.get("damaging_pass_count") or 0)<min_damaging:continue
            if x.get("mapping_artifact_review"):continue
        rows.append({
          "gene":gene,"role":role,"genome_rank":rank,
          "rank_score_not_diagnostic":x.get("score"),
          "phenotype_score":x.get("phenotype_score"),
          "damaging_variant_count":x.get("damaging_pass_count"),
          "best_variant_severity":x.get("best_severity"),
          "clinvar_weight":x.get("clinvar_weight"),
          "homozygous_damaging_count":x.get("homozygous_damaging_count"),
          "heterozygous_damaging_count":x.get("heterozygous_damaging_count"),
          "phase_state":(x.get("allelic_phase_evidence") or {}).get("state"),
          "directional_pharmacology_allowed":role=="PRIMARY",
          "diagnosis_asserted":False,
        })
    # preserve primary genes even if absent from genome-wide ledger
    known={x["gene"] for x in rows}
    for gene in sorted(primary-known):
        rows.append({"gene":gene,"role":"PRIMARY","genome_rank":None,
                     "directional_pharmacology_allowed":True,
                     "diagnosis_asserted":False})
    rows.sort(key=lambda x:(0 if x["role"]=="PRIMARY" else 1,
                            x["genome_rank"] or 99999,x["gene"]))
    return rows

def build(private,patient_state,mechanisms,objectives,signor_edges=None):
    genes=select_genes(private)
    by_gene={x["gene"]:x for x in genes}
    nodes=[];edges=[]
    for g in genes:
        nodes.append({"id":"GENE:"+g["gene"],"type":"gene",**g})
    for o in objectives.get("objectives",[]):
        oid="OBJECTIVE:"+o["objective_id"]
        nodes.append({"id":oid,"type":"rescue_objective",
                      "gene":o.get("gene"),"objective_type":o.get("objective_type"),
                      "priority_tier":o.get("priority_tier"),
                      "desired_action":o.get("desired_action"),
                      "variable":o.get("variable")})
        edges.append({"source":"GENE:"+o["gene"],"target":oid,
                      "type":"has_rescue_objective",
                      "directional":True})
    for m in mechanisms.get("hypotheses",[]):
        if m.get("status")=="NOT_INFERRED_WITHOUT_DIRECT_EVIDENCE":continue
        mid="MECH:"+m["mechanism_id"]
        nodes.append({"id":mid,"type":"mechanism","gene":m.get("gene"),
                      "mechanism_type":m.get("mechanism_type"),
                      "defect_direction":m.get("defect_direction"),
                      "status":m.get("status")})
        edges.append({"source":"GENE:"+m["gene"],"target":mid,
                      "type":"has_mechanism_hypothesis"})
    for e in signor_edges or []:
        a=str(e.get("source_gene") or "").upper()
        b=str(e.get("target_gene") or "").upper()
        if a in by_gene and b in by_gene:
            edges.append({"source":"GENE:"+a,"target":"GENE:"+b,
                          "type":"signed_interaction","sign":e.get("sign"),
                          "direct":bool(e.get("direct")),
                          "source_id":e.get("signor_id"),"pmid":e.get("pmid")})
    primary=[x["gene"] for x in genes if x["role"]=="PRIMARY"]
    secondary=[x["gene"] for x in genes if x["role"]=="SECONDARY"]
    return {
      "schema":SCHEMA,
      "graph_id":"DSG:"+_hid({"primary":primary,"secondary":secondary}),
      "primary_genes":primary,"secondary_genes":secondary,
      "network_seed_genes":primary+secondary,
      "directional_rescue_genes":primary,
      "nodes":nodes,"edges":edges,
      "secondary_genes_are_diagnostic_claims":False,
      "secondary_genes_can_expand_network_search":True,
      "secondary_genes_can_create_directional_drug_claims":False,
      "clinical_recommendation":False,
    }
