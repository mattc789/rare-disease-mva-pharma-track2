"""Signed causal network expansion from SIGNOR for V3.1.

Generates multi-hop upstream regulator and downstream compensation routes.
Multi-hop sign propagation is mechanistic hypothesis generation, not validated
disease rescue; every >1-hop route is explicitly held to source verification.
"""
from __future__ import annotations
from collections import defaultdict,deque
from v2810 import routes as v28_routes

SCHEMA="RDMVA_SIGNED_NETWORK_EXPANSION_V31_1"

def _graphs():
    fwd=defaultdict(list);rev=defaultdict(list)
    for e in v28_routes.signor_edges():
        if not e.get("direct"):continue
        a=e.get("source_acc");b=e.get("target_acc")
        if not a or not b or e.get("sign") not in (-1,1):continue
        fwd[a].append(e);rev[b].append(e)
    return fwd,rev

def _walk(start,adj,direction,max_hops):
    out=[];q=deque([(start,1,[],{start})])
    while q:
        node,sign,path,seen=q.popleft()
        if len(path)>=max_hops:continue
        for e in adj.get(node,[]):
            nxt=e["source_acc"] if direction=="reverse" else e["target_acc"]
            if nxt in seen:continue
            ns=sign*e["sign"];np=path+[e]
            out.append((nxt,ns,np))
            q.append((nxt,ns,np,seen|{nxt}))
    return out

def build(private,objectives,cache_dir,max_hops=2):
    fwd,rev=_graphs();rows=[]
    objectives_by_gene=defaultdict(list)
    for o in objectives.get("objectives",[]):
        objectives_by_gene[o.get("gene")].append(o)
    for gene in sorted(private.get("genes") or {}):
        try:acc=v28_routes.accession_for_gene(gene,cache_dir)
        except Exception:acc=None
        if not acc:continue
        # Upstream: choose regulator action that tends to increase/support primary gene.
        for target_acc,sgn,path in _walk(acc,rev,"reverse",max_hops):
            e0=path[-1]  # farthest source is source of final reverse edge
            target_gene=e0.get("source_gene")
            if not target_gene:continue
            hops=len(path);desired="increase" if sgn>0 else "decrease"
            for o in objectives_by_gene.get(gene,[]):
                # A generic signed activity path can support an activity-rescue
                # hypothesis only. It cannot prove abundance, folding, splicing,
                # NMD, translation, complex assembly or downstream phenotype.
                if o.get("objective_type")!="increase_residual_activity":continue
                rows.append({
                  "family":"v31_signed_upstream_multihop",
                  "case_gene":gene,"target_concept":target_gene,
                  "targets":[{"gene":target_gene,"accession":target_acc}],
                  "desired_action":desired,"variable":"activity",
                  "mechanism_concept":f"{desired} {target_gene} through signed {hops}-hop path to support {gene}",
                  "score":max(6.0,13.0-2.0*(hops-1)),
                  "state":"conditional","source":"SIGNOR_signed_path",
                  "network_hops":hops,"network_path_sign":sgn,
                  "network_path":[{"source_gene":e.get("source_gene"),
                                   "target_gene":e.get("target_gene"),
                                   "sign":e.get("sign"),"pmid":e.get("pmid"),
                                   "signor_id":e.get("signor_id")} for e in reversed(path)],
                  "rescue_objective_id":o.get("objective_id"),
                  "rescue_objective_ids":[o.get("objective_id")],
                  "rescue_objective_type":o.get("objective_type"),
                  "rescue_objective_types":[o.get("objective_type")],
                  "case_gate_reasons":(["multi_hop_signed_path_requires_context_validation"] if hops>1 else []),
                })
        # Downstream: mimic the direction normally imposed by the primary gene.
        for target_acc,sgn,path in _walk(acc,fwd,"forward",max_hops):
            elast=path[-1];target_gene=elast.get("target_gene")
            if not target_gene:continue
            hops=len(path);desired="increase" if sgn>0 else "decrease"
            for o in objectives_by_gene.get(gene,[]):
                if o.get("objective_type") not in (
                    "compensate_downstream_deficit","bypass_or_compensate_downstream"):continue
                rows.append({
                  "family":"v31_signed_downstream_multihop",
                  "case_gene":gene,"target_concept":target_gene,
                  "targets":[{"gene":target_gene,"accession":target_acc}],
                  "desired_action":desired,"variable":"activity",
                  "mechanism_concept":f"{desired} {target_gene} to mimic {gene} downstream output over {hops} hop(s)",
                  "score":max(5.0,11.0-2.0*(hops-1)),
                  "state":"conditional","source":"SIGNOR_signed_path",
                  "network_hops":hops,"network_path_sign":sgn,
                  "network_path":[{"source_gene":e.get("source_gene"),
                                   "target_gene":e.get("target_gene"),
                                   "sign":e.get("sign"),"pmid":e.get("pmid"),
                                   "signor_id":e.get("signor_id")} for e in path],
                  "rescue_objective_id":o.get("objective_id"),
                  "rescue_objective_ids":[o.get("objective_id")],
                  "rescue_objective_type":o.get("objective_type"),
                  "rescue_objective_types":[o.get("objective_type")],
                  "case_gate_reasons":["downstream_compensation_requires_disease_endpoint"]+
                    (["multi_hop_signed_path_requires_context_validation"] if hops>1 else []),
                })
    # dedupe by exact pharmacologic action + objective
    d={}
    for r in rows:
        k=(r["case_gene"],r["target_concept"],r["desired_action"],
           r["variable"],r["rescue_objective_id"])
        if k not in d or r["score"]>d[k]["score"]:d[k]=r
    out=sorted(d.values(),key=lambda x:(-x["score"],x["network_hops"],
                                        x["case_gene"],x["target_concept"]))
    for i,r in enumerate(out,1):
        r["route_id"]=f"V31S{i:05d}";r["rank"]=i
        hard=False
        issues=list(r.get("case_gate_reasons") or [])
        r["scientific_gate"]={
          "schema":"V31_SIGNED_NETWORK_GATE_1",
          "primary_gene_target_identity":False,
          "allele_specific_rescue_verified":False,
          "eligible_for_exploratory_screen":not hard,
          "functional_rescue_validated":False,
          "causal_direction_validated_for_disease":False,
          "requires_independent_endpoint":True,
          "hard_hold":hard,"issues":issues}
    return {"schema":SCHEMA,"routes":out,"route_count":len(out),
            "max_hops":max_hops,"multi_hop_is_validated_rescue":False,
            "clinical_recommendation":False}
