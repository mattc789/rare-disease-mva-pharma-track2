"""Rescue Objective Generator V3.0.

Transforms mechanistic hypotheses into directional research objectives BEFORE
drug search. Objectives are hypotheses and never prescriptions.
"""
from __future__ import annotations
import hashlib,json

SCHEMA="RDMVA_RESCUE_OBJECTIVES_V30_1"

def _oid(gene,kind,action,variable):
    blob=json.dumps([gene,kind,action,variable],separators=(",",":"))
    return "OBJ:"+hashlib.sha256(blob.encode()).hexdigest()[:16]

def _axis(kind):
    return {
      "restore_functional_gene_product":"primary_function",
      "increase_residual_activity":"primary_function",
      "restore_complex_or_interaction_function":"primary_function",
      "restore_full_length_translation":"translation",
      "restore_productive_splicing":"splicing",
      "stabilize_residual_protein":"proteostasis",
      "reduce_causal_degradation":"proteostasis",
      "compensate_downstream_deficit":"downstream_compensation",
      "bypass_or_compensate_downstream":"downstream_compensation",
      "reduce_nmd_conditionally":"rna_surveillance",
      "measurement_first":"diagnostic_research",
    }.get(kind,kind)

def _obj(out,gene,mech,kind,action,variable,layer,tier,conditional,requirements):
    out.append({
      "objective_id":_oid(gene,kind,action,variable),
      "gene":gene,"mechanism_id":mech,"objective_type":kind,
      "objective_axis":_axis(kind),
      "desired_action":action,"variable":variable,"target_layer":layer,
      "priority_tier":tier,"conditional":bool(conditional),
      "requirements_before_promotion":requirements,
      "drug_search_eligible":kind!="measurement_first",
      "clinical_recommendation":False,
    })

def build(mechanisms):
    out=[]
    for m in mechanisms.get("hypotheses",[]):
        if m.get("status")=="NOT_INFERRED_WITHOUT_DIRECT_EVIDENCE":continue
        gene=m["gene"];mid=m["mechanism_id"];kind=m["mechanism_type"]
        if kind=="loss_of_function_or_truncated_product":
            _obj(out,gene,mid,"restore_functional_gene_product","increase","functional_abundance",
                 "primary_gene","O1",True,
                 ["confirm residual/productive allele potential","demonstrate functional rescue"])
            _obj(out,gene,mid,"compensate_downstream_deficit","restore","disease_relevant_pathway_output",
                 "downstream_network","O2",True,
                 ["signed causal path","disease-relevant functional endpoint"])
        elif kind=="premature_stop_with_possible_nmd":
            _obj(out,gene,mid,"restore_full_length_translation","increase","full_length_translation",
                 "RNA_translation","O1",True,
                 ["confirm premature-stop context","show full-length functional product"])
            _obj(out,gene,mid,"reduce_nmd_conditionally","decrease","nmd_of_implicated_transcript",
                 "RNA_surveillance","O3",True,
                 ["prove transcript is NMD-sensitive","exclude harmful truncated product"])
        elif kind=="aberrant_splicing":
            _obj(out,gene,mid,"restore_productive_splicing","increase","productive_transcript",
                 "RNA_splicing","O1",True,
                 ["quantify abnormal isoform","show corrected transcript/protein/function"])
        elif kind=="protein_stability_or_folding_defect":
            _obj(out,gene,mid,"stabilize_residual_protein","increase","folding_stability",
                 "proteostasis","O1",True,
                 ["confirm residual protein","show stability plus functional rescue"])
            _obj(out,gene,mid,"reduce_causal_degradation","decrease","degradation_pressure",
                 "proteostasis","O2",True,
                 ["identify causal degradation machinery","exclude global proteostasis toxicity"])
        elif kind=="protein_activity_or_interaction_defect":
            _obj(out,gene,mid,"increase_residual_activity","increase","protein_specific_activity",
                 "primary_gene","O1",True,
                 ["measure residual activity","show direct functional rescue"])
            _obj(out,gene,mid,"restore_complex_or_interaction_function","restore","complex_function",
                 "protein_network","O2",True,
                 ["identify disrupted interaction","rescue disease-relevant endpoint"])
        elif kind=="downstream_function_deficit":
            _obj(out,gene,mid,"bypass_or_compensate_downstream","restore","disease_relevant_pathway_output",
                 "downstream_network","O2",True,
                 ["signed causal path","same-drug functional rescue"])
        elif kind=="mechanism_unresolved":
            _obj(out,gene,mid,"measurement_first","measure","causal_molecular_defect",
                 "diagnostic_research","O0",False,
                 ["resolve mechanism before directional pharmacology"])
    # de-duplicate exact objectives produced by overlapping mechanisms
    dedup={}
    for x in out:
        dedup.setdefault((x["gene"],x["objective_type"],x["desired_action"],x["variable"]),x)
    rows=list(dedup.values())
    rows.sort(key=lambda x:(x["priority_tier"],x["gene"],x["objective_type"]))
    return {
      "schema":SCHEMA,"mechanism_schema":mechanisms.get("schema"),
      "objectives":rows,
      "drug_search_objectives":sum(x["drug_search_eligible"] for x in rows),
      "objectives_are_patient_specific_hypotheses_not_treatments":True,
    }
