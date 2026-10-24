"""Mechanism Engine V3.0: variant class -> explicit falsifiable mechanism hypotheses."""
from __future__ import annotations
import hashlib,json

SCHEMA="RDMVA_MECHANISM_ENGINE_V30_1"

def _id(payload):
    return "MECH:"+hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()[:16]

def _add(out,gene,kind,defect,status,evidence,requires):
    payload={"gene":gene,"mechanism_type":kind,"defect_direction":defect}
    out.append({
      "mechanism_id":_id(payload),"gene":gene,"mechanism_type":kind,
      "defect_direction":defect,"status":status,
      "evidence_basis":sorted(set(evidence)),
      "requirements_to_confirm":requires,
      "clinical_mechanism_confirmed":False,
      "falsifiable":True,
    })

def build(patient_state):
    hypotheses=[]
    for g in patient_state.get("genes",[]):
        gene=g["gene"];classes=set(g.get("effect_classes") or [])
        refs=["variant_class:"+x for x in sorted(classes)]
        if classes & {"premature_stop","frameshift"}:
            _add(hypotheses,gene,"loss_of_function_or_truncated_product","decrease_function",
                 "PLAUSIBLE_FROM_VARIANT_CLASS",refs,
                 ["RNA abundance/NMD assay","full-length protein abundance","gene-specific functional assay"])
            _add(hypotheses,gene,"downstream_function_deficit","decrease_pathway_output",
                 "PLAUSIBLE_SECONDARY_MECHANISM",refs,
                 ["disease-relevant cellular phenotype","signed downstream pathway measurements"])
        if "premature_stop" in classes:
            _add(hypotheses,gene,"premature_stop_with_possible_nmd","decrease_productive_transcript",
                 "PLAUSIBLE_FROM_VARIANT_CLASS",refs,
                 ["allele-specific transcript abundance","NMD dependence","full-length translation assay"])
        if "splice_altering" in classes:
            _add(hypotheses,gene,"aberrant_splicing","decrease_productive_transcript",
                 "PLAUSIBLE_FROM_VARIANT_CLASS",refs,
                 ["RNA isoform quantification","junction validation","protein/function rescue assay"])
        if classes & {"missense","inframe_indel"}:
            _add(hypotheses,gene,"protein_stability_or_folding_defect","decrease_function",
                 "PLAUSIBLE_NOT_PROVEN",refs,
                 ["protein abundance/half-life","folding/localization assay","variant-specific functional assay"])
            _add(hypotheses,gene,"protein_activity_or_interaction_defect","decrease_function",
                 "PLAUSIBLE_NOT_PROVEN",refs,
                 ["biochemical activity","interaction/complex assembly","cellular rescue assay"])
        if not classes or classes=={"uncertain_effect"}:
            _add(hypotheses,gene,"mechanism_unresolved","unknown",
                 "MEASUREMENT_REQUIRED",refs or ["variant_effect_unknown"],
                 ["orthogonal variant interpretation","RNA/protein/function measurements"])
        # Never infer gain-of-function merely from a coding variant.
        hypotheses.append({
          "mechanism_id":_id({"gene":gene,"mechanism_type":"gain_of_function"}),
          "gene":gene,"mechanism_type":"gain_of_function",
          "defect_direction":"increase_function",
          "status":"NOT_INFERRED_WITHOUT_DIRECT_EVIDENCE",
          "evidence_basis":[],"requirements_to_confirm":["direct variant-specific gain-of-function evidence"],
          "clinical_mechanism_confirmed":False,"falsifiable":True,
        })
    return {
      "schema":SCHEMA,"patient_state_id":patient_state.get("state_id"),
      "hypotheses":hypotheses,
      "confirmed_mechanisms":sum(x["clinical_mechanism_confirmed"] for x in hypotheses),
      "gain_of_function_not_assumed":True,
      "clinical_recommendation":False,
    }
