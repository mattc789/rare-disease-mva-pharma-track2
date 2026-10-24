"""Reproducible quantitative PK / PBPK evidence qualification.

Supports source-qualified PK-Sim exports with independent observed PK, plus
a reference-tested analytic IV bolus one-compartment model ONLY when complete
literature-qualified parameters and validation data are provided.

No default drug dosing, pediatric extrapolation or efficacy probability.
"""
from __future__ import annotations
import csv,hashlib,json,math,statistics
from collections import Counter
from pathlib import Path

SCHEMA="V289_QUANT_PK_1"
MAX_POINTS=2000
def fnum(x,*,positive=True):
    try:v=float(x)
    except (ValueError,TypeError,OverflowError):return None
    if not math.isfinite(v) or (positive and v<=0):return None
    return v

def _sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        while True:
            b=f.read(1024*1024)
            if not b:break
            h.update(b)
    return h.hexdigest()

def analytical_iv_profile(dose_mg,cl_l_h,vd_l,mw_g_mol,fu,times_h):
    """Closed-form 1C IV bolus: C_free[nM] = (dose/V)*1e6/MW*fu*exp[-(CL/V)t]."""
    vals=[fnum(x) for x in (dose_mg,cl_l_h,vd_l,mw_g_mol,fu)]
    if any(v is None for v in vals):raise ValueError("V289_PK_NONPOSITIVE_INPUT")
    dose,cl,vd,mw,unbound=vals
    if unbound>1:raise ValueError("V289_PK_UNBOUND_FRACTION_GT_1")
    if len(times_h)>MAX_POINTS:raise ValueError("V289_PK_TIMEGRID_TOO_LONG")
    k=cl/vd
    c0=(dose/vd)*1_000_000/mw*unbound
    result=[]
    for t in times_h:
        tv=fnum(t,positive=False)
        if tv is None or tv<0:raise ValueError("V289_PK_INVALID_TIME")
        result.append({"time_h":tv,"predicted_free_plasma_nM":c0*math.exp(-k*tv)})
    return result

def _load_curves(file):
    rows=[]
    with Path(file).open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f)
        fields=set(rd.fieldnames or [])
        if not {"time_h","unbound_concentration_nM"}<=fields:
            raise ValueError("V289_PK_CURVE_HEADER_MISMATCH")
        for row in rd:
            t=fnum(row["time_h"],positive=False)
            conc=fnum(row["unbound_concentration_nM"])
            if t is None or t<0 or conc is None:continue
            rows.append((t,conc))
    if not rows or len(rows)>MAX_POINTS:raise ValueError("V289_PK_CURVE_EMPTY_OR_TOO_LONG")
    if len({t for t,c in rows})!=len(rows):raise ValueError("V289_PK_DUPLICATE_TIMES")
    return sorted(rows)

def _interpolate(curve,t):
    if not curve or t<curve[0][0] or t>curve[-1][0]:return None
    if t==curve[0][0]:return curve[0][1]
    for i in range(1,len(curve)):
        t0,y0=curve[i-1];t1,y1=curve[i]
        if t0<=t<=t1:
            if t==t1:return y1
            f=(t-t0)/(t1-t0)
            return y0+(y1-y0)*f
    return None

def validate_against_observed(predicted,observed,require_min=3):
    """External qualification: geometric fold error is a descriptive fit statistic."""
    ratios=[];matched=[]
    for t,y in observed:
        guess=_interpolate(predicted,t)
        if guess is None or guess<=0 or y<=0:continue
        ratios.append(guess/y);matched.append((t,guess,y))
    if len(matched)<require_min:
        return {"state":"INSUFFICIENT_INDEPENDENT_OBSERVATIONS",
         "matched_observations":len(matched),"independent_clinical_model_validated":False}
    logs=[abs(math.log(r)) for r in ratios]
    aafe=math.exp(sum(logs)/len(logs))
    bias=math.exp(sum(math.log(r) for r in ratios)/len(ratios))
    all_twofold=all(0.5<=r<=2 for r in ratios)
    return {"state":"EXTERNAL_PROFILE_CHECK_COMPLETED",
       "matched_observations":len(matched),
       "absolute_average_fold_error":round(aafe,4),
       "geometric_mean_prediction_observation_ratio":round(bias,4),
       "all_evaluation_points_within_twofold":all_twofold,
       "predeclared_twofold_screen_met":bool(aafe<=2 and all_twofold),
       "clinical_individual_safety_or_efficacy_validated":False}

def _trial_has_evidence(row):
    return bool(
        row.get("source_reference") and row.get("independent_validation_reference")
        and row.get("compound") and row.get("analyte")=="unbound_parent_drug"
        and row.get("population") and row.get("dose_regimen_reference")
        and row.get("formulation") and row.get("route")
        and row.get("source_checked") is True
        and row.get("independent_evaluation_cohort") is True
        and row.get("study_reports_dose_and_sampling") is True)

def assess(row):
    base={"compound":str(row.get("compound") or ""),"model_kind":row.get("model_kind"),
      "state":"INSUFFICIENT_MODEL_PROVENANCE","independent_pk_validation":False,
      "model_is_patient_specific":False,"licensed_for_patient_dose_decisions":False,
      "predicted_therapeutic_benefit":None,
      "observed_validation_count":0}
    if not _trial_has_evidence(row):return base
    source=Path(str(row.get("observations_csv") or ""))
    if not source.is_file():return {**base,"state":"QUALIFICATION_OBSERVATIONS_MISSING"}
    try:observations=_load_curves(source)
    except (ValueError,OSError):return {**base,"state":"INVALID_QUALIFICATION_OBSERVATIONS"}
    kind=row.get("model_kind")
    if kind=="OSP_PK_SIM_PREVALIDATED":
        f=Path(str(row.get("simulation_export_csv") or ""))
        p=Path(str(row.get("pksim_project_file") or ""))
        report=Path(str(row.get("model_evaluation_report") or ""))
        if not f.is_file() or not p.is_file() or not report.is_file():
            return {**base,"state":"PK_SIM_MODEL_EXPORT_PROJECT_OR_EVALUATION_REPORT_MISSING"}
        if _sha(p)!=row.get("project_sha256") or _sha(report)!=row.get("report_sha256"):
            return {**base,"state":"PK_SIM_MODEL_PROVENANCE_HASH_MISMATCH"}
        try:predicted=_load_curves(f)
        except (ValueError,OSError):return {**base,"state":"INVALID_PBPK_SIMULATION_EXPORT"}
        tool="Open Systems Pharmacology PK-Sim"
    elif kind=="LITERATURE_QUALIFIED_ONE_COMPARTMENT_IV":
        if row.get("route") not in ("intravenous_bolus", "intravenous"):
            return {**base,"state":"ONE_COMPARTMENT_IV_ROUTE_UNSUPPORTED"}
        fields=["dose_mg","clearance_L_h","volume_L","molecular_weight_g_mol","unbound_fraction"]
        if any(fnum(row.get(k)) is None for k in fields):
            return {**base,"state":"ONE_COMPARTMENT_PARAMETERS_MISSING"}
        if not row.get("parameter_source_reference"):
            return {**base,"state":"ONE_COMPARTMENT_PARAMETER_SOURCE_MISSING"}
        try:
            grid=sorted({0.0,*[t for t,y in observations]})
            calculated=analytical_iv_profile(row["dose_mg"],row["clearance_L_h"],row["volume_L"],
               row["molecular_weight_g_mol"],row["unbound_fraction"],grid)
            predicted=[(x["time_h"],x["predicted_free_plasma_nM"]) for x in calculated]
        except (ValueError,OverflowError):
            return {**base,"state":"ONE_COMPARTMENT_SIMULATION_FAILED"}
        tool="closed_form_1_compartment_IV_bolus"
    else:return {**base,"state":"UNSUPPORTED_PK_MODEL_TYPE"}
    validation=validate_against_observed(predicted,observations)
    # The criteria are a transparent ENGINEERING screen, not regulatory model qualification.
    qualified=(validation.get("predeclared_twofold_screen_met") is True
        and row.get("model_applicability_for_this_population") is True
        and row.get("clinical_equations_qualification_documented") is True)
    state=("SOURCE_QUALIFIED_POPULATION_PK_CONTEXT" if qualified
           else "MODEL_EVALUATED_NOT_QUALIFIED_FOR_TARGET_POPULATION")
    # Exposure vs target potency is a research-context filter only; no treatment recommendation.
    v={**base,"state":state,
       "simulation_engine":tool,
       "observed_validation_count":validation["matched_observations"],
       "quantitative_validation":validation,
       "independent_pk_validation":qualified,
       "model_is_patient_specific":False,
       "qualification_is_not_target_engagement_proof":True,
       "extrapolation_to_different_age_or_tissue_disallowed":True,
       "source_reference":row["source_reference"],
       "independent_validation_reference":row["independent_validation_reference"],
       "population":row["population"],
       "dose_regimen_reference":row["dose_regimen_reference"],
       "formulation":row["formulation"],
       "simulated_timepoints":len(predicted),
       "concentration_unit":"unbound_plasma_nM",
    }
    if qualified:
        v["estimated_simulated_max_free_plasma_nM"]=round(max(y for t,y in predicted),3)
        # No specific clinical dosages are emitted in reports.
    return v

def run(outdir,manifest_path=None,approved_drug_names=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    approved=set(str(x).upper() for x in (approved_drug_names or []))
    if not manifest_path or not Path(manifest_path).is_file():
        summary={"schema":SCHEMA,"state":"NO_QUALIFIED_QUANTITATIVE_DRUG_MODELS_INSTALLED",
          "models_evaluated":0,"models_source_qualified":0,
          "patients_with_validated_pk_predictions":0,
          "no_default_pk_parameters_invented":True,
          "validation_data_required":True}
        (out/"QUANT_PK_STATUS.json").write_text(json.dumps(summary,indent=2)+"\n")
        return summary,[]
    data=json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    if data.get("schema")!=SCHEMA:
        raise ValueError("V289_QUANT_PK_MANIFEST_SCHEMA_MISMATCH")
    items=[]
    for row in data.get("models",[]):
        if approved and str(row.get("compound") or "").upper() not in approved:continue
        items.append(assess(row))
    summary={"schema":SCHEMA,"state":"QUANTITATIVE_MODELS_AUDITED",
       "models_evaluated":len(items),
       "models_source_qualified":sum(x["independent_pk_validation"] for x in items),
       "qualification_statuses":dict(Counter(x["state"] for x in items)),
       "clinical_efficacy_probability_calculated":False,
       "models_are_valid_for_same_population_only":True,
       "no_patient_specific_treatment_recommendations":True}
    (out/"QUANT_PK_STATUS.json").write_text(json.dumps(summary,indent=2)+"\n")
    (out/"QUANT_PK_MODELS_AUDIT.json").write_text(json.dumps(items,indent=2)+"\n")
    return summary,items
