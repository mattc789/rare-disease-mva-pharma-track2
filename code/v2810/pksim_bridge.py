"""Read actual PK-Sim CLI wide-format CSV, with strict units and source roles.

The public Midazolam OSP reference control can be audited at total-plasma level,
but not used as patient unbound exposure without measured f_u and other data.
"""
from __future__ import annotations
from pathlib import Path
import csv,json

SCHEMA="V289_PK_SIM_EXPORT_BRIDGE_1"

def available_series(csv_file):
    with Path(csv_file).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.reader(f)
        fields=next(reader,[])
    return {
      "schema":SCHEMA,"time_column":"Time [min]" if "Time [min]" in fields else None,
      "concentration_columns":[c for c in fields
          if ("concentration" in c.lower() or "plasma" in c.lower())
          and ("umol/l" in c.lower() or "µmol/l" in c.lower() or
               "nmol/l" in c.lower() or "\u00b5mol/l" in c.lower())],
      "unit_conversion_requires_explicit_column_selection":True,
      "unbound_pk_not_inferred_from_total_plasma":True,
    }

def _unit(col):
    low=col.lower().replace("μ","µ")
    if "nmol/l" in low:return 1.0
    if "µmol/l" in low or "umol/l" in low:return 1000.0
    if "mmol/l" in low:return 1e6
    raise ValueError("V289_PKSIM_UNRECOGNIZED_UNIT_NEVER_GUESS")

def extract_curve(csv_file,series_column,*,selected_individual_id="0",
       unbound_fraction=None,unbound_fraction_source_reference=None):
    """Converts PKSim time in min to hours and µmol/L to nM.

    Without externally sourced unbound fraction, keeps total_plasma_nM;
    never invents patient free concentration or a validated exposure model.
    """
    from .quant_pk import fnum
    factor=_unit(series_column)
    fu=fnum(unbound_fraction) if unbound_fraction is not None else None
    if fu is not None and fu>1:raise ValueError("V289_PK_FU_GT_ONE")
    if fu is not None and not str(unbound_fraction_source_reference or "").strip():
        raise ValueError("V289_PK_UNBOUND_FRACTION_SOURCE_MISSING")
    if unbound_fraction is not None and fu is None:raise ValueError("V289_PK_FU_INVALID")
    with Path(csv_file).open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f)
        if series_column not in (rd.fieldnames or []) or "Time [min]" not in (rd.fieldnames or []):
            raise ValueError("V289_PKSIM_UNVERIFIED_TIME_OR_SERIES_COLUMN")
        if "IndividualId" not in (rd.fieldnames or []):
            raise ValueError("V289_PKSIM_INDIVIDUAL_ID_MISSING")
        out=[]
        for r in rd:
            if str(r["IndividualId"])!=str(selected_individual_id):continue
            t=fnum(r["Time [min]"],positive=False)
            c=fnum(r.get(series_column),positive=False)
            if t is None or c is None or t<0 or c<0:continue
            row={"time_h":round(t/60,8),"total_plasma_nM":c*factor}
            if fu is not None:row["calculated_free_plasma_nM"]=round(c*factor*fu,8)
            out.append(row)
    if not out:raise ValueError("V289_PKSIM_SERIES_HAS_NO_USABLE_OBSERVATIONS")
    times=[p["time_h"] for p in out]
    if len(times)!=len(set(times)):raise ValueError("V289_PKSIM_DUPLICATE_TIMEPOINTS")
    return {"schema":SCHEMA,"series_column":series_column,
       "individual_id":str(selected_individual_id),
       "total_plasma_concentration_available":True,
       "unbound_fraction_source_qualification_required":fu is None,
       "free_plasma_calculation_if_curated":fu is not None,
       "unbound_fraction_source_reference":unbound_fraction_source_reference if fu is not None else None,
       "prediction_is_not_experimentally_validated_by_import":True,
       "patient_safety_or_treatment_inference_allowed":False,
       "observations":out}

def audit_control(run_folder,output):
    """Checks an existing exported source PBPK control; never executes dosing."""
    folder=Path(run_folder)
    files=sorted(folder.rglob("*-Results.csv")) if folder.is_dir() else []
    rows=[]
    for f in files:
        try:
            src=available_series(f)
            rows.append({"simulation_file":f.name,
              "series_columns_count":len(src["concentration_columns"]),
              "valid_time_column":bool(src["time_column"]),
              "quantitative_profiles_present":bool(src["concentration_columns"])})
        except (OSError,ValueError,UnicodeError):
            rows.append({"simulation_file":f.name,"error":"PKSIM_SERIES_HEADER_NOT_READABLE"})
    summary={"schema":SCHEMA,
      "model_kind":"PUBLIC_MIDAZOLAM_ENGINEERING_CONTROL_ONLY",
      "simulations_exported":len(files),
      "simulations_with_parsable_concentrations":sum(bool(x.get("quantitative_profiles_present")) for x in rows),
      "independent_pk_population_qualification_from_just_csv":False,
      "clinical_patient_benefit_or_dose_inferred":False,
      "control_is_not_a_candidate_for_rare_disease":True}
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps({"summary":summary,"simulations":rows},indent=2)+"\n",encoding="utf-8")
    return summary
