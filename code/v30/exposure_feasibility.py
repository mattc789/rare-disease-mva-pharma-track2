"""Exposure-feasibility gate for V3.1 research hypotheses.

Compares source-bound concentrations needed for target modulation with
observed human exposure from public regulatory labels. It never recommends a
dose and never treats total plasma concentration as free target exposure.
"""
from __future__ import annotations
from pathlib import Path
import html,json,os,re,urllib.request
from . import openai_reasoner

SCHEMA="RDMVA_EXPOSURE_FEASIBILITY_V31_1"
UA="RareDiseaseMVAPharma-V31/1.0"

def _html_text(pmcid):
    if not pmcid:return ""
    url=f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
    req=urllib.request.Request(url,headers={"User-Agent":UA})
    with urllib.request.urlopen(req,timeout=35) as r:
        raw=r.read(2_000_000).decode("utf-8","replace")
    raw=re.sub(r"(?is)<script.*?</script>|<style.*?</style>"," ",raw)
    text=html.unescape(re.sub(r"(?s)<[^>]+>"," ",raw))
    return " ".join(text.split())

def _relevant_context(text,candidate,target,target_aliases=None,max_chars=30000):
    """Prioritize numeric assay windows so PK-relevant evidence is not crowded out."""
    lo=text.lower();spans=[]
    aliases=[str(target or "").lower()]+[
      str(x or "").lower() for x in (target_aliases or []) if x]
    aliases=list(dict.fromkeys(x for x in aliases if x))
    # Highest priority: every explicit concentration window.
    conc_rx=re.compile(r"\b\d+(?:\.\d+)?\s*(?:µm|μm|um|nm|µg/ml|μg/ml|ug/ml|mg/l)\b",re.I)
    for m in conc_rx.finditer(text):
        a=max(0,m.start()-950);b=min(len(text),m.end()+1600)
        window=text[a:b]
        wlo=window.lower()
        priority=0
        if any(a in wlo for a in aliases):priority+=5
        if str(candidate or "").lower() in wlo:priority+=4
        if any(k in wlo for k in ("kinase","activity","binding","inhibit","assay","profil")):priority+=3
        spans.append((priority,window))
    # Next: target-centric windows, then candidate-centric windows.
    terms=[(a,5) for a in aliases]
    terms += [(str(candidate or "").lower(),2),("kinase profiling",4),
              ("directly bound",4),("activity",2),("binding",2)]
    for term,base in terms:
        start=0
        while term:
            i=lo.find(term,start)
            if i<0:break
            a=max(0,i-850);b=min(len(text),i+1500)
            spans.append((base,text[a:b]));start=i+max(1,len(term))
            if start>=len(text):break
    spans.sort(key=lambda x:-x[0])
    out=[];seen=set();n=0
    for _,s in spans:
        key=hash(s[:500])
        if key in seen:continue
        seen.add(key)
        if n+len(s)>max_chars:continue
        out.append(s);n+=len(s)
    return "\n---\n".join(out)

def _model():
    cfg=openai_reasoner.config()
    return (os.environ.get(cfg.get("environment_literature_model_override",
                                   "RDMVA_OPENAI_LITERATURE_MODEL"))
            or cfg.get("literature_extraction_model") or cfg["default_model"])

def extract_target_concentrations(candidate,target,source,outfile,budget_root,target_aliases=None):
    pmcid=(source or {}).get("pmcid")
    text=_html_text(pmcid)
    ctx=_relevant_context(text,candidate,target,target_aliases=target_aliases)
    schema={
      "type":"object","additionalProperties":False,
      "properties":{
        "records":{"type":"array","items":{
          "type":"object","additionalProperties":False,
          "properties":{
            "value":{"type":"number"},"unit":{"type":"string","enum":["nM","uM","ug/mL","mg/L","unknown"]},
            "assay_context":{"type":"string","enum":[
              "DIRECT_BINDING","BIOCHEMICAL_TARGET_ACTIVITY","TARGET_SPECIFIC_CELLULAR_ASSAY",
              "CELLULAR_PHENOTYPE","OTHER"]},
            "target_action":{"type":"string","enum":["increase","decrease","unknown"]},
            "effect_fraction_of_control":{"type":["number","null"]},
            "effect_description":{"type":"string","maxLength":500},
            "explicit_numeric_in_source":{"type":"boolean"}
          },
          "required":["value","unit","assay_context","target_action",
                      "effect_fraction_of_control","effect_description","explicit_numeric_in_source"]}},
        "lowest_direct_target_concentration_uM":{"type":["number","null"]},
        "highest_cellular_effect_concentration_uM":{"type":["number","null"]},
        "limitations":{"type":"array","items":{"type":"string"}}
      },
      "required":["records","lowest_direct_target_concentration_uM",
                  "highest_cellular_effect_concentration_uM","limitations"]
    }
    payload={"candidate":candidate,"target":target,
             "target_aliases":list(target_aliases or []),"pmcid":pmcid,
             "source_title":(source or {}).get("title"),"source_context":ctx,
             "patient_data_included":False}
    instructions=(
      "Extract ONLY explicit numeric concentrations from the supplied source context for the named candidate. "
      "For lowest_direct_target_concentration_uM use the lowest concentration at which direct binding, biochemical "
      "target activity, or a target-specific assay actually demonstrates modulation of the named target. "
      "IMPORTANT: when a kinase-profiling table reports percent-of-control for the named target at an explicit "
      "compound concentration, treat that as BIOCHEMICAL_TARGET_ACTIVITY. A percent-of-control below 100 means "
      "decreased activity; for example 45 means effect_fraction_of_control=0.45 and target_action=decrease. "
      "Do not use docking-only values, unrelated assays, concentrations of other compounds, or inferred IC50 values. "
      "Convert explicit nM to uM when needed. For cellular phenotype concentrations keep them separate. "
      "Do not recommend a dose."
    )
    body={"model":_model(),"input":json.dumps(payload,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],"max_output_tokens":2200,
          "text":{"format":{"type":"json_schema","name":"v31_exposure_target",
                           "strict":True,"schema":schema}}}
    rec=openai_reasoner.request_json(body,outfile,budget_root=budget_root)
    return rec.get("value") or {}

def extract_label_exposure(candidate,label,outfile,budget_root):
    label=label or {}
    sections=[]
    for key in ("clinical_pharmacology","pharmacokinetics","clinical_pharmacology_table",
                "pediatric_use","dosage_and_administration"):
        v=label.get(key)
        if isinstance(v,list):sections.extend(str(x) for x in v[:5])
        elif isinstance(v,str):sections.append(v)
    text="\n".join(sections)[:30000]
    schema={
      "type":"object","additionalProperties":False,
      "properties":{
        "cmax_records":{"type":"array","items":{
          "type":"object","additionalProperties":False,
          "properties":{
            "value_ug_per_mL":{"type":"number"},
            "route":{"type":"string"},"dose_description":{"type":"string"},
            "population":{"type":"string"},"time_description":{"type":"string"},
            "observed_label_value":{"type":"boolean"}
          },
          "required":["value_ug_per_mL","route","dose_description","population",
                      "time_description","observed_label_value"]}},
        "max_observed_cmax_ug_per_mL":{"type":["number","null"]},
        "protein_binding_low_pct":{"type":["number","null"]},
        "protein_binding_high_pct":{"type":["number","null"]},
        "binding_is_concentration_dependent":{"type":"boolean"},
        "limitations":{"type":"array","items":{"type":"string"}}
      },
      "required":["cmax_records","max_observed_cmax_ug_per_mL",
                  "protein_binding_low_pct","protein_binding_high_pct",
                  "binding_is_concentration_dependent","limitations"]
    }
    payload={"candidate":candidate,"public_label_sections":text,
             "instructional_context":"Observed regulatory-label PK only; no dose recommendation."}
    instructions=(
      "Extract observed or tabulated plasma Cmax/maximum plasma concentrations from the supplied public label. "
      "Do not invent values or recommend doses. Record route, the dose exactly as context text, and population. "
      "Extract the stated plasma protein binding percentage range. max_observed_cmax_ug_per_mL is the highest "
      "explicit observed/tabulated value in the supplied label, not a proposed dose or extrapolation."
    )
    body={"model":_model(),"input":json.dumps(payload,ensure_ascii=False),
          "instructions":instructions,"store":False,"tools":[],"max_output_tokens":2200,
          "text":{"format":{"type":"json_schema","name":"v31_exposure_label",
                           "strict":True,"schema":schema}}}
    rec=openai_reasoner.request_json(body,outfile,budget_root=budget_root)
    return rec.get("value") or {}

def molecular_weight_from_chembl(entity):
    ids=list((entity or {}).get("chembl_parent_ids") or [])
    if not ids:return None,None
    cid=ids[0]
    url=f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
    try:
        from v2810.chembl_sqlite_v31 import configured,local_json,PROVENANCE
        if configured() or PROVENANCE.is_file():
            j=local_json(url)
        else:
            req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=30) as r:j=json.loads(r.read(2_000_000))
        p=j.get("molecule_properties") or {}
        mw=float(p.get("full_mwt") or p.get("mw_freebase"))
        return mw,cid
    except Exception:return None,cid

def assess(candidate,entity,bridge,label,outdir,quantitative_override=None):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    target=(bridge or {}).get("target_symbol")
    source=(bridge or {}).get("drug_target_source") or {}
    if quantitative_override is None:
        project_root=Path(__file__).resolve().parents[2]
        quantitative_override=load_curated_record(project_root,candidate,target)
    if quantitative_override:
        target_exp=dict(quantitative_override)
        target_exp["source_mode"]="CURATED_SOURCE_AUDITED_NUMERIC_FACTS"
        target_exp["raw_source_text_redistributed"]=False
        (out/"TARGET_CONCENTRATIONS.json").write_text(
          json.dumps(target_exp,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    else:
        target_exp=extract_target_concentrations(
            candidate,target,source,out/"TARGET_CONCENTRATIONS.json",out)
    label_exp=extract_label_exposure(
        candidate,label,out/"LABEL_EXPOSURE.json",out)
    label_mw=molecular_weight_from_label(label)
    chembl_mw,cid=molecular_weight_from_chembl(entity)
    mw=label_mw if label_mw is not None else chembl_mw
    mw_source="REGULATORY_LABEL_FORMULATION" if label_mw is not None else "CHEMBL_PARENT"
    cmax=label_exp.get("max_observed_cmax_ug_per_mL")
    total_uM=(float(cmax)*1000.0/mw) if (cmax is not None and mw) else None
    bind_hi=label_exp.get("protein_binding_high_pct")
    bind_lo=label_exp.get("protein_binding_low_pct")
    # Most permissive plausible unbound upper bound uses the LOWEST stated binding.
    unbound_upper=None
    if total_uM is not None and bind_lo is not None:
        unbound_upper=total_uM*max(0.0,1.0-float(bind_lo)/100.0)
    required=target_exp.get("lowest_direct_target_concentration_uM")
    ratio_total=ratio_unbound=None
    if required is not None and total_uM:
        ratio_total=float(required)/total_uM
    if required is not None and unbound_upper:
        ratio_unbound=float(required)/unbound_upper

    if required is None or total_uM is None:
        state="UNKNOWN"
    elif unbound_upper is not None and float(required)<=unbound_upper:
        state="EXPOSURE_PLAUSIBLE"
    elif float(required)<=total_uM:
        state="EXPOSURE_BORDERLINE_PROTEIN_BINDING"
    elif float(required)<=3.0*total_uM:
        state="EXPOSURE_BORDERLINE_ABOVE_TOTAL_CMAX"
    else:
        state="EXPOSURE_IMPLAUSIBLE_FROM_LABEL_CMAX"

    result={
      "schema":SCHEMA,"candidate":candidate,"target":target,
      "state":state,"molecular_weight_g_per_mol":mw,
      "molecular_weight_source":mw_source,
      "chembl_parent_molecular_weight_g_per_mol":chembl_mw,
      "chembl_id":cid,
      "lowest_direct_target_concentration_uM":required,
      "max_observed_label_total_cmax_ug_per_mL":cmax,
      "max_observed_label_total_cmax_uM":round(total_uM,3) if total_uM is not None else None,
      "protein_binding_low_pct":bind_lo,"protein_binding_high_pct":bind_hi,
      "estimated_unbound_cmax_upper_uM":round(unbound_upper,3) if unbound_upper is not None else None,
      "required_to_total_cmax_ratio":round(ratio_total,3) if ratio_total is not None else None,
      "required_to_unbound_upper_ratio":round(ratio_unbound,3) if ratio_unbound is not None else None,
      "target_evidence":target_exp,"label_exposure":label_exp,
      "total_plasma_is_not_free_target_exposure":True,
      "tissue_exposure_not_established":True,
      "dose_recommendation":False,"clinical_recommendation":False,
      "interpretation":"Research exposure plausibility only; not a dose, safety clearance, or efficacy estimate."}
    (out/"EXPOSURE_FEASIBILITY_V31.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return result

def annotate_track2(track2,assessments):
    amap={x.get("candidate"):x for x in assessments or []}
    def ann(x0):
        x=dict(x0);a=amap.get(x.get("candidate"))
        if a:
            x["exposure_feasibility"]=a
            state=a.get("state")
            x["exposure_submission_hold"]=state=="EXPOSURE_IMPLAUSIBLE_FROM_LABEL_CMAX"
            if state and state.startswith("EXPOSURE_BORDERLINE"):
                x["exposure_uncertainty_requires_resolution"]=True
        return x
    out=dict(track2)
    out["single_agent_hypotheses"]=[ann(x) for x in track2.get("single_agent_hypotheses",[])]
    out["submission_level_single_agents"]=[
      ann(x) for x in track2.get("submission_level_single_agents",[])
      if not (amap.get(x.get("candidate"),{}).get("state")=="EXPOSURE_IMPLAUSIBLE_FROM_LABEL_CMAX")]
    out["exposure_feasibility_applied"]=True
    out["has_submission_level_hypothesis"]=bool(
      out["submission_level_single_agents"] or out.get("combination_hypotheses"))
    return out


def load_curated_record(project_root,candidate,target):
    p=Path(project_root)/"config/v31_curated_quantitative_evidence.json"
    if not p.is_file():return None
    try:data=json.loads(p.read_text(encoding="utf-8"))
    except Exception:return None
    key=f"{str(candidate).upper()}|{str(target).upper()}"
    rec=(data.get("records") or {}).get(key)
    return dict(rec) if rec else None

def _exposure_adjustment(state):
    return {
      "EXPOSURE_PLAUSIBLE":2.0,
      "EXPOSURE_BORDERLINE_PROTEIN_BINDING":-5.0,
      "EXPOSURE_BORDERLINE_ABOVE_TOTAL_CMAX":-10.0,
      "EXPOSURE_IMPLAUSIBLE_FROM_LABEL_CMAX":-30.0,
      "UNKNOWN":-2.0,
      None:0.0,
    }.get(state,-2.0)

def annotate_track2(track2,assessments):
    amap={x.get("candidate"):x for x in assessments or []}
    def ann(x0):
        x=dict(x0);a=amap.get(x.get("candidate"))
        if a:
            x["exposure_feasibility"]=a
            state=a.get("state")
            x["exposure_submission_hold"]=state=="EXPOSURE_IMPLAUSIBLE_FROM_LABEL_CMAX"
            if state and state.startswith("EXPOSURE_BORDERLINE"):
                x["exposure_uncertainty_requires_resolution"]=True
            base=(x.get("challenge_priority_score_source_adjusted")
                  if x.get("challenge_priority_score_source_adjusted") is not None
                  else x.get("challenge_priority_score") or 0)
            x["challenge_priority_score_exposure_adjusted"]=round(
                float(base)+_exposure_adjustment(state),3)
        else:
            x["exposure_feasibility"]={"state":"NOT_ASSESSED"}
            x["exposure_submission_hold"]=False
        return x
    out=dict(track2)
    out["single_agent_hypotheses"]=[ann(x) for x in track2.get("single_agent_hypotheses",[])]
    sub=[ann(x) for x in track2.get("submission_level_single_agents",[])]
    sub=[x for x in sub if not x.get("exposure_submission_hold")]
    sub.sort(key=lambda x:(
      -float(x.get("challenge_priority_score_exposure_adjusted",
                   x.get("challenge_priority_score_source_adjusted",
                         x.get("challenge_priority_score",0))) or 0),
      str(x.get("candidate") or "")))
    out["submission_level_single_agents"]=sub
    out["exposure_feasibility_applied"]=True
    out["has_submission_level_hypothesis"]=bool(
      sub or out.get("combination_hypotheses"))
    out["clinical_safety_clearance"]=False
    return out


def _alias_pattern(alias):
    parts=[re.escape(x) for x in re.split(r"[\s_-]+",str(alias or "").strip()) if x]
    if not parts:return None
    return r"(?<![A-Za-z0-9])"+r"[-\s_]*".join(parts)+r"(?![A-Za-z0-9])"

def _deterministic_percent_control_records(text,candidate,target_aliases):
    """Extract kinase-panel percent-of-control rows without an LLM."""
    out=[];aliases=[x for x in (target_aliases or []) if x]
    if not text or not aliases:return out
    lo=text.lower()
    for alias in aliases:
        ap=_alias_pattern(alias)
        if not ap:continue
        for m in re.finditer(ap,text,re.I):
            a=max(0,m.start()-2200);b=min(len(text),m.end()+900)
            window=text[a:b]
            wlo=window.lower()
            if "percent of control" not in wlo and "% of control" not in wlo:
                continue
            if "kinase" not in wlo:
                continue
            if str(candidate or "").lower() not in wlo:
                continue
            # Prefer an explicit candidate concentration in the same table/context.
            concs=[]
            for cm in re.finditer(
                re.escape(str(candidate or ""))+r".{0,80}?(\d+(?:\.\d+)?)\s*(µm|μm|um|nm)",
                window,re.I|re.S):
                val=float(cm.group(1));unit=cm.group(2).lower()
                if unit=="nm":val/=1000.0
                concs.append(val)
            if not concs:
                for cm in re.finditer(r"(\d+(?:\.\d+)?)\s*(µm|μm|um|nm)",window,re.I):
                    val=float(cm.group(1));unit=cm.group(2).lower()
                    if unit=="nm":val/=1000.0
                    concs.append(val)
            if not concs:continue
            # Parse the numeric score immediately after the target row label.
            tail=window[m.start()-a+len(m.group(0)):m.start()-a+len(m.group(0))+80]
            pm=re.search(r"(?:\([A-Za-z]+\))?\s*([0-9]+(?:\.[0-9]+)?)",tail)
            if not pm:continue
            pct=float(pm.group(1))
            if pct<0 or pct>250:continue
            action="decrease" if pct<100 else "increase" if pct>100 else "unknown"
            rec={
              "value":min(concs),"unit":"uM",
              "assay_context":"BIOCHEMICAL_TARGET_ACTIVITY",
              "target_action":action,
              "effect_fraction_of_control":pct/100.0,
              "effect_description":f"Deterministic kinase-panel extraction: target activity {pct}% of control.",
              "explicit_numeric_in_source":True,
              "deterministic_source_parser":True}
            key=(rec["value"],rec["target_action"],rec["effect_fraction_of_control"])
            if not any((x.get("value"),x.get("target_action"),x.get("effect_fraction_of_control"))==key for x in out):
                out.append(rec)
    return out


def assess_track2(track2,entities,raw_labels,outdir):
    """Assess exposure only for submission hypotheses backed by a cross-source bridge."""
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    emap={e.get("canonical_label"):e for e in entities or []}
    seen=set();assessments=[]
    rows=(track2.get("submission_level_single_agents") or [])
    for x in rows:
        if not x.get("source_bound_cross_source_bridge"):continue
        cand=x.get("candidate")
        if not cand or cand in seen:continue
        seen.add(cand)
        entity=emap.get(cand)
        label=(raw_labels or {}).get(cand)
        source=x.get("cross_source_bridge_drug_target_source") or {}
        bridge={
          "target_symbol":x.get("route_target"),
          "target_aliases":x.get("cross_source_bridge_target_aliases") or [],
          "drug_target_source":source}
        if not entity or not label:
            assessments.append({
              "schema":SCHEMA,"candidate":cand,"target":x.get("route_target"),
              "state":"UNKNOWN","reason":"approved entity or public label unavailable",
              "dose_recommendation":False,"clinical_recommendation":False})
            continue
        assessments.append(assess(
            cand,entity,bridge,label,out/cand))
    annotated=annotate_track2(track2,assessments)
    (out/"EXPOSURE_FEASIBILITY_SUMMARY_V31.json").write_text(
      json.dumps({"schema":SCHEMA,"assessments":assessments,
                  "clinical_recommendation":False,"dose_recommendation":False},
                 indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return annotated,assessments


def molecular_weight_from_label(label):
    """Prefer the molecular weight explicitly stated for the regulated formulation."""
    vals=[]
    for key in ("description","clinical_pharmacology"):
        v=(label or {}).get(key)
        if isinstance(v,list):vals.extend(str(x) for x in v)
        elif isinstance(v,str):vals.append(v)
    text=" ".join(vals)
    pats=[
      r"(?:calculated\s+)?molecular\s+weight\s+(?:of|is|=|:)?\s*([0-9]+(?:\.[0-9]+)?)",
      r"\bMW\s*(?:=|:)?\s*([0-9]+(?:\.[0-9]+)?)"]
    for pat in pats:
        m=re.search(pat,text,re.I)
        if m:
            try:
                mw=float(m.group(1))
                if 50<=mw<=5000:return mw
            except Exception:pass
    return None
