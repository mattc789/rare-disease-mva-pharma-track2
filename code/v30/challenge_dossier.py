"""Track-2 submission dossier for Rare Disease MVA Pharma V3.1.

Transforms structured research hypotheses into a concise, falsifiable,
source-aware submission artifact. It never converts a research hypothesis into
a treatment recommendation or dosing instruction.
"""
from __future__ import annotations
from pathlib import Path
import json

SCHEMA="RDMVA_TRACK2_DOSSIER_V31_1"

def _gene_context(gene,patient_state,disease_graph):
    ps=next((x for x in patient_state.get("genes",[]) if x.get("gene")==gene),{})
    node=next((x for x in disease_graph.get("nodes",[])
               if x.get("id")=="GENE:"+str(gene)),{})
    return {
      "gene":gene,
      "genome_rank":node.get("genome_rank"),
      "rank_score_not_diagnostic":node.get("rank_score_not_diagnostic"),
      "phenotype_score":node.get("phenotype_score"),
      "variant_count":ps.get("variant_count"),
      "effect_classes":ps.get("effect_classes") or [],
      "clinvar_reviewed_pathogenic_present":
          (ps.get("evidence_flags") or {}).get("clinvar_reviewed_pathogenic_present"),
      "phase_state":ps.get("phase_state"),
      "mechanism_confirmed":bool(ps.get("mechanism_confirmed")),
      "diagnosis_asserted":False,
    }

def _sources(x):
    rows=[]
    if x.get("source_bound_cross_source_bridge"):
        for s in x.get("cross_source_bridge_rescue_sources") or []:
            rows.append({
              "role":"PRIMARY_GENE_OR_DISEASE_RESCUE",
              "pmid":s.get("pmid"),"pmcid":s.get("pmcid"),
              "title":s.get("title"),"license":s.get("license"),
              "text_redistributed":False})
        s=x.get("cross_source_bridge_drug_target_source") or {}
        if s:
            rows.append({
              "role":"APPROVED_DRUG_TARGET_MODULATION",
              "pmid":s.get("pmid"),"pmcid":s.get("pmcid"),
              "title":s.get("title"),"license":s.get("license"),
              "text_redistributed":False})
    ev=x.get("source_bound_public_evidence") or {}
    return rows,ev

def _tests(x):
    gene=x.get("gene");target=x.get("route_target")
    tests=[
      "Verify the implicated patient variant/mechanism with an orthogonal molecular assay before interpreting drug rescue.",
      "Test the candidate in a disease-relevant patient-derived or engineered cellular model and include vehicle and mechanistic controls.",
      "Measure target engagement at the same concentrations used for phenotypic rescue.",
      "Measure both intended rescue and plausible worsening endpoints; kill the hypothesis if the disease-relevant phenotype worsens.",
      "Replicate the result independently before any translational interpretation.",
    ]
    if gene=="BUB1B":
        tests += [
          "Quantify BUBR1 protein abundance and spindle-assembly/chromosome-congression phenotypes.",
          "Measure chromosome alignment, micronuclei/aneuploidy and cell viability in parallel so apparent rescue is not caused by selective loss of abnormal cells.",
        ]
        if str(target).upper() in ("AURKB","AURORA B","AURORA KINASE B"):
            tests += [
              "Measure Aurora-B target engagement (for example phospho-histone-H3/kinetochore signaling) together with chromosome-congression rescue.",
              "Test a clinically plausible free-exposure range and the biochemical reference concentration; kill the hypothesis if rescue occurs only above plausible exposure.",
              "Use a research Aurora-B inhibitor as a mechanistic positive control, not as a proposed treatment.",
            ]
    elif gene=="LZTR1":
        tests += [
          "Measure LZTR1 protein abundance/function and RAS-MAPK output in the same experimental system.",
          "Confirm that rescue is not restricted to an unrelated oncologic or species-specific context.",
        ]
    # preserve order, remove duplicates
    out=[];seen=set()
    for t in tests:
        if t not in seen:out.append(t);seen.add(t)
    return out

def _entry(x,patient_state,disease_graph):
    sources,source_state=_sources(x)
    return {
      "candidate":x.get("candidate"),
      "status":x.get("status"),
      "gene_context":_gene_context(x.get("gene"),patient_state,disease_graph),
      "causal_chain":{
        "gene":x.get("gene"),
        "rescue_objective":x.get("objective_type"),
        "route_target":x.get("route_target"),
        "desired_action":x.get("route_desired_action"),
        "candidate":x.get("candidate"),
      },
      "evidence_class":(
        x.get("direct_source_rescue_evidence_class")
        or ("CROSS_SOURCE_MECHANISTIC_RESCUE"
            if x.get("source_bound_cross_source_bridge") else None)),
      "source_bound_evidence_state":x.get("source_bound_evidence_state"),
      "sources":sources,
      "source_bound_public_evidence":source_state,
      "challenge_priority_score":x.get("challenge_priority_score"),
      "source_adjusted_score":x.get("challenge_priority_score_source_adjusted"),
      "exposure_adjusted_score":x.get("challenge_priority_score_exposure_adjusted"),
      "exposure_feasibility":x.get("exposure_feasibility") or {"state":"NOT_ASSESSED"},
      "public_label_safety_review":x.get("public_label_safety_review") or {"state":"NOT_REVIEWED"},
      "main_blockers":x.get("main_blockers") or [],
      "context_conflict_requires_resolution":bool(
          x.get("context_conflict_requires_resolution")),
      "decisive_falsification_tests":_tests(x),
      "hypothesis_for_follow_up_not_evidence_it_works":True,
      "clinical_recommendation":False,
      "dose_recommendation":False,
    }

def build(track2,patient_state,disease_graph,source_policy=None):
    subs=list(track2.get("submission_level_single_agents") or [])
    primary=_entry(subs[0],patient_state,disease_graph) if subs else None
    alternatives=[_entry(x,patient_state,disease_graph) for x in subs[1:5]]
    return {
      "schema":SCHEMA,
      "track2_scope":"market-approved medicines proposed as hypotheses for further investigation",
      "primary_hypothesis":primary,
      "alternative_hypotheses":alternatives,
      "combination_hypotheses":track2.get("combination_hypotheses") or [],
      "combination_proposed":bool(track2.get("combination_hypotheses")),
      "source_policy_version":(source_policy or {}).get("version"),
      "clinical_safety_clearance":False,
      "clinical_efficacy_claim":False,
      "dose_recommendation":False,
      "clinical_recommendation":False,
      "interpretation":"Submission-oriented research hypothesis only; experimental validation is required.",
    }

def write(outdir,dossier):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    (out/"TRACK2_RESEARCH_DOSSIER_V31.json").write_text(
      json.dumps(dossier,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    p=dossier.get("primary_hypothesis")
    lines=[
      "# Rare Disease MVA Pharma — Track 2 Research Dossier","",
      "**Research hypothesis only — not a prescription, dosing recommendation, efficacy claim, or safety clearance.**","",
    ]
    if not p:
        lines += ["## Result","No submission-level hypothesis passed the current gates."]
    else:
        chain=p["causal_chain"];exp=p.get("exposure_feasibility") or {};safety=p.get("public_label_safety_review") or {}
        lines += [
          "## Primary hypothesis",
          f"**Candidate:** {p['candidate']}",
          f"**Evidence class:** {p.get('evidence_class')}",
          f"**Gene / rescue chain:** {chain.get('gene')} → {chain.get('rescue_objective')} → {chain.get('route_target')} ({chain.get('desired_action')}) → {p['candidate']}",
          f"**Source-bound state:** {p.get('source_bound_evidence_state')}",
          f"**Exposure state:** {exp.get('state','NOT_ASSESSED')}",
          f"**Public-label safety state:** {safety.get('state','NOT_REVIEWED')}","",
          "### Why it is scientifically interesting",
          "- It is a market-approved medicine linked to an explicit patient-first rescue objective.",
        ]
        if p.get("evidence_class")=="CROSS_SOURCE_MECHANISTIC_RESCUE":
            lines += [
              "- One independent source supports perturbing the intermediate target to rescue disease-relevant biology.",
              "- A second independent source supports direct modulation of that target by the approved medicine.",
              "- The same medicine has not yet been shown to rescue this patient's disease model; that is the central testable hypothesis.",
            ]
        lines += ["","### Material limitations"]
        lines += ["- "+str(x).replace("_"," ") for x in p.get("main_blockers",[])]
        if exp.get("state") and exp.get("state")!="NOT_ASSESSED":
            lines += [
              f"- Exposure feasibility: {exp.get('state')}.",
              "- Total plasma concentration is not equivalent to free intracellular target exposure.",
            ]
        lines += ["","### Decisive experiments"]
        lines += ["- "+t for t in p.get("decisive_falsification_tests",[])]
        lines += ["","### Supporting sources"]
        for s in p.get("sources",[]):
            lines += [f"- {s.get('role')}: {s.get('title')} (PMID {s.get('pmid')}; PMCID {s.get('pmcid')}; source terms: {s.get('license')})."]
    alts=dossier.get("alternative_hypotheses") or []
    if alts:
        lines += ["","## Alternative hypotheses"]
        for a in alts:
            lines += [
              f"- **{a.get('candidate')}** — {a.get('gene_context',{}).get('gene')} / {a.get('evidence_class')}; "
              f"source state {a.get('source_bound_evidence_state')}; exposure {(a.get('exposure_feasibility') or {}).get('state','NOT_ASSESSED')}."
            ]
    if not dossier.get("combination_proposed"):
        lines += ["","## Combination result","No A+B or A+B+C combination passed the current complementary-mechanism gates; no combination is proposed."]
    lines += ["","## Rights / redistribution guardrails",
      "- Third-party source text is not redistributed in this dossier; only citation metadata, structured facts and original analysis are retained.",
      "- AlphaFold DB context, when used, is attribution-bearing prediction data and is not treated as clinical validation.",
      "- AlphaMissense is disabled in the default Challenge profile because official notices conflict on reuse terms.",
    ]
    (out/"TRACK2_RESEARCH_DOSSIER_V31.md").write_text(
      "\n".join(lines)+"\n",encoding="utf-8")
    return out/"TRACK2_RESEARCH_DOSSIER_V31.json"
