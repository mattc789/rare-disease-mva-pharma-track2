#!/usr/bin/env python3
"""Validate autonomous V3.1 execution using current run evidence and receipts.

A technical PASS never depends on the name of a particular drug. Full clinical
validity, the truth of external sources, and zero leakage cannot be established
merely by inspecting a run's output flags.
"""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json

SCHEMA="RDMVA_CLEANROOM_ATTESTATION_V31_2"

def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8-sig"))

def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def ai_receipts(run):
    """Audit direct API receipts and require summary references to be backed by one.

    V3.0/3.1 critic reports copy a response_id but do not issue an additional
    API request. They omit the store field even though their corresponding
    request receipt explicitly records store=False. Never treat an unverified
    report as proof of a compliant API request. No response IDs are exported.
    """
    ids=set();models=set();summary_refs=[];without_store_false=0
    direct_count=0
    for p in run.rglob("*.json"):
        if p.stat().st_size>6_000_000 or p.name in (
                "CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json",
                "RUN_REPRODUCIBILITY_CONTEXT_V31.json"):
            continue
        try:
            v=load(p)
        except (OSError,ValueError,UnicodeError):
            continue
        if not isinstance(v,dict):continue
        resp=v.get("response_id")
        if not isinstance(resp,str) or not resp.startswith("resp_"):continue
        # Direct stored response receipts include the original request hash
        # and API usage. Critic reports are references, not request receipts.
        direct=(isinstance(v.get("request_sha256"),str) and
                isinstance(v.get("usage"),dict))
        if direct:
            direct_count+=1
            ids.add(resp)
            if isinstance(v.get("model"),str):models.add(v["model"])
            if v.get("store") is not False:without_store_false+=1
        else:
            # A report with store=True is inconsistent even if the request
            # receipt claims False. A report with no store field is valid ONLY
            # when its response ID matches an independently checked receipt.
            summary_refs.append((resp,v.get("store")))
    orphan_references=sum(resp not in ids for resp,_ in summary_refs)
    contradictory_references=sum(store is True for _,store in summary_refs)
    without_store_false+=orphan_references+contradictory_references
    return {"unique_api_response_ids":len(ids),
            "models_observed":sorted(models),
            "direct_receipt_records_checked":direct_count,
            "summary_references_verified":len(summary_refs)-orphan_references,
            "orphan_summary_response_references":orphan_references,
            "response_records_without_store_false":without_store_false,
            "response_ids_not_exported":True}

def verify(run,expected_path,require_ai):
    final_path=run/"FINAL_RESULT_V31_NATIVE.json"
    if not final_path.is_file():
        raise FileNotFoundError("Final V3.1 result was not generated")
    final=load(final_path)
    name=final.get("track2_submission_hypotheses")
    if not isinstance(name,str) or not name:
        raise ValueError("Final result lacks Track-2 hypotheses path")
    track_path=(run/name).resolve()
    if run.resolve() not in track_path.parents:
        raise ValueError("Track-2 path escapes run directory")
    track=load(track_path)
    stat_path=run/"STATUS_V30_NATIVE.json"
    status=load(stat_path) if stat_path.is_file() else {}
    ai=ai_receipts(run)
    supported=int(final.get("all_approved_entities_screened") or 0)
    min_entities=7000
    if expected_path:
        expected=load(expected_path)
        min_entities=int(expected.get("expected_minimum_approved_entities_screened") or 7000)
    else:
        expected={}
    provenance=run/"E2E_EXECUTION_V31.json"
    run_receipt=load(provenance) if provenance.is_file() else {}
    is_current=run_receipt.get("public_data_profile")=="CURRENT_OFFICIAL_SOURCES"
    # For current official source snapshots, require actual complete screening
    # of the source-bound medicine universe and nonzero multi-region coverage.
    master_path=run/"master/summary.json"
    universe=load(master_path) if master_path.is_file() else {}
    us_count=int(universe.get("us_snapshot_entity_count") or 0)
    extension_count=int(universe.get("new_entity_labels") or 0)
    live_coverage=(supported>0 and
                   supported==int(universe.get("entities") or 0) and
                   us_count>=2000 and extension_count>0 and
                   universe.get("dataset_profile")=="CURRENT_OFFICIAL_SOURCES" and
                   bool(run_receipt.get("current_public_data_receipt_sha256")))
    genomics=run/"GENOMIC_REQUIRED_EVIDENCE_V31.json"
    genomic_audit=load(genomics) if genomics.is_file() else {}
    checks={
        "required_genomic_evidence_full_profile":genomic_audit.get("required_science_pass") is True,
        "v31_final_result_present":final_path.is_file(),
        "v31_stage_completed":status.get("completed") is True and status.get("final_result")==final_path.name,
        "approved_entities_screened_minimum":live_coverage if is_current else supported>=min_entities,
        "track2_file_present":track_path.is_file(),
        "track2_hypotheses_object_valid":isinstance(track.get("submission_level_single_agents"),list),
        "raw_genome_not_declared_as_transmitted":final.get("patient_raw_data_sent_to_ai") is False,
        "clinical_recommendation_disabled":final.get("clinical_recommendation") is False,
        "source_policy_exists":(run/"V31_THIRD_PARTY_RIGHTS_AUDIT.json").is_file(),
        "online_ai_flag_correct":final.get("openai_enabled") is require_ai,
        "real_ai_response_receipts_present":ai["unique_api_response_ids"]>0 if require_ai else True,
        "ai_response_store_false":ai["response_records_without_store_false"]==0,
    }
    technical_pass=all(checks.values())
    agents=list(track.get("submission_level_single_agents") or [])
    lead=agents[0] if agents else {}
    evidence=(lead.get("direct_source_rescue_evidence_class") or
              ("CROSS_SOURCE_MECHANISTIC_RESCUE"
               if lead.get("source_bound_cross_source_bridge") else None))
    run_context={
        "schema":"RDMVA_RUN_REPRODUCIBILITY_CONTEXT_V31_1",
        "note":"Observed evidence from this run. No particular drug is required for technical success.",
        "current_public_data_profile":is_current,
        "primary_hypothesis":{"candidate":lead.get("candidate"),
                              "gene":lead.get("gene"),
                              "target":lead.get("route_target"),
                              "evidence_class":evidence},
        "approved_entities":{"observed":supported,
                             "observed_us_ingredients":us_count,
                             "observed_new_jurisdiction_ingredients":extension_count},
        "submission_level_single_agents":len(agents),
        "combination_hypotheses":len(track.get("combination_hypotheses") or []),
        "genomic_required_science_pass":genomic_audit.get("required_science_pass") is True,
        "ai_response_receipts_verified":ai["unique_api_response_ids"],
    }
    (run/"RUN_REPRODUCIBILITY_CONTEXT_V31.json").write_text(
        json.dumps(run_context,indent=2)+"\n",encoding="utf-8")
    attest={
        "schema":SCHEMA,
        "passed":technical_pass,
        "checks":checks,
        "ai_provenance_metadata":ai,
        "run_context_file":"RUN_REPRODUCIBILITY_CONTEXT_V31.json",
        "technical_success_is_not_clinical_validation":True,
        "privacy_field_is_a_program_declaration_not_network_capture":True,
        "final_result_file":final_path.name,
        "final_result_sha256":sha256(final_path),
        "track2_file":track_path.relative_to(run).as_posix(),
        "track2_sha256":sha256(track_path),
        "contains_patient_raw_data":False,
    }
    output=run/"CLEANROOM_REPRODUCIBILITY_ATTESTATION_V31.json"
    output.write_text(json.dumps(attest,indent=2)+"\n",encoding="utf-8")
    return technical_pass,output,ai["unique_api_response_ids"]

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--run-dir",required=True)
    p.add_argument("--expected",default=None)
    g=p.add_mutually_exclusive_group()
    g.add_argument("--require-ai",action="store_true",default=True)
    g.add_argument("--no-require-ai",dest="require_ai",action="store_false")
    a=p.parse_args()
    try:
        ok,out,num=verify(Path(a.run_dir).resolve(),a.expected,a.require_ai)
    except (ValueError,OSError,KeyError,TypeError) as e:
        print("CLEANROOM_REPRODUCTION_FAIL",type(e).__name__,str(e)[:300],flush=True)
        raise SystemExit(2)
    print("CLEANROOM_REPRODUCTION_PASS",ok,flush=True)
    print("AI_RESPONSE_RECEIPTS",num,flush=True)
    print("ATTESTATION",out,flush=True)
    raise SystemExit(0 if ok else 2)

if __name__=="__main__":main()
