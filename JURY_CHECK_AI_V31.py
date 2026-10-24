#!/usr/bin/env python3
"""Private, optional LIVE model readiness check. Does not load patient data.

Use --verify-calls only with the reviewer's explicit consent to tiny billable
synthetic OpenAI requests. Does not print, cache or return the API key.
"""
from __future__ import annotations
import argparse,datetime,json,os,sys,urllib.request,urllib.error
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"code"))
from v2810.ai import api_key,MODELS

MODELS_CONFIG=ROOT/"config/openai_v30.json"
CHECK_PATH=(Path(os.environ.get("LOCALAPPDATA",Path.home()/"AppData/Local"))/
            "Examedix"/"RareDiseaseMVAPharma"/"JURY_AI_MODEL_CHECK_V31.json")

def required_models():
    cfg=json.loads(MODELS_CONFIG.read_text(encoding="utf-8"))
    roles={
       "v30_default":cfg["default_model"],
       "v30_critic":cfg["critic_model"],
       "v30_literature":cfg["literature_extraction_model"],
    }
    for role,definition in MODELS.items():
        roles["v2810_"+role]=definition["model"]
    return dict(sorted(roles.items()))

def request_json(method,url,key,body=None,timeout=60):
    headers={"Authorization":"Bearer "+key}
    data=None
    if body is not None:
        headers["Content-Type"]="application/json"
        data=json.dumps(body,separators=(",",":"),sort_keys=True).encode("utf-8")
    req=urllib.request.Request(url,data=data,headers=headers,method=method)
    with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req,timeout=timeout) as res:
        payload=res.read(4_000_000)
        return json.loads(payload)

def validate(*,test_responses=False,report_path=CHECK_PATH):
    key=api_key() # Windows-user DPAPI or explicit environment variable
    roles=required_models()
    result={"schema":"RDMVA_JURY_MODEL_VALIDATION_V31_1",
       "utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
       "tested_using_patient_data":False,"contains_secret":False,
       "billable_synthetic_calls":bool(test_responses),"roles":roles,"models":{}}
    try:
        discovery=request_json("GET","https://api.openai.com/v1/models",key,timeout=35)
        known={x["id"] for x in discovery.get("data",[]) if isinstance(x,dict) and isinstance(x.get("id"),str)}
    except (OSError,ValueError,KeyError,urllib.error.HTTPError) as exc:
        print("AI_MODEL_DISCOVERY_UNAVAILABLE",type(exc).__name__,flush=True)
        return False,result
    targets=sorted(set(roles.values()))
    for model in targets:
        observed={"listed_by_api":model in known,"live_response_tested":False,
                  "live_response_completed":None}
        if test_responses:
            # The only content sent is identical, non-clinical synthetic text.
            body={"model":model,"store":False,"tools":[],
                "input":"Synthetic public connectivity test. Respond only with OK.",
                "max_output_tokens":256}
            try:
                response=request_json("POST","https://api.openai.com/v1/responses",key,body,timeout=90)
                observed["live_response_tested"]=True
                observed["live_response_completed"]=(
                    response.get("status")=="completed" and
                    str(response.get("id") or "").startswith("resp_"))
                observed["returned_model"]=response.get("model")
                observed["response_id"]=response.get("id")
                usage=response.get("usage") or {}
                observed["usage_tokens"]={
                    "input":usage.get("input_tokens"),"output":usage.get("output_tokens")}
            except urllib.error.HTTPError as exc:
                observed["live_response_tested"]=True
                observed["live_response_completed"]=False
                observed["http_status"]=exc.code
            except (OSError,ValueError,KeyError,TimeoutError) as exc:
                observed["live_response_tested"]=True
                observed["live_response_completed"]=False
                observed["error_class"]=type(exc).__name__
        result["models"][model]=observed
        status=("PASS" if observed["live_response_completed"] else
                "FAIL" if test_responses else
                "LISTED" if observed["listed_by_api"] else "NOT_LISTED")
        print("AI_MODEL",model,status,flush=True)
    # GET model listing is informative; live Responses calls determine whether
    # the API can actually execute each configured model under this account.
    ok=all(x["live_response_completed"] is True for x in result["models"].values()) if test_responses else all(
        x["listed_by_api"] for x in result["models"].values())
    result["passed"]=bool(ok)
    if test_responses:
        report_path=Path(report_path)
        report_path.parent.mkdir(parents=True,exist_ok=True)
        report_path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return ok,result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--verify-calls",action="store_true",
                   help="Run small billable synthetic Responses API calls for each required model")
    args=p.parse_args()
    try: ok,_=validate(test_responses=args.verify_calls)
    except (OSError,ValueError,RuntimeError) as exc:
        print("AI_MODEL_VALIDATION_FAILED",type(exc).__name__,str(exc)[:100])
        return 2
    print("AI_MODEL_VALIDATION_PASS" if ok else "AI_MODEL_VALIDATION_INCOMPLETE")
    return 0 if ok else 2

if __name__=="__main__":
    raise SystemExit(main())
