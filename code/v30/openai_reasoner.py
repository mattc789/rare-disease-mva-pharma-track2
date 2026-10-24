"""Optional OpenAI research critic for V3.0.

Only a manually constructed sanitized envelope is transmitted. The model may
identify gaps and falsification priorities; it cannot upgrade deterministic
evidence tiers or create a clinical recommendation.
"""
from __future__ import annotations
from pathlib import Path
import json,os,re,urllib.request,urllib.error,hashlib,time,socket
from v2810 import ai as v28_ai
from . import api_budget

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/"config/openai_v30.json"
SCHEMA="RDMVA_OPENAI_CRITIC_V30_1"

FORBIDDEN_KEY_PARTS=("chrom","position","pos","ref","alt","hgv","genotype","phase_tag",
                     "variant_ref","case_id","patient_id","phenotype","hpo","vcf","bam","cram")
FORBIDDEN_TEXT=[
 re.compile(r"(?i)\b(?:chr)?(?:[1-9]|1\d|2[0-2]|X|Y|MT):\d+:[ACGTN]+:[ACGTN]+\b"),
 re.compile(r"(?i)\b(?:c|g|p|n|r|m)\.[A-Za-z0-9_*+>\-]+"),
 re.compile(r"\b[012][|/][012]\b"),
 re.compile(r"\bHP:\d{7}\b",re.I),
 re.compile(r"#CHROM\b",re.I),
]

def config():
    d=json.loads(CONFIG.read_text(encoding="utf-8"))
    if d.get("schema")!="RDMVA_OPENAI_V30_1":raise ValueError("V30_OPENAI_CONFIG_SCHEMA")
    return d

def _assert_safe(obj,path="$"):
    if isinstance(obj,dict):
        for k,v in obj.items():
            low=str(k).lower()
            if any(part in low for part in FORBIDDEN_KEY_PARTS):
                raise ValueError("V30_API_PRIVACY_FORBIDDEN_KEY:"+path+"."+str(k))
            _assert_safe(v,path+"."+str(k))
    elif isinstance(obj,list):
        for i,v in enumerate(obj):_assert_safe(v,path+f"[{i}]")
    elif isinstance(obj,str):
        for rx in FORBIDDEN_TEXT:
            if rx.search(obj):raise ValueError("V30_API_PRIVACY_FORBIDDEN_TEXT:"+path)

def safe_envelope(patient_state,mechanisms,objectives,falsification,max_candidates=8):
    gene_state={x["gene"]:x for x in patient_state.get("genes",[])}
    mech_by={x["gene"]:[] for x in patient_state.get("genes",[])}
    for m in mechanisms.get("hypotheses",[]):
        if m.get("status")!="NOT_INFERRED_WITHOUT_DIRECT_EVIDENCE":
            mech_by.setdefault(m.get("gene"),[]).append({
              "mechanism_type":m.get("mechanism_type"),"status":m.get("status"),
              "defect_direction":m.get("defect_direction")})
    obj_by={}
    for o in objectives.get("objectives",[]):
        obj_by.setdefault(o.get("gene"),[]).append({
          "objective_type":o.get("objective_type"),"desired_action":o.get("desired_action"),
          "variable":o.get("variable"),"priority_tier":o.get("priority_tier")})
    genes=[]
    for gene,g in gene_state.items():
        genes.append({"gene":gene,"effect_classes":g.get("effect_classes"),
          "evidence_completeness_tier":g.get("evidence_completeness_tier"),
          "mechanisms":mech_by.get(gene,[]),"rescue_objectives":obj_by.get(gene,[])})
    candidates=[]
    for c in falsification.get("candidates",[])[:max_candidates]:
        candidates.append({
          "candidate":c.get("candidate"),"gene":c.get("gene"),
          "route_target":c.get("route_target"),"objective_type":c.get("objective_type"),
          "decision":c.get("decision"),"supporting_factors":c.get("supporting_factors"),
          "hold_reasons":c.get("hold_reasons"),"kill_reasons":c.get("kill_reasons"),
        })
    payload={"schema":"RDMVA_V30_SANITIZED_AI_ENVELOPE_1","genes":genes,"candidates":candidates,
             "instructional_context":"Research hypothesis critique only; no treatment, dose, efficacy probability or patient identification."}
    _assert_safe(payload)
    return payload

def _response_text(raw):
    chunks=[]
    for item in raw.get("output",[]):
        for c in item.get("content",[]) if isinstance(item,dict) else []:
            if c.get("type") in ("output_text","text") and isinstance(c.get("text"),str):
                chunks.append(c["text"])
    if not chunks and isinstance(raw.get("output_text"),str):chunks.append(raw["output_text"])
    return "".join(chunks)

def request_json(body,outfile,budget_root=None):
    """Budgeted, cached Responses API call returning parsed structured JSON."""
    cfg=config();outfile=Path(outfile);outfile.parent.mkdir(parents=True,exist_ok=True)
    model=str(body.get("model") or cfg["default_model"])
    pricing=(cfg.get("model_pricing_usd_per_million") or {}).get(model)
    if not pricing:raise RuntimeError("V30_OPENAI_MODEL_PRICING_NOT_CONFIGURED:"+model)
    root=Path(os.environ.get('RDMVA_E2E_GLOBAL_AI_BUDGET_ROOT') or budget_root or outfile.parent)
    rawbody=json.dumps(body,sort_keys=True,ensure_ascii=False).encode("utf-8")
    h=hashlib.sha256(rawbody).hexdigest()
    cache=root/"V30_API_CACHE"/(h+".json")
    if cache.is_file():
        rec=json.loads(cache.read_text(encoding="utf-8"))
        replay=dict(rec);replay["replayed_from_cache"]=True
        outfile.write_text(json.dumps(replay,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
        return replay
    maxout=int(body.get("max_output_tokens") or cfg.get("max_output_tokens",4500))
    api_budget.reserve(root,rawbody,maxout,pricing["input"],pricing["output"],model)
    req=urllib.request.Request(cfg["endpoint"],data=rawbody,method="POST",
        headers={"Authorization":"Bearer "+v28_ai.api_key(),"Content-Type":"application/json"})
    raw=None;last_exc=None
    for attempt in range(4):
        try:
            with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(req,timeout=120) as r:
                raw=json.loads(r.read(2_000_000))
            break
        except urllib.error.HTTPError as exc:
            last_exc=exc
            if exc.code not in (408,409,429,500,502,503,504):
                api_budget.fail(root,h);raise
        except (urllib.error.URLError,TimeoutError,socket.timeout) as exc:
            last_exc=exc
        except Exception:
            api_budget.fail(root,h);raise
        if attempt<3:
            time.sleep(2**attempt)
    if raw is None:
        api_budget.fail(root,h)
        raise RuntimeError("V30_OPENAI_TRANSIENT_FAILURE_AFTER_RETRIES") from last_exc
    usage=raw.get("usage") or {}
    if raw.get("status") not in (None,"completed"):
        api_budget.finish(root,h,usage,pricing["input"],pricing["output"],"incomplete")
        raise RuntimeError("V30_OPENAI_INCOMPLETE_RESPONSE")
    try:value=json.loads(_response_text(raw))
    except Exception:
        api_budget.finish(root,h,usage,pricing["input"],pricing["output"],"invalid_json")
        raise RuntimeError("V30_OPENAI_INVALID_STRUCTURED_JSON")
    api_budget.finish(root,h,usage,pricing["input"],pricing["output"],"completed")
    rec={"request_sha256":h,"model":raw.get("model") or model,"response_id":raw.get("id"),
         "usage":usage,"value":value,"store":False,"patient_raw_data_transmitted":False}
    cache.parent.mkdir(parents=True,exist_ok=True)
    cache.write_text(json.dumps(rec,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    outfile.write_text(json.dumps(rec,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return rec

def review(patient_state,mechanisms,objectives,falsification,outfile):
    cfg=config();payload=safe_envelope(patient_state,mechanisms,objectives,falsification)
    model=os.environ.get(cfg.get("environment_model_override","RDMVA_OPENAI_MODEL"),
                         cfg.get("critic_model") or cfg["default_model"])
    schema={
      "type":"object","additionalProperties":False,
      "properties":{
        "candidate_reviews":{"type":"array","items":{
          "type":"object","additionalProperties":False,
          "properties":{
            "candidate":{"type":"string"},"causal_gaps":{"type":"array","items":{"type":"string"}},
            "counter_hypotheses":{"type":"array","items":{"type":"string"}},
            "highest_information_falsification_tests":{"type":"array","items":{"type":"string"}},
            "reason_to_hold_or_kill":{"type":"string"}
          },
          "required":["candidate","causal_gaps","counter_hypotheses",
                      "highest_information_falsification_tests","reason_to_hold_or_kill"]}},
        "global_missing_evidence":{"type":"array","items":{"type":"string"}},
        "new_research_questions":{"type":"array","items":{"type":"string"}}
      },
      "required":["candidate_reviews","global_missing_evidence","new_research_questions"]
    }
    instructions=(
      "Act as a skeptical molecular pharmacology research critic. Use only the supplied sanitized abstract evidence. "
      "Do not infer patient details, treatment suitability, dose, clinical benefit, or efficacy probabilities. "
      "Do not upgrade evidence tiers. Identify causal gaps, plausible counter-hypotheses and the smallest high-information "
      "experiments that could falsify each candidate mechanism. If evidence is weak, say HOLD or KILL rather than invent support."
    )
    body={"model":model,"input":json.dumps(payload,sort_keys=True),"instructions":instructions,
          "store":False,"tools":[],"max_output_tokens":int(cfg.get("max_output_tokens",4500)),
          "text":{"format":{"type":"json_schema","name":"v30_research_critic","strict":True,"schema":schema}}}
    rec=request_json(body,outfile,budget_root=Path(outfile).parent)
    result={"schema":SCHEMA,"model":rec.get("model") or model,"response_id":rec.get("response_id"),
            "value":rec["value"],"patient_raw_data_transmitted":False,
            "deterministic_evidence_tiers_modified":False,"clinical_recommendation":False,
            "api_budget":api_budget.status(Path(outfile).parent)}
    Path(outfile).write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return result
