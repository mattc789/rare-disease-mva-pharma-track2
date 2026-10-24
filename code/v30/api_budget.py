"""Per-run OpenAI budget guard for V3.0.

Tracks conservative reservations and actual token-estimated cost in a local
JSON ledger. It is not included in source releases or research results.
"""
from __future__ import annotations
from pathlib import Path
import json,os,hashlib,decimal

MICRO=decimal.Decimal("1000000")

def _money(x):return decimal.Decimal(str(x))

def _ledger_path(root):return Path(root)/"V30_API_BUDGET_LEDGER.json"

def _load(root):
    p=_ledger_path(root)
    if p.is_file():return json.loads(p.read_text(encoding="utf-8"))
    return {"schema":"RDMVA_V30_API_BUDGET_1","requests":{}}

def _save(root,j):
    p=_ledger_path(root);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(".tmp")
    tmp.write_text(json.dumps(j,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    tmp.replace(p)

def cap_usd():
    raw=os.environ.get("RDMVA_OPENAI_BUDGET_USD","3.0")
    cap=_money(raw)
    if cap<=0 or cap>1000:raise ValueError("V30_API_BUDGET_INVALID")
    return cap

def reserve(root,body_bytes,max_output_tokens,input_rate,output_rate,model):
    j=_load(root);h=hashlib.sha256(body_bytes).hexdigest()
    if h in j["requests"]:return h,j["requests"][h].get("reserved_usd",0),True
    # deliberately conservative input upper bound, same principle as V2.8.10
    # Token count cannot exceed UTF-8 byte count for normal text payloads; add a fixed safety margin.
    input_token_upper=max(1,len(body_bytes)+8192)
    est=(_money(input_token_upper)*_money(input_rate)+
         _money(max_output_tokens)*_money(output_rate))/MICRO
    reserved=sum(_money(x.get("reserved_usd",0)) for x in j["requests"].values())
    if reserved+est>cap_usd():
        raise RuntimeError("V30_OPENAI_BUDGET_HOLD")
    rec={"model":model,"reserved_usd":float(est),"actual_estimate_usd":None,"status":"reserved"}
    j["requests"][h]=rec;_save(root,j)
    return h,float(est),False

def finish(root,h,usage,input_rate,output_rate,status="completed"):
    j=_load(root);rec=j["requests"].get(h)
    if not rec:return
    actual=(_money((usage or {}).get("input_tokens",0))*_money(input_rate)+
            _money((usage or {}).get("output_tokens",0))*_money(output_rate))/MICRO
    rec["actual_estimate_usd"]=float(actual);rec["status"]=status
    _save(root,j)

def fail(root,h,status="transport_or_output_failure"):
    j=_load(root);rec=j["requests"].get(h)
    if rec:rec["status"]=status;_save(root,j)

def status(root):
    j=_load(root);rows=list(j["requests"].values())
    return {"cap_usd":float(cap_usd()),"requests":len(rows),
      "reserved_usd":round(sum(float(x.get("reserved_usd") or 0) for x in rows),6),
      "actual_estimate_usd":round(sum(float(x.get("actual_estimate_usd") or 0) for x in rows),6)}
