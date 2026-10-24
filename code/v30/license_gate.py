"""Machine-readable third-party licensing gate for V3.0."""
from __future__ import annotations
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[2]
POLICY_PATH=ROOT/"config/v30_source_policy.json"

def load_policy(path=POLICY_PATH):
    p=Path(path)
    data=json.loads(p.read_text(encoding="utf-8"))
    if data.get("schema")!="RDMVA_V30_SOURCE_POLICY_1":
        raise ValueError("V30_SOURCE_POLICY_SCHEMA_UNSUPPORTED")
    return data

def source_status(name,policy=None):
    data=policy or load_policy()
    rec=(data.get("sources") or {}).get(name)
    if not rec: return {"status":"UNREGISTERED","allowed_default":False}
    return dict(rec)

def assert_default_profile(policy=None):
    data=policy or load_policy()
    bad=[]
    for name,rec in (data.get("sources") or {}).items():
        cls=str(rec.get("class") or "")
        if rec.get("enabled_by_default") and (cls=="RED" or cls.startswith("RED_") or cls=="YELLOW_LICENSE_CONFLICT"):
            bad.append(name)
        if rec.get("enabled_by_default") and not rec.get("redistribution_strategy"):
            bad.append(name+":missing_strategy")
    if bad:
        raise RuntimeError("V30_LICENSE_GATE_FAILED:"+",".join(sorted(bad)))
    return {
      "schema":"RDMVA_V30_LICENSE_GATE_RESULT_1",
      "policy_version":data.get("version"),
      "default_profile_valid":True,
      "red_sources_enabled":False,
      "license_conflict_sources_enabled":False,
      "source_count":len(data.get("sources") or {}),
    }
