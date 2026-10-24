"""V3.1 third-party rights audit and notices generator.

Conservative engineering compliance gate; not legal advice.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,json
from . import license_gate

SCHEMA="RDMVA_THIRD_PARTY_RIGHTS_AUDIT_V31_1"

def _sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

def build(policy_path=None):
    path=Path(policy_path or license_gate.POLICY_PATH)
    policy=license_gate.load_policy(path)
    gate=license_gate.assert_default_profile(policy)
    sources=policy.get("sources") or {}
    rows=[]
    for name,rec in sorted(sources.items()):
        cls=str(rec.get("class") or "")
        default=bool(rec.get("enabled_by_default"))
        rows.append({
          "source":name,
          "class":cls,
          "enabled_by_default":default,
          "license_or_terms":rec.get("license_or_terms"),
          "commercial_use":rec.get("commercial_use"),
          "clinical_use_validated":rec.get("clinical_use_validated"),
          "redistribution_strategy":rec.get("redistribution_strategy"),
          "default_profile_allowed":bool(
              default and cls!="RED" and not cls.startswith("RED_")
              and cls!="YELLOW_LICENSE_CONFLICT"),
          "requires_attribution":(
              "ATTRIBUTION" in cls or
              " BY " in (" "+str(rec.get("license_or_terms") or "").upper()+" "))
        })
    af=sources.get("AlphaFold_DB") or {}
    am=sources.get("AlphaMissense") or {}
    findings={
      "alphafold_db":{
        "status":"ALLOWED_WITH_ATTRIBUTION" if af.get("enabled_by_default") else "DISABLED",
        "commercial_use_permitted":bool(af.get("commercial_use")),
        "default_enabled":bool(af.get("enabled_by_default")),
        "clinical_validation_claimed":False,
        "required_controls":[
          "retain source/accession/retrieval-date attribution",
          "do not present predicted structure as clinical validation",
          "do not imply third-party-rights warranty",
          "prefer runtime retrieval over bundling large mirrored collections"]},
      "alphamissense":{
        "status":"DISABLED_LICENSE_CONFLICT",
        "default_enabled":bool(am.get("enabled_by_default")),
        "reason":am.get("license_or_terms"),
        "commercial_profile_allowed":False,
        "challenge_baseline_allowed":False,
        "required_resolution":"rights-holder clarification or unambiguous upstream artifact terms before enabling"}
    }
    failures=[]
    if findings["alphamissense"]["default_enabled"]:
        failures.append("AlphaMissense_must_be_disabled")
    if not findings["alphafold_db"]["commercial_use_permitted"]:
        failures.append("AlphaFold_DB_commercial_status_not_green")
    if any(x["enabled_by_default"] and
           (x["class"]=="RED" or x["class"].startswith("RED_") or
            x["class"]=="YELLOW_LICENSE_CONFLICT") for x in rows):
        failures.append("restricted_source_enabled_by_default")
    return {
      "schema":SCHEMA,
      "policy_version":policy.get("version"),
      "policy_sha256":_sha256(path),
      "default_license_gate":gate,
      "source_count":len(rows),
      "sources":rows,
      "special_findings":findings,
      "challenge_public_profile_pass":not failures,
      "failures":failures,
      "legal_advice":False,
      "interpretation":"Conservative engineering rights audit; consult qualified counsel for legal conclusions."
    }

def write(audit_path,notices_path,policy_path=None):
    audit=build(policy_path)
    ap=Path(audit_path);ap.parent.mkdir(parents=True,exist_ok=True)
    ap.write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=[
      "# Rare Disease MVA Pharma — Third-Party Notices (V3.1)","",
      "This file documents third-party source/tool terms used by the research pipeline. It is an engineering compliance record, not legal advice.","",
      "Policy version: "+str(audit.get("policy_version")),"",
      "## Default Challenge/public profile","",
      "- Rights gate: "+("PASS" if audit.get("challenge_public_profile_pass") else "FAIL"),
      "- Third-party bulk datasets are not silently relicensed under the project's own licence.",
      "- Gated/patient data and API secrets are never redistribution assets.","",
      "## Sources",""]
    for x in audit.get("sources",[]):
        lines += [
          "### "+str(x["source"]),
          "- Class: "+str(x["class"]),
          "- Enabled by default: "+str(x["enabled_by_default"]).lower(),
          "- Terms: "+str(x.get("license_or_terms") or "See source-specific policy record."),
          "- Redistribution strategy: "+str(x.get("redistribution_strategy") or "Not specified."),""]
    lines += [
      "## AlphaFold DB","",
      "- Default profile: allowed with attribution.",
      "- Prediction data are treated as theoretical research context, not clinically validated evidence.",
      "- Selected structures/metadata should be retrieved at runtime with accession/source attribution.","",
      "## AlphaMissense","",
      "- Disabled by default.",
      "- Current project policy records conflicting official notices for the precomputed prediction artifact.",
      "- It is not used for the Challenge/public/commercial baseline until the rights conflict is resolved.","",
      "## ChEMBL / DGIdb","",
      "- ChEMBL remains external/source-specific data with its own share-alike obligations; the source package does not bundle an adapted ChEMBL database.",
      "- DGIdb is discovery/index only; an aggregator-only interaction cannot become submission evidence without primary-source verification."
    ]
    np=Path(notices_path);np.parent.mkdir(parents=True,exist_ok=True)
    np.write_text("\n".join(lines)+"\n",encoding="utf-8")
    return audit
