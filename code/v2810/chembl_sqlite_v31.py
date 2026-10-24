"""Official ChEMBL 37 SQLite adapter for REST-compatible, fully local research.
Only source code is redistributed. Database CC BY-SA 3.0 (ChEMBL/EMBL-EBI).
Queries are read-only and scoped to the immutable verified official archive.
"""
from __future__ import annotations
from pathlib import Path
import json,sqlite3,urllib.parse,os,contextlib

PROJECT=Path(__file__).resolve().parents[2]
ROOT=PROJECT/"public_inputs/third_party/chembl37"
DATABASE=ROOT/"chembl_37.db"
PROVENANCE=ROOT/"CHEMBL_SQLITE_PROVENANCE_V31.json"
BASE="https://www.ebi.ac.uk/chembl/api/data/"
RELEASE=37

def configured():
    return DATABASE.is_file() and PROVENANCE.is_file()

def connection():
    if not PROVENANCE.is_file():
        raise RuntimeError("V31_CHEMBL_SQLITE_UNVERIFIED_PROVENANCE")
    record=json.loads(PROVENANCE.read_text(encoding="utf-8"))
    if (record.get("schema")!="RDMVA_CHEMBL_SQLITE_OFFICIAL_V31_1" or
        record.get("chembl_release")!=RELEASE or
        record.get("archive_sha256")!="33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281"):
        raise RuntimeError("V31_CHEMBL_SQLITE_PROVENANCE_MISMATCH")
    if not DATABASE.is_file() or DATABASE.stat().st_size<100_000:
        raise RuntimeError("V31_CHEMBL_SQLITE_DATABASE_MISSING")
    uri=DATABASE.as_uri()+"?mode=ro&immutable=1"
    db=sqlite3.connect(uri,uri=True,timeout=90)
    db.row_factory=sqlite3.Row
    db.execute("PRAGMA query_only=ON")
    db.execute("PRAGMA busy_timeout=90000")
    return db

def _in(values):
    if not values:return "(NULL)",[]
    if len(values)>900:raise ValueError("V31_CHEMBL_SQLITE_TOO_MANY_IDS")
    return "("+",".join("?" for _ in values)+")",values

def _metadata(url,params,offset,limit,more):
    p=dict(params)
    p.pop("offset",None)
    p.pop("limit",None)
    p["limit"]=limit
    nxt=None
    if more:
        p["offset"]=offset+limit
        nxt="/chembl/api/data/"+url+"?"+urllib.parse.urlencode(p,safe=",")
    previous=None
    if offset:
        p["offset"]=max(0,offset-limit)
        previous="/chembl/api/data/"+url+"?"+urllib.parse.urlencode(p,safe=",")
    return {"limit":limit,"next":nxt,"previous":previous,"total_count":None}

def _molecules(db,params,limit,offset):
    conditions=[]
    args=[]
    if params.get("max_phase") is not None:
        conditions.append("md.max_phase=?")
        args.append(int(params["max_phase"]))
    ids=[x for x in params.get("molecule_chembl_id__in","").split(",") if x]
    if ids:
        place,vals=_in(ids)
        conditions.append("md.chembl_id IN "+place)
        args.extend(vals)
    where=" WHERE "+" AND ".join(conditions) if conditions else ""
    sql=(
        "SELECT md.molregno,md.chembl_id,md.pref_name,md.max_phase,"
        "md.molecule_type,md.withdrawn_flag,md.availability_type,"
        "md.veterinary,md.first_approval,md.oral,md.parenteral,"
        "md.topical,md.prodrug,mh.parent_molregno,mh.active_molregno,"
        "parent.chembl_id AS parent_id,active.chembl_id AS active_id,"
        "cs.canonical_smiles,cp.full_mwt,cp.mw_freebase "
        "FROM molecule_dictionary md "
        "LEFT JOIN molecule_hierarchy mh ON mh.molregno=md.molregno "
        "LEFT JOIN molecule_dictionary parent ON parent.molregno=mh.parent_molregno "
        "LEFT JOIN molecule_dictionary active ON active.molregno=mh.active_molregno "
        "LEFT JOIN compound_structures cs ON cs.molregno=md.molregno "
        "LEFT JOIN compound_properties cp ON cp.molregno=md.molregno "
        +where+" ORDER BY md.molregno LIMIT ? OFFSET ?")
    rows=[dict(x) for x in db.execute(sql,[*args,limit+1,offset])]
    more=len(rows)>limit
    rows=rows[:limit]
    numbers=[x["molregno"] for x in rows]
    syn={}
    atc={}
    if numbers:
        # A REST-compatible page can contain 1000 molecules. SQLite queries
        # use conservative batches of <=900 bound parameters: do not reject
        # an entire legitimate page or silently omit its synonyms/ATC codes.
        for start in range(0,len(numbers),900):
            placeholder,params_=_in(numbers[start:start+900])
            for r in db.execute("SELECT molregno,synonyms,syn_type FROM molecule_synonyms WHERE molregno IN "+placeholder,params_):
                syn.setdefault(r["molregno"],[]).append({
                    "molecule_synonym":r["synonyms"],"synonyms":r["synonyms"],"syn_type":r["syn_type"]})
            for r in db.execute("SELECT molregno,level5 FROM molecule_atc_classification WHERE molregno IN "+placeholder,params_):
                atc.setdefault(r["molregno"],[]).append(r["level5"])
    out=[]
    for r in rows:
        ident=r["chembl_id"]
        parent=r["parent_id"] or ident
        active=r["active_id"] or parent
        out.append({
            "molecule_chembl_id":ident,"pref_name":r["pref_name"],"max_phase":r["max_phase"],
            "molecule_type":r["molecule_type"],"withdrawn_flag":r["withdrawn_flag"],
            "availability_type":r["availability_type"],"veterinary":r["veterinary"],
            "first_approval":r["first_approval"],"oral":r["oral"],
            "parenteral":r["parenteral"],"topical":r["topical"],"prodrug":r["prodrug"],
            "molecule_hierarchy":{"molecule_chembl_id":ident,
                "parent_chembl_id":parent,"active_chembl_id":active},
            "molecule_synonyms":syn.get(r["molregno"],[]),
            "atc_classifications":sorted(set(atc.get(r["molregno"],[]))),
            "molecule_structures":{"canonical_smiles":r["canonical_smiles"]} if r["canonical_smiles"] else None,
            "molecule_properties":{"full_mwt":r["full_mwt"],"mw_freebase":r["mw_freebase"]}
        })
    return out,more

def _mechanisms(db,limit,offset):
    sql=(
        "SELECT dm.mec_id,dm.record_id,dm.site_id,dm.action_type,"
        "dm.mechanism_of_action,dm.direct_interaction,dm.disease_efficacy,"
        "dm.molecular_mechanism,dm.mechanism_comment,dm.selectivity_comment,"
        "dm.binding_site_comment,md.chembl_id AS molecule_chembl_id,"
        "parent.chembl_id AS parent_chembl_id,md.max_phase,"
        "td.chembl_id AS target_chembl_id "
        "FROM drug_mechanism dm "
        "LEFT JOIN compound_records cr ON dm.record_id=cr.record_id "
        "LEFT JOIN molecule_dictionary md ON md.molregno=COALESCE(dm.molregno,cr.molregno) "
        "LEFT JOIN molecule_hierarchy mh ON mh.molregno=md.molregno "
        "LEFT JOIN molecule_dictionary parent ON parent.molregno=mh.parent_molregno "
        "LEFT JOIN target_dictionary td ON td.tid=dm.tid "
        "ORDER BY dm.mec_id LIMIT ? OFFSET ?")
    rows=[dict(r) for r in db.execute(sql,(limit+1,offset))]
    more=len(rows)>limit
    out=[]
    for r in rows[:limit]:
        r["parent_molecule_chembl_id"]=r.pop("parent_chembl_id") or r["molecule_chembl_id"]
        r["mechanism_refs"]=[]
        r["variant_sequence"]=None
        out.append(r)
    return out,more

def _targets(db,params,limit,offset):
    cond=[]
    args=[]
    ids=[x for x in params.get("target_chembl_id__in","").split(",") if x]
    if ids:
        ph,v=_in(ids);cond.append("td.chembl_id IN "+ph);args.extend(v)
    accession=params.get("target_components__accession")
    if accession:
        cond.append("EXISTS (SELECT 1 FROM target_components tc2 JOIN component_sequences cs2 "
                    "ON cs2.component_id=tc2.component_id "
                    "WHERE tc2.tid=td.tid AND cs2.accession=?)")
        args.append(accession)
    sql=("SELECT td.tid,td.chembl_id,td.pref_name,td.organism,td.tax_id,"
         "td.target_type,td.species_group_flag FROM target_dictionary td")
    if cond:sql+=" WHERE "+" AND ".join(cond)
    sql+=" ORDER BY td.tid LIMIT ? OFFSET ?"
    rows=[dict(x) for x in db.execute(sql,[*args,limit+1,offset])]
    more=len(rows)>limit
    out=[]
    for r in rows[:limit]:
        comps=[dict(x) for x in db.execute(
            "SELECT cs.component_id,cs.accession,cs.component_type,"
            "cs.description AS component_description FROM target_components tc "
            "JOIN component_sequences cs ON cs.component_id=tc.component_id "
            "WHERE tc.tid=? ORDER BY cs.component_id",(r["tid"],))]
        out.append({
          "target_chembl_id":r["chembl_id"],"pref_name":r["pref_name"],
          "organism":r["organism"],"tax_id":r["tax_id"],
          "target_type":r["target_type"],"species_group_flag":r["species_group_flag"],
          "target_components":comps,"cross_references":[]})
    return out,more

def _activities(db,params,limit,offset):
    cond=[]
    args=[]
    tid=params.get("target_chembl_id")
    if tid:
        cond.append("td.chembl_id=?");args.append(tid)
    types=[x for x in params.get("standard_type__in","").split(",") if x]
    if types:
        ph,v=_in(types);cond.append("act.standard_type IN "+ph);args.extend(v)
    if "standard_units" in params:
        cond.append("act.standard_units=?");args.append(params["standard_units"])
    if "standard_value__lte" in params:
        cond.append("act.standard_value<=?");args.append(float(params["standard_value__lte"]))
    where=" WHERE "+" AND ".join(cond) if cond else ""
    if not cond:raise RuntimeError("V31_SQLITE_UNBOUNDED_ACTIVITY_SCAN_REFUSED")
    sql=(
        "SELECT act.activity_id,ass.chembl_id AS assay_chembl_id,"
        "md.chembl_id AS molecule_chembl_id,"
        "parent.chembl_id AS parent_molecule_chembl_id,"
        "td.chembl_id AS target_chembl_id,"
        "act.standard_type,act.standard_value,act.standard_units,"
        "act.standard_relation,act.data_validity_comment,"
        "cs.canonical_smiles,docs.chembl_id AS document_chembl_id "
        "FROM activities act "
        "JOIN assays ass ON ass.assay_id=act.assay_id "
        "LEFT JOIN target_dictionary td ON td.tid=ass.tid "
        "LEFT JOIN molecule_dictionary md ON md.molregno=act.molregno "
        "LEFT JOIN molecule_hierarchy mh ON mh.molregno=md.molregno "
        "LEFT JOIN molecule_dictionary parent ON parent.molregno=mh.parent_molregno "
        "LEFT JOIN compound_structures cs ON cs.molregno=md.molregno "
        "LEFT JOIN docs ON docs.doc_id=act.doc_id "
        +where+" ORDER BY act.activity_id LIMIT ? OFFSET ?")
    rows=[dict(x) for x in db.execute(sql,[*args,limit+1,offset])]
    more=len(rows)>limit
    for r in rows[:limit]:
        if not r.get("parent_molecule_chembl_id"):
            r["parent_molecule_chembl_id"]=r.get("molecule_chembl_id")
    return rows[:limit],more

def _assays(db,params,limit,offset):
    ids=[x for x in params.get("assay_chembl_id__in","").split(",") if x]
    if not ids:raise RuntimeError("V31_SQLITE_UNBOUNDED_ASSAY_SCAN_REFUSED")
    ph,v=_in(ids)
    sql=("SELECT chembl_id AS assay_chembl_id,description,confidence_score,assay_type "
         "FROM assays WHERE chembl_id IN "+ph+
         " ORDER BY assay_id LIMIT ? OFFSET ?")
    rows=[dict(x) for x in db.execute(sql,[*v,limit+1,offset])]
    return rows[:limit],len(rows)>limit

def local_json(url):
    parsed=urllib.parse.urlsplit(url)
    if parsed.scheme!="https" or parsed.hostname!="www.ebi.ac.uk" or not parsed.path.startswith("/chembl/api/data/"):
        raise RuntimeError("V31_CHEMBL_SQLITE_URL_NOT_OFFICIAL")
    path=parsed.path.split("/chembl/api/data/",1)[1]
    parts=path.split("/")
    resource=parts[0].split(".")[0]
    if resource not in ("molecule","mechanism","target","activity","assay"):
        raise RuntimeError("V31_CHEMBL_SQLITE_RESOURCE_UNSUPPORTED:"+resource)
    params={k:v[-1] for k,v in urllib.parse.parse_qs(parsed.query,keep_blank_values=True).items()}
    limit=int(params.get("limit","1000"))
    offset=int(params.get("offset","0"))
    if not (1<=limit<=1000 and 0<=offset<=10000000):
        raise ValueError("V31_CHEMBL_SQLITE_PAGINATION_OUT_OF_BOUNDS")
    with contextlib.closing(connection()) as db:
        if resource=="molecule":
            if len(parts)==2:
                if not parts[1].endswith(".json"):
                    raise ValueError("V31_SQLITE_MOLECULE_DETAIL_FORMAT_INVALID")
                molecule=parts[1][:-5]
                params={"molecule_chembl_id__in":molecule}
                rows,_=_molecules(db,params,1,0)
                if len(rows)!=1:raise LookupError("V31_SQLITE_MOLECULE_NOT_FOUND")
                return rows[0]
            rows,more=_molecules(db,params,limit,offset)
        elif resource=="mechanism": rows,more=_mechanisms(db,limit,offset)
        elif resource=="target": rows,more=_targets(db,params,limit,offset)
        elif resource=="activity": rows,more=_activities(db,params,limit,offset)
        elif resource=="assay": rows,more=_assays(db,params,limit,offset)
    key={"molecule":"molecules","mechanism":"mechanisms","target":"targets",
         "activity":"activities","assay":"assays"}[resource]
    return {key:rows,"page_meta":_metadata(resource+".json",params,offset,limit,more)}
