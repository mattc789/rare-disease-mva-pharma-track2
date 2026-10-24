"""V3.1 accelerated ChEMBL activity acquisition.

Semantics are intentionally identical to v2810.screen.route_activities; only
independent target-page retrieval is parallelised. Existing per-target caches
are reused, so interrupted runs resume without repeated downloads.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import concurrent.futures
import urllib.parse
import hashlib,json,os,re,shutil
from v2810 import screen
from v2810.core import load,save

SCHEMA="RDMVA_FAST_SCREEN_V31_1"


def target_index_parallel(cache_dir,routes,max_workers=8):
    cache=Path(cache_dir)/"target_index.json"
    if cache.exists():return load(cache)
    idx={}
    for t in screen.load(screen.CHEMBL_TARGETS):
        if t.get("organism")!="Homo sapiens":continue
        comps=t.get("target_components") or []
        accs={x.get("accession") for x in comps if x.get("accession")}
        if len(accs)==1:
            idx.setdefault(next(iter(accs)),[]).append(t["target_chembl_id"])
    needed={t.get("accession") for r in routes for t in r.get("targets",[]) if t.get("accession")}
    missing=sorted(acc for acc in needed if acc not in idx)

    def resolve(acc):
        u="https://www.ebi.ac.uk/chembl/api/data/target.json?"+urllib.parse.urlencode(
            {"target_components__accession":acc,"limit":100})
        j=screen.fetch_json(u,20,8_000_000)
        ids=[]
        for t in j.get("targets",[]):
            comps={x.get("accession") for x in t.get("target_components",[]) if x.get("accession")}
            if t.get("tax_id")==9606 and t.get("target_type")=="SINGLE PROTEIN" and comps=={acc}:
                ids.append(t["target_chembl_id"])
        return acc,sorted(set(ids))

    if missing:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs={ex.submit(resolve,acc):acc for acc in missing}
            for fut in concurrent.futures.as_completed(futs):
                acc=futs[fut]
                try:
                    _,ids=fut.result()
                except Exception as exc:
                    raise RuntimeError("V31_CHEMBL_TARGET_UNAVAILABLE:"+str(acc)) from exc
                if ids:idx[acc]=ids
    for acc in list(idx):idx[acc]=sorted(set(idx[acc]))
    save(cache,idx)
    save(Path(cache_dir)/"V31_TARGET_INDEX_FETCH_AUDIT.json",{
      "schema":"RDMVA_FAST_TARGET_INDEX_V31_1",
      "needed_accessions":len(needed),
      "resolved_accessions":sum(acc in idx for acc in needed),
      "api_missing_accessions_attempted":len(missing),
      "max_workers":max_workers})
    return idx

def route_activities_parallel(routes,target_idx,cache_dir,max_workers=8):
    cache=Path(cache_dir)/"route_activities.json"
    if cache.exists():return load(cache)
    tid_routes=defaultdict(list)
    for r in routes:
        for t in r.get("targets",[]):
            for tid in target_idx.get(t.get("accession"),[]):
                tid_routes[tid].append(r["route_id"])

    def fetch_target(tid):
        cp=Path(cache_dir)/f"act_{tid}.json"
        if cp.exists():
            rr=load(cp)
            return tid,rr,True
        rr=screen.chembl_paged(
          "activity",{
            "target_chembl_id":tid,
            "standard_type__in":"IC50,EC50,AC50,Ki,Kd",
            "standard_units":"nM","standard_value__lte":"10000"},
          "activities")
        save(cp,rr)
        return tid,rr,False

    raw={}
    # Moderate concurrency: enough to remove serial latency without abusive
    # request fan-out to the public ChEMBL service.
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs={ex.submit(fetch_target,tid):tid for tid in sorted(tid_routes)}
        for fut in concurrent.futures.as_completed(futs):
            tid=futs[fut]
            try:
                _,rr,_=fut.result();raw[tid]=rr
            except Exception as exc:
                raise RuntimeError("V31_CHEMBL_ACTIVITY_UNAVAILABLE:"+str(tid)) from exc

    acts=[];assays={}
    for tid in sorted(raw):
        for a in raw[tid]:
            try:v=float(a.get("standard_value"))
            except:continue
            if a.get("standard_relation")!="=" or a.get("standard_units")!="nM" or not 0<v<=10000:
                continue
            acts.append({k:a.get(k) for k in [
              "activity_id","assay_chembl_id","molecule_chembl_id",
              "parent_molecule_chembl_id","target_chembl_id","standard_type",
              "standard_value","canonical_smiles","document_chembl_id",
              "data_validity_comment"]})
    aids=sorted({a["assay_chembl_id"] for a in acts if a.get("assay_chembl_id")})
    batches=[aids[off:off+100] for off in range(0,len(aids),100)]
    def fetch_assay_batch(ids):
        return screen.chembl_paged(
            "assay",{"assay_chembl_id__in":",".join(ids)},"assays",5000)
    if batches:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs=[ex.submit(fetch_assay_batch,ids) for ids in batches]
            for fut in concurrent.futures.as_completed(futs):
                try:
                    rr=fut.result()
                except Exception as exc:
                    raise RuntimeError("V31_CHEMBL_ASSAY_UNAVAILABLE") from exc
                for a in rr:
                    assays[a["assay_chembl_id"]]={
                      "description":a.get("description"),
                      "confidence_score":a.get("confidence_score"),
                      "assay_type":a.get("assay_type")}

    route_by={r["route_id"]:r for r in routes};rows=[]
    for a in acts:
        ass=assays.get(a.get("assay_chembl_id"),{})
        if ass.get("confidence_score") not in (9,"9") or ass.get("assay_type") not in ("B","F"):
            continue
        lo=(ass.get("description") or "").lower()
        typ=a.get("standard_type")
        if typ=="IC50" and ("agonist" not in lo and "activate" not in lo):
            direction=-1
        elif typ in ("EC50","AC50") and any(k in lo for k in ("agonist","stimulat","activat","induc")):
            direction=1
        elif typ in ("EC50","AC50") and any(k in lo for k in ("inhibit","antagon","repress","suppress")):
            direction=-1
        else:
            direction=0
        direct_functional=(
          ass.get("assay_type")=="F" and direction in (-1,1) and
          typ in ("IC50","EC50","AC50") and
          any(token in lo for token in
              ("activat","agonist","inhibit","antagon","stimulat","suppress")))
        for rid in tid_routes[a["target_chembl_id"]]:
            desired=1 if route_by[rid]["desired_action"]=="increase" else -1 if route_by[rid]["desired_action"]=="decrease" else 0
            role="compatible" if direction and direction==desired else "opposite" if direction and direction!=desired else "unknown"
            rows.append({**a,"route_id":rid,"value_nM":float(a["standard_value"]),
              "direction":direction,"role":role,
              "assay_description":ass.get("description"),
              "assay_type":ass.get("assay_type"),
              "direct_target_assay_verified":bool(direct_functional and role=="compatible"),
              "assay_target_confidence":ass.get("confidence_score")})
    save(cache,rows)
    save(Path(cache_dir)/"V31_ACTIVITY_FETCH_AUDIT.json",{
      "schema":SCHEMA,"target_ids":len(tid_routes),"activity_rows_raw":len(acts),
      "qualified_route_activity_rows":len(rows),"max_workers":max_workers,"assay_batches":len(batches),
      "semantic_equivalence_target":"v2810.screen.route_activities"})
    return rows

def _file_sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as source:
        for block in iter(lambda:source.read(4*1024*1024),b""):h.update(block)
    return h.hexdigest()

def restore_public_chembl_cache(previous_run,new_run):
    # Public ChEMBL cache, and only when the approved-drug universe was
    # independently rebuilt from byte-identical source snapshots.
    source=Path(previous_run)
    target=Path(new_run)
    audit={"schema":"RDMVA_VERIFIED_CH_EMBL_CACHE_REUSE_V31_1",
           "source_run":str(source),"reused_public_target_files":0,
           "routes_identical":False,"fresh_chembl_downloads_still_allowed":True,
           "patient_case_or_genomic_evidence_reused":False}
    old_master=source/"master/summary.json"
    new_master=target/"master/summary.json"
    if not old_master.is_file() or not new_master.is_file():
        raise RuntimeError("V31_REUSE_MASTER_UNIVERSE_MISSING")
    x=load(old_master);y=load(new_master)
    for key in ("entities","products","baseline_us_entities","dataset_profile","source_hashes"):
        if x.get(key)!=y.get(key):
            raise RuntimeError("V31_REUSE_UNIVERSE_SNAPSHOT_MISMATCH:"+key)
    old_ledger=source/"master/entities.jsonl"
    new_ledger=target/"master/entities.jsonl"
    if not old_ledger.is_file() or not new_ledger.is_file() or _file_sha(old_ledger)!=_file_sha(new_ledger):
        raise RuntimeError("V31_REUSE_ENTITY_LEDGER_MISMATCH")
    old_cache=source/"screen/sources"
    new_cache=target/"screen/sources"
    if not old_cache.is_dir():return audit
    new_cache.mkdir(parents=True,exist_ok=True)
    def transfer(name):
        f=old_cache/name
        if not f.is_file():return False
        if f.stat().st_size>150*1024*1024:
            raise RuntimeError("V31_REUSE_UNEXPECTED_CACHE_FILE_SIZE:"+name)
        dest=new_cache/name
        if dest.exists():raise RuntimeError("V31_REUSE_DESTINATION_EXISTS:"+name)
        shutil.copy2(f,dest)
        if _file_sha(f)!=_file_sha(dest):
            raise RuntimeError("V31_REUSE_CACHE_SHA256_MISMATCH:"+name)
        return True
    for f in sorted(old_cache.glob("act_CHEMBL*.json")):
        m=re.fullmatch(r"act_(CHEMBL[0-9]+)\.json",f.name)
        if not m:continue
        rows=load(f)
        if (not isinstance(rows,list) or
            any(not isinstance(x,dict) or x.get("target_chembl_id")!=m.group(1) for x in rows)):
            raise RuntimeError("V31_REUSE_PUBLIC_TARGET_PAYLOAD_INVALID:"+f.name)
        if transfer(f.name):audit["reused_public_target_files"]+=1
    # Structures are universe-dependent; identical approved-entity hashes
    # and source snapshots allow reuse of public structure responses.
    structures=old_cache/"structures.json"
    if structures.is_file():
        entries=load(structures)
        if not isinstance(entries,dict):
            raise RuntimeError("V31_REUSE_STRUCTURES_FORMAT_INVALID")
        audit["structures_reused"]=transfer("structures.json")
    # Target-index and route-activities are route-specific. Reuse them ONLY
    # if both unmodified full route JSON graphs are exactly identical.
    old_routes=source/"NATIVE_ROUTES_V31.json"
    new_routes=target/"NATIVE_ROUTES_V31.json"
    if old_routes.is_file() and new_routes.is_file():
        x=load(old_routes);y=load(new_routes)
        audit["routes_identical"]=x.get("routes")==y.get("routes")
    if audit["routes_identical"]:
        for name in ("target_index.json","route_activities.json"):
            if (old_cache/name).is_file():
                obj=load(old_cache/name)
                if not isinstance(obj,(dict,list)):
                    raise RuntimeError("V31_REUSE_ROUTE_CACHE_FORMAT_INVALID:"+name)
                audit[name+"_reused"]=transfer(name)
    audit["verified_master_sha256"]=_file_sha(new_master)
    audit["verified_entities_sha256"]=_file_sha(new_ledger)
    audit["cache_reuse_not_fresh_external_retrieval"]=True
    save(target/"V31_CHEMBL_CACHE_REUSE_AUDIT.json",audit)
    print("V31_VERIFIED_PUBLIC_CHEMBL_CACHES",
          audit["reused_public_target_files"],
          "ROUTE_MATCH",audit["routes_identical"],flush=True)
    return audit

def universal_screen(entities,routes,outdir,max_workers=8):
    original_activity=screen.route_activities
    original_target=screen.target_index
    def _fast_activities(routes_,target_idx_,cache_dir_):
        return route_activities_parallel(routes_,target_idx_,cache_dir_,max_workers=max_workers)
    def _fast_targets(cache_dir_,routes_):
        return target_index_parallel(cache_dir_,routes_,max_workers=max_workers)
    screen.route_activities=_fast_activities
    screen.target_index=_fast_targets
    try:
        print("V31_CH_EMBL_SCREEN_START",len(entities),"ENTITIES",len(routes),"ROUTES",flush=True)
        return screen.universal_screen(entities,routes,outdir)
    finally:
        screen.route_activities=original_activity
        screen.target_index=original_target
