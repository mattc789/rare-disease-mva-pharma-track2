#!/usr/bin/env python3
"""Download and verify the pinned 2026-10-10 global regulatory snapshot.

Rare Disease MVA Pharma — V2.8.3
Developed by Matteo Caprani | ExaMedix Technologies | www.examedix.ch

No patient data are read, shared or downloaded. The V2.8.2 snapshot is untouched.
"""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,os,tempfile,urllib.request

ROOT=Path(__file__).resolve().parent
DIR=ROOT/'public_inputs/third_party/regulatory_v283'
SPEC=(ROOT/'public_inputs/current_v31/regulatory_snapshot_v283_current.json'
  if os.environ.get('RDMVA_CURRENT_PUBLIC_DATA')=='1'
  else ROOT/'config/regulatory_snapshot_v283.json')
AGENT='RareDiseaseMVA-V283-RegulatorySnapshot/1.0'

def checked(path,details):
 if not path.is_file():return False
 if path.stat().st_size!=details['bytes']:return False
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()==details['sha256']

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--verify-only',action='store_true')
 args=p.parse_args()
 m=json.loads(SPEC.read_text(encoding='utf-8'))
 if m['schema']!='V283_REGULATORY_SCHEMA_1':raise RuntimeError('regulatory schema mismatch')
 missing=[]
 for fname,entry in sorted(m['files'].items()):
  dest=DIR/fname
  if not checked(dest,entry) and not args.verify_only:
   dest.parent.mkdir(parents=True,exist_ok=True)
   # Keep a stable partial filename for checkpoint-based resume.
   tmp=dest.with_name(dest.name+'.rdmva-download')
   try:
    # Hash-pinned regulatory public datasets share the accelerated source
    # transport. An upstream snapshot change remains a hard failure.
    import sys
    code=str(ROOT/'code')
    if code not in sys.path:sys.path.insert(0,code)
    from v2810.scientific_download_v31 import download_public_data
    download_public_data(entry['url'],tmp,
                         expected_bytes=int(entry['bytes']),
                         max_bytes=int(entry['bytes'])+1048576,
                         expected_sha256=entry['sha256'],label=fname)
    if not checked(tmp,entry):
     raise RuntimeError('Pinned upstream snapshot changed for '+fname+'; do not silently replace frozen inputs')
    os.replace(tmp,dest)
   finally:
    if tmp.exists():tmp.unlink()
  if checked(dest,entry):print('V283_REGULATORY_SHA_OK',fname,flush=True)
  else:missing.append(fname)
 if missing:raise RuntimeError('Missing / changed pinned regulatory inputs: '+', '.join(missing))
 print('V283_REGULATORY_SNAPSHOT_OK',len(m['files']),'sources',flush=True)
if __name__=='__main__':main()
