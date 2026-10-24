#!/usr/bin/env python3
"""GPL-3.0 compatibility adaptation for independent Pangolin source tree.

PyVCF3 added a seventh _Info(..., type_code) argument; original Pangolin
uses legacy six-argument _Info. The adapter retains original GPL source and
upstream copyright; it does not alter model weights or mix Pangolin into
the application's CC BY 4.0 code archive.
"""
import pathlib,sys
if len(sys.argv)!=2:raise SystemExit("Usage: patch_pangolin_pyvcf3.py /path/to/Pangolin")
repo=pathlib.Path(sys.argv[1]).resolve()
file=repo/"pangolin/pangolin.py"
license=repo/"LICENSE"
if not file.is_file() or not license.is_file():
    raise RuntimeError("PANGOLIN_SOURCE_MISSING")
text=license.read_text(encoding="utf-8")
if "GNU GENERAL PUBLIC LICENSE" not in text or "Version 3" not in text:
    raise RuntimeError("PANGOLIN_LICENSE_NOT_GPL3")
s=file.read_text(encoding="utf-8")
old='"Format: gene|pos:score_change|pos:score_change|warnings,...",\'.\',\'.\')'
new='"Format: gene|pos:score_change|pos:score_change|warnings,...",\'.\',\'.\',vcf.parser._encode_type("String"))'
if s.count(new)==1:
    print("PANGOLIN_PYVCF3_SHIM_ALREADY_APPLIED")
elif s.count(old)==1:
    file.write_text(s.replace(old,new,1),encoding="utf-8")
    print("PANGOLIN_PYVCF3_SHIM_APPLIED")
else:
    raise RuntimeError("PANGOLIN_UPSTREAM_LAYOUT_UNRECOGNIZED_REFUSE_MODIFY")
