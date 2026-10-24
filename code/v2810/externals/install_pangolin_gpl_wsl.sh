#!/usr/bin/env bash
set -euo pipefail
BASE="$HOME/rdmva_tools"
BOOT="$BASE/.bootstrap_venv"
RUNTIME="$BASE/pangolin_py312"
mkdir -p "$BASE"
if [[ ! -x "$BOOT/bin/uv" ]]; then
  python3 -m venv "$BOOT"
  "$BOOT/bin/python" -m pip install --no-cache-dir uv
fi
UV="$BOOT/bin/uv"
"$UV" python install 3.12
if [[ ! -x "$RUNTIME/bin/python" ]]; then
  "$UV" venv --python 3.12 "$RUNTIME"
fi
"$UV" pip install --python "$RUNTIME/bin/python" torch numpy pandas gffutils pyfastx biopython PyVCF3 setuptools==80.9.0
"$UV" pip install --python "$RUNTIME/bin/python" --no-deps -e "$BASE/Pangolin"
"$RUNTIME/bin/python" -c 'import torch,vcf,pangolin,gffutils,pyfastx,pandas; print("PANGOLIN_RUNTIME_READY"); print("TORCH",torch.__version__)'
echo PANGOLIN_RUNTIME_SETUP_COMPLETE
