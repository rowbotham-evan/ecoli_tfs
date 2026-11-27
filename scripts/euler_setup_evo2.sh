#!/usr/bin/env bash
# Build the Evo 2 inference environment on Euler, entirely inside $SCRATCH.
#
# Euler is the staging and de-risking host: it has outbound internet, which
# bs-slurm does not, and spare GPUs for a smoke test. The faithful run happens
# later on bs-slurm H100s, using artefacts copied from here.
#
# HuggingFace transformers, not vLLM: every public Evo 2 port declares
# model_type "evo2" with trust_remote_code custom code, and vLLM only runs
# architectures in its own registry. That costs us prefix caching, which is why
# the scorer subsamples negatives instead.
#
# Follows the conventions of ecoli_tfs_dbd/euler_setup.sh: home is at its 50 GB
# quota, so venv, pip cache and HF cache all live in $SCRATCH, and weights are
# pre-downloaded on the login node because compute nodes have no outbound
# internet.
#
# Run as:  ssh eth "bash -lc '/cluster/scratch/\$USER/ecoli_tfs_evo2/euler_setup_evo2.sh'"
set -euo pipefail

module load stack/2024-06 python/3.11.6 >/dev/null 2>&1 || {
  echo "could not load python module; run this under: bash -lc" >&2
  exit 1
}
echo ">>> python: $(python3 -V) at $(command -v python3)"

BASE="${SCRATCH:?SCRATCH is not set}/ecoli_tfs_evo2"
VENV="$BASE/venv"
MODEL="${EVO2_MODEL:-Taykhoom/Evo2-7B-1M}"

mkdir -p "$BASE"/{data,logs,hf_cache,pip_cache,tmp,wheels}
export PIP_CACHE_DIR="$BASE/pip_cache"
export HF_HOME="$BASE/hf_cache"
export TORCH_HOME="$BASE/hf_cache/torch"
export TMPDIR="$BASE/tmp"

# An earlier attempt left a half-installed vLLM venv here; vLLM is not used.
if [ -d "$VENV" ] && [ ! -f "$BASE/.venv_is_transformers" ]; then
  echo ">>> removing the earlier partial venv (half-installed vLLM, unusable)"
  rm -rf "$VENV"
fi

if [ ! -x "$VENV/bin/python" ]; then
  echo ">>> creating venv at $VENV"
  python3 -m venv "$VENV"
fi
PY="$VENV/bin/python"
"$PY" -m pip install --quiet --upgrade pip wheel setuptools

echo
echo ">>> STAGE 1: torch + transformers"
# No flash-attn and no Transformer Engine here. The port falls back to sdpa,
# and FP8 input projections are disabled at load time on non-Hopper hardware.
"$PY" -m pip install --quiet torch transformers accelerate safetensors huggingface_hub numpy
"$PY" - <<'PY'
import torch, transformers
print("torch       ", torch.__version__, "| cuda build", torch.version.cuda)
print("transformers", transformers.__version__)
PY
touch "$BASE/.venv_is_transformers"

echo
echo ">>> STAGE 2: checkpoint $MODEL (login node has internet; compute nodes do not)"
EVO2_MODEL="$MODEL" "$PY" - <<'PY'
import os, pathlib, sys
from huggingface_hub import snapshot_download
repo = os.environ["EVO2_MODEL"]
try:
    p = snapshot_download(repo_id=repo)
except Exception as e:
    print("download FAILED:", type(e).__name__, e, file=sys.stderr)
    sys.exit(2)
total = sum(f.stat().st_size for f in pathlib.Path(p).rglob("*") if f.is_file())
print("cached at:", p)
print("size: %.1f GB" % (total / 1e9))
print("files:", sorted(f.name for f in pathlib.Path(p).iterdir())[:20])
PY

echo
echo ">>> STAGE 3: wheels for offline transfer to bs-slurm (python 3.12, linux x86_64)"
# bs-slurm has python 3.12 and no internet, so its wheels are collected here.
"$PY" -m pip download --quiet --dest "$BASE/wheels" \
    --python-version 312 --only-binary=:all: \
    --platform manylinux2014_x86_64 \
    torch transformers accelerate safetensors huggingface_hub numpy 2>&1 | tail -3 || \
    echo "NOTE: wheel pre-download for 3.12 had problems; revisit before staging"
echo "wheels: $(ls "$BASE/wheels" 2>/dev/null | wc -l) files, $(du -sh "$BASE/wheels" 2>/dev/null | cut -f1)"

echo
echo ">>> setup complete"
echo "    venv:    $VENV"
echo "    HF_HOME: $HF_HOME"
echo "    model:   $MODEL"
