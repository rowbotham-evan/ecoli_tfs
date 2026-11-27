#!/usr/bin/env bash
# Stage the Evo 2 run from Euler to bs-slurm. RUNS ON EULER.
#
# bs-slurm has no outbound network on login or compute nodes, so everything --
# wheels, checkpoint, prompts, scripts -- has to arrive over ssh from here.
# Euler can reach euler->basel on port 22.
#
# It also has no pip and no ensurepip, so a venv cannot be built there the
# normal way. get-pip.py is fetched here (it embeds pip) and the venv is
# created on basel with --without-pip, then bootstrapped offline.
#
#   ssh eth "bash -lc '/cluster/scratch/\$USER/ecoli_tfs_evo2/scripts/stage_euler_to_basel.sh'"
#
# Re-runnable: rsync skips what is already there.
set -euo pipefail

BASE="${SCRATCH:?}/ecoli_tfs_evo2"
REMOTE_HOST="${REMOTE_HOST:-bs-slurm}"
REMOTE_USER="${REMOTE_USER:-$USER}"
REMOTE_BASE="${REMOTE_BASE:-ecoli_tfs_evo2}"   # relative to basel's $HOME
SSH="ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new"
DEST="$REMOTE_USER@$REMOTE_HOST"

echo ">>> checking reachability"
$SSH "$DEST" 'echo "connected to $(hostname)"'

echo
echo ">>> get-pip.py (basel has no ensurepip)"
if [ ! -f "$BASE/wheels/get-pip.py" ]; then
    curl -sSL -o "$BASE/wheels/get-pip.py" https://bootstrap.pypa.io/get-pip.py
fi
echo "    $(wc -c < "$BASE/wheels/get-pip.py") bytes"

# The model snapshot is a tree of symlinks into hf_cache/hub/.../blobs, so it
# must be copied with -L to send real files rather than dangling links.
SNAP="$(find "$BASE/hf_cache/hub" -maxdepth 3 -type d -name 'snapshots' -print -quit)"
MODEL_DIR="$(find "$SNAP" -maxdepth 1 -mindepth 1 -type d -print -quit)"
echo
echo ">>> model snapshot: $MODEL_DIR"
du -shL "$MODEL_DIR"

echo
echo ">>> creating remote layout"
$SSH "$DEST" "mkdir -p '$REMOTE_BASE'/{scripts,data/evo2_icl,wheels,logs,tmp,model}"

echo
echo ">>> 1/4 scripts + prompts (small)"
rsync -az --info=stats2 -e "$SSH" \
    "$BASE/scripts/" "$DEST:$REMOTE_BASE/scripts/"
rsync -az --info=stats2 -e "$SSH" \
    "$BASE/data/evo2_icl/prompts.jsonl" "$BASE/data/evo2_icl/candidates.tsv" \
    "$DEST:$REMOTE_BASE/data/evo2_icl/"

echo
echo ">>> 2/4 wheels (~2.8 GB)"
rsync -az --partial --info=progress2 -e "$SSH" \
    "$BASE/wheels/" "$DEST:$REMOTE_BASE/wheels/"

echo
echo ">>> 3/4 checkpoint (~13.2 GB, dereferencing symlinks)"
rsync -azL --partial --info=progress2 -e "$SSH" \
    "$MODEL_DIR/" "$DEST:$REMOTE_BASE/model/"

echo
echo ">>> 4/4 building the venv on basel, offline"
$SSH "$DEST" "bash -s" <<REMOTE
set -euo pipefail
cd "\$HOME/$REMOTE_BASE"
if [ ! -x venv/bin/python ]; then
    python3 -m venv --without-pip venv
    ./venv/bin/python wheels/get-pip.py --no-index --find-links wheels
fi
./venv/bin/python -m pip install --no-index --find-links wheels \
    torch transformers accelerate safetensors huggingface_hub numpy
./venv/bin/python - <<'PY'
import torch, transformers
print("basel torch       ", torch.__version__)
print("basel transformers", transformers.__version__)
PY
df -h "\$HOME" | tail -1
REMOTE

echo
echo ">>> staged. On basel, submit with:"
echo "    cd ~/$REMOTE_BASE && EVO2_MODEL=\$HOME/$REMOTE_BASE/model sbatch scripts/basel_score_evo2.sbatch"
echo "    (EVO2_MODEL is a local path: there is no hub access there)"
