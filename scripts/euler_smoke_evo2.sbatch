#!/bin/bash
#SBATCH --job-name=evo2_smoke
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --gpus=1
# No gpumem constraint here: it selects on memory alone and happily lands on a
# TITAN RTX (sm_7.5, Turing), which has no native bfloat16. Pass a GPU type on
# the command line instead, e.g. --gpus=rtx_4090:1 or --gpus=rtx_3090:1.
#
# Smoke test on Euler before staging anything to bs-slurm.
#
# The point is NOT the numbers -- Euler has no Hopper, so FP8 input projections
# are off and the model card warns outputs then differ noticeably from the
# reference. The point is to prove the pipeline runs end to end: the checkpoint
# loads, the tokenizer really is one token per nucleotide, the candidate slice
# lines up, and a 12 kb prompt fits in memory. Discovering a shape bug here
# costs an hour; discovering it after a 3-day bs-slurm queue costs three days.
#
# Submit from $SCRATCH/ecoli_tfs_evo2:
#   sbatch scripts/euler_smoke_evo2.sbatch

set -euo pipefail

BASE="${SLURM_SUBMIT_DIR:-$PWD}"
VENV="$BASE/venv"
export HF_HOME="$BASE/hf_cache"
export TORCH_HOME="$BASE/hf_cache/torch"
export TMPDIR="$BASE/tmp"
export HF_HUB_OFFLINE=1          # weights are already cached; never dial out
mkdir -p "$BASE/logs" "$TMPDIR"

echo "host: $(hostname)"
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader

source "$VENV/bin/activate"
python -c "import torch, transformers; print('torch', torch.__version__, '| transformers', transformers.__version__, '| cuda', torch.cuda.is_available())"

# 2 prompts (one real/swap pair) x 12 candidates: enough to prove alignment and
# to confirm the pair differs only where it should.
#
# Run python directly, NOT under srun. On a shared Euler GPU node srun tries to
# create a job step and fails with "Requested nodes are busy", retrying until
# the walltime expires. A single-task batch job needs no step.
python "$BASE/scripts/score_evo2_icl.py" \
    --indir "$BASE/data/evo2_icl" \
    --model "${EVO2_MODEL:-Taykhoom/Evo2-7B-1M}" \
    --batch-size "${BATCH_SIZE:-2}" \
    --out "$BASE/data/evo2_icl/${SMOKE_OUT:-scores_smoke.tsv}" \
    --limit-prompts 2 \
    --limit-candidates 12 \
    --no-fp8

echo
echo "=== smoke output ==="
head -5 "$BASE/data/evo2_icl/${SMOKE_OUT:-scores_smoke.tsv}"
wc -l "$BASE/data/evo2_icl/${SMOKE_OUT:-scores_smoke.tsv}"
