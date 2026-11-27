#!/bin/bash
#SBATCH --job-name=evo2_icl
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=bsse-gpu
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --gres=gpu:GH100:4
#
# The headline run: Evo 2 ICL swap control on bs-slurm H100s.
#
# Hopper matters here. The model card puts FP8 input projections on by default
# and warns that disabling them changes outputs noticeably, so this is the only
# place the numbers are trustworthy. The Euler smoke test proves the plumbing;
# this produces the result.
#
# --time=02:00:00 is deliberate. 67 jobs are queued for GH100 and the ones
# running hold GPUs for the full 3-day limit, so a short job can backfill into
# the gaps between them. Recent history: long requests waited ~3.4 days, one
# short job started in 4.4 hours. Do not pad this.
#
# No internet on this cluster, on login or compute nodes. Everything -- venv,
# checkpoint, prompts -- must already be staged from Euler.
#
# Submit from ~/ecoli_tfs_evo2:
#   sbatch scripts/basel_score_evo2.sbatch

set -euo pipefail

BASE="${SLURM_SUBMIT_DIR:-$PWD}"
VENV="$BASE/venv"
NGPU="${NGPU:-4}"

export HF_HOME="$BASE/hf_cache"
export TORCH_HOME="$BASE/hf_cache/torch"
export TMPDIR="${TMPDIR:-$BASE/tmp}"
export HF_HUB_OFFLINE=1      # not optional here: there is no outbound network
export TRANSFORMERS_OFFLINE=1
mkdir -p "$BASE/logs" "$TMPDIR"

echo "host:  $(hostname)"
echo "gpus:  $NGPU"
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader

source "$VENV/bin/activate"
python -c "import torch; print('torch', torch.__version__, '| cuda', torch.cuda.is_available(), '| devices', torch.cuda.device_count())"

# One worker per GPU. Prompts are split round-robin by --shard, and each worker
# appends to its own file so there is no write contention; the shards are
# concatenated afterwards. Each worker resumes from its own file if requeued.
pids=()
for i in $(seq 0 $((NGPU - 1))); do
    CUDA_VISIBLE_DEVICES=$i python "$BASE/scripts/score_evo2_icl.py" \
        --indir "$BASE/data/evo2_icl" \
        --out "$BASE/data/evo2_icl/scores_shard${i}.tsv" \
        --model "${EVO2_MODEL:-Taykhoom/Evo2-7B-1M}" \
        --batch-size "${BATCH_SIZE:-4}" \
        --shard "$i" --num-shards "$NGPU" \
        > "$BASE/logs/shard${i}.log" 2>&1 &
    pids+=($!)
done

fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done

OUT="$BASE/data/evo2_icl/scores.tsv"
head -1 "$BASE/data/evo2_icl/scores_shard0.tsv" > "$OUT"
for i in $(seq 0 $((NGPU - 1))); do
    tail -n +2 "$BASE/data/evo2_icl/scores_shard${i}.tsv" >> "$OUT"
done
echo "merged $(($(wc -l < "$OUT") - 1)) rows -> $OUT"

if [ "$fail" -ne 0 ]; then
    echo "WARNING: at least one shard exited non-zero; check logs/shard*.log" >&2
    exit 1
fi

python "$BASE/scripts/analyze_evo2_icl.py" \
    --scores "$OUT" \
    --prompts "$BASE/data/evo2_icl/prompts.jsonl" \
    --out "$BASE/data/evo2_icl/results.tsv"
