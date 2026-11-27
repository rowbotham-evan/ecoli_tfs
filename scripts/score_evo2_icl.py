#!/usr/bin/env python3
"""Score Evo 2 ICL candidate windows with HuggingFace transformers.

For every (prompt, candidate) pair, records the summed log-probability of the
candidate's nucleotides conditioned on the prompt.

Why not vLLM: every public Evo 2 port declares model_type "evo2" with
trust_remote_code custom modelling code, and vLLM only serves architectures in
its own registry. So there is no prefix caching, and each candidate costs a
full forward pass over the ~12 kb prompt. The builder compensates by sampling
~200 negatives per prompt, bringing the job to ~316M tokens.

Every sequence in a batch is the same length -- one fixed prompt plus a
fixed-width candidate -- so no padding and no attention-mask bookkeeping.

Output is appended per prompt so a job killed by the scheduler can resume.

Usage:
  python score_evo2_icl.py --indir data/evo2_icl --out scores.tsv \
      --model Taykhoom/Evo2-7B-1M --batch-size 4
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

csv.field_size_limit(10_000_000)

FIELDS = ["prompt_id", "pair_id", "condition", "query_tf", "dbd_tf", "split",
          "protein_half", "n_shot", "seed", "site_id", "cand_tf", "label",
          "logprob_sum", "logprob_mean"]


def load_jsonl(path: Path) -> list[dict]:
    with path.open() as fh:
        return [json.loads(line) for line in fh if line.strip()]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/evo2_icl"))
    ap.add_argument("--out", type=Path, default=Path("data/evo2_icl/scores.tsv"))
    ap.add_argument("--model", default="Taykhoom/Evo2-7B-1M")
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--limit-prompts", type=int, default=None,
                    help="process only the first N prompts (smoke test)")
    ap.add_argument("--limit-negatives", type=int, default=None,
                    help="keep every positive but cap negatives per prompt "
                         "(pilot runs). Taking the first N candidates instead "
                         "would yield positives only and no AUROC.")
    ap.add_argument("--limit-candidates", type=int, default=None,
                    help="hard cap on candidates per prompt (smoke test only; "
                         "prefer --limit-negatives, which keeps both classes)")
    ap.add_argument("--no-flank", action="store_true",
                    help="score the bare core with no genomic flank (the "
                         "original design; measured at uniform entropy, kept "
                         "only for comparison)")
    ap.add_argument("--shard", type=int, default=0,
                    help="this worker's index, for splitting prompts across GPUs")
    ap.add_argument("--num-shards", type=int, default=1)
    ap.add_argument("--force-fp8", dest="fp8", action="store_true", default=None,
                    help="keep FP8 input projections (needs Hopper + Transformer Engine)")
    ap.add_argument("--no-fp8", dest="fp8", action="store_false",
                    help="disable FP8 input projections (required off Hopper; "
                         "the model card warns outputs then differ noticeably)")
    args = ap.parse_args()

    import torch
    from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer

    if not torch.cuda.is_available():
        print("error: no CUDA device visible", file=sys.stderr)
        return 1
    cap = torch.cuda.get_device_capability()
    name = torch.cuda.get_device_name()
    is_hopper = cap[0] >= 9
    # FP8 needs Transformer Engine, which is a loud stub in this environment,
    # so it can never be the default -- cap[0] >= 9 also catches sm_120
    # Blackwell, which killed a pilot job on an RTX PRO 6000 at model load.
    # The sanity check showed the model is fine without it.
    use_fp8 = False if args.fp8 is None else args.fp8
    print(f"GPU: {name} (sm_{cap[0]}{cap[1]}) | hopper={is_hopper} | fp8_projections={use_fp8}",
          flush=True)
    if use_fp8 and not is_hopper:
        print("warning: FP8 requested on non-Hopper hardware; this will likely fail",
              file=sys.stderr)
    if not use_fp8:
        print("note: FP8 input projections are OFF. The model card states outputs "
              "then differ noticeably from the reference implementation. Fine for "
              "a smoke test, not for the headline run.", flush=True)

    prompts = load_jsonl(args.indir / "prompts.jsonl")
    cands = {}
    with (args.indir / "candidates.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cands[row["site_id"]] = row

    if args.num_shards > 1:
        prompts = [p for i, p in enumerate(prompts) if i % args.num_shards == args.shard]
        print(f"shard {args.shard}/{args.num_shards}: {len(prompts)} prompts", flush=True)
    if args.limit_prompts:
        prompts = prompts[: args.limit_prompts]

    done = set()
    if args.out.is_file():
        with args.out.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                done.add(row["prompt_id"])
    todo = [p for p in prompts if p["prompt_id"] not in done]
    print(f"prompts {len(prompts)} | already done {len(done)} | to do {len(todo)}",
          flush=True)
    if not todo:
        print("nothing to do")
        return 0

    print(f"loading {args.model} ...", flush=True)
    cfg = AutoConfig.from_pretrained(args.model, trust_remote_code=True)
    if hasattr(cfg, "use_fp8_input_projections"):
        cfg.use_fp8_input_projections = bool(use_fp8)
    tok = AutoTokenizer.from_pretrained(args.model, trust_remote_code=True)
    # Pin sdpa explicitly. The port sets _supports_flash_attn_2 = True, so
    # transformers may otherwise auto-select flash_attention_2 and hit the
    # lazy flash_attn import. sdpa is not bit-exact with FA2, but the swap
    # control is a paired contrast -- both arms share this implementation, so
    # the drift cancels in the delta.
    model = AutoModelForCausalLM.from_pretrained(
        args.model, config=cfg, trust_remote_code=True,
        torch_dtype=torch.bfloat16, device_map="cuda",
        attn_implementation="sdpa",
    )
    model.eval()

    def encode(s: str) -> list[int]:
        ids = tok(s, add_special_tokens=False)["input_ids"]
        return ids

    # Evo 2's tokenizer is byte-level: one token per nucleotide. If that does
    # not hold, the candidate slice below would be misaligned and every score
    # silently wrong, so check it rather than assume it.
    probe = "ACGTACGTAC"
    if len(encode(probe)) != len(probe):
        print(f"error: tokenizer is not 1 token per nucleotide "
              f"({len(encode(probe))} tokens for {len(probe)} nt); "
              "candidate alignment would be wrong", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    new_file = not args.out.is_file()
    fh = args.out.open("a", newline="")
    w = csv.writer(fh, delimiter="\t")
    if new_file:
        w.writerow(FIELDS)

    t0 = time.time()
    for i, p in enumerate(todo, 1):
        pid = encode(p["prompt"])
        negs = p["negatives"]
        if args.limit_negatives:
            negs = negs[: args.limit_negatives]
        items = [(s, 1) for s in p["positives"]] + [(s, 0) for s in negs]
        if args.limit_candidates:
            items = items[: args.limit_candidates]

        widths = {len(cands[s]["seq"]) for s, _ in items}
        if len(widths) != 1:
            print(f"error: mixed candidate widths {sorted(widths)}", file=sys.stderr)
            return 1
        width = widths.pop()

        # Each candidate is prompt + its true genomic left flank + the scored
        # core. The flank is what makes the core scoreable at all: without it
        # the model sits at uniform entropy (-1.44/nt vs -0.88 with 800 bp).
        built = []
        for sid, label in items:
            flank = "" if args.no_flank else cands[sid].get("flank", "")
            fids = encode(flank) if flank else []
            cids = encode(cands[sid]["seq"])
            built.append((sid, label, pid + fids + cids, len(pid) + len(fids)))

        # Candidates near the genome start have a truncated flank, so group by
        # (total length, prefix length) to keep every batch rectangular rather
        # than padding and tracking a mask.
        groups: dict[tuple[int, int], list] = {}
        for item in built:
            groups.setdefault((len(item[2]), item[3]), []).append(item)

        rows = []
        for (_tot, plen), members in groups.items():
            for b0 in range(0, len(members), args.batch_size):
                chunk = members[b0 : b0 + args.batch_size]
                ids = torch.tensor([m[2] for m in chunk], dtype=torch.long,
                                   device="cuda")

                with torch.inference_mode():
                    out = model(ids)
                    logits = out.logits if hasattr(out, "logits") else out[0]
                    # Position t's logits predict token t+1, so the logits that
                    # predict the core start one before it.
                    sl = logits[:, plen - 1 : plen + width - 1, :].float()
                    lp = torch.log_softmax(sl, dim=-1)
                    tgt = ids[:, plen : plen + width]
                    tok_lp = lp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1)
                    totals = tok_lp.sum(dim=-1).tolist()

                for (sid, label, _ids, _pl), total in zip(chunk, totals):
                    rows.append([
                        p["prompt_id"], p["pair_id"], p["condition"],
                        p["query_tf"], p["dbd_tf"], p["split"],
                        p["protein_half"], p["n_shot"], p["seed"],
                        sid, cands[sid]["tf_name"], label,
                        f"{total:.5f}", f"{total / width:.6f}",
                    ])

        w.writerows(rows)
        fh.flush()
        os.fsync(fh.fileno())
        el = time.time() - t0
        print(f"[{i}/{len(todo)}] {p['prompt_id']}  {len(rows)} cands  "
              f"{el/60:.1f} min elapsed, ~{el/i*(len(todo)-i)/60:.1f} min left",
              flush=True)

    fh.close()
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
