#!/usr/bin/env python3
"""Does the ICL prompt format help the model, or wreck it?

The sanity check established that Evo 2 scores contiguous E. coli at about
-0.53 log-prob/nt against a uniform floor of -1.386. The ICL smoke test scored
40 bp windows at about -1.40 -- i.e. the model was no better than guessing once
it sat behind our spacer-separated concatenation of coding sequences and site
windows.

That concatenation is far out of distribution for a model trained on contiguous
genomes. If a spacer leaves the model maximally uncertain, there is nothing for
the swap control to measure, and the full run would produce a null result that
says more about the prompt format than about Evo 2.

So score the SAME candidate windows under several contexts:

  cold          the 40 bp alone, no context at all (worst case)
  genomic       the window preceded by its true upstream genomic sequence
                (best case: what the model was trained on)
  icl_full      the real ICL prompt: DBD CDS + sites per block, spacers
  icl_sites     sites only, same spacers, no coding sequences
                (isolates whether the CDS blocks are the problem)
  icl_nospacer  blocks concatenated with no spacer at all

If icl_* sits near uniform while genomic sits near -0.5, the format is the
problem and the experiment needs redesigning before any large run.

Usage:
  python evo2_prompt_format_check.py --indir data/evo2_icl --genome genome.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import statistics as st
from pathlib import Path

csv.field_size_limit(10_000_000)
SEP_SITE = "GGGGGGGG"
SEP_PROT_SITE = "ATGCATGA"
SEP_BLOCK = "ATGCATGAATGCATGA"


def read_genome(path: Path) -> str:
    if path.suffix == ".txt":
        return path.read_text().strip().upper()
    return "".join(l.strip() for l in path.open() if not l.startswith(">")).upper()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/evo2_icl"))
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--model", default="Taykhoom/Evo2-7B-1M")
    ap.add_argument("--n-candidates", type=int, default=10)
    ap.add_argument("--genomic-context", type=int, default=2000)
    ap.add_argument("--query-tf", default=None, help="default: first prompt's TF")
    args = ap.parse_args()

    import torch
    from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer

    genome = read_genome(args.genome)
    cands = {r["site_id"]: r for r in
             csv.DictReader((args.indir / "candidates.tsv").open(), delimiter="\t")}
    prompts = [json.loads(l) for l in (args.indir / "prompts.jsonl").open() if l.strip()]
    rec = next(p for p in prompts if p["condition"] == "real"
               and (args.query_tf is None or p["query_tf"] == args.query_tf))
    print(f"query TF: {rec['query_tf']} | prompt {rec['prompt_nt']} nt | "
          f"{rec['n_shot']} shots", flush=True)

    # Rebuild the sites-only and no-spacer variants from the real prompt by
    # splitting it on the block separator.
    blocks = rec["prompt"].split(SEP_BLOCK)
    sites_only_blocks, nospacer_blocks = [], []
    for b in blocks[:-1]:
        if SEP_PROT_SITE in b:
            _, _, sites = b.partition(SEP_PROT_SITE)
            sites_only_blocks.append(sites)
            nospacer_blocks.append(sites.replace(SEP_SITE, ""))
    icl_sites = SEP_BLOCK.join(sites_only_blocks) + SEP_BLOCK + SEP_PROT_SITE
    icl_nospacer = "".join(nospacer_blocks)

    contexts = {
        "cold": "",
        "icl_full": rec["prompt"],
        "icl_sites": icl_sites,
        "icl_nospacer": icl_nospacer,
    }
    for k, v in contexts.items():
        print(f"  context {k:<14} {len(v):>7} nt")

    cfg = AutoConfig.from_pretrained(args.model, trust_remote_code=True)
    cap = torch.cuda.get_device_capability()
    if hasattr(cfg, "use_fp8_input_projections"):
        cfg.use_fp8_input_projections = cap[0] >= 9
    tok = AutoTokenizer.from_pretrained(args.model, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        args.model, config=cfg, trust_remote_code=True,
        torch_dtype=torch.bfloat16, device_map="cuda",
        attn_implementation="sdpa",
    ).eval()

    def score_tail(context: str, tail: str) -> float:
        """Mean log-prob per nucleotide of `tail`, conditioned on `context`."""
        ctx_ids = tok(context, add_special_tokens=False)["input_ids"] if context else []
        tail_ids = tok(tail, add_special_tokens=False)["input_ids"]
        ids = torch.tensor([ctx_ids + tail_ids], dtype=torch.long, device="cuda")
        n = len(ctx_ids)
        with torch.inference_mode():
            out = model(ids)
            lg = (out.logits if hasattr(out, "logits") else out[0]).float()
            # Position t predicts t+1. With no context, the first tail token has
            # no predictor, so start one in.
            lo = max(n - 1, 0)
            hi = n + len(tail_ids) - 1
            lp = torch.log_softmax(lg[:, lo:hi, :], dim=-1)
            tgt = ids[:, lo + 1 : hi + 1]
            return lp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1).mean().item()

    rng = random.Random(0)
    picks = rec["positives"][: args.n_candidates]
    results: dict[str, list[float]] = {k: [] for k in contexts}
    results["genomic"] = []

    for i, sid in enumerate(picks, 1):
        c = cands[sid]
        seq = c["seq"]
        for k, ctx in contexts.items():
            results[k].append(score_tail(ctx, seq))
        # True upstream genomic context for this exact window.
        mid = int(c["mid"])
        start = mid - len(seq) // 2
        up0 = max(0, start - args.genomic_context)
        results["genomic"].append(score_tail(genome[up0:start], seq))
        print(f"  candidate {i}/{len(picks)}", flush=True)

    print(f"\n{'context':<16}{'log-prob/nt':>13}{'vs cold':>10}")
    cold = st.mean(results["cold"])
    for k in ("genomic", "icl_full", "icl_sites", "icl_nospacer", "cold"):
        v = st.mean(results[k])
        print(f"{k:<16}{v:>13.4f}{v - cold:>+10.4f}")
    print(f"{'ln(1/4)':<16}{-1.3863:>13.4f}")

    g, f = st.mean(results["genomic"]), st.mean(results["icl_full"])
    print()
    print(f"ICL prompt vs true genomic context: {f - g:+.4f} log-prob/nt")
    if f < -1.30:
        print("=> The ICL prompt leaves the model at ~uniform. The format, not the")
        print("   model, is the blocker. Redesign before scaling up.")
    elif f < g - 0.2:
        print("=> The ICL prompt is much worse than natural context, but not dead.")
        print("   Compare icl_sites / icl_nospacer to see which part costs most.")
    else:
        print("=> The ICL prompt holds up. Proceed to the full run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
