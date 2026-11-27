#!/usr/bin/env python3
"""Is this Evo 2 checkpoint actually working? Gate the full run on the answer.

The smoke test scored real E. coli DNA at about -1.40 log-prob per nucleotide.
Uniform over {A,C,G,T} is ln(1/4) = -1.3863. A genomic language model that has
seen prokaryotic genomes should find E. coli far more predictable than a coin
flip, so scoring at uniform means something is broken -- most likely the FP8
input projections running in bf16 without the checkpoint's calibration scales.

This distinguishes "model works" from "model is noise" in about a minute, with
no dependence on the ICL setup at all:

  real        contiguous E. coli genomic sequence
  dinuc       the same sequence, dinucleotide-shuffled (same 2-mer stats,
              destroyed long-range structure)
  mono        the same sequence, mononucleotide-shuffled (same base composition)
  uniform     i.i.d. uniform ACGT

A working model: real < dinuc < mono < uniform in entropy, with real clearly
below uniform. A broken model: all four land on top of each other near -1.386.

Run it both ways (--no-fp8 and, on Hopper with Transformer Engine, --force-fp8)
to find out whether FP8 is the culprit.

Usage:
  python evo2_sanity_check.py --genome genome.fna --model Taykhoom/Evo2-7B-1M --no-fp8
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path


def read_genome(path: Path) -> str:
    if path.suffix == ".txt":
        return path.read_text().strip().upper()
    return "".join(l.strip() for l in path.open() if not l.startswith(">")).upper()


def mono_shuffle(s: str, rng: random.Random) -> str:
    x = list(s)
    rng.shuffle(x)
    return "".join(x)


def dinuc_shuffle(s: str, rng: random.Random) -> str:
    """Altschul-Erikson: Eulerian walk preserving exact dinucleotide counts."""
    if len(s) < 3:
        return s
    edges: dict[str, list[str]] = {}
    for a, b in zip(s, s[1:]):
        edges.setdefault(a, []).append(b)
    last = s[-1]
    # Pick a random last-edge per vertex forming a tree rooted at `last`.
    tree: dict[str, str] = {}
    for v in list(edges):
        if v == last:
            continue
        while True:
            cand = rng.choice(edges[v])
            seen, u, ok = {v}, cand, False
            while True:
                if u == last:
                    ok = True
                    break
                if u in seen or u not in tree:
                    break
                seen.add(u)
                u = tree[u]
            if ok or cand == last:
                tree[v] = cand
                break
            if all(e == cand for e in edges[v]):
                tree[v] = cand
                break
    for v, t in tree.items():
        edges[v].remove(t)
    for v in edges:
        rng.shuffle(edges[v])
    for v, t in tree.items():
        edges[v].append(t)
    out, cur = [s[0]], s[0]
    for _ in range(len(s) - 1):
        if not edges.get(cur):
            break
        nxt = edges[cur].pop(0)
        out.append(nxt)
        cur = nxt
    res = "".join(out)
    return res if len(res) == len(s) else s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--model", default="Taykhoom/Evo2-7B-1M")
    ap.add_argument("--length", type=int, default=2000)
    ap.add_argument("--n-windows", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--force-fp8", dest="fp8", action="store_true", default=None)
    ap.add_argument("--no-fp8", dest="fp8", action="store_false")
    args = ap.parse_args()

    import torch
    from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer

    cap = torch.cuda.get_device_capability()
    is_hopper = cap[0] >= 9
    use_fp8 = is_hopper if args.fp8 is None else args.fp8
    print(f"GPU {torch.cuda.get_device_name()} sm_{cap[0]}{cap[1]} | fp8={use_fp8}", flush=True)

    g = read_genome(args.genome)
    rng = random.Random(args.seed)
    starts = [rng.randrange(0, len(g) - args.length) for _ in range(args.n_windows)]

    cfg = AutoConfig.from_pretrained(args.model, trust_remote_code=True)
    if hasattr(cfg, "use_fp8_input_projections"):
        cfg.use_fp8_input_projections = bool(use_fp8)
    tok = AutoTokenizer.from_pretrained(args.model, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        args.model, config=cfg, trust_remote_code=True,
        torch_dtype=torch.bfloat16, device_map="cuda",
        attn_implementation="sdpa",
    ).eval()

    def score(seq: str) -> float:
        ids = torch.tensor([tok(seq, add_special_tokens=False)["input_ids"]],
                           dtype=torch.long, device="cuda")
        with torch.inference_mode():
            out = model(ids)
            lg = (out.logits if hasattr(out, "logits") else out[0]).float()
            lp = torch.log_softmax(lg[:, :-1, :], dim=-1)
            tgt = ids[:, 1:]
            return lp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1).mean().item()

    kinds = {"real": [], "dinuc": [], "mono": [], "uniform": []}
    for i, st0 in enumerate(starts, 1):
        real = g[st0 : st0 + args.length]
        if set(real) - set("ACGT"):
            continue
        kinds["real"].append(score(real))
        kinds["dinuc"].append(score(dinuc_shuffle(real, rng)))
        kinds["mono"].append(score(mono_shuffle(real, rng)))
        kinds["uniform"].append(score("".join(rng.choice("ACGT") for _ in range(args.length))))
        print(f"  window {i}/{len(starts)} done", flush=True)

    print(f"\n{'sequence':<10}{'log-prob/nt':>14}{'vs uniform':>13}")
    uni = sum(kinds["uniform"]) / len(kinds["uniform"])
    for k in ("real", "dinuc", "mono", "uniform"):
        v = sum(kinds[k]) / len(kinds[k])
        print(f"{k:<10}{v:>14.4f}{v - uni:>+13.4f}")
    print(f"{'ln(1/4)':<10}{-1.3863:>14.4f}")

    real_m = sum(kinds["real"]) / len(kinds["real"])
    print()
    if real_m - uni > 0.15:
        print("=> WORKING. Real E. coli is clearly more predictable than random.")
        return 0
    print("=> BROKEN. Real genomic DNA scores no better than random sequence;")
    print("   this checkpoint is not usable for scoring in this configuration.")
    print("   Do not spend GPU hours on the full run until this is resolved.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
