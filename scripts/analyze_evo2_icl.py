#!/usr/bin/env python3
"""Analyse Evo 2 ICL scores: per-TF AUROC and the swap-control contrast.

For each prompt, positives are the query TF's held-out sites and negatives are
other TFs' real sites. Sites shown verbatim in the prompt, and sites overlapping
a query-TF site, were marked excluded at build time and are dropped here.

The headline number is not the AUROC itself. Evo 2 has memorised the E. coli
genome, so a high AUROC under the real prompt proves little on its own. What
matters is the paired contrast:

    delta = AUROC(real query DBD) - AUROC(out-of-family DBD)

Both prompts share identical demonstrations and identical candidates, so delta
isolates the contribution of the query TF's own coding sequence. A delta
indistinguishable from zero means the DBD is doing no work.

Usage:
  python3 -I scripts/analyze_evo2_icl.py --scores data/evo2_icl/scores.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)


def auroc(pos: list[float], neg: list[float]) -> float | None:
    """Mann-Whitney U / |pos||neg|, with ties at half credit."""
    if not pos or not neg:
        return None
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks, i, n = {}, 0, len(merged)
    rank_sum_pos = 0.0
    while i < n:
        j = i
        while j + 1 < n and merged[j + 1][0] == merged[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0          # 1-based average rank over the tie
        for k in range(i, j + 1):
            if merged[k][1] == 1:
                rank_sum_pos += avg_rank
        i = j + 1
    npos, nneg = len(pos), len(neg)
    u = rank_sum_pos - npos * (npos + 1) / 2.0
    return u / (npos * nneg)


def wilcoxon_signed_rank_p(deltas: list[float]) -> float | None:
    """Two-sided p for the paired real-minus-swap deltas.

    Exact by enumeration for n <= 20 and no tied |d|; the normal
    approximation only above that. At the pilot's n=6 the approximation is
    badly anticonservative -- it returned 0.046 where the exact test gives
    0.0625, which is in fact the SMALLEST p a 5-of-6 split can produce.
    """
    d = [x for x in deltas if x != 0.0]
    n = len(d)
    if n < 6:
        return None
    absd = [abs(x) for x in d]
    if n <= 20 and len(set(absd)) == n:
        order = sorted(range(n), key=lambda i: absd[i])
        rank = [0] * n
        for pos, i in enumerate(order, 1):
            rank[i] = pos
        w_plus_obs = sum(r for r, x in zip(rank, d) if x > 0)
        total = n * (n + 1) // 2
        obs = min(w_plus_obs, total - w_plus_obs)
        hits = 0
        for mask in range(1 << n):
            w = sum(k + 1 for k in range(n) if mask >> k & 1)
            if min(w, total - w) <= obs:
                hits += 1
        return hits / float(1 << n)
    order = sorted(range(n), key=lambda i: abs(d[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(d[order[j + 1]]) == abs(d[order[i]]):
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    w_plus = sum(r for r, x in zip(ranks, d) if x > 0)
    mean = n * (n + 1) / 4.0
    sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    if sd == 0:
        return None
    z = (w_plus - mean) / sd
    return 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / math.sqrt(2.0))))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", type=Path, default=Path("data/evo2_icl/scores.tsv"))
    ap.add_argument("--prompts", type=Path, default=Path("data/evo2_icl/prompts.jsonl"))
    ap.add_argument("--metric", default="logprob_sum",
                    choices=["logprob_sum", "logprob_mean"])
    ap.add_argument("--out", type=Path, default=None, help="write per-TF table here")
    args = ap.parse_args()

    meta = {}
    with args.prompts.open() as fh:
        for line in fh:
            if line.strip():
                r = json.loads(line)
                meta[r["prompt_id"]] = r

    # prompt_id -> ([positive scores], [negative scores]), split by the scorer's
    # own label column rather than re-deriving membership here.
    scores: dict[str, tuple[list[float], list[float]]] = {}
    with args.scores.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "label" not in (reader.fieldnames or []):
            print("error: scores file has no 'label' column; it predates the "
                  "current scorer. Re-run scoring.", file=sys.stderr)
            return 1
        for row in reader:
            pos, neg = scores.setdefault(row["prompt_id"], ([], []))
            (pos if row["label"] == "1" else neg).append(float(row[args.metric]))
    print(f"prompts scored: {len(scores)}")

    per_prompt = {}
    for pid, (pos, neg) in scores.items():
        m = meta.get(pid)
        if m is None:
            continue
        a = auroc(pos, neg)
        if a is not None:
            per_prompt[pid] = {"auroc": a, "n_pos": len(pos), "n_neg": len(neg), **m}

    # Pair real against swap.
    pairs = {}
    for pid, r in per_prompt.items():
        pairs.setdefault(r["pair_id"], {})[r["condition"]] = r

    rows = []
    for pair_id, p in sorted(pairs.items()):
        if "real" not in p:
            continue
        real = p["real"]
        swap = p.get("swap")
        rows.append({
            "pair_id": pair_id,
            "query_tf": real["query_tf"],
            "split": real["split"],
            "protein_half": real["protein_half"],
            "n_shot": real["n_shot"],
            "seed": real["seed"],
            "n_pos": real["n_pos"],
            "n_neg": real["n_neg"],
            "auroc_real": round(real["auroc"], 4),
            "auroc_swap": round(swap["auroc"], 4) if swap else "",
            "swap_tf": swap["dbd_tf"] if swap else "",
            "delta": round(real["auroc"] - swap["auroc"], 4) if swap else "",
        })

    if not rows:
        print("no scored pairs found")
        return 1

    print(f"\n{'TF':<10} {'n+':>4} {'n-':>5} {'AUROC real':>11} {'AUROC swap':>11} {'delta':>8}")
    for r in sorted(rows, key=lambda x: -(x["delta"] if x["delta"] != "" else -9)):
        print(f"{r['query_tf']:<10} {r['n_pos']:>4} {r['n_neg']:>5} "
              f"{r['auroc_real']:>11} {str(r['auroc_swap']):>11} {str(r['delta']):>8}")

    reals = [r["auroc_real"] for r in rows]
    deltas = [r["delta"] for r in rows if r["delta"] != ""]
    print(f"\nAUROC (real DBD):  median {st.median(reals):.4f}  "
          f"mean {st.mean(reals):.4f}  n={len(reals)}")
    if deltas:
        pos_n = sum(1 for d in deltas if d > 0)
        print(f"AUROC (swap DBD):  median "
              f"{st.median([r['auroc_swap'] for r in rows if r['auroc_swap'] != '']):.4f}")
        print(f"delta real-swap:   median {st.median(deltas):+.4f}  "
              f"mean {st.mean(deltas):+.4f}")
        print(f"                   positive in {pos_n}/{len(deltas)} TFs")
        p = wilcoxon_signed_rank_p(deltas)
        if p is not None:
            print(f"                   Wilcoxon signed-rank p = {p:.3g}")
        print()
        if p is not None and p < 0.05 and st.median(deltas) > 0:
            print("=> The query TF's own coding sequence shifts the ranking. Signal present.")
        else:
            print("=> No detectable effect of the query DBD. The AUROC above, however")
            print("   high, is explained by genome memorisation rather than by the TF.")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
