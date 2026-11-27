#!/usr/bin/env python3
"""Per-TF PWM baseline on exactly the candidate set Evo 2 is scored against.

Without this, an Evo 2 AUROC has no scale. This supplies two reference points,
computed over the identical candidates, positives and negatives that
build_evo2_icl_prompts.py produced:

  ceiling  A PWM fitted on half of the query TF's OWN sites, scoring the other
           half. This is not a fair competitor to in-context learning -- it
           sees the query TF's sites, which the ICL prompts deliberately never
           do -- it is the practical upper bound for this candidate set. If a
           TF's own PWM cannot separate its sites from other TFs' sites, no
           method will, and that TF tells you nothing.

  floor    A PWM fitted on a DIFFERENT TF's sites, scoring the query TF's.
           The PWM analogue of the swap control.

Evo 2 operating between floor and ceiling is the interesting outcome. Evo 2 at
the floor means the DBD is doing nothing.

Usage:
  python3 -I scripts/pwm_baseline_evo2_icl.py --indir data/evo2_icl
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics as st
from pathlib import Path

csv.field_size_limit(10_000_000)
BASES = "ACGT"
IDX = {b: i for i, b in enumerate(BASES)}
COMP = str.maketrans("ACGT", "TGCA")


def auroc(pos: list[float], neg: list[float]) -> float | None:
    if not pos or not neg:
        return None
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    n = len(merged)
    rank_sum_pos, i = 0.0, 0
    while i < n:
        j = i
        while j + 1 < n and merged[j + 1][0] == merged[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            if merged[k][1] == 1:
                rank_sum_pos += avg
        i = j + 1
    npos, nneg = len(pos), len(neg)
    return (rank_sum_pos - npos * (npos + 1) / 2.0) / (npos * nneg)


def fit_pwm(seqs: list[str], width: int, bg: list[float], pseudo: float = 0.25):
    """Log-odds matrix from centre-aligned sequences."""
    counts = [[pseudo] * 4 for _ in range(width)]
    used = 0
    for s in seqs:
        if len(s) < width:
            continue
        off = (len(s) - width) // 2
        sub = s[off : off + width]
        if any(c not in IDX for c in sub):
            continue
        for i, c in enumerate(sub):
            counts[i][IDX[c]] += 1
        used += 1
    if used == 0:
        return None
    pwm = []
    for col in counts:
        tot = sum(col)
        pwm.append([math.log((col[k] / tot) / bg[k]) for k in range(4)])
    return pwm


def scan(pwm, seq: str) -> float:
    """Best log-odds over all offsets and both strands."""
    w = len(pwm)
    best = -1e9
    for s in (seq, seq.translate(COMP)[::-1]):
        for off in range(len(s) - w + 1):
            tot = 0.0
            ok = True
            for i in range(w):
                c = s[off + i]
                if c not in IDX:
                    ok = False
                    break
                tot += pwm[i][IDX[c]]
            if ok and tot > best:
                best = tot
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/evo2_icl"))
    ap.add_argument("--pwm-width", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("data/evo2_icl/pwm_baseline.tsv"))
    args = ap.parse_args()

    cands = list(csv.DictReader((args.indir / "candidates.tsv").open(), delimiter="\t"))
    by_id = {c["site_id"]: c for c in cands}
    by_tf: dict[str, list[dict]] = {}
    for c in cands:
        by_tf.setdefault(c["tf_name"], []).append(c)

    # Background composition from the candidate pool itself.
    tot = [0] * 4
    for c in cands:
        for ch in c["seq"]:
            if ch in IDX:
                tot[IDX[ch]] += 1
    n = sum(tot)
    bg = [t / n for t in tot]
    print(f"candidates {len(cands)} | TFs {len(by_tf)} | "
          f"background A{bg[0]:.3f} C{bg[1]:.3f} G{bg[2]:.3f} T{bg[3]:.3f}")

    # Use the 'real' prompt records purely for their exclusion bookkeeping, so
    # positives and negatives match the Evo 2 run exactly.
    meta = {}
    with (args.indir / "prompts.jsonl").open() as fh:
        for line in fh:
            r = json.loads(line)
            if r["condition"] == "real":
                meta[r["query_tf"]] = r

    rows = []
    tfs = sorted(meta)
    for q in tfs:
        m = meta[q]
        excluded = set(m["excluded"])
        own = [c for c in by_tf.get(q, []) if c["site_id"] not in excluded]
        if len(own) < 6:
            continue

        rng = random.Random(f"{args.seed}|{q}")
        shuffled = own[:]
        rng.shuffle(shuffled)
        half = len(shuffled) // 2
        train, test = shuffled[:half], shuffled[half:]

        negs = [c for c in cands
                if c["tf_name"] != q and c["site_id"] not in excluded]

        pwm = fit_pwm([c["seq"] for c in train], args.pwm_width, bg)
        if pwm is None:
            continue
        pos_scores = [scan(pwm, c["seq"]) for c in test]
        neg_scores = [scan(pwm, c["seq"]) for c in negs]
        ceiling = auroc(pos_scores, neg_scores)

        # Floor: a PWM from an unrelated TF, scoring this TF's sites.
        other = rng.choice([t for t in tfs if t != q])
        other_sites = [c["seq"] for c in by_tf.get(other, [])]
        pwm_o = fit_pwm(other_sites, args.pwm_width, bg)
        floor = None
        if pwm_o is not None:
            floor = auroc([scan(pwm_o, c["seq"]) for c in test],
                          [scan(pwm_o, c["seq"]) for c in negs])

        rows.append({
            "query_tf": q,
            "n_train": len(train),
            "n_test": len(test),
            "n_neg": len(negs),
            "auroc_own_pwm": round(ceiling, 4) if ceiling is not None else "",
            "auroc_other_pwm": round(floor, 4) if floor is not None else "",
            "other_tf": other,
        })
        print(f"  {q:<10} train {len(train):>3} test {len(test):>3}  "
              f"own {ceiling:.3f}  other({other}) "
              f"{floor:.3f}" if floor is not None else "", flush=True)

    if not rows:
        print("no TFs had enough sites")
        return 1

    ceil_vals = [r["auroc_own_pwm"] for r in rows if r["auroc_own_pwm"] != ""]
    floor_vals = [r["auroc_other_pwm"] for r in rows if r["auroc_other_pwm"] != ""]
    print(f"\n{'':<22}{'median':>8}{'mean':>8}{'min':>8}{'max':>8}")
    print(f"{'PWM, own sites':<22}{st.median(ceil_vals):>8.4f}{st.mean(ceil_vals):>8.4f}"
          f"{min(ceil_vals):>8.4f}{max(ceil_vals):>8.4f}")
    print(f"{'PWM, other TF sites':<22}{st.median(floor_vals):>8.4f}"
          f"{st.mean(floor_vals):>8.4f}{min(floor_vals):>8.4f}{max(floor_vals):>8.4f}")
    print(f"\nn TFs: {len(rows)}")
    print("Evo 2 should be read against this range, not against 0.5.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
