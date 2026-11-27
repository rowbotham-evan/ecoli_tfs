#!/usr/bin/env python3
"""How much natural flank does a candidate need to recover predictive confidence?

Established so far, scoring 40 bp windows:

  true genomic context   -0.870 log-prob/nt
  cold, no context       -1.321
  after the ICL prompt   -1.444   (worse than cold)

The candidates in the ICL prompt had NO natural context -- a bare 40 bp window
appended after a spacer. A genome model needs contiguous sequence to re-establish
state after a discontinuity, so that is the worst possible way to present one.

This measures the recovery curve: give the candidate F nucleotides of its TRUE
genomic left flank inside the prompt, score only the central 40 bp core, and
sweep F. Two prefixes are compared at every F:

  icl       the full ICL prompt, then the flank, then the core
  genomic   the flank alone, then the core (no ICL prefix)

That separates two different problems:
  * if icl improves with F and converges on genomic, the discontinuity was the
    only issue and the fix is simply to give candidates natural flanks
  * if icl stays flat while genomic improves, the ICL prefix is itself toxic
    and the demonstrations need redesigning, not just the candidates

Usage:
  python evo2_flank_recovery_check.py --indir data/evo2_icl --genome genome.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics as st
from pathlib import Path

csv.field_size_limit(10_000_000)


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
    ap.add_argument("--n-candidates", type=int, default=8)
    ap.add_argument("--core", type=int, default=40, help="scored width")
    ap.add_argument("--flanks", type=int, nargs="+",
                    default=[0, 50, 100, 200, 400, 800, 1600])
    ap.add_argument("--query-tf", default=None)
    args = ap.parse_args()

    import torch
    from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer

    genome = read_genome(args.genome)
    cands = {r["site_id"]: r for r in
             csv.DictReader((args.indir / "candidates.tsv").open(), delimiter="\t")}
    prompts = [json.loads(l) for l in (args.indir / "prompts.jsonl").open() if l.strip()]
    rec = next(p for p in prompts if p["condition"] == "real"
               and (args.query_tf is None or p["query_tf"] == args.query_tf))
    print(f"query TF {rec['query_tf']} | ICL prefix {rec['prompt_nt']} nt | "
          f"core {args.core} nt | flanks {args.flanks}", flush=True)

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

    def score_core(prefix: str, core: str) -> float:
        pre = tok(prefix, add_special_tokens=False)["input_ids"] if prefix else []
        cor = tok(core, add_special_tokens=False)["input_ids"]
        ids = torch.tensor([pre + cor], dtype=torch.long, device="cuda")
        n = len(pre)
        with torch.inference_mode():
            out = model(ids)
            lg = (out.logits if hasattr(out, "logits") else out[0]).float()
            lo = max(n - 1, 0)
            hi = n + len(cor) - 1
            lp = torch.log_softmax(lg[:, lo:hi, :], dim=-1)
            tgt = ids[:, lo + 1 : hi + 1]
            return lp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1).mean().item()

    picks = rec["positives"][: args.n_candidates]
    res: dict[tuple[str, int], list[float]] = {}

    for i, sid in enumerate(picks, 1):
        c = cands[sid]
        mid = int(c["mid"])
        cs = mid - args.core // 2
        core = genome[cs : cs + args.core]
        if set(core) - set("ACGT"):
            continue
        for F in args.flanks:
            flank = genome[max(0, cs - F) : cs]
            res.setdefault(("genomic", F), []).append(score_core(flank, core))
            res.setdefault(("icl", F), []).append(
                score_core(rec["prompt"] + flank, core))
        print(f"  candidate {i}/{len(picks)}", flush=True)

    print(f"\n{'flank nt':>9}{'icl':>10}{'genomic':>10}{'gap':>9}")
    for F in args.flanks:
        a = st.mean(res[("icl", F)])
        b = st.mean(res[("genomic", F)])
        print(f"{F:>9}{a:>10.4f}{b:>10.4f}{a - b:>+9.4f}")
    print(f"\nreference: uniform -1.3863 | cold -1.3209 | full genomic context -0.8703")

    icl0, iclmax = st.mean(res[("icl", args.flanks[0])]), st.mean(res[("icl", args.flanks[-1])])
    gmax = st.mean(res[("genomic", args.flanks[-1])])
    print()
    print(f"ICL recovery with flank: {icl0:.4f} -> {iclmax:.4f}  ({iclmax - icl0:+.4f})")
    print(f"gap to genomic-only at max flank: {iclmax - gmax:+.4f}")
    if iclmax > -1.0 and abs(iclmax - gmax) < 0.15:
        print("=> Flanks fix it. Candidates just needed natural context; the ICL")
        print("   prefix is not toxic once they have it. Rebuild with flanked candidates.")
    elif iclmax - icl0 > 0.2:
        print("=> Flanks help substantially but a gap to genomic remains. Usable,")
        print("   but the demonstrations likely need lengthening too.")
    else:
        print("=> Flanks do not rescue it. The ICL prefix itself is the problem;")
        print("   lengthening demonstrations is the only remaining option.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
