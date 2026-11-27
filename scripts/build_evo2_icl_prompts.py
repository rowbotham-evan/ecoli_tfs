#!/usr/bin/env python3
"""Build Evo 2 in-context-learning prompts for TF binding-site recognition.

The question: does Evo 2, with no weight updates, already encode a mapping from
a transcription factor's coding sequence to the DNA it binds?

Each prompt is a list of worked examples, one per demonstration TF:

    <TF1 DBD nucleotides> ATGCATGA <site window> GGGGGGGG <site window> ...
    ATGCATGAATGCATGA
    <TF2 DBD nucleotides> ATGCATGA <site window> GGGGGGGG ...
    ATGCATGAATGCATGA
    <QUERY TF DBD nucleotides> ATGCATGA          <- prompt ends here

A scorer then appends each candidate window and sums the log-probability of the
candidate's tokens. Positives are the query TF's held-out sites; negatives are
*other TFs' real sites*, so every candidate is genuine E. coli regulatory DNA
and the model cannot win by recognising regulatory sequence in general.

Swap control (--swap-control, the experiment to run first): for each query TF,
emit a second prompt that is byte-identical except that the final DBD belongs to
a TF from a different family. Same demonstrations, same candidates. If the two
score distributions do not differ, the DBD is doing no work and the result is
pure genome memorisation.

Outputs
  <outdir>/candidates.tsv  every candidate window, emitted once
  <outdir>/prompts.jsonl   one record per (query TF, split, protein half, k, seed)

Usage
  python3 -I scripts/build_evo2_icl_prompts.py <genome.fna> --swap-control
  python3 -I scripts/build_evo2_icl_prompts.py <genome.fna> --sweep
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import statistics as st
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SITES_CSV = REPO / "data" / "regulondb_sites.csv"
DBD_CDS_TSV = REPO / "data" / "dbd_cds.tsv"
DBD_UNIPROT_TSV = REPO / "data" / "sources" / "uniprot" / "dbd_uniprot.tsv"
DEFAULT_OUT = REPO / "data" / "evo2_icl"

# Spacers follow Brixi et al. 2026. They are ordinary nucleotide strings, not
# special tokens; the inner/outer split is this script's extension, since the
# paper only ever separates one level of example.
SEP_SITE = "GGGGGGGG"          # between sites inside one TF block
SEP_PROT_SITE = "ATGCATGA"     # between a TF's coding sequence and its sites
SEP_BLOCK = "ATGCATGAATGCATGA" # between TF blocks

csv.field_size_limit(10_000_000)


def read_genome(path: Path) -> str:
    if path.suffix == ".txt":
        return path.read_text().strip().upper()
    seq = []
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def load_tfs() -> dict:
    """TF -> {dbd_nt, cds_nt, family}, verified entries only."""
    fam = {}
    with DBD_UNIPROT_TSV.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            note = (r.get("dbd_note") or "").strip()
            # "no annotated DBD; compact chain" is not a family; make each such
            # TF its own singleton so leave-one-family-out stays honest.
            if not note or note.startswith("no annotated DBD"):
                note = f"__singleton__{r['tf_name']}"
            fam[r["tf_name"]] = note

    tfs = {}
    with DBD_CDS_TSV.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] != "ok":
                continue
            tfs[r["tf_name"]] = {
                "dbd_nt": r["dbd_cds_nt"],
                "cds_nt": r["cds_nt"],
                "family": fam.get(r["tf_name"], f"__singleton__{r['tf_name']}"),
            }
    return tfs


def load_sites(genome: str, window: int, flank: int, region: int) -> list[dict]:
    """Experimental RegulonDB sites as fixed-width windows centred on the site.

    Each site carries three views of the same locus:
      seq     the `window`-bp core that gets scored
      flank   `flank` bp of true genomic sequence immediately left of the core
      region  a `region`-bp contiguous window centred on the site, used as a
              demonstration so that examples are in-distribution too
    """
    out = []
    with SITES_CSV.open() as fh:
        for r in csv.DictReader(fh):
            if r["experimental"] != "yes":
                continue
            try:
                left, right = int(r["left"]), int(r["right"])
            except ValueError:
                continue
            g0, g1 = left - 1, right                      # 0-based half-open
            if genome[g0:g1] != r["sequence"].upper():    # verified upstream
                continue
            mid = (g0 + g1) // 2
            start = mid - window // 2
            if start < 0 or start + window > len(genome):
                continue
            # Natural context on both sides. A bare 40 bp window dropped after a
            # spacer scores at uniform entropy (-1.44/nt); with 800 bp of true
            # left flank it recovers to -0.88, matching contiguous genomic DNA.
            # So candidates carry a left flank, and demonstration sites are
            # shown as contiguous regions rather than excerpts.
            fl_start = max(0, start - flank)
            rg_start = max(0, mid - region // 2)
            out.append(
                {
                    "site_id": r["site_id"],
                    "tf_name": r["tf_name"],
                    "left": left,
                    "right": right,
                    "mid": mid,
                    "core_len": len(r["sequence"]),
                    "seq": genome[start : start + window],
                    "flank": genome[fl_start:start],
                    "region": genome[rg_start : rg_start + region],
                }
            )
    return out


def build_prompt(demo_tfs, demo_sites, query_tf, tfs, protein_half, demo_key="region"):
    """Assemble one prompt. Ends at the spacer, ready for flank + candidate.

    demo_key selects what represents each demonstration site: "region" for the
    contiguous genomic window (default) or "seq" for the bare core, which is
    what the first design used and which scores at uniform entropy.
    """
    key = "dbd_nt" if protein_half == "dbd" else "cds_nt"
    blocks = []
    for tf in demo_tfs:
        sites = SEP_SITE.join(s[demo_key] for s in demo_sites[tf])
        blocks.append(tfs[tf][key] + SEP_PROT_SITE + sites)
    blocks.append(tfs[query_tf][key] + SEP_PROT_SITE)
    return SEP_BLOCK.join(blocks)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("genome", type=Path, help="GCF_000005845.2_ASM584v2_genomic.fna")
    ap.add_argument("-o", "--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--window", type=int, default=40, help="scored candidate width in bp")
    ap.add_argument("--candidate-flank", type=int, default=800,
                    help="true genomic sequence placed left of each candidate. "
                         "Measured recovery: 0 bp scores -1.44/nt (uniform), "
                         "800 bp -0.88, matching contiguous genomic DNA.")
    ap.add_argument("--demo-region", type=int, default=400,
                    help="contiguous genomic window representing each "
                         "demonstration site, so examples are in-distribution too")
    ap.add_argument("--demo-as-core", action="store_true",
                    help="show demonstrations as bare cores (the original design, "
                         "kept for comparison; scores at uniform entropy)")
    ap.add_argument("--min-sites", type=int, default=10,
                    help="minimum experimental sites for a TF to be eligible")
    ap.add_argument("--sites-per-demo", type=int, default=4,
                    help="sites shown per demonstration TF")
    ap.add_argument("--overlap-margin", type=int, default=None,
                    help="a negative whose midpoint is within this many bp of any "
                         "query-TF site is dropped as ambiguous (default: window)")
    ap.add_argument("--n-shot", type=int, nargs="+", default=[8],
                    help="number of demonstration TFs")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0])
    ap.add_argument("--split", nargs="+", default=["lofo"],
                    choices=["loto", "lofo"])
    ap.add_argument("--protein-half", nargs="+", default=["dbd"],
                    choices=["dbd", "cds"])
    ap.add_argument("--max-negatives", type=int, default=200,
                    help="negatives sampled per prompt. Without vLLM there is no "
                         "prefix caching, so the full 2.4k pool would mean "
                         "recomputing a 12 kb prompt per candidate (~3e9 tokens). "
                         "0 means use every eligible negative.")
    ap.add_argument("--swap-control", action="store_true",
                    help="also emit a paired prompt with the query DBD replaced "
                         "by an out-of-family TF's")
    ap.add_argument("--sweep", action="store_true",
                    help="shorthand for the full grid")
    args = ap.parse_args()

    if args.sweep:
        args.n_shot = [2, 4, 8, 16]
        args.seeds = [0, 1, 2]
        args.split = ["loto", "lofo"]
        args.protein_half = ["dbd", "cds"]
    margin = args.window if args.overlap_margin is None else args.overlap_margin
    demo_key = "seq" if args.demo_as_core else "region"

    if not args.genome.is_file():
        print(f"error: no such file: {args.genome}", file=sys.stderr)
        return 1

    genome = read_genome(args.genome)
    tfs = load_tfs()
    sites = load_sites(genome, args.window, args.candidate_flank, args.demo_region)
    print(f"genome {len(genome):,} bp | verified TFs {len(tfs)} | windowed sites {len(sites)}")

    by_tf = {}
    for s in sites:
        by_tf.setdefault(s["tf_name"], []).append(s)
    eligible = sorted(t for t, ss in by_tf.items() if len(ss) >= args.min_sites and t in tfs)
    print(f"eligible query TFs (>={args.min_sites} sites, verified CDS): {len(eligible)}")

    fams = {}
    for t in eligible:
        fams.setdefault(tfs[t]["family"], []).append(t)
    print(f"families among eligible TFs: {len(fams)}")

    # Candidate pool: every site of every eligible TF, emitted once.
    pool = [s for t in eligible for s in by_tf[t]]
    args.outdir.mkdir(parents=True, exist_ok=True)
    cand_path = args.outdir / "candidates.tsv"
    with cand_path.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["site_id", "tf_name", "left", "right", "mid", "core_len", "seq", "flank"],
            delimiter="\t",
            # 'region' is only used to build demonstration blocks; it never
            # appears in the candidate table.
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(pool)
    print(f"candidate pool: {len(pool)} windows of {args.window} bp -> {cand_path}")

    records, prompt_lens = [], []
    for split in args.split:
        for half in args.protein_half:
            for k in args.n_shot:
                for seed in args.seeds:
                    for q in eligible:
                        rng = random.Random(f"{split}|{half}|{k}|{seed}|{q}")

                        if split == "loto":
                            pool_tfs = [t for t in eligible if t != q]
                        else:
                            qf = tfs[q]["family"]
                            pool_tfs = [t for t in eligible if tfs[t]["family"] != qf]
                        if len(pool_tfs) < k:
                            continue
                        demo_tfs = rng.sample(pool_tfs, k)
                        rng.shuffle(demo_tfs)

                        demo_sites, used = {}, set()
                        for t in demo_tfs:
                            pick = rng.sample(
                                by_tf[t], min(args.sites_per_demo, len(by_tf[t]))
                            )
                            demo_sites[t] = pick
                            used.update(s["site_id"] for s in pick)

                        # Sites shown verbatim in the prompt cannot be scored fairly.
                        # Sites overlapping a query-TF site are ambiguous, not negative.
                        qmids = [s["mid"] for s in by_tf[q]]
                        excluded = set(used)
                        for s in pool:
                            if s["tf_name"] == q:
                                continue
                            if any(abs(s["mid"] - m) < margin for m in qmids):
                                excluded.add(s["site_id"])

                        positives = [
                            s["site_id"] for s in by_tf[q] if s["site_id"] not in excluded
                        ]
                        if not positives:
                            continue

                        negatives = [
                            s["site_id"] for s in pool
                            if s["tf_name"] != q and s["site_id"] not in excluded
                        ]
                        if args.max_negatives and len(negatives) > args.max_negatives:
                            # Seeded on the pair, not the condition, so the real
                            # and swap prompts always share a candidate set.
                            negatives = random.Random(
                                f"neg|{split}|{half}|{k}|{seed}|{q}"
                            ).sample(negatives, args.max_negatives)

                        prompt = build_prompt(demo_tfs, demo_sites, q, tfs, half, demo_key)
                        prompt_lens.append(len(prompt))
                        base = {
                            "split": split,
                            "protein_half": half,
                            "n_shot": k,
                            "seed": seed,
                            "query_tf": q,
                            "demo_tfs": demo_tfs,
                            "positives": positives,
                            "negatives": negatives,
                            "excluded": sorted(excluded),
                            "prompt_nt": len(prompt),
                        }
                        records.append(
                            {
                                **base,
                                "prompt_id": f"{split}|{half}|k{k}|s{seed}|{q}|real",
                                "pair_id": f"{split}|{half}|k{k}|s{seed}|{q}",
                                "condition": "real",
                                "dbd_tf": q,
                                "prompt": prompt,
                            }
                        )

                        if args.swap_control:
                            qf = tfs[q]["family"]
                            others = [t for t in eligible if tfs[t]["family"] != qf]
                            if not others:
                                continue
                            swap_tf = random.Random(f"swap|{q}|{seed}").choice(others)
                            swapped = build_prompt(
                                demo_tfs, demo_sites, swap_tf, tfs, half, demo_key
                            )
                            records.append(
                                {
                                    **base,
                                    "prompt_id": f"{split}|{half}|k{k}|s{seed}|{q}|swap",
                                    "pair_id": f"{split}|{half}|k{k}|s{seed}|{q}",
                                    "condition": "swap",
                                    "dbd_tf": swap_tf,
                                    "prompt": swapped,
                                    "prompt_nt": len(swapped),
                                }
                            )

    prompts_path = args.outdir / "prompts.jsonl"
    with prompts_path.open("w") as fh:
        for r in records:
            fh.write(json.dumps(r) + "\n")

    print(f"\nprompts: {len(records)} -> {prompts_path}")
    if prompt_lens:
        print(
            "prompt length nt: min %d  median %d  max %d"
            % (min(prompt_lens), st.median(prompt_lens), max(prompt_lens))
        )
        n_cand = [len(r["positives"]) + len(r["negatives"]) for r in records]
        print("candidates per prompt: median %d  total %d"
              % (st.median(n_cand), sum(n_cand)))
        # No prefix caching, so each candidate costs a full forward over
        # prompt + its flank + the scored core.
        per = args.prompt_nt_overhead = args.candidate_flank + args.window
        fwd = sum((r["prompt_nt"] + per) *
                  (len(r["positives"]) + len(r["negatives"])) for r in records)
        print("tokens to forward (no prefix cache): %.0fM" % (fwd / 1e6))
        print("  = %d candidates x (prompt + %d flank + %d core)"
              % (sum(n_cand), args.candidate_flank, args.window))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
