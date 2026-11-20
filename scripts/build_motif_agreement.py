#!/usr/bin/env python3
"""Build data/motif_agreement.json: how well each TF's two published motifs agree.

For every TF with a DAP-seq motif and a second motif, the two are turned into
position probability matrices (PPMs), aligned with Tomtom (MEME Suite), allowing
shifts and reverse complements, and scored by the Pearson correlation between
the aligned PPM columns.

The motifs:
  - DAP-seq (Baumgart et al. 2021): a frequency matrix, used as the PPM as is.
  - ChIP-seq (Lally et al. 2025): a BoltzNet energy matrix, where a weight above
    zero favors binding. It is turned into a PPM with the Boltzmann relation,
    p(base) = exp(w) / sum over the four bases of exp(w), at every position.
  - PhoB only (Fitzgerald et al. 2023): the aligned ChIP-seq sites stored in
    data/dataset.csv, counted per position.

Tomtom is run once, every DAP-seq motif against every second motif, with
`-dist pearson -min-overlap MIN_OVERLAP -thresh 1.0`, so it reports its best
alignment (offset and orientation) for all pairs. The Pearson correlation is then
computed here from the untouched PPMs over the aligned columns. A TF's own pair
gives its score; the pairs of different TFs give the "mismatched" reference
distribution, since choosing the best alignment raises every score.

The same all-pairs scores answer a second question: does a TF's identity survive the
change of method? The output's "identity" block holds the full cross-method matrix
(DAP-seq motif of TF i against the second motif of TF j) and, for every TF, the rank
of its own partner among all candidates in both directions (rank 1 = retrieved).
Because related TFs recognize similar sequences, each TF's DNA-binding family (the
UniProt domain note in data/sources/uniprot/dbd_uniprot.tsv) is recorded too, and
the retrieval is repeated among same-family candidates only.

Needs the `tomtom` program: `conda create -n meme -c conda-forge -c bioconda meme`,
then `conda run -n meme python3 scripts/build_motif_agreement.py`, or pass
`--tomtom PATH`. The output is rewritten on every run.
"""

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATRICES = ROOT / "data" / "motif_matrices.json"
DATASET = ROOT / "data" / "dataset.csv"
OUTPUT = ROOT / "data" / "motif_agreement.json"
UNIPROT = ROOT / "data" / "sources" / "uniprot" / "dbd_uniprot.tsv"
# UniProt domain notes that describe a motif or the lack of one, not a family.
NOT_A_FAMILY = {"", "H-T-H motif", "no annotated DBD; compact chain"}
MIN_OVERLAP = 8       # aligned positions Tomtom must keep
MIN_SITES = 10        # aligned sites needed to build a PPM from dataset.csv
BASES = "ACGT"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def normalize(matrix):
    return [[v / (sum(row) or 1) for v in row] for row in matrix]


def boltzmann(matrix):
    """Energy matrix -> PPM: a base's probability is proportional to exp(weight)."""
    return normalize([[math.exp(w) for w in row] for row in matrix])


def reverse_complement(matrix):
    return [row[::-1] for row in matrix[::-1]]


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx, syy = sum((x - mx) ** 2 for x in xs), sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0


def load_pairs():
    """One entry per TF with two motifs: its DAP-seq PPM and its second PPM."""
    matrices = json.loads(MATRICES.read_text())
    sites = {}
    with DATASET.open(newline="") as handle:
        for row in csv.DictReader(handle):
            seq = (row.get("dna_binding_sequence") or "").strip().upper()
            if seq:
                sites.setdefault(row["tf_name"], []).append(seq)
    pairs = []
    for name, papers in matrices.items():
        entries = list(papers.values())
        freq = next((e for e in entries if e["kind"] == "frequency"), None)
        energy = next((e for e in entries if e["kind"] == "energy"), None)
        if not freq:
            continue
        if energy:
            second, other = boltzmann(energy["matrix"]), "ChIP-seq energy matrix (Lally et al.)"
        else:
            seqs = sites.get(name, [])
            if len(seqs) < MIN_SITES or len({len(s) for s in seqs}) != 1 or any(set(s) - set(BASES) for s in seqs):
                continue
            second = normalize([[sum(s[i] == b for s in seqs) for b in BASES] for i in range(len(seqs[0]))])
            other = "ChIP-seq aligned sites (Fitzgerald et al.)"
        pairs.append({"name": name, "other": other, "dap": normalize(freq["matrix"]), "second": second})
    require(pairs, "No TF has two motifs to compare.")
    return pairs


def load_families():
    """TF name -> DNA-binding family, for single-protein TFs whose UniProt domain names one."""
    notes = {}
    with UNIPROT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            notes.setdefault(row["tf_name"], []).append(row["dbd_note"].strip())
    return {name: found[0].replace("HTH ", "").replace("-type", "") for name, found in notes.items()
            if len(found) == 1 and found[0] not in NOT_A_FAMILY}


def identity(pairs, scores):
    """Cross-method retrieval: is a TF's own partner its most similar motif?

    scores[(q, t)] is the similarity of TF q's DAP-seq motif to TF t's second motif."""
    families = load_families()
    family = [families.get(p["name"]) for p in pairs]
    sizes = {f: family.count(f) for f in set(family) if f}
    # Families with several members first (largest first), so they form blocks in the matrix.
    order = sorted(range(len(pairs)), key=lambda i: (
        -sizes.get(family[i], 0) if sizes.get(family[i], 0) > 1 else 1, family[i] or "~", pairs[i]["name"].lower()))
    n = len(pairs)

    def rank(own, candidates, score):
        """1 + the number of candidates scoring above the TF's own partner."""
        return 1 + sum(score(c) > score(own) for c in candidates if c != own)

    tfs = []
    for i in order:
        others = [j for j in range(n) if j != i]
        mates = [j for j in others if family[i] and family[j] == family[i]]
        best_forward = max(range(n), key=lambda j: scores[(i, j)])    # which second motif the DAP-seq motif retrieves
        best_backward = max(range(n), key=lambda j: scores[(j, i)])   # which DAP-seq motif the second motif retrieves
        tfs.append({
            "name": pairs[i]["name"], "family": family[i] if sizes.get(family[i], 0) > 1 else None,
            "rank_forward": rank(i, range(n), lambda j: scores[(i, j)]),
            "rank_backward": rank(i, range(n), lambda j: scores[(j, i)]),
            "best_forward": pairs[best_forward]["name"], "best_backward": pairs[best_backward]["name"],
            "family_mates": len(mates),
            # The same question among same-family TFs only (None when the TF has no family-mate here).
            "family_rank_forward": rank(i, [i] + mates, lambda j: scores[(i, j)]) if mates else None,
            "family_rank_backward": rank(i, [i] + mates, lambda j: scores[(j, i)]) if mates else None,
        })
    mean = lambda values: round(sum(values) / len(values), 4) if values else None
    same = [scores[(i, j)] for i in range(n) for j in range(n) if i != j and family[i] and family[i] == family[j]]
    different = [scores[(i, j)] for i in range(n) for j in range(n) if i != j and not (family[i] and family[i] == family[j])]
    return {
        "tfs": tfs,
        "matrix": [[round(scores[(i, j)], 3) for j in order] for i in order],   # rows: DAP-seq motif; columns: second motif
        "mean_own": mean([scores[(i, i)] for i in range(n)]),
        "mean_same_family": mean(same), "mean_other_family": mean(different),
    }


def write_meme(path, motifs):
    lines = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25", ""]
    for motif_id, ppm in motifs:
        lines += [f"MOTIF {motif_id}", f"letter-probability matrix: alength= 4 w= {len(ppm)} nsites= 20 E= 0"]
        lines += [" ".join(f"{v:.6f}" for v in row) for row in ppm] + [""]
    path.write_text("\n".join(lines))


def run_tomtom(tomtom, pairs):
    """Tomtom's best alignment for every (DAP-seq motif, second motif) pair."""
    args = ["-dist", "pearson", "-min-overlap", str(MIN_OVERLAP), "-thresh", "1.0", "-text", "-verbosity", "1"]
    with tempfile.TemporaryDirectory() as tmp:
        query, target = Path(tmp) / "dap.meme", Path(tmp) / "second.meme"
        # Tomtom ids cannot hold spaces or parentheses, so motifs are numbered.
        write_meme(query, [(f"q{i}", p["dap"]) for i, p in enumerate(pairs)])
        write_meme(target, [(f"t{i}", p["second"]) for i, p in enumerate(pairs)])
        out = subprocess.run([tomtom, *args, str(query), str(target)], check=True, capture_output=True, text=True).stdout
    rows = [r for r in csv.DictReader((line for line in out.splitlines() if line and not line.startswith("#")), delimiter="\t")]
    hits = {}
    for r in rows:
        hits[(int(r["Query_ID"][1:]), int(r["Target_ID"][1:]))] = {
            "offset": int(r["Optimal_offset"]), "overlap": int(r["Overlap"]),
            "reversed": r["Orientation"] == "-", "p": float(r["p-value"]),
        }
    require(len(hits) == len(pairs) ** 2, f"Tomtom reported {len(hits)} of {len(pairs) ** 2} pairs.")
    return hits, args


def overlap_columns(query_len, target_len, shift):
    """Query columns i that sit over target column i + shift."""
    return [i for i in range(query_len) if 0 <= i + shift < target_len]


def offset_sign(pairs, hits):
    """Whether Tomtom's offset is added to or subtracted from the query column.

    Settled from the data: only one sign reproduces Tomtom's own overlap counts."""
    fits = [sign for sign in (1, -1) if all(
        len(overlap_columns(len(pairs[q]["dap"]), len(pairs[t]["second"]), sign * h["offset"])) == h["overlap"]
        for (q, t), h in hits.items())]
    require(len(fits) == 1, f"Could not settle Tomtom's offset convention (candidates: {fits}).")
    return fits[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tomtom", default=shutil.which("tomtom"), help="path to the tomtom program (default: from PATH)")
    tomtom = parser.parse_args().tomtom
    require(tomtom, "tomtom not found: install the MEME Suite or pass --tomtom PATH.")
    version = subprocess.run([tomtom, "-version"], check=True, capture_output=True, text=True).stdout.strip()

    pairs = load_pairs()
    hits, args = run_tomtom(tomtom, pairs)
    sign = offset_sign(pairs, hits)

    tfs, mismatched, scores = [], [], {}
    for (q, t), hit in sorted(hits.items()):
        dap = pairs[q]["dap"]
        second = reverse_complement(pairs[t]["second"]) if hit["reversed"] else pairs[t]["second"]
        shift = sign * hit["offset"]
        columns = overlap_columns(len(dap), len(second), shift)
        r = pearson([v for i in columns for v in dap[i]], [v for i in columns for v in second[i + shift]])
        scores[(q, t)] = r
        if q != t:
            mismatched.append(round(r, 4))
            continue
        tfs.append({
            "name": pairs[q]["name"], "other": pairs[q]["other"], "r": round(r, 4),
            "overlap": hit["overlap"], "reversed": hit["reversed"], "tomtom_p": hit["p"],
            # The aligned cells: [DAP-seq probability, other motif's probability, base index, DAP-seq position].
            "points": [[round(dap[i][k], 4), round(second[i + shift][k], 4), k, i + 1] for i in columns for k in range(4)],
        })
    mismatched.sort()
    for tf in tfs:
        tf["p"] = round(sum(v >= tf["r"] for v in mismatched) / len(mismatched), 4) if mismatched else 1
    tfs.sort(key=lambda tf: -tf["r"])
    result = {
        "method": {
            "aligner": f"Tomtom {version} (MEME Suite)", "tomtom_args": " ".join(args),
            "score": "Pearson correlation between the aligned position probability matrices",
            "energy_to_ppm": "p = exp(w) / sum(exp(w)) at each position",
            "min_overlap": MIN_OVERLAP,
        },
        "threshold": mismatched[int(0.95 * len(mismatched))] if mismatched else 1,
        "tfs": tfs, "mismatched": mismatched,
        "identity": identity(pairs, scores),
    }
    OUTPUT.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    rs = [tf["r"] for tf in tfs]
    print(f"{version and 'Tomtom ' + version}: {len(tfs)} TFs, median r {sorted(rs)[len(rs) // 2]:.2f}, "
          f"{sum(r > 0.5 for r in rs)} above 0.5; wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    sys.exit(main())
