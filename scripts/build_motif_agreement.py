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

Tomtom is run once, every two-motif TF's DAP-seq motif against every two-motif TF's
second motif, with
`-dist pearson -min-overlap MIN_OVERLAP -thresh 1.0`, so it reports its best
alignment (offset and orientation) for all pairs. The Pearson correlation is then
computed here from the untouched PPMs over the aligned columns. A TF's own pair
gives its score; the pairs of different two-motif TFs give the "mismatched"
reference distribution, since choosing the best alignment raises every score.

The same all-pairs scores answer a second question: does a TF's identity survive the
change of method? For every two-motif TF, the output's "identity" block gives the
rank of its own partner among the two-motif TFs' motifs of the other method (rank 1 =
retrieved), in both directions. Because related TFs recognize similar sequences, each TF's DNA-binding
family (the UniProt domain note in data/sources/uniprot/dbd_uniprot.tsv) is recorded
too, and the retrieval is repeated among same-family candidates only. The block
also holds the cross-method matrix for the heatmap.

A second file, data/motif_grid.json, is for the "all motifs" page: every DAP-seq motif
against every ChIP-seq motif, including the TFs that only one method covers. It
comes from a second Tomtom run over all motifs. Tomtom picks each pair's alignment
using statistics of the whole target set, so a few weak pairs score differently
there than in the two-motif run; every TF's own pair scores the same.

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
GRID = ROOT / "data" / "motif_grid.json"
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


def load_motifs():
    """Every TF's motifs as PPMs: (daps, seconds), two lists of {"name", "ppm", ...}.

    daps are the DAP-seq motifs; seconds are the ChIP-seq motifs (with "other" naming
    their kind). A TF with both appears in both lists under the same name."""
    matrices = json.loads(MATRICES.read_text())
    sites = {}
    with DATASET.open(newline="") as handle:
        for row in csv.DictReader(handle):
            seq = (row.get("dna_binding_sequence") or "").strip().upper()
            if seq:
                sites.setdefault(row["tf_name"], []).append(seq)
    daps, seconds = [], []
    for name, papers in matrices.items():
        entries = list(papers.values())
        freq = next((e for e in entries if e["kind"] == "frequency"), None)
        energy = next((e for e in entries if e["kind"] == "energy"), None)
        if freq:
            daps.append({"name": name, "ppm": normalize(freq["matrix"])})
        if energy:
            seconds.append({"name": name, "ppm": boltzmann(energy["matrix"]), "other": "ChIP-seq energy matrix (Lally et al.)"})
        elif freq:
            # A TF whose second motif is a set of aligned sites in the dataset (PhoB).
            seqs = sites.get(name, [])
            if len(seqs) >= MIN_SITES and len({len(s) for s in seqs}) == 1 and not any(set(s) - set(BASES) for s in seqs):
                ppm = normalize([[sum(s[i] == b for s in seqs) for b in BASES] for i in range(len(seqs[0]))])
                seconds.append({"name": name, "ppm": ppm, "other": "ChIP-seq aligned sites (Fitzgerald et al.)"})
    require(daps and seconds, "No motifs to compare.")
    return daps, seconds


def load_families():
    """TF name -> DNA-binding family, for single-protein TFs whose UniProt domain names one."""
    notes = {}
    with UNIPROT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            notes.setdefault(row["tf_name"], []).append(row["dbd_note"].strip())
    return {name: found[0].replace("HTH ", "").replace("-type", "") for name, found in notes.items()
            if len(found) == 1 and found[0] not in NOT_A_FAMILY}


def identity(daps, seconds, paired, scores):
    """Cross-method retrieval: is a TF's own partner its most similar motif?

    scores[(q, t)] is the similarity of daps[q] to seconds[t]; paired lists the
    (q, t) of every TF that has both. Each motif is ranked against all motifs of the
    other method, not only those of the paired TFs."""
    families = load_families()
    family = lambda motif: families.get(motif["name"])
    names = [daps[q]["name"] for q, _ in paired]
    own_family = [family(daps[q]) for q, _ in paired]
    sizes = {f: own_family.count(f) for f in set(own_family) if f}
    # Families with several paired members first (largest first), so they form blocks in the matrix.
    order = sorted(range(len(paired)), key=lambda i: (
        -sizes.get(own_family[i], 0) if sizes.get(own_family[i], 0) > 1 else 1, own_family[i] or "~", names[i].lower()))

    def rank(own, candidates, score):
        """1 + the number of candidates scoring above the TF's own partner."""
        return 1 + sum(score(c) > score(own) for c in candidates if c != own)

    tfs = []
    for i in order:
        q, t = paired[i]
        forward = lambda c: scores[(q, c)]      # this TF's DAP-seq motif against a second motif
        backward = lambda c: scores[(c, t)]     # a DAP-seq motif against this TF's second motif
        kin_seconds = [c for c in range(len(seconds)) if c != t and own_family[i] and family(seconds[c]) == own_family[i]]
        kin_daps = [c for c in range(len(daps)) if c != q and own_family[i] and family(daps[c]) == own_family[i]]
        best_second, best_dap = max(range(len(seconds)), key=forward), max(range(len(daps)), key=backward)
        tfs.append({
            "name": names[i], "family": own_family[i] if sizes.get(own_family[i], 0) > 1 else None,
            "rank_forward": rank(t, range(len(seconds)), forward),
            "rank_backward": rank(q, range(len(daps)), backward),
            "best_forward": seconds[best_second]["name"], "best_backward": daps[best_dap]["name"],
            # Whether the motif retrieved instead of the TF's own belongs to a relative.
            "best_forward_related": best_second in kin_seconds, "best_backward_related": best_dap in kin_daps,
            # The same question among same-family TFs only (None when the TF has no relative in that pool).
            "family_mates": len(kin_seconds),
            "family_rank_forward": rank(t, [t] + kin_seconds, forward) if kin_seconds else None,
            "family_mates_backward": len(kin_daps),
            "family_rank_backward": rank(q, [q] + kin_daps, backward) if kin_daps else None,
        })
    block = [paired[i] for i in order]
    mean = lambda values: round(sum(values) / len(values), 4) if values else None
    related = lambda i, j: i != j and own_family[i] and own_family[i] == own_family[j]
    cross = [(i, j) for i in range(len(paired)) for j in range(len(paired)) if i != j]
    return {
        "pool_forward": len(seconds), "pool_backward": len(daps),   # candidates each motif is ranked against
        "tfs": tfs,
        # The paired TFs only (rows: DAP-seq motif; columns: second motif), for the heatmap.
        "matrix": [[round(scores[(q, t)], 3) for _, t in block] for q, _ in block],
        "mean_own": mean([scores[pair] for pair in paired]),
        "mean_same_family": mean([scores[(paired[i][0], paired[j][1])] for i, j in cross if related(i, j)]),
        "mean_other_family": mean([scores[(paired[i][0], paired[j][1])] for i, j in cross if not related(i, j)]),
    }


def write_meme(path, motifs):
    lines = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25", ""]
    for motif_id, ppm in motifs:
        lines += [f"MOTIF {motif_id}", f"letter-probability matrix: alength= 4 w= {len(ppm)} nsites= 20 E= 0"]
        lines += [" ".join(f"{v:.6f}" for v in row) for row in ppm] + [""]
    path.write_text("\n".join(lines))


def run_tomtom(tomtom, daps, seconds):
    """Tomtom's best alignment for every (DAP-seq motif, second motif) pair."""
    args = ["-dist", "pearson", "-min-overlap", str(MIN_OVERLAP), "-thresh", "1.0", "-text", "-verbosity", "1"]
    with tempfile.TemporaryDirectory() as tmp:
        query, target = Path(tmp) / "dap.meme", Path(tmp) / "second.meme"
        # Tomtom ids cannot hold spaces or parentheses, so motifs are numbered.
        write_meme(query, [(f"q{i}", m["ppm"]) for i, m in enumerate(daps)])
        write_meme(target, [(f"t{i}", m["ppm"]) for i, m in enumerate(seconds)])
        out = subprocess.run([tomtom, *args, str(query), str(target)], check=True, capture_output=True, text=True).stdout
    rows = [r for r in csv.DictReader((line for line in out.splitlines() if line and not line.startswith("#")), delimiter="\t")]
    hits = {}
    for r in rows:
        hits[(int(r["Query_ID"][1:]), int(r["Target_ID"][1:]))] = {
            "offset": int(r["Optimal_offset"]), "overlap": int(r["Overlap"]),
            "reversed": r["Orientation"] == "-", "p": float(r["p-value"]),
        }
    require(len(hits) == len(daps) * len(seconds), f"Tomtom reported {len(hits)} of {len(daps) * len(seconds)} pairs.")
    return hits, args


def overlap_columns(query_len, target_len, shift):
    """Query columns i that sit over target column i + shift."""
    return [i for i in range(query_len) if 0 <= i + shift < target_len]


def offset_sign(daps, seconds, hits):
    """Whether Tomtom's offset is added to or subtracted from the query column.

    Settled from the data: only one sign reproduces Tomtom's own overlap counts."""
    fits = [sign for sign in (1, -1) if all(
        len(overlap_columns(len(daps[q]["ppm"]), len(seconds[t]["ppm"]), sign * h["offset"])) == h["overlap"]
        for (q, t), h in hits.items())]
    require(len(fits) == 1, f"Could not settle Tomtom's offset convention (candidates: {fits}).")
    return fits[0]


def score_pairs(tomtom, daps, seconds):
    """Align every DAP-seq motif with every second motif and score each pair.

    Returns Tomtom's hits, the Pearson score of every (q, t), the aligned PPMs behind
    it (dap, second, shift, columns) and the Tomtom options used."""
    hits, args = run_tomtom(tomtom, daps, seconds)
    sign = offset_sign(daps, seconds, hits)
    scores, aligned = {}, {}
    for (q, t), hit in hits.items():
        dap = daps[q]["ppm"]
        second = reverse_complement(seconds[t]["ppm"]) if hit["reversed"] else seconds[t]["ppm"]
        shift = sign * hit["offset"]
        columns = overlap_columns(len(dap), len(second), shift)
        scores[(q, t)] = pearson([v for i in columns for v in dap[i]], [v for i in columns for v in second[i + shift]])
        aligned[(q, t)] = (dap, second, shift, columns)
    return hits, scores, aligned, args


def write_grid(tomtom, every_dap, every_second, block):
    """data/motif_grid.json: every DAP-seq motif (rows) against every second motif (columns).

    block lists the two-motif TFs in the order of the comparison's matrix (name and
    family); they come first on both axes, so their own pairs are the diagonal of the
    top-left square. The TFs with one motif follow, by name."""
    _, scores, _, _ = score_pairs(tomtom, every_dap, every_second)
    first = [tf["name"] for tf in block]
    rest = lambda motifs: sorted((m["name"] for m in motifs if m["name"] not in first), key=str.lower)
    rows, cols = first + rest(every_dap), first + rest(every_second)
    dap_at = {m["name"]: q for q, m in enumerate(every_dap)}
    second_at = {m["name"]: t for t, m in enumerate(every_second)}
    GRID.write_text(json.dumps({
        "rows": rows, "cols": cols, "paired": len(first),
        "families": [tf["family"] for tf in block],   # of the paired TFs, in order; None when the TF has no relative among them
        "matrix": [[round(scores[(dap_at[r], second_at[c])], 2) for c in cols] for r in rows],
    }, separators=(",", ":")) + "\n")
    return len(rows), len(cols)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--tomtom", default=shutil.which("tomtom"), help="path to the tomtom program (default: from PATH)")
    tomtom = parser.parse_args().tomtom
    require(tomtom, "tomtom not found: install the MEME Suite or pass --tomtom PATH.")
    version = subprocess.run([tomtom, "-version"], check=True, capture_output=True, text=True).stdout.strip()

    every_dap, every_second = load_motifs()
    # Only the TFs with both motifs take part in the comparison: each one's motif is compared
    # with its own partner and with the other two-motif TFs' motifs.
    both = {m["name"] for m in every_dap} & {m["name"] for m in every_second}
    daps, seconds = [m for m in every_dap if m["name"] in both], [m for m in every_second if m["name"] in both]
    hits, scores, aligned, args = score_pairs(tomtom, daps, seconds)

    # The TFs with both motifs, in the order of the matrices file.
    second_of = {m["name"]: t for t, m in enumerate(seconds)}
    paired = [(q, second_of[m["name"]]) for q, m in enumerate(daps) if m["name"] in second_of]
    require(paired, "No TF has two motifs to compare.")

    mismatched = sorted(round(scores[(q, t)], 4) for q, _ in paired for _, t in paired if daps[q]["name"] != seconds[t]["name"])
    tfs = []
    for q, t in paired:
        dap, second, shift, columns = aligned[(q, t)]
        r = round(scores[(q, t)], 4)
        tfs.append({
            "name": daps[q]["name"], "other": seconds[t]["other"], "r": r,
            "overlap": hits[(q, t)]["overlap"], "reversed": hits[(q, t)]["reversed"], "tomtom_p": hits[(q, t)]["p"],
            "p": round(sum(v >= r for v in mismatched) / len(mismatched), 4) if mismatched else 1,
            # The aligned cells: [DAP-seq probability, other motif's probability, base index, DAP-seq position].
            "points": [[round(dap[i][k], 4), round(second[i + shift][k], 4), k, i + 1] for i in columns for k in range(4)],
        })
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
        "identity": identity(daps, seconds, paired, scores),
    }
    OUTPUT.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    grid_rows, grid_cols = write_grid(tomtom, every_dap, every_second, result["identity"]["tfs"])
    rs = [tf["r"] for tf in tfs]
    found = result["identity"]["tfs"]
    print(f"Tomtom {version}: {len(tfs)} TFs with two motifs, median r {sorted(rs)[len(rs) // 2]:.2f}, {sum(r > 0.5 for r in rs)} above 0.5. "
          f"Own partner retrieved for {sum(t['rank_forward'] == 1 for t in found)} (among {len(seconds)} ChIP-seq motifs) and "
          f"{sum(t['rank_backward'] == 1 for t in found)} (among {len(daps)} DAP-seq motifs); wrote {OUTPUT.relative_to(ROOT)} "
          f"and {GRID.relative_to(ROOT)} ({grid_rows} x {grid_cols})")


if __name__ == "__main__":
    sys.exit(main())
