#!/usr/bin/env python3
"""Export the position frequency matrix of every TF motif published as one.

    python3 scripts/export_frequency_matrices.py

Writes two files, rewritten in full on each run:

  data/motif_frequency_matrices.csv   one row per motif position: A, C, G, T probabilities
                                      and the base counts behind them
  data/motif_frequency_matrices.meme  the same matrices in MEME format, for Tomtom, FIMO etc.

Only matrices a paper computed from its experiment and published as numbers are exported,
as is: the DAP-seq MEME frequency matrices of Baumgart et al. 2021. These are the only
motifs scripts/build_motifs.py keeps; see its docstring for the ones left out.

The CSV also carries the base counts behind each frequency (frequency x the motif's number
of sites, in count_A..count_T).

Uses only the standard library.
"""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MATRICES = ROOT / "data" / "motif_matrices.json"
MOTIFS = ROOT / "data" / "motifs.csv"
OUT_CSV = ROOT / "data" / "motif_frequency_matrices.csv"
OUT_MEME = ROOT / "data" / "motif_frequency_matrices.meme"
BASES = "ACGT"
BASIS = "published frequency matrix"


def main():
    matrices = json.loads(MATRICES.read_text())
    with MOTIFS.open(newline="") as handle:
        info = {(r["tf_name"], r["paper"]): r for r in csv.DictReader(handle)}
    table, meme = [], ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
                        "Background letter frequencies", "A 0.25 C 0.25 G 0.25 T 0.25", ""]
    motifs, tfs, skipped = 0, set(), len(info)
    for tf in sorted(matrices, key=str.lower):
        for paper, entry in matrices[tf].items():
            if entry["kind"] != "frequency":
                continue  # an energy matrix: its frequencies would be derived here, not published
            meta = info.get((tf, paper), {})
            ppm = [[v / sum(row) for v in row] for row in entry["matrix"]]
            sites = int(meta["sites"]) if meta.get("sites") else None
            for position, row in enumerate(ppm, 1):
                counts = [round(v * sites) for v in row] if sites else [""] * 4
                table.append([tf, meta.get("method", ""), paper, BASIS, len(ppm), meta.get("sites", ""), position]
                             + [f"{v:.6f}" for v in row] + counts)
            extra = f" nsites= {meta['sites']} E= {meta.get('e_value') or 0}" if meta.get("sites") else ""
            meme += [f"MOTIF {tf}_{meta.get('method', 'motif')} {tf} ({paper}; {BASIS})",
                     f"letter-probability matrix: alength= 4 w= {len(ppm)}{extra}"]
            meme += [" ".join(f"{v:.6f}" for v in row) for row in ppm] + [""]
            motifs += 1
            tfs.add(tf)
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["tf_name", "method", "paper", "basis", "width", "sites", "position", *BASES,
                         *(f"count_{b}" for b in BASES)])
        writer.writerows(table)
    OUT_MEME.write_text("\n".join(meme))
    print(f"{motifs} matrices for {len(tfs)} TFs, {len(table)} positions "
          f"({skipped - motifs} motifs without a published frequency matrix left out)")
    print(f"  {OUT_CSV.relative_to(ROOT)}\n  {OUT_MEME.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
