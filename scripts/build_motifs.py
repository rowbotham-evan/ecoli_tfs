#!/usr/bin/env python3
"""Build data/motifs.csv and data/motif_matrices.json from the motif snapshots.

data/motifs.csv has one row per (TF, paper): which TFs a source paper computed
a binding motif for from its experiment and published as a numeric frequency
matrix, the motif's type, where it is published, its width, site count and
E-value. data/motif_matrices.json holds the numeric matrices, keyed by TF name
and paper, for drawing logos: {"kind": "frequency"}.

Source (regenerate with scripts/snapshot_motif_supplements.py):
  - data/sources/dap/: Baumgart et al. 2021, Supplementary Data 2. One MEME
    motif per E. coli TF (124): dap_motifs.csv and dap_motif_matrices.json.

Left out, because the paper gives no numeric frequency matrix:
  - Lally et al. 2025: BoltzNet energy matrices for 124 ChIP-seq TFs, with the
    frequency matrices only as logo images. Its snapshot stays in
    data/sources/chip_compendium/.
  - Fitzgerald et al. 2023: the PhoB MEME motif (Fig. 2B), published as a logo
    and its 59 sites. The sites stay in dataset.csv.

Every snapshot is checked against the SHA-256 in its manifest. TF names follow
the importers' alias map, so motif rows use the same table names as the DAP
records. Both outputs are rewritten on every run.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.dont_write_bytecode = True  # importing the alias maps should not leave __pycache__ in scripts/
sys.path.insert(0, str(ROOT / "scripts"))
from import_dap import ALIASES as DAP_ALIASES  # noqa: E402

DATASET = ROOT / "data" / "dataset.csv"
OUTPUT = ROOT / "data" / "motifs.csv"
MATRICES = ROOT / "data" / "motif_matrices.json"
DAP_DIR = ROOT / "data" / "sources" / "dap"
FIELDS = ["tf_name", "method", "paper", "source", "motif_type", "published_as",
          "paper_gene", "width", "sites", "e_value", "supplement_url", "supplement_file",
          "logo_images", "matrix_url"]
PAPER_ORDER = {"Baumgart et al. 2021": 0}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest(source_dir):
    return json.loads((source_dir / "manifest.json").read_text())


def checked(source_dir, filename, checksum_key):
    path = source_dir / filename
    require(sha256(path) == manifest(source_dir)[checksum_key],
            f"{path.relative_to(ROOT)} differs from the reviewed snapshot.")
    return path


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def table_name(gene, aliases):
    return aliases.get(gene, gene[0].upper() + gene[1:])


def main():
    out, matrices = [], {}

    dap = manifest(DAP_DIR)
    dap_rows = read_csv(checked(DAP_DIR, "dap_motifs.csv", "dap_motifs_sha256"))
    require(len(dap_rows) == dap["motif_count_ecoli"], "Unexpected DAP-seq motif count.")
    dap_matrices = json.loads(checked(DAP_DIR, "dap_motif_matrices.json", "dap_motif_matrices_sha256").read_text())
    require(set(dap_matrices) == {r["tf"] for r in dap_rows}, "DAP-seq matrices and motif list disagree.")
    for r in dap_rows:
        name = table_name(r["tf"], DAP_ALIASES)
        out.append({
            "tf_name": name, "method": "DAP", "paper": "Baumgart et al. 2021",
            "source": "https://doi.org/10.1038/s41592-021-01312-2",
            "motif_type": "MEME frequency matrix",
            "published_as": "MEME motif file (Supplementary Data 2)",
            "paper_gene": r["tf"], "width": r["width"], "sites": r["sites"], "e_value": r["e_value"],
            "supplement_url": dap["motifs_url"], "supplement_file": r["motif_file"], "logo_images": "",
            "matrix_url": "",
        })
        matrices.setdefault(name, {})["Baumgart et al. 2021"] = {
            "kind": "frequency", "matrix": dap_matrices[r["tf"]]["matrix"]}

    rows = read_csv(DATASET)

    keys = [(r["paper"], r["tf_name"].casefold()) for r in out]
    require(len(keys) == len(set(keys)), "Two motifs map to the same TF within one paper.")
    out.sort(key=lambda r: (PAPER_ORDER[r["paper"]], r["tf_name"].casefold()))
    with OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)
    ordered = {n: matrices[n] for n in sorted(matrices, key=str.casefold)}
    MATRICES.write_text(json.dumps(ordered, separators=(",", ":")) + "\n")

    base = {r["tf_name"].casefold() for r in rows if r["method"] in ("", "SELEX")}
    tfs = {r["tf_name"] for r in out}
    outside = sorted(n for n in tfs if n.casefold() not in base)
    print(f"Wrote {len(out)} motif rows for {len(tfs)} distinct TF names to {OUTPUT.relative_to(ROOT)}, "
          f"and matrices for {len(ordered)} TFs ({sum(len(v) for v in ordered.values())} motifs) to {MATRICES.relative_to(ROOT)}.")
    print(f"{len(tfs) - len(outside)} are TF-table entries; {len(outside)} are not: {', '.join(outside)}")


if __name__ == "__main__":
    main()
