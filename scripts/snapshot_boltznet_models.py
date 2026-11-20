#!/usr/bin/env python3
"""Snapshot the Lally et al. 2025 BoltzNet weight matrices from a download folder.

The paper's supplement has its models only as logo images. The authors serve
the numeric models on their BoltzNet site (the files their `boltznet` package
downloads):

  https://boltznet.bu.edu/static/data/tf_list.txt               the 124 TF names
  https://boltznet.bu.edu/static/data/models/<tf>_wm_bias.json  one model per TF

Download them into one folder, e.g.

  mkdir -p DIR && cd DIR
  curl -sO https://boltznet.bu.edu/static/data/tf_list.txt
  while read tf; do
    curl -sO "https://boltznet.bu.edu/static/data/models/${tf}_wm_bias.json"
    printf '%s\\t%s\\n' "$tf" "$(curl -sI "https://boltznet.bu.edu/static/data/models/${tf}_wm_bias.json" \\
      | grep -i '^last-modified' | cut -d' ' -f2- | tr -d '\\r')" >> last_modified.tsv
  done < tf_list.txt

then run

  python3 scripts/snapshot_boltznet_models.py DIR

It writes data/sources/chip_compendium/boltznet_weight_matrices.json: for each
TF, the 25 x (A, C, G, T) weight matrix exactly as served, its bias, the
SHA-256 of the downloaded file and the server's Last-Modified date. It checks
that the TF set equals the 124 modelled TFs in boltznet_models.csv. It prints
the output's SHA-256; data/sources/chip_compendium/manifest.json keeps the
reviewed value, which scripts/build_motifs.py checks.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CHIP_DIR = ROOT / "data" / "sources" / "chip_compendium"
OUTPUT = CHIP_DIR / "boltznet_weight_matrices.json"
WIDTH = 25


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    require(len(sys.argv) == 2, __doc__)
    folder = Path(sys.argv[1])
    tfs = [line.strip() for line in (folder / "tf_list.txt").read_text().splitlines() if line.strip()]
    with (CHIP_DIR / "boltznet_models.csv").open(newline="") as f:
        modelled = {r["tf"] for r in csv.DictReader(f) if r["has_model"] == "yes"}
    require(set(tfs) == modelled and len(tfs) == len(modelled),
            "The site's TF list differs from the modelled TFs in boltznet_models.csv.")
    dates = {}
    if (folder / "last_modified.tsv").exists():
        for line in (folder / "last_modified.tsv").read_text().splitlines():
            if line.strip():
                tf, date = line.split("\t", 1)
                dates[tf] = date.strip()

    out = {}
    for tf in sorted(tfs, key=str.casefold):
        raw = (folder / f"{tf}_wm_bias.json").read_bytes()
        model = json.loads(raw)
        wm = model["weight_matrix"]
        require(set(wm) == set("ACGT"), f"{tf}: unexpected bases {sorted(wm)}.")
        require(all(set(wm[b]) == {str(i) for i in range(WIDTH)} for b in "ACGT"),
                f"{tf}: expected positions 0-{WIDTH - 1}.")
        require(isinstance(model.get("bias"), (int, float)), f"{tf}: missing bias.")
        out[tf] = {
            "matrix": [[wm[b][str(i)] for b in "ACGT"] for i in range(WIDTH)],
            "bias": model["bias"],
            "source_sha256": hashlib.sha256(raw).hexdigest(),
            "last_modified": dates.get(tf, ""),
        }
    OUTPUT.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    print(f"Wrote {len(out)} weight matrices to {OUTPUT.relative_to(ROOT)}.")
    print(f"{hashlib.sha256((folder / 'tf_list.txt').read_bytes()).hexdigest()}  tf_list.txt")
    print(f"{hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}  {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
