#!/usr/bin/env python3
"""Import the Lally et al. 2025 ChIP-seq compendium TF list into dataset.csv.

Source: data/sources/chip_compendium/tf_chip_stats.csv — Supplementary Data 2
of Lally et al., Nat Commun 16, 4255 (2025), the 139 E. coli TFs whose
ChIP-seq experiments passed the paper's QC. One row is added per TF
(method CHIP). Binding-region coordinates are not imported here; per-TF
union-region counts stay in the source snapshot, and the full region data
live in GEO GSE268698 and RegulonDB.

Repeated runs are idempotent: an existing (tf_name, CHIP, this source) row
is never duplicated or overwritten. New record IDs continue above the
current maximum.
"""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = ROOT / "data" / "sources" / "chip_compendium"
DATASET = ROOT / "data" / "dataset.csv"
PAPER = "https://doi.org/10.1038/s41467-025-58862-8"

# Paper gene name -> existing table entry. Verified one-to-one synonyms
# (RegulonDB canonical names / renames); same policy as import_pptp.py.
ALIASES = {
    "dgsA": "Mlc", "glnG": "NtrC", "hdfR": "YifA", "lldR": "LctR",
    "ycgE": "BluR", "yciT": "YciT(DeoT)", "ydcN": "YdcN(SutR)",
    "yeaM": "YeaM(NimR)", "yihW": "CsqR", "ykgD": "YkgD(RclR)",
    "yncC": "McbR", "ygiT": "MqsA",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    manifest = json.loads((SOURCE_DIR / "manifest.json").read_text())
    source = SOURCE_DIR / "tf_chip_stats.csv"
    require(hashlib.sha256(source.read_bytes()).hexdigest() == manifest["tf_chip_stats_sha256"],
            "Source snapshot checksum differs from the reviewed import.")
    with source.open(newline="") as f:
        incoming = list(csv.DictReader(f))
    require(len(incoming) == manifest["tf_count"], "Unexpected TF count in snapshot.")

    with DATASET.open(newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    names = {r["tf_name"].casefold() for r in rows if r["method"] in ("", "SELEX")}
    require(all(name.casefold() in names for name in ALIASES.values()),
            "An existing alias target has changed; review the mapping.")

    existing = {r["tf_name"].casefold() for r in rows
                if r["method"] == "CHIP" and r["source"] == PAPER}
    ids = [int(r["record_id"]) for r in rows if r["record_id"].isdigit()]
    next_id = max(ids) + 1

    additions = []
    for raw in incoming:
        gene = raw["tf"]
        name = ALIASES.get(gene, gene[0].upper() + gene[1:])
        if name.casefold() in existing:
            continue
        modes = [label for key, label in [
            ("native_spa", "native SPA"), ("native_flag", "native 3xFLAG"),
            ("native_v5", "native V5"), ("inducible", "inducible pBAD24-3xFLAG")]
            if raw[key] == "yes"]
        require(modes, f"No expression mode marked for {gene}.")
        row = {f: "" for f in fields}
        row.update({
            "tf_name": name,
            "record_id": str(next_id),
            "method": "CHIP",
            "sequence_strain": "K-12 MG1655",
            "strain": f"K-12 MG1655 ({', '.join(modes)})",
            "growth_medium": "M9 glycerol (0.4%)",
            "source": PAPER,
        })
        additions.append(row)
        next_id += 1

    if not additions:
        print("Nothing to import; all compendium TFs are already present.")
        return
    with DATASET.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writerows(additions)
    in_table = sum(1 for r in additions if r["tf_name"].casefold() in names)
    print(f"Imported {len(additions)} CHIP records "
          f"({in_table} match existing TF-table entries, "
          f"{len(additions) - in_table} are new names outside the table).")


if __name__ == "__main__":
    main()
