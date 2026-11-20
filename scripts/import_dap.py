#!/usr/bin/env python3
"""Import the Baumgart et al. 2021 E. coli DAP-seq TF list into dataset.csv.

Source: data/sources/dap/dap_tf_list.csv — the E. coli subset of
Supplementary Table 2 of Baumgart et al., Nat Methods 18, 1499-1510 (2021):
214 E. coli K-12 MG1655 TFs assayed by in vitro biotin DAP-seq, 107 of them
successful. One row is added per successful TF (method DAP). Binding-region
coordinates are not imported here; they live in the paper's Supplementary
Data 1 and in RegulonDB.

Repeated runs are idempotent: an existing (tf_name, DAP, this source) row
is never duplicated or overwritten. New record IDs continue above the
current maximum.
"""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = ROOT / "data" / "sources" / "dap"
DATASET = ROOT / "data" / "dataset.csv"
PAPER = "https://doi.org/10.1038/s41592-021-01312-2"

# Paper gene name -> existing table entry. Verified one-to-one synonyms
# (RegulonDB canonical names / renames); same policy as import_pptp.py
# and import_chip_compendium.py.
ALIASES = {
    "alsR": "RpiR", "cecR": "YbiH", "dpiA": "CitB", "glnG": "NtrC",
    "hypT": "YjiE", "lldR": "LctR", "nimR": "YeaM(NimR)", "ttdR": "Dan",
    "yiaJ": "YiaJ(BdcA)",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    manifest = json.loads((SOURCE_DIR / "manifest.json").read_text())
    source = SOURCE_DIR / "dap_tf_list.csv"
    require(hashlib.sha256(source.read_bytes()).hexdigest() == manifest["dap_tf_list_sha256"],
            "Source snapshot checksum differs from the reviewed import.")
    with source.open(newline="") as f:
        incoming = list(csv.DictReader(f))
    require(len(incoming) == manifest["tf_count_attempted"], "Unexpected TF count in snapshot.")
    successful = [r for r in incoming if r["successful"] == "yes"]
    require(len(successful) == manifest["tf_count_successful"], "Unexpected successful-TF count.")

    with DATASET.open(newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    names = {r["tf_name"].casefold() for r in rows if r["method"] in ("", "SELEX")}
    require(all(name.casefold() in names for name in ALIASES.values()),
            "An existing alias target has changed; review the mapping.")

    existing = {r["tf_name"].casefold() for r in rows
                if r["method"] == "DAP" and r["source"] == PAPER}
    ids = [int(r["record_id"]) for r in rows if r["record_id"].isdigit()]
    next_id = max(ids) + 1

    additions = []
    for raw in successful:
        gene = raw["tf"]
        name = ALIASES.get(gene, gene[0].upper() + gene[1:])
        if name.casefold() in existing:
            continue
        row = {f: "" for f in fields}
        row.update({
            "tf_name": name,
            "record_id": str(next_id),
            "method": "DAP",
            "sequence_strain": "K-12 MG1655",
            "strain": "K-12 MG1655",
            "source": PAPER,
        })
        additions.append(row)
        next_id += 1

    if not additions:
        print("Nothing to import; all DAP-seq TFs are already present.")
        return
    with DATASET.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writerows(additions)
    in_table = sum(1 for r in additions if r["tf_name"].casefold() in names)
    print(f"Imported {len(additions)} DAP records "
          f"({in_table} match existing TF-table entries, "
          f"{len(additions) - in_table} are new names outside the table).")


if __name__ == "__main__":
    main()
