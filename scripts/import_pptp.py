#!/usr/bin/env python3
"""Import the pinned Han et al. PPTP responses without changing existing IDs.

Run from any directory: python3 scripts/import_pptp.py
The source snapshot retains the authors' values, classifications and line numbers.
"""
import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = ROOT / "data/sources/pptp"
DATASET = ROOT / "data/dataset.csv"
PAPER = "https://doi.org/10.1038/s41467-023-41572-4"
FIELDS = ["perturbed_gene", "growth_medium", "response", "log2_fold_change",
          "log2_fold_change_vs_median", "neg_log10_q_value", "replicates"]
# Single-gene identities checked against RegulonDB TFSet (14.0). Keep the
# existing display names, and retain the paper's exact targeted gene separately.
# Do not merge individual subunit/toxin knockdowns into protein-complex rows.
ALIASES = {
    "cecR": "YbiH", "comR": "YcfQ", "dpiA": "CitB", "glnG": "NtrC",
    "hdfR": "YifA", "hprR": "YedW", "hypT": "YjiE", "lldR": "LctR",
    "nimR": "YeaM(NimR)", "pdeL": "YahA", "pyrR": "YpdB",
    "rclR": "YkgD(RclR)", "rob": "CbpB", "sutR": "YdcN(SutR)",
    "ttdR": "Dan",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def key(row):
    return row["perturbed_gene"], row["associated_genes"], row["growth_medium"]


def main():
    manifest = json.loads((SOURCE_DIR / "manifest.json").read_text())
    source = SOURCE_DIR / "significant_responses.csv"
    require(hashlib.sha256(source.read_bytes()).hexdigest() == manifest["selected_sha256"],
            "Source snapshot checksum differs from the reviewed import.")
    with source.open(newline="") as f:
        incoming = list(csv.DictReader(f))
    with DATASET.open(newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    ids = [r["record_id"] for r in rows]
    require(all(i.isdecimal() for i in ids) and len(ids) == len(set(ids)),
            "Existing record IDs must be unique positive integers.")
    names = {r["tf_name"].casefold(): r["tf_name"] for r in rows}
    require(all(name.casefold() in names for name in ALIASES.values()),
            "An existing alias target has changed; review the mapping.")
    existing = {}
    for row in rows:
        if row["method"] == "PPTP" and row["source"] == PAPER:
            require(key(row) not in existing, "Duplicate PPTP observation in dataset.")
            existing[key(row)] = row
    next_id = max([328] + [int(i) for i in ids]) + 1
    output_fields = fields + [f for f in FIELDS if f not in fields]
    counts = Counter()
    seen = set()
    additions = []
    for raw in incoming:
        gene = raw["tf_gene"]
        require(raw["class"] in {"1", "-1"} and float(raw["std"]) < 0.7,
                "Unexpected classification or replicate variability.")
        direction = int(raw["class"])
        # The upstream CSV is rounded to three decimals. Retain its class;
        # validate with rounding tolerance instead of reclassifying boundary rows.
        for field in ["FC", "FCOM"]:
            value = float(raw[field])
            require(math.isfinite(value) and direction * value > math.log2(1.7) - 0.00051,
                    "Response direction or fold-change threshold is inconsistent.")
        require(float(raw["#NAME?"]) >= 2 - 0.00051 and int(raw["n_rep"]) >= 2,
                "Unexpected significance score or replicate count.")
        display_name = names.get(gene.casefold(), ALIASES.get(gene, gene[0].upper() + gene[1:]))
        row = {f: "" for f in output_fields}
        row.update(tf_name=display_name, method="PPTP", strain="K-12 MG1655 derivative FR-E01",
                   associated_genes=raw["operon"], source=PAPER, perturbed_gene=gene,
                   growth_medium=raw["growth_medium"],
                   response="up after knockdown" if direction == 1 else "down after knockdown",
                   log2_fold_change=raw["FC"], log2_fold_change_vs_median=raw["FCOM"],
                   neg_log10_q_value=raw["#NAME?"], replicates=raw["n_rep"])
        identity = key(row)
        require(identity not in seen, "Duplicate source TF/promoter/medium observation.")
        seen.add(identity)
        counts[raw["source_file"]] += 1
        if identity in existing:
            old = existing[identity]
            require(all(old.get(f, "") == row[f] for f in output_fields if f != "record_id"),
                    "Existing PPTP observation differs; review instead of overwriting it.")
        else:
            row["record_id"] = str(next_id)
            next_id += 1
            additions.append(row)
    require(set(existing) <= seen, "An existing PPTP record is absent from this source.")
    require(dict(counts) == {f["file"]: f["selected_rows"] for f in manifest["files"]},
            "Imported counts differ from the source manifest.")
    if additions:
        buffer = io.StringIO(newline="")
        writer = csv.DictWriter(buffer, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(rows + additions)
        DATASET.write_bytes(buffer.getvalue().encode("utf-8"))
    print(f"Added {len(additions)} PPTP records; {len(incoming)} verified source observations; "
          f"{len(rows) + len(additions)} total records.")


if __name__ == "__main__":
    main()
