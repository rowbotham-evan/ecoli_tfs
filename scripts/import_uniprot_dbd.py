#!/usr/bin/env python3
"""Write UniProt DNA-binding-domain sequences into the TF table.

Fills `tf_binding_domain` (and `sequence_strain`) on each TF-table base row from
the reviewed UniProt snapshot in data/sources/uniprot/.

Only spans that are a real domain unit are imported - a DNA_BIND feature at
least 40 aa, or the annotated HTH DOMAIN. Short recognition-helix motifs and
full-chain fallbacks are left out, and composite entries (IHF, FlhDC, ...) are
skipped because one cell cannot hold two subunit sequences; all of these are
listed in the manifest.

Repeated runs are idempotent. An existing sequence is only replaced when the new
span strictly contains it - the table held bare ~20 aa recognition motifs, and
UniProt gives the full domain around them. Anything else is reported as a
conflict and left untouched.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = ROOT / "data" / "sources" / "uniprot"
DATASET = ROOT / "data" / "dataset.csv"

IMPORT_SPANS = {
    "DNA_BIND",            # annotated DNA-binding region, >= 40 aa
    "DOMAIN(HTH)",         # annotated helix-turn-helix domain
    "DNA_BIND(windowed)",  # bare H-T-H motif widened to a 60 aa window
    "FULL_CHAIN(compact)", # unannotated but <= 200 aa: a single-domain binder
}
STRAIN = "K-12 MG1655"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    manifest = json.loads((SOURCE_DIR / "manifest.json").read_text())
    for name, meta in manifest["files"].items():
        digest = hashlib.sha256((SOURCE_DIR / name).read_bytes()).hexdigest()
        require(digest == meta["sha256"],
                f"{name} checksum does not match the manifest; re-snapshot the source")
    require(sorted(manifest["imported_span_types"]) == sorted(IMPORT_SPANS),
            "manifest span types disagree with this importer")

    with (SOURCE_DIR / "dbd_uniprot.tsv").open(newline="", encoding="utf-8") as fh:
        dbd_rows = list(csv.DictReader(fh, delimiter="\t"))
    require(len(dbd_rows) == manifest["files"]["dbd_uniprot.tsv"]["records"],
            "dbd_uniprot.tsv row count disagrees with the manifest")

    # tf_name -> (sequence, accession, span) for the importable spans only.
    chosen = {}
    for r in dbd_rows:
        if r["dbd_source"] not in IMPORT_SPANS:
            continue
        if ":" in r["label"]:            # composite subunit
            continue
        seq = r["dbd_sequence"].strip().upper()
        if not seq:
            continue
        require(r["tf_name"] not in chosen,
                f"two importable spans for {r['tf_name']}; the source is ambiguous")
        chosen[r["tf_name"]] = (seq, r["accession"], f"{r['dbd_start']}-{r['dbd_end']}")

    with DATASET.open(newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = list(reader)
    col = {name: i for i, name in enumerate(header)}
    for needed in ("tf_name", "method", "tf_binding_domain", "sequence_strain"):
        require(needed in col, f"dataset.csv is missing the {needed} column")

    # One row per TF name carries the sequence: its TF-table base row where one
    # exists, otherwise the first record bearing that name (composite subunits
    # and paper-only names have no base row, but the site shows their rows and
    # the domain is a property of the protein either way).
    target = {}
    for i, row in enumerate(rows):
        tf = row[col["tf_name"]].strip()
        if tf not in chosen:
            continue
        is_base = row[col["method"]].strip() in ("", "SELEX")
        cur = target.get(tf)
        if cur is None or (is_base and not cur[1]):
            target[tf] = (i, is_base)

    filled, extended, unchanged, conflicts = [], [], [], []
    for tf, (i, _is_base) in target.items():
        row = rows[i]
        new_seq, accession, span = chosen[tf]
        current = row[col["tf_binding_domain"]].strip().upper()
        if current == new_seq:
            unchanged.append(tf)
            continue
        if current and current not in new_seq:
            conflicts.append((tf, accession, len(current), len(new_seq)))
            continue
        row[col["tf_binding_domain"]] = new_seq
        row[col["sequence_strain"]] = STRAIN
        (extended if current else filled).append(tf)

    require(not conflicts,
            "existing sequences do not fit inside the UniProt span: "
            + ", ".join(f"{t} ({a}, {o}aa vs {n}aa)" for t, a, o, n in conflicts))

    if not filled and not extended:
        print(f"No change. {len(unchanged)} TF entries already carry the UniProt DBD.")
        return 0

    with DATASET.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\r\n")
        writer.writerow(header)
        writer.writerows(rows)

    print(f"UniProt release {manifest['uniprot_release']} "
          f"({manifest['uniprot_release_date']})")
    print(f"  filled empty tf_binding_domain : {len(filled)}")
    print(f"  replaced a shorter fragment    : {len(extended)}")
    print(f"  already correct                : {len(unchanged)}")
    print(f"  TF entries with a DBD now      : {len(filled) + len(extended) + len(unchanged)}")
    skipped = len(dbd_rows) - len(chosen)
    print(f"  source rows not imported       : {skipped} "
          f"(short motifs, full-chain fallbacks, composite subunits)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as exc:
        sys.exit(f"error: {exc}")
