#!/usr/bin/env python3
"""Build data/regulondb_sites.csv from the RegulonDB TF-RISet snapshot.

Source: data/sources/regulondb/TF-RISet.tsv — RegulonDB's curated set of
regulatory interactions between a transcription factor and the promoter,
gene or transcription unit it regulates (release and checksum in
manifest.json). Each interaction may name a TF binding site with genome
coordinates (U00096.3) and its sequence.

Output: one row per TF and binding site that has a sequence. A site that
takes part in several interactions (several promoters or genes) is listed
once, with its functions, targets, evidence codes and PMIDs merged.

dataset.csv is not modified. The output is rewritten in full on every run,
so repeated runs give the same file.
"""

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIR = ROOT / "data" / "sources" / "regulondb"
DATASET = ROOT / "data" / "dataset.csv"
OUTPUT = ROOT / "data" / "regulondb_sites.csv"

# RegulonDB TF name -> TF-table entry, for names that differ. These are the
# synonyms verified against RegulonDB's TFSet download. RegulonDB's YciT is
# listed in the table twice (DeoT and YciT(DeoT)); its sites go to the entry
# matched by name.
ALIASES = {
    "Rob": "CbpB", "DpiA": "CitB", "YciT": "YciT(DeoT)", "UidR": "GusR",
    "H-NS": "HNS", "LldR": "LctR", "AaeR": "QseA", "AlsR": "RpiR",
    "XynR": "YagI", "PdeL": "YahA", "DecR": "YbaO", "CecR": "YbiH",
    "ComR": "YcfQ", "SutR": "YdcN(SutR)", "NimR": "YeaM(NimR)",
    "HprR": "YedW", "SrsR": "YgfI", "PlaR": "YiaJ(BdcA)", "HdfR": "YifA",
    "HypT": "YjiE", "RclR": "YkgD(RclR)", "PtrR": "YneJ", "PyrR": "YpdB",
}

FIELDS = [
    "tf_name", "regulondb_tf_name", "regulondb_tf_id", "site_id",
    "left", "right", "length", "sequence", "sequence_with_flanks",
    "function", "target_genes", "promoters", "confidence", "experimental",
    "site_evidence", "site_pmids", "interaction_ids",
]
CONFIDENCE = {"C": "confirmed", "S": "strong", "W": "weak", "?": "unknown"}
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def merged(values):
    """Distinct non-empty values of '|'-separated cells, in first-seen order."""
    seen = []
    for value in values:
        for part in value.split("|"):
            part = part.strip()
            if part and part not in seen:
                seen.append(part)
    return seen


def main():
    manifest = json.loads((SOURCE_DIR / "manifest.json").read_text())
    source = SOURCE_DIR / "TF-RISet.tsv"
    require(hashlib.sha256(source.read_bytes()).hexdigest() == manifest["tf_riset_sha256"],
            "Source snapshot checksum differs from the reviewed import.")
    lines = source.read_text().splitlines()
    header = [h.split(")", 1)[1] for h in lines[0].split("\t")]
    interactions = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]
    require(len(interactions) == manifest["interaction_count"], "Unexpected interaction count in snapshot.")

    with DATASET.open(newline="") as f:
        dataset = list(csv.DictReader(f))
    table = {r["tf_name"].casefold(): r["tf_name"] for r in dataset if r["method"] in ("", "SELEX")}
    other = {r["tf_name"].casefold(): r["tf_name"] for r in dataset}
    require(all(name.casefold() in table for name in ALIASES.values()),
            "An existing alias target has changed; review the mapping.")

    def local_name(regulondb_name):
        if regulondb_name in ALIASES:
            return ALIASES[regulondb_name]
        key = regulondb_name.casefold()
        return table.get(key) or other.get(key) or regulondb_name

    sites = {}
    for ri in interactions:
        if not ri["tfrsSeq"]:
            continue                       # interaction without a located binding site
        require(re.fullmatch(r"[acgt]*[ACGT]+[acgt]*", ri["tfrsSeq"]), f"Unexpected sequence format in {ri['id']}.")
        sites.setdefault((ri["regulatorName"], ri["tfrsID"]), []).append(ri)

    rows = []
    for (regulator, site_id), group in sites.items():
        require(len({(ri["tfrsLeft"], ri["tfrsRight"]) for ri in group}) == 1, f"Site {site_id} has conflicting coordinates.")
        # RegulonDB prints the sequence on the strand of each interaction; store the
        # forward-strand reading (site in capitals, 10 bp flanks in lower case).
        forward = [ri["tfrsSeq"] for ri in group if ri["strand"] == "forward"]
        flanked = forward[0] if forward else group[0]["tfrsSeq"][::-1].translate(COMPLEMENT)
        core = re.sub(r"[a-z]", "", flanked)
        left, right = int(group[0]["tfrsLeft"]), int(group[0]["tfrsRight"])
        require(right - left + 1 == len(core), f"Site {site_id} length differs from its coordinates.")
        levels = {ri["confidenceLevel"] for ri in group}
        evidence = merged(ri["tfrsEvidence"] for ri in group)
        rows.append({
            "tf_name": local_name(regulator),
            "regulondb_tf_name": regulator,
            "regulondb_tf_id": group[0]["regulatorId"],
            "site_id": site_id,
            "left": left,
            "right": right,
            "length": len(core),
            "sequence": core,
            "sequence_with_flanks": flanked,
            "function": "; ".join(sorted(merged(ri["riFunction"] for ri in group))),
            "target_genes": "; ".join(merged(ri["firstGene"] for ri in group)),
            "promoters": "; ".join(merged(ri["promoterName"] for ri in group)),
            "confidence": next(CONFIDENCE[c] for c in "CSW?" if c in levels),
            "experimental": "yes" if any(code.startswith("EXP-") for code in evidence) else "no",
            "site_evidence": "; ".join(evidence),
            "site_pmids": "; ".join(merged(ri["tfrsPMIDS"] for ri in group)),
            "interaction_ids": "; ".join(ri["id"] for ri in group),
        })
    rows.sort(key=lambda r: (r["tf_name"].casefold(), r["left"], r["right"], r["site_id"]))

    with OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    per_tf = Counter(r["tf_name"] for r in rows)
    in_table = [name for name in per_tf if name.casefold() in table]
    outside = sorted(name for name in per_tf if name.casefold() not in table)
    print(f"Wrote {len(rows)} binding sites for {len(per_tf)} TFs to {OUTPUT.relative_to(ROOT)} "
          f"(RegulonDB {manifest['regulondb_release']}).")
    print(f"{len(in_table)} of the {len(table)} TF-table entries have at least one site; "
          f"{sum(1 for n in in_table if per_tf[n] >= 10)} have ten or more.")
    if outside:
        print("Not TF-table entries: " + ", ".join(outside))


if __name__ == "__main__":
    main()
