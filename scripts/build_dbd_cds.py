#!/usr/bin/env python3
"""Map each TF's DNA-binding domain onto its native nucleotide coding sequence.

Evo 2 reads nucleotides, not amino acids, so prompting it with a TF's DBD means
supplying the DNA that encodes the DBD. Every TF here is E. coli K-12, so the
real codons exist in the genome and no reverse-translation is needed.

Input:
  data/sources/uniprot/dbd_uniprot.tsv   (accession, dbd_start, dbd_end, dbd_sequence)
  <cds_fasta>  NCBI GCF_000005845.2_ASM584v2_cds_from_genomic.fna

The CDS FASTA headers carry [db_xref=UniProtKB/Swiss-Prot:<acc>], so the join to
UniProt is exact. Sequences in that file are already in coding orientation, so
minus-strand genes need no special handling.

Every slice is verified by translating it back and comparing to the UniProt DBD
sequence; rows that fail are written out with the reason rather than dropped.

Usage:
  python3 -I scripts/build_dbd_cds.py <cds_fasta> [-o data/dbd_cds.tsv]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DBD_TSV = REPO / "data" / "sources" / "uniprot" / "dbd_uniprot.tsv"
SUPPLEMENT_TSV = REPO / "data" / "sources" / "uniprot" / "dbd_cds_supplement.tsv"
DEFAULT_OUT = REPO / "data" / "dbd_cds.tsv"

CODON_TABLE = {}
_BASES = "TCAG"
_AAS = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
for _i, _aa in enumerate(_AAS):
    CODON_TABLE[_BASES[_i >> 4] + _BASES[(_i >> 2) & 3] + _BASES[_i & 3]] = _aa

# Bacterial initiation codons translate as Met at position 1 (NCBI table 11).
ALT_STARTS = {"GTG", "TTG", "ATT", "CTG"}

ACC_RE = re.compile(r"db_xref=UniProtKB/(?:Swiss-Prot|TrEMBL):([A-Z0-9]+)")
GENE_RE = re.compile(r"\[gene=([^\]]+)\]")
LOCUS_RE = re.compile(r"\[locus_tag=([^\]]+)\]")
PROTID_RE = re.compile(r"\[protein_id=([^\]]+)\]")


def translate(nt: str, first_codon_is_start: bool = False) -> str:
    aas = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        codon = nt[i : i + 3]
        if i == 0 and first_codon_is_start and codon in ALT_STARTS:
            aas.append("M")
            continue
        aas.append(CODON_TABLE.get(codon, "X"))
    return "".join(aas)


def read_cds_fasta(path: Path) -> list[dict]:
    records, header, chunks = [], None, []

    def flush():
        if header is None:
            return
        seq = "".join(chunks).upper().replace("U", "T")
        records.append(
            {
                "accession": (ACC_RE.search(header).group(1) if ACC_RE.search(header) else ""),
                "gene": (GENE_RE.search(header).group(1) if GENE_RE.search(header) else ""),
                "locus_tag": (LOCUS_RE.search(header).group(1) if LOCUS_RE.search(header) else ""),
                "protein_id": (PROTID_RE.search(header).group(1) if PROTID_RE.search(header) else ""),
                "cds_nt": seq,
            }
        )

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header, chunks = line, []
            elif header is not None:
                chunks.append(line.strip())
    flush()
    return records


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cds_fasta", type=Path, help="GCF_000005845.2_ASM584v2_cds_from_genomic.fna")
    ap.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    if not args.cds_fasta.is_file():
        print(f"error: no such file: {args.cds_fasta}", file=sys.stderr)
        return 1

    cds = read_cds_fasta(args.cds_fasta)
    by_acc, by_gene = {}, {}
    for rec in cds:
        if rec["accession"]:
            by_acc.setdefault(rec["accession"], rec)
        if rec["gene"]:
            by_gene.setdefault(rec["gene"].lower(), rec)
    print(f"CDS records: {len(cds)}  with UniProt xref: {len(by_acc)}")

    # A few K-12 genes have no CDS feature in U00096.3 (e.g. glpR) and so never
    # appear in the RefSeq extraction. Those are supplied from ENA, once, here.
    n_supp = 0
    if SUPPLEMENT_TSV.is_file():
        with SUPPLEMENT_TSV.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                seq = row.get("cds_nt", "").strip().upper()
                if not seq:
                    continue
                rec = {
                    "accession": row.get("accession", ""),
                    "gene": row.get("gene", ""),
                    "locus_tag": row.get("locus_tag", ""),
                    "protein_id": row.get("protein_id", ""),
                    "cds_nt": seq,
                }
                if rec["accession"]:
                    by_acc.setdefault(rec["accession"], rec)
                if rec["gene"]:
                    by_gene.setdefault(rec["gene"].lower(), rec)
                n_supp += 1
        print(f"supplemented from {SUPPLEMENT_TSV.name}: {n_supp}")

    with DBD_TSV.open() as fh:
        tf_rows = list(csv.DictReader(fh, delimiter="\t"))

    out_rows, stats = [], {"ok": 0, "no_cds": 0, "bad_range": 0, "mismatch": 0, "no_dbd": 0}

    for row in tf_rows:
        start, end = row.get("dbd_start", ""), row.get("dbd_end", "")
        rec = by_acc.get(row.get("accession", "")) or by_gene.get(row.get("uniprot_gene", "").lower())
        joined_on = (
            "accession" if row.get("accession") in by_acc
            else ("gene" if rec else "")
        )

        out = {
            "tf_name": row.get("tf_name", ""),
            "accession": row.get("accession", ""),
            "uniprot_gene": row.get("uniprot_gene", ""),
            "locus_tag": rec["locus_tag"] if rec else "",
            "protein_id": rec["protein_id"] if rec else "",
            "joined_on": joined_on,
            "dbd_start": start,
            "dbd_end": end,
            "dbd_sequence": row.get("dbd_sequence", ""),
            "dbd_cds_nt": "",
            "cds_nt": "",
            "cds_nt_length": len(rec["cds_nt"]) if rec else "",
            "status": "",
        }

        if not (start.isdigit() and end.isdigit()):
            out["status"] = "no_dbd_annotation"
            stats["no_dbd"] += 1
            out_rows.append(out)
            continue
        if rec is None:
            out["status"] = "no_cds_match"
            stats["no_cds"] += 1
            out_rows.append(out)
            continue

        # Drop the terminal stop codon so protein positions map 1:1 onto codons.
        coding = rec["cds_nt"]
        if len(coding) >= 3 and CODON_TABLE.get(coding[-3:]) == "*":
            coding = coding[:-3]
        # Stop-stripped CDS: what a whole-gene prompt variant would use.
        out["cds_nt"] = coding

        s, e = int(start), int(end)
        if s < 1 or e * 3 > len(coding) or s > e:
            out["status"] = f"range_outside_cds (cds {len(coding)//3} aa, dbd {s}-{e})"
            stats["bad_range"] += 1
            out_rows.append(out)
            continue

        dbd_nt = coding[(s - 1) * 3 : e * 3]
        back = translate(dbd_nt, first_codon_is_start=(s == 1))
        expected = row.get("dbd_sequence", "")

        out["dbd_cds_nt"] = dbd_nt
        if expected and back != expected:
            n_diff = sum(1 for a, b in zip(back, expected) if a != b) + abs(len(back) - len(expected))
            out["status"] = f"translation_mismatch ({n_diff} aa)"
            stats["mismatch"] += 1
        else:
            out["status"] = "ok"
            stats["ok"] += 1
        out_rows.append(out)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    total = len(out_rows)
    print(f"\nTFs in dbd_uniprot.tsv: {total}")
    print(f"  verified (translates back to the UniProt DBD): {stats['ok']}")
    print(f"  translation mismatch:                          {stats['mismatch']}")
    print(f"  no CDS match:                                  {stats['no_cds']}")
    print(f"  DBD range outside CDS:                         {stats['bad_range']}")
    print(f"  no DBD annotation:                             {stats['no_dbd']}")
    print(f"\nwrote {args.out}")

    bad = [r for r in out_rows if r["status"] not in ("ok", "no_dbd_annotation")]
    if bad:
        print("\nrows needing attention:")
        for r in bad[:25]:
            print(f"  {r['tf_name']:12s} {r['accession']:10s} {r['status']}")
        if len(bad) > 25:
            print(f"  ... and {len(bad) - 25} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
