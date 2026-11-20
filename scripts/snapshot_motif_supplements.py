#!/usr/bin/env python3
"""Regenerate the motif snapshots from the two downloaded supplement files.

Usage:
  python3 scripts/snapshot_motif_supplements.py DAP_ZIP LALLY_XLSX

  DAP_ZIP     Baumgart et al. 2021 Supplementary Data 2
              (41592_2021_1312_MOESM5_ESM.zip, `motifs_url` in data/sources/dap/manifest.json)
  LALLY_XLSX  Lally et al. 2025 Supplementary Data 3
              (41467_2025_58862_MOESM5_ESM.xlsx, `models_url` in data/sources/chip_compendium/manifest.json)

Both inputs are checked against the SHA-256 in their manifest before anything
is written. Outputs:

  data/sources/dap/dap_motifs.csv               one row per E. coli motif: width, sites, E-value
  data/sources/dap/dap_motif_matrices.json      the MEME letter-probability matrices (A, C, G, T per position)
  data/sources/chip_compendium/boltznet_models.csv   which of the 139 TFs have a model
  data/sources/chip_compendium/logos/*.png      the published logo images, byte-for-byte unmodified
  data/sources/chip_compendium/logos/index.csv  gene, kind, file, original image name, SHA-256

The script prints the SHA-256 of each output. data/sources/*/manifest.json keep
the reviewed values, and scripts/build_motifs.py refuses to run if a snapshot
no longer matches them, so update the manifests deliberately after reviewing.
"""

import csv
import hashlib
import json
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DAP_DIR = ROOT / "data" / "sources" / "dap"
CHIP_DIR = ROOT / "data" / "sources" / "chip_compendium"
LOGO_DIR = CHIP_DIR / "logos"
LOGO_KINDS = {2: ("weight_matrix", "_WM"), 3: ("pfm", "_PFM")}   # sheet column index -> (kind, name suffix)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def snapshot_dap(zip_path):
    manifest = json.loads((DAP_DIR / "manifest.json").read_text())
    require(sha256(zip_path.read_bytes()) == manifest["motifs_zip_sha256"],
            "DAP-seq motif zip differs from the reviewed download.")
    with (DAP_DIR / "dap_tf_list.csv").open(newline="") as f:
        successful = {r["tf"]: r["successful"] for r in csv.DictReader(f)}
    rows, matrices = [], {}
    with zipfile.ZipFile(zip_path) as z:
        for name in sorted(z.namelist()):
            m = re.fullmatch(r"motifs/escherichia_coli_motifs/(.+)\.txt", name)
            if not m:
                continue
            gene, text = m.group(1), z.read(name).decode().replace("\r\n", "\n")   # files use CRLF
            blocks = re.findall(r"letter-probability matrix: alength= 4 w= (\d+) nsites= (\d+) E= (\S+)\s*\n"
                                r"((?:[ \t]*[\d.eE+-]+[ \t]+[\d.eE+-]+[ \t]+[\d.eE+-]+[ \t]+[\d.eE+-]+[ \t]*\n)+)", text)
            require(len(blocks) == 1, f"{name}: expected exactly one motif.")
            width, sites, evalue, body = blocks[0]
            matrix = [[round(float(v), 6) for v in line.split()] for line in body.strip().splitlines()]
            require(len(matrix) == int(width), f"{name}: matrix has {len(matrix)} rows, header says {width}.")
            require(all(abs(sum(r) - 1) < 0.01 for r in matrix), f"{name}: a matrix row does not sum to 1.")
            rows.append({"tf": gene, "motif_file": name, "width": width, "sites": sites,
                         "e_value": evalue, "dap_successful": successful.get(gene, "not assayed")})
            matrices[gene] = {"width": int(width), "sites": int(sites), "e_value": evalue,
                              "alphabet": "ACGT", "matrix": matrix}
    require(len(rows) == manifest["motif_count_ecoli"], f"Expected {manifest['motif_count_ecoli']} E. coli motifs.")
    rows.sort(key=lambda r: r["tf"].casefold())
    write_csv(DAP_DIR / "dap_motifs.csv", list(rows[0]), rows)
    ordered = {g: matrices[g] for g in sorted(matrices, key=str.casefold)}
    (DAP_DIR / "dap_motif_matrices.json").write_text(json.dumps(ordered, separators=(",", ":")) + "\n")
    print(f"DAP-seq: {len(rows)} motifs.")


def cell_text(strings, attrs, inner):
    v = re.search(r"<v>(.*?)</v>", inner or "")
    if not v:
        return ""
    return strings[int(v.group(1))] if 't="s"' in attrs else v.group(1)


def snapshot_lally(xlsx_path):
    manifest = json.loads((CHIP_DIR / "manifest.json").read_text())
    require(sha256(xlsx_path.read_bytes()) == manifest["models_xlsx_sha256"],
            "Lally Supplementary Data 3 differs from the reviewed download.")
    with zipfile.ZipFile(xlsx_path) as z:
        strings = [re.sub(r"<[^>]+>", "", s) for s in
                   re.findall(r"<si>(.*?)</si>", z.read("xl/sharedStrings.xml").decode(), re.S)]
        table = {}   # 0-based row index -> {column letter: text}
        for row in re.findall(r"<row [^>]*>.*?</row>", z.read("xl/worksheets/sheet1.xml").decode(), re.S):
            index = int(re.search(r'<row [^>]*r="(\d+)"', row).group(1)) - 1
            table[index] = {col: cell_text(strings, attrs, inner) for col, attrs, inner in
                            re.findall(r'<c r="([A-Z]+)\d+"([^>]*?)(?:/>|>(.*?)</c>)', row, re.S)}
        require([table[0].get(c) for c in "ABCD"] == ["TF", "Bnumber", "Model PWM", "Model PFM"],
                "Unexpected Supplementary Data 3 header.")
        rels = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="\.\./media/([^"]+)"',
                               z.read("xl/drawings/_rels/drawing1.xml.rels").decode()))
        anchors = re.findall(r"<xdr:twoCellAnchor\b.*?</xdr:twoCellAnchor>",
                             z.read("xl/drawings/drawing1.xml").decode(), re.S)
        models, logos = [], []
        LOGO_DIR.mkdir(exist_ok=True)
        for old in LOGO_DIR.glob("*.png"):
            old.unlink()
        seen = set()
        for a in anchors:
            col = int(re.search(r"<xdr:from><xdr:col>(\d+)</xdr:col>", a).group(1))
            row = int(re.search(r"<xdr:from>.*?<xdr:row>(\d+)</xdr:row>", a, re.S).group(1))
            original = re.search(r'<xdr:cNvPr id="\d+" name="([^"]+)"', a).group(1)
            media = rels[re.search(r'r:embed="(rId\d+)"', a).group(1)]
            gene = table[row]["A"]
            kind, suffix = LOGO_KINDS[col]
            require(original.casefold() == f"{gene}{suffix}.png".casefold(),
                    f"Image {original} is anchored in the row for {gene}, column {col}.")
            require((gene, kind) not in seen, f"Two {kind} images for {gene}.")
            seen.add((gene, kind))
            data = z.read(f"xl/media/{media}")
            out = LOGO_DIR / f"{gene}_{kind}.png"
            out.write_bytes(data)
            logos.append({"tf": gene, "kind": kind, "file": out.name, "original_name": original,
                          "sha256": sha256(data)})
        for index in sorted(i for i in table if i > 0):
            r = table[index]
            has_model = r.get("C") != "None"
            require(has_model == (r.get("D") != "None"), f"{r['A']}: one logo marked None, the other not.")
            require(has_model == ((r["A"], "weight_matrix") in seen and (r["A"], "pfm") in seen),
                    f"{r['A']}: logo images do not match the None markers.")
            models.append({"tf": r["A"], "bnumber": r["B"], "has_model": "yes" if has_model else "no"})
    require(len(models) == manifest["tf_count"], "Unexpected TF count in Supplementary Data 3.")
    require(sum(m["has_model"] == "yes" for m in models) == manifest["model_count"], "Unexpected model count.")
    write_csv(CHIP_DIR / "boltznet_models.csv", ["tf", "bnumber", "has_model"], models)
    logos.sort(key=lambda r: (r["tf"].casefold(), r["kind"]))
    write_csv(LOGO_DIR / "index.csv", list(logos[0]), logos)
    print(f"Lally: {len(models)} TFs, {sum(m['has_model'] == 'yes' for m in models)} with a model, {len(logos)} logo images.")


def main():
    require(len(sys.argv) == 3, __doc__)
    snapshot_dap(Path(sys.argv[1]))
    snapshot_lally(Path(sys.argv[2]))
    for path in [DAP_DIR / "dap_motifs.csv", DAP_DIR / "dap_motif_matrices.json",
                 CHIP_DIR / "boltznet_models.csv", LOGO_DIR / "index.csv"]:
        print(f"{sha256(path.read_bytes())}  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
