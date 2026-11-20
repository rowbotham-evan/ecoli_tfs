#!/usr/bin/env python3
"""Build data/motifs.csv and data/motif_matrices.json from the motif snapshots.

data/motifs.csv has one row per (TF, paper): which TFs each source paper
computed a binding motif for, the motif's type, where it is published, its
width, site count and E-value, and (for Lally et al.) the published logo images
and the address of the TF's numeric weight matrix on the authors' BoltzNet site.
data/motif_matrices.json holds the numeric matrices, keyed by TF name and
paper, for drawing logos: {"kind": "frequency"} for the DAP-seq MEME motifs and
{"kind": "energy", "bias", "published_logo"} for the BoltzNet weight matrices.

Sources (regenerate the first two with scripts/snapshot_motif_supplements.py):
  - data/sources/dap/: Baumgart et al. 2021, Supplementary Data 2. One MEME
    motif per E. coli TF (124): dap_motifs.csv and dap_motif_matrices.json.
  - data/sources/chip_compendium/: Lally et al. 2025, Supplementary Data 3.
    A BoltzNet energy matrix and its derived PFM for 124 of the 139 ChIP-seq
    TFs, published only as logo images: boltznet_models.csv and logos/. The
    numeric weight matrices come from the authors' BoltzNet site:
    boltznet_weight_matrices.json (scripts/snapshot_boltznet_models.py), and
    boltznet_logo_check.csv records whether each one matches the published logo.
  - Fitzgerald et al. 2023: the PhoB MEME motif (Fig. 2B). Its 59 sites are
    already in dataset.csv; motif.html draws the logo from them.

Every snapshot is checked against the SHA-256 in its manifest. TF names follow
the importers' alias maps, so motif rows use the same table names as the DAP
and CHIP records. Both outputs are rewritten on every run.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.dont_write_bytecode = True  # importing the alias maps should not leave __pycache__ in scripts/
sys.path.insert(0, str(ROOT / "scripts"))
from import_chip_compendium import ALIASES as CHIP_ALIASES  # noqa: E402
from import_dap import ALIASES as DAP_ALIASES  # noqa: E402

DATASET = ROOT / "data" / "dataset.csv"
OUTPUT = ROOT / "data" / "motifs.csv"
MATRICES = ROOT / "data" / "motif_matrices.json"
DAP_DIR = ROOT / "data" / "sources" / "dap"
CHIP_DIR = ROOT / "data" / "sources" / "chip_compendium"
FIELDS = ["tf_name", "method", "paper", "source", "motif_type", "published_as",
          "paper_gene", "width", "sites", "e_value", "supplement_url", "supplement_file",
          "logo_images", "matrix_url"]
PAPER_ORDER = {"Baumgart et al. 2021": 0, "Lally et al. 2025": 1, "Fitzgerald et al. 2023": 2}


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

    chip = manifest(CHIP_DIR)
    models = [r for r in read_csv(checked(CHIP_DIR, "boltznet_models.csv", "boltznet_models_sha256"))
              if r["has_model"] == "yes"]
    require(len(models) == chip["model_count"], "Unexpected BoltzNet model count.")
    logo_dir = CHIP_DIR / "logos"
    require(sha256(logo_dir / "index.csv") == chip["logos_index_sha256"],
            "logos/index.csv differs from the reviewed snapshot.")
    logos = read_csv(logo_dir / "index.csv")
    require(len(logos) == chip["logo_count"], "Unexpected logo count.")
    images = {}
    for l in logos:
        require(sha256(logo_dir / l["file"]) == l["sha256"], f"{l['file']} differs from logos/index.csv.")
        images.setdefault(l["tf"], {})[l["kind"]] = (logo_dir / l["file"]).relative_to(ROOT).as_posix()
    weights = json.loads(checked(CHIP_DIR, "boltznet_weight_matrices.json", "boltznet_weight_matrices_sha256").read_text())
    check = {r["tf"]: r["published_logo"] for r in
             read_csv(checked(CHIP_DIR, "boltznet_logo_check.csv", "boltznet_logo_check_sha256"))}
    require(set(weights) == {r["tf"] for r in models} == set(check),
            "BoltzNet weight matrices, logo check and model list cover different TFs.")
    require(set(check.values()) <= {"same", "differs"}, "Unexpected value in boltznet_logo_check.csv.")
    for r in models:
        pair = images.get(r["tf"], {})
        require(set(pair) == {"weight_matrix", "pfm"}, f"Missing logo images for {r['tf']}.")
        out.append({
            "tf_name": table_name(r["tf"], CHIP_ALIASES), "method": "CHIP", "paper": "Lally et al. 2025",
            "source": "https://doi.org/10.1038/s41467-025-58862-8",
            "motif_type": "BoltzNet energy matrix + derived frequency matrix",
            "published_as": "Logo images (Supplementary Data 3); numeric weight matrix on the authors' BoltzNet site",
            "paper_gene": r["tf"], "width": "25", "sites": "", "e_value": "",
            "supplement_url": chip["models_url"], "supplement_file": "",
            "logo_images": f"{pair['weight_matrix']};{pair['pfm']}",
            "matrix_url": chip["boltznet_models_url_pattern"].replace("<tf>", r["tf"]),
        })
        matrices.setdefault(table_name(r["tf"], CHIP_ALIASES), {})["Lally et al. 2025"] = {
            "kind": "energy", "matrix": weights[r["tf"]]["matrix"], "bias": weights[r["tf"]]["bias"],
            "published_logo": check[r["tf"]]}

    rows = read_csv(DATASET)
    phob_sites = [r for r in rows if r["tf_name"] == "PhoB" and r["dna_binding_sequence"]]
    require(len(phob_sites) == 59, "Expected the 59 PhoB motif sites from Fitzgerald et al. 2023.")
    out.append({
        "tf_name": "PhoB", "method": "CHIP", "paper": "Fitzgerald et al. 2023",
        "source": "https://doi.org/10.1128/mbio.02535-22",
        "motif_type": "MEME frequency matrix",
        "published_as": "Figure 2B; the 59 sites are in Table 1",
        "paper_gene": "phoB", "width": "18", "sites": "59", "e_value": "3.0e-72",
        "supplement_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10294691/#tab1",
        "supplement_file": "", "logo_images": "", "matrix_url": "",
    })

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
