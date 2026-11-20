# DAP-seq source provenance (Baumgart et al. 2021)

Baumgart LA, Lee JE, Salamov A, Dilworth DJ, Na H, Mingay M, Blow MJ, Zhang Y, Yoshinaga Y, Daum CG, O'Malley RC. *Nature Methods* 18, 1499–1510 (2021).
[Paper](https://doi.org/10.1038/s41592-021-01312-2), raw reads in [NCBI BioProject PRJNA758434](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA758434). The method follows the DAP-seq protocol of [Bartlett et al., *Nat Protoc* 12, 1659–1672 (2017)](https://doi.org/10.1038/nprot.2017.055) (biotin/multiDAP variant). This is the single publication RegulonDB's *E. coli* DAP-seq collection is built from (RegulonDB 11.0, Tierrafría et al. 2022).

`supplementary_tables.xlsx` is the paper's combined Supplementary Tables 1–6, retrieved 2026-11-06 from the publisher (URL and SHA-256 checksums in `manifest.json`). Supplementary Table 2 lists the 214 *E. coli* K-12 MG1655 TFs cloned and assayed by in vitro biotin DAP-seq, with locus tag, protein ID, and a success flag; 107 are marked successful. Supplementary Table 4 is the 92-TF multiDAP benchmark (replicate peak counts), a subset of the 107. `dap_tf_list.csv` is the Table 2 *E. coli* subset converted to CSV (1/0 → yes/no, plus an `in_multidap_benchmark` column from Table 4), and is the input the importer validates by checksum.

DAP-seq is an in vitro assay: in vitro–expressed, biotinylated TFs are bound to fragmented naked genomic DNA, so records carry no growth condition. The assay pool contained 48 bacterial genomes, including *E. coli* K-12 MG1655 (GCF_000005845.2, i.e. U00096.3) and BW25113 (GCF_004355105.1); the TFs themselves are MG1655 proteins.

To import into the website:

```sh
python3 scripts/import_dap.py
```

The importer adds **one `DAP` record per successful TF** (107) — a statement that the TF has DAP-seq binding data in this study — not the binding regions themselves. Peak coordinates for all 48 genomes are in the paper's Supplementary Data 1 (ZIP, ~164 MB, not snapshotted here) and per-TF datasets in RegulonDB, and could be imported later. Repeated runs are idempotent: existing records are never duplicated or overwritten, and new record IDs continue above the current maximum.

`tf_name` uses an existing table name for verified one-to-one synonyms (alsR → RpiR, cecR → YbiH, dpiA → CitB, glnG → NtrC, hypT → YjiE, lldR → LctR, nimR → YeaM(NimR), ttdR → Dan, yiaJ → YiaJ(BdcA); same policy as `scripts/import_pptp.py` and `scripts/import_chip_compendium.py`); otherwise the paper's gene name is kept with the first letter capitalized. 4 of the 107 TFs (CsiR, PdeL, Rob, UidR) are not TF-table entries and appear in the coverage page's excluded-names footnote. `regulondb` is left blank, consistent with the other paper-derived imports; a blank does not assert absence from RegulonDB.

Stored here, with attribution: the factual TF list and per-TF statistics from the supplementary tables, and the 124 *E. coli* motif matrices from Supplementary Data 2 (numbers only; see below). The paper's text and figures are not reproduced.

## Motifs (Supplementary Data 2)

The paper's Supplementary Data 2 (`41592_2021_1312_MOESM5_ESM.zip`, retrieved 2026-11-12; URL and SHA-256 in `manifest.json`) holds one MEME v4.12.0 motif file per TF: 124 for *E. coli* and 118 for *Pseudomonas simiae*, the paper's 242 motifs. Each file has exactly one motif. MEME was run on peak summits pooled across all 48 genomes in the multiDAP assay, so most of the aligned sites are not *E. coli* DNA.

`dap_motifs.csv` lists the 124 *E. coli* motifs: the paper's gene name, the file inside the zip, and the motif's width, site count and E-value as read from its `letter-probability matrix` line. `dap_successful` is the Supplementary Table 2 flag. 20 TFs with a motif are not marked successful there, and three successful TFs (hns, hypT, stpA) have no motif. `dap_motif_matrices.json` holds each of those motifs' letter-probability matrix (one A, C, G, T row per position, as printed by MEME), keyed by the paper's gene name, so the motif pages can draw a logo. The rest of each MEME file (the aligned site instances, run log and so on) is not stored.

Both files are regenerated from the downloaded zip by `python3 scripts/snapshot_motif_supplements.py DAP_ZIP LALLY_XLSX`, which first checks the zip against `motifs_zip_sha256`. `scripts/build_motifs.py` checks them against `dap_motifs_sha256` and `dap_motif_matrices_sha256` and builds `data/motifs.csv` and `data/motif_matrices.json`.
