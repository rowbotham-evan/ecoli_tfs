# ChIP-seq compendium source provenance (Lally et al. 2025)

Lally P, Gómez-Romero L, Tierrafría VH, Aquino P, Rioualen C, Zhang X, Kim S, Baniulyte G, Plitnick J, Smith C, Babu M, Collado-Vides J, Wade JT, Galagan JE. *Nature Communications* 16, 4255 (2025).
[Paper](https://doi.org/10.1038/s41467-025-58862-8), [GEO GSE268698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268698).

`supplementary_data_2.xlsx` is the paper's Supplementary Data 2 ("TF ChIP-Seq Stats"), retrieved 2026-11-06 from the publisher (URL and SHA-256 checksums in `manifest.json`). It lists the 139 *E. coli* K-12 MG1655 TFs whose in vivo ChIP-seq experiments passed the paper's QC (at least two passing experiments per TF), with each TF's b-number, expression modes used (native SPA / 3xFLAG / V5 tags, or arabinose-inducible pBAD24-3xFLAG), and union binding-region counts. `tf_chip_stats.csv` is the same table converted to CSV (✔/✘ → yes/no), and is the input the importer validates by checksum.

All experiments used one standardized condition: M9 minimal medium + 0.4% glycerol, log phase, 30 °C (inducible strains: + 0.2% arabinose for 1 h before crosslinking). The reference genome is U00096.3.

To import into the website:

```sh
python3 scripts/import_chip_compendium.py
```

The importer adds **one `CHIP` record per TF** — a statement that the TF has in vivo ChIP-seq binding data in this compendium — not the binding regions themselves. Per-TF union-region counts remain in the snapshot; region coordinates are available per experiment in GEO GSE268698 (`*_filtered_regions.gff`) and as per-TF union peaks in RegulonDB, and could be imported later. Repeated runs are idempotent: existing records are never duplicated or overwritten, and new record IDs continue above the current maximum.

`tf_name` uses an existing table name for verified one-to-one synonyms (e.g. dgsA → Mlc, ycgE → BluR, yihW → CsqR; full map in the importer, following the same policy as `scripts/import_pptp.py`); otherwise the paper's gene name is kept with the first letter capitalized. `sxy` is deliberately **not** mapped to the table's composite `CRP-Sxy` entry. 36 of the 139 TFs are not TF-table entries and appear in the coverage page's excluded-names footnote. `regulondb` is left blank, consistent with the other paper-derived imports; a blank does not assert absence from RegulonDB.

The article is published under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/). This import stores the factual TF list and per-TF statistics from Supplementary Data 2, and the logo images from Supplementary Data 3 reproduced unmodified (see below), with attribution; the paper's text and other figures are not reproduced.

## BoltzNet models (Supplementary Data 3)

The paper's Supplementary Data 3 ("Model Weight Matrices and Corresponding PFMs", `41467_2025_58862_MOESM5_ESM.xlsx`, retrieved 2026-11-12; URL and SHA-256 in `manifest.json`) lists all 139 TFs. For each TF with a model it shows two logo images: the 25-bp BoltzNet energy matrix, and a frequency matrix (PFM). The PFM is built by counting bases in the TF's own peak sequences that score at least 0.8 of the matrix's maximum. The 15 TFs without a model are marked `None`. The workbook contains no numeric matrix values (248 embedded PNGs, no number cells). The paper's data-availability statement says the matrices are deposited in RegulonDB, but RegulonDB's GraphQL API does not expose them.

`boltznet_models.csv` records, for each of the 139 TFs, its b-number and whether it has a model (124 do). The 15 without a model (allS, csgD, csiR, gadX, gcvR, hipA, hipB, mazE, nanR, sgrR, stpA, yafC, yedW, yfiE, yjgJ) differ from the failure list in the paper's Methods text: the text names aaeR and rutR, which have logos here, and omits gcvR, hipA, sgrR and yedW. This supplement is used as the authoritative list.

`logos/` holds the 248 logo images from that workbook, extracted byte-for-byte unmodified: `<gene>_weight_matrix.png` (the energy matrix; letters above zero favor binding) and `<gene>_pfm.png` (the derived frequency matrix, in bits) for each modelled TF. Each image is matched to its TF by the row it is anchored in, and its original name inside the workbook (e.g. `aaeR_WM.png`) is checked against that row. `logos/index.csv` lists gene, kind, file, original name and SHA-256. The images are shared under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/), which allows unmodified copies with credit for non-commercial use; the motif pages show them with that credit.

`boltznet_models.csv` and `logos/` are regenerated from the downloaded workbook by `python3 scripts/snapshot_motif_supplements.py DAP_ZIP LALLY_XLSX`, which first checks it against `models_xlsx_sha256`. `scripts/build_motifs.py` checks them against `boltznet_models_sha256`, `logos_index_sha256` and each image's SHA-256, and builds `data/motifs.csv` from them.

## BoltzNet weight matrices (authors' site)

The supplement has the models only as images, but the authors publish the numeric models on their BoltzNet site: `https://boltznet.bu.edu/static/data/models/<tf>_wm_bias.json` for each TF in `https://boltznet.bu.edu/static/data/tf_list.txt`. These are the files their [`boltznet`](https://pypi.org/project/boltznet/) Python package downloads. Retrieved 2026-11-12; the site dates the files 29–30 April 2025, except pdhR (6 August 2025). The TF list matches the 124 modelled TFs in Supplementary Data 3 exactly. The model files carry no licence of their own; the package is CC BY-NC 4.0 (James Galagan, 2025), and the site says the code is free for academic use.

`boltznet_weight_matrices.json` stores, for each TF, the 25-position weight matrix exactly as served (one A, C, G, T row per position), its bias term, the SHA-256 of the downloaded file and its Last-Modified date. It is produced by `python3 scripts/snapshot_boltznet_models.py DIR` from a folder of downloads (the script's docstring has the download commands).

`boltznet_logo_check.csv` records, for each TF, whether the site's matrix matches the energy logo published in Supplementary Data 3. Checked 2026-11-12 by drawing every site matrix as an energy logo next to its published image. 119 are identical (`same`). 5 differ (`differs`): nac, pdhR, glnG, allR and ulaR, the five TFs the paper studies most closely and calibrates with binding measurements. Their published logos probably come from separately trained models; which version the authors consider definitive is not stated. The motif pages show both for those five.

`scripts/build_motifs.py` checks both files against `boltznet_weight_matrices_sha256` and `boltznet_logo_check_sha256` and writes the matrices into `data/motif_matrices.json`.
