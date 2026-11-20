# E. Coli

A database of *Escherichia coli* transcription factors, published with GitHub Pages.


The database can be accessed [here](https://rowbotham-evan.github.io/ecoli_tfs/html/).

All HTML pages live in `html/`. For local development, run `python3 -m http.server 8000` from the repository root and open `http://localhost:8000/html/`. Shared assets, data, and `notes.md` remain at the repository root.

The Coverage tab (`html/coverage.html`) charts TF coverage per dataset — overall and restricted to RegulonDB-marked TFs — computed in the browser from `data/dataset.csv`, so it stays current as records are added. Its Motifs section reads `data/motifs.csv` and lists every TF with a motif, ranked by number of motifs. Each TF links to `html/motif.html?tf=<name>`, a page showing where that TF's motifs come from (paper, method, supplementary file, width, sites, E-value) and any binding sites stored for it in `data/dataset.csv`.

The Coverage tab's "Motif correlation between methods" card compares each TF's DAP-seq motif (a frequency matrix) with its ChIP-seq motif (a BoltzNet energy matrix), for the TFs that have both (plus PhoB, whose second motif is the frequency matrix of its aligned sites in `data/dataset.csv`). Both are compared as position probability matrices (the energy matrix is converted with the Boltzmann relation, p ∝ exp(weight)): the pair is aligned with Tomtom, allowing shifts and reverse complements, and scored by the Pearson correlation between the aligned matrices. TFs are listed in two groups, correlation above and below 0.5; the notes give the scores of mismatched pairs for reference. The scores are precomputed: run `conda run -n meme python3 scripts/build_motif_agreement.py` to rebuild `data/motif_agreement.json` (it needs Tomtom from the MEME Suite: `conda create -n meme -c conda-forge -c bioconda meme`), which `assets/motif-agreement.js` loads. The plot is one scatter of every TF's aligned cells; clicking a TF shows it on its own. Each two-motif TF's motif page shows its own correlation and scatter. Below it, the "Does TF identity dominate experimental method?" card uses the same all-pairs scores: a heatmap of every TF's DAP-seq motif against every TF's ChIP-seq motif, the top-1 retrieval accuracy in both directions (is a TF's own partner its most similar motif?), and the same test against TFs of the same DNA-binding family only (families from the UniProt domain notes). These are in the `identity` block of `data/motif_agreement.json`.

Picking a TF on the Database page (or opening `html/index.html?tf=<name>`) shows a dashboard above its records: one tile per dataset (gSELEX, PPTP-seq, DAP-seq, ChIP-seq) saying whether it covers the TF and with how many records, and one card per published motif with its logo, grouped by method. A covered tile filters the table to that method; several tiles can be selected at once to show their records together (`method=PPTP,CHIP`). The selected TF and the method, medium and strain filters are kept in the address (`?tf=&method=&medium=&strain=`), so the browser's Back button, a reload or a shared link returns to the same view; the motif page's back link returns to the TF's dashboard when reached from it. The logos are drawn by `assets/logo.js`, shared with `html/motif.html`.

The table's `#` column numbers the filtered results from 1, continues across pages, and restarts when filters change. `Record ID` identifies a specific data row and stays fixed across filtering and reloads; multiple records for the same TF have different IDs. The CSV stores these IDs in `record_id`. Preserve existing IDs when editing or reordering records, assign new records IDs above the highest previously assigned value (currently 13475), and never reuse retired IDs. The view's row numbers are not exported.

PPTP import: run `python3 scripts/import_pptp.py` to verify or import the pinned significant-response snapshot. Repeated runs preserve IDs and do not duplicate records. See [source provenance](data/sources/pptp/README.md).

ChIP-seq compendium import: run `python3 scripts/import_chip_compendium.py` to import the 139-TF Lally et al. 2025 ChIP-seq compendium (one record per TF). See [source provenance](data/sources/chip_compendium/README.md).

DAP-seq import: run `python3 scripts/import_dap.py` to import the 107 successful *E. coli* TFs from the Baumgart et al. 2021 biotin DAP-seq study (one record per TF). See [source provenance](data/sources/dap/README.md).

RegulonDB binding sites: run `python3 scripts/import_regulondb_sites.py` to rebuild `data/regulondb_sites.csv`, the 3,620 curated TF binding sites (genome position and sequence) for 206 TFs from RegulonDB's TF-RISet, release 14.5.1. It is a separate file; `data/dataset.csv` is not changed. The Coverage tab's "DNA Binding Domain and Transcription Factor Binding Site" card reads it to count the TFs that have both a `tf_binding_domain` sequence and at least one binding-site sequence. RegulonDB's terms restrict integrating its data into other databases without consent, so read the [source provenance](data/sources/regulondb/README.md) before publishing this file.

Motifs: run `python3 scripts/build_motifs.py` to rebuild `data/motifs.csv`, which lists the TFs each source paper computed a binding motif for: the Baumgart et al. 2021 DAP-seq MEME motifs (124 TFs), the Lally et al. 2025 BoltzNet models (124 TFs) and the Fitzgerald et al. 2023 PhoB motif. One row per TF and paper. It records each motif's type, where it is published, (where given) its width, site count and E-value, and for Lally et al. the paths of their published logo images. It also writes `data/motif_matrices.json`, the numeric matrices the motif pages draw logos from: the DAP-seq frequency matrices and the Lally et al. BoltzNet weight matrices from the authors' site (PhoB's second logo is drawn from its 59 stored sites). To refresh the BoltzNet matrices, run `python3 scripts/snapshot_boltznet_models.py DIR` on a folder of downloads. To regenerate the snapshots it reads from the downloaded supplements, run `python3 scripts/snapshot_motif_supplements.py DAP_ZIP LALLY_XLSX`. The build script verifies the snapshot checksums and maps names with the importers' alias lists, so motif rows use the same TF names as the DAP and CHIP records. Provenance is in the [DAP-seq](data/sources/dap/README.md) and [ChIP-seq compendium](data/sources/chip_compendium/README.md) READMEs.


1. Find the proteome for Escherichia coli (strain K12)
2. Then under components click "View proteins"
3. Searched for desired transcription factor
4. Under Entry tab for that protein there is a section called "DNA binding" 
5. Click tools and then copy sequence. 


https://www.uniprot.org/proteomes/UP000000625


