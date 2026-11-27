# E. coli Transcription Factor Database

What is known about how each *Escherichia coli* transcription factor (TF) binds DNA, gathered from four kinds of experiment. [Open the site](https://rowbotham-evan.github.io/ecoli_tfs/html/).

## The idea

A TF recognizes a short DNA pattern, its **motif**. Each experiment sees it from a different angle:

| Dataset | What it measures |
| --- | --- |
| gSELEX (TEC database) | Binding to genomic DNA in a test tube |
| DAP-seq (Baumgart et al. 2021) | Binding in a test tube, with a motif per TF |
| ChIP-seq (Lally et al. 2025; Fitzgerald et al. 2023 for PhoB) | Binding inside living cells |
| PPTP-seq (Han et al. 2023) | Which promoters respond when the TF is knocked down (an effect, not binding) |

Only motifs a paper computed from its experiment and published as a numeric frequency matrix are kept: the 124 DAP-seq motifs. The ChIP-seq papers publish theirs as energy matrices, logos or aligned sites, so those are left out.

## Pages

- **Database**: every record. Pick a TF, or search by TF or gene, to see its datasets and motifs.
- **Coverage**: what each dataset covers, and which TFs have a motif.
- **Motif page** (one per TF): the motif's logo and source.
- **Ingestion**: the source papers behind each method, scraped for their Methods sections.
  - A table of methods against the repositories their data went to (GEO, GenBank, ...). Click a count to list its papers below the table.
  - Pick a method above the paper cards for a summary of its papers: the RegulonDB TFs they cover between them, the TFs with a motif, the growth conditions and strains they name, and the TFs not yet in the database.
  - Select a scraped paper for a card of what it covered.

## Run it

```sh
python3 -m http.server 8000    # then open http://localhost:8000/html/
```

Pages are in `html/`, scripts and styles in `assets/`, data in `data/`.

## Rebuild the data

| To rebuild | Run |
| --- | --- |
| Records from a source | `scripts/import_pptp.py`, `import_dap.py`, `import_chip_compendium.py`, `import_uniprot_dbd.py` |
| RegulonDB binding sites | `scripts/import_regulondb_sites.py` |
| A paper's entry on the Ingestion page | `scripts/scrape_paper.py DOI --method PPTP` (open-access papers in Europe PMC; writes only `data/ingestion/papers.json`, never `data/dataset.csv`, so a paper with no records shows as awaiting approval) |
| Motif list and matrices (DAP-seq only) | `scripts/build_motifs.py` (writes `data/motifs.csv` and `data/motif_matrices.json`) |
| The same matrices as a table and in MEME format | `scripts/export_frequency_matrices.py` (writes `data/motif_frequency_matrices.csv` and `.meme`) |

Each source's provenance is in `data/sources/<source>/README.md`.

## How coverage is counted

- A **TF entry** is a distinct `tf_name` with a base-table record (a row whose method is blank or SELEX). DeoT and YciT(DeoT) are listed separately although RegulonDB treats them as one TF, and RacR is absent.
- The **RegulonDB** figures keep entries marked ✓ in the `regulondb` column.
- PPTP-seq rows are expression responses, which may be indirect; SELEX, DAP-seq and ChIP-seq rows are DNA-binding evidence.
- 30 TF names in the datasets are not TF-table entries (complex subunits, toxin–antitoxin genes, uncharacterized y-genes) and are left out of the coverage figures: AaeR, BglJ, CsiR, DicC, DinJ, FlhC, FlhD, GcvR, HigB, HipA, Hns, HupA, HupB, IhfA, IhfB, MazF, PdeL, RacR, RcsA, RelE, Rob, SfsB, SgcR, Sxy, Ttk, UidR, YafN, YafQ, YdaS, YoeB.

## How the Ingestion summary is counted

- A paper's **TFs** are the ones it characterized, curated by hand, not every TF its text mentions. TFs are matched to the database, RegulonDB and the motif list by name, ignoring case.
- **RegulonDB TFs covered** is the union over the method's papers, out of the TFs marked ✓ in `data/dataset.csv`.
- **Growth conditions and strains** are names matched in each paper's text, so a cloning host counts as a strain and a paper that was not scraped names none. M9 variants with the same carbon source count as one condition.

## Rules to keep

- **Record IDs are permanent.** Keep existing `record_id` values, number new records above 13475, and never reuse one.
- **RegulonDB has terms.** Read `data/sources/regulondb/README.md` before publishing `data/regulondb_sites.csv`.
