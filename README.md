# E. coli Transcription Factor Database

What is known about how each *Escherichia coli* transcription factor (TF) binds DNA, gathered from four kinds of experiment, with checks on whether the experiments agree. [Open the site](https://rowbotham-evan.github.io/ecoli_tfs/html/).

## The idea

A TF recognizes a short DNA pattern, its **motif**. Each experiment sees it from a different angle:

| Dataset | What it measures |
| --- | --- |
| gSELEX (TEC database) | Binding to genomic DNA in a test tube |
| DAP-seq (Baumgart et al. 2021) | Binding in a test tube, with a motif per TF |
| ChIP-seq (Lally et al. 2025; Fitzgerald et al. 2023 for PhoB) | Binding inside living cells, with a motif per TF |
| PPTP-seq (Han et al. 2023) | Which promoters respond when the TF is knocked down (an effect, not binding) |

If two independent methods give the same motif for a TF, the motif is probably real. The site is built around that check.

## Pages

- **Database**: every record. Pick a TF, or search by TF or gene, to see its datasets and motifs.
- **Coverage**: what each dataset covers, plus two motif comparisons:
  - *Motif correlation*: how alike a TF's DAP-seq and ChIP-seq motifs are.
  - *TF identity*: whether a TF's motif matches its own partner better than any other TF's. This is the control for the first.
- **Motif page** (one per TF): each motif's logo and source.

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
| Motif list and matrices | `scripts/build_motifs.py` |
| Motif comparison | `conda run -n meme python3 scripts/build_motif_agreement.py` (needs Tomtom: `conda create -n meme -c conda-forge -c bioconda meme`) |

Each source's provenance is in `data/sources/<source>/README.md`.

## How coverage is counted

- A **TF entry** is a distinct `tf_name` with a base-table record (a row whose method is blank or SELEX). DeoT and YciT(DeoT) are listed separately although RegulonDB treats them as one TF, and RacR is absent.
- The **RegulonDB** figures keep entries marked ✓ in the `regulondb` column.
- PPTP-seq rows are expression responses, which may be indirect; SELEX, DAP-seq and ChIP-seq rows are DNA-binding evidence.
- 30 TF names in the datasets are not TF-table entries (complex subunits, toxin–antitoxin genes, uncharacterized y-genes) and are left out of the coverage figures: AaeR, BglJ, CsiR, DicC, DinJ, FlhC, FlhD, GcvR, HigB, HipA, Hns, HupA, HupB, IhfA, IhfB, MazF, PdeL, RacR, RcsA, RelE, Rob, SfsB, SgcR, Sxy, Ttk, UidR, YafN, YafQ, YdaS, YoeB.

## Rules to keep

- **Record IDs are permanent.** Keep existing `record_id` values, number new records above 13475, and never reuse one.
- **RegulonDB has terms.** Read `data/sources/regulondb/README.md` before publishing `data/regulondb_sites.csv`.
