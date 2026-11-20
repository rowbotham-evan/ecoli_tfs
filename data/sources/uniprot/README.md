# UniProt DNA-binding-domain source provenance

UniProtKB, release **2026_04** (released 4 November 2026), retrieved 2026-11-09
from the REST `stream` endpoint. Query, fields, URL, release headers and SHA-256
checksums are in `manifest.json`; both files below are verified by checksum
before import.

The query is the reviewed (Swiss-Prot) subset of the *E. coli* K-12 MG1655
reference proteome **UP000000625** — the same strain and assembly (U00096.3) the
rest of the dataset uses. Restricting to reviewed entries is what makes the
feature coordinates stable: TrEMBL carries predicted features that change
between releases.

`uniprot_k12_reviewed.tsv` is the raw download (4,403 entries, with the
`DNA_BIND` and `DOMAIN` feature columns and the chain sequence).
`dbd_uniprot.tsv` is the derived table of 296 DBD spans — one row per TF-table
entry, or one per subunit for composite entries. Both are snapshotted here so the
import is reproducible offline.

To import into the website:

```sh
python3 scripts/import_uniprot_dbd.py
```

The importer fills `tf_binding_domain` and sets `sequence_strain` to
`K-12 MG1655` on the first base row of each TF entry. Repeated runs are
idempotent. **214 of the 288 TF entries** received a domain: 124 cells that were
empty, 80 that held a shorter fragment, and 10 that already matched.

## Which span is taken as the DBD

A `DNA_BIND` feature is not always a domain. CRP's is annotated
`DNA_BIND 180..186; /note="H-T-H motif"` — the recognition helix alone — while
its usable unit is `DOMAIN 138..210; /note="HTH crp-type"`. The span is therefore
chosen by priority, recorded per row in `dbd_source`:

| Priority | Rule | Rows | Imported |
|---|---|---|---|
| 1 | `DNA_BIND` ≥ 40 aa | 17 | yes |
| 2 | longest `DOMAIN` whose note matches HTH / helix-turn-helix | 204 | yes |
| 3 | short `DNA_BIND` motif | 22 | **no** |
| 4 | full chain (no DNA_BIND, no HTH domain) | 53 | **no** |

Coordinates are 1-based inclusive as UniProt gives them. Imported domains are
46–135 aa (median 61).

Priorities 3 and 4 are **not** imported because they are not comparable units.
The full-chain group mixes genuine non-HTH DNA binders whose fold UniProt
annotates under other feature types (HU, IHF, MetJ, ArgR, Fur, Zur, NikR, Dps,
and the ribbon-helix-helix antitoxins) with large moonlighting enzymes where the
chain is emphatically not a DBD — PutA (1320 aa), AlaS, PepA, AidB, PaaX.
Importing either would put a misleading sequence in the table.

Composite entries (FlhDC, IHF, HU, RcsAB, RcsB-BglJ, CRP-Sxy, and the
toxin–antitoxin operons DinJ-YafQ, HigBA, HipAB, MazEF, RelBE, YefMB) are
resolved to their subunits in `dbd_uniprot.tsv` but are not written to the table,
because one cell cannot hold two subunit sequences.

## Relationship to the sequences that were already there

Of the 95 DBD sequences the table held before this import:

- **10 were byte-identical** to the UniProt span.
- **80 were strict substrings** of it — the table carried the bare ~20 aa
  recognition motif, UniProt gives the full domain around it.
- **0 conflicted.** The importer refuses to replace a sequence that is not
  contained in the new span, so a future release that genuinely disagrees will
  stop the import rather than overwrite curation.

Five entries keep their original fragment because UniProt offers only a short
`DNA_BIND` motif for them, which would be no improvement: **CitB** (20 aa),
**DcuR** (20 aa), **HNS** (6 aa), **Mlc** (10 aa), **SfsA** (20 aa). These retain
`sequence_strain = K-12 W3110`, their original provenance.

## Name resolution

Table entries are matched to UniProt by gene name, trying in order: the
established one-to-one aliases already recorded in `scripts/import_*.py`
(Mlc→dgsA, NtrC→glnG, LctR→lldR, CsqR→yihW, …), parenthetical alternates such as
`YciT(DeoT)`, then the table name itself. Every gene synonym UniProt lists is
indexed, so b-numbers resolve too.

Three table names were found to be wrong, by searching the table's own DBD
fragment as an exact peptide against all 4,403 chains — the fragment identifies
the protein unambiguously:

| Table name | Actually | Note |
|---|---|---|
| **FlhA** | `fhlA` / FhlA (P19323), formate hydrogenlyase activator | UniProt `flhA` is the 692 aa flagellar biosynthesis membrane protein with no DNA_BIND feature. Transposed letters. |
| **MrlA** | `mlrA` / MlrA (P33358), MerR-like regulator | No `mrlA` gene exists in UniProt. Transposed letters. |
| **DeoT** | `yciT` (P76034) | The same protein as the table's separate `YciT(DeoT)` entry, confirmed by sequence; RegulonDB also treats these as one TF. |

These are mapped in the importer so the right domain is fetched, but the **table
names themselves are left unchanged** — renaming entries or merging DeoT with
YciT(DeoT) moves the TF counts and every coverage figure, so it belongs in a
separate, deliberate edit.

Four entries have no UniProt protein record at all, checked reviewed and
unreviewed under taxonomy 83333 and 562, by gene and protein name: **CitR**,
**DicD**, **GadRcs**, **NfeR**. They are RegulonDB-only or obsolete names and get
no domain.

## Licence

UniProtKB is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Only factual sequence and feature-coordinate data are stored here, with the
release recorded for citation.
