# RegulonDB TF-RISet source provenance

RegulonDB release **14.5.1**, file `TF-RISet` (`RDBECOLIDLF00002`, file version 1.0, created 11-11-2026), retrieved 2026-11-18 from the public read-only [GraphQL API](https://regulondb.ccg.unam.mx/graphql). The query, release fields, column descriptions and SHA-256 checksum are in `manifest.json`.

Cite: Salgado H, Gama-Castro S, et al. RegulonDB v12.0: a comprehensive resource of transcriptional regulation in *E. coli* K-12. *Nucleic Acids Research*, 2023, gkad1072. https://doi.org/10.1093/nar/gkad1072

## License: read before publishing

RegulonDB is free for academic and noncommercial use. Its terms also say that a user may not "integrate RegulonDB partly or as a whole into other databank systems, without prior written consent from CCG-UNAM" ([terms and conditions](https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions)). `TF-RISet.tsv` and `data/regulondb_sites.csv` are copies of RegulonDB content, so get that consent before pushing them to the public site, or keep both files out of the published repository.

## Files

`TF-RISet.tsv` is the download's `content` field, unchanged: 5,785 regulatory interactions between a transcription factor and the promoter, gene or transcription unit it regulates, for 242 TFs. 4,733 interactions name a binding site with a sequence; 21 name a site with coordinates but no sequence, and the rest have no located site.

To rebuild the site table:

```sh
python3 scripts/import_regulondb_sites.py
```

The script checks the snapshot's checksum and writes `data/regulondb_sites.csv`: **3,620 binding sites for 206 TFs**, one row per TF and site. `data/dataset.csv` is not changed. The output is rewritten in full on each run, so repeated runs give the same file.

| Column | Meaning |
| --- | --- |
| `tf_name` | The TF's name in `data/dataset.csv` |
| `regulondb_tf_name`, `regulondb_tf_id` | RegulonDB's name and identifier for the TF |
| `site_id` | RegulonDB binding-site identifier |
| `left`, `right`, `length` | Site position in U00096.3 (1-based, inclusive) and its length in bp |
| `sequence` | The site's sequence on the forward strand |
| `sequence_with_flanks` | The same with RegulonDB's 10 bp flanks in lower case |
| `function` | activator, repressor or both, over the site's interactions |
| `target_genes`, `promoters` | First transcribed genes and promoters the site regulates |
| `confidence` | Highest RegulonDB confidence level among the site's interactions (confirmed, strong, weak, unknown) |
| `experimental` | yes if any site evidence code is experimental (`EXP-…`) |
| `site_evidence`, `site_pmids` | RegulonDB evidence codes and PubMed IDs for the site |
| `interaction_ids` | The TF-RISet interactions the row was merged from |

RegulonDB prints each site on the strand of its interaction, so a site shared by two divergent promoters appears in both orientations. The script stores the forward-strand reading. All 3,620 sequences were checked against U00096.3 at their coordinates when the import was written.

## Names

`tf_name` uses the TF-table entry for the 23 RegulonDB names that differ (for example Rob → CbpB, DpiA → CitB, H-NS → HNS; the full list is `ALIASES` in the script). 205 of the 288 TF-table entries have at least one site and 60 have ten or more. RacR has 4 sites but is not a TF-table entry. The table lists RegulonDB's YciT twice, as DeoT and YciT(DeoT); its 4 sites are stored under YciT(DeoT) only, so DeoT counts as having none.
