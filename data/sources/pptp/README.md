# PPTP source provenance

Han Y, Li W, Filko A, Li J, Zhang F. *Nature Communications* 14, 5757 (2023).
[Paper](https://doi.org/10.1038/s41467-023-41572-4), [GEO GSE213624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213624).

Retrieved 2026-11-03 from the authors' [PPTP-seq_ECOLI repository](https://github.com/yhan0410/PPTP-seq_ECOLI/tree/b8a5f5d7e694de8938bf51a2d905a34326504316), pinned to commit `b8a5f5d7e694de8938bf51a2d905a34326504316`. Supplementary Data 1 and 11 were checked for the target library and experimental strain. The supplement archive is available from [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10505187/supplementaryFiles).

`manifest.json` records source URLs, SHA-256 checksums, input and retained counts, and the selection rule. `significant_responses.csv` retains the original source columns and values for every imported observation, prefixed with its growth medium, source filename, and original CSV line number (header is line 1). Missing upstream values remain blank. This is a selected snapshot, not the complete collection of nonsignificant and unmeasured pairs.

Selection follows the authors' `Notebooks/Figure2_promoter_response.ipynb`: retain class ±1, require natural-log activity standard deviation <0.7, and exclude `(operon, tf_gene)` pairs in `promoter-targeting_sgRNAs.txt`. The downloadable processed tables already omit the listed pairs; keeping the explicit exclusion documents the rule. Apply the consistency filter to each medium. The retained counts are 4058 glucose, 5279 LB, and 3810 glycerol. The glucose count matches Figure 2a exactly. Source values are rounded; classifications are retained as published.

The authors' README calls the `#NAME?` column `-logP`. Their classification notebook defines it as `-log10(q + 1e-10)`, not an unadjusted P-value. Do not convert rounded scores at the ceiling into purported exact q values. `mean` and `std` are natural-log activity statistics; `std_linear` uses a linear scale. `FC`, `FCOM`, `FC_global`, and `FC_specific` are retained as supplied in this snapshot.

To reproduce the source selection, download the three `Processed_data/PPTP-seq_{Glu,LB,Gly}.csv` URLs in the manifest, verify each SHA-256, read in Glu/LB/Gly order, enumerate original data lines starting at 2, and apply the rule above. Medium labels are `M9 glucose (0.4%)`, `LB`, and `M9 glycerol (0.5%)`. Preserve row order and source column order; add the three provenance columns as shown in the snapshot. Use CSV with CRLF line endings to reproduce its checksum.

To import the reviewed snapshot into the website:

```sh
python3 scripts/import_pptp.py
```

The importer verifies its checksum, source counts, classifications, unique observation keys, and existing IDs. It appends IDs above the current maximum and refuses to overwrite an existing PPTP observation whose content differs. Repeated runs are idempotent. Existing observations retain every original field value. The one-to-one alias mappings in the importer were checked against RegulonDB TFSet 14.0, preferring canonical names and single-gene assignments; they do not confer RegulonDB membership on these records. `perturbed_gene` retains the original target even when the table uses an older TF alias.

The raw GEO summary was also downloaded for validation; its higher-precision fold changes and replicate counts were compared against the retained processed observations. It does not expose `FCOM` or the response classes, so it was not used to independently reclassify responses.

The [MIT license](LICENSE.txt) is copied from the authors' repository. The paper and supplementary information are [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Local modifications are selection, metadata/column mapping, TF alias alignment, and assignment of stable record IDs.
