# ChIP-seq (PhoB) source provenance

Fitzgerald DM, Stringer AM, Smith C, Lapierre P, Wade JT. *mBio* 14(3):e02535-22 (2023), [Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10294691/#tab1), [doi:10.1128/mbio.02535-22](https://doi.org/10.1128/mbio.02535-22). Table data adapted under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The `CHIP` entries were extracted directly from Table 1 of this paper. The table reports 65 PhoB-bound regions in *E. coli* K-12 MG1655 phoB-FLAG3 (DMF34), grown under low-phosphate conditions (0.2 mM K₂HPO₄). Each database entry represents one ChIP-seq peak for PhoB. The associated genes are genomic neighbors or overlapping genes; binding does not by itself establish transcriptional regulation.

Table 1 provides 59 individual DNA sequences matching a PhoB motif identified by MEME from the ChIP-seq regions. These are motif instances, not 59 distinct motif models. The six peaks marked "NA" in the paper retain a blank DNA sequence in the database. The protein binding-domain and effector fields are blank because Table 1 does not provide them. `characterized` describes PhoB, not the regulatory function of each binding site.

`peak_center` is the center of the enriched region in the MG1655 reference genome **NC_000913.3**, not the start of the listed motif. `chip_score` is the paper's relative ChIP enrichment score. Parentheses in `associated_genes` indicate an intragenic site, as in the original table. The `source` column links to the paper's table. No RegulonDB records were used for this import; the RegulonDB membership field is left blank for these entries. Raw ChIP-seq data are deposited under [E-MTAB-9293](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-9293).

## Motif analysis

The paper's "PhoB motif discovery and analysis" section reports **MEME v5.1.0** analysis of 100-bp regions surrounding ChIP-seq peaks, with default parameters. **Figure 2B** shows the PhoB motif: matches occur in 59 of 65 regions, with MEME E-value **3.0 × 10⁻⁷²**. **Figure 2C** uses **CentriMo v5.1.0** through MEME-ChIP to test central enrichment relative to peak centers (E-value **1.2 × 10⁻¹²**). **Figure 3** repeats motif discovery for peaks shared with, or unique relative to, an earlier ChIP-chip study. The authors also use **FIMO v5.1.0** to scan the genome with a PWM built from scrambled binding sites as a control. [Paper](https://doi.org/10.1128/mbio.02535-22).
