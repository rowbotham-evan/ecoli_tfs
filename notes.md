<!-- New entries go at the top. -->

## 2026-10-27

Started the database.



## Methods Compilation 

The following methods are catalogued in the database. 



### DAP-seq



### gSelex-seq

Method Procedure: 
1. Extract DNA from E. Coli
2. Mechnically break up the DNA using sonication or other techniques
3. Select DNA fragments that are 200 - 300 bp long. There's usually a lot of fragments, because redunancy with fragments improves genome coverage. 
4. 



\[1\] [TEC Database](https://shigen.nig.ac.jp/ecoli/tec/top/aboutSelex)


### CHIP-seq 

#### Procedure

1. **Tag the TF with an epitope.** Native-locus tagging: SPA or 3×FLAG tags introduced by recombineering in strain DY330, moved into MG1655 by P1 transduction, KanR marker removed with FLP recombinase; V5 tags added with the FRUIT method. Alternatively, inducible expression: the TF coding sequence cloned into pBAD24 fused to 3×FLAG (arabinose-inducible), validated by Sanger sequencing.
2. **Grow under one standardized condition.** M9 minimal medium + 0.4% glycerol (plus 100 µg/mL ampicillin for plasmid strains), 30 °C with shaking; subculture to OD ~0.5. For inducible strains, add arabinose to 0.2% and grow 1 h to induce the TF before crosslinking. A single medium is used for all TFs so binding profiles are comparable across the compendium.
3. **Crosslink.** Add formaldehyde to 1% final, shake 30 min at room temperature to fix TFs onto bound DNA; quench with glycine (250 mM final, 15 min). Pellet, wash twice in PBS, and resuspend in HEPES lysis buffer with protease inhibitor.
4. **Shear the DNA.** Sonicate on a Covaris S2 (12 min, 20% amplitude, intensity 5, 200 cycles/burst), giving a tight fragment distribution centered on ~250 bp — important for downstream resolution. Clear the lysate by centrifugation and salt-adjust the supernatant (10 mM Tris-HCl pH 8, 150 mM NaCl, 0.1% NP-40).
5. **Immunoprecipitate.** Add 5 µg of a ChIP-grade monoclonal antibody against the tag (anti-FLAG, Sigma #F1804, for SPA/3×FLAG; anti-V5, Sigma #V8012) and rock overnight at 4 °C. Capture complexes on Protein G agarose beads; wash five times with IPP150 and twice with TE to remove unbound lysate.
6. **Purify the DNA.** Elute TF–DNA complexes from the beads in SDS elution buffer at 65 °C; digest protein with Proteinase K (37 °C, 1 h); reverse the crosslinks by overnight incubation at 65 °C; column-purify the DNA.
7. **Prepare sequencing libraries** (NEBNext Ultra II): end repair + dA-tailing, adapter ligation with USER treatment, AMPure XP cleanup, 14-cycle indexing PCR, and double-sided size selection to ~350 bp library fragments. Verify size and concentration on a Bioanalyzer, normalize to 4 nM, and pool.
8. **Sequence** on Illumina NextSeq 500/2000, 75 bp single-end reads (paired-end in a few cases).
9. **Analyze.** Trim adapters, align to the U00096.3 reference, and compute strand-specific coverage. Call enriched regions (contiguous stretches > 150 bp above background), then filter on the expected forward/reverse strand coverage signature, enrichment in no-tag/control experiments, and known artifact regions. An experiment passes QC only if ≥ 90% of genome positions have coverage and the strongest region is ≥ 10-fold enriched; a TF is retained with at least two passing experiments. Overlapping regions are merged and kept when they appear in ≥ 2 experiments, or in one experiment at a RegulonDB-known site.

\[2\] [Genome-Wide Mapping of the Escherichia coli PhoB Regulon Reveals Many Transcriptionally Inert, Intragenic Binding Sites.](https://pubmed.ncbi.nlm.nih.gov/37067422/)

\[3\] [Predictive biophysical neural network modeling of a compendium of in vivo transcription factor DNA binding profiles for Escherichia coli.](https://www.nature.com/articles/s41467-025-58862-8)


### PPTP-Seq

\[4\] [Genome-wide promoter responses to CRISPR perturbations of regulators reveal regulatory networks in Escherichia coli](https://www.nature.com/articles/s41467-023-41572-4)


### IPOD-HR

\[5\] [Genomic footprinting uncovers global transcription factor responses to amino acids in Escherichia coli](https://pubmed.ncbi.nlm.nih.gov/37820729/)
