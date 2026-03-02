# Chapter 3: The metabolic diversity of endolithic and soil crust microbial communities

Code for Chapter 3 of my doctoral thesis: "Traces of Microbial Life and Activity in Arctic Endolithic Habitats".

## Contents

| File/Directory | Description |
|:---:|:---:|
| svalbard_map.R | R script for general map of Svalbard used in Figure 2.1. |
| QIIME2_16S.sh | QIIME2 workflow for processing of raw 16S rRNA sequences. |
| QIIME2_18S.sh | QIIME2 workflow for processing of raw 18S rRNA sequences. |
| microeco.R | R script for microbial community analysis as shown in Figures 2.3 - 2.6 using QIIME2 outputs. |
| phyloseq_data_prep.R | R script for creation of phyloseq objects from QIIME2 outputs to use as input for ancombc script. |
| ancombc_16S.R | R script for differential abundance analysis in Figure 2.7 using 16S phyloseq object. |
| ancombc_18S.R | R script for differential abundance analysis in Figure 2.8 using 18S phyloseq object. |
| appendix_a_archaea.R | R script for mean relative abundances of archaea as shown in table A.2 in Appendix A. |
| `elemental analyses` | R scripts for additional elemental analyses that didn't make it into the thesis. |





SP 16S microeco phyloflash - mine, creates rel ab plot

SP 18S '' - mine, creates rel ab plot

mtg-qc.bash - creates 1_Raw and 2_trimmed folders, based on NEOF course. "3_" was taxonomic profiling which i skipped 

humann.bash - runs humann, then? (need to create 1_ and 2_ scripts)

scycdb barplots


FIGURES 
3.1 - rel ab 16S (SP 16S microeco phyloflash, any prerequisite files? needs a phyloseq object)

3.2 - same but 18S

3.3, 3.4 - produced by Laura 

3.5 - scycdb bar plots


