# Pangenome analysis reveals the evolutionary dynamics of repeat-based holocentromeres

Scripts used for the preparation of manuscript available at:

https://www.biorxiv.org/content/10.64898/2026.01.17.700053v1

## Directories

**Repeat_analysis/** - R scripts for processing and analyzing repeat elements across multiple genomes. Scripts include: data filtering and quality control (TYBA, EDTA transposable element annotations), comparative synteny analysis of repeat sequences, identification of higher-order repeat (HOR) arrays and islands, gap analysis, and generation of statistical summaries and publication-quality visualizations. Helper functions are provided in `aux_fun.R`.

**Simulations/** - Python and Fortran scripts implementing molecular dynamics simulations of chromatin structure and topology. Includes generation of simulated higher-order repeat arrays with specified sizes and spacings, analysis of array structures, and computational prediction of chromatin properties such as loop size equilibration and three-dimensional fiber width.