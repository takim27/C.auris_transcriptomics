# Novel lncRNA responses to *Candida auris* stimulation

This repository documents the computational workflow used to discover and prioritize previously unannotated polyadenylated long noncoding RNAs (lncRNAs) in human peripheral blood mononuclear cells (PBMCs) exposed to *Candida auris* or fungal cell-wall components.

The analysis reuses 24 QuantSeq 3′ mRNA-sequencing libraries from NCBI BioProject [PRJNA647871](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA647871): three donors, four conditions (RPMI control, live *C. auris*, mannan, and β-glucan), and two time points (4 h and 24 h). The source study is Bruno et al., *Nature Microbiology* (2020), [doi:10.1038/s41564-020-0780-3](https://doi.org/10.1038/s41564-020-0780-3).

## Repository contents

```text
C.auris_transcriptomics/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── docs/
│   ├── detailed_method.md
│   └── detailed_method.pdf
│
├── environments/
│   └── Conda environment YAML files
│
├── scripts/
│   └── Bash and SLURM scripts for sequence processing and lncRNA identification
│
└── R/
    ├── R_session_Info.txt
    │
    ├── expression_screening/
    │   ├── expression_screen.R
    │   ├── candidate_gene_ids.txt
    │   ├── screening_gene_counts.txt.gz
    │   └── sample_metadata.csv
    │
    └── differential_expression_and_downStream_analysis/
        ├── differential_expression_and_downStream_analysis.R
        ├── final_gene_annotation.tsv
        ├── final_gene_counts.txt.gz
        └── sample_metadata.csv
```

The files in the [docs/](docs/) folder describe the complete workflow from public FASTQ files through transcript discovery, candidate validation, expression analysis, functional analysis, cis-proximity screening, and network construction.


## Analysis overview

1. Raw reads were assessed with FastQC/MultiQC, trimmed with Cutadapt, and aligned to GRCh38 with STAR.
2. Reference-guided StringTie assemblies were merged and compared with GENCODE Release 48 using GffCompare.
3. Intergenic (`u`) and antisense exonic-overlap (`x`) transcripts were screened for length, splice-junction support, internal priming, noncoding RNA families, coding potential, and Pfam domains.
4. Donor-supported expression and manual IGV review produced a final catalog of 608 transcript models representing 607 novel lncRNA loci.
5. The novel catalog was appended to GENCODE Release 48 and quantified with featureCounts.
6. edgeR was used for expression filtering, TMM normalization, donor-adjusted differential-expression testing, and stimulus-by-time interaction testing.
7. Significant novel lncRNAs were evaluated by co-expression, pathway enrichment, genomic proximity, and an integrated network.

## Principal checkpoints

| Checkpoint | Expected result |
|---|---:|
| Initial `u`/`x` transcript models | 1,682 |
| Nuclear candidates after removal of one chrM transcript | 1,681 |
| Multi-exon / single-exon candidates | 196 / 1,485 |
| Structurally supported candidates after junction and internal-priming filters | 676 |
| Candidates after Rfam, CPC2, and Pfam screening | 658 transcripts / 657 loci |
| Expression-supported candidates | 614 transcripts / 613 loci |
| Final manually reviewed catalog | 608 transcripts / 607 loci |
| Features in the final featureCounts matrix | 79,293 |
| Genes retained by `filterByExpr()` | 13,321 |
| Retained novel lncRNA loci | 449 |
| Novel lncRNAs significant in at least one stimulus-versus-RPMI contrast | 102 |

These values are validation checkpoints, not parameters to force. A different reference release, software version, or downloaded dataset revision may change them.


## Reproducing the R analyses
Two .R files in the [R/](R/) folders contain the expression screening analysis, differential expression analysis, UpSetR overlap, enrichment, co-expression, cis-proximity, and co-expression-based network analysis.
### Requirements

- R 4.5.1 was used for the archived analysis.
- Bioconductor or CRAN packages are listed in [R_session_Info.txt](R/R_session_Info.txt)
- The main packages include edgeR, limma, ggplot2, UpSetR, pheatmap, ComplexHeatmap, clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, igraph, and writexl.

### Expression-support screening

The files required for expression-support screening are grouped under [`R/expression_screening/`](R/expression_screening/). The [`expression_screen.R`](R/expression_screening/expression_screen.R) script applies the criterion CPM ≥0.5 in at least two distinct donors within at least one identical stimulus–time group.

Before running the script, change `base_dir` near the beginning of the file to the local path of the `R/expression_screening/` directory. From the repository root, run:

```bash
Rscript R/expression_screening/expression_screen.R
```

This directory contains:

* [`candidate_gene_ids.txt`](R/expression_screening/candidate_gene_ids.txt)
* [`screening_gene_counts.txt.gz`](R/expression_screening/screening_gene_counts.txt.gz)
* [`sample_metadata.csv`](R/expression_screening/sample_metadata.csv)
* [`expression_screen.R`](R/expression_screening/expression_screen.R)

### Differential-expression and downstream analyses

The files required for the final statistical and functional analyses are grouped under [`R/differential_expression_and_downStream_analysis/`](R/differential_expression_and_downStream_analysis/).

Before running [`differential_expression_and_downStream_analysis.R`](R/differential_expression_and_downStream_analysis/differential_expression_and_downStream_analysis.R), replace its computer-specific `D:/a/final_analysis` paths with the local path of the `R/differential_expression_and_downStream_analysis/` directory. 

* [`differential_expression_and_downStream_analysis.R`](R/differential_expression_and_downStream_analysis/differential_expression_and_downStream_analysis.R)
* [`final_gene_annotation.tsv`](R/differential_expression_and_downStream_analysis/final_gene_annotation.tsv)
* [`final_gene_counts.txt.gz`](R/differential_expression_and_downStream_analysis/final_gene_counts.txt.gz)
* [`sample_metadata.csv`](R/differential_expression_and_downStream_analysis/sample_metadata.csv)

The principal statistical settings are:

* `filterByExpr(dge, group = metadata$group)` using the edgeR defaults;
* TMM normalization using `calcNormFactors()`;
* donor-adjusted design `~0 + group + donor`;
* robust dispersion estimation and quasi-likelihood fitting;
* six stimulus-versus-time-matched-RPMI contrasts;
* differential expression at FDR <0.05 and |log2 fold change| >1;
* Pearson co-expression at |r| ≥0.70 and BH FDR <0.05;
* stimulus-by-time interactions at FDR <0.05 and |interaction log2 fold change| >1;
* cis proximity defined as ≤100 kb.

The R version and package information used for these analyses are provided in [`R/R_session_Info.txt`](R/R_session_Info.txt).

## Important distinction between the two count matrices

[`R/expression_screening/screening_gene_counts.txt.gz`](R/expression_screening/screening_gene_counts.txt.gz) is the provisional featureCounts matrix generated using GENCODE Release 48 together with all 658 Pfam-pass candidate transcripts. It was used only to identify donor-supported candidate loci before manual IGV review.

[`R/differential_expression_and_downStream_analysis/final_gene_counts.txt.gz`](R/differential_expression_and_downStream_analysis/final_gene_counts.txt.gz) is the final raw gene-level featureCounts matrix generated using GENCODE Release 48 together with the manually retained 608 transcripts representing 607 loci. It contains 79,293 feature rows before edgeR expression filtering.

The 449 expressed novel loci were retained after `filterByExpr()` and were not identified through differential-expression testing.

## Current reproducibility status

The repository currently contains:

* the expression-screening R script and its required count matrix, candidate IDs, and sample metadata;
* the differential-expression and downstream-analysis R script with its required final count matrix, annotation table, and sample metadata;
* the R session information;
* detailed documentation under [`docs/`](docs/);
* computational-environment materials under [`environments/`](environments/);
* sequence-processing and lncRNA-identification scripts under [`scripts/`](scripts/);


The scripts were originally developed for specific Windows and SLURM/HPC directory structures. Although the scientific procedures are documented, the repository should not be described as fully one-command reproducible until the remaining absolute paths and interactive commands have been addressed.























