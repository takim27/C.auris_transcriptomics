# Novel lncRNA responses to *Candida auris* stimulation

This repository documents the computational workflow used to discover and prioritize previously unannotated polyadenylated long noncoding RNAs (lncRNAs) in human peripheral blood mononuclear cells (PBMCs) exposed to *Candida auris* or fungal cell-wall components.

The analysis reuses 24 QuantSeq 3′ mRNA-sequencing libraries from NCBI BioProject [PRJNA647871](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA647871): three donors, four conditions (RPMI control, live *C. auris*, mannan, and β-glucan), and two time points (4 h and 24 h). The source study is Bruno et al., *Nature Microbiology* (2020), [doi:10.1038/s41564-020-0780-3](https://doi.org/10.1038/s41564-020-0780-3).

## Repository contents

```text
C.auris_transcriptomics/
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

The [detailed_methods.md](documentation/detailed_methods.md) describes the complete workflow from public FASTQ files through transcript discovery, candidate validation, expression analysis, functional analysis, cis-proximity screening, and network construction.


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

### Requirements

- R 4.5.1 was used for the archived analysis.
- Bioconductor/CRAN packages are listed in [`final_analysis/R_session_Info.txt`](final_analysis/R_session_Info.txt).
- The main packages include edgeR, limma, ggplot2, UpSetR, pheatmap, ComplexHeatmap, clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, igraph, and writexl.

### Expression-support screen

Inputs are already grouped under [`expression_screen/`](expression_screen/). The script applies the rule CPM ≥0.5 in at least two distinct donors within at least one identical stimulus-time group.

Before running, change `base_dir` near the top of [`expression_screen.R`](expression_screen/expression_screen.R) to the local `expression_screen` directory. Then run:

```bash
Rscript expression_screen/expression_screen.R
```

### Final statistical and functional analysis

Inputs are under [`final_analysis/`](final_analysis/). Before running [`analysis_part_1.R`](final_analysis/analysis_part_1.R), replace its computer-specific `D:/a/final_analysis` paths with the local `final_analysis` directory. The script contains `View()` calls and is currently best run section by section in an interactive R session.

The principal statistical settings are:

- `filterByExpr(dge, group=metadata$group)` using the edgeR defaults;
- TMM normalization with `calcNormFactors()`;
- donor-adjusted design `~0 + group + donor`;
- robust dispersion estimation and quasi-likelihood fitting;
- six stimulus-versus-time-matched-RPMI contrasts;
- differential expression at FDR <0.05 and |log2 fold change| >1;
- Pearson co-expression at |r| ≥0.70 and BH FDR <0.05;
- stimulus-by-time interactions at FDR <0.05 and |interaction log2 fold change| >1;
- cis proximity defined as ≤100 kb.

## Important distinction between the two count matrices

`expression_screen/screening_gene_counts.txt` is a provisional featureCounts matrix generated from GENCODE Release 48 plus all 658 Pfam-pass candidate transcripts. It was used only to identify donor-supported candidate loci before IGV review.

`final_analysis/final_gene_counts.txt` is the final raw gene-level featureCounts matrix generated from GENCODE Release 48 plus the manually retained 608 transcripts representing 607 loci. It contains 79,293 feature rows before edgeR expression filtering. The 449 novel loci were obtained after `filterByExpr()`; they were not produced by differential-expression testing.

## Current reproducibility status

The count matrices, metadata, final annotation table, R scripts, outputs, and detailed protocol are present. Before a public GitHub/Zenodo release, the following items should also be added:

- path-independent shell scripts for FASTQ retrieval, trimming, alignment, assembly, candidate filtering, and featureCounts;
- the Conda environment YAML files used for the core and lncRNA-specific tools;
- an accession list and checksums for downloaded reference/input files;
- a small configuration file defining project and reference paths;
- a non-interactive version of the final R analysis with `View()` calls removed or guarded;
- a license and a citation file.

The existing scripts were originally written for a specific SLURM/HPC directory structure. Their scientific commands are documented in the detailed protocol, but the repository should not be described as one-command reproducible until those portability items are completed.

## Data availability and archiving

Raw sequencing data should not be duplicated in this repository; they are available from NCBI under PRJNA647871. Large derived files may be deposited in Zenodo or another stable repository and linked here. For manuscript submission, archive a tagged GitHub release in Zenodo and cite the resulting DOI so the exact code version remains accessible.

## Citation

If this workflow is used before the associated manuscript is published, cite the source RNA-seq study and this repository/Zenodo record. Replace this section with the final manuscript citation and repository DOI after acceptance or public archiving.

## License

No license has yet been specified. Add an appropriate open-source license before public release; without one, reuse rights remain restricted by default.
























