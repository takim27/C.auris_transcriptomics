# Detailed computational protocol

## Purpose and scope

This protocol describes the computational analysis used to discover, validate, quantify, and prioritize putative novel polyadenylated lncRNAs in human PBMCs exposed to live *Candida auris*, purified mannan, purified β-glucan, or RPMI control for 4 h or 24 h.

The experiment contains 24 QuantSeq 3′ mRNA-seq libraries: three donors × four conditions × two time points. Raw data are available from NCBI BioProject [PRJNA647871](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA647871). The reference annotation was GENCODE Release 48 on GRCh38.

This document records the successful workflow. Failed exploratory commands that did not contribute to an output, including the unavailable `CPC2_output_peptide.py` helper, are intentionally omitted. Pfam input proteins were generated successfully with `orfipy`, not with a CPC2 peptide-conversion script.

## Reproducibility status

The original workflow was run on a SLURM cluster and used site-specific absolute paths. Commands below use symbolic project paths. Before public release, place the shell wrappers and Conda YAML files in the repository and replace all cluster account, email, partition, and absolute-path values with configuration variables.

The current R scripts also contain `D:/a/...` paths and interactive `View()` calls. They reproduce the archived analysis in the original environment but require path editing and are not yet a one-command portable pipeline. Do not use `.RData` workspaces as required inputs; every analysis should start from the documented files in a fresh R session.

## Study design

| Factor | Levels |
|---|---|
| Donor | A, B, C |
| Condition | RPMI, live *C. auris*, mannan, β-glucan |
| Time | 4 h, 24 h |
| Libraries | 24 single-end QuantSeq 3′ mRNA-seq libraries |
| Strandedness | Forward-stranded |

The metadata file must contain at least these columns:

```text
sample_id,donor,stimulus,time
```

`sample_id` values must match the cleaned BAM/count-matrix sample names exactly.

## Recommended directory structure

```text
project/
├── reference/
├── metadata/
├── raw_fastq/
├── trimmed_fastq/
├── qc/
│   ├── raw/
│   └── trimmed/
├── star_index/
├── alignments/
├── assembly/
├── candidate_filtering/
│   ├── 01_structure/
│   ├── 02_junction_support/
│   ├── 03_internal_priming/
│   ├── 04_rfam/
│   ├── 05_cpc2/
│   ├── 06_pfam/
│   ├── 07_genomic_context/
│   ├── 08_expression_screen/
│   └── 09_final_catalog/
└── final_analysis/
```

For the examples below, define paths without using a personal home-directory variable:

```bash
export PROJECT_DIR=/path/to/project
export REF_DIR="$PROJECT_DIR/reference"
export RAW_DIR="$PROJECT_DIR/raw_fastq"
export TRIM_DIR="$PROJECT_DIR/trimmed_fastq"
export ALIGN_DIR="$PROJECT_DIR/alignments"
export ASSEMBLY_DIR="$PROJECT_DIR/assembly"
export FILTER_DIR="$PROJECT_DIR/candidate_filtering"
export GENOME="$REF_DIR/GRCh38.gencode_v48.chromosomes.fa"
export GTF="$REF_DIR/gencode.v48.annotation.gtf"
```

## Software

The workflow used FastQC, MultiQC, Cutadapt, STAR, SAMtools, StringTie, GffCompare, gffread, BEDTools, SeqKit, Infernal, Rfam, CPC2, orfipy, HMMER/Pfam, Subread featureCounts, IGV, R, edgeR, limma, UpSetR, clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, ggplot2, pheatmap, ComplexHeatmap, igraph, and writexl.

R 4.5.1 was used for the archived downstream analysis. See [`../final_analysis/R_session_Info.txt`](../final_analysis/R_session_Info.txt) for the recorded R environment. The exact non-R versions should be captured in the final Conda YAML files and a `software_versions.txt` file before archiving.

## Phase I: reference preparation and read processing

### 1. Download GENCODE Release 48

Download and decompress the GRCh38 primary-assembly FASTA and GENCODE Release 48 primary-assembly annotation:

```bash
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz -O "$REF_DIR/GRCh38.primary_assembly.genome.fa.gz"
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.annotation.gtf.gz -O "$REF_DIR/gencode.v48.primary_assembly.annotation.gtf.gz"

gzip -dc "$REF_DIR/GRCh38.primary_assembly.genome.fa.gz" > "$REF_DIR/GRCh38.primary_assembly.genome.fa"
gzip -dc "$REF_DIR/gencode.v48.primary_assembly.annotation.gtf.gz" > "$GTF"
```

Record file dates, sizes, and SHA-256 checksums. The original analysis retained the 25 chromosome sequences represented in the annotation and called the resulting file `GRCh38.gencode_v48.chromosomes.fa`. Verify FASTA/GTF chromosome-name compatibility before indexing.

```bash
cut -f1 "$GTF" | grep -v '^#' | sort -u > "$REF_DIR/gtf_contigs.txt"
samtools faidx "$REF_DIR/GRCh38.primary_assembly.genome.fa"
cut -f1 "$REF_DIR/GRCh38.primary_assembly.genome.fa.fai" | sort -u > "$REF_DIR/fasta_contigs.txt"
comm -23 "$REF_DIR/gtf_contigs.txt" "$REF_DIR/fasta_contigs.txt"
```

The final command should return no GTF contig absent from the FASTA.

### 2. Obtain metadata and FASTQ files

Use the 24 accessions belonging to PRJNA647871 and preserve a tabular accession-to-sample mapping. Download with the SRA Toolkit or ENA, then verify file integrity and sample count.

Example SRA Toolkit pattern:

```bash
while read -r accession; do
  prefetch "$accession"
  fasterq-dump "$accession" --outdir "$RAW_DIR" --threads 8
  gzip "$RAW_DIR/${accession}.fastq"
done < metadata/accessions.txt
```

The public repository should include `metadata/accessions.txt` and the complete sample metadata, but not duplicate the raw FASTQ files.

### 3. Build the STAR genome index

The reads were predominantly 75 nt, so `sjdbOverhang` was set to 74.

```bash
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir "$PROJECT_DIR/star_index" \
  --genomeFastaFiles "$GENOME" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 74
```

### 4. Assess raw-read quality

Run FastQC for all 24 untrimmed libraries and aggregate the reports with MultiQC.

```bash
mkdir -p "$PROJECT_DIR/qc/raw/fastqc"
fastqc --threads 8 --outdir "$PROJECT_DIR/qc/raw/fastqc" "$RAW_DIR"/*.fastq.gz
multiqc "$PROJECT_DIR/qc/raw/fastqc" --outdir "$PROJECT_DIR/qc/raw"
```

Review read counts, length distribution, per-base quality, ambiguous bases, duplication, adapter sequence, and poly(A) signal.

### 5. Trim QuantSeq reads

The original `cutadapt_quantseq.sh` used three sequential Cutadapt operations:

1. Remove 20-nt poly(A) and poly(G) tracts, requiring a 20-nt overlap and retaining reads ≥20 nt.
2. Apply NextSeq quality trimming and remove the 3′ adapter with a minimum 3-nt overlap and maximum error rate 0.1.
3. Remove reads containing a full 5′ adapter match with a 20-nt overlap.

Equivalent parameterization:

```bash
cutadapt -m 20 -O 20 -a 'polyA=A{20}' -a 'QUALITY=G{20}' -n 2 \
  -o stage1.fastq.gz input.fastq.gz

cutadapt -m 20 -O 3 --nextseq-trim=10 \
  -a 'r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.1' \
  -o stage2.fastq.gz stage1.fastq.gz

cutadapt -m 20 -O 20 \
  -g 'r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20' \
  --discard-trimmed \
  -o final_trimmed.fastq.gz stage2.fastq.gz
```

Run this operation independently for each library and retain Cutadapt logs.

### 6. Assess trimmed-read quality

Repeat FastQC and MultiQC on the final trimmed FASTQ files. Confirm the retained read counts and check that adapter/poly(A) content is reduced without an unexpected collapse in read length or quality.

### 7. Align reads with STAR

The original alignment used two-pass STAR and produced coordinate-sorted BAM files, gene-count files, and splice-junction tables.

```bash
STAR \
  --runThreadN 8 \
  --genomeDir "$PROJECT_DIR/star_index" \
  --readFilesIn sample.trimmed.fastq.gz \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --quantMode GeneCounts \
  --outFileNamePrefix "$ALIGN_DIR/sample_"
```

For every sample retain:

- `*_Aligned.sortedByCoord.out.bam`
- `*_Log.final.out`
- `*_ReadsPerGene.out.tab`
- `*_SJ.out.tab`

Verify mapping summaries with MultiQC. Index every BAM:

```bash
samtools index "$ALIGN_DIR/sample_Aligned.sortedByCoord.out.bam"
samtools quickcheck -v "$ALIGN_DIR"/*_Aligned.sortedByCoord.out.bam
```

No output from `samtools quickcheck -v` indicates that the checked BAMs passed the basic integrity test.

### 8. Confirm library orientation

Use the STAR `ReadsPerGene.out.tab` columns to compare unstranded, forward-strand, and reverse-strand assignments. In this dataset, forward-strand counts were substantially greater than reverse-strand counts; therefore, StringTie used `--fr` and featureCounts used `-s 1`.

This is an empirical orientation check, not a universal property of QuantSeq libraries. Reconfirm it whenever the library protocol or STAR settings change.

## Phase II: reference-guided transcript discovery

### 9. Assemble each library with StringTie

For every indexed STAR BAM:

```bash
stringtie sample_Aligned.sortedByCoord.out.bam \
  -p 8 \
  --fr \
  -G "$GTF" \
  -m 200 \
  -o "$ASSEMBLY_DIR/sample.gtf" \
  -A "$ASSEMBLY_DIR/sample.gene_abundance.tsv"
```

Create a text file containing the 24 sample GTF paths and merge them:

```bash
stringtie --merge \
  -p 8 \
  -G "$GTF" \
  -m 200 \
  -o "$ASSEMBLY_DIR/stringtie_merged.gtf" \
  "$ASSEMBLY_DIR/mergelist.txt"
```

### 10. Compare the merged assembly with GENCODE

```bash
gffcompare \
  -V \
  -r "$GTF" \
  -o "$ASSEMBLY_DIR/gffcmp" \
  "$ASSEMBLY_DIR/stringtie_merged.gtf"
```

Retain transcript models carrying GffCompare class code:

- `u`: intergenic relative to the reference annotation;
- `x`: exonic overlap with a reference transcript on the opposite strand.

The initial set contained 1,682 models: 388 `u` and 1,294 `x`. One chrM transcript was removed, leaving 1,681 nuclear candidates.

### 11. Require spliced transcript length ≥200 nt

Use gffread to extract transcript sequences or a transcript table and calculate the sum of exon lengths rather than genomic span. Retain candidates with spliced length ≥200 nt. All 1,681 nuclear candidates passed this requirement.

Separate candidates by exon number:

- 196 multi-exon transcript models;
- 1,485 single-exon transcript models.

The branches were evaluated differently because only multi-exon transcripts have splice junctions that can be validated.

## Phase III: structural and noncoding validation

### 12. Validate multi-exon splice junctions

Confirm that all 24 STAR junction tables are available:

```bash
find "$ALIGN_DIR" -maxdepth 1 -name '*_SJ.out.tab' | sort | wc -l
```

Expected value: 24.

Create candidate intron coordinates:

```bash
JDIR="$FILTER_DIR/02_junction_support"
mkdir -p "$JDIR"

gffread "$FILTER_DIR/01_structure/primary_ux_multiexon.gtf" \
  --table @id,@chr,@strand,@introns \
  > "$JDIR/multiexon_candidate_introns.tsv"
```

Aggregate evidence for every junction across the 24 STAR `SJ.out.tab` files. STAR splice-motif codes 1–6 represent recognized canonical or semi-canonical motifs recorded by STAR; they are not a separate splice-junction database.

An intron passed when all four conditions were met:

- STAR splice-motif code 1–6;
- at least one uniquely mapped supporting read in at least two libraries;
- at least three uniquely mapped supporting reads in total;
- maximum junction overhang ≥12 nt.

The overhang is the aligned portion of a splice-spanning read on one side of the exon–exon junction. Requiring ≥12 nt means that at least one supporting alignment extended 12 or more bases into an adjacent exon, reducing the chance that the junction arose from a very short, ambiguous match.

Retain a multi-exon transcript only when every intron passed. This retained 41 of 196 models and excluded 155.

### 13. Screen for internal poly(A) priming

For each junction-supported multi-exon model and every single-exon model, extract the strand-oriented 20 genomic nucleotides immediately downstream of the inferred transcript 3′ end:

```bash
bedtools getfasta \
  -fi "$GENOME" \
  -bed downstream20.bed \
  -s \
  -nameOnly \
  -fo downstream20.fa
```

Flag a candidate as A-rich when the downstream 20 nt contain either:

- at least six consecutive adenines; or
- at least 12 adenines in total.

Example classification logic:

```awk
seq=toupper($0)
temp=seq
n=gsub(/A/,"",temp)
status=(seq~/AAAAAA/ || n>=12)?"A_RICH":"PASS"
```

After removing A-rich candidates, retain:

- 22 of the 41 junction-supported multi-exon models;
- 654 of the 1,485 single-exon models.

Combining both branches produced 676 structurally supported candidates. These were still provisional candidates, not confirmed lncRNAs.

### 14. Remove structured noncoding RNA families with Rfam

Download the Rfam covariance models and clan information, record the database release/checksum, decompress `Rfam.cm`, and run `cmpress`.

The original `rfam_scan.sh` used Infernal `cmscan` with:

```bash
cmscan \
  --cpu 8 \
  -Z <database_size_Mb> \
  --toponly \
  --cut_ga \
  --rfam \
  --nohmmonly \
  --fmt 2 \
  --clanin Rfam.clanin \
  --tblout rfam.tblout \
  Rfam.cm \
  primary_ux_structural_pass.fa
```

Exclude every candidate with a reported Rfam hit. Fifteen of 676 models were excluded, leaving 661.

`Rfam/CURRENT` is mutable. A public reproducibility release should state the exact Rfam version and checksum rather than relying only on the `CURRENT` URL.

### 15. Evaluate coding potential with CPC2

Run CPC2 on the 661 Rfam-pass transcript sequences:

```bash
CPC2.py \
  -i primary_ux_rfam_pass.fa \
  -o cpc2_results
```

CPC2 appends `.txt` to the requested output prefix. Retain candidates classified as `noncoding` and remove those classified as `coding`. CPC2 retained 659 and excluded two.

### 16. Predict ORFs and search Pfam domains

Pfam searches protein sequences, so ORFs were predicted from the 659 CPC2-noncoding transcript sequences with orfipy. The transcript FASTA sequences were already oriented in the direction of transcription; therefore, search only the forward strand.

```bash
orfipy primary_ux_cpc2_noncoding.fa \
  --pep primary_ux_candidate_orfs.pep.fa \
  --min 30 \
  --strand f \
  --start ATG \
  --partial-3 \
  --table 1 \
  --procs 8 \
  --outdir pfam_results
```

Parameters:

- minimum ORF length: 30 nt;
- start codon: ATG;
- translation table: genetic code 1;
- 3′-partial ORFs permitted because the libraries were 3′ enriched.

For each transcript containing one or more qualifying ORFs, keep the longest translated ORF. Candidates with no qualifying ORF remain in the noncoding candidate set.

Search the longest proteins against Pfam-A release 38.2 using HMMER `hmmscan` and family-specific gathering thresholds:

```bash
hmmpress Pfam-A.hmm

hmmscan \
  --cpu 8 \
  --cut_ga \
  --domtblout pfam.domtblout \
  Pfam-A.hmm \
  primary_ux_longest_orfs.pep.fa \
  > pfam.hmmscan.txt
```

Remove any transcript with a significant Pfam-A hit. One model was excluded, leaving 658 transcripts representing 657 loci.

For exact reproducibility, archive the Pfam-A 38.2 file or its official archived URL and checksum; do not substitute a later `current_release` silently.

## Phase IV: genomic context, expression support, and manual review

### 17. Flag possible coding-gene 3′ extensions

Because QuantSeq preferentially samples transcript 3′ ends, a candidate immediately downstream of a protein-coding gene on the same strand may represent an unannotated coding-gene extension or readthrough product. Construct strand-aware 2-kb windows downstream of annotated protein-coding genes and intersect them with candidate coordinates:

```bash
bedtools intersect \
  -s \
  -wa \
  -wb \
  -a pfam_pass_candidates.bed \
  -b coding_gene_3prime_2kb_windows.bed \
  > coding_3prime_risk_pairs.tsv
```

The 2-kb window is a practical review boundary, not a universal biological cutoff. Flagged candidates were not removed automatically.

This screen flagged 34 transcripts representing 33 loci; 624 candidates were considered lower risk.

### 18. Generate the provisional expression-screen matrix

Create a provisional annotation containing the complete GENCODE Release 48 GTF plus all 658 Pfam-pass candidate transcripts:

```bash
cat "$GTF" primary_ux_pfam_pass.gtf > provisional_annotation_658.gtf
```

Count all 24 BAMs by `gene_id` with forward-strand featureCounts:

```bash
featureCounts \
  -T 8 \
  -s 1 \
  -t exon \
  -g gene_id \
  -a provisional_annotation_658.gtf \
  -o screening_gene_counts.txt \
  "$ALIGN_DIR"/*_Aligned.sortedByCoord.out.bam
```

This command is the origin of `expression_screen/screening_gene_counts.txt`. It is a provisional raw-count matrix and is distinct from the final count matrix used for differential expression.

### 19. Apply donor-supported expression filtering

Run [`../expression_screen/expression_screen.R`](../expression_screen/expression_screen.R) with:

- CPM threshold: ≥0.5;
- biological support: at least two distinct donors;
- grouping: the same stimulus–time combination.

A candidate locus passed if it met the CPM threshold in at least two donors within at least one of the eight stimulus–time groups. Of 657 loci, 613 passed and 44 failed. The 613 loci corresponded to 614 transcript models because one locus contained two retained isoforms.

### 20. Select the expression-supported risk subset

Intersect the expression-pass locus IDs with the 33 risk loci. One risk locus failed the expression filter, leaving 33 transcripts representing 32 loci for IGV review. The review table contained 34 transcript–coding-gene relationships because one transcript had two coding-gene partners.

Prepare:

- a GTF containing the risk candidates and nearby coding genes;
- a table of searchable IGV regions;
- a BED file of merged review regions;
- 24 reduced BAM files containing reads in those regions;
- a BAI index for every reduced BAM.

Use `samtools view -bh -L review_regions.bed` to create the reduced BAMs, then index and validate them with `samtools quickcheck`.

### 21. Perform manual IGV review

IGV 2.19.8 was used with the hg38 genome. For each flagged relationship, inspect all 24 reduced BAM tracks and the candidate/coding-gene annotation.

Record one of four decisions:

- `RETAIN`: a separate, reproducible signal with an observable gap from the coding gene;
- `EXCLUDE_READTHROUGH`: continuous coverage from the coding gene into the candidate;
- `EXCLUDE_EXTENSION`: the candidate appears to extend the terminal coding exon;
- `AMBIGUOUS`: insufficient evidence for a confident decision.

The final strict exclusion list contained six unique transcripts:

```text
MSTRG.1803.1
MSTRG.17016.1
MSTRG.22578.15
MSTRG.34694.1
MSTRG.36939.1
MSTRG.13391.1
```

MSTRG.36939.1 appeared twice in the relationship table because it was compared with two coding genes, but it was excluded only once.

For reproducibility, retain the review table, decision for every row, reviewer identity, IGV version, displayed tracks, and representative screenshots. Manual review is a judgment-based step and should remain auditable.

### 22. Build the final novel-lncRNA catalog

Remove the six IGV-excluded transcripts from the 614 expression-supported models. The final catalog contains 608 transcript models representing 607 loci.

Extract the final GTF with gffread and add explicit attributes when absent:

```text
gene_type "novel_lncRNA";
transcript_type "novel_lncRNA";
```

Append the labeled novel GTF to the complete GENCODE Release 48 annotation to create:

```text
final_annotation_GENCODEv48_plus_novel_lncRNA.gtf
```

Verify that the combined annotation contains 608 novel transcript models, 607 novel loci, and no chrM novel models.

## Phase V: final quantification and expression analysis

### 23. Generate the final count matrix

Run featureCounts against the combined final annotation:

```bash
featureCounts \
  -T 8 \
  -s 1 \
  -t exon \
  -g gene_id \
  -a final_annotation_GENCODEv48_plus_novel_lncRNA.gtf \
  -o final_gene_counts.txt \
  "$ALIGN_DIR"/*_Aligned.sortedByCoord.out.bam
```

Confirm 24 sample-count columns and generate a gene-annotation table containing `gene_id`, `gene_name`, `gene_type`, chromosome, and strand. The final raw matrix contained 79,293 gene-level features, including 607 novel loci.

### 24. Filter and normalize expression values

In [`../final_analysis/analysis_part_1.R`](../final_analysis/analysis_part_1.R):

```r
metadata$group <- interaction(
  metadata$stimulus,
  metadata$time,
  sep = "_",
  drop = TRUE
)

dge <- edgeR::DGEList(counts = counts)
keep <- edgeR::filterByExpr(dge, group = metadata$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge)
```

`filterByExpr()` used its edgeR defaults with the eight-level stimulus–time group supplied. It retained 13,321 genes, including 449 novel lncRNA loci. These 449 are expression-filtered loci, not differentially expressed loci.

Use normalized log2-CPM values with prior count 2 for PCA and exploratory expression visualization:

```r
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
```

### 25. Fit the donor-adjusted edgeR model

```r
design <- model.matrix(~ 0 + group + donor, data = metadata)
colnames(design) <- make.names(colnames(design))

dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
```

Define six stimulus-versus-time-matched-RPMI contrasts:

```r
contrasts <- limma::makeContrasts(
  KTCbg_4h    = groupcaur_KTCbg_4h    - groupRPMI_4h,
  KTClive_4h  = groupcaur_KTClive_4h  - groupRPMI_4h,
  KTCman_4h   = groupcaur_KTCman_4h   - groupRPMI_4h,
  KTCbg_24h   = groupcaur_KTCbg_24h   - groupRPMI_24h,
  KTClive_24h = groupcaur_KTClive_24h - groupRPMI_24h,
  KTCman_24h  = groupcaur_KTCman_24h  - groupRPMI_24h,
  levels = design
)
```

For each contrast, run `glmQLFTest()`. Adjust P values across all tested genes in that contrast using the Benjamini–Hochberg procedure returned by edgeR. A gene was classified as differentially expressed when:

- FDR <0.05; and
- |log2 fold change| >1.

Across the six contrasts, 102 unique novel lncRNAs met these criteria in at least one contrast.

### 26. Compare DE sets with UpSet analysis

Create membership sets for significant genes in all six contrasts and visualize their intersections with UpSet plots. The R script includes analyses of all DE genes and novel-lncRNA subsets, including direction-specific sets. If only selected UpSet panels are reported in the manuscript, retain the remaining outputs in the repository and label them as supplementary or unreported exploratory results.

## Phase VI: exploratory functional analyses

### 27. Perform donor-adjusted co-expression analysis

Construct normalized expression matrices for:

- 102 unique DE novel lncRNAs;
- 2,125 unique DE protein-coding genes.

Remove donor effects while preserving the eight stimulus–time groups:

```r
logCPM_coexpression <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
group_design <- model.matrix(~ 0 + group, data = metadata)

logCPM_donor_adjusted <- limma::removeBatchEffect(
  logCPM_coexpression,
  batch = metadata$donor,
  design = group_design
)
```

Calculate all 102 × 2,125 = 216,750 Pearson correlations. Because two donor indicator parameters were removed, the original script used:

```r
correlation_df <- 24 - (nlevels(metadata$donor) - 1) - 2
```

Convert correlations to two-sided t-test P values, adjust all 216,750 P values together by Benjamini–Hochberg, and retain pairs with:

- |r| ≥0.70;
- FDR <0.05.

The checkpoint results were 10,866 selected pairs involving 94 novel lncRNAs and 1,774 protein-coding genes: 5,880 positive and 4,986 negative pairs.

These correlations are exploratory associations and do not establish regulation or causality.

### 28. Perform functional enrichment

Rank lncRNAs by the number of significant coding-gene partners and select the 15 most highly connected lncRNAs for enrichment analyses. For the combined analysis, use their 1,370 unique coding partners as the target set and the 2,125 DE coding genes as the background universe.

Map Ensembl IDs through `org.Hs.eg.db`, then run:

- GO Biological Process enrichment with `clusterProfiler::enrichGO()`;
- Reactome enrichment with `ReactomePA::enrichPathway()`;
- KEGG enrichment with `clusterProfiler::enrichKEGG()`.

The combined GO settings were:

```r
enrichGO(
  gene = target_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)
```

Reactome and KEGG used BH adjustment, `pvalueCutoff=0.05`, `qvalueCutoff=0.05`, and gene-set sizes of 10–500. KEGG used `organism="hsa"`, NCBI Gene IDs, and `use_internal_data=FALSE`. The script also performs lncRNA-specific enrichment for the 15 highly connected candidates.

### 29. Test stimulus-by-time interactions

Use the same fitted donor-adjusted edgeR model and define, for each stimulus:

```text
(stimulus24 − RPMI24) − (stimulus4 − RPMI4)
```

In R:

```r
time_contrasts <- limma::makeContrasts(
  KTCbg_time_change =
    (groupcaur_KTCbg_24h - groupRPMI_24h) -
    (groupcaur_KTCbg_4h - groupRPMI_4h),
  KTClive_time_change =
    (groupcaur_KTClive_24h - groupRPMI_24h) -
    (groupcaur_KTClive_4h - groupRPMI_4h),
  KTCman_time_change =
    (groupcaur_KTCman_24h - groupRPMI_24h) -
    (groupcaur_KTCman_4h - groupRPMI_4h),
  levels = design
)
```

Retain interaction results with FDR <0.05 and |interaction log2 fold change| >1. The analysis identified 19 unique novel lncRNAs with a significant temporal interaction; 14 were also among the 102 stimulus-versus-RPMI DE novel lncRNAs and were carried into integrated prioritization. Four of those 14—MSTRG.6308, MSTRG.32011, MSTRG.25627, and MSTRG.22879—were also among the 15 most highly connected lncRNAs.

## Phase VII: integrated prioritization

### 30. Screen for cis-proximal coding genes

Use the 14 DE-and-temporally responsive lncRNAs. Compare their genomic coordinates with protein-coding genes on the same chromosome and calculate interval-to-interval distance, assigning zero to overlaps. Retain pairs separated by ≤100,000 bp.

Checkpoints:

- 14,930 same-chromosome lncRNA–coding-gene combinations;
- 60 pairs within 100 kb;
- those 60 pairs involved 13 lncRNAs and 56 unique coding genes;
- 13 proximal pairs contained a coding gene belonging to the 2,125 DE coding-gene set.

Retrieve the co-expression statistics for those 13 analyzable proximal pairs and retain pairs meeting |r| ≥0.70 and FDR <0.05. Four candidate cis-associated pairs passed:

```text
MSTRG.32011–TLR2
MSTRG.5763–DUSP5
MSTRG.25627–HCK
MSTRG.18159–ITGB3
```

Proximity plus co-expression supports candidate cis association; it does not demonstrate direct cis regulation.

### 31. Construct the integrated co-expression network

From the global significant co-expression table, extract edges involving the 14 DE-and-temporal candidate lncRNAs. Annotate edges by correlation direction, coding-gene name, 100-kb proximity, distance, strand relationship, and membership in the four high-confidence cis pairs.

The full network checkpoint is:

- 1,744 significant co-expression edges;
- 13 connected novel lncRNAs;
- 1,149 unique coding genes;
- 896 positive and 848 negative edges;
- MSTRG.6052 had no qualifying edge.

For the reduced display network, sort coding partners for each lncRNA by FDR and then absolute correlation, retain the top five partners per connected lncRNA, and add any of the four high-confidence cis edges not already selected. Remove duplicate edges.

The reduced display checkpoint is 77 nodes (13 lncRNAs and 64 coding genes) and 66 edges.

Export full and reduced node/edge tables for Cytoscape:

```text
candidate_14_Cytoscape_edges.csv
candidate_14_Cytoscape_nodes.csv
candidate_14_display_edges.csv
candidate_14_display_nodes.csv
```

## Expected end-to-end checkpoints

| Stage | Retained | Excluded at stage |
|---|---:|---:|
| GffCompare `u`/`x` models | 1,682 | — |
| Nuclear candidates | 1,681 | 1 chrM transcript |
| Multi-exon splice-junction support | 41 of 196 | 155 multi-exon models |
| Internal-priming screen | 676: 22 multi-exon + 654 single-exon | 850 models |
| Rfam | 661 | 15 |
| CPC2 noncoding classification | 659 | 2 |
| Pfam | 658 transcripts / 657 loci | 1 |
| Donor-supported expression | 614 transcripts / 613 loci | 44 loci |
| Manual IGV review | 608 transcripts / 607 loci | 6 transcripts |

## Quality-control and audit requirements

For every public or archived run, retain:

- input URLs, release numbers, file sizes, and SHA-256 checksums;
- the 24-accession/sample metadata table;
- Conda YAML files and `software_versions.txt`;
- command logs and SLURM logs;
- FastQC and MultiQC reports before and after trimming;
- STAR `Log.final.out`, `ReadsPerGene.out.tab`, and `SJ.out.tab` files;
- intermediate candidate ID lists and GTFs at every exclusion step;
- the IGV review table, decisions, version, and representative screenshots;
- raw featureCounts outputs and summaries;
- R scripts, `sessionInfo()`, and generated tables/figures;
- a machine-readable manifest mapping every manuscript figure/table to the generating script and source file.

## Known limitations of the workflow

- The experiment contains three donors, so correlation and interaction results are exploratory.
- QuantSeq is 3′ enriched and is not optimized for reconstructing full-length transcript structures; single-exon candidates require especially cautious interpretation.
- Manual IGV review introduces reviewer judgment and should be independently checked when feasible.
- The 2-kb 3′-extension flag and 100-kb cis window are operational thresholds, not universal biological boundaries.
- Co-expression and proximity do not establish direct regulatory activity.
- Rfam, Pfam, annotation, and pathway databases change over time; releases and checksums must be frozen for exact reproduction.

## Minimum steps before GitHub/Zenodo release

1. Add all successful shell scripts under `scripts/` and replace absolute paths with a shared configuration file.
2. Add the Conda environment YAML files.
3. Add the accession list and reference/database checksum manifest.
4. Refactor `analysis_part_1.R` so it runs non-interactively from a clean R session and uses repository-relative paths.
5. Add a license, `CITATION.cff`, and tagged release.
6. Archive that release in Zenodo and cite its DOI in the manuscript’s availability statement.
