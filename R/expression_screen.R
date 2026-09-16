# ============================================================
# Novel-lncRNA expression screening
# Study design: 3 donors × 4 stimuli × 2 timepoints = 24 samples
# Expression rule: CPM >= 0.5 in at least 2 distinct donors
# within at least one identical stimulus-time group
# ============================================================

# ------------------------------------------------------------
# 0. Install and load edgeR
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("edgeR", quietly=TRUE)) {
  BiocManager::install("edgeR")
}

library(edgeR)

# ------------------------------------------------------------
# 1. Define directories and filtering thresholds
# ------------------------------------------------------------

base_dir <- "D:/a/expression_screen"
output_dir <- file.path(base_dir, "results")

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

cpm_threshold <- 0.5
minimum_donors <- 2

# ------------------------------------------------------------
# 2. Define and verify input files
# ------------------------------------------------------------

count_file <- file.path(base_dir, "screening_gene_counts.txt")
candidate_file <- file.path(base_dir, "candidate_gene_ids.txt")
metadata_file <- file.path(base_dir, "sample_metadata.csv")

if (!file.exists(count_file)) {
  stop("screening_gene_counts.txt was not found")
}

if (!file.exists(candidate_file)) {
  stop("candidate_gene_ids.txt was not found")
}

if (!file.exists(metadata_file)) {
  stop("sample_metadata.csv was not found")
}

# ------------------------------------------------------------
# 3. Import featureCounts output
# ------------------------------------------------------------

fc <- read.delim(
  count_file,
  comment.char="#",
  check.names=FALSE,
  stringsAsFactors=FALSE
)

cat("FeatureCounts dimensions:", nrow(fc), "rows and", ncol(fc), "columns\n")

required_fc_columns <- c(
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

if (!all(required_fc_columns %in% colnames(fc))) {
  stop("The featureCounts annotation columns are incomplete")
}

if (ncol(fc) <= 6) {
  stop("No sample-count columns were detected")
}

head(fc[, 1:7])

# ------------------------------------------------------------
# 4. Create the count matrix
# ------------------------------------------------------------

counts <- as.matrix(fc[, 7:ncol(fc)])
storage.mode(counts) <- "numeric"
rownames(counts) <- fc$Geneid

# Removed the HPC paths and BAM suffixes from sample column names

sample_names <- colnames(counts)
sample_names <- sub("^.*/", "", sample_names)
sample_names <- sub("^.*\\\\", "", sample_names)
sample_names <- sub("_Aligned.sortedByCoord.out.bam$", "", sample_names)

colnames(counts) <- sample_names

cat("Count-matrix dimensions:", nrow(counts), "genes and", ncol(counts), "samples\n")
print(colnames(counts))

# Confirmed 24 samples

if (ncol(counts) != 24) {
  stop(paste("Expected 24 samples but found", ncol(counts)))
}

# Checked count values

if (anyNA(counts)) {
  stop("Missing values were detected in the count matrix")
}

if (any(counts < 0)) {
  stop("Negative values were detected in the count matrix")
}

# Checked duplicated gene IDs

duplicated_gene_number <- sum(duplicated(rownames(counts)))

cat("Duplicated gene IDs:", duplicated_gene_number, "\n")

if (anyDuplicated(rownames(counts))) {
  stop("Duplicated gene IDs were found in the featureCounts file")
}

# ------------------------------------------------------------
# 5. Import and validate sample metadata
# ------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_metadata_columns <- c(
  "sample_id",
  "donor",
  "stimulus",
  "time"
)

if (!all(required_metadata_columns %in% colnames(metadata))) {
  stop("Metadata must contain: sample_id, donor, stimulus and time")
}

# Keeping only required columns

metadata <- metadata[, required_metadata_columns]

# Removing accidental spaces

metadata$sample_id <- trimws(metadata$sample_id)
metadata$donor <- trimws(metadata$donor)
metadata$stimulus <- trimws(metadata$stimulus)
metadata$time <- trimws(as.character(metadata$time))

# Checking for blank metadata entries

if (any(metadata$sample_id == "")) {
  stop("Blank sample IDs were found")
}

if (any(metadata$donor == "")) {
  stop("Blank donor labels were found")
}

if (any(metadata$stimulus == "")) {
  stop("Blank stimulus labels were found")
}

if (any(metadata$time == "")) {
  stop("Blank time labels were found")
}

# Confirmed 24 metadata rows

if (nrow(metadata) != 24) {
  stop(paste("Expected 24 metadata rows but found", nrow(metadata)))
}

# Checked duplicated sample IDs

if (anyDuplicated(metadata$sample_id)) {
  stop("Duplicated sample IDs were found in sample_metadata.csv")
}

# ------------------------------------------------------------
# 6. Match metadata to count-matrix columns
# ------------------------------------------------------------

missing_from_metadata <- setdiff(
  colnames(counts),
  metadata$sample_id
)

extra_in_metadata <- setdiff(
  metadata$sample_id,
  colnames(counts)
)

if (length(missing_from_metadata) > 0) {
  stop(
    paste(
      "Samples missing from metadata:",
      paste(missing_from_metadata, collapse=", ")
    )
  )
}

if (length(extra_in_metadata) > 0) {
  stop(
    paste(
      "Extra samples in metadata:",
      paste(extra_in_metadata, collapse=", ")
    )
  )
}

# Reordering metadata once to match the count matrix

sample_index <- match(
  colnames(counts),
  metadata$sample_id
)

if (anyNA(sample_index)) {
  stop("Some count-matrix samples could not be matched to metadata")
}

metadata <- metadata[sample_index, , drop=FALSE]

if (!identical(metadata$sample_id, colnames(counts))) {
  stop("Metadata and count-matrix sample order do not match")
}

cat("Metadata and count-matrix sample order: MATCHED\n")

# ------------------------------------------------------------
# 7. Validating the experimental design
# ------------------------------------------------------------

# We Expect:
# 3 donors
# 4 stimuli
# 2 timepoints
# 8 stimulus-time groups
# 3 samples in each group

if (length(unique(metadata$donor)) != 3) {
  stop("Expected exactly three donors")
}

if (length(unique(metadata$stimulus)) != 4) {
  stop("Expected exactly four stimuli")
}

if (length(unique(metadata$time)) != 2) {
  stop("Expected exactly two timepoints")
}

metadata$group <- interaction(
  metadata$stimulus,
  metadata$time,
  drop=TRUE,
  sep="_"
)

if (nlevels(metadata$group) != 8) {
  stop("Expected eight stimulus-time groups")
}

cat("\nSamples per stimulus-time group:\n")
print(table(metadata$group))

cat("\nDonors within each stimulus-time group:\n")
print(table(metadata$group, metadata$donor))

if (!all(table(metadata$group) == 3)) {
  stop("Every stimulus-time group must contain exactly three samples")
}

if (!all(table(metadata$group, metadata$donor) == 1)) {
  stop("Every stimulus-time group must contain one sample from each donor")
}

# Saved the correctly ordered metadata

write.csv(
  metadata,
  file.path(output_dir, "ordered_sample_metadata.csv"),
  row.names=FALSE
)

# ------------------------------------------------------------
# 8. Import candidate novel-lncRNA gene IDs
# ------------------------------------------------------------

candidate_ids <- scan(
  candidate_file,
  what="character",
  quiet=TRUE
)

candidate_ids <- trimws(candidate_ids)
candidate_ids <- candidate_ids[candidate_ids != ""]
candidate_ids <- unique(candidate_ids)

cat("\nCandidate gene IDs supplied:", length(candidate_ids), "\n")
print(head(candidate_ids))

if (length(candidate_ids) == 0) {
  stop("candidate_gene_ids.txt is empty")
}

# ------------------------------------------------------------
# 9. Matched candidate IDs to the count matrix
# ------------------------------------------------------------

candidate_ids_found <- intersect(
  candidate_ids,
  rownames(counts)
)

candidate_ids_missing <- setdiff(
  candidate_ids,
  rownames(counts)
)

cat("Candidate loci found in count matrix:", length(candidate_ids_found), "\n")
cat("Candidate IDs missing from count matrix:", length(candidate_ids_missing), "\n")

if (length(candidate_ids_found) == 0) {
  stop("None of the candidate gene IDs were found in the count matrix")
}

write.table(
  candidate_ids_found,
  file.path(output_dir, "candidate_gene_ids_found_in_counts.txt"),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)

write.table(
  candidate_ids_missing,
  file.path(output_dir, "candidate_gene_ids_missing_from_counts.txt"),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)

# ------------------------------------------------------------
# 10. Confirming that no novel candidate remains on chrM
# ------------------------------------------------------------

gene_chr <- setNames(fc$Chr, fc$Geneid)

candidate_chrM <- candidate_ids_found[
  grepl(
    "(^|;)chrM(;|$)",
    gene_chr[candidate_ids_found]
  )
]

cat("Novel candidate loci found on chrM:", length(candidate_chrM), "\n")

write.table(
  candidate_chrM,
  file.path(output_dir, "candidate_chrM_gene_ids.txt"),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)

if (length(candidate_chrM) > 0) {
  warning("One or more candidate MSTRG loci remain on chrM")
}

# ------------------------------------------------------------
# 11. TMM normalization using the complete count matrix
# ------------------------------------------------------------

dge_all <- DGEList(counts=counts)

dge_all <- calcNormFactors(
  dge_all,
  method="TMM"
)

cat("\nTMM normalization factors:\n")
print(dge_all$samples)

all_cpm <- cpm(
  dge_all,
  normalized.lib.sizes=TRUE
)

write.csv(
  all_cpm,
  file.path(output_dir, "all_gene_normalized_cpm.csv")
)

# ------------------------------------------------------------
# 12. Extracting candidate novel-lncRNA CPM values
# ------------------------------------------------------------

novel_cpm <- all_cpm[
  candidate_ids_found,
  ,
  drop=FALSE
]

cat(
  "\nNovel CPM matrix:",
  nrow(novel_cpm),
  "candidate loci and",
  ncol(novel_cpm),
  "samples\n"
)

write.csv(
  novel_cpm,
  file.path(output_dir, "novel_candidate_normalized_cpm.csv")
)

# ------------------------------------------------------------
# 13. Count supporting donors within each group
# ------------------------------------------------------------

group_levels <- levels(metadata$group)

donor_support_matrix <- sapply(
  group_levels,
  function(current_group) {
    group_samples <- metadata$group == current_group
    
    apply(
      novel_cpm[, group_samples, drop=FALSE],
      1,
      function(x) {
        passing_donors <- metadata$donor[group_samples][
          x >= cpm_threshold
        ]
        
        length(unique(passing_donors))
      }
    )
  }
)

donor_support_matrix <- as.matrix(donor_support_matrix)
rownames(donor_support_matrix) <- rownames(novel_cpm)
colnames(donor_support_matrix) <- group_levels

write.csv(
  donor_support_matrix,
  file.path(output_dir, "candidate_donor_support_by_group.csv")
)

# ------------------------------------------------------------
# 14. Applying the expression-reproducibility filter
# ------------------------------------------------------------

maximum_donor_support <- apply(
  donor_support_matrix,
  1,
  max
)

expression_pass <- maximum_donor_support >= minimum_donors

expression_pass_ids <- rownames(novel_cpm)[expression_pass]
expression_failed_ids <- rownames(novel_cpm)[!expression_pass]

cat("\nExpression filtering results:\n")
print(table(expression_pass))

cat("Expression-pass loci:", length(expression_pass_ids), "\n")
cat("Expression-failed loci:", length(expression_failed_ids), "\n")

# ------------------------------------------------------------
# 15. Identify the best-supported group for each candidate
# ------------------------------------------------------------

best_supported_group <- apply(
  donor_support_matrix,
  1,
  function(x) {
    paste(names(x)[x == max(x)], collapse=";")
  }
)

candidate_screening_results <- data.frame(
  gene_id=rownames(novel_cpm),
  maximum_CPM=apply(novel_cpm, 1, max),
  mean_CPM=apply(novel_cpm, 1, mean),
  maximum_donor_support=maximum_donor_support,
  best_supported_group=best_supported_group,
  expression_pass=expression_pass,
  stringsAsFactors=FALSE
)

candidate_screening_results <- candidate_screening_results[
  order(
    candidate_screening_results$expression_pass,
    candidate_screening_results$maximum_donor_support,
    candidate_screening_results$maximum_CPM,
    decreasing=TRUE
  ),
]

write.csv(
  candidate_screening_results,
  file.path(output_dir, "candidate_expression_screening_results.csv"),
  row.names=FALSE
)

# ------------------------------------------------------------
# 16. Save expression-pass and expression-failed IDs
# ------------------------------------------------------------

write.table(
  expression_pass_ids,
  file.path(output_dir, "expression_pass_gene_ids.txt"),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)

write.table(
  expression_failed_ids,
  file.path(output_dir, "expression_failed_gene_ids.txt"),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)

# ------------------------------------------------------------
# 17. Create expression-filter summary
# ------------------------------------------------------------

expression_summary <- data.frame(
  category=c(
    "Candidate IDs supplied",
    "Candidate loci found in counts",
    "Candidate IDs missing from counts",
    "Candidate loci on chrM",
    "Expression-pass loci",
    "Expression-failed loci"
  ),
  number=c(
    length(candidate_ids),
    length(candidate_ids_found),
    length(candidate_ids_missing),
    length(candidate_chrM),
    length(expression_pass_ids),
    length(expression_failed_ids)
  )
)

print(expression_summary)

write.csv(
  expression_summary,
  file.path(output_dir, "expression_filter_summary.csv"),
  row.names=FALSE
)

# ------------------------------------------------------------
# 18. Prepared expressed non-mitochondrial genes for MDS QC
# ------------------------------------------------------------

is_mitochondrial <- grepl(
  "(^|;)chrM(;|$)",
  fc$Chr
)

keep_for_mds <- filterByExpr(
  dge_all,
  group=metadata$group
)

keep_for_mds <- keep_for_mds & !is_mitochondrial

dge_mds <- dge_all[
  keep_for_mds,
  ,
  keep.lib.sizes=FALSE
]

dge_mds <- calcNormFactors(
  dge_mds,
  method="TMM"
)

cat(
  "\nGenes retained for non-mitochondrial MDS:",
  nrow(dge_mds),
  "\n"
)

# ------------------------------------------------------------
# 19. Generate sample-level MDS plot
# ------------------------------------------------------------

group_palette <- setNames(
  hcl.colors(
    nlevels(metadata$group),
    palette="Dark 3"
  ),
  levels(metadata$group)
)

sample_colors <- group_palette[
  as.character(metadata$group)
]

pdf(
  file.path(output_dir, "sample_MDS_plot.pdf"),
  width=11,
  height=8
)

plotMDS(
  dge_mds,
  col=sample_colors,
  pch=19,
  labels=metadata$sample_id,
  main="MDS plot of 24 RNA-seq samples"
)

legend(
  "topright",
  legend=names(group_palette),
  col=group_palette,
  pch=19,
  cex=0.75
)

dev.off()

plotMDS(dge_mds, col=sample_colors, pch=19, labels=metadata$sample_id, main="MDS plot of 24 RNA-seq samples")


# ------------------------------------------------------------
# 20. Printing final results
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Novel-lncRNA expression screening complete\n")
cat("========================================\n")
cat("Candidate IDs supplied:", length(candidate_ids), "\n")
cat("Candidate loci found:", length(candidate_ids_found), "\n")
cat("Candidate IDs missing:", length(candidate_ids_missing), "\n")
cat("Novel candidate loci on chrM:", length(candidate_chrM), "\n")
cat("Expression-pass loci:", length(expression_pass_ids), "\n")
cat("Expression-failed loci:", length(expression_failed_ids), "\n")
cat("CPM threshold:", cpm_threshold, "\n")
cat("Required supporting donors:", minimum_donors, "\n")
cat("Results saved in:", output_dir, "\n")


# ------------------------------------------------------------
# 21. Expression-screening bar plot
# ------------------------------------------------------------


if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
library(ggplot2)

screening_bar_data <- data.frame(
  status=factor(
    c("Expression pass", "Expression failed"),
    levels=c("Expression pass", "Expression failed")
  ),
  loci=c(
    length(expression_pass_ids),
    length(expression_failed_ids)
  )
)

screening_bar_data$percentage <- 100 * screening_bar_data$loci / sum(screening_bar_data$loci)

screening_bar_data$label <- paste0(
  screening_bar_data$loci,
  " (",
  sprintf("%.1f", screening_bar_data$percentage),
  "%)"
)

expression_bar_plot <- ggplot(
  screening_bar_data,
  aes(x=status, y=loci, fill=status)
) +
  geom_col(width=0.65) +
  geom_text(
    aes(label=label),
    vjust=-0.5,
    size=5
  ) +
  scale_fill_manual(
    values=c(
      "Expression pass"="#2E8B57",
      "Expression failed"="#D95F5F"
    )
  ) +
  scale_y_continuous(
    limits=c(0, max(screening_bar_data$loci) * 1.12),
    expand=c(0, 0)
  ) +
  labs(
    title="Novel-lncRNA expression screening",
    subtitle="CPM ≥0.5 in at least two donors within the same stimulus–time group",
    x=NULL,
    y="Number of candidate loci"
  ) +
  theme_classic(base_size=14) +
  theme(
    legend.position="none",
    plot.title=element_text(face="bold"),
    axis.text.x=element_text(face="bold")
  )

expression_bar_plot

ggsave(file.path(output_dir, "expression_screen_pass_fail_bar.pdf"), expression_bar_plot, width=8, height=6)
ggsave(file.path(output_dir, "expression_screen_pass_fail_bar.png"), expression_bar_plot, width=8, height=6, dpi=300)




# ------------------------------------------------------------
# 21. Donor-support heatmap
# ------------------------------------------------------------


if (!requireNamespace("pheatmap", quietly=TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# ------------------------------------------------------------
# 1. Create pass/failed row annotation
# ------------------------------------------------------------

heatmap_status <- data.frame(
  Status=factor(
    ifelse(expression_pass, "Pass", "Failed"),
    levels=c("Pass", "Failed")
  )
)

rownames(heatmap_status) <- rownames(donor_support_matrix)

# ------------------------------------------------------------
# 2. Ordering candidates
# Pass candidates first, followed by failed candidates
# ------------------------------------------------------------

heatmap_order <- order(
  expression_pass,
  maximum_donor_support,
  decreasing=TRUE
)

ordered_support <- donor_support_matrix[
  heatmap_order,
  ,
  drop=FALSE
]

ordered_annotation <- heatmap_status[
  heatmap_order,
  ,
  drop=FALSE
]

# ------------------------------------------------------------
# 3. Define pass/failed annotation colors
# ------------------------------------------------------------

annotation_colors <- list(
  Status=c(
    "Pass"="#2E8B57",
    "Failed"="#D95F5F"
  )
)

# ------------------------------------------------------------
# 4. Display the heatmap 
# ------------------------------------------------------------


if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("ComplexHeatmap", quietly=TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

library(ComplexHeatmap)
library(grid)

# Convert donor-support values to discrete categories

support_matrix <- matrix(
  as.character(ordered_support),
  nrow=nrow(ordered_support),
  ncol=ncol(ordered_support),
  dimnames=dimnames(ordered_support)
)

# Define colors

support_colors <- c(
  "0"="#F2F2F2",
  "1"="#56B4E9",
  "2"="#E69F00",
  "3"="#6A3D9A"
)

status_colors <- c(
  "Pass"="#2E8B57",
  "Failed"="#D95F5F"
)

status_vector <- ordered_annotation$Status

# Row annotation without its automatic legend

status_annotation <- rowAnnotation(
  Status=status_vector,
  col=list(Status=status_colors),
  show_annotation_name=FALSE,
  show_legend=FALSE,
  width=unit(4, "mm")
)

# Create heatmap without its automatic legend

support_ht <- Heatmap(
  support_matrix,
  name="Donor support",
  col=support_colors,
  
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  
  show_row_names=FALSE,
  show_column_names=TRUE,
  column_names_rot=90,
  column_names_gp=gpar(fontsize=9),
  
  row_split=status_vector,
  cluster_row_slices=FALSE,
  row_gap=unit(1.5, "mm"),
  row_title=NULL,
  
  show_heatmap_legend=FALSE,
  use_raster=TRUE,
  
  column_title=paste0(
    "Donor support for candidate novel-lncRNA loci\n",
    "Pass: ",
    length(expression_pass_ids),
    " | Failed: ",
    length(expression_failed_ids)
  )
)

# Donor-support legend

support_legend <- Legend(
  title="Donor support",
  labels=c(
    "3 donors",
    "2 donors",
    "1 donor",
    "0 donors"
  ),
  legend_gp=gpar(
    fill=c(
      support_colors["3"],
      support_colors["2"],
      support_colors["1"],
      support_colors["0"]
    )
  )
)

# Status legend

status_legend <- Legend(
  title="Status",
  labels=c("Pass", "Failed"),
  legend_gp=gpar(
    fill=c(
      status_colors["Pass"],
      status_colors["Failed"]
    )
  )
)

# Drawing Donor support plot

draw(
  status_annotation + support_ht,
  heatmap_legend_side="right",
  heatmap_legend_list=list(
    support_legend,
    status_legend
  ),
  merge_legends=TRUE
)
 






