##################################################
#Importing counts, metadata and gene annotation

base_dir <- "D:/a/final_analysis"
output_dir <- file.path(base_dir, "results")
dir.create(output_dir, showWarnings=FALSE)

# Importing featureCounts output
fc <- read.delim(
  file.path(base_dir, "final_gene_counts.txt"),
  comment.char="#",
  check.names=FALSE
)

counts <- as.matrix(fc[, 7:ncol(fc)])
storage.mode(counts) <- "numeric"
rownames(counts) <- fc$Geneid

# Cleaning sample names
sample_names <- gsub("\\\\", "/", colnames(counts))
sample_names <- basename(sample_names)
sample_names <- sub("_Aligned.sortedByCoord.out.bam$", "", sample_names)
colnames(counts) <- sample_names

# Keeping featureCounts gene information
gene_info <- fc[, 1:6]

# Importing metadata
metadata <- read.csv(
  file.path(base_dir, "sample_metadata.csv"),
  stringsAsFactors=FALSE
)

metadata$sample_id <- trimws(metadata$sample_id)
metadata$donor <- trimws(metadata$donor)
metadata$stimulus <- trimws(metadata$stimulus)
metadata$time <- trimws(as.character(metadata$time))

# Arrange metadata in count-matrix order
metadata <- metadata[
  match(colnames(counts), metadata$sample_id),
  ,
  drop=FALSE
]

metadata$donor <- factor(metadata$donor)
metadata$stimulus <- factor(metadata$stimulus)
metadata$time <- factor(metadata$time)

# Importing gene annotation
annotation <- read.delim(
  file.path(base_dir, "final_gene_annotation.tsv"),
  stringsAsFactors=FALSE
)

annotation$gene_name[
  is.na(annotation$gene_name) | annotation$gene_name==""
] <- annotation$gene_id[
  is.na(annotation$gene_name) | annotation$gene_name==""
]

annotation <- annotation[
  match(rownames(counts), annotation$gene_id),
  ,
  drop=FALSE
]

# Classify genes
known_lnc_types <- c(
  "lncRNA",
  "lincRNA",
  "antisense",
  "processed_transcript",
  "sense_intronic",
  "sense_overlapping"
)

annotation$analysis_class <- "other"

annotation$analysis_class[
  annotation$gene_type=="protein_coding"
] <- "protein_coding"

annotation$analysis_class[
  annotation$gene_type %in% known_lnc_types
] <- "known_lncRNA"

annotation$analysis_class[
  annotation$gene_type=="novel_lncRNA" |
    grepl("^MSTRG\\.", annotation$gene_id)
] <- "novel_lncRNA"

# One essential verification
stopifnot(
  ncol(counts)==24,
  !anyDuplicated(rownames(counts)),
  identical(metadata$sample_id, colnames(counts)),
  !anyNA(annotation$gene_id),
  sum(annotation$analysis_class=="novel_lncRNA")==607
)

# Saved imported dataset
saveRDS(
  list(
    counts=counts,
    metadata=metadata,
    annotation=annotation,
    gene_info=gene_info
  ),
  file.path(output_dir, "imported_final_dataset.rds")
)

# Important summary
dim(counts)
table(annotation$analysis_class)
table(metadata$stimulus, metadata$time)



##############################################################
#Generatng Sequencing-library bar plot


library(edgeR)

# Create stimulus-time groups
metadata$group <- interaction(
  metadata$stimulus,
  metadata$time,
  sep="_",
  drop=TRUE
)

# Create edgeR object
dge <- DGEList(counts=counts)

# Remove genes with insufficient expression for DE testing
keep <- filterByExpr(
  dge,
  group=metadata$group
)

table(keep)

# Apply filtering and TMM normalization
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# Annotation corresponding to retained genes
annotation_filtered <- annotation[keep, , drop=FALSE]

# Number of retained genes by class
table(annotation_filtered$analysis_class)

# Sequencing-library sizes
barplot(
  dge$samples$lib.size / 1000000,
  names.arg=metadata$sample_id,
  las=2,
  cex.names=0.6,
  ylab="Library size (million reads)",
  main="RNA-seq library sizes"
)

#############################################################

#Generating the PCA

library(edgeR)
library(ggplot2)

# Log-normalized expression for PCA
logCPM <- cpm(dge, log=TRUE, prior.count=2)

# Global PCA using all 24 samples
pca <- prcomp(t(logCPM), scale.=FALSE)

# Percentage of variance explained
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Set panel order only
metadata$time <- factor(
  as.character(metadata$time),
  levels=c("4h", "24h")
)

# Combine PCA coordinates with metadata
pca_data <- data.frame(
  PC1=pca$x[, 1],
  PC2=pca$x[, 2],
  metadata
)

# Create boundaries around the three donors
hulls <- do.call(
  rbind,
  lapply(
    split(
      pca_data,
      interaction(pca_data$time, pca_data$stimulus)
    ),
    function(x) {
      x[chull(x$PC1, x$PC2), , drop=FALSE]
    }
  )
)

# Create the PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x=PC1,
    y=PC2,
    color=stimulus,
    shape=donor
  )
) +
  geom_polygon(
    data=hulls,
    aes(
      x=PC1,
      y=PC2,
      group=interaction(time, stimulus),
      fill=stimulus
    ),
    inherit.aes=FALSE,
    alpha=0.12,
    color=NA
  ) +
  geom_point(size=3.5) +
  geom_text(
    aes(label=donor),
    vjust=-0.8,
    show.legend=FALSE,
    size=3.5
  ) +
  scale_shape_manual(
    values=c(
      "A"=16,
      "B"=17,
      "C"=15
    )
  ) +
  facet_wrap(
    ~time,
    nrow=1
  ) +
  labs(
    title="PCA of 24 RNA-seq samples",
    subtitle="Color = stimulus; shape and label = donor",
    x=paste0(
      "PC1 (",
      round(variance[1], 1),
      "% variance)"
    ),
    y=paste0(
      "PC2 (",
      round(variance[2], 1),
      "% variance)"
    ),
    color="Stimulus",
    fill="Stimulus",
    shape="Donor"
  ) +
  theme_bw() +
  theme(
    legend.position="right",
    plot.title=element_text(hjust=0.5),
    plot.subtitle=element_text(hjust=0.5)
  )

pca_plot


##############################################################
#Creating the donor-adjusted edgeR model

# Create factors
metadata$donor <- factor(metadata$donor)
metadata$group <- factor(metadata$group)

# Donor-adjusted design matrix
design <- model.matrix(
  ~ 0 + group + donor,
  data=metadata
)

colnames(design) <- make.names(colnames(design))

# Essential check
stopifnot(qr(design)$rank == ncol(design))

# View coefficient names
colnames(design)


#estimating the dispersion and fit the edgeR quasi-likelihood model

dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmQLFit(dge, design, robust=TRUE)
plotQLDisp(fit)

#defining the six stimulus-versus-RPMI contrasts

library(limma)

contrasts <- makeContrasts(
  KTCbg_4h   = groupcaur_KTCbg_4h   - groupRPMI_4h,
  KTClive_4h = groupcaur_KTClive_4h - groupRPMI_4h,
  KTCman_4h  = groupcaur_KTCman_4h  - groupRPMI_4h,
  KTCbg_24h   = groupcaur_KTCbg_24h   - groupRPMI_24h,
  KTClive_24h = groupcaur_KTClive_24h - groupRPMI_24h,
  KTCman_24h  = groupcaur_KTCman_24h  - groupRPMI_24h,
  levels=design
)

colnames(contrasts)




##############################################################
# Test all six contrasts and save annotated DEG results

de_dir <- file.path(output_dir, "DE_results")
dir.create(de_dir, showWarnings=FALSE, recursive=TRUE)

# Confirm that everything remains aligned

stopifnot(
  identical(colnames(dge), metadata$sample_id),
  identical(rownames(dge), annotation_filtered$gene_id),
  identical(rownames(contrasts), colnames(design))
)

# Objects for storing results

qlf_tests <- setNames(
  vector("list", ncol(contrasts)),
  colnames(contrasts)
)

de_results <- setNames(
  vector("list", ncol(contrasts)),
  colnames(contrasts)
)

# Testing each contrast

for (contrast_name in colnames(contrasts)) {
  
# Quasi-likelihood F-test
  qlf <- glmQLFTest(
    fit,
    contrast=contrasts[, contrast_name]
  )
  
  qlf_tests[[contrast_name]] <- qlf
  
# Extract all tested genes
  result_table <- topTags(
    qlf,
    n=Inf,
    sort.by="PValue"
  )$table
  
  gene_ids <- rownames(result_table)
  
# Match gene annotation
  annotation_index <- match(
    gene_ids,
    annotation_filtered$gene_id
  )
  
  stopifnot(!anyNA(annotation_index))
  
# Combine statistics and annotation
  result_table <- data.frame(
    contrast=contrast_name,
    gene_id=gene_ids,
    gene_name=annotation_filtered$gene_name[annotation_index],
    gene_type=annotation_filtered$gene_type[annotation_index],
    analysis_class=annotation_filtered$analysis_class[annotation_index],
    result_table,
    row.names=NULL,
    check.names=FALSE
  )
  
# Classify DE genes
  result_table$DE_status <- "Not_significant"
  
  result_table$DE_status[
    result_table$FDR < 0.05 &
      result_table$logFC > 1
  ] <- "Up"
  
  result_table$DE_status[
    result_table$FDR < 0.05 &
      result_table$logFC < -1
  ] <- "Down"
  
  de_results[[contrast_name]] <- result_table
  
# Save all tested genes
  write.csv(
    result_table,
    file.path(
      de_dir,
      paste0(contrast_name, "_all_genes.csv")
    ),
    row.names=FALSE
  )
  
# Save only significant DE genes
  significant_table <- result_table[
    result_table$DE_status != "Not_significant",
    ,
    drop=FALSE
  ]
  
  write.csv(
    significant_table,
    file.path(
      de_dir,
      paste0(contrast_name, "_significant_DEGs.csv")
    ),
    row.names=FALSE
  )
}



##############################################################
# Summarize DE results

de_summary <- do.call(
  rbind,
  lapply(names(de_results), function(contrast_name) {
    
    result_table <- de_results[[contrast_name]]
    
    significant_table <- result_table[
      result_table$DE_status != "Not_significant",
      ,
      drop=FALSE
    ]
    
    data.frame(
      contrast=contrast_name,
      genes_tested=nrow(result_table),
      upregulated=sum(result_table$DE_status == "Up"),
      downregulated=sum(result_table$DE_status == "Down"),
      total_DE=nrow(significant_table),
      
      protein_coding_DE=sum(
        significant_table$analysis_class == "protein_coding"
      ),
      
      known_lncRNA_DE=sum(
        significant_table$analysis_class == "known_lncRNA"
      ),
      
      novel_lncRNA_DE=sum(
        significant_table$analysis_class == "novel_lncRNA"
      )
    )
  })
)

rownames(de_summary) <- NULL

de_summary

write.csv(
  de_summary,
  file.path(de_dir, "DE_summary_all_contrasts.csv"),
  row.names=FALSE
)

saveRDS(
  list(
    dge=dge,
    design=design,
    contrasts=contrasts,
    fit=fit,
    qlf_tests=qlf_tests,
    de_results=de_results,
    de_summary=de_summary
  ),
  file.path(output_dir, "edgeR_DE_analysis.rds")
)



##############################################################
# MD plots for all six contrasts

md_dir <- file.path(output_dir, "MD_plots")
dir.create(md_dir, showWarnings=FALSE, recursive=TRUE)

png(
  file.path(md_dir, "MD_plots_all_contrasts.png"),
  width=2600,
  height=1700,
  res=220
)

par(
  mfrow=c(2, 3),
  mar=c(4.5, 4.5, 3.5, 1)
)

for (contrast_name in names(de_results)) {
  
  result_table <- de_results[[contrast_name]]
  
  point_colors <- rep(
    adjustcolor("grey60", alpha.f=0.45),
    nrow(result_table)
  )
  
  point_colors[
    result_table$DE_status == "Up"
  ] <- adjustcolor("#D73027", alpha.f=0.75)
  
  point_colors[
    result_table$DE_status == "Down"
  ] <- adjustcolor("#2166AC", alpha.f=0.75)
  
  plot(
    result_table$logCPM,
    result_table$logFC,
    pch=16,
    cex=0.45,
    col=point_colors,
    xlab="Average logCPM",
    ylab="log2 fold change",
    main=contrast_name
  )
  
  abline(
    h=c(-1, 0, 1),
    lty=c(2, 1, 2),
    col=c("grey40", "black", "grey40")
  )
  
  legend(
    "topright",
    legend=c("Up", "Down", "Not significant"),
    col=c("#D73027", "#2166AC", "grey60"),
    pch=16,
    cex=0.7,
    bty="n"
  )
}

dev.off()




##############################################################
# Six-panel volcano plot

library(ggplot2)

volcano_dir <- file.path(output_dir, "volcano_plots")
dir.create(volcano_dir, showWarnings=FALSE, recursive=TRUE)

# Combine the six result tables

volcano_data <- do.call(
  rbind,
  de_results
)

rownames(volcano_data) <- NULL

# Set contrast order

contrast_order <- c(
  "KTCbg_4h",
  "KTClive_4h",
  "KTCman_4h",
  "KTCbg_24h",
  "KTClive_24h",
  "KTCman_24h"
)

volcano_data$contrast <- factor(
  volcano_data$contrast,
  levels=contrast_order
)

volcano_data$DE_status <- factor(
  volcano_data$DE_status,
  levels=c("Down", "Not_significant", "Up")
)

# Prevent infinite values if an FDR is numerically zero

volcano_data$minus_log10_FDR <- -log10(
  pmax(volcano_data$FDR, .Machine$double.xmin)
)

# Plot

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x=logFC,
    y=minus_log10_FDR
  )
) +
  geom_point(
    data=subset(
      volcano_data,
      DE_status == "Not_significant"
    ),
    color="grey70",
    alpha=0.45,
    size=0.65
  ) +
  geom_point(
    data=subset(
      volcano_data,
      DE_status != "Not_significant"
    ),
    aes(color=DE_status),
    alpha=0.80,
    size=0.85
  ) +
  geom_vline(
    xintercept=c(-1, 1),
    linetype="dashed",
    color="grey35"
  ) +
  geom_hline(
    yintercept=-log10(0.05),
    linetype="dashed",
    color="grey35"
  ) +
  scale_color_manual(
    values=c(
      "Down"="#2166AC",
      "Up"="#D73027"
    ),
    labels=c(
      "Downregulated",
      "Upregulated"
    )
  ) +
  facet_wrap(
    ~contrast,
    ncol=3,
    scales="free_y"
  ) +
  labs(
    title="Differential expression across C. auris stimuli",
    subtitle="Donor-adjusted edgeR quasi-likelihood analysis",
    x=expression(log[2]~fold~change),
    y=expression(-log[10]~FDR),
    color="DE status"
  ) +
  theme_bw(base_size=12) +
  theme(
    strip.background=element_rect(
      fill="grey95",
      color="grey30"
    ),
    strip.text=element_text(
      face="bold",
      size=11
    ),
    plot.title=element_text(
      face="bold",
      hjust=0.5
    ),
    plot.subtitle=element_text(hjust=0.5),
    legend.position="right",
    panel.grid.minor=element_blank()
  )

volcano_plot

ggsave(
  file.path(volcano_dir, "volcano_all_contrasts.png"),
  volcano_plot,
  width=13,
  height=8,
  dpi=600,
  bg="white"
)

ggsave(
  file.path(volcano_dir, "volcano_all_contrasts.pdf"),
  volcano_plot,
  width=13,
  height=8
)




##############################################################
# Volcano plots colored by gene class

# Assigning readable class names

volcano_data$gene_class_plot <- "Other"

volcano_data$gene_class_plot[
  volcano_data$analysis_class == "protein_coding"
] <- "Protein-coding"

volcano_data$gene_class_plot[
  volcano_data$analysis_class == "known_lncRNA"
] <- "Known lncRNA"

volcano_data$gene_class_plot[
  volcano_data$analysis_class == "novel_lncRNA"
] <- "Novel lncRNA"

# Separate significant and non-significant genes

nonsignificant_data <- subset(
  volcano_data,
  DE_status == "Not_significant"
)

significant_data <- subset(
  volcano_data,
  DE_status != "Not_significant"
)

significant_data$gene_class_plot <- factor(
  significant_data$gene_class_plot,
  levels=c(
    "Protein-coding",
    "Known lncRNA",
    "Novel lncRNA",
    "Other"
  )
)

# Create plot

class_volcano_plot <- ggplot(
  volcano_data,
  aes(
    x=logFC,
    y=minus_log10_FDR
  )
) +
 
# Nonsignificant genes
  geom_point(
    data=nonsignificant_data,
    color="grey65",
    alpha=0.20,
    size=0.60
  ) +
  
# Protein-coding and other genes
  geom_point(
    data=subset(
      significant_data,
      gene_class_plot %in% c("Protein-coding", "Other")
    ),
    aes(color=gene_class_plot),
    alpha=0.40,
    size=0.80
  ) +
  
# Known lncRNAs
  geom_point(
    data=subset(
      significant_data,
      gene_class_plot == "Known lncRNA"
    ),
    aes(color=gene_class_plot),
    alpha=0.50,
    size=1.10
  ) +
  
# Novel lncRNAs
  geom_point(
    data=subset(
      significant_data,
      gene_class_plot == "Novel lncRNA"
    ),
    aes(color=gene_class_plot),
    alpha=0.60,
    size=1.50
  ) +
  
  geom_vline(
    xintercept=c(-1, 1),
    linetype="dashed",
    color="grey35"
  ) +
  geom_hline(
    yintercept=-log10(0.05),
    linetype="dashed",
    color="grey35"
  ) +
  
  scale_color_manual(
    values=c(
      "Protein-coding"="#E69F00",
      "Known lncRNA"="#0072B2",
      "Novel lncRNA"="#CC79A7",
      "Other"="#000000"
    ),
    limits=c(
      "Protein-coding",
      "Known lncRNA",
      "Novel lncRNA",
      "Other"
    ),
    drop=FALSE
  ) +
  
  facet_wrap(
    ~contrast,
    ncol=3,
    scales="free_y"
  ) +
  
  labs(
    title="Differential expression across C. auris stimuli",
    subtitle=paste(
      "Significant genes: FDR < 0.05 and |log2FC| > 1;",
      "right = upregulated, left = downregulated"
    ),
    x=expression(log[2]~fold~change),
    y=expression(-log[10]~FDR),
    color="Gene class"
  ) +
  
  theme_bw(base_size=12) +
  theme(
    strip.background=element_rect(
      fill="grey95",
      color="grey30"
    ),
    strip.text=element_text(
      face="bold",
      size=11
    ),
    plot.title=element_text(
      face="bold",
      hjust=0.5
    ),
    plot.subtitle=element_text(hjust=0.5),
    legend.position="right",
    panel.grid.minor=element_blank()
  ) +
  
  guides(
    color=guide_legend(
      override.aes=list(
        size=3,
        alpha=1
      )
    )
  )

class_volcano_plot


ggsave(
  file.path(volcano_dir, "volcano_by_gene_class.png"),
  class_volcano_plot,
  width=13,
  height=8,
  dpi=600,
  bg="white"
)

ggsave(
  file.path(volcano_dir, "volcano_by_gene_class.pdf"),
  class_volcano_plot,
  width=13,
  height=8
)


#####################################################################
#UpSet overlap analysis
#It will determine how many of the DEGs are shared or unique across the six contrasts
#this part includes all gene type in the list (coding, known and novel lncRNA and others)

 

library(UpSetR)

overlap_dir <- file.path(
  output_dir,
  "overlap_analysis"
)

dir.create(
  overlap_dir,
  showWarnings=FALSE,
  recursive=TRUE
)

##############################################################
# Extract significant gene IDs separately

KTCbg_4h_table <- de_results[["KTCbg_4h"]]

KTCbg_4h_genes <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$DE_status != "Not_significant"
]


KTClive_4h_table <- de_results[["KTClive_4h"]]

KTClive_4h_genes <- KTClive_4h_table$gene_id[
  KTClive_4h_table$DE_status != "Not_significant"
]


KTCman_4h_table <- de_results[["KTCman_4h"]]

KTCman_4h_genes <- KTCman_4h_table$gene_id[
  KTCman_4h_table$DE_status != "Not_significant"
]


KTCbg_24h_table <- de_results[["KTCbg_24h"]]

KTCbg_24h_genes <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$DE_status != "Not_significant"
]


KTClive_24h_table <- de_results[["KTClive_24h"]]

KTClive_24h_genes <- KTClive_24h_table$gene_id[
  KTClive_24h_table$DE_status != "Not_significant"
]


KTCman_24h_table <- de_results[["KTCman_24h"]]

KTCman_24h_genes <- KTCman_24h_table$gene_id[
  KTCman_24h_table$DE_status != "Not_significant"
]

#Verify that the numbers match the summary:
c(
  KTCbg_4h=length(KTCbg_4h_genes),
  KTClive_4h=length(KTClive_4h_genes),
  KTCman_4h=length(KTCman_4h_genes),
  KTCbg_24h=length(KTCbg_24h_genes),
  KTClive_24h=length(KTClive_24h_genes),
  KTCman_24h=length(KTCman_24h_genes)
)


#Create the named list

DE_gene_sets <- list(
  KTCbg_4h=KTCbg_4h_genes,
  KTClive_4h=KTClive_4h_genes,
  KTCman_4h=KTCman_4h_genes,
  KTCbg_24h=KTCbg_24h_genes,
  KTClive_24h=KTClive_24h_genes,
  KTCman_24h=KTCman_24h_genes
)

upset_data <- fromList(DE_gene_sets)


#Create and save the UpSet plot

upset(
  upset_data,
  sets=c(
    "KTCbg_4h",
    "KTClive_4h",
    "KTCman_4h",
    "KTCbg_24h",
    "KTClive_24h",
    "KTCman_24h"
  ),
  nsets=6,
  nintersects=30,
  keep.order=TRUE,
  order.by="freq",
  show.numbers="yes",
  point.size=3,
  line.size=1,
  main.bar.color="#E69F00",
  sets.bar.color="#4D4D4D",
  matrix.color="#000000",
  mainbar.y.label="Number of shared DE genes",
  sets.x.label="DE genes per contrast",
  text.scale=1.2
)

#################################################################

#separate overlaps into upregulated and downregulated genes

##############################################################
# Extract upregulated genes

KTCbg_4h_up <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$DE_status == "Up"
]

KTClive_4h_up <- KTClive_4h_table$gene_id[
  KTClive_4h_table$DE_status == "Up"
]

KTCman_4h_up <- KTCman_4h_table$gene_id[
  KTCman_4h_table$DE_status == "Up"
]

KTCbg_24h_up <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$DE_status == "Up"
]

KTClive_24h_up <- KTClive_24h_table$gene_id[
  KTClive_24h_table$DE_status == "Up"
]

KTCman_24h_up <- KTCman_24h_table$gene_id[
  KTCman_24h_table$DE_status == "Up"
]


#Create the upregulated-gene UpSet plot:


up_gene_sets <- list(
  KTCbg_4h=KTCbg_4h_up,
  KTClive_4h=KTClive_4h_up,
  KTCman_4h=KTCman_4h_up,
  KTCbg_24h=KTCbg_24h_up,
  KTClive_24h=KTClive_24h_up,
  KTCman_24h=KTCman_24h_up
)

up_upset_data <- fromList(up_gene_sets)

graphics.off()

upset(
  up_upset_data,
  sets=names(up_gene_sets),
  nsets=6,
  nintersects=30,
  keep.order=TRUE,
  order.by="freq",
  show.numbers="yes",
  main.bar.color="#D73027",
  sets.bar.color="#D73027",
  matrix.color="#000000",
  mainbar.y.label="Shared upregulated genes",
  sets.x.label="Upregulated genes per contrast",
  text.scale=1.2
)



##############################################################
# Extract downregulated genes

KTCbg_4h_down <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$DE_status == "Down"
]

KTClive_4h_down <- KTClive_4h_table$gene_id[
  KTClive_4h_table$DE_status == "Down"
]

KTCman_4h_down <- KTCman_4h_table$gene_id[
  KTCman_4h_table$DE_status == "Down"
]

KTCbg_24h_down <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$DE_status == "Down"
]

KTClive_24h_down <- KTClive_24h_table$gene_id[
  KTClive_24h_table$DE_status == "Down"
]

KTCman_24h_down <- KTCman_24h_table$gene_id[
  KTCman_24h_table$DE_status == "Down"
]


#Create the downregulated-gene UpSet plot

down_gene_sets <- list(
  KTCbg_4h=KTCbg_4h_down,
  KTClive_4h=KTClive_4h_down,
  KTCman_4h=KTCman_4h_down,
  KTCbg_24h=KTCbg_24h_down,
  KTClive_24h=KTClive_24h_down,
  KTCman_24h=KTCman_24h_down
)

down_upset_data <- fromList(down_gene_sets)

graphics.off()

upset(
  down_upset_data,
  sets=names(down_gene_sets),
  nsets=6,
  nintersects=30,
  keep.order=TRUE,
  order.by="freq",
  show.numbers="yes",
  main.bar.color="#2166AC",
  sets.bar.color="#2166AC",
  matrix.color="#000000",
  mainbar.y.label="Shared downregulated genes",
  sets.x.label="Downregulated genes per contrast",
  text.scale=1.2
)



####################################################################
#this part includes only the novel lncRNA 

##############################################################
# Significant novel lncRNAs-all novel lncRNAs

novel_KTCbg_4h <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$DE_status != "Not_significant" &
    KTCbg_4h_table$analysis_class == "novel_lncRNA"
]

novel_KTClive_4h <- KTClive_4h_table$gene_id[
  KTClive_4h_table$DE_status != "Not_significant" &
    KTClive_4h_table$analysis_class == "novel_lncRNA"
]

novel_KTCman_4h <- KTCman_4h_table$gene_id[
  KTCman_4h_table$DE_status != "Not_significant" &
    KTCman_4h_table$analysis_class == "novel_lncRNA"
]

novel_KTCbg_24h <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$DE_status != "Not_significant" &
    KTCbg_24h_table$analysis_class == "novel_lncRNA"
]

novel_KTClive_24h <- KTClive_24h_table$gene_id[
  KTClive_24h_table$DE_status != "Not_significant" &
    KTClive_24h_table$analysis_class == "novel_lncRNA"
]

novel_KTCman_24h <- KTCman_24h_table$gene_id[
  KTCman_24h_table$DE_status != "Not_significant" &
    KTCman_24h_table$analysis_class == "novel_lncRNA"
]


#Create the novel-lncRNA plot

novel_lnc_sets <- list(
  KTCbg_4h=novel_KTCbg_4h,
  KTClive_4h=novel_KTClive_4h,
  KTCman_4h=novel_KTCman_4h,
  KTCbg_24h=novel_KTCbg_24h,
  KTClive_24h=novel_KTClive_24h,
  KTCman_24h=novel_KTCman_24h
)

novel_lnc_upset_data <- fromList(novel_lnc_sets)

graphics.off()

upset(
  novel_lnc_upset_data,
  sets=names(novel_lnc_sets),
  nsets=6,
  nintersects=20,
  keep.order=TRUE,
  order.by="freq",
  show.numbers="yes",
  main.bar.color="#E69F00",
  sets.bar.color="#4D4D4D",
  matrix.color="#000000",
  mainbar.y.label="Shared novel lncRNAs",
  sets.x.label="Novel lncRNAs per contrast",
  text.scale=1.2
)



############################################################
# Extract upregulated novel lncRNAs

novel_up_KTCbg_4h <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$analysis_class == "novel_lncRNA" &
    KTCbg_4h_table$DE_status == "Up"
]

novel_up_KTClive_4h <- KTClive_4h_table$gene_id[
  KTClive_4h_table$analysis_class == "novel_lncRNA" &
  KTClive_4h_table$DE_status == "Up"
]

novel_up_KTCman_4h <- KTCman_4h_table$gene_id[
  KTCman_4h_table$analysis_class == "novel_lncRNA" &
  KTCman_4h_table$DE_status == "Up"
]

novel_up_KTCbg_24h <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$analysis_class == "novel_lncRNA" &
  KTCbg_24h_table$DE_status == "Up"
]

novel_up_KTClive_24h <- KTClive_24h_table$gene_id[
  KTClive_24h_table$analysis_class == "novel_lncRNA" &
  KTClive_24h_table$DE_status == "Up"
]

novel_up_KTCman_24h <- KTCman_24h_table$gene_id[
  KTCman_24h_table$analysis_class == "novel_lncRNA" &
  KTCman_24h_table$DE_status == "Up"
]




############################################################
# Extract downregulated novel lncRNAs

novel_down_KTCbg_4h <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$analysis_class == "novel_lncRNA" &
    KTCbg_4h_table$DE_status == "Down"
]

novel_down_KTClive_4h <- KTClive_4h_table$gene_id[
  KTClive_4h_table$analysis_class == "novel_lncRNA" &
    KTClive_4h_table$DE_status == "Down"
]

novel_down_KTCman_4h <- KTCman_4h_table$gene_id[
  KTCman_4h_table$analysis_class == "novel_lncRNA" &
    KTCman_4h_table$DE_status == "Down"
]

novel_down_KTCbg_24h <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$analysis_class == "novel_lncRNA" &
    KTCbg_24h_table$DE_status == "Down"
]

novel_down_KTClive_24h <- KTClive_24h_table$gene_id[
  KTClive_24h_table$analysis_class == "novel_lncRNA" &
    KTClive_24h_table$DE_status == "Down"
]

novel_down_KTCman_24h <- KTCman_24h_table$gene_id[
  KTCman_24h_table$analysis_class == "novel_lncRNA" &
    KTCman_24h_table$DE_status == "Down"
]




############################################################
# Upregulated novel-lncRNA overlap

library(UpSetR)

novel_up_sets <- list(
  KTCbg_4h = novel_up_KTCbg_4h,
  KTClive_4h = novel_up_KTClive_4h,
  KTCman_4h = novel_up_KTCman_4h,
  KTCbg_24h = novel_up_KTCbg_24h,
  KTClive_24h = novel_up_KTClive_24h,
  KTCman_24h = novel_up_KTCman_24h
)

novel_up_upset_data <- fromList(novel_up_sets)

upset(
  novel_up_upset_data,
  sets = names(novel_up_sets),
  nsets = 6,
  nintersects = 20,
  keep.order = TRUE,
  order.by = "freq",
  show.numbers = "yes",
  main.bar.color = "#D73027",
  sets.bar.color = "#D73027",
  matrix.color = "#000000",
  mainbar.y.label = "Shared upregulated novel lncRNAs",
  sets.x.label = "Upregulated novel lncRNAs per contrast",
  text.scale = 1.2
)


############################################################
# Downregulated novel-lncRNA overlap

novel_down_sets <- list(
  KTCbg_4h = novel_down_KTCbg_4h,
  KTClive_4h = novel_down_KTClive_4h,
  KTCman_4h = novel_down_KTCman_4h,
  KTCbg_24h = novel_down_KTCbg_24h,
  KTClive_24h = novel_down_KTClive_24h,
  KTCman_24h = novel_down_KTCman_24h
)

novel_down_upset_data <- fromList(novel_down_sets)

upset(
  novel_down_upset_data,
  sets = names(novel_down_sets),
  nsets = 6,
  nintersects = 20,
  keep.order = TRUE,
  order.by = "freq",
  show.numbers = "yes",
  main.bar.color = "#0072B2",
  sets.bar.color = "#0072B2",
  matrix.color = "#000000",
  mainbar.y.label = "Shared downregulated novel lncRNAs",
  sets.x.label = "Downregulated novel lncRNAs per contrast",
  text.scale = 1.2
)



#total number of significantly expressed up-down novel lncRNA

novel_DE_genes <- unique(
  c(
    novel_KTCbg_4h,
    novel_KTClive_4h,
    novel_KTCman_4h,
    novel_KTCbg_24h,
    novel_KTClive_24h,
    novel_KTCman_24h
  )
)

length(novel_DE_genes)



#Extract unique DE protein-coding genes
############################################################
# Extract DE protein-coding genes from each contrast

pc_KTCbg_4h <- KTCbg_4h_table$gene_id[
  KTCbg_4h_table$analysis_class == "protein_coding" &
    KTCbg_4h_table$DE_status %in% c("Up", "Down")
]

pc_KTClive_4h <- KTClive_4h_table$gene_id[
  KTClive_4h_table$analysis_class == "protein_coding" &
    KTClive_4h_table$DE_status %in% c("Up", "Down")
]

pc_KTCman_4h <- KTCman_4h_table$gene_id[
  KTCman_4h_table$analysis_class == "protein_coding" &
    KTCman_4h_table$DE_status %in% c("Up", "Down")
]

pc_KTCbg_24h <- KTCbg_24h_table$gene_id[
  KTCbg_24h_table$analysis_class == "protein_coding" &
    KTCbg_24h_table$DE_status %in% c("Up", "Down")
]

pc_KTClive_24h <- KTClive_24h_table$gene_id[
  KTClive_24h_table$analysis_class == "protein_coding" &
    KTClive_24h_table$DE_status %in% c("Up", "Down")
]

pc_KTCman_24h <- KTCman_24h_table$gene_id[
  KTCman_24h_table$analysis_class == "protein_coding" &
    KTCman_24h_table$DE_status %in% c("Up", "Down")
]



#Combine them while removing repeated genes

protein_coding_DE_unique <- unique(c(
  pc_KTCbg_4h,
  pc_KTClive_4h,
  pc_KTCman_4h,
  pc_KTCbg_24h,
  pc_KTClive_24h,
  pc_KTCman_24h
))

length(protein_coding_DE_unique)


############################################################
# Combine upregulated novel lncRNAs across all contrasts

novel_up_unique <- unique(c(
  novel_up_KTCbg_4h,
  novel_up_KTClive_4h,
  novel_up_KTCman_4h,
  novel_up_KTCbg_24h,
  novel_up_KTClive_24h,
  novel_up_KTCman_24h
))

length(novel_up_unique)

############################################################
# Combine downregulated novel lncRNAs across all contrasts

novel_down_unique <- unique(c(
  novel_down_KTCbg_4h,
  novel_down_KTClive_4h,
  novel_down_KTCman_4h,
  novel_down_KTCbg_24h,
  novel_down_KTClive_24h,
  novel_down_KTCman_24h
))

length(novel_down_unique)


############################################################
# All unique DE novel lncRNAs

novel_DE_unique <- unique(c(
  novel_up_unique,
  novel_down_unique
))

length(novel_DE_unique)


#the two direction-switching genes

direction_switching_novel <- intersect(
  novel_up_unique,
  novel_down_unique
)

length(direction_switching_novel)

direction_switching_novel



############################################################
# Prepare expression data for coexpression analysis

library(edgeR)
library(limma)

logCPM_coexpression <- cpm(
  dge,
  log=TRUE,
  prior.count=2
)

stopifnot(
  identical(
    colnames(logCPM_coexpression),
    metadata$sample_id
  )
)



#Remove donor effects while preserving experimental groups
#This removes systematic differences among donors A, B and C 
#while retaining stimulus-time responses.
group_design_coexpression <- model.matrix(
  ~ 0 + group,
  data=metadata
)

logCPM_donor_adjusted <- removeBatchEffect(
  logCPM_coexpression,
  batch=metadata$donor,
  design=group_design_coexpression
)



#Confirm the selected gene IDs
protein_coding_ids_present <- protein_coding_DE_unique[
  protein_coding_DE_unique %in% rownames(logCPM_donor_adjusted)
]

novel_lncRNA_ids_present <- novel_DE_unique[
  novel_DE_unique %in% rownames(logCPM_donor_adjusted)
]

length(protein_coding_ids_present)


length(novel_lncRNA_ids_present)


#Create the two expression matrices

protein_coding_expression <- logCPM_donor_adjusted[
  protein_coding_ids_present,
  ,
  drop=FALSE
]

novel_lncRNA_expression <- logCPM_donor_adjusted[
  novel_lncRNA_ids_present,
  ,
  drop=FALSE
]


#Check their dimensions
dim(protein_coding_expression)


dim(novel_lncRNA_expression)


stopifnot(
  identical(
    colnames(protein_coding_expression),
    colnames(novel_lncRNA_expression)
  )
)



#####################################################
#Calculate the correlation matrix
############################################################
# Protein coding - novel lncRNA correlations

coexpression_cor <- cor(
  t(novel_lncRNA_expression),
  t(protein_coding_expression),
  method="pearson"
)

dim(coexpression_cor)

#number of possible correlations

length(coexpression_cor)     

#Calculate correlation P-values

number_samples <- ncol(novel_lncRNA_expression)

number_donor_parameters <- nlevels(metadata$donor) - 1

correlation_df <- number_samples -
  number_donor_parameters -
  2

correlation_df



#Calculate the P-values
correlation_t <- coexpression_cor *
  sqrt(
    correlation_df /
      pmax(
        1 - coexpression_cor^2,
        .Machine$double.eps
      )
  )

coexpression_pvalue <- 2 * pt(
  -abs(correlation_t),
  df=correlation_df
)



#Correct for all 216,750 tests
coexpression_FDR <- matrix(
  p.adjust(
    as.vector(coexpression_pvalue),
    method="BH"
  ),
  nrow=nrow(coexpression_cor),
  ncol=ncol(coexpression_cor),
  dimnames=dimnames(coexpression_cor)
)



#Identify strong coexpression pairs
selected_pairs <- (
  abs(coexpression_cor) >= 0.70 &
    coexpression_FDR < 0.05
)



#Count the results

sum(selected_pairs)

#Positively correlated pairs
sum(
  selected_pairs &
    coexpression_cor > 0
)

#negatively correlated pairs
sum(
  selected_pairs &
    coexpression_cor < 0
)



#Making a results table

selected_indices <- which(
  selected_pairs,
  arr.ind=TRUE
)

coexpression_results <- data.frame(
  novel_lncRNA=rownames(coexpression_cor)[
    selected_indices[, "row"]
  ],
  
  protein_coding_gene=colnames(coexpression_cor)[
    selected_indices[, "col"]
  ],
  
  correlation=coexpression_cor[selected_indices],
  
  PValue=coexpression_pvalue[selected_indices],
  
  FDR=coexpression_FDR[selected_indices],
  
  stringsAsFactors=FALSE
)



# Determine the sorting order

sort_position <- base::order(
  coexpression_results$FDR,
  -abs(coexpression_results$correlation)
)


# Sort the results table

coexpression_results <- coexpression_results[
  sort_position,
  ,
  drop=FALSE
]

rownames(coexpression_results) <- NULL



# Inspect the results

dim(coexpression_results)

head(coexpression_results, 20)



############################################################
# Add protein-coding gene names

coding_gene_position <- match(
  coexpression_results$protein_coding_gene,
  annotation$gene_id
)

coexpression_results$protein_coding_gene_name <-
  annotation$gene_name[coding_gene_position]


#Replace missing names with their Ensembl IDs
missing_gene_name <- (
  is.na(coexpression_results$protein_coding_gene_name) |
    coexpression_results$protein_coding_gene_name == ""
)

coexpression_results$protein_coding_gene_name[
  missing_gene_name
] <- coexpression_results$protein_coding_gene[
  missing_gene_name
]



#Arrange the columns
coexpression_results <- coexpression_results[
  ,
  c(
    "novel_lncRNA",
    "protein_coding_gene",
    "protein_coding_gene_name",
    "correlation",
    "PValue",
    "FDR"
  )
]


#Check the annotated result

head(coexpression_results, 20)

sum(is.na(
  coexpression_results$protein_coding_gene_name           
))


#Summarize the selected pairs

coexpression_summary <- data.frame(
  possible_pairs=length(coexpression_cor),
  
  selected_pairs=nrow(coexpression_results),
  
  unique_novel_lncRNAs=length(
    unique(coexpression_results$novel_lncRNA)
  ),
  
  unique_protein_coding_genes=length(
    unique(coexpression_results$protein_coding_gene)
  ),
  
  positive_pairs=sum(
    coexpression_results$correlation > 0
  ),
  
  negative_pairs=sum(
    coexpression_results$correlation < 0
  )
)

coexpression_summary



#Calculate the percentage selected

100 * nrow(coexpression_results) /
  length(coexpression_cor)


#Count coding partners for each novel lncRNA


partners_per_lncRNA <- table(
  coexpression_results$novel_lncRNA
)

partners_per_lncRNA <- sort(
  partners_per_lncRNA,
  decreasing=TRUE
)

head(partners_per_lncRNA, 20)



#Check whether any of the 102 lncRNAs had no selected partner
lncRNAs_with_partners <- unique(
  coexpression_results$novel_lncRNA
)

lncRNAs_without_partners <- setdiff(
  novel_lncRNA_ids_present,
  lncRNAs_with_partners
)

length(lncRNAs_without_partners)

lncRNAs_without_partners



############################################################
# Select novel lncRNAs for the heatmap (15 to highly connected novel lncRNA)

top_lncRNAs <- names(
  partners_per_lncRNA
)[1:15]

heatmap_lncRNAs <- top_lncRNAs

heatmap_lncRNAs
length(heatmap_lncRNAs)    


#Count coding-gene links to these 15 lncRNAs
############################################################
# Correlations involving the selected lncRNAs

selected_lncRNA_correlations <- coexpression_cor[
  heatmap_lncRNAs,
  ,
  drop=FALSE
]

selected_lncRNA_pairs <- selected_pairs[
  heatmap_lncRNAs,
  ,
  drop=FALSE
]

coding_link_counts <- colSums(
  selected_lncRNA_pairs
)




#Calculating the mean absolute correlation 
#For example, genes connected to 12 of the 15 lncRNAs rank above genes connected to only 10. 
#If two genes both have 12 links, the one with the stronger mean absolute correlation ranks first.

coding_mean_absolute_correlation <- colMeans(
  abs(selected_lncRNA_correlations)
)



#Rank and selecting 25 protein-coding genes

coding_rank_position <- base::order(
  -coding_link_counts,
  -coding_mean_absolute_correlation
)

heatmap_coding_ids <- colnames(
  selected_lncRNA_correlations
)[
  coding_rank_position[1:25]
]

heatmap_coding_ids
length(heatmap_coding_ids)         


#Reason why these genes were selected

heatmap_coding_selection <- data.frame(
  gene_id=heatmap_coding_ids,
  
  significant_links=coding_link_counts[
    heatmap_coding_ids
  ],
  
  mean_absolute_correlation=
    coding_mean_absolute_correlation[
      heatmap_coding_ids
    ],
  
  stringsAsFactors=FALSE
)

heatmap_coding_selection


#Add protein-coding gene names


heatmap_coding_names <- annotation$gene_name[
  match(
    heatmap_coding_ids,
    annotation$gene_id
  )
]

missing_heatmap_names <- (
  is.na(heatmap_coding_names) |
    heatmap_coding_names == ""
)

heatmap_coding_names[
  missing_heatmap_names
] <- heatmap_coding_ids[
  missing_heatmap_names
]

heatmap_coding_labels <- make.unique(
  heatmap_coding_names
)


#Construct the heatmap matrix

############################################################
# Correlation matrix for the selected genes

correlation_heatmap_matrix <- coexpression_cor[
  heatmap_lncRNAs,
  heatmap_coding_ids,
  drop=FALSE
]

dim(correlation_heatmap_matrix)         


#Markinng correlations that passed both criteria

heatmap_selected_pairs <- selected_pairs[
  heatmap_lncRNAs,
  heatmap_coding_ids,
  drop=FALSE
]

heatmap_markers <- matrix(
  "",
  nrow=nrow(correlation_heatmap_matrix),
  ncol=ncol(correlation_heatmap_matrix),
  dimnames=dimnames(correlation_heatmap_matrix)
)

heatmap_markers[
  heatmap_selected_pairs
] <- "*"

sum(heatmap_selected_pairs)



############################################################
# Plot the targeted correlation heatmap

library(pheatmap)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "white", "#B2182B")
)(101)

heatmap_breaks <- seq(
  -1,
  1,
  length.out=102
)

pheatmap(
  correlation_heatmap_matrix,
  
  color=heatmap_colors,
  breaks=heatmap_breaks,
  
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  
  scale="none",
  border_color=NA,
  
  labels_col=heatmap_coding_labels,
  
  display_numbers=heatmap_markers,
  number_color="black",
  fontsize_number=9,
  
  fontsize_row=9,
  fontsize_col=8,
  angle_col=45,
  
  main=paste0(
    "Protein-coding–novel lncRNA coexpression\n",
    "* = |r| ≥ 0.70 and FDR < 0.05"
  )
)



#Extract the pairs displayed in the heatmap
heatmap_pair_results <- coexpression_results[
  coexpression_results$novel_lncRNA %in%
    heatmap_lncRNAs &
    
    coexpression_results$protein_coding_gene %in%
    heatmap_coding_ids,
  ,
  drop=FALSE
]

dim(heatmap_pair_results)

head(
  heatmap_pair_results,
  20
)



#Extract all pairs involving the 15 selected lncRNAs
############################################################
# All strong pairs involving the selected 15 lncRNAs

top_lncRNA_all_pairs <- coexpression_results[
  coexpression_results$novel_lncRNA %in%
    heatmap_lncRNAs,
  ,
  drop=FALSE
]

rownames(top_lncRNA_all_pairs) <- NULL

dim(top_lncRNA_all_pairs)

head(top_lncRNA_all_pairs, 20)



#Extract unique protein-coding partners

all_coding_partners <- unique(
  top_lncRNA_all_pairs$protein_coding_gene
)

positive_coding_partners <- unique(
  top_lncRNA_all_pairs$protein_coding_gene[
    top_lncRNA_all_pairs$correlation > 0
  ]
)

negative_coding_partners <- unique(
  top_lncRNA_all_pairs$protein_coding_gene[
    top_lncRNA_all_pairs$correlation < 0
  ]
)



#Count them

coding_partner_summary <- data.frame(
  category=c(
    "All partners",
    "Positive partners",
    "Negative partners"
  ),
  
  number=c(
    length(all_coding_partners),
    length(positive_coding_partners),
    length(negative_coding_partners)
  )
)

coding_partner_summary




############################################################
# Prepare enrichment target and background IDs

enrichment_gene_ids <- unique(
  sub(
    "\\..*$",
    "",
    all_coding_partners
  )
)

enrichment_background_ids <- unique(
  sub(
    "\\..*$",
    "",
    protein_coding_DE_unique
  )
)



#validation

length(enrichment_gene_ids)

#1370 unique pprotein coding ids were connected to those 15 lncRNAs

length(enrichment_background_ids)

all(
  enrichment_gene_ids %in%
    enrichment_background_ids
)



############################################
#GO BP and KEGG analysis 
#using 1370 target proteins against the 2125 background proteins
############################################


library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)


#Map the background genes

background_mapping <- bitr(
  enrichment_background_ids,
  fromType="ENSEMBL",
  toType=c("ENTREZID", "SYMBOL"),
  OrgDb=org.Hs.eg.db
)

#Derive the target mapping from the same table
target_mapping <- background_mapping[
  background_mapping$ENSEMBL %in%
    enrichment_gene_ids,
  ,
  drop=FALSE
]

background_entrez <- unique(
  background_mapping$ENTREZID
)

target_entrez <- unique(
  target_mapping$ENTREZID
)


#validation

c(
  target_Ensembl_input=length(enrichment_gene_ids),
  target_Ensembl_mapped=length(unique(target_mapping$ENSEMBL)),
  target_Entrez=length(target_entrez),
  background_Ensembl_input=length(enrichment_background_ids),
  background_Ensembl_mapped=length(unique(background_mapping$ENSEMBL)),
  background_Entrez=length(background_entrez)
)

stopifnot(
  all(target_entrez %in% background_entrez)
)



#GO Biological Process enrichment

go_bp <- enrichGO(
  gene=target_entrez,
  universe=background_entrez,
  OrgDb=org.Hs.eg.db,
  keyType="ENTREZID",
  ont="BP",
  pAdjustMethod="BH",
  pvalueCutoff=0.05,
  qvalueCutoff=0.20,
  minGSSize=10,
  maxGSSize=500,
  readable=TRUE
)

go_bp_results <- as.data.frame(go_bp)

go_bp_results <- go_bp_results[
  order(go_bp_results$p.adjust),
  ,
  drop=FALSE
]

dim(go_bp_results)

head(go_bp_results, 20)


#Display dot plots

dotplot(
  go_bp,
  showCategory=20,
  font.size=10
) +
  ggtitle(
    "GO biological processes associated with\nnovel lncRNA coexpression partners"
  )



#Reactome enrichment

reactome_result <- enrichPathway(
  gene=target_entrez,
  universe=background_entrez,
  organism="human",
  pvalueCutoff=0.05,
  pAdjustMethod="BH",
  qvalueCutoff=0.05,
  minGSSize=10,
  maxGSSize=500,
  readable=TRUE
)


# Convert results to a table

reactome_results <- as.data.frame(
  reactome_result
)


# Arrange by adjusted P-value

reactome_results <- reactome_results[
  order(reactome_results$p.adjust),
  ,
  drop=FALSE
]


# Number of significant Reactome pathways

number_reactome_significant <- nrow(
  reactome_results
)

number_reactome_significant


############################################################

# results and plot 

if (number_reactome_significant > 0) {
  
  print(
    head(
      reactome_results[
        ,
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "FoldEnrichment",
          "pvalue",
          "p.adjust",
          "qvalue"
        )
      ],
      20
    )
  )
  
  reactome_plot <- dotplot(
    reactome_result,
    showCategory=min(
      20,
      number_reactome_significant
    ),
    font.size=10,
    color="p.adjust"
  ) +
    ggtitle(
      "Significantly enriched Reactome pathways associated with\nnovel lncRNA coexpression partners"
    )
  
  print(reactome_plot)
  
} else {
  
  message(
    "No Reactome pathways were significant at BH FDR < 0.05."
  )
  
}





#KEGG
#############################################
#global strict KEGG enrichment using the combined coding partners of the 15 lncRNAs.
############################################################

# Prepare Entrez IDs for global KEGG enrichment

target_entrez <- unique(
  as.character(target_entrez)
)

target_entrez <- target_entrez[
  !is.na(target_entrez) &
    target_entrez != ""
]


background_entrez <- unique(
  as.character(background_entrez)
)

background_entrez <- background_entrez[
  !is.na(background_entrez) &
    background_entrez != ""
]


# Ensure every target gene belongs to the background

target_entrez <- intersect(
  target_entrez,
  background_entrez
)


# Check input numbers

length(target_entrez)

length(background_entrez)

all(
  target_entrez %in% background_entrez
)


############################################################

# Strict global KEGG enrichment

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

kegg_result <- enrichKEGG(
  gene=target_entrez,
  universe=background_entrez,
  organism="hsa",
  keyType="ncbi-geneid",
  pvalueCutoff=0.05,
  pAdjustMethod="BH",
  qvalueCutoff=0.05,
  minGSSize=10,
  maxGSSize=500,
  use_internal_data=FALSE
)


############################################################

# Convert KEGG results to a table

kegg_results <- as.data.frame(
  kegg_result
)


# Arrange significant pathways by adjusted P-value

if (nrow(kegg_results) > 0) {
  
  kegg_results <- kegg_results[
    order(kegg_results$p.adjust),
    ,
    drop=FALSE
  ]
  
}


# Number of significant KEGG pathways

number_kegg_significant <- nrow(
  kegg_results
)

number_kegg_significant


############################################################

# Display and plot significant KEGG pathways

if (number_kegg_significant > 0) {
  
  print(
    head(
      kegg_results[
        ,
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "pvalue",
          "p.adjust",
          "qvalue",
          "Count"
        )
      ],
      20
    )
  )
  
  
  kegg_plot <- dotplot(
    kegg_result,
    showCategory=min(
      20,
      number_kegg_significant
    ),
    font.size=10,
    color="p.adjust"
  ) +
    ggtitle(
      "Significantly enriched KEGG pathways associated with\nnovel lncRNA coexpression partners"
    )
  
  print(kegg_plot)
  
} else {
  
  message(
    "No KEGG pathways were significant at BH FDR and q-value < 0.05."
  )
  
}





########################################################

#GO BP and KEGG analysis of the 15 lncRNAs individually

############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


############################################################
# Step 1: Select the 15 novel lncRNAs

selected_lncRNAs <- unique(
  as.character(top_lncRNAs)
)


# Check selected lncRNAs

length(selected_lncRNAs)

selected_lncRNAs

############################################################
# Step 2: Prepare background Ensembl IDs

background_ensembl <- unique(
  sub(
    "\\.[0-9]+$",
    "",
    as.character(enrichment_background_ids)
  )
)


# Check the background

length(enrichment_background_ids)

length(background_ensembl)

head(background_ensembl)

length(
  setdiff(
    enrichment_background_ids,
    background_ensembl
  )
)


############################################################
# Step 3: Map background Ensembl IDs to Entrez IDs

background_mapping_all <- bitr(
  background_ensembl,
  fromType="ENSEMBL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)


# Check the initial mapping

dim(background_mapping_all)

head(background_mapping_all)


# Number of input Ensembl genes mapped at least once

length(
  unique(
    background_mapping_all$ENSEMBL
  )
)


# Number of Ensembl genes that did not map

length(
  setdiff(
    background_ensembl,
    background_mapping_all$ENSEMBL
  )
)


#The extra mapping rows occur because some Ensembl genes mapped to multiple Entrez IDs


############################################################
# Step 4: Identify one-to-many mappings

# Remove any completely duplicated mapping rows

background_mapping_all <- unique(
  background_mapping_all[
    ,
    c(
      "ENSEMBL",
      "ENTREZID"
    ),
    drop=FALSE
  ]
)


# Count the Entrez mappings for each Ensembl gene

mapping_frequency <- table(
  background_mapping_all$ENSEMBL
)


# Identify Ensembl genes with more than one Entrez mapping

ambiguous_ensembl <- names(
  mapping_frequency[
    mapping_frequency > 1
  ]
)


# Number of ambiguous Ensembl genes

length(ambiguous_ensembl)


# Display the ambiguous Ensembl IDs

ambiguous_ensembl


# Display their Ensembl-to-Entrez mappings

ambiguous_mapping_rows <- background_mapping_all[
  which(
    background_mapping_all$ENSEMBL %in%
      ambiguous_ensembl
  ),
  c(
    "ENSEMBL",
    "ENTREZID"
  ),
  drop=FALSE
]


ambiguous_mapping_rows


############################################################
# Step 5: Remove ambiguous background mappings

keep_unambiguous_mapping <- (
  !background_mapping_all$ENSEMBL %in%
    ambiguous_ensembl
)


background_mapping <- background_mapping_all[
  which(keep_unambiguous_mapping),
  c(
    "ENSEMBL",
    "ENTREZID"
  ),
  drop=FALSE
]


# Remove any duplicate mapping rows

background_mapping <- unique(
  background_mapping
)


############################################################
# Create the final Entrez background

background_entrez <- unique(
  as.character(
    background_mapping$ENTREZID
  )
)


background_entrez <- background_entrez[
  !is.na(background_entrez) &
    background_entrez != ""
]


############################################################
# Check the final background

dim(background_mapping)


length(
  unique(
    background_mapping$ENSEMBL
  )
)


length(background_entrez)


all(
  !background_mapping$ENSEMBL %in%
    ambiguous_ensembl
)



############################################################
# Step 6: Test partner mapping for the first lncRNA

test_lncRNA <- selected_lncRNAs[1]

test_lncRNA


############################################################
# Extract its protein-coding partners

test_partner_ensembl <- unique(
  as.character(
    coexpression_results$protein_coding_gene[
      coexpression_results$novel_lncRNA ==
        test_lncRNA
    ]
  )
)


# Remove Ensembl version numbers

test_partner_ensembl <- unique(
  sub(
    "\\.[0-9]+$",
    "",
    test_partner_ensembl
  )
)


# Original number of coding partners

length(test_partner_ensembl)


############################################################
# Map partners using the cleaned background mapping

keep_test_partner_mapping <- (
  background_mapping$ENSEMBL %in%
    test_partner_ensembl
)


test_partner_mapping <- background_mapping[
  which(keep_test_partner_mapping),
  c(
    "ENSEMBL",
    "ENTREZID"
  ),
  drop=FALSE
]


# Create the Entrez target list

test_partner_entrez <- unique(
  as.character(
    test_partner_mapping$ENTREZID
  )
)


############################################################
# Check the partner mapping

length(
  unique(
    test_partner_mapping$ENSEMBL
  )
)


length(test_partner_entrez)


all(
  test_partner_entrez %in%
    background_entrez
)



############################################################
# Step 7: Prepare partner Entrez lists for all 15 lncRNAs

partner_entrez_lists <- list()


partner_mapping_summary <- data.frame(
  novel_lncRNA=character(),
  coding_partners=integer(),
  mapped_partners=integer(),
  Entrez_IDs_used=integer(),
  stringsAsFactors=FALSE
)


############################################################
# Process the partners of all 15 lncRNAs

for (lncRNA_id in selected_lncRNAs) {
  
  # Extract protein-coding partners
  
  partner_ensembl <- unique(
    as.character(
      coexpression_results$protein_coding_gene[
        coexpression_results$novel_lncRNA ==
          lncRNA_id
      ]
    )
  )
  
  
  # Remove Ensembl version numbers
  
  partner_ensembl <- unique(
    sub(
      "\\.[0-9]+$",
      "",
      partner_ensembl
    )
  )
  
  
  # Map through the cleaned background
  
  keep_partner_mapping <- (
    background_mapping$ENSEMBL %in%
      partner_ensembl
  )
  
  
  partner_mapping <- background_mapping[
    which(keep_partner_mapping),
    c(
      "ENSEMBL",
      "ENTREZID"
    ),
    drop=FALSE
  ]
  
  
  # Create the Entrez list
  
  partner_entrez <- unique(
    as.character(
      partner_mapping$ENTREZID
    )
  )
  
  
  partner_entrez <- intersect(
    partner_entrez,
    background_entrez
  )
  
  
  # Save the Entrez list for this lncRNA
  
  partner_entrez_lists[[lncRNA_id]] <-
    partner_entrez
  
  
  # Create one summary row
  
  one_summary_row <- data.frame(
    
    novel_lncRNA=lncRNA_id,
    
    coding_partners=
      length(partner_ensembl),
    
    mapped_partners=
      length(
        unique(
          partner_mapping$ENSEMBL
        )
      ),
    
    Entrez_IDs_used=
      length(partner_entrez),
    
    stringsAsFactors=FALSE
  )
  
  
  # Add the row to the summary
  
  partner_mapping_summary <- rbind(
    partner_mapping_summary,
    one_summary_row
  )
  
}


############################################################
# Verify results for all 15 lncRNAs

rownames(partner_mapping_summary) <- NULL


partner_mapping_summary


length(partner_entrez_lists)


stopifnot(
  nrow(partner_mapping_summary) == 15,
  
  length(partner_entrez_lists) == 15,
  
  all(
    partner_mapping_summary$mapped_partners <=
      partner_mapping_summary$coding_partners
  ),
  
  all(
    partner_mapping_summary$Entrez_IDs_used <=
      partner_mapping_summary$mapped_partners
  )
)


View(partner_mapping_summary)



############################################################
# Step 8: GO BP enrichment for all 15 lncRNAs

go_results_list <- list()

significant_go_list <- list()


go_enrichment_summary <- data.frame(
  novel_lncRNA=character(),
  Entrez_IDs_used=integer(),
  significant_GO_BP=integer(),
  stringsAsFactors=FALSE
)


############################################################
# Run GO BP enrichment

for (lncRNA_id in selected_lncRNAs) {
  
  message(
    "GO BP analysis: ",
    lncRNA_id
  )
  
  
  # Retrieve this lncRNA's coding-partner Entrez IDs
  
  partner_entrez <- partner_entrez_lists[[lncRNA_id]]
  
  # Run strict GO BP enrichment
  
  go_result <- enrichGO(
    gene=partner_entrez,
    universe=background_entrez,
    OrgDb=org.Hs.eg.db,
    keyType="ENTREZID",
    ont="BP",
    pvalueCutoff=0.05,
    pAdjustMethod="BH",
    qvalueCutoff=0.05,
    minGSSize=10,
    maxGSSize=500,
    readable=TRUE
  )
  
  
  # Convert the result to a data frame
  
  go_table <- as.data.frame(
    go_result
  )
  
  
  # Arrange significant terms
  
  if (nrow(go_table) > 0) {
    
    go_table$novel_lncRNA <-
      lncRNA_id
    
    go_table$Entrez_IDs_used <-
      length(partner_entrez)
    
    
    go_order <- base::order(
      go_table$p.adjust
    )
    
    
    go_table <- go_table[
      go_order,
      colnames(go_table),
      drop=FALSE
    ]
    
    
    significant_go_list[[lncRNA_id]] <-
      go_table
    
  }
  
  
  # Store results, including empty results
  
  go_results_list[[lncRNA_id]] <-
    go_table
  
  
  # Add one row to the GO summary
  
  one_go_summary <- data.frame(
    novel_lncRNA=lncRNA_id,
    Entrez_IDs_used=length(partner_entrez),
    significant_GO_BP=nrow(go_table),
    stringsAsFactors=FALSE
  )
  
  
  go_enrichment_summary <- rbind(
    go_enrichment_summary,
    one_go_summary
  )
  
}


############################################################
# Check the GO summary

rownames(go_enrichment_summary) <- NULL


go_enrichment_summary


length(go_results_list)


sum(
  go_enrichment_summary$significant_GO_BP > 0
)


stopifnot(
  nrow(go_enrichment_summary) == 15,
  length(go_results_list) == 15
)


View(go_enrichment_summary)



nrow(go_enrichment_summary)

length(go_results_list)


############################################################
# Create GO counts required for the following Step 

stored_GO_counts <- vapply(
  go_results_list,
  nrow,
  integer(1)
)


# Check stored counts

stored_GO_counts


stopifnot(
  length(stored_GO_counts) == 15,
  
  all(
    unname(stored_GO_counts) ==
      go_enrichment_summary$significant_GO_BP
  )
)





############################################################
# Step 9: Combine significant GO BP results

# Identify lncRNAs with significant GO BP terms

go_lncRNAs_with_results <- names(
  stored_GO_counts[
    stored_GO_counts > 0
  ]
)


go_lncRNAs_with_results

length(go_lncRNAs_with_results)




############################################################
# Combine significant GO BP tables

lncRNA_GO_BP_results <- do.call(
  rbind,
  go_results_list[
    go_lncRNAs_with_results
  ]
)


# Remove inherited row names

rownames(
  lncRNA_GO_BP_results
) <- NULL


############################################################
# Arrange all results by adjusted P-value

go_all_order <- base::order(
  lncRNA_GO_BP_results$p.adjust
)


lncRNA_GO_BP_results <- lncRNA_GO_BP_results[
  go_all_order,
  colnames(lncRNA_GO_BP_results),
  drop=FALSE
]


############################################################
# Verify the combined table

dim(
  lncRNA_GO_BP_results
)


table(
  lncRNA_GO_BP_results$novel_lncRNA
)


stopifnot(
  
  nrow(lncRNA_GO_BP_results) ==
    sum(stored_GO_counts),
  
  length(
    unique(
      lncRNA_GO_BP_results$novel_lncRNA
    )
  ) ==
    length(go_lncRNAs_with_results)
)


############################################################
# Display the top GO BP results

head(
  lncRNA_GO_BP_results[
    ,
    c(
      "novel_lncRNA",
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "qvalue",
      "Count"
    ),
    drop=FALSE
  ],
  30
)


View(
  lncRNA_GO_BP_results
)




############################################################
# Create compareCluster GO BP object for plotting

go_partner_lists <- partner_entrez_lists[
  selected_lncRNAs
]


go_compare_result <- compareCluster(
  geneCluster=go_partner_lists,
  fun="enrichGO",
  universe=background_entrez,
  OrgDb=org.Hs.eg.db,
  keyType="ENTREZID",
  ont="BP",
  pvalueCutoff=0.05,
  pAdjustMethod="BH",
  qvalueCutoff=0.05,
  minGSSize=10,
  maxGSSize=500,
  readable=TRUE
)


# Convert to a table for verification

go_compare_results <- as.data.frame(
  go_compare_result
)


# Check the results

c(
  significant_associations=
    nrow(go_compare_results),
  
  lncRNAs_with_GO_BP=
    length(
      unique(
        go_compare_results$Cluster
      )
    )
)





############################################################
# GO BP dot plot

go_bp_dotplot <- enrichplot::dotplot(
  go_compare_result,
  x="Cluster",
  color="p.adjust",
  showCategory=3,
  by="geneRatio",
  includeAll=FALSE,
  font.size=9,
  title="GO biological processes associated with novel lncRNA coding partners"
) +
  labs(
    x="Novel lncRNA",
    y="GO biological process",
    color="BH-adjusted P value",
    size="Gene ratio"
  ) +
  theme_bw() +
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    ),
    plot.title=element_text(
      hjust=0.5,
      face="bold"
    )
  )


print(go_bp_dotplot)




############################################################
# Step 10: Strict KEGG enrichment for all 15 lncRNAs

kegg_partner_lists <- partner_entrez_lists[
  selected_lncRNAs
]


# Verify the input

stopifnot(
  length(kegg_partner_lists) == 15,
  identical(
    names(kegg_partner_lists),
    selected_lncRNAs
  ),
  all(
    lengths(kegg_partner_lists) > 0
  )
)


############################################################
# Run KEGG enrichment for all 15 lncRNAs

kegg_compare_result <- compareCluster(
  geneCluster=kegg_partner_lists,
  fun="enrichKEGG",
  universe=background_entrez,
  organism="hsa",
  keyType="ncbi-geneid",
  pvalueCutoff=0.05,
  pAdjustMethod="BH",
  qvalueCutoff=0.05,
  minGSSize=10,
  maxGSSize=500,
  use_internal_data=FALSE
)


############################################################
# Convert results to a table

if (is.null(kegg_compare_result)) {
  
  lncRNA_KEGG_results <- data.frame()
  
  stored_KEGG_counts <- setNames(
    rep(
      0L,
      length(selected_lncRNAs)
    ),
    selected_lncRNAs
  )
  
} else {
  
  lncRNA_KEGG_results <- as.data.frame(
    kegg_compare_result
  )
  
  
# Arrange pathways by lncRNA and adjusted P-value
  
  kegg_order <- base::order(
    lncRNA_KEGG_results$Cluster,
    lncRNA_KEGG_results$p.adjust
  )
  
  
  lncRNA_KEGG_results <- lncRNA_KEGG_results[
    kegg_order,
    colnames(lncRNA_KEGG_results),
    drop=FALSE
  ]
  
  
  rownames(
    lncRNA_KEGG_results
  ) <- NULL
  
  
# Count significant pathways for each lncRNA
  
  stored_KEGG_counts <- table(
    factor(
      lncRNA_KEGG_results$Cluster,
      levels=selected_lncRNAs
    )
  )
  
}


############################################################
# Create a summary for all 15 lncRNAs

kegg_enrichment_summary <- data.frame(
  
  novel_lncRNA=
    selected_lncRNAs,
  
  Entrez_IDs_used=
    as.integer(
      lengths(kegg_partner_lists)
    ),
  
  significant_KEGG=
    as.integer(
      stored_KEGG_counts
    ),
  
  stringsAsFactors=FALSE
)


############################################################
# Verify the results

stopifnot(
  nrow(kegg_enrichment_summary) == 15,
  
  sum(
    kegg_enrichment_summary$significant_KEGG
  ) ==
    nrow(lncRNA_KEGG_results)
)


############################################################
# Display summary

kegg_enrichment_summary

View(
  kegg_enrichment_summary
)


############################################################
# Number of lncRNAs with significant KEGG pathways

number_lncRNAs_with_KEGG <- sum(
  kegg_enrichment_summary$significant_KEGG > 0
)

number_lncRNAs_with_KEGG


############################################################
# Display significant KEGG results

if (nrow(lncRNA_KEGG_results) > 0) {
  
  print(
    lncRNA_KEGG_results[
      ,
      c(
        "Cluster",
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "qvalue",
        "Count"
      ),
      drop=FALSE
    ]
  )
  
  
  View(
    lncRNA_KEGG_results
  )
  
} else {
  
  message(
    "No significant KEGG pathways were found for any of the 15 lncRNAs."
  )
  
}


############################################################
# Verify KEGG results

c(
  significant_associations=
    nrow(lncRNA_KEGG_results),
  
  lncRNAs_with_KEGG=
    length(
      unique(
        lncRNA_KEGG_results$Cluster
      )
    ),
  
  maximum_FDR=
    max(
      lncRNA_KEGG_results$p.adjust
    ),
  
  maximum_qvalue=
    max(
      lncRNA_KEGG_results$qvalue
    )
)



# Number of unique significant KEGG pathways

length(
  unique(
    lncRNA_KEGG_results$ID
  )
)


# Number of significant pathways for each lncRNA

sort(
  table(
    lncRNA_KEGG_results$Cluster
  ),
  decreasing=TRUE
)


# Identify the four lncRNAs without significant KEGG pathways

setdiff(
  selected_lncRNAs,
  unique(
    lncRNA_KEGG_results$Cluster
  )
)


#which pathways recur across multiple lncRNAs

pathway_frequency <- sort(
  table(
    lncRNA_KEGG_results$Description
  ),
  decreasing=TRUE
)

pathway_frequency



#plot

library(enrichplot)
library(ggplot2)

kegg_dotplot <- enrichplot::dotplot(
  kegg_compare_result,
  x="Cluster",
  color="p.adjust",
  showCategory=10,
  by="geneRatio",
  includeAll=TRUE,
  font.size=9,
  title="KEGG pathways associated with novel lncRNA coding partners"
) +
  labs(
    x="Novel lncRNA",
    y="KEGG pathway",
    color="BH-adjusted P value",
    size="Gene ratio"
  ) +
  theme_bw() +
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    ),
    plot.title=element_text(
      hjust=0.5,
      face="bold"
    )
  )

print(kegg_dotplot)





######################################### Stimulus-by-time interaction analysis  ####################


library(edgeR)
library(limma)

############################################################
# Step 1: Confirm the required design coefficients

required_coefficients <- c(
  "groupcaur_KTCbg_4h",
  "groupcaur_KTCbg_24h",
  "groupcaur_KTClive_4h",
  "groupcaur_KTClive_24h",
  "groupcaur_KTCman_4h",
  "groupcaur_KTCman_24h",
  "groupRPMI_4h",
  "groupRPMI_24h"
)

required_coefficients %in% colnames(design)

stopifnot(
  all(required_coefficients %in% colnames(design))
)


############################################################
# Step 2: Define the three time-interaction contrasts

time_contrasts <- makeContrasts(
  
  KTCbg_time_change =
    (groupcaur_KTCbg_24h - groupRPMI_24h) -
    (groupcaur_KTCbg_4h - groupRPMI_4h),
  
  KTClive_time_change =
    (groupcaur_KTClive_24h - groupRPMI_24h) -
    (groupcaur_KTClive_4h - groupRPMI_4h),
  
  KTCman_time_change =
    (groupcaur_KTCman_24h - groupRPMI_24h) -
    (groupcaur_KTCman_4h - groupRPMI_4h),
  
  levels=design
)

colnames(time_contrasts)


############################################################
# Step 3: Test the beta-glucan time interaction

KTCbg_time_test <- glmQLFTest(
  fit,
  contrast=time_contrasts[, "KTCbg_time_change"]
)

KTCbg_time_table <- topTags(
  KTCbg_time_test,
  n=Inf,
  sort.by="PValue"
)$table

KTCbg_time_table <- data.frame(
  gene_id=rownames(KTCbg_time_table),
  KTCbg_time_table,
  row.names=NULL,
  check.names=FALSE
)


############################################################
# Step 4: Test the live C. auris time interaction

KTClive_time_test <- glmQLFTest(
  fit,
  contrast=time_contrasts[, "KTClive_time_change"]
)

KTClive_time_table <- topTags(
  KTClive_time_test,
  n=Inf,
  sort.by="PValue"
)$table

KTClive_time_table <- data.frame(
  gene_id=rownames(KTClive_time_table),
  KTClive_time_table,
  row.names=NULL,
  check.names=FALSE
)


############################################################
# Step 5: Test the mannan time interaction

KTCman_time_test <- glmQLFTest(
  fit,
  contrast=time_contrasts[, "KTCman_time_change"]
)

KTCman_time_table <- topTags(
  KTCman_time_test,
  n=Inf,
  sort.by="PValue"
)$table

KTCman_time_table <- data.frame(
  gene_id=rownames(KTCman_time_table),
  KTCman_time_table,
  row.names=NULL,
  check.names=FALSE
)


############################################################
# Step 6: Add gene annotation

KTCbg_time_table$gene_name <- annotation_filtered$gene_name[
  match(
    KTCbg_time_table$gene_id,
    annotation_filtered$gene_id
  )
]

KTCbg_time_table$analysis_class <- annotation_filtered$analysis_class[
  match(
    KTCbg_time_table$gene_id,
    annotation_filtered$gene_id
  )
]


KTClive_time_table$gene_name <- annotation_filtered$gene_name[
  match(
    KTClive_time_table$gene_id,
    annotation_filtered$gene_id
  )
]

KTClive_time_table$analysis_class <- annotation_filtered$analysis_class[
  match(
    KTClive_time_table$gene_id,
    annotation_filtered$gene_id
  )
]


KTCman_time_table$gene_name <- annotation_filtered$gene_name[
  match(
    KTCman_time_table$gene_id,
    annotation_filtered$gene_id
  )
]

KTCman_time_table$analysis_class <- annotation_filtered$analysis_class[
  match(
    KTCman_time_table$gene_id,
    annotation_filtered$gene_id
  )
]


############################################################
# Step 7: Extract substantial time-interaction results
# Using the same reporting criteria used in the main DE analysis
# FDR < 0.05 and absolute interaction logFC > 1

KTCbg_time_significant <- KTCbg_time_table[
  KTCbg_time_table$FDR < 0.05 &
    abs(KTCbg_time_table$logFC) > 1,
  ,
  drop=FALSE
]

KTClive_time_significant <- KTClive_time_table[
  KTClive_time_table$FDR < 0.05 &
    abs(KTClive_time_table$logFC) > 1,
  ,
  drop=FALSE
]

KTCman_time_significant <- KTCman_time_table[
  KTCman_time_table$FDR < 0.05 &
    abs(KTCman_time_table$logFC) > 1,
  ,
  drop=FALSE
]


############################################################
# Step 8: Extract time-dependent novel lncRNAs

KTCbg_time_novel <- KTCbg_time_significant[
  KTCbg_time_significant$analysis_class == "novel_lncRNA",
  ,
  drop=FALSE
]

KTClive_time_novel <- KTClive_time_significant[
  KTClive_time_significant$analysis_class == "novel_lncRNA",
  ,
  drop=FALSE
]

KTCman_time_novel <- KTCman_time_significant[
  KTCman_time_significant$analysis_class == "novel_lncRNA",
  ,
  drop=FALSE
]


############################################################
# Step 9: Create the interaction summary

time_interaction_summary <- data.frame(
  
  stimulus=c(
    "KTCbg",
    "KTClive",
    "KTCman"
  ),
  
  genes_FDR_below_0.05=c(
    sum(KTCbg_time_table$FDR < 0.05, na.rm=TRUE),
    sum(KTClive_time_table$FDR < 0.05, na.rm=TRUE),
    sum(KTCman_time_table$FDR < 0.05, na.rm=TRUE)
  ),
  
  genes_FDR_and_abs_logFC_above_1=c(
    nrow(KTCbg_time_significant),
    nrow(KTClive_time_significant),
    nrow(KTCman_time_significant)
  ),
  
  novel_lncRNAs_FDR_and_abs_logFC_above_1=c(
    nrow(KTCbg_time_novel),
    nrow(KTClive_time_novel),
    nrow(KTCman_time_novel)
  ),
  
  stringsAsFactors=FALSE
)

time_interaction_summary

View(time_interaction_summary)

View(KTCbg_time_novel)
View(KTClive_time_novel)
View(KTCman_time_novel)





############################################################
# Combine time-dependent novel lncRNA results

KTCbg_time_novel_combined <- KTCbg_time_novel
KTCbg_time_novel_combined$stimulus <- "Beta-glucan"

KTClive_time_novel_combined <- KTClive_time_novel
KTClive_time_novel_combined$stimulus <- "Live_C_auris"

KTCman_time_novel_combined <- KTCman_time_novel
KTCman_time_novel_combined$stimulus <- "Mannan"


############################################################
# Combine all three tables

all_time_novel_results <- rbind(
  KTCbg_time_novel_combined,
  KTClive_time_novel_combined,
  KTCman_time_novel_combined
)

rownames(all_time_novel_results) <- NULL


############################################################
# Add an interpretable interaction direction

all_time_novel_results$temporal_response <- ifelse(
  all_time_novel_results$logFC > 0,
  "Stronger_or_more_positive_at_24h",
  "Weaker_or_more_negative_at_24h"
)


############################################################
# Check total stimulus-specific calls

nrow(all_time_novel_results)



############################################################
# Extract unique time-dependent novel lncRNAs

time_dependent_novel_ids <- unique(
  all_time_novel_results$gene_id
)

length(time_dependent_novel_ids)



############################################################
# Find lncRNAs significant for more than one stimulus

time_interaction_frequency <- table(
  all_time_novel_results$gene_id
)

shared_time_dependent_lncRNAs <- names(
  time_interaction_frequency[
    time_interaction_frequency > 1
  ]
)

shared_time_dependent_lncRNAs


############################################################
# Display their stimulus-specific results

shared_time_interaction_results <- all_time_novel_results[
  all_time_novel_results$gene_id %in%
    shared_time_dependent_lncRNAs,
  ,
  drop=FALSE
]

shared_time_interaction_results[
  ,
  c(
    "gene_id",
    "stimulus",
    "logFC",
    "FDR",
    "temporal_response"
  ),
  drop=FALSE
]


############################################################
# Check how many belong to the selected top 15 lncRNAs

top15_time_dependent_lncRNAs <- intersect(
  selected_lncRNAs,
  time_dependent_novel_ids
)

top15_time_dependent_lncRNAs

length(top15_time_dependent_lncRNAs)


############################################################
# Reconstruct the 102 unique DE novel lncRNAs

novel_DE_102_ids <- unique(c(
  novel_up_KTCbg_4h,
  novel_down_KTCbg_4h,
  
  novel_up_KTClive_4h,
  novel_down_KTClive_4h,
  
  novel_up_KTCman_4h,
  novel_down_KTCman_4h,
  
  novel_up_KTCbg_24h,
  novel_down_KTCbg_24h,
  
  novel_up_KTClive_24h,
  novel_down_KTClive_24h,
  
  novel_up_KTCman_24h,
  novel_down_KTCman_24h
))

length(novel_DE_102_ids)




############################################################
# Interaction lncRNAs that also belong to the 102 DE set

time_and_DE_lncRNAs <- intersect(
  time_dependent_novel_ids,
  novel_DE_102_ids
)

length(time_and_DE_lncRNAs)

time_and_DE_lncRNAs


############################################################
# Interaction lncRNAs not included in the original 102

time_interaction_only_lncRNAs <- setdiff(
  time_dependent_novel_ids,
  novel_DE_102_ids
)

length(time_interaction_only_lncRNAs)

time_interaction_only_lncRNAs


############################################################
# Final count summary

stimulus_time_count_summary <- data.frame(
  category=c(
    "Expressed novel lncRNAs",
    "DE novel lncRNAs",
    "Time-interacting novel lncRNAs",
    "Both DE and time-interacting",
    "Time-interacting but not original DE",
    "Time-interacting among top 15"
  ),
  
  number=c(
    449,
    length(novel_DE_102_ids),
    length(time_dependent_novel_ids),
    length(time_and_DE_lncRNAs),
    length(time_interaction_only_lncRNAs),
    length(top15_time_dependent_lncRNAs)
  )
)

stimulus_time_count_summary



############################################################
# Temporal profiles for beta-glucan interaction lncRNAs

KTCbg_temporal_profile <- data.frame(
  
  gene_id=KTCbg_time_novel$gene_id,
  
  stimulus="Beta-glucan",
  
  logFC_4h=KTCbg_4h_table$logFC[
    match(
      KTCbg_time_novel$gene_id,
      KTCbg_4h_table$gene_id
    )
  ],
  
  FDR_4h=KTCbg_4h_table$FDR[
    match(
      KTCbg_time_novel$gene_id,
      KTCbg_4h_table$gene_id
    )
  ],
  
  logFC_24h=KTCbg_24h_table$logFC[
    match(
      KTCbg_time_novel$gene_id,
      KTCbg_24h_table$gene_id
    )
  ],
  
  FDR_24h=KTCbg_24h_table$FDR[
    match(
      KTCbg_time_novel$gene_id,
      KTCbg_24h_table$gene_id
    )
  ],
  
  interaction_logFC=KTCbg_time_novel$logFC,
  
  interaction_FDR=KTCbg_time_novel$FDR,
  
  stringsAsFactors=FALSE
)


############################################################
# Temporal profiles for live C. auris interaction lncRNAs

KTClive_temporal_profile <- data.frame(
  
  gene_id=KTClive_time_novel$gene_id,
  
  stimulus="Live_C_auris",
  
  logFC_4h=KTClive_4h_table$logFC[
    match(
      KTClive_time_novel$gene_id,
      KTClive_4h_table$gene_id
    )
  ],
  
  FDR_4h=KTClive_4h_table$FDR[
    match(
      KTClive_time_novel$gene_id,
      KTClive_4h_table$gene_id
    )
  ],
  
  logFC_24h=KTClive_24h_table$logFC[
    match(
      KTClive_time_novel$gene_id,
      KTClive_24h_table$gene_id
    )
  ],
  
  FDR_24h=KTClive_24h_table$FDR[
    match(
      KTClive_time_novel$gene_id,
      KTClive_24h_table$gene_id
    )
  ],
  
  interaction_logFC=KTClive_time_novel$logFC,
  
  interaction_FDR=KTClive_time_novel$FDR,
  
  stringsAsFactors=FALSE
)


############################################################
# Temporal profiles for mannan interaction lncRNAs

KTCman_temporal_profile <- data.frame(
  
  gene_id=KTCman_time_novel$gene_id,
  
  stimulus="Mannan",
  
  logFC_4h=KTCman_4h_table$logFC[
    match(
      KTCman_time_novel$gene_id,
      KTCman_4h_table$gene_id
    )
  ],
  
  FDR_4h=KTCman_4h_table$FDR[
    match(
      KTCman_time_novel$gene_id,
      KTCman_4h_table$gene_id
    )
  ],
  
  logFC_24h=KTCman_24h_table$logFC[
    match(
      KTCman_time_novel$gene_id,
      KTCman_24h_table$gene_id
    )
  ],
  
  FDR_24h=KTCman_24h_table$FDR[
    match(
      KTCman_time_novel$gene_id,
      KTCman_24h_table$gene_id
    )
  ],
  
  interaction_logFC=KTCman_time_novel$logFC,
  
  interaction_FDR=KTCman_time_novel$FDR,
  
  stringsAsFactors=FALSE
)


############################################################
# Combine the three stimuli

novel_lncRNA_temporal_profiles <- rbind(
  KTCbg_temporal_profile,
  KTClive_temporal_profile,
  KTCman_temporal_profile
)

rownames(novel_lncRNA_temporal_profiles) <- NULL


############################################################
# Assign DE status at each timepoint

novel_lncRNA_temporal_profiles$status_4h <- "Not_DE"

novel_lncRNA_temporal_profiles$status_4h[
  novel_lncRNA_temporal_profiles$FDR_4h < 0.05 &
    novel_lncRNA_temporal_profiles$logFC_4h > 1
] <- "Up"

novel_lncRNA_temporal_profiles$status_4h[
  novel_lncRNA_temporal_profiles$FDR_4h < 0.05 &
    novel_lncRNA_temporal_profiles$logFC_4h < -1
] <- "Down"


novel_lncRNA_temporal_profiles$status_24h <- "Not_DE"

novel_lncRNA_temporal_profiles$status_24h[
  novel_lncRNA_temporal_profiles$FDR_24h < 0.05 &
    novel_lncRNA_temporal_profiles$logFC_24h > 1
] <- "Up"

novel_lncRNA_temporal_profiles$status_24h[
  novel_lncRNA_temporal_profiles$FDR_24h < 0.05 &
    novel_lncRNA_temporal_profiles$logFC_24h < -1
] <- "Down"


############################################################
# Add previous evidence

novel_lncRNA_temporal_profiles$in_original_102 <-
  novel_lncRNA_temporal_profiles$gene_id %in%
  novel_DE_102_ids

novel_lncRNA_temporal_profiles$in_top15 <-
  novel_lncRNA_temporal_profiles$gene_id %in%
  selected_lncRNAs


############################################################
# Arrange and examine

novel_lncRNA_temporal_profiles <-
  novel_lncRNA_temporal_profiles[
    order(
      novel_lncRNA_temporal_profiles$stimulus,
      novel_lncRNA_temporal_profiles$interaction_FDR
    ),
    ,
    drop=FALSE
  ]

View(novel_lncRNA_temporal_profiles)

table(
  novel_lncRNA_temporal_profiles$status_4h,
  novel_lncRNA_temporal_profiles$status_24h
)


############################################################
# Save temporal-analysis results

temporal_dir <- file.path(
  output_dir,
  "temporal_analysis"
)

dir.create(
  temporal_dir,
  showWarnings=FALSE,
  recursive=TRUE
)


############################################################
# Save the complete temporal-profile table

write.csv(
  novel_lncRNA_temporal_profiles,
  file.path(
    temporal_dir,
    "novel_lncRNA_temporal_profiles.csv"
  ),
  row.names=FALSE
)


############################################################
# Save the 4h vs 24h status summary

temporal_status_summary <- as.data.frame(
  table(
    status_4h=novel_lncRNA_temporal_profiles$status_4h,
    status_24h=novel_lncRNA_temporal_profiles$status_24h
  )
)

write.csv(
  temporal_status_summary,
  file.path(
    temporal_dir,
    "novel_lncRNA_temporal_status_summary.csv"
  ),
  row.names=FALSE
)


############################################################
# Assign temporal classifications

novel_lncRNA_temporal_profiles$temporal_class <-
  "Interaction_only_subthreshold"

novel_lncRNA_temporal_profiles$temporal_class[
  novel_lncRNA_temporal_profiles$status_4h == "Not_DE" &
    novel_lncRNA_temporal_profiles$status_24h == "Up"
] <- "Late_induced"

novel_lncRNA_temporal_profiles$temporal_class[
  novel_lncRNA_temporal_profiles$status_4h == "Not_DE" &
    novel_lncRNA_temporal_profiles$status_24h == "Down"
] <- "Late_repressed"


############################################################
# Extract unique lncRNAs in each temporal category

late_induced_lncRNAs <- unique(
  novel_lncRNA_temporal_profiles$gene_id[
    novel_lncRNA_temporal_profiles$temporal_class ==
      "Late_induced"
  ]
)

late_repressed_lncRNAs <- unique(
  novel_lncRNA_temporal_profiles$gene_id[
    novel_lncRNA_temporal_profiles$temporal_class ==
      "Late_repressed"
  ]
)

subthreshold_interaction_lncRNAs <- unique(
  novel_lncRNA_temporal_profiles$gene_id[
    novel_lncRNA_temporal_profiles$temporal_class ==
      "Interaction_only_subthreshold"
  ]
)


############################################################
# Create unique-lncRNA summary

temporal_classification_summary <- data.frame(
  
  temporal_class=c(
    "Late induced",
    "Late repressed",
    "Interaction only, subthreshold"
  ),
  
  unique_novel_lncRNAs=c(
    length(late_induced_lncRNAs),
    length(late_repressed_lncRNAs),
    length(subthreshold_interaction_lncRNAs)
  )
)

temporal_classification_summary


############################################################
# Save the classified results

write.csv(
  novel_lncRNA_temporal_profiles,
  file.path(
    temporal_dir,
    "novel_lncRNA_temporal_profiles_classified.csv"
  ),
  row.names=FALSE
)

write.csv(
  temporal_classification_summary,
  file.path(
    temporal_dir,
    "temporal_classification_summary.csv"
  ),
  row.names=FALSE
)






############################################################
# Final 14 priority novel lncRNAs

candidate_14_lncRNAs <- intersect(
  novel_DE_102_ids,
  time_dependent_novel_ids
)

candidate_14_lncRNAs

length(candidate_14_lncRNAs)

stopifnot(
  length(candidate_14_lncRNAs) == 14
)




################################### Cis-annalysis   ##########################




# Extract coordinates of the 14 candidate lncRNAs
#############################################################
candidate_14_coordinates <- gene_info[
  match(
    candidate_14_lncRNAs,
    gene_info$Geneid
  ),
  c(
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
  ),
  drop=FALSE
]

############################################################
# Examine the coordinates

candidate_14_coordinates

View(candidate_14_coordinates)

############################################################
# Check for missing IDs or coordinates

setdiff(
  candidate_14_lncRNAs,
  gene_info$Geneid
)

colSums(
  is.na(candidate_14_coordinates)
)

############################################################
# Essential checks

stopifnot(
  nrow(candidate_14_coordinates) == 14,
  !anyNA(candidate_14_coordinates$Geneid),
  !anyDuplicated(candidate_14_coordinates$Geneid),
  !anyNA(candidate_14_coordinates$Chr),
  !anyNA(candidate_14_coordinates$Start),
  !anyNA(candidate_14_coordinates$End)
)



############################################################
# Confirm that annotation and gene_info are aligned

stopifnot(
  all(
    as.character(gene_info$Geneid) ==
      as.character(annotation$gene_id)
  )
)


############################################################
# Functions for collapsing featureCounts coordinates

get_single_value <- function(x) {
  
  values <- strsplit(
    as.character(x),
    ";",
    fixed=TRUE
  )[[1]]
  
  values <- unique(values)
  
  if (length(values) == 1) {
    return(values)
  }
  
  return(NA_character_)
}


get_minimum_coordinate <- function(x) {
  
  values <- strsplit(
    as.character(x),
    ";",
    fixed=TRUE
  )[[1]]
  
  values <- as.numeric(values)
  
  min(
    values,
    na.rm=TRUE
  )
}


get_maximum_coordinate <- function(x) {
  
  values <- strsplit(
    as.character(x),
    ";",
    fixed=TRUE
  )[[1]]
  
  values <- as.numeric(values)
  
  max(
    values,
    na.rm=TRUE
  )
}



############################################################
# Identify protein-coding rows

protein_coding_rows <- which(
  annotation$analysis_class == "protein_coding"
)

############################################################
# Create gene-level protein-coding coordinates

protein_coding_coordinates <- data.frame(
  
  coding_gene_id=
    as.character(
      annotation$gene_id[protein_coding_rows]
    ),
  
  coding_gene_name=
    as.character(
      annotation$gene_name[protein_coding_rows]
    ),
  
  Chr=vapply(
    gene_info$Chr[protein_coding_rows],
    get_single_value,
    character(1)
  ),
  
  coding_Start=vapply(
    gene_info$Start[protein_coding_rows],
    get_minimum_coordinate,
    numeric(1)
  ),
  
  coding_End=vapply(
    gene_info$End[protein_coding_rows],
    get_maximum_coordinate,
    numeric(1)
  ),
  
  coding_Strand=vapply(
    gene_info$Strand[protein_coding_rows],
    get_single_value,
    character(1)
  ),
  
  stringsAsFactors=FALSE
)



############################################################
# Identify valid coordinate rows

valid_coordinate_rows <- (
  !is.na(protein_coding_coordinates$Chr) &
    !is.na(protein_coding_coordinates$coding_Start) &
    !is.na(protein_coding_coordinates$coding_End) &
    !is.na(protein_coding_coordinates$coding_Strand) &
    is.finite(protein_coding_coordinates$coding_Start) &
    is.finite(protein_coding_coordinates$coding_End)
)

############################################################
# Number removed because of ambiguous coordinates

table(valid_coordinate_rows)

############################################################
# Keep valid rows

protein_coding_coordinates <-
  protein_coding_coordinates[
    valid_coordinate_rows,
    ,
    drop=FALSE
  ]

rownames(protein_coding_coordinates) <- NULL


#verify


dim(protein_coding_coordinates)

head(protein_coding_coordinates)

colSums(
  is.na(protein_coding_coordinates)
)

any(
  grepl(
    ";",
    protein_coding_coordinates$Chr,
    fixed=TRUE
  )
)

all(
  protein_coding_coordinates$coding_Start <=
    protein_coding_coordinates$coding_End
)


############################################################
# Prepare gene-level coordinates for the 14 lncRNAs

lncRNA_coordinates <- data.frame(
  
  lncRNA_id=
    as.character(
      candidate_14_coordinates$Geneid
    ),
  
  Chr=vapply(
    candidate_14_coordinates$Chr,
    get_single_value,
    character(1)
  ),
  
  lncRNA_Start=vapply(
    candidate_14_coordinates$Start,
    get_minimum_coordinate,
    numeric(1)
  ),
  
  lncRNA_End=vapply(
    candidate_14_coordinates$End,
    get_maximum_coordinate,
    numeric(1)
  ),
  
  lncRNA_Strand=vapply(
    candidate_14_coordinates$Strand,
    get_single_value,
    character(1)
  ),
  
  stringsAsFactors=FALSE
)


############################################################
# Check the lncRNA coordinates

dim(lncRNA_coordinates)

lncRNA_coordinates

colSums(
  is.na(lncRNA_coordinates)
)

any(
  grepl(
    ";",
    lncRNA_coordinates$Chr,
    fixed=TRUE
  )
)

all(
  lncRNA_coordinates$lncRNA_Start <=
    lncRNA_coordinates$lncRNA_End
)

stopifnot(
  nrow(lncRNA_coordinates) == 14,
  !anyNA(lncRNA_coordinates),
  !anyDuplicated(lncRNA_coordinates$lncRNA_id)
)


############################################################
# Match lncRNAs and coding genes on the same chromosome

same_chromosome_pairs <- merge(
  lncRNA_coordinates,
  protein_coding_coordinates,
  by="Chr",
  all=FALSE,
  sort=FALSE
)

dim(same_chromosome_pairs)


#The 14,930 rows are not 14,930 cis neighbors. 
#They are all possible lncRNA–protein-coding-gene pairs located on the same chromosome.


############################################################
# Calculate distance in base pairs

same_chromosome_pairs$distance_bp <- ifelse(
  
  same_chromosome_pairs$coding_End <
    same_chromosome_pairs$lncRNA_Start,
  
  same_chromosome_pairs$lncRNA_Start -
    same_chromosome_pairs$coding_End,
  
  ifelse(
    
    same_chromosome_pairs$coding_Start >
      same_chromosome_pairs$lncRNA_End,
    
    same_chromosome_pairs$coding_Start -
      same_chromosome_pairs$lncRNA_End,
    
    0
  )
)

same_chromosome_pairs$distance_kb <-
  same_chromosome_pairs$distance_bp / 1000



# Coding-gene position relative to genomic coordinates

same_chromosome_pairs$genomic_position <- ifelse(
  
  same_chromosome_pairs$coding_End <
    same_chromosome_pairs$lncRNA_Start,
  
  "Left",
  
  ifelse(
    
    same_chromosome_pairs$coding_Start >
      same_chromosome_pairs$lncRNA_End,
    
    "Right",
    
    "Overlapping"
  )
)



############################################################
# Select cis-proximal genes within 100 kb

cis_window_bp <- 100000

cis_neighbors_100kb <- same_chromosome_pairs[
  same_chromosome_pairs$distance_bp <= cis_window_bp,
  ,
  drop=FALSE
]

############################################################
# Arrange by lncRNA and distance

cis_order <- order(
  cis_neighbors_100kb$lncRNA_id,
  cis_neighbors_100kb$distance_bp
)

cis_neighbors_100kb <- cis_neighbors_100kb[
  cis_order,
  ,
  drop=FALSE
]

rownames(cis_neighbors_100kb) <- NULL


#verify
nrow(cis_neighbors_100kb)

length(
  unique(
    cis_neighbors_100kb$lncRNA_id
  )
)

table(
  cis_neighbors_100kb$lncRNA_id
)

setdiff(
  candidate_14_lncRNAs,
  unique(cis_neighbors_100kb$lncRNA_id)
)

View(cis_neighbors_100kb)




############################################################
#save the files

library(writexl)

results_dir <- "D:/a/final_analysis/results"

dir.create(
  results_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

write_xlsx(
  same_chromosome_pairs,
  path=file.path(
    results_dir,
    "same_chromosome_pairs.xlsx"
  )
)



write_xlsx(
  cis_neighbors_100kb,
  path=file.path(
    results_dir,
    "candidate_14_cis_neighbors_100kb.xlsx"
  )
)





############################################################
# Copy the cis-neighbor table

cis_neighbors_evidence <- cis_neighbors_100kb

############################################################
# Match lncRNAs and coding genes to the correlation matrices

lncRNA_matrix_position <- match(
  cis_neighbors_evidence$lncRNA_id,
  rownames(coexpression_cor)
)

coding_matrix_position <- match(
  cis_neighbors_evidence$coding_gene_id,
  colnames(coexpression_cor)
)

############################################################
# Mark coding genes belonging to the 2,125 DE coding genes

cis_neighbors_evidence$in_DE_coding_2125 <-
  !is.na(coding_matrix_position)

############################################################
# Initialize correlation columns

cis_neighbors_evidence$correlation <- NA_real_

cis_neighbors_evidence$coexpression_PValue <- NA_real_

cis_neighbors_evidence$coexpression_FDR <- NA_real_

############################################################
# Identify pairs available in the correlation matrices

available_correlation_pair <- (
  !is.na(lncRNA_matrix_position) &
    !is.na(coding_matrix_position)
)

############################################################
# Matrix positions for available pairs

matrix_positions <- cbind(
  lncRNA_matrix_position[available_correlation_pair],
  coding_matrix_position[available_correlation_pair]
)

############################################################
# Extract correlation statistics

cis_neighbors_evidence$correlation[
  available_correlation_pair
] <- coexpression_cor[
  matrix_positions
]

cis_neighbors_evidence$coexpression_PValue[
  available_correlation_pair
] <- coexpression_pvalue[
  matrix_positions
]

cis_neighbors_evidence$coexpression_FDR[
  available_correlation_pair
] <- coexpression_FDR[
  matrix_positions
]



############################################################
# Significant coexpression criteria

cis_neighbors_evidence$significant_coexpression <- (
  !is.na(cis_neighbors_evidence$correlation) &
    !is.na(cis_neighbors_evidence$coexpression_FDR) &
    abs(cis_neighbors_evidence$correlation) >= 0.70 &
    cis_neighbors_evidence$coexpression_FDR < 0.05
)

############################################################
# Add strand relationship

cis_neighbors_evidence$strand_relationship <- ifelse(
  cis_neighbors_evidence$lncRNA_Strand ==
    cis_neighbors_evidence$coding_Strand,
  "Same_strand",
  "Opposite_strand"
)

############################################################
# Examine evidence counts

table(
  cis_neighbors_evidence$in_DE_coding_2125
)

table(
  cis_neighbors_evidence$significant_coexpression
)


############################################################
# Nearby, DE and significantly coexpressed pairs

high_confidence_cis_pairs <- cis_neighbors_evidence[
  cis_neighbors_evidence$in_DE_coding_2125 &
    cis_neighbors_evidence$significant_coexpression,
  ,
  drop=FALSE
]

high_confidence_order <- base::order(
  high_confidence_cis_pairs$coexpression_FDR,
  high_confidence_cis_pairs$distance_bp
)

high_confidence_cis_pairs <- high_confidence_cis_pairs[
  high_confidence_order,
  ,
  drop=FALSE
]

rownames(high_confidence_cis_pairs) <- NULL

dim(high_confidence_cis_pairs)

View(high_confidence_cis_pairs)




writexl::write_xlsx(
  high_confidence_cis_pairs,
  path="D:/a/final_analysis/results/high_confidence_cis_pairs.xlsx"
)





####################Network Analysi###############################




############################################################
# Extract significant coexpression edges for the 14 lncRNAs

candidate_14_network_edges <- coexpression_results[
  coexpression_results$novel_lncRNA %in% candidate_14_lncRNAs,
  ,
  drop=FALSE
]

############################################################
# Check the extracted network

dim(candidate_14_network_edges)

length(
  unique(candidate_14_network_edges$novel_lncRNA)
)

length(
  unique(candidate_14_network_edges$protein_coding_gene)
)

table(
  candidate_14_network_edges$novel_lncRNA
)

############################################################
# Check whether any of the 14 lncRNAs has no significant partner

setdiff(
  candidate_14_lncRNAs,
  unique(candidate_14_network_edges$novel_lncRNA)
)




############################################################
# Create the complete Cytoscape edge table

network_edge_table <- data.frame(
  
  source=
    as.character(
      candidate_14_network_edges$novel_lncRNA
    ),
  
  interaction=
    ifelse(
      candidate_14_network_edges$correlation > 0,
      "Positive_coexpression",
      "Negative_coexpression"
    ),
  
  target=
    as.character(
      candidate_14_network_edges$protein_coding_gene
    ),
  
  correlation=
    candidate_14_network_edges$correlation,
  
  absolute_correlation=
    abs(candidate_14_network_edges$correlation),
  
  PValue=
    candidate_14_network_edges$PValue,
  
  FDR=
    candidate_14_network_edges$FDR,
  
  stringsAsFactors=FALSE
)

############################################################
# Add protein-coding gene names

network_target_clean <- sub(
  "\\.[0-9]+$",
  "",
  network_edge_table$target
)

annotation_gene_clean <- sub(
  "\\.[0-9]+$",
  "",
  as.character(annotation$gene_id)
)

target_annotation_position <- match(
  network_target_clean,
  annotation_gene_clean
)

network_edge_table$target_name <- as.character(
  annotation$gene_name[target_annotation_position]
)

# Use Ensembl ID when a gene name is unavailable

missing_target_name <- (
  is.na(network_edge_table$target_name) |
    network_edge_table$target_name == ""
)

network_edge_table$target_name[missing_target_name] <-
  network_edge_table$target[missing_target_name]

############################################################
# Create identifiers for matching cis pairs

network_edge_table$pair_id <- paste(
  network_edge_table$source,
  network_edge_table$target,
  sep="__"
)

cis_neighbors_evidence$pair_id <- paste(
  cis_neighbors_evidence$lncRNA_id,
  cis_neighbors_evidence$coding_gene_id,
  sep="__"
)

cis_match_position <- match(
  network_edge_table$pair_id,
  cis_neighbors_evidence$pair_id
)

############################################################
# Add 100-kb cis-neighbor evidence

network_edge_table$is_cis_100kb <- !is.na(
  cis_match_position
)

network_edge_table$distance_kb <-
  cis_neighbors_evidence$distance_kb[
    cis_match_position
  ]

network_edge_table$genomic_position <-
  cis_neighbors_evidence$genomic_position[
    cis_match_position
  ]

network_edge_table$strand_relationship <-
  cis_neighbors_evidence$strand_relationship[
    cis_match_position
  ]

############################################################
# Mark the four high-confidence cis pairs

high_confidence_pair_id <- paste(
  high_confidence_cis_pairs$lncRNA_id,
  high_confidence_cis_pairs$coding_gene_id,
  sep="__"
)

network_edge_table$high_confidence_cis <-
  network_edge_table$pair_id %in%
  high_confidence_pair_id

############################################################
# Mark the previous four high-priority lncRNAs

previous_priority_lncRNAs <- c(
  "MSTRG.6308",
  "MSTRG.32011",
  "MSTRG.25627",
  "MSTRG.22879"
)

network_edge_table$source_previous_priority <-
  network_edge_table$source %in%
  previous_priority_lncRNAs

############################################################
# Mark the two lncRNAs with maximum convergent evidence

maximum_evidence_lncRNAs <- c(
  "MSTRG.32011",
  "MSTRG.25627"
)

network_edge_table$source_maximum_evidence <-
  network_edge_table$source %in%
  maximum_evidence_lncRNAs

############################################################
# Arrange the columns

network_edge_table <- network_edge_table[
  ,
  c(
    "source",
    "interaction",
    "target",
    "target_name",
    "correlation",
    "absolute_correlation",
    "PValue",
    "FDR",
    "is_cis_100kb",
    "distance_kb",
    "genomic_position",
    "strand_relationship",
    "high_confidence_cis",
    "source_previous_priority",
    "source_maximum_evidence",
    "pair_id"
  ),
  drop=FALSE
]

############################################################
# Put high-confidence cis edges first

edge_order <- base::order(
  -as.integer(network_edge_table$high_confidence_cis),
  -network_edge_table$absolute_correlation
)

network_edge_table <- network_edge_table[
  edge_order,
  ,
  drop=FALSE
]

rownames(network_edge_table) <- NULL

############################################################
# Check the completed edge table

dim(network_edge_table)

length(
  unique(network_edge_table$source)
)

length(
  unique(network_edge_table$target)
)

table(network_edge_table$interaction)

table(network_edge_table$is_cis_100kb)

table(network_edge_table$high_confidence_cis)

View(network_edge_table)






############################################################
# Create the novel-lncRNA nodes

lncRNA_nodes <- data.frame(
  id=as.character(candidate_14_lncRNAs),
  label=as.character(candidate_14_lncRNAs),
  node_type="Novel_lncRNA",
  stringsAsFactors=FALSE
)

############################################################
# Create the protein-coding nodes

coding_node_ids <- unique(
  as.character(network_edge_table$target)
)

coding_name_position <- match(
  coding_node_ids,
  network_edge_table$target
)

coding_nodes <- data.frame(
  id=coding_node_ids,
  label=as.character(
    network_edge_table$target_name[
      coding_name_position
    ]
  ),
  node_type="Protein_coding",
  stringsAsFactors=FALSE
)

############################################################
# Combine both node types

network_node_table <- rbind(
  lncRNA_nodes,
  coding_nodes
)

############################################################
# Calculate node degree

degree_counts <- table(
  c(
    network_edge_table$source,
    network_edge_table$target
  )
)

degree_position <- match(
  network_node_table$id,
  names(degree_counts)
)

network_node_table$degree <- 0L

nodes_with_edges <- !is.na(
  degree_position
)

network_node_table$degree[nodes_with_edges] <-
  as.integer(
    degree_counts[
      degree_position[nodes_with_edges]
    ]
  )

############################################################
# Mark the previous four priority lncRNAs

previous_priority_lncRNAs <- c(
  "MSTRG.6308",
  "MSTRG.32011",
  "MSTRG.25627",
  "MSTRG.22879"
)

network_node_table$previous_priority <-
  network_node_table$id %in%
  previous_priority_lncRNAs

############################################################
# Mark the two maximum-evidence lncRNAs

maximum_evidence_lncRNAs <- c(
  "MSTRG.32011",
  "MSTRG.25627"
)

network_node_table$maximum_evidence <-
  network_node_table$id %in%
  maximum_evidence_lncRNAs

############################################################
# Mark nodes in the four high-confidence cis pairs

high_confidence_cis_nodes <- unique(
  c(
    as.character(
      high_confidence_cis_pairs$lncRNA_id
    ),
    as.character(
      high_confidence_cis_pairs$coding_gene_id
    )
  )
)

network_node_table$high_confidence_cis_node <-
  network_node_table$id %in%
  high_confidence_cis_nodes

############################################################
# Indicate whether each node is connected

network_node_table$connected_in_network <- (
  network_node_table$degree > 0
)

############################################################
# Put lncRNAs first and arrange by connectivity

node_order <- base::order(
  network_node_table$node_type != "Novel_lncRNA",
  -network_node_table$degree,
  network_node_table$label
)

network_node_table <- network_node_table[
  node_order,
  ,
  drop=FALSE
]

rownames(network_node_table) <- NULL

############################################################
# Check the completed node table

dim(network_node_table)

table(network_node_table$node_type)

table(network_node_table$connected_in_network)

anyDuplicated(network_node_table$id)

all(
  network_edge_table$source %in%
    network_node_table$id
)

all(
  network_edge_table$target %in%
    network_node_table$id
)

############################################################
# Examine the 14 novel lncRNA nodes

network_node_table[
  network_node_table$node_type == "Novel_lncRNA",
  ,
  drop=FALSE
]

View(network_node_table)



############################################################
# Create the network-results directory

network_output_dir <- file.path(
  "D:/a/final_analysis/results",
  "network_analysis"
)

dir.create(
  network_output_dir,
  recursive=TRUE,
  showWarnings=FALSE
)



############################################################
# Save Cytoscape edge table without quotation marks

write.csv(
  network_edge_table,
  file=file.path(
    network_output_dir,
    "candidate_14_Cytoscape_edges.csv"
  ),
  row.names=FALSE,
  quote=FALSE,
  na=""
)

############################################################
# Save Cytoscape node table without quotation marks

write.csv(
  network_node_table,
  file=file.path(
    network_output_dir,
    "candidate_14_Cytoscape_nodes.csv"
  ),
  row.names=FALSE,
  quote=FALSE,
  na=""
)




############################################################
# Create a smaller network for visualization

display_edge_candidates <- network_edge_table

############################################################
# Sort partners within each lncRNA:
# first by FDR, then by absolute correlation

display_order <- base::order(
  display_edge_candidates$source,
  display_edge_candidates$FDR,
  -display_edge_candidates$absolute_correlation
)

display_edge_candidates <- display_edge_candidates[
  display_order,
  ,
  drop=FALSE
]

############################################################
# Rank coding partners separately for each lncRNA

display_edge_candidates$partner_rank <- ave(
  seq_len(nrow(display_edge_candidates)),
  display_edge_candidates$source,
  FUN=seq_along
)

############################################################
# Keep the top five partners per lncRNA
# and retain every high-confidence cis edge

display_network_edges <- display_edge_candidates[
  display_edge_candidates$partner_rank <= 5 |
    display_edge_candidates$high_confidence_cis,
  ,
  drop=FALSE
]

rownames(display_network_edges) <- NULL

############################################################
# Check the reduced network

dim(display_network_edges)

table(display_network_edges$source)

sum(display_network_edges$high_confidence_cis)

length(
  unique(display_network_edges$source)
)

length(
  unique(display_network_edges$target)
)



############################################################
# Identify nodes present in the display network

display_node_ids <- unique(
  c(
    display_network_edges$source,
    display_network_edges$target
  )
)

############################################################
# Extract their node information

display_network_nodes <- network_node_table[
  network_node_table$id %in%
    display_node_ids,
  ,
  drop=FALSE
]

############################################################
# Calculate degree within the reduced display network

display_degree_counts <- table(
  c(
    display_network_edges$source,
    display_network_edges$target
  )
)

display_degree_position <- match(
  display_network_nodes$id,
  names(display_degree_counts)
)

display_network_nodes$display_degree <- as.integer(
  display_degree_counts[
    display_degree_position
  ]
)

############################################################
# Create a selective label column

display_network_nodes$key_label <- ""

label_key_nodes <- (
  display_network_nodes$node_type == "Novel_lncRNA" |
    display_network_nodes$high_confidence_cis_node
)

display_network_nodes$key_label[label_key_nodes] <-
  display_network_nodes$label[label_key_nodes]

rownames(display_network_nodes) <- NULL

############################################################
# Check display nodes

dim(display_network_nodes)

table(display_network_nodes$node_type)

View(display_network_nodes)





############################################################
# Save the reduced Cytoscape display network

write.csv(
  display_network_edges,
  file=file.path(
    network_output_dir,
    "candidate_14_display_edges.csv"
  ),
  row.names=FALSE,
  quote=FALSE,
  na=""
)

write.csv(
  display_network_nodes,
  file=file.path(
    network_output_dir,
    "candidate_14_display_nodes.csv"
  ),
  row.names=FALSE,
  quote=FALSE,
  na=""
)








