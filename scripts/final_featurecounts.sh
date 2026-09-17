#!/bin/bash
#SBATCH --job-name=final_featureCounts
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --output=/project/bishalab/msarker/c.auris_new/lncrna_filtering/09_final_catalog/final_counts/featureCounts_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/lncrna_filtering/09_final_catalog/final_counts/featureCounts_%j.err

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

featureCounts -T "$SLURM_CPUS_PER_TASK" -s 1 -t exon -g gene_id -a /project/bishalab/msarker/c.auris_new/lncrna_filtering/09_final_catalog/final_annotation_GENCODEv48_plus_novel_lncRNA.gtf -o /project/bishalab/msarker/c.auris_new/lncrna_filtering/09_final_catalog/final_counts/final_gene_counts.txt /project/bishalab/msarker/c.auris_new/STAR_alignment/*Aligned.sortedByCoord.out.bam