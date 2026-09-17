#!/bin/bash
#SBATCH --job-name=lnc_screen_counts
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/project/bishalab/msarker/c.auris_new/lncrna_filtering/08_expression_screen/logs/featureCounts_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/lncrna_filtering/08_expression_screen/logs/featureCounts_%j.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

BAM_DIR=/project/bishalab/msarker/c.auris_new/STAR_alignment
OUT_DIR=/project/bishalab/msarker/c.auris_new/lncrna_filtering/08_expression_screen/results
ANNOTATION=/project/bishalab/msarker/c.auris_new/lncrna_filtering/08_expression_screen/provisional_annotation_658.gtf

shopt -s nullglob
BAMS=("$BAM_DIR"/*Aligned.sortedByCoord.out.bam)

if [ "${#BAMS[@]}" -ne 24 ]; then
    echo "ERROR: Expected 24 BAM files but found ${#BAMS[@]}"
    exit 1
fi

featureCounts -T "$SLURM_CPUS_PER_TASK" -s 1 -t exon -g gene_id -a "$ANNOTATION" -o "$OUT_DIR/screening_gene_counts.txt" "${BAMS[@]}"

echo "featureCounts completed for ${#BAMS[@]} BAM files"