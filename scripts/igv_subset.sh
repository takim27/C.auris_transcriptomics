#!/bin/bash
#SBATCH --job-name=IGV_subset
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-23%6
#SBATCH --output=/project/bishalab/msarker/c.auris_new/lncrna_filtering/07_genomic_context/igv_bams/logs/IGV_subset_%A_%a.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/lncrna_filtering/07_genomic_context/igv_bams/logs/IGV_subset_%A_%a.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

BAM_DIR="/project/bishalab/msarker/c.auris_new/STAR_alignment"
REGIONS="/project/bishalab/msarker/c.auris_new/lncrna_filtering/07_genomic_context/igv_review_regions.bed"
OUT_DIR="/project/bishalab/msarker/c.auris_new/lncrna_filtering/07_genomic_context/igv_bams"

mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*_Aligned.sortedByCoord.out.bam" | sort)

if [[ "${#BAMS[@]}" -ne 24 ]]; then
    echo "ERROR: Expected 24 BAM files but found ${#BAMS[@]}"
    exit 1
fi

BAM="${BAMS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
OUTPUT_BAM="${OUT_DIR}/${SAMPLE}_IGV_regions.bam"

echo "Sample: $SAMPLE"
echo "Input BAM: $BAM"
echo "Output BAM: $OUTPUT_BAM"

samtools view -@ "$SLURM_CPUS_PER_TASK" -bh -L "$REGIONS" "$BAM" -o "$OUTPUT_BAM"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$OUTPUT_BAM"
samtools quickcheck -v "$OUTPUT_BAM"

echo "Completed: $SAMPLE"