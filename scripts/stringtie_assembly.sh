#!/bin/bash
#SBATCH --job-name=StringTie
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-23%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --chdir=/project/bishalab/msarker/c.auris_new
#SBATCH --output=/project/bishalab/msarker/c.auris_new/stringtie_assembly/logs/StringTie_%A_%a.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/stringtie_assembly/logs/StringTie_%A_%a.err
#SBATCH --mail-user=msarker@uwyo.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
shopt -s nullglob

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

BAM_DIR="/project/bishalab/msarker/c.auris_new/STAR_alignment"
GTF_DIR="/project/bishalab/msarker/c.auris_new/stringtie_assembly/gtf"
ABUND_DIR="/project/bishalab/msarker/c.auris_new/stringtie_assembly/abundance"
REFERENCE_GTF="/project/bishalab/msarker/c.auris_new/gencode.v48.annotation.gtf"

BAMS=("$BAM_DIR"/*Aligned.sortedByCoord.out.bam)

if [[ "${#BAMS[@]}" -ne 24 ]]; then
    echo "ERROR: Expected 24 BAM files, but found ${#BAMS[@]}"
    exit 1
fi

BAM="${BAMS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")

echo "Sample: $SAMPLE"
echo "BAM: $BAM"
echo "Start: $(date)"

stringtie "$BAM" \
    -p "$SLURM_CPUS_PER_TASK" \
    --fr \
    -G "$REFERENCE_GTF" \
    -m 200 \
    -o "$GTF_DIR/${SAMPLE}.gtf" \
    -A "$ABUND_DIR/${SAMPLE}.gene_abundance.tsv"

echo "Completed: $SAMPLE"
echo "End: $(date)"