#!/bin/bash
#SBATCH --job-name=BAM_index
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/project/bishalab/msarker/c.auris_new
#SBATCH --output=/project/bishalab/msarker/c.auris_new/STAR_alignment/BAM_index_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/STAR_alignment/BAM_index_%j.err
#SBATCH --mail-user=msarker@uwyo.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
shopt -s nullglob

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

BAM_DIR="/project/bishalab/msarker/c.auris_new/STAR_alignment"
BAMS=("$BAM_DIR"/*Aligned.sortedByCoord.out.bam)

if [[ "${#BAMS[@]}" -ne 24 ]]; then
    echo "ERROR: Expected 24 BAM files, but found ${#BAMS[@]}"
    exit 1
fi

for bam in "${BAMS[@]}"; do
    echo "Indexing: $(basename "$bam")"
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$bam"
done

echo "BAM indexing completed: $(date)"