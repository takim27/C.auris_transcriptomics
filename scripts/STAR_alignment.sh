#!/bin/bash
#SBATCH --job-name=STAR_align
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-23%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --chdir=/project/bishalab/msarker/c.auris_new
#SBATCH --output=/project/bishalab/msarker/c.auris_new/STAR_alignment/STAR_%A_%a.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/STAR_alignment/STAR_%A_%a.err
#SBATCH --mail-user=msarker@uwyo.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
shopt -s nullglob

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

READ_DIR="/project/bishalab/msarker/c.auris_new/cutadapt/trimmed"
INDEX_DIR="/project/bishalab/msarker/c.auris_new/STAR_index"
OUTPUT_DIR="/project/bishalab/msarker/c.auris_new/STAR_alignment"

FASTQS=("$READ_DIR"/SRR*.trimmed.fastq.gz)

if [[ "${#FASTQS[@]}" -ne 24 ]]; then
    echo "ERROR: Expected 24 trimmed FASTQ files, but found ${#FASTQS[@]}"
    exit 1
fi

if [[ ! -s "$INDEX_DIR/Genome" || ! -s "$INDEX_DIR/SA" ]]; then
    echo "ERROR: STAR genome index is missing or incomplete"
    exit 1
fi

FASTQ="${FASTQS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "$FASTQ" ".trimmed.fastq.gz")

echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Sample: $SAMPLE"
echo "Input: $FASTQ"
echo "Start: $(date)"

STAR \
  --runThreadN "$SLURM_CPUS_PER_TASK" \
  --genomeDir "$INDEX_DIR" \
  --readFilesIn "$FASTQ" \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --quantMode GeneCounts

echo "Completed: $SAMPLE"
echo "End: $(date)"