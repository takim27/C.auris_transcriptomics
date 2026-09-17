#!/bin/bash
#SBATCH --job-name=STARidx_v48chr
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --chdir=/project/bishalab/msarker/c.auris_new
#SBATCH --mail-user=msarker@uwyo.edu
#SBATCH --mail-type=ALL
#SBATCH --output=STAR_index_%j.out
#SBATCH --error=STAR_index_%j.err

set -euo pipefail

echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_JOB_NAME: $SLURM_JOB_NAME"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date +'%D %T')"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

FASTA="GRCh38.gencode_v48.chromosomes.fa"
GTF="gencode.v48.annotation.gtf"
INDEX_DIR="STAR_index"

if [[ ! -s "$FASTA" ]]; then
    echo "ERROR: FASTA not found or empty: $FASTA"
    exit 1
fi

if [[ ! -s "$GTF" ]]; then
    echo "ERROR: GTF not found or empty: $GTF"
    exit 1
fi

mkdir -p "$INDEX_DIR"

echo "STAR version:"
STAR --version

STAR --runThreadN "$SLURM_CPUS_PER_TASK" --runMode genomeGenerate --genomeDir "$INDEX_DIR" --genomeFastaFiles "$FASTA" --sjdbGTFfile "$GTF" --sjdbOverhang 74

echo "STAR index generated in: $INDEX_DIR"
echo "End time: $(date +'%D %T')"