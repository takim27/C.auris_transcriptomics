#!/bin/bash
#SBATCH --job-name=cutadapt_QS
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --array=0-23%6
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --chdir=/project/bishalab/msarker/c.auris_new/cutadapt
#SBATCH --mail-user=msarker@uwyo.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=cutadapt_%A_%a.out
#SBATCH --error=cutadapt_%A_%a.err

set -euo pipefail
shopt -s nullglob

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

mkdir -p trimmed cutadapt_logs

FASTQS=(SRR*.fastq.gz)

if [[ "${#FASTQS[@]}" -ne 24 ]]; then
    echo "ERROR: Expected 24 FASTQ files but found ${#FASTQS[@]}"
    exit 1
fi

R1="${FASTQS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE="${R1%.fastq.gz}"
OUTPUT="trimmed/${SAMPLE}.trimmed.fastq.gz"

echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Sample: $SAMPLE"
echo "Input: $R1"
echo "Output: $OUTPUT"
echo "Start: $(date)"

cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 "$R1" 2> "cutadapt_logs/${SAMPLE}.stage1.log" | cutadapt -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.1" - 2> "cutadapt_logs/${SAMPLE}.stage2.log" | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o "$OUTPUT" - 2> "cutadapt_logs/${SAMPLE}.stage3.log"

echo "Finished: $SAMPLE"
echo "End: $(date)"