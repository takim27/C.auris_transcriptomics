#!/bin/bash
#SBATCH --job-name=GFFCompare
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --chdir=/project/bishalab/msarker/c.auris_new
#SBATCH --output=/project/bishalab/msarker/c.auris_new/gffcompare/GFFCompare_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/gffcompare/GFFCompare_%j.err
#SBATCH --mail-user=msarker@uwyo.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_core_v2

REFERENCE_GTF="/project/bishalab/msarker/c.auris_new/gencode.v48.annotation.gtf"
MERGED_GTF="/project/bishalab/msarker/c.auris_new/stringtie_merge/stringtie_merged.gtf"
OUTPUT_PREFIX="/project/bishalab/msarker/c.auris_new/gffcompare/gffcmp"

gffcompare -V -r "$REFERENCE_GTF" -o "$OUTPUT_PREFIX" "$MERGED_GTF"

echo "GFFCompare completed: $(date)"