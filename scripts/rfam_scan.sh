#!/bin/bash
#SBATCH --job-name=Rfam_scan
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/project/bishalab/msarker/c.auris_new/lncrna_filtering/04_rfam/logs/Rfam_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/lncrna_filtering/04_rfam/logs/Rfam_%j.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_lnc_v2

INPUT=/project/bishalab/msarker/c.auris_new/lncrna_filtering/03_internal_priming/primary_ux_structural_pass.fa
DBDIR=/project/bishalab/msarker/c.auris_new/lncrna_filtering/04_rfam/database
OUTDIR=/project/bishalab/msarker/c.auris_new/lncrna_filtering/04_rfam/results

Z=$(seqkit stats -T "$INPUT" | awk 'NR==2 {printf "%.6f",$5/1000000}')

echo "Transcript sequence database size: ${Z} Mb"

cmscan \
  --cpu "$SLURM_CPUS_PER_TASK" \
  -Z "$Z" \
  --toponly \
  --cut_ga \
  --rfam \
  --nohmmonly \
  --fmt 2 \
  --clanin "$DBDIR/Rfam.clanin" \
  --tblout "$OUTDIR/rfam.tblout" \
  "$DBDIR/Rfam.cm" \
  "$INPUT" \
  > "$OUTDIR/rfam.cmscan"