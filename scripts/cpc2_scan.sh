#!/bin/bash
#SBATCH --job-name=CPC2_scan
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/project/bishalab/msarker/c.auris_new/lncrna_filtering/05_cpc2/logs/CPC2_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/lncrna_filtering/05_cpc2/logs/CPC2_%j.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_lnc_v2

INPUT=/project/bishalab/msarker/c.auris_new/lncrna_filtering/04_rfam/results/primary_ux_rfam_pass.fa
OUTBASE=/project/bishalab/msarker/c.auris_new/lncrna_filtering/05_cpc2/results/cpc2_results

CPC2_BIN=$(command -v CPC2.py)

if [ -z "$CPC2_BIN" ]; then
    echo "ERROR: CPC2.py not found"
    exit 1
fi

"$CPC2_BIN" -i "$INPUT" -o "$OUTBASE"