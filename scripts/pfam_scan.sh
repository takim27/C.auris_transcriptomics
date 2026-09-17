#!/bin/bash
#SBATCH --job-name=Pfam_scan
#SBATCH --account=bishalab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/project/bishalab/msarker/c.auris_new/lncrna_filtering/06_pfam/logs/Pfam_%j.out
#SBATCH --error=/project/bishalab/msarker/c.auris_new/lncrna_filtering/06_pfam/logs/Pfam_%j.err

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate auris_lnc_v2

DB_GZ=/project/bishalab/msarker/c.auris_new/lncrna_filtering/06_pfam/database/Pfam-A.hmm.gz
DB=/project/bishalab/msarker/c.auris_new/lncrna_filtering/06_pfam/database/Pfam-A.hmm
PEP=/project/bishalab/msarker/c.auris_new/lncrna_filtering/06_pfam/results/primary_ux_longest_orfs.pep.fa
OUT=/project/bishalab/msarker/c.auris_new/lncrna_filtering/06_pfam/results

if [ ! -s "${DB}.h3m" ]; then
    gzip -dc "$DB_GZ" > "$DB"
    hmmpress -f "$DB"
fi

hmmscan --cpu "$SLURM_CPUS_PER_TASK" --cut_ga --domtblout "$OUT/pfam.domtblout" "$DB" "$PEP" > "$OUT/pfam_scan.txt"

echo "Pfam scan completed: $(date)"