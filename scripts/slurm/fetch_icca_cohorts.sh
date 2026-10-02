#!/bin/bash
# Download and inventory the iCCA cohorts (config/icca_cohorts.yaml).
# Needs internet: run on a login/transfer node from the project root, not via sbatch on a compute node.
#   bash scripts/slurm/fetch_icca_cohorts.sh                # all cohorts
#   bash scripts/slurm/fetch_icca_cohorts.sh --dry-run
#   bash scripts/slurm/fetch_icca_cohorts.sh --only GSE244807 TCGA_CHOL
# Long downloads: wrap in tmux/screen or nohup.
set -euo pipefail
source "${CONDA_BASE:-/users/sturkars/mambaforge}/etc/profile.d/conda.sh"
export PYTHONUTF8=1
conda activate "$PWD/envs/hcc-prep"
python scripts/tools/fetch_icca_cohorts.py "$@"
