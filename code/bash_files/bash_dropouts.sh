#!/bin/bash
#SBATCH --job-name=MICOM_dropouts
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150gb
#SBATCH --time=24:00:00 
#SBATCH --account=ag-toepfer
#SBATCH --output=/home/freichha/CarveMe/%j/slurm-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=freichha@smail.uni-koeln.de

module load lang/Miniconda3/23.9.0-0
module load math/CPLEX/22.1.1 
conda activate micom

set -e

MEDIUM_FILE="/home/emma/Dokumente/thesis/media_creation/created_media/af_diff_growth/combined_af7_c1i127.csv"
DROP_DIR_SC1="thesis/communities/simulations/dropouts/HvSC1"
DROP_DIR_SC2="thesis/communities/simulations/dropouts/HvSC2"

python3 run_dropouts.py \
    --medium "$MEDIUM_FILE" \
    --dropout_dirs "$DROP_DIR_SC1" "$DROP_DIR_SC2"