#!/bin/bash
#SBATCH --job-name=MICOM_dropouts_l2pfba
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150gb
#SBATCH --time=24:00:00 
#SBATCH --account=ag-toepfer
#SBATCH --output=/home/freichha/MICOM/%j/slurm-%j.out
#SBATCH --mail-type=END

module load lang/Miniconda3/23.9.0-0
module load math/CPLEX/22.1.1 
conda activate micom
export PYTHONPATH=$PYTHONPATH:$CPLEX_HOME/python/3.10/x86-64_linux #global variable  

MEDIUM_FILE="/home/freichha/MICOM/inputs/dropout_coms/combined_af7_c1i127.csv"
'''DROP_DIR_SC1="/home/freichha/MICOM/inputs/dropout_coms/HvSC1"'''
DROP_DIR_SC2="/home/freichha/MICOM/inputs/dropout_coms/HvSC2"

python3 /home/freichha/MICOM/inputs/dropout_coms/l2pfba_dropouts.py \
    --medium "$MEDIUM_FILE" \
    --dropout_dirs "$DROP_DIR_SC2" #"$DROP_DIR_SC1"
