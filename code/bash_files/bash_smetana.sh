#!/bin/bash
#SBATCH --job-name=smetana_scores
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150gb
#SBATCH --time=24:00:00 
#SBATCH --account=ag-toepfer
#SBATCH --output=/home/freichha/smetana/%j/slurm-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=freichha@smail.uni-koeln.de

module load lang/Miniconda3/23.9.0-0
module load math/CPLEX/22.1.1 
conda activate smetana

set -e
working_dir = "/home/freichha/smetana"
set working dir 
MODELS_DIR = "/home/freichha/MICOM/inputs/syncom_models"

smetana /MODELS_SC1/*.xml -m AF_complete --mediadb af7_complete_smetanaformat.tsv -c community_list_smetana.tsv -o MODELS_SC1/smetanaout_sc1.csv
