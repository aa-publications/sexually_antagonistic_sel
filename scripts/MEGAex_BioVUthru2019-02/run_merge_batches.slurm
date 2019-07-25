#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=18:00:00
#SBATCH --mem=26G
#SBATCH --output=pca_%j
#SBATCH --job-name=no_batch4_pca
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL


module load Anaconda3
source activate py36_r_ml

# Pass your Python script to the Anaconda3 python intepreter for execution
python merge_batches.py
