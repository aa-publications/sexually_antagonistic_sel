#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=23:00:00
#SBATCH --mem=90G
#SBATCH --output=merge_%j
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL


module load Anaconda3
source activate py36_r_ml

# Pass your Python script to the Anaconda3 python intepreter for execution
python merge_batches.py
