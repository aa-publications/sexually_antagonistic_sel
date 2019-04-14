#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=06:10:00
#SBATCH --mem=18G
#SBATCH --output=assoc_%J
#SBATCH --job-name=assoc

module load Anaconda3
source activate py36_r_ml

# Pass your Python script to the Anaconda3 python intepreter for execution
python test_for_assoc.py
