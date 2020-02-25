#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:60:00
#SBATCH --mem=24GB
#SBATCH --output=python_job_slurm.out


module load Anaconda3
source activate py36
# Pass your Python script to the Anaconda3 python intepreter for execution
python summarize_alignments.py