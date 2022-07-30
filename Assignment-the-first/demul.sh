#!/bin/bash
#SBATCH --account=bgmp         ### SLURM account which will be charged for the job
#SBATCH --job-name=Demultiplex    ### Job Name
#SBATCH --output=Demultiplex_%j.out         ### File in which to store job output
#SBATCH --error=Demultiplex-%j.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --cpus-per-task=16   ### Number of cpus (cores) per task
#SBATCH --partition=bgmp          ### partition to run things

conda activate base
#base is the environment you are on

# ./demultiplex_1.py -l 101 -r 363246735
# ./demultiplex_2.py -l 8 -r 363246735
# ./demultiplex_3.py -l 8 -r 363246735
./demultiplex_4.py -l 101 -r 363246735
