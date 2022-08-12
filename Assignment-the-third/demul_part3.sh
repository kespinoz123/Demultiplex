#!/bin/bash
#SBATCH --account=bgmp         ### SLURM account which will be charged for the job
#SBATCH --job-name=Demultiplex_P3   ### Job Name
#SBATCH --output=Demultiplex_%j.out         ### File in which to store job output
#SBATCH --error=Demultiplex-%j.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --cpus-per-task=1  ### Number of cpus (cores) per task
#SBATCH --partition=bgmp          ### partition to run things

R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" 
R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" 
R3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" 
R4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" 
Index="/projects/bgmp/shared/2017_sequencing/indexes.txt"

conda activate base

#base is the environment you are on
./Demul_Parse.py -R1 $R1 -R2 $R2 -R3 $R3 -R4 $R4 -i $Index 


