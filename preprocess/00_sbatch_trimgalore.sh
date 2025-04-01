#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=100Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load trim-galore/0.6.10
#trim_galore --paired ../20240530_Joshua_Bowen/Unaligned_20240530_0060X_Joshua_Bowen/JB*R*.fastq.gz
#trim_galore --paired ../20240312_Joshua_Bowen/Unaligned_20240312_Joshua_Bowen-8055/240308*R*.fastq.gz
trim_galore --paired ../../241126-0144_Joshua_Bowen-8326/241017*R*.fastq.gz

