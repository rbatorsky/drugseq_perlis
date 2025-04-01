#!/bin/bash -l
#SBATCH -J dmso_2
#SBATCH --time=7-0:00:00 
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=200Gb
#SBATCH --exclude=s1cmp006,s1cmp007
#SBATCH --partition=preempt,patralab,largemem,batch
#SBATCH --output=dmso_2_%j.out #saving standard output to file
#SBATCH --error=dmso_2_%j.err #saving standard error to file

module purge
export SINGULARITY_BIND="/cluster/tufts"

/cluster/tufts/biocontainers/tools/r-bioinformatics/4.4.0/bin/Rscript --no-save 08_de_dmso_vs_best_dmso_step2.R


