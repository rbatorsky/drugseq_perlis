#!/bin/bash -l
#SBATCH -J wgcna
#SBATCH --time=7-0:00:00 
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=300Gb
#SBATCH --partition=patralab,largemem,batch,preempt
#SBATCH --output=err/sample_%j.out #saving standard output to file
#SBATCH --error=err/sample_%j.err #saving standard error to file
 
module purge
export SINGULARITY_BIND="/cluster/tufts"

/cluster/tufts/biocontainers/tools/r-bioinformatics/4.4.0/bin/Rscript --no-save wgcna.R

