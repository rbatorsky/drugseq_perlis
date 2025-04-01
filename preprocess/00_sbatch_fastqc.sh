#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=64Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load fastqc/0.11.8

fastqc --outdir /cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/cspr_jan25_analysis/align_trim/ /cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/cspr_jan25_analysis/trim/*val*fq*gz
#fastqc --outdir ../20240530_Joshua_Bowen/analysis/fastqc/ ../20240530_Joshua_Bowen/Unaligned_20240530_0060X_Joshua_Bowen/*.fastq.gz
#fastqc --outdir ../20240312_Joshua_Bowen/analysis/fastqc/ ../20240312_Joshua_Bowen/Unaligned_20240312_Joshua_Bowen-8055/240308_Drug-seq_S1_L002_R2_001.fastq.gz
