#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=64Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

export PATH="/cluster/tufts/patralab/rbator01/software/STAR-2.7.11b/source/:$PATH"
STAR --runMode genomeGenerate \
--genomeDir  /cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/STAR/ \
--genomeFastaFiles  /cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/Homo_sapiens.GRCh38.108.gtf \
--runThreadN 8
