#!/bin/bash
#SBATCH --job-name=qualimap
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=100Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load qualimap/2.3
mkdir -p ../analysis/qualimap/$1

bam=../analysis/star/demultiplex/${1}.bam
qualimap rnaseq \
-outdir ../analysis/qualimap/$1 \
-a proportional \
-bam $bam \
-gtf ~/patralab/rbator01/perlis_lab/drugseq_mar24/ref/Homo_sapiens.GRCh38.108.gtf \
--java-mem-size=8G

