#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=100Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load subread/2.0.1

#genomeDir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/STAR/
#barcode_file=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/barcodes.txt
#bamdir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/analysis/star/
#R2=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/Unaligned_20240312_Joshua_Bowen-8055/240308_Drug-seq_S1_L002_R2_001.fastq.gz
#R1=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/Unaligned_20240312_Joshua_Bowen-8055/240308_Drug-seq_S1_L002_R1_001.fastq.gz

gtf=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/Homo_sapiens.GRCh38.108.gtf
output_name=read_counts.txt
bam_dir=~/patralab/rbator01/perlis_lab/drugseq_mar24/analysis/star/Aligned.sortedByCoord.out.bam

featureCounts -t exon -g gene_id -s 1 -a $gtf -o $output_name
 
