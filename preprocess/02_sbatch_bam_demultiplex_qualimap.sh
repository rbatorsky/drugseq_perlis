#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH --partition=preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --exclude=s1cmp006,s1cmp007,p1cmp072
#SBATCH --mem=64Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

bam_dir=$1
input_bam=$2
sample_id=$3
tag_value=$4

echo "------- making bam -------"
echo $bam_dir
module load picard/2.26.10
picard FilterSamReads I=${bam_dir}/${input_bam}  O=${bam_dir}/demultiplex/${sample_id}.bam TAG=CR TAG_VALUE=${tag_value} FILTER=includeTagValues 

module load qualimap/2.3
mkdir -p ${bam_dir}/qualimap/$sample_id
 
output_bam=${bam_dir}/demultiplex/${sample_id}.bam

ls $output_bam

echo "------- running qualimap -------"

qualimap rnaseq \
-outdir ${bam_dir}/qualimap/${sample_id} \
-a proportional \
-bam $output_bam \
-gtf ~/patralab/rbator01/perlis_lab/drugseq_mar24/ref/Homo_sapiens.GRCh38.108.gtf \
--java-mem-size=8G
