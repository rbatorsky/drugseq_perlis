#!/bin/bash
#SBATCH --job-name=star
#SBATCH --partition=patralab,preempt,largemem,batch
#SBATCH -N 1		
#SBATCH -n 10
#SBATCH --mem=100Gb
#SBATCH --time=0-24:00:00
#SBATCH --output=star_%j.out
#SBATCH --error=star_%j.err

export PATH="/cluster/tufts/patralab/rbator01/software/STAR-2.7.11b/source/:$PATH"

#rd1
genomeDir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/STAR/
barcode_file=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/barcodes.txt
bamdir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/star/align_trim/
R1=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/trim/240308_Drug-seq_S1_L002_R1_001_val_1.fq.gz
R2=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/trim/240308_Drug-seq_S1_L002_R2_001_val_2.fq.gz
#rd2
genomeDir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/STAR/
barcode_file=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/rd2_barcodefile.txt
bamdir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/star_rd2/align_trim/
R1=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/trim/JB01_S1_L001_R1_001_val_1.fq.gz
R2=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/trim/JB01_S1_L001_R2_001_val_2.fq.gz
#rd3
genomeDir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/ref/STAR/
barcode_file=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/rd3_barcodefile.txt
R1=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/trim/JB02_S2_L002_R1_001_val_1.fq.gz
R2=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/trim/JB02_S2_L002_R2_001_val_2.fq.gz


mkdir $bam_dir

STAR --runMode alignReads \
  --outSAMmapqUnique 60 \
  --runThreadN 8 \
  --outSAMunmapped Within \
  --soloStrand Forward \
  --quantMode GeneCounts \
  --outBAMsortingThreadN 8 \
  --genomeDir $genomeDir \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 14 \
  --soloUMIstart 15 \
  --soloUMIlen 14 \
  --soloUMIdedup NoDedup 1MM_All \
  --soloCellFilter None \
  --soloCBwhitelist $barcode_file \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --outFilterMultimapNmax 1 \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix $bamdir \
  --readFilesIn $R2 $R1
