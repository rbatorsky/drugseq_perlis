#!/bin/bash

## rd1
demultiplexed_bam_out_dir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/star/align_trim/
input_bam=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/analysis/star/align_trim/Aligned.sortedByCoord.out.bam
barcode_info=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240312_Joshua_Bowen/sample_barcodes.txt

##rd2
#bam_dir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/star_rd2/align_trim/
#input_bam=Aligned.sortedByCoord.out.bam
#barcode_info=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/rd2_sample_barcodefile.txt

##rd3
#bam_dir=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/analysis/star_rd3/align_trim/
#input_bam=Aligned.sortedByCoord.out.bam
#barcode_info=/cluster/tufts/patralab/rbator01/perlis_lab/drugseq_mar24/20240530_Joshua_Bowen/rd3_sample_barcodefile.txt


mkdir ${bam_dir}/demultiplex/

while IFS=$'\t' read -r -a line;
do
  sample_id=${line[0]}
  tag_value=${line[1]}
  echo $sample_id
  echo $tag_value
  #ls ${demultiplexed_bam_out_dir}/../../qualimap/$sample_id/qualimapReport.html
  sbatch 02_sbatch_bam_demultiplex_qualimap.sh ${bam_dir} ${input_bam} ${sample_id}  ${tag_value}  
  
done < $barcode_info

