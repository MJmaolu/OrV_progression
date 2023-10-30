#!/bin/bash
#alignN2.sh
#SBATCH --job-name=alignN2.sh 
#SBATCH --partition=short
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=12 
#SBATCH --mem=20gb 
#SBATCH --time=1-00:00:00 
#SBATCH --output=alignN2_%j.log

<< "alignN2"
For each sample in samples.txt:
- Align paired end reads to the complete genome with STAR
- Mark duplicates with gatk
- remove duplicates

At the end of the process, the sorted by coordinate alignments are stored
as ${dir_sorted}/_N2_${sample}_Aligned.sortedByCoord.out.bam

MJ OLMO-UCEDA
2022/06/06
alignN2

samples=$1

dir="/storage/evsysvir/OrvProgressionInfection"
REF="/storage/evsysvir/REFERENCES"
dir_host=$REF"/N2_bristol"

dir_clean_fq=$dir"/clean_data"

dir_alignments="/storage/evsysvir/TimeSeries/PRJNA636173/alignments"

dir_qc1=$dir"/QC1"
dir_alignment=$dir"/alignments"
dir_sorted=$dir_alignment"/sorted"

STAR="/home/maolu/programs/STAR-2.7.9a/bin/Linux_x86_64/STAR"

while read sample
do
    # Align vs complete genome hg19
    $STAR --runThreadN 12 \
    --runMode alignReads \
    --genomeDir $dir_host \
    --readFilesCommand gunzip -c\
    --quantMode GeneCounts \
    --readFilesIn ${dir_clean_fq}/${sample}_clean_1.fq.gz ${dir_clean_fq}/${sample}_clean_2.fq.gz \
    --outFileNamePrefix ${dir_sorted}/${sample}_N2_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes NH HI NM MD AS

done < $samples
