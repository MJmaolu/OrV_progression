#!/bin/bash 
#preprocessOrV.sh
#SBATCH --job-name=preprocessOrV 
#SBATCH --partition=short
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=12 
#SBATCH --mem=10gb 
#SBATCH --time=1-00:00:00 
#SBATCH --output=preprocessOrV_%j.log 

<< "preprocessOrV.sh"
Procesado con cada muestra:
- Clean fq.gz with bbduck.sh
 
MJ OLMO-UCEDA
2022/06/06
preprocessOrV.sh

# rutas
samples=$1

dir="/storage/evsysvir/OrvProgressionInfection"
REF="/storage/evsysvir/REFERENCES"
adapters=$REF"/novogen_adapters_orvProgression.fasta"
ref_vir=$REF"/OrV_lab/OrV_complete_genome_mut_evolsysvir.fasta"

dir_fq=$dir"/raw_data"
dir_clean_fq=$dir"/clean_data"
dir_discarded=$dir"/discarded_data"


dir_qc1=$dir"/QC1"

# parámetros 
minlength=60
threads=20

start=`date +%s`

# PROCESO POR CADA MUESTRA DEL FICHERO $SAMPLES
while read sample
do
    echo "****************** PROCESANDO MUESTRA $sample ******************"
    echo

    # 1. Preprocesado de los fastqs con bbduck.SH
    echo "········································"
    echo "@--> Preprocesando los raw_fastqs"
    echo "········································"

    bbduk.sh \
    in=${dir_fq}/${sample}_1.fq.gz \
    in2=${dir_fq}/${sample}_2.fq.gz \
    out=${dir_clean_fq}/${sample}_clean_1.fq.gz \
    out2=${dir_clean_fq}/${sample}_clean_2.fq.gz \
    outm=${dir_discarded}/${sample}_dis.fq \
    ref=${adapters} \
    ktrim=r k=31 mink=11 \
    qtrim=rl trimq=10 maq=5 minlength=$minlength \
    forcetrimleft=10
    echo

done < $samples

# Tiempo total de ejecución
end=`date +%s`
runtime=$((end-start))
echo "Total time: $runtime s"
