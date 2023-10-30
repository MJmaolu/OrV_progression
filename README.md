# OrV progression

Code generated for analyze the Dynamics of Orsay virus accumulation and the __C. elegans__ responses within the course of an infection process

SCRIPTS:

A) PREPROCESS AND ALIGN READS

1) Preprocess fastq files with `bbduk.sh`: `preprocessOrV.sh`
2) Align and quantify gene expression with `STAR`: `alignN2.sh`

B) Differential Expression Analysis and related
1) Create count matrix from individual STAR count files: `prepareCountMatrix.R`
2) DEA: `DEA.R`
3) 

C) VISUALIZATION
- Fig1.R
- 
