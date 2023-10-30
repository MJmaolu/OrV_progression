# OrV progression

Code generated for analyze the Dynamics of Orsay virus accumulation and the _C. elegans_ responses within the course of an infection

SCRIPTS:
---

A) PREPROCESS AND ALIGN READS (preDEA)

1) Preprocess fastq files with `bbduk.sh`: `preprocessOrV.sh`
2) Align and quantify gene expression with `STAR`: `alignN2.sh`

B) Differential Expression Analysis and related
1) Create count matrix from individual STAR count files: `prepareCountMatrix.R`
2) DEA: `DEA.R`
3) Correlated/anticorrelated genes: `correlationLog2FCwVirus.R`
4) Comparison analysis: `comparisonAnalysis.R`

C) VISUALIZATION: code to create the plots
- Main figures: `fig1.R`, `fig2.R`, `fig3.R`, `fig4.R`, `fig5.R` and `fig6.R`
