# 20220531_Oscars_tool
## General information
This repository contains scripts used to prepare and deconvolve cell type proportions from several public bulk RNAseq datasets involving ALS.

## BisqueRNA.R
The main script written by Oscar Bedoya Reina to perform deconvolution of bulk RNAseq data based on a single cell RNAseq dataset. The package used to perform deconvolution here is BisqueRNA.

To run the script two libraries are needed:
- BisqueRNA
- Biobase

## Data files
Data used as input for deconvolution were extracted from the following repositories:
- https://github.com/NathanSkene/ALS_Human_EWCE
- https://github.com/NathanSkene/ALS_Mouse_EWCE

Single cell as well as bulk RNAseq data was extracted by running the scripts inside the repositories and perform additional filtering.  
For more details see the documentation:  
<i>Deconvolution of bulk RNAseq – PRJNA512012.docx</i> - Available upon request</i>
## Resources
The folder called resources contains several resources files used for correct conversion between Ensembl IDs and gene names. Also, it contains a file listing one-to-one orthologs between human and mouse.

## READOUT
Folder containing deconvolved cell proportions of bulk RNAseq data.

## Test_files
The folder called test_files contains pseudo-data to test the BisqueRNA script.

## Contact
For questions don't hesitate to send me a message at:
q.c.lin@students.uu.nl