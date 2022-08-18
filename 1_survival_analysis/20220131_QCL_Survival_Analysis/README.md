# 20220131_QCL_Survival_Analysis
## General information
This repository contains adjusted scripts and data files used in:  
<i><b>Altered perivascular fibroblast activity precedes ALS disease onset</i></b>

Several adjustment made to the original scripts include:  
- Updating used packages
- Excluding redundant code
- Inclusion of novel data
- Enabling bulk analyses

For the repository with the original scripts see:  
https://github.com/lewandowskilab/PVF_Manuscript

## Data
All data files used in the scripts are located in this these folders:
- generated_data_R
- R_data
(Access upon request)

## 2020_02_24_Cutpoint_Analyses_Publ.R
<i>The script is an adjusted version of the original script 2020_02_24_Cutpoint_Analyses_Publ.R</i>  
As an addition to the script, the tidyverse package was added. Seed was set to '747' to enable reproducibility of the results. Furthermore, a results naming variable was added to the custom plot function to enable bulk analyses of all proteins interested.  

To see if we could see any difference in survival in genotype stratified samples with non-stratified samples, we included a p-val comparison section in the results.

## 2022_02_17_Cox_bulk.R
<i>The script is an adjusted version of the original script 2020_02_24_CoxOutput.R</i>  
As an addition to the script, the tidyverse package was added. Seed was set to '747' to enable reproducibility of the results. Furthermore, a results naming variable was added to the custom plot function to enable bulk analyses of all proteins interested.

Multivariate Cox Proportional Hazard Model functions were rewritten in this new script to enable analysis of genotype interactions with the proteins of interest. All proteins can be analysed all together with the same function.

## Additional notes
To access the encrypted data, a specific package is needed. 

Specific R package:
- SBA 1.37 (Suspension Bead Array, Version 1.37)

Additional information on the scripts can be found in the original script depository.

## Contact
For questions don't hesitate to send me a message at:
q.c.lin@students.uu.nl