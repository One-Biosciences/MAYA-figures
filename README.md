# Multi-modal quantification of pathway activity with MAYA 

Scripts to analyze single-cell datasets with MAYA, a computational method that enables the detection and scoring of the diverse modes of activation of biological pathways across cell populations. 

## 1- Checking package dependencies

The first step is to check that all required packages are installed. The list is found at the top of the ./scripts/utils.R file.

## 2- Retrieving/Generating datasets

The second step is to generate all datasets required for the analysis. 
You can find the code to generate all datasets in ./scripts/generate_datasets.R
Some require to download data where indicated and store them in the indicated folder.
Some datasets are already provided (Kidney, Colon, Larynx) as SingleCellExperiment objects in corresponding folders in ./datasets/.

## 3- Reproducing analysis and figures

All other scripts in the directory ./scripts can then be run in the indicated order. 
Scripts 1 to 3 reproduce all analyses and generate all files necessary to generate the plots, that are then used in script 4 to 8 to reproduce the figures from the manuscript.



### Manuscript

Landais, Y. & Vallot, C. Multi-modal quantification of pathway activity with MAYA. bioRxiv (2022) doi:10.1101/2022.07.19.500633