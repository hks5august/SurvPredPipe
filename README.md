# SurvPredPipe
The SurvPredPipe package is for a computational pipeline for the prediction Survival Probability of Cancer Patients. It performs various steps: Data Processing, Split data into training and test subset, Data Normalization, Select Significant features based on Univariate survival,  Generate LASSO PI Score, Develop Prediction model for survival probability based on different features  and draw survival curve based on predicted survival probability values and barplots for predicted mean and median survival time of patients.



![SurvPredPipe_workflow](https://github.com/hks5august/SurvPredPipe/assets/46756503/795d1d62-9d8a-487f-88b5-7c76d5fdaf6a)

Figure: The Workflow of SurvPredPipe representing different steps performed by various functions of SurvPredPipe package.


# Follow the Steps to Install SurvPredPipe package on your Local R system:
```r
#Step1: First Install remote package
install.packages("remotes") 
#load remortes package
library("remotes")
#Step2: install SurvPredPipe package
remotes::install_github("hks5august/SurvPredPipe", local = TRUE)
# load package
library("SurvPredPipe") 
```



# SurvPredPipe_Package_Description.pdf:
This file explains description, installation detail, details of each function, variable.

# SurvPredPipe_Package_workflow.pdf: 
This file explains workflow (Step by Step) of SurvPredPipe package with an example for user to follow.

# SurvPredPipe_workflow.Rmd: 
This file explains workflow with an example for user to follow.

# R Folder:
It folder/directory Rscript for all the functions

# "inst" Folder: 
This folder/directory contains two folders: (1) "extdata" subfolder contains example input data files, user can run example and prepare their input files accoding to example files; (2) output_data contains all the output after running the pipeline with example input data.

# The vignettes: 
This directory contains Figures, extdata, workflow file.


