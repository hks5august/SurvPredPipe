# SurvPredPipe
SurvPredPipe: A R-Package for the Pipeline to Predict Survival Probability of Cancer Patients

Follow the Steps to Install package 
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

### SurvPredPipe_Manual.pdf file explains description of each function, variable, installation detail, and workflow of the pipeline with an example.

### SurvPredPipe_workflow.Rmd file explains workflow with an example for user to follow.

### The "R" folder  contains Rscript for all the functions

### The "inst" folder/directory contains two folders: (1) "extdata" subfolder contains example input data files, user can run example and prepare their input files accoding to example files; (2) output_data contains all the output after running the pipeline with example input data.

### The vignettes directory contains Figures, extdata, workflow file.
