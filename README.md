# Severity Stratifed Analysis of Pre-chemotherapy TB Studies

This directory contains the code and data to run analyses and produce results for 
Natural history of tuberculosis disease according to disease severity at presentation
by Rodriguez, CA and Leavitt SV, et al. This paper describes a meta-analysis of disease
prognosis (including survival and self-cure rates) of pre-chemotherapy tuberculosis
studies stratified by disease severity.

## Data

### cure_data.csv

This is the individual-level self-cure data for studies used in CITE PAPER 1

### cure_data_all.csv

LAURA REMOVE OR DESCRIBE

### cureDataSummary.csv

LAURA REMOVE OR DESCRIBE (if keeping, change name to cure_data_summary.csv to be consistent)

### mortality_data.csv

This is the individual-level mortality data used in CITE PAPER 1

### NaturalRecovaryData.csv

LAURA REMOVE OR DESCRIBE (if keeping, change name to natural_recovery_data.csv to be consistent)

### study_id.csv

This table details the concordance between the numeric study IDs and the papers they 
refer to (first author and year) as well as various characteristics of the studies.




***

## Scripts

### cure_analysis.R

This script runs the Bayesian logistic regression model for self-cure and extracts the 
data from the results for the table used in the manuscript.

### mortality_functions.R

This script contains the functions to run all of the Bayesian mortality survival
models: TB-specific mortality both with (stratified model) and without (complete model)
a fixed effect for disease severity.

### plot_all_cure_data.R

LAURA REMOVE OR DESCRIBE

### plotAllCureData.html

LAURA REMOVE OR DESCRIBE (if keeping, change name to plot_all_cure_data.html for consistency).

### plotAllCureData.Rmd

LAURA REMOVE OR DESCRIBE (if keeping, change name to plot_all_cure_data.Rmd for consistency).

### Ragonnet calculations.R

LAURA REMOVE OR DESCRIBE (if keeping, change name to ragonnet_calculations.R for consitency)

### severity_analysis.R

This script runs and saves the results of the Bayesian mortality TB-survival analysis
for the complete and stratified model and saves the output in an R workspace
(bayesian_mortality.RMD) 

### severity_results.R

This script takes the results of the analyses and formats them, creating tables and 
figures for the main text and supplement. 

### utils.R

This script contains many functions called by the other programs to format the data
for analysis and results for presentation.
