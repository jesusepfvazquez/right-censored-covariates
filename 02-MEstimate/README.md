Right censored covariates
================

This GitHub repository contains the code and materials neccesary to 
to reproduce the simulation study, analyses, visualizations, and
tables of the paper ["Establishing the Parallels and Differences Between 
Right-Censored and Missing Covariates"](google.com). The template provided by the Journal of the 
American Statistical Association was used in the creation of this repository.

## Structure of the repository

This repository is divided into separate sections that contain the R-code neccesary to 
generate simulation data, R functions to estimate the parameters of interest, and the 
neccesary code to create the figures. The folders include the following:

1.  01-Data. This folder contains all required data genrating functions to generate simulation data under non-informative and informative covariate censoring.
2.  02-MEstimate. This folder contains all functions used to calculate the paramaters of the linear regression model using the geex package [1].
3.  03-Simulation-Example. This folder contains a simulation example for a sample size of n=300. 

For more details of data simulation procress, please reffer to the Supplementary Material of ["Establishing the Parallels and Differences Between 
Right-Censored and Missing Covariates"](google.com). 

## References

[1] Saul, B.C. and Hudgens, M.G., 2020. The calculus of M-estimation in R with geex. Journal of statistical software, 92(2).
