Right censored covariates
================

This GitHub repository contains the materials necessary to 
to reproduce the simulation study, analyses, visualizations, and
tables of the paper ["Establishing the Parallels and Differences Between 
Right-Censored and Missing Covariates"](google.com). Additionally, users can use this code to implement the proposed estimators in their real data. The template provided by the Journal of the 
American Statistical Association was used in the creation of this repository.

## Structure of the repository

This repository is divided into separate sections that contain the R-code necessary to 
generate simulation data, R-functions to estimate the parameters of interest, and the 
necessary code to create the figures presented in the paper. The folders include the following:

1.  01-Data. This folder contains all helper functions used to simulate data, as well as all simulation data used  to create the tables and figures of the published paper. 

2.  02-MEstimate. This folder contains all helper functions used to calculate the parameters of the linear regression model using the geex package [1].

3.  03-Simulation-Example. This folder contains a **Example.Rmd** which is a tutorial of how to replicate all analyses of the published paper. 

For more details of data simulation process, please refer to ["Establishing the Parallels and Differences Between 
Right-Censored and Missing Covariates"](google.com). 

## References

[1] Saul, B.C. and Hudgens, M.G., 2020. The calculus of M-estimation in R with geex. Journal of statistical software, 92(2).
