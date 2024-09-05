Right censored covariates
================

This GitHub repository contains the materials necessary to 
to reproduce the simulation study, analyses, visualizations, and
tables of the paper ["Establishing the Parallels and Differences Between 
Right-Censored and Missing Covariates"](google.com). Additionally, users can use this code to implement the proposed estimators in their real data. The template provided by the Journal of the 
American Statistical Association was used in the creation of this repository.

## Structure of the repository

This repository is organized into several sections, each containing the R code necessary to generate simulation data, estimate parameters, and create the figures presented in the paper. The folders are structured as follows:

1.	  01-Data: Contains all helper functions for simulating data, as well as the simulation data used to create the tables and figures in the paper.
2.	  02-MEstimate: Includes all helper functions for calculating the parameters of the linear regression model using the geex package [1].
3.	  03-Simulation-Example: Contains the Example.Rmd, a tutorial that demonstrates how to replicate all the analyses from the published paper.
4.    docs: This folder includes the necessary files for hosting the Example.Rmd tutorial online [link](https://jesusepfvazquez.github.io/right-censored-covariates/#1_Introduction).

For more details on the data simulation process, refer to [“Establishing the Parallels and Differences Between Right-Censored and Missing Covariates”](google).

## References

[1] Saul, B.C. and Hudgens, M.G., 2020. The calculus of M-estimation in R with geex. Journal of statistical software, 92(2).
