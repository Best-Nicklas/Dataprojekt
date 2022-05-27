# RyouSick <img src="man/figures/logo.png" align="right" width="120"/>
<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-brightgreen)](https://opensource.org/licenses/MIT/)
[![reposize](https://img.shields.io/github/repo-size/Best-Nicklas/Dataprojekt)](https://github.com/Best-Nicklas/Dataprojekt)

<!-- badges: end -->

## Overview
RyouSick is an R package designed to simulate both small and large data sets containing genotype and phenotype information on individuals (read more in `vignette("Simulation_Basis")`). It furthermore provides tools to analyse such datasets to find out which genetic compositions might be causing various inherited illnesses (read more in the articles `vignette("GWAS")`, `vignette("GWAX")` and `vignette("LTFH")`). 

As such RyouSick provides a dual functionality. Firstly, users will be able input real genotype and phenotype data, provided that it is correctly formatted and loaded (read more in `vignette("Writing_efficient_R_code_and_working_with_large_datasets")`), and analyze it using one of the provide analysis methods to find out which SNPs might be causing the illness. The output of the analysis can then be passed on to the prediction function provided in the package (read more in  `vignette("Prediction")`), which will generate a model than can be used to predict the disease status of future genotype cases.    
Secondly, users can simulate genetic data using the in-package functions and perform their own analyses on it and then compare their performance with more established analysis methods that are provided in the package. 

Reading the articles in the presented order will present a clear picture of the subject matter, but if you are just dying to begin using the package head over to the "get started" page to dive right into the package and its functionality. 

## Installation

You can install the development version of RyouSick from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools") if devtools not already installed
devtools::install_github("Best-Nicklas/Dataprojekt")
```

## Dependencies
To  store large quantities of data outside the RAM and perform calculations on this data RyouSick uses the bigsnpr and bigstatr packages. Additionally RyouSick allows speedup of simulation by the use of paralization which is underpinned by the future.apply package and the future backend. Besides these a number of tidyverse packages, such tibble, dplyr and purrr are used in the package and will also be installed along with RyouSick.



