# FLSelex
*Analyzing impacts of alternative selectivity patterns in FLR* 


### Author: Henning Winker (EC-JRC) 

# Installation
Installing ss3diags requires the librabry(devtools), which can be install by 'install.packages('devtools')' and a R version >= 3.5. `FLSRTBMbeta` also requires the latest version of `library(FLCore)` and suggests using `library(ggplotFL)` for plotting. All can be installed from github.

`devtools::install_github("flr/FLCore")`

`devtools::install_github("flr/ggplotFL")`

`devtools::install_github("henning-winker/FLSRTMBbeta")`

`library(FLSRTMBbeta)`

Compiling C++ in windows can be troublesome. As an alternative to installing from github, a windows package binary zip file can be downloaded [here](https://github.com/Henning-Winker/FLSRTMBbeta/tree/main/BinaryPackage/win).

# User Handbook

Documentation of the available applications and plotting functions is provided in the preliminary [`FLSelex` Handbook](https://github.com/Henning-Winker/FLSelex/blob/main/vignette/FLSelex_handbook.pdf)

# Licence

European Commission Joint Research Centre D.02. Released under the EUPL 1.1.
