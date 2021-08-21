# FLSelex
*Analyzing impacts of alternative selectivity patterns in FLR* 


### Author: Henning Winker (EC-JRC) 

# Installation

`FLSelex` requires very recent versions of `FLR` libraries `FLCore`, `FLBRP`, `FLasher` and `ggplotFL`. This can be installed together with `FLSelex` from gihtub using library(devtools):

```{r, eval=FALSE}
installed.packages("devtools")

devtools::install_github("flr/FLCore")

devtools::install_github("flr/FLBRP")

devtools::install_github("flr/FLasher")

devtools::install_github("flr/ggplotFL")

devtools::install_github("henning-winker/FLSelex")

```
However, due to increasing difficulties of compiling C++ code with Rtools for Windows systems, these are also provided a binary package zip files [here](https://github.com/Henning-Winker/FLSelex/tree/main/binary_package/win). Not a dependency, but a very useful to explore selectivity pattern under alternative stock recruitment relationships is the new [`FLSRTMBbeta`](https://github.com/Henning-Winker/FLSRTMBbeta/blob/main/README.md) (Winker and Mosqueira, 2021), for which the latest binary package zip for Windows can be found [here](https://github.com/Henning-Winker/FLSRTMBbeta/tree/main/BinaryPackage/win). 


# User Handbook

Documentation of the available applications and plotting functions is provided in the preliminary [`FLSelex` Handbook](https://github.com/Henning-Winker/FLSelex/blob/main/vignette/FLSelex_handbook.pdf)

# Licence

European Commission Joint Research Centre D.02. Released under the EUPL 1.1.
