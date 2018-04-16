# Shiny_DMirNet
Shiny-based Web application for exploring direct microRNA-mRNA associations from expression profiles
## Getting Started
These instructions will get you a copy of the web application up and running on your local machine for development and testing purposes.
### Prerequisites
The following application must be installed to run DMirNet.
* [R version 3.4.3](https://cran.r-project.org/) 
* [R studio version 1.1.419](https://www.rstudio.com/products/rstudio/download/)
### Installing
The following R packages must be install to run Shiny_DMirNet. Run the following scripts to install the packages.
```R
install.packages("shiny")
install.packages("checkpoint")
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")
```
## Running Shiny_DMirNet
There are two ways of installing DMirNet
* By running Script in R studio: To run DMirNet locally from GitHub, run the following script
```R
shiny::runGitHub('Shiny_DMirNet','dmirnet')
```
* By downloading the source code: Download the zip file of Shiny_DMirNet and extract the source file. Then run global.R, ui.R, or server. R source files. 
