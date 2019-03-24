# Shiny_DMirNet
Shiny-based Web application for exploring direct microRNA-mRNA associations from expression profiles
## Getting Started
These instructions will get you a copy of the web application up and running on your local machine for development and testing purposes.
### Prerequisites
The following application must be installed to run DMirNet.
* [R Software](https://cran.r-project.org/) 
* The following R packages must be install to run Shiny_DMirNet. Run the following scripts to install the packages.
```R
install.packages("shiny")
install.packages("checkpoint")
```
## Running Shiny_DMirNet
There are two ways of installing DMirNet
* By running Script in R studio: To run DMirNet locally from GitHub, run the following script
```R
shiny::runGitHub('Shiny_DMirNet','dmirnet')
```
* By downloading the source code: Download the zip file of Shiny_DMirNet and extract the source file. Then run global.R, ui.R, or server. R source files. To run the source file, run the following script
```R
setwd("<path to Shiny_DMirNet root directory>/")
shiny::runApp()
```
> Note: A request for upgrading packages might appear when running the program for the first time. Since upgrading the packages might cause version incompatability issues, please DO NOT upgrade the packages by selecting 'n' as shown below. 
```R
Old packages: 'assertthat', 'BH', 'cli', 'clue', 'colorspace', 'curl',
  'digest', 'DT', 'ggplot2', 'glue', 'gmp', 'htmlwidgets', 'httpuv', 'igraph',
  'jsonlite', 'later', 'lazyeval', 'mime', 'munsell', 'pcalg', 'pillar',
  'pkgconfig', 'R6', 'Rcpp', 'RcppArmadillo', 'rlang', 'robustbase', 'scales',
  'sfsmisc', 'shiny', 'shinythemes', 'stringi', 'stringr', 'tibble', 'utf8',
  'V8', 'xtable', 'yaml', 'zoo'
Update all/some/none? [a/s/n]: n
```
## Citation
DMirNet is provided under a free-of-charge, open-source license (A-GPL3). All we require is that you cite/attribute the following in any work that benefits from this code or application.
### Citing the Web Application
Minsu Lee, Anene Bekuma Terefe, HyungJune Lee, and Sangyoon Oh, "DMirNet: Shiny-based application for exploring direct microRNA-mRNA associations", Bioinformatics (Under Review)
### Citing the DMirNet framework 
Minsu Lee and HyungJune Lee, (2016) "DMirNet: Inferring direct microRNA-mRNA association networks", BMC Systems Biology, 10(suppl5):125. 
