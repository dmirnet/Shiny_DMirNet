if(Sys.info()['sysname']=='Windows'){
  Sys.setlocale(category = 'LC_ALL', 'English')
 }
#Libararies 
library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(corpcor)
library(space)
library(pcalg)
library(rlist)
library(parallel)
library(ParallelPC)
library(DT)
library(checkpoint)

#R Files 
source('ui/directoryInput.R') #source https://github.com/wleepang/shiny-directory-input/
source("server/Common_Functions.R")
source("server/Direct_Correlation_Methods_and_Bootstrapping.R")
source("server/Ensemble_Analysis_Vaidation.R")
source("server/Setup_Working_Directory.R")
source("server/File_read_and_Write.R")
source("ui.R")
source("server.R")


#Global variables 
#working environment
root_dir=getwd()
working_dir=data=NULL
rerun=FALSE
rerun_iter=ensemble_iter=1
Result_Corpcor=NULL

dir_experment=dir_direct_bootstrap=dir_ensemble=dir_analysis=dir_validation=dir_direct_bootstrap_uppertri=NULL

#expermental parameters
cores=iterations=miRNAs=dataset=sample.percentage=NULL 

#description for experment 
readme_view1=readme_view2=readme_view3=readme_view4=NULL

checkpoint::checkpoint("2018-01-01")
