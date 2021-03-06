instal_chk=0
if(Sys.info()['sysname']=='Windows'){
  Sys.setlocale(category = 'LC_ALL', 'English')
  plat=version["platform"]
  plat=substring(plat,1)
  ver=version["version.string"]
  ver=substring(ver,11,15)
  lib_dir=paste0("~/.checkpoint/2018-04-29/lib/",plat)
  lib_dir=paste0(lib_dir,"/")
  lib_dir=paste0(lib_dir,ver)
  packg_rbg=paste0(lib_dir,"/RBGL")
  packg_graph=paste0(lib_dir,"/graph")
  if(!file.exists(packg_rbg) || !file.exists(packg_graph)){
	instal_chk=1
    }
}else{
  if(!file.exists("~/.checkpoint/2018-04-29")){
	instal_chk<<-1
  }
}
checkpoint::checkpoint("2018-04-29")
if(instal_chk==1){
	source("https://bioconductor.org/biocLite.R")
	biocLite("RBGL")
	biocLite("graph")
}
#Libararies 
library(shiny)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(corpcor)
library(space)
library(pcalg)
library(parallel)
library(ParallelPC)
library(DT)

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

#seed number
seed_num=NULL



