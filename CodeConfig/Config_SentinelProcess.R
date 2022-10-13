print("##SET UP PATH for LABRIRIES ")
mainLib<-"./lib/"
# mainLib<-"F:/PROGETTI_AGRONOMICO/LAMPvm/R_IrrWatB_Process/lib/"
print(mainLib)
###############################################################################
print("##LOAD LABRIRIES from CLOUD.........")
# Load Libraries
list.of.packages <- c("raster", "rgdal", "stringr", "lubridate", "ggplot2",   "scales",  "sp",  "zoo", "reshape2", "devtools", "plyr","hsdar", "dplyr", "data.table", "randomForest", "C50", "caTools", "MASS", "data.table", "raster", "rlang", "parallel", "httr", "XML", "doParallel", "geosphere", "rgeos")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos="https://cran.stat.unipd.it/")
lapply(list.of.packages, library, character.only = TRUE)
# print("##/*/*/*/*/set memory usage:.....")
## MEMORY USAGE
## #linux
# devtools::install_github("krlmlr/ulimit", force = TRUE, dependencies = TRUE) # set memeory usage in linux
# set_config(config( ssl_verifypeer = 0L ))
# ulimit::memory_limit(1000000)
## #windows
# memory.limit(size = 10000000) # limit per windows
###############################################################################
print("##LOAD LABRIRIES from LOCAL.........")
# Load Libraries
# Rpackage_LocalLIB_path<-"//data/LAMPvm/R_IrrWatB_Process/lib/"
# localRSTOOLBOX<-paste0(mainLib, 'RStoolbox_0.2.6.zip')
# localRSTOOLBOX<-paste0(Rpackage_LocalLIB_path, 'RStoolbox_0.2.6.zip')
localRSTOOLBOX<-paste0(mainLib, 'RStoolbox_0.2.6.tar.gz')
install.packages(localRSTOOLBOX, repos = NULL, type = "source", dependencies = TRUE)
library(RStoolbox)
#source(localRSTOOLBOX)
###############################################################################
print("##LOAD PRIVATE LABRIRIES from LOCAL.........")
source(paste0(mainLib, 'R_RasterStackTimeSeriesPLOT.R'))## carica libreria sviluppata
source(paste0(mainLib, 'R_PostPROCESS_ML.R'))## carica libreria sviluppata
###################################################
# ###################################################
# SET TEMP DIR for RASTER and LARG FILE
print("##SET TEMP FOLDER.........")
tempFolder <- "./DataFOLDER/tempProcess"
# tempFolder <- "F:/PROGETTI_AGRONOMICO/temp"
print(tempFolder)
rasterOptions(tmpdir = tempFolder)
writeTempPath<-paste ("R_USER = ",  tempFolder, sep = "")
write(writeTempPath, file=file.path(Sys.getenv('R_USER'), '.Renviron'))
