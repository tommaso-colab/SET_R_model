#! /usr/bin/Rscript
# ##############################################################
# ##############################################################
# ##load libraries and CONFIG file
mainfolder = '//data/LAMPvm/R_IrrWatB_Process'
setwd(mainfolder)
source('./CodeConfig/Config_SentinelProcess.R')
# ## dayRUN
# todayRUNdate <- ('20210504')
todayRUNdate <- format(Sys.time(), "%Y%m%d")
# ##############################################################
# ##Define folder input for RASTER: 1 folder each image
inputFolder<-"//data/LAMPvm/temp"
copyInputToOTHERFolder<-"//data/LAMPvm/tempTQ"
# #
listFolder<-list.files(inputFolder, full = TRUE, pattern = ".SAFE",)
# ##############################################################
# logger run
outputfile <- "./logfile.txt"
sink(file = outputfile, append = TRUE, type = c("output", "message"),
     split = FALSE)
sink.number(type = c("output", "message"))
print(paste0("#/*/*start R process: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
# ##############################################################
# ##PARAMS input for raster process
cutArea<-"./DataFOLDER/PlotCUT_Shape/CutArea/StudyAREA.shp" #shape file dove fare il CUT
CutYorN <- "N" #definisci se eseguire il CUT Y/N
SaveNew<- "N"  # deinifsci se salvare il file o no  Y/N SE FAI CUT deve ESSERE SALVATO
CUTRasterfolderpath<-"./DataFOLDER/rasterCUT/tileTQ" #output folder dove salvare il risultato del CUT
# ##PARAMS input FOR MODEL ML
filenamepathBrickSI<-paste0("./ModelOUT/tempSIbrick_", todayRUNdate, "/") ## FOLDER DOVE VENGONO SALVATI TUTTI I BRICK CON GLI INDICI CALCOLATI che poi saranno RILETTI per il calcolo del modello ML
filenamepathOUTPUT<-paste0("./ModelOUT/MLout_", todayRUNdate, "/") # cartella dove saranno creati i raster 1-0 del modello
dir.create(filenamepathBrickSI, showWarnings = FALSE)
dir.create(filenamepathOUTPUT, showWarnings = FALSE)
modelpath<-"./modelML/model20190302_DVI_GNDVI_MNDWI_NDVI_NDWI2_SATVI.rds" # ptah del modello RandomFOREST da leggere
# ##############################################################
# ##MAIN
GetFileN<-GetRASTERFilePathANDArray(listFolder, cutArea, CUTRasterfolderpath, CutYorN, SaveNew) 
# creo la l'array con un vettore di path dei file per ogni SUBFOLDER
ArrayRaster<-GetFileN$GetRasterArray()# creo un array di raster per ogni folderdss se Y salvo e taglio
# ArrayRaster<-GetFileN$GetRasterArray_onlyREAD()
ArrayRasterDate<-GetFileN$GetRasterDATEvector() # creo una lista di date dei rilievi
# ArrayRasterDate<-GetFileN$GetRasterDATEvectorTIFF() # creo una lista di date dei rilievi da lettura cartella TIFF
ArrayRasterNAMEandDate<-GetFileN$GetRasterNAMEvector() # creo una lista di date dei rilievi
ArrayDateVector<-as.Date(unlist(ArrayRasterDate)) # converto la lista in vettore
# # ##############################################################
# # selezioni INDICI da CALCOLARE
IndexToCalc<-c("DVI", "GNDVI", "MNDWI", "NDVI", "NDWI2" , "SATVI")
# # ###################################################
# # # # RUN PROCESS
# vectorIDraster <- c(1,2,3,4,25,21,8,9,10,11,12)## vettore con ID su BRICK di rasterINPUT con .safe FILE
vectorIDraster <- c(1,2,3,4,5,6,7,8,9,10,11)## vettore con ID su BRICK di rasterINPUT con .safe FILE
# # VECTOR ID with FILENAME	1:B02_10m, 2:B03_10m, 3:B04_10m, 4:B08_10m, 25:SCL_20m , 21:B09_60m, 8:B05_20m,9:B06_20m,10:B07_20m, 11:_B11_20m, 12:B12_20m
# # VECTOR ID with FILENAME	1:B02_10m, 2:B03_10m, 3:B04_10m, 4:B08_10m, 5:SCL_20m , 6:B09_60m, 7:B05_20m,8:B06_20m,9:B07_20m, 10:_B11_20m, 11:B12_20m
# # vectorIDraster <- c(2,5,8,17,25,18,11,13,15,19,21) ## posizione in BRICK su ritagli in RasterCUT
shapefolder<-"./DataFOLDER/PlotCUT_Shape/in" # poligono sui cui estrarre e calcolare le statistiche dei valori
shapeOUTArray<-c()
# shapeOUTArray<-fileListShape(shapefolder)
# shapeOUTArrayFileName<-fileListShapeFileName(shapefolder)
# VECTOR ID with FILENAME	1:B02_10m, 2:B03_10m, 3:B04_10m, 4:B08_10m, 25:SCL_20m , 21:B09_60m, 8:B05_20m,9:B06_20m,10:B07_20m, 11:_B11_20m, 12:B12_20m
IndexValues3<-lapply(ArrayRaster, function(x){
											print("##RUN Cycle##")
                      # browser()
											CalcIndexVector<-CalcIndex(x, shapeOUTArray, IndexToCalc, vectorIDraster, filenamepathBrickSI)
											VIindexes<-CalcIndexVector$RStoolCalcINDEX_NOT_polygon()
											print("################ADD file name to BRICK")
											RasterFileID<-paste(unlist(strsplit(x[[1]]@data@names, "_"))[1], unlist(strsplit(x[[1]]@data@names, "_"))[2], sep="_")#create ID raster from bricks 
											# print(paste0("####VIindexes: ",VIindexes))
											VIindexes@file@name<-as.character(RasterFileID)# add ID filename TO BRICKS
											return(VIindexes)
												})#calcolo gli indici vegetazionali e maschero con SCL raster	


# r <-"/data/LAMPvm/R_IrrWatB_Process/ModelOUT/tempSIbrick_20220512/T32TPQ_20220511T101601.gri"
# b <- brick(r)
# outfile<- "/data/LAMPvm/R_IrrWatB_Process/ModelOUT/temp/hdr_sat_T32TPQ_20220511T101601.tiff"
# writeRaster(b,outfile, format="GTiff", bylayer =TRUE, overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

# plot(dvi)
# plot(dvi@file)


######################################################################
##APPLY MODEL ML to ALL RASTER BRICKS##################################
# input 1: model ML in rds
# input 2: FOLDER PATH per i risultati
# input 3: RASTER BRICK con tutti gli VI che SERVONO ottenuti da CICLO 
# inoput 4: SOGlIE per THRESHOLDING di NVDI e GNDVI
# PROCESS: calcolo del RISULTATO di ML per ogni elemento
# output: file con 0-1 per ogni raster brick di input in FOLDER 
###################################
# # # Threshold for acceptance
# # NDVI min 0.192
# # GNDVI min 0.1945
###################################
print("/*/*/*/*###RUN ML alg:")
gc()
# input
ndviThre<-0.192
gndvithre<-0.1945
# main			
# memory.limit(size=10000000)
# ulimit::memory_limit(10000000)# limit per linux	
# memory.limit(size = 10000000) # limit per windows
VectorFile<-list.files(filenamepathBrickSI, pattern = ".grd", recursive = FALSE, full = FALSE)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

list.of.packages <- c("foreach", "raster", "rgdal", "stringr", "lubridate", "ggplot2",   "scales",  "sp",  "zoo", "reshape2", "devtools", "plyr","hsdar", "dplyr", "data.table", "randomForest", "C50", "caTools", "MASS", "data.table", "raster", "rlang", "parallel", "httr", "XML", "doParallel", "geosphere", "rgeos")

data = foreach(i = 1:length(VectorFile),
                .packages = list.of.packages) %dopar% {
                 try({
                   sink(file = outputfile, append = TRUE, type = c("output", "message"), split = FALSE)
                   print("###START ML:")
                   filein <- VectorFile[i]
                   filenameInput<-paste0(filenamepathBrickSI, filein)
                   print(filenameInput)
                   rasterBrick<-brick(filenameInput)
                   filenamepathRasterML<-paste0(filenamepathOUTPUT, unlist(strsplit(filein, '.g'))[1])
                   print(filenamepathRasterML)
                   modeloutCalc<-PredictRasterWithRMmodelAndSAVE_byROW(rasterBrick, filenamepathRasterML, modelpath, ndviThre, gndvithre)
                   Calc_model<-modeloutCalc$CalculateModelML()
                 })
               }
stopCluster(cl)

# # ###################################################
##RUN POST PROCESS
# print('###RUN POST PROCESS of ML: ')
# IniPostProcessofML_output<-PostProcess_ML_modelFOLDER(filenamepathOUTPUT)
# MergeCreate<-IniPostProcessofML_output$CalculateWGS84_AND_MERGE()
# print('###POST PROCESS of ML COMPLETED.')

# # ###################################################
##CLEAN WORK DIRECTORY
print('###file sentinel move to other folder .tempTQ/temp/')
file.copy(from=inputFolder, to=copyInputToOTHERFolder,
          overwrite = FALSE, recursive = TRUE,
          copy.mode = TRUE)
print('###file sentinel DELETE from working dir .temp/')
unlink(listFolder, recursive=TRUE)
print('###file sentinel DELETE from working dir ./DataFOLDER/tempProcess')
tempFolderR <- "./DataFOLDER/tempProcess/*"
unlink(tempFolderR, recursive=TRUE)
# print('###file BRICK SI DELETE from working dir ./tempBRICK/')
# unlink(filenamepathBrickSI, recursive=TRUE)


################################################################################################################################
################################################################################################################################
################################
##TEST VELOCITA PROCESSORE######
# start 17.38 15/03/2021 - 2 file su GISserver senza multicore
# stop 19.31 15/03/2021

# start 11.39 16/03/2021 - 9 file su GISserver con multicore
# stop 14.29 16/03/2021 

# start 13.36 18/03/2021 - 9 file su CENTO7server con multicore
# stop  15.21 18/03/2021   
#################################
# PARALLEL COMPUTATION in UNIX
# Raster_ArrayML_Parallel<-lapply(VectorFile, function(x) mcparallel({
# 								# browser()
# 								filenameInput<-paste0(filenamepathBrickSI, x)
# 								rasterBrick<-brick(filenameInput)
# 								filenamepathRasterML<-paste0(filenamepathOUTPUT, unlist(strsplit(x, '.g'))[1])
# 								print(filenamepathRasterML)
# 								print(x)
# 								modeloutCalc<-PredictRasterWithRMmodelAndSAVE_byROW(rasterBrick, filenamepathRasterML, modelpath, ndviThre, gndvithre)
# 								Calc_model<-modeloutCalc$CalculateModelML()
# 								})
# 								)
# mccollect(Raster_ArrayML_Parallel)
