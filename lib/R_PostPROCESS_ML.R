PostProcess_ML_modelFOLDER<- function(filenamepathOUTPUT) {
    structure(class = "PostProcess_ML_modelFOLDER", list(
        # attributes
        filenamepathOUTPUT = filenamepathOUTPUT,  
        
        # methods
        CalculateWGS84_AND_MERGE = function() {	
            
            print("##postPROCESS_ML output: ")
            listFolder<-list.files(filenamepathOUTPUT, full = TRUE, pattern = "\\.tif$")
            # a partire dalla cartelal di input del modello ML "MLout_data" creo altra cartella con WGS84 converto il CRS ed unifico file per stessa data
            # creo la cartella di LAVORO 
            FolderMainToModified<-strsplit(filenamepathOUTPUT, "/")[[1]][(length(strsplit(filenamepathOUTPUT, "/")[[1]]))]
            PATHmain<-strsplit(filenamepathOUTPUT, "/")[[1]][0:(length(strsplit(filenamepathOUTPUT, "/")[[1]])-1)]
            pathMain_slash<-paste0(PATHmain, "/")
            pathMain_slashJoin<-paste0(pathMain_slash[1],pathMain_slash[2])
            newFolderToCreate_wgs84<-paste0(pathMain_slashJoin, FolderMainToModified,"_wgs84")
            dir.create(newFolderToCreate_wgs84, showWarnings = FALSE)
            # eseguo conversione di CRS    
            ExportRastesr_wgs84<-lapply(listFolder, function(x){
                
                listTofilename <- as.character(x)
                readRasterF<-raster(listTofilename)
                pr1 <- projectRaster(readRasterF, crs="+proj=longlat +datum=WGS84")
                filname<-(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
                outputFilePath <- paste0(newFolderToCreate_wgs84 ,"/", filname)
                print("##WGS84 filename path: ")
                print(outputFilePath)
                rc <- writeRaster(pr1, filename=outputFilePath,  fprj = TRUE,format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
                writeLines(showWKT(crs(pr1, asText=TRUE)), extension(outputFilePath, 'prj') )
                print("##end SAVE")
                return(listTofilename)
            })
            # faccio ciclo su CARTELLA di WGS84 per LISTA FILE
            listFolder_wgs84<-list.files(newFolderToCreate_wgs84, full = TRUE, pattern = "\\.tif$")
            # creo vettore con unique DATe nel caso per 1 giorno di download abbia più date
            dateFromFile<-lapply(listFolder_wgs84, function(x){
                fullpath<-x
                filname<-(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
                all_dateFromFile<-strsplit(filname, "_")[[1]][2]
                dateFromFile<- substr(all_dateFromFile, 1, 8)
                export<- c( x, dateFromFile)
                return(export)})
            df_FielList<-as.data.frame(do.call(rbind, dateFromFile))
            print(df_FielList)
            # browser()
            datedownloaded<-unique(df_FielList[,2])
            print(datedownloaded)
            # eseguo ciclo per ogni data unica prendo ed unisco in unico FILE
            Selected_RasterbyDATE_AND_PROCESS<-lapply(datedownloaded, function(x){
                SelectedITEM <- df_FielList[df_FielList$V2 == x,1]
                print("/*/*/* select ITEM:")
                # print(x)
                filevector<-c(as.character(unlist(SelectedITEM)))
                print(SelectedITEM)
                outputFILE_path<-paste0(newFolderToCreate_wgs84,"/MERGED_wgs84_" , x ,".tif")
                print(outputFILE_path)
                e <- extent(7, 17, 40, 50)
                template <- raster(e)
                projection(template) <- CRS('+init=epsg:4326')
                writeRaster(template, file=outputFILE_path, format="GTiff")
                mosaic_rasters(gdalfile=filevector,dst_dataset=outputFILE_path, of="GTiff")
                
                return(filevector)
            })
            
            
            # return(out)
        }
    ))}




# # #! /usr/bin/Rscript
# # # ##############################################################
# # # ##############################################################
# # # ##load libraries and CONFIG file
# # mainfolder = 'F:/PROGETTI_AGRONOMICO/LAMPvm/R_IrrWatB_Process/'
# # mainfolder = '//mnt/LAMPvm/R_IrrWatB_Process/'
# # 
# # setwd(mainfolder)
# # source('./CodeConfig/Config_SentinelProcess.R')
# # # ##############################################################
# # # ##Define folder input for RASTER: 1 folder each image
# # 
# # inputFolder<-"//mnt/LAMPvm/temp"
# # copyInputToOTHERFolder<-"//mnt/LAMPvm/tempTQ"
# # # #
# # listFolder<-list.files(inputFolder, full = TRUE)
# # 
# # # # # ###################################################
# # # # # ###################################################
# # # # # ###################################################
# # # # # # ###################################################
# # # # # # # # RUN POST PROCESS
# # # # # convert into WGS84 in other FOLDER
# # filenamepathOUTPUT<-"./ModelOUT/MLout_TEST/"# cartella dove saranno creati i raster 1-0 del modello
# # filenamepathOUTPUT_wgs84<-"./ModelOUT/MLout_wgs84/"# cartella dove saranno creati i raster 1-0 del modello
# 
# # input = filenamepathOUTPUT
# getwd()
# filenamepathOUTPUT
# filenamepathOUTPUT<-"./ModelOUT/MLout/"
# listFolder<-list.files(filenamepathOUTPUT, full = TRUE, pattern = "\\.tif$")
# # a partire dalla cartelal di input del modello ML "MLout_data" creo altra cartella con WGS84 converto il CRS ed unifico file per stessa data
# # creo la cartella di LAVORO 
# FolderMainToModified<-strsplit(filenamepathOUTPUT, "/")[[1]][(length(strsplit(filenamepathOUTPUT, "/")[[1]]))]
# PATHmain<-strsplit(filenamepathOUTPUT, "/")[[1]][0:(length(strsplit(filenamepathOUTPUT, "/")[[1]])-1)]
# pathMain_slash<-paste0(PATHmain, "/")
# pathMain_slashJoin<-paste0(pathMain_slash[1],pathMain_slash[2])
# newFolderToCreate_wgs84<-paste0(pathMain_slashJoin, FolderMainToModified,"_wgs84")
# dir.create(newFolderToCreate_wgs84, showWarnings = FALSE)
# # eseguo conversione di CRS    
# ExportRastesr_wgs84<-lapply(listFolder, function(x){
# 
#     listTofilename <- as.character(x)
#     readRasterF<-raster(listTofilename)
#     pr1 <- projectRaster(readRasterF, crs="+proj=longlat +datum=WGS84")
#     filname<-(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
#     outputFilePath <- paste0( filenamepathOUTPUT_wgs84 , filname)
#     print("##WGS84 filename path: ")
#     print(outputFilePath)
#     rc <- writeRaster(pr1, filename=outputFilePath,  fprj = TRUE,format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
#     writeLines(showWKT(crs(pr1, asText=TRUE)), extension(outputFilePath, 'prj') )
#     print("##end SAVE")
#     return(listTofilename)
#     })
# # faccio ciclo su CARTELLA di WGS84 per LISTA FILE
# listFolder_wgs84<-list.files(newFolderToCreate_wgs84, full = TRUE, pattern = "\\.tif$")
# # creo vettore con unique DATe nel caso per 1 giorno di download abbia più date
# dateFromFile<-lapply(listFolder_wgs84, function(x){
#                               fullpath<-x
#                               filname<-(strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])])
#                               all_dateFromFile<-strsplit(filname, "_")[[1]][2]
#                               dateFromFile<- substr(all_dateFromFile, 1, 8)
#                               export<- c( x, dateFromFile)
#                               return(export)})
# df_FielList<-as.data.frame(do.call(rbind, dateFromFile))
# datedownloaded<-unique(df_FielList[,2])
# # eseguo ciclo per ogni data unica prendo ed unisco in unico FILE
# Selected_RasterbyDATE_AND_PROCESS<-lapply(datedownloaded, function(x){
#                                                       SelectedITEM <- df_FielList[df_FielList$V2 == x,1]
#                                                       print("/*/*/* select ITEM:")
#                                                       # print(x)
#                                                       filevector<-c(as.character(unlist(SelectedITEM)))
#                                                       print(SelectedITEM)
# 
#                                                       outputFILE_path<-paste0(filenamepathOUTPUT_wgs84,"MERGED_wgs84_" , x, ".tif")
#                                                       print(outputFILE_path)
#                                                       # rast.mosaic <- do.call(mosaic,filevector)
#                                                       
#                                                       e <- extent(7, 17, 40, 50)
#                                                       template <- raster(e)
#                                                       projection(template) <- CRS('+init=epsg:4326')
#                                                       writeRaster(template, file=outputFILE_path, format="GTiff")
#                                                       mosaic_rasters(gdalfile=filevector,dst_dataset=outputFILE_path, of="GTiff")
# 
#                                                       return(filevector)
#                                                       })
# 
# 




