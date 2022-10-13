#/*/*!!!! IL POLIGONO DI TAGLIO DELL'AREA STUDIO DEVE ESSERE fatto RETTANGOLARE preciso sul 60m
#memory.limit()
##############################################################################
# Load Libraries
list.of.packages <- c("gdalUtils", "lubridate", "ggplot2", "raster",  "scales", "rgdal", "sp",  "zoo", "reshape2", "devtools",  "plyr", "IDPmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages, repos="https://cran.stat.unipd.it/")
lapply(list.of.packages, library, character.only = TRUE)
library(rgdal)
# ###################################################
# ###################################################
# # # INPUT  FILE PATH SETTINGS reading
# Mainpath <-"F:/PROGETTI_AGRONOMICO/LAMPvm/"
# Mainpath <-"C:/dbGIS/"
# shapefolder<-paste(Mainpath, "satCentos7/AcquaCAMPUS_shape", sep ="") # folder con shapefile per calcolo
# folderpath<- paste(Mainpath, "2A_new", sep ="") # folder con le cartelle del raster sentinel2
# pathOUTplot<-paste(Mainpath, "satCentos7/sentinelScript/outscript/plot/", sep ="")   # folder di output
# pathOUTtable <- paste(Mainpath, "satCentos7/sentinelScript/outscript/plot/outTable1.csv", sep ="") # file csv
# tempFolder<-paste(Mainpath, "temp", sep ="") # file csv  # folder dei file temp
# tempFolderClean<-paste(tempFolder, "/", sep = "")
######POLYGON AREAS
# shapefolder<-"//172.17.0.30/DocumentiCER/Agronomicoamb/Letterio/PSR2014/SOIA/Aziende/"
# shapefolder<-"/mnt/LAMPvm/satCentos7/AcquaCAMPUS_shape/"
# shapefolder<-"F:/PROGETTI_AGRONOMICO/LAMPvm/satCentos7/AcquaCAMPUS_shape/"
# shapefolder<-"C:/Users/tonpo/Documents/scriptLocal/R_sentinel_input/AcquaCAMPUS_shape"
######RASTER FOLDER
# folderpath<- "//172.17.0.30/giscer/PROGETTI_AGRONOMICO/LAMPvm/satCentos7/sentinelScript/outscript/2A/"
# folderpath<- "/mnt/LAMPvm/2A/"
# folderpath<- "C:/Users/tonpo/Documents/scriptLocal/R_sentinel_input/2a_input/"
######OUTPUT DIR
# pathOUTplot<- "/mnt/LAMPvm/satCentos7/sentinelScript/outscript/plot/"
# pathOUTplot<- "C:/Users/tonpo/Documents/scriptLocal/R_sentinel_input/out/"
# pathOUTtable <- "C:/Users/tonpo/Documents/scriptLocal/R_sentinel_input/out/outTable1.csv"
# pathOUTtable <- "/mnt/LAMPvm/satCentos7/sentinelScript/outscript/plot/outTable1.csv"
######TEMPFOLDER
# tempFolderClean<- "/mnt/LAMPvm/temp/"
# tempFolderClean<- "C:/Users/tonpo/Documents/scriptLocal/R_sentinel_input/temp/"
###################################################
# ###################################################
# SET TEMP DIR for RASTER and LARG FILE
# rasterOptions(tmpdir = "//172.17.0.30/giscer/PROGETTI_AGRONOMICO/temp")
# write("R_USER = //172.17.0.30/giscer/PROGETTI_AGRONOMICO/temp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
# rasterOptions(tmpdir = "/mnt/LAMPvm/temp")
# rasterOptions(tmpdir = tempFolder)
# rasterOptions(tmpdir = "F:/PROGETTI_AGRONOMICO/LAMPvm/temp")
# write("R_USER = /mnt/LAMPvm/temp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
# writeTempPath<-paste ("R_USER = ",  tempFolder, sep = "")
# write(writeTempPath, file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# write("R_USER = F:/PROGETTI_AGRONOMICO/LAMPvm/temp", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
# ###################################################
# ###################################################
#*#*# Process for folder analysis of sentinel2 raster data
##input : folder with multiple subfolder eachone as a image date as sentinel2
# ###################################################
# ###################################################
# #CLASSI:
# /*/GetRASTERFilePathANDArray:
# input: un folder con 1 subfolder per ogni sentinel2 image già bottom-of-atm quindi livello 2A
# metodi:
# 1) GetRasterArray: un oggetto contente tutti i raster del folder con
# rilfettanza delle bande e il Rastes SCL che servirà per MASK l'indice che
# indica la presenza di vegetazione e suolo nudo è 4-5, tutto il resto
# viene mascherato
# 2) GetRasterDATEvector: un vettore con le date delle immagini del sentinel
# 3) GetMASK: extract raster mask
# 4) Add raster Mask to raster brick
# ###################################################
# /*/CalcIndex:
# input un array con tutti i file raster di una data del sentinel
# 1) RStoolCalcINDEX:  un array dove sono calcolati tutti i gli indici vegetazionali escludendo i valori dei MASK
# ###################################################
# /*/ListToDf:
# input :
# ListVal :list of value as CalcIndex Obj output
# dateVec  :vector of date of Raster
# VectorShapeName :VectorOf Shape file Name where we extract value
# indexToCalc  :vector of index To calc from raster brick
# out: formatted Table
# method:
# 1)GiveFormattebTable: elabora l'output di CalcIndex formattando una tabella interpretabile
# ###################################################
# /*/PolygonValueSTAT, plot seasonal graph and generate a table output
# #input:
# 1)Array di rasterbrick: ogni rasterbrick è composto da 1 data e dal vari raster di indici vegetazionali
# 2) vettore con date dei raster bricks
# 3) polygon dove calcolare i valori ripeilogativi
# #output: tabella con statistiche per ogni indice vegetazionale
# ###################################################
# /*/calcoloStatistichesuPolygon, calcolo statistiche su polyg
# #input:
# 1)Array di rasterbrick: ogni rasterbrick è composto da 1 data e dal vari raster di indici vegetazionali
# 2) polygon dove calcolare i valori ripeilogativi
# #output: tabella con statistiche per ogni indice vegetazionale
# ###################################################
# /*/fileListShape, lista di shape file polgon dove calcolare valori
# #input:
# 1)Folde path:
# #output: vettore di shape file
# ###################################################
# /*/fileListShapeFileName, lista di  NOMI shape file polgon dove calcolare valori
# #input:
# 1)Folde path:
# #output: vettore di NOMI di shape file
# ###################################################
# /*/fplotTable, plotTable
# #input:
# 1)tabella
# #output: plot
# ###################################################
# /*/funT
# #input:
# 1)funzione per mapply per elaborazione output di CalcIndex
# #output: List di table
# ###################################################
# /*/ResElab
# #input:
# 1)output di mapply 
# #output: Tabella formattata
# ###################################################
# /*/PredictRasterWithRMmodelAnsSAVE: (_byRow)
# input: file raster brick con 5 indici CALCOLATI - path del model ML randomFOREST - e percorso folder e FILE TIFF per CREARE il file e Salvare l'output 
# metodi:
# 1) modelResultANDsave in folder 
# output:
# FileSALVATO in TIFF con risulati 1 IRR 2 NoIRR 
# ###################################################
GetRASTERFilePathANDArray <- function(input, cutAREA, outFolder, CutYorN, SaveNew) {
    structure(class = "GetRASTERFilePathANDArray", list(
    # attributes
    input = input,## path of folder contente tanti folder quante sono le date delle immagini da analizzare
		cutAREA =  cutAREA, #Polygon TOCUT polugono dell'area di studio
		outFolder = outFolder, #outputfolder dove ricreo una cartealla ogni cartella in input e copio i file crop
		CutYorN = CutYorN, #stringa con Y o N per fare ritaglio
		SaveNew = SaveNew, #stringa per sapere se salvare o no file raster elaborati 
    # methods #create raster ARRAY for each FOLDER 
	GetRasterArray_onlyREAD= function() {
				lapply(input, function (x) {
						# 
				  # browser()
						print(paste0("print x: ", x))
						# print(paste0("print input: ", input))
						FileList<-list.files(x,pattern = "*.grd", recursive = TRUE, full = TRUE)
						FileRaster1<-lapply(FileList, function(x){return(raster(x))})
						return(FileRaster1)})
				},
    # methods #create raster ARRAY for each FOLDER anc cut with STUDY AREA Polygon and SAVE in output PATH
    GetRasterArray	 = function() {
				myshp <- readOGR(cutAREA)
				e <- extent(myshp)
				ExtentValue <- e
				print(paste0("print ExtentValue: ", ExtentValue))
				lapply(input, function (x) { ## faccio ciclo per ogni sotto cartella trovata nel path main - un cartella per ogni immagine come scaricata dal sito esa
				  # VECTOR ID with FILENAME	1:B02_10m, 2:B03_10m, 3:B04_10m, 4:B08_10m, 25:SCL_20m , 21:B09_60m, 8:B05_20m,9:B06_20m,10:B07_20m, 11:_B11_20m, 12:B12_20m
				  MainFolder <- list.dirs(x, full.names=TRUE)
				  FileList_img_data <- MainFolder[ grepl("IMG_DATA", MainFolder) ]
				  B02_10m <-list.files(x,pattern = "*B02_10m(.*)jp2$", recursive = TRUE, full = TRUE)
				  B03_10m <-list.files(x,pattern = "*B03_10m(.*)jp2$", recursive = TRUE, full = TRUE)
				  B04_10m <-list.files(x,pattern = "*B04_10m(.*)jp2$", recursive = TRUE, full = TRUE) 
				  B08_10m <-list.files(x,pattern = "*B08_10m(.*)jp2$", recursive = TRUE, full = TRUE)  
				  SCL<-list.files(x,pattern = "*_SCL_20m(.*)jp2$", recursive = TRUE, full = TRUE)
				  B09_60m <-list.files(x,pattern = "*B09_60m(.*)jp2$", recursive = TRUE, full = TRUE)  
				  B05_20m <-list.files(x,pattern = "*B05_20m(.*)jp2$", recursive = TRUE, full = TRUE)   
				  B06_20m <-list.files(x,pattern = "*B06_20m(.*)jp2$", recursive = TRUE, full = TRUE) 
				  B07_20m <-list.files(x,pattern = "*B07_20m(.*)jp2$", recursive = TRUE, full = TRUE) 
				  B11_20m <-list.files(x,pattern = "*B11_20m(.*)jp2$", recursive = TRUE, full = TRUE) 
				  B12_20m <-list.files(x,pattern = "*B12_20m(.*)jp2$", recursive = TRUE, full = TRUE)  
				  FileList<-c(B02_10m, B03_10m, B04_10m, B08_10m, SCL, B09_60m, B05_20m, B06_20m, B07_20m, B11_20m, B12_20m)
				  # browser()
# 				  FileList<-list.files(x,pattern = "*_B(.*)jp2$", recursive = TRUE, full = TRUE) # trovo gruppi di file _B per ogni cartella di data
#   				FileList1<-list.files(x,pattern = "*_SCL_20m(.*)jp2$", recursive = TRUE, full = TRUE)
#   				FileList<-c(FileList,FileList1)
  				outputFOLDERS<-""
  				if (CutYorN == "N") {} # se non SALVO FILE non CREO nemmeno il FOLDER
  				if (SaveNew == "Y") {
  									print(paste0("print FileList: ", FileList[1]))
  									folderString<-length(head(unlist(strsplit(FileList[1], "/") )))
  									folderName<-head(unlist(strsplit(FileList[1], "/") ))[(folderString-1)]
  									print(paste0("print folderName: ", folderName))
  									outputFOLDERS<-paste(outFolder, folderName, sep="/")
  									dir.create(outputFOLDERS)}#creo un folder nel outputpath per ogni cartella di input
  				FileRaster1<-lapply(FileList, function(rasterF){ # faccio ciclo per tutti i file identificati nella cartella della DATA e TAGLIO area di studio e la salvo in nuovo percorso
  								if (CutYorN == "N") {myraster.crop <- raster(rasterF)}
  								else { myraster.crop <- crop(raster(rasterF), ExtentValue)} # ritaglio raster per area di studio
  								if (SaveNew == "Y") {
  													filenamePATH_split<-strsplit(rasterF, "/")
  													filenameShort<-tail(unlist(filenamePATH_split), n=1) # identifico nome file da filename path
  													outputFilePath<-paste(outputFOLDERS, str_replace(filenameShort, ".jp2", ".grd"), sep="/") # identifico filepath del file che voglio salvare
  													rc <- writeRaster(myraster.crop, filename=outputFilePath, prj = TRUE, format = "raster") # savlo file in directory più generale
  													writeLines(showWKT(crs(myraster.crop, asText=TRUE)), extension(outputFilePath, 'prj') )
  													# # unlink(paste0(filePATHdaCancellare,"/*"), recursive = TRUE)
  													# print(paste0("####print extent out: ", extent(raster(outputFilePath))))
  													}
  								return(myraster.crop)}
  								)
  				return(FileRaster1)
  				})
				# unlink(paste0(filePATHdaCancellare,"/*"), recursive = TRUE)

				},
    GetMASK	 = function() {
				lapply(input, function (x) {
				FileList<-list.files(x,pattern = "*_SCL_20m(.*)jp2$", recursive = TRUE, full = TRUE)
				FileRaster<-lapply(FileList, function(x){return(raster(x))})
				return(unlist(FileRaster))})
				},
    # methods #create date vector from raster ARRAY for each FOLDER 
	GetRasterDATEvector	 = function() {
				lapply(input, function (x) {
				FileList<-list.files(x,pattern = "*_B(.*)jp2$", recursive = TRUE, full = FALSE)
				FileN<-strsplit(FileList[1], "/")
				print(FileN)
				LatItem = (length(FileN[[1]]))
				outName<-lapply(FileN, function(x){unlist(x)[LatItem]})
				print(outName)
				DateFromFileName<-unlist(strsplit(outName[[1]], "_"))[2]
				print(DateFromFileName)
				year<-substr(DateFromFileName,1,4)
				month<-substr(DateFromFileName,5,6)
				day<-substr(DateFromFileName,7,8)
				ndate1 <- as.Date(paste(year, month, day, sep="-"), "%Y-%m-%d")
				print(DateFromFileName)
				print(ndate1)
				return(unlist(ndate1))})
				},
	# methods #create date vector from raster ARRAY for each FOLDER  READING CUT RESULTS folders
	GetRasterDATEvectorTIFF	 = function() {
				lapply(input, function (x) {
				FileList<-list.files(x,pattern = "*.grd*", recursive = FALSE, full = FALSE)[1]
				# FileN<-strsplit(FileList[1], "/")
				# print(FileN)
				# outName<-lapply(FileN, function(x){unlist(x)[1]})
				# print(outName)
				DateFromFileName<-unlist(strsplit(FileList[[1]], "_"))[2]
				print(paste0("### DateFromFileName: ",DateFromFileName))
				year<-substr(DateFromFileName,1,4)
				month<-substr(DateFromFileName,5,6)
				day<-substr(DateFromFileName,7,8)
				ndate1 <- as.Date(paste(year, month, day, sep="-"), "%Y-%m-%d")
				print(DateFromFileName)
				print(ndate1)
				# browser()
				return(unlist(ndate1))})
				},
	# methods #create NAMEandDATE VECTOR from TIFF reading methods
	GetRasterNAMEvectorTIFF	 = function() {
				lapply(input, function (x) {
				FileList<-list.files(x,pattern = "*.grd*", recursive = FALSE, full = FALSE)[1]
				# browser()
				NAME_DateFromFileName<-unlist(strsplit(FileList[[1]], "_"))[1:2]
				NAME_DateFromFileNameJOIN<-paste(NAME_DateFromFileName[1],NAME_DateFromFileName[2], sep = "_" )	
				print(paste0("### DateFromFileName: ",NAME_DateFromFileNameJOIN))
				return(unlist(NAME_DateFromFileNameJOIN))})
				},
	# methods #create NAMEandDATE VECTOR from TIFF reading methods
	GetRasterNAMEvector	 = function() {
				lapply(input, function (x) {
					FileList<-unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]
					# browser()
					NAME_DateFromFileName<-unlist(strsplit(FileList[[1]], "_"))[6:7]
					NAME_DateFromFileNameJOIN<-paste(NAME_DateFromFileName[1],NAME_DateFromFileName[2], sep = "_" )	
					print(paste0("### DateFromFileName: ",NAME_DateFromFileNameJOIN))
					return(unlist(NAME_DateFromFileNameJOIN))})
				}
		))}
####################################################
PredictRasterWithRMmodelAndSAVE <- function(raster, filenamepath, modelpath, threNDVI, threGNDVI) {
    structure(class = "PredictRasterWithRMmodel", list(
    # attributes
	raster = raster,  
	filenamepath = filenamepath,
	modelpath = modelpath,
	threNDVI = threNDVI,  # threshold NDVI for exclunding pixel with lower values
	threGNDVI = threGNDVI, # threshold GNDVI for exclunding pixel with lower values
	# methods
	CalculateModelML = function() {	
		print("in calc method")
		modelML <- readRDS(modelpath)
		outRaster<-raster(raster)
		out <- writeStart(outRaster, filenamepath, format = 'GTiff', overwrite=TRUE, options=c("COMPRESS=NONE", "TFW=YES"))
		# out <- writeStart(outRaster, filenamepath, overwrite=TRUE)
		bs <- blockSize(out)
			for (i in 1:bs$n) {
				# browser()
				v1 <- getValues(raster[[1]], row=bs$row[i], nrows=bs$nrows[i] )
				v2 <- getValues(raster[[2]], row=bs$row[i], nrows=bs$nrows[i] )
				v3 <- getValues(raster[[3]], row=bs$row[i], nrows=bs$nrows[i] )
				v4 <- getValues(raster[[4]], row=bs$row[i], nrows=bs$nrows[i] )
				v5 <- getValues(raster[[5]], row=bs$row[i], nrows=bs$nrows[i] )
				v6 <- getValues(raster[[6]], row=bs$row[i], nrows=bs$nrows[i] )
				# print(paste0("##v1: ", v1))
				# print(paste0("##Block number: ", str(i)))
				db<-as.data.frame(cbind(v1,v2,v3,v4,v5,v6))
				colnames(db)<-c("DVI", "GNDVI", "MNDWI", "NDVI", "NDWI2" , "SATVI")
				# print("##DB")
				# print(db)
				head(db)
				y_pred <- predict(modelML, typre = 'response', newdata = db, na.rm=TRUE, inf.rm=TRUE)
				vectPRED<-as.numeric(as.character(y_pred))
				# print(paste0("##length y pred: ", length(vectPRED)))
				print(str(y_pred))
				Vectoroutput<-Filteringoutput_thre(y_pred, v2, v4, threNDVI, threGNDVI)
				# print("##str filter predict: ")
				print(str(Vectoroutput))
				# browser()
				# print("##stop filter predict")
				out <- writeValues(out, Vectoroutput, bs$row[i])
				# out <- writeValues(out, y_pred, bs$row[i])
				gc()
				}
		out <- writeStop(out)
		return(out)
		}
	))}
####################################################
####################################################
PredictRasterWithRMmodelAndSAVE_Complete <- function(rasterfile, filenamepath, modelpath, threNDVI, threGNDVI) {
    structure(class = "PredictRasterWithRMmodel", list(
    # attributes
	rasterfile = rasterfile,  
	filenamepath = filenamepath,
	modelpath = modelpath,
	threNDVI = threNDVI,  # threshold NDVI for exclunding pixel with lower values
	threGNDVI = threGNDVI, # threshold GNDVI for exclunding pixel with lower values
	# methods
	CalculateModelML = function() {	
		print("in calc method")
		modelML <- readRDS(modelpath)
		outRaster<-raster(rasterfile)
		# browser()
		big <- !canProcessInMemory(outRaster, 2)
		filenamepath <- trim(filenamepath)
		# out <- writeStart(outRaster, filenamepath, format = 'GTiff', overwrite=TRUE, options=c("COMPRESS=NONE", "TFW=YES"))
		if (big & filenamepath == '') {
					filenamepath <- rasterTmpFile()
					}
		if (filenamepath != '') {
						# out <- writeStart(outRaster, filenamepath, overwrite=TRUE)
						out <- writeStart(outRaster, filenamepath, format = 'GTiff', overwrite=TRUE, options=c("COMPRESS=NONE", "TFW=YES"))
						todisk <- TRUE
					} else {
						print("##please insert valid raster filename path for output")
					}
		bs <- blockSize(outRaster)
		pb <- pbCreate(bs$n)
		if (todisk) {
			for (i in 1:bs$n) {
					# browser()
					v1 <- getValues(rasterfile[[1]], row=bs$row[i], nrows=bs$nrows[i] )
					v2 <- getValues(rasterfile[[2]], row=bs$row[i], nrows=bs$nrows[i] )
					v3 <- getValues(rasterfile[[3]], row=bs$row[i], nrows=bs$nrows[i] )
					v4 <- getValues(rasterfile[[4]], row=bs$row[i], nrows=bs$nrows[i] )
					v5 <- getValues(rasterfile[[5]], row=bs$row[i], nrows=bs$nrows[i] )
					v6 <- getValues(rasterfile[[6]], row=bs$row[i], nrows=bs$nrows[i] )
					# print(paste0("##v1: ", v1))
					print(paste0("##Block number: ", str(i)))
					db<-as.data.frame(cbind(v1,v2,v3,v4,v5,v6))
					colnames(db)<-c("DVI", "GNDVI", "MNDWI", "NDVI", "NDWI2" , "SATVI")
					# print("##DB")
					# print(db)
					head(db)
					y_pred <- predict(modelML, typre = 'response', newdata = db, na.rm=TRUE, inf.rm=TRUE)
					vectPRED<-as.numeric(as.character(y_pred))
					print(paste0("##length y pred: ", length(vectPRED)))
					print(str(y_pred))
					Vectoroutput<-Filteringoutput_thre(y_pred, v2, v4, threNDVI, threGNDVI)
					print("##str filter predict: ")
					print(str(Vectoroutput))
					# browser()
					print("##stop filter predict")
					out <- writeValues(out, Vectoroutput, bs$row[i])
					pbStep(pb, i)
					gc()
					# out <- writeValues(out, y_pred, bs$row[i])
					}
			out <- writeStop(out)
			}
		pbClose(pb)
		return(out)
		}
	)
	)
	}
####################################################
PredictRasterWithRMmodelAndSAVE_byROW <- function(raster, filenamepath, modelpath, threNDVI, threGNDVI) {
    structure(class = "PredictRasterWithRMmodel", list(
    # attributes
	raster = raster,  
	filenamepath = filenamepath,
	modelpath = modelpath,
	threNDVI = threNDVI,  # threshold NDVI for exclunding pixel with lower values
	threGNDVI = threGNDVI, # threshold GNDVI for exclunding pixel with lower values
	# methods
	CalculateModelML = function() {	
		# print("in calc method")
		modelML <- readRDS(modelpath)
		outRaster<-raster(raster)
		out <- writeStart(outRaster, filenamepath, format = 'GTiff', overwrite=TRUE, options=c("COMPRESS=NONE", "TFW=YES"))
		# out <- writeStart(outRaster, filenamepath, overwrite=TRUE)
		for (i in 1:nrow(out)) {
				# browser()
				v1 <- getValues(raster[[1]], i)
				v2 <- getValues(raster[[2]], i)
				v3 <- getValues(raster[[3]], i)
				v4 <- getValues(raster[[4]], i)
				v5 <- getValues(raster[[5]], i)
				v6 <- getValues(raster[[6]], i)
				# print(paste0("##v1: ", v1))
				# print(paste0("##Block number: ", str(i)))
				db<-as.data.frame(cbind(v1,v2,v3,v4,v5,v6))
				colnames(db)<-c("DVI", "GNDVI", "MNDWI", "NDVI", "NDWI2" , "SATVI")
				# print("##DB")
				# print(db)
				# head(db)
				y_pred <- predict(modelML, typre = 'response', newdata = db, na.rm=TRUE, inf.rm=TRUE)
				vectPRED<-as.numeric(as.character(y_pred))
				# print(paste0("##length y pred: ", length(vectPRED)))
				# print(str(y_pred))
				Vectoroutput<-Filteringoutput_thre(y_pred, v2, v4, threNDVI, threGNDVI)
				# print("##str filter predict: ")
				# print(str(Vectoroutput))
				# browser()
				# print("##stop filter predict")
				out <- writeValues(out, Vectoroutput, i)
				# out <- writeValues(out, y_pred, bs$row[i])
				gc()
				}
		out <- writeStop(out)
		return(out)
		}
	))}
####################################################
####################################################
CalcIndex <- function(x, pol, indexVector, idraster, filenamepath) {
				structure(class = "CalcIndex", list(
				# attributes
					x = x, #raster bricks
					pol = pol, #polygon to extract
					indexVector = indexVector, #index of vector Name to Calc as for spectralIndices Function
					idraster = idraster, #index of vector Name to Calc as for spectralIndices Function
					filenamepath = filenamepath, #filename path dove salvare il brick con SI 
				# methods 1: read and calculate INDEX after reading SENTINEL raw data
		RStoolCalcINDEX = function(){
					blu <- idraster[1]
					green <- idraster[2]
					red <- idraster[3]
					nir <- idraster[4]
					maskSCL <- idraster[5]
					swir1 <- idraster[6]
					re1<- idraster[7]
					re2<- idraster[8]
					re3<- idraster[9]
					swir2<- idraster[10]
					swir3<- idraster[11]
					# print(paste0("###briks_raster: ", x))
					# browser()
					maskLayer <- raster::disaggregate(x[maskSCL][[1]], fact=2) #mask layer
					print(paste0("###mask layer: ", maskLayer))
					m <- c(-Inf, 3, 1,  3.1, 5, 2,  5.1, 12, 1) #reclassify matrix Metto tutto 1 Dove faccio Maschea, Lascio 2 solo su vegetazione o terreno nudo
					# m <- c(-Inf, 3, 1,  3.1, 10, 2,  10.1, 12, 1) #reclassify matrix Metto tutto 1 Dove faccio Maschea, Lascio 2 solo su vegetazione o terreno nudo
					rclmat <- matrix(m, ncol=3, byrow=TRUE)
					rc_main <- reclassify(maskLayer, rclmat) #reclassify values 1 NOdata, 2 OkValue
					# plot(rc)
					SWIRc <- raster::disaggregate(x[swir1][[1]], fact=6)
					SWIR2band <- raster::disaggregate(x[swir2][[1]], fact=2)
					SWIR3band <- raster::disaggregate(x[swir3][[1]], fact=2)
					rededge1Raster <- raster::disaggregate(x[re1][[1]], fact=2)
					rededge2Raster <- raster::disaggregate(x[re2][[1]], fact=2)
					rededge3Raster <- raster::disaggregate(x[re3][[1]], fact=2)
					# reample masked layer to same extent (to avoid compareraster ERROR)
					rc <- projectRaster(rc_main,SWIRc,method = 'ngb')
					# plot(rededge3Raster)
					print(paste0("###rededge1Raster: ", rededge1Raster))
					print(paste0("###rededge2Raster: ", rededge2Raster))
					print(paste0("###rededge3Raster: ", rededge3Raster))
					print(paste0("###SWIR: ", SWIRc))
					print(paste0("###SWIR: ", SWIR3band))
					print(paste0("###mask layer RC: ", rc))
					# browser()
					redR<-x[red][[1]] ## remove problem on RSTOOLBox li quando trova valori di RED < 1.5 non fa calcolo di EVI
					redR[redR>15000] <- 15000 #all values > di 1.5 diventano 1.5
					Rstack <- stack(x[blu][[1]], x[green][[1]], redR, x[nir][[1]], rededge1Raster, rededge2Raster, rededge3Raster, SWIRc, SWIR2band, SWIR3band)
					# SI <- spectralIndices(Rstack, blu=1, green=2, red=3, redEdge1=4, redEdge2=5, redEdge3=6, nir=7,  swir2=5, swir3=9 ,scaleFactor = 10000, maskLayer = rc, maskValue = 1,  index = indexVector)
					SI<- spectralIndices(Rstack, blue = 1, green = 2, red = 3, nir = 4, redEdge1 = 5, redEdge2 = 6, redEdge3 = 7, swir1 = NULL, swir2 = 9, swir3 = 10, scaleFactor = 10000, index = indexVector, maskLayer = rc, maskValue = 1, coefs = list(L = 0.5, G = 2.5,  L_evi = 1, C1 = 6, C2 = 7.5, s = 1, swir2ccc = NULL, swir2coc = NULL))
					# SItest <- spectralIndices(Rstack,blue = 1, green = 2, red = 3, nir = 4, redEdge1 = 5, redEdge2 = 6, redEdge3 = 7, swir1 = NULL, swir2 = 9, swir3 = 10, index = ('EVI'), maskLayer = rc, maskValue = 1 , coefs = list(L = 0.5, G = 2.5,  L_evi = 1, C1 = 6, C2 = 7.5, s = 1, swir2ccc = NULL, swir2coc = NULL))
					# SI <- spectralIndices(Rstack,blue = 1, green = 2, red = 3, nir = 4, redEdge1 = 5, redEdge2 = 6, redEdge3 = 7, swir1 = NULL, swir2 = 9, swir3 = 10, scaleFactor = 10000, index = indexVector, maskLayer = rc, maskValue = 1 , coefs = list(L = 0.5, G = 2.5,  L_evi = 1, C1 = 6, C2 = 7.5, s = 1, swir2ccc = NULL, swir2coc = NULL))
					print(paste0("###Rstack: ", Rstack))
					print("#######SI###")
					print(paste0("###SI stack: ", SI))
					# browser()
					# print(SI[1])
					# print(SI[2])
					ExtractEachPol<-DefinStat(SI, pol) #estrai statistiche ZONALI per singolo shaoe con + feature o multiple polygion
					# ExtractEachPol<-lapply(pol, function(feat){
															# print("#############")
															# print(paste0("###SI: ",SI))
															# # print(paste0("###feat: ", feat ))
															# rextract <- raster::extract(SI, feat, df = TRUE)
															# print(paste0("###rextract: ", rextract))
															# rextractOUT<-apply(rextract,2, my_summary)
															# rextractOUTTab<-cbind.data.frame(Val = c("avg", "med", "min", "max", "p25", "p75", "NA"), rextractOUT)
															# # print(rextractOUT)
															# # browser()
															# return (rextractOUTTab)})
					# rextractt <- raster::extract(SI, pol[[1]], df = TRUE)
					# summary(rextractt)	
					print(paste0("###ExtractEachPol: ",ExtractEachPol))
					return(ExtractEachPol)		
					},
# methods 2: read and calculate INDEX after cutting DATA in grd FORMAT	
		RStoolCalcINDEX_NOT_polygon = function(){
		  # VECTOR ID with FILENAME	1:B02_10m, 2:B03_10m, 3:B04_10m, 4:B08_10m, 5:SCL_20m , 6:B09_60m, 7:B05_20m,8:B06_20m,9:B07_20m, 10:_B11_20m, 11:B12_20m
					blu <- idraster[1]
					green <- idraster[2]
					red <- idraster[3]
					nir <- idraster[4]
					maskSCL <- idraster[5]
					swir1 <- idraster[6]
					re1<- idraster[7]
					re2<- idraster[8]
					re3<- idraster[9]
					swir2<- idraster[10]
					swir3<- idraster[11]
					# browser()
					# print(paste0("###briks_raster: ", x))
					maskLayer <- raster::disaggregate(x[maskSCL][[1]], fact=2) #mask layer
					print(paste0("###mask layer: ", maskLayer))
					m <- c(-Inf, 3, 1,  3.1, 5, 2,  5.1, 12, 1) #reclassify matrix Metto tutto 1 Dove faccio Maschea, Lascio 2 solo su vegetazione o terreno nudo
					# m <- c(-Inf, 3, 1,  3.1, 10, 2,  10.1, 12, 1) #reclassify matrix Metto tutto 1 Dove faccio Maschea, Lascio 2 solo su vegetazione o terreno nudo
					rclmat <- matrix(m, ncol=3, byrow=TRUE)
					rc_main <- reclassify(maskLayer, rclmat) #reclassify values 1 NOdata, 2 OkValue
					# plot(rc)
					SWIRc <- raster::disaggregate(x[swir1][[1]], fact=6)
					SWIR2band <- raster::disaggregate(x[swir2][[1]], fact=2)
					SWIR3band <- raster::disaggregate(x[swir3][[1]], fact=2)
					rededge1Raster <- raster::disaggregate(x[re1][[1]], fact=2)
					rededge2Raster <- raster::disaggregate(x[re2][[1]], fact=2)
					rededge3Raster <- raster::disaggregate(x[re3][[1]], fact=2)
					# resample masked layer to same extent (to avoid compare raster ERROR)
					rc <- projectRaster(rc_main,SWIRc,method = 'ngb')
					# plot(rededge3Raster)
					print(paste0("###rededge1Raster: ", rededge1Raster))
					print(paste0("###rededge2Raster: ", rededge2Raster))
					print(paste0("###rededge3Raster: ", rededge3Raster))
					print(paste0("###SWIR: ", SWIRc))
					print(paste0("###SWIR: ", SWIR3band))
					print(paste0("###mask layer RC: ", rc))
					# redR<-x[red][[1]] ## remove problem on RSTOOLBox li quando trova valori di RED < 1.5 non fa calcolo di EVI
					# redR[redR>15000] <- 15000 #all values > di 1.5 diventano 1.5  #RIMOSSO 10/05/2022 funziona lo stesso con X[red]
					# browser()
					# Rstack <- stack(x[blu][[1]], x[green][[1]], redR, x[nir][[1]], rededge1Raster, rededge2Raster, rededge3Raster, SWIRc, SWIR2band, SWIR3band)
					Rstack <- stack(x[blu][[1]], x[green][[1]], x[red][[1]], x[nir][[1]], rededge1Raster, rededge2Raster, rededge3Raster, SWIRc, SWIR2band, SWIR3band)
          # SI <- spectralIndices(Rstack, blue = 1, green = 2, red = 3, nir = 4, redEdge1 = 5, redEdge2 = 6, redEdge3 = 7, swir1 = NULL, swir2 = 9, swir3 = 10, scaleFactor = 10000, index = c('EVI'), maskLayer = rc, maskValue = 1, coefs = list(L = 0.5, G = 2.5,  L_evi = 1, C1 = 6, C2 = 7.5, s = 1, swir2ccc = NULL, swir2coc = NULL))
					# SI <- spectralIndices(Rstack, blu=1, green=2, red=3, redEdge1=4, redEdge2=5, redEdge3=6, nir=7,  swir2=5, swir3=9 ,scaleFactor = 10000, maskLayer = rc, maskValue = 1,  index = indexVector)
					SI <- spectralIndices(Rstack, blue = 1, green = 2, red = 3, nir = 4, redEdge1 = 5, redEdge2 = 6, redEdge3 = 7, swir1 = NULL, swir2 = 9, swir3 = 10, scaleFactor = 10000, index = indexVector, maskLayer = rc, maskValue = 1, coefs = list(L = 0.5, G = 2.5,  L_evi = 1, C1 = 6, C2 = 7.5, s = 1, swir2ccc = NULL, swir2coc = NULL))
					print(paste0("###Rstack: ", Rstack))
					print("#######SI###")
					print(paste0("###SI stack: ", SI))
					# browser()
					#save brick raster to folder
					filenameOUTPUT<-paste(unlist(strsplit(x[[1]]@data@names, '_'))[1], unlist(strsplit(x[[1]]@data@names, '_'))[2], sep="_")
					fileOutputComplete<-paste0(filenamepath, filenameOUTPUT)
					print(paste0("###FileToSave: ", fileOutputComplete))
					# rc <- writeRaster(SI, filename=fileOutputComplete, prj = TRUE, format = "raster", bylayer=TRUE, suffix='names') 			
					rc <- writeRaster(SI, filename=fileOutputComplete,  format = "raster", bylayer=FALSE, suffix='names') #correzione 09-03-2021 rimosso prj =TRUE perch? dava errore emsso bylaer = false x avere un BRICK
					# rc <- writeRaster(SI, filename=fileOutputComplete, prj = TRUE,format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
					# writeLines(showWKT(crs(x[[1]], asText=TRUE)), extension(fileOutputComplete, 'prj') )
					print("###file SAVED!!!")
					return(SI)		
					}
				   ))}
####################################################
ListToDf <- function(ListVal, dateVec, VectorShapeName, indexToCalc) {
    structure(class = "ListToDf", list(
    # attributes
        ListVal = ListVal, #list of value as CalcIndex Obj output
		dateVec = dateVec, #vector of date of Raster
		VectorShapeName = VectorShapeName, #VectorOf Shape file Name where we extract value
		indexToCalc = indexToCalc, #vector of index To calc from raster brick
	# methods
		GiveFormattebTable = function(){
				ret<-sapply(ListVal, function(x){
															Vector<-unlist(data.frame((strsplit(unlist(x), ":")))[2,], use.names=FALSE)
															vector1<-as.numeric(as.character(Vector))
															return(vector1)
															})
				print("##PRINT-head(ret)")
				print(head(ret))
				VecIndex<-sapply(IndexToCalc, function(x){rep(x,7)})
				VecInd<-c(rep("ID", 7), VecIndex) #creo vettore con ID degli INDICI
				vectoLoc<-as.vector(sapply(VectorShapeName, function(x){rep(x,length(VectorShapeName)*7)}))# creo vettore con NomideiFILE degli shape di input
				vec<-c("min", "q25", "med", "avg", "q75", "Max", "NaC")
				VecVal<-rep(vec, length(VectorShapeName))#vettore con valori calcolati
				VecIndTo<-rep(VecInd, length(VectorShapeName))
				Dbdef<-as.data.frame(cbind(vectoLoc,VecIndTo, VecVal, ret))
				print("##PRINT-head(Dbdef)")
				print(head(Dbdef))
				mdata <- melt(Dbdef, id.vars=c("vectoLoc","VecIndTo", "VecVal"))
				mdata_Red<-mdata[mdata$VecIndTo != "ID", ]
				VectorDate<-as.Date(sapply(ArrayDateVector, function(x){
													return(rep(x, length(IndexToCalc)*7))
													}))
				TabDEf<-cbind(as.Date(VectorDate), mdata_Red)
		return(TabDEf)
		}
	   ))}
####################################################
PolygonValueSTAT <- function(arrayBricks, datevector,  pol) {
    structure(class = "PolygonValueSTAT", list(
    # attributes
        arrayBricks = arrayBricks,
		datevector = datevector,
		pol=pol,
		VectorIndex<-labels(arrayBricks[[1]]), #estraggo dal primo raster brick il numero di indici
		NumeroCicli<-length(VectorIndex), #conto per il ciclo FOR
		print(NumeroCicli),
		a <- data.frame (),# data frame for output value
    # methods
		StatPolygValue = function(){
		for (i in 1:NumeroCicli){
			index<-sapply(arrayBricks, function(x){
										rasterFile<-x[[i]]
										OUTdt<-calcoloStatistichesuPolygon(rasterFile, pol) #estraggo statistiche zonali
										print("#####OUTdt")
										print(OUTdt)
										return(OUTdt[1:6])})
			print("t(index)")
			print(t(index))
			indextbOUT<-data.frame(t(index),VectorIndex[i], datevector)
			print("indextbOUT")
			print(indextbOUT)
			a <- rbind (a, indextbOUT)
			print(a)
			}
		return(a)}
	   ))}
###############################################################
###############################################################
###############################################################
###############################################################
# #FUNZIONI:
# /*/calcoloStatistichesuPolygon:
# # input: 
# 1) BRICK dalla funzione SpectralIndices
# 2) poligono o poligonI su cui calcolare le statistiche
# # process:
# per il multi poligono fa un ciclo per ogni poligono ed estrare le statistiche da ogni raster 
# per il singolo poligono con + feature prende media
# # ouput:
# nel caso di Singolo Poligono avrai altro shapefile con MEDIA dei valori dei pixel per i poligono
# nel caso di MULTI POLIGONI un DATAFRAME con per ogni poligono e per ogni raster le statistiche
DefinStat<-function(SI, pol) {if(is.list(pol)){
											print("#############")
											print(paste0("###SI: ",SI))
											ExtractEachPol<-lapply(pol, function(feat){
												# browser()
												# str(feat@polygons[[1]])
												# feat@polygons[[1]]@ID
												# if (feat@polygons[[1]]@ID == '10'){browser()}
												rextract <- raster::extract(SI, feat, df = TRUE)
												print(paste0("###rextract: ", rextract))
												# rextractOUT<-my_summary(rextract)
												rextractOUT<-apply(rextract,2, my_summary)
												print(paste0("###rextractOUT: ", rextractOUT))
												rextractOUTTab<-cbind.data.frame(Val = c("avg", "med", "min", "max", "p25", "p75", "NA"), rextractOUT)
												return(rextractOUTTab)
												}
												)
										}else{
											ExtractEachPol<-raster::extract(SI, pol, df=TRUE, fun=mean, sp=TRUE)
											}
										return(ExtractEachPol)
										}
###############################################################
# /*/calcolo EVI:
# # input: scale factor + 3 bande + coefficienti
# # ouput: indice EVI
CalculateEVI<-function(ScaleFactor, blue, red, nir, G, C1, C2, L_evi){
						index <- G * (((nir/ScaleFactor) - (red/ScaleFactor))/((nir/ScaleFactor) + C1 * red/ScaleFactor - C2 * (blue/ScaleFactor) + L_evi))
						return(index)
						}	
###############################################################
# /*/calcoloStatistichesuPolygon:
# # input:
# 1) Raster input dove ho calcolato INDICE dal quale estrarre informazioni su certa AREA
# 2) polygono dove calcolare valori
# # ouput:
# vettore con 6 elementi media, mediana e devst.
# ###################################################
calcoloStatistichesuPolygon<-function(x, polys){
											rextract <- raster::extract(x[[1]], polys, df = TRUE)
											# str(rextract)
											Vector<-unlist(data.frame((strsplit(summary(rextract[2], na.rm=TRUE), ":")))[2,], use.names=FALSE)
											vector1<-as.numeric(as.character(Vector))
											return(vector1)
										}
# /*/fileListShape:
# # input:
# 1) Folder path
# # ouput:
# vettore shapefile nei quali calcolare il le statistiche sui raster
# ###################################################
fileListShape<-function(vect){
							FileList<-list.files(vect,pattern = "*.shp", recursive = TRUE, full = TRUE)
							rout<-lapply(FileList, function(x){
										polys1 <-shapefile(x)
										print(polys1)
										polys <- spTransform(polys1, CRS=CRS("+init=epsg:32632"))
										return(polys)})
							return(rout)}
# /*/fileListShapeFileName:
# # input:
# 1) Folder path
# # ouput:
# vettore shapefile nei quali calcolare il le statistiche sui raster
# ###################################################
fileListShapeFileName<-function(vect){
							FileList<-list.files(vect,pattern = "*.shp", recursive = F, full = FALSE)
							return(FileList)}
# ###################################################
# /*/mysummary:
# # input:
# 1) vettore
# # ouput:
# summary con valore di Na
# ###################################################
my_summary <- function(v){
						meanV<-mean(v,  na.rm=TRUE)
						medianV<-median(v,  na.rm=TRUE)
						minV<-min(v,  na.rm=TRUE)
						maxV<-max(v,  na.rm=TRUE)
						per25<-quantile(v, probs = c(0.25),  na.rm=TRUE)
						per75<-quantile(v, probs = c(0.75),  na.rm=TRUE)
						NAvec<-sum(is.na(v))
						result<-c(meanV, medianV, minV, maxV, per25, per75, NAvec)
						return(result)
						}					
 ###################################################
# /*/plotTimeSPolygon:
# # input:
# 1) tabella con indice come output da PolygonValueSTAT
# # ouput:
# plot  + tabella riepilogativa
# ###################################################

plotTimeSPolygon<-function(res1){
									colnames(res1)<-c("min", "Q1", "Med", "Mean", "Q3","Max","vi", "date")
									plot1<-ggplot(data = res1, aes(x = date, y = Mean))+
									 geom_line(aes(color = vi), lwd = 1)+
									 geom_point(aes(color = vi))  +
									 geom_errorbar(aes(ymin = Q1, ymax = Q3, color = vi), width=0.8, position=position_dodge(0.2))
									print(plot1)
									return(res1)
									}
# ###################################################
# /*/plotTable:
# # input:
# 1) tabella come output da ListToDf
# # ouput:
# plot
# ###################################################
plotTable<-function(Tabout){
						colnames(Tabout)[1]<-"VectorDate"
						Tabout$value<-as.numeric(as.character(Tabout$value))
						TabCast<-dcast(Tabout,  VectorDate+vectoLoc+VecIndTo~VecVal)
						names(TabCast)
						PlotDBtotRain <- ggplot(TabCast, aes(VectorDate, med)) +
										facet_grid(vectoLoc ~ .) +
										geom_line(aes(color = factor(VecIndTo))) +
										geom_line(aes(color = factor(VecIndTo)))+
										theme(legend.text=element_text(size=8))+
										theme(legend.position = "bottom")+
										ylim(0,1)+
										geom_errorbar(aes(ymin=q25, ymax=q75), width=.2,
										 position=position_dodge(0.05))
						PlotDBtotRain
						}
# ###################################################
# /*/funT:
# # input:
#  CalcIndexForEachFolder, ArrayDateVector
# # ouput:
# tabella formattata
# ###################################################						
funT<-function(x,y){print(x)
					print(y)
					res1<-lapply(x, function (z){
												Vector<-data.frame(z)
												print("##DbVal")
												DbVal<-as.data.frame(t(as.data.frame(strsplit(as.character(Vector$Freq), ":"))))
												print(DbVal)
												# print(colnames(DbVal))
												Vec1<-as.vector(DbVal$V1)
												Vec2<-as.vector(DbVal$V2)
												VecDef<-as.data.frame(cbind(Vec1, Vec2))
												print(Vec1)
												Sel<-VecDef[VecDef$Vec1 != c("NA's   "),]
												Sel<-Sel[!is.na(Sel$Vec1),]
												print("######Sel")
												print(Sel)
												VectorNAME<-as.data.frame(cbind(y, Sel))
												colnames(VectorNAME)[1]<-c("Date")
												# browser()
												return(VectorNAME)})
					return(res1)
					}	
# ###################################################
# /*/ResElab:
# # input:
#  res di MAPPLY
# # ouput:
# tabella formattata
# ###################################################	
ResElab<-function(res, indexTC, ArrayDateVector){
print(indexTC)
# browser()
resMedian<-apply(res,1,function(x){
						# print(x)
						meanV<-lapply(x, function(y){	
													# print(y$Vec1)
													outMean<-y[y$Vec1 == "Median ",]
													return(outMean)
													})
						df <- do.call(rbind, meanV)
						VecIn<- c("i", indexTC)
						RepV<-rep(VecIn, length(indexTC))
						RepV1<-rep(RepV, length(ArrayDateVector))
						# browser()
						dfdef<-cbind(RepV1, df)
						return(dfdef)
						})
resMean<-apply(res,1,function(x){
						print(x)
						meanV<-lapply(x, function(y){	
													# print(y$Vec1)
													outMean<-y[y$Vec1 == "Mean   ",]
													return(outMean)
													})
						df <- do.call(rbind, meanV)
						VecIn<- c("i", indexTC)
						RepV<-rep(VecIn, length(indexTC))
						RepV1<-rep(RepV, length(ArrayDateVector))
						dfdef<-cbind(RepV1, df)
						return(dfdef)
						})
resMin<-apply(res,1,function(x){
						print(x)
						meanV<-lapply(x, function(y){	
													# print(y$Vec1)
													outMean<-y[y$Vec1 == "Min.   ",]
													return(outMean)
													})
						df <- do.call(rbind, meanV)
						VecIn<- c("i", indexTC)
						RepV<-rep(VecIn, length(indexTC))
						RepV1<-rep(RepV, length(ArrayDateVector))
						dfdef<-cbind(RepV1, df)
						return(dfdef)
						})
 res1Q<-apply(res,1,function(x){
						print(x)
						meanV<-lapply(x, function(y){	
													# print(y$Vec1)
													outMean<-y[y$Vec1 == "1st Qu.",]
													return(outMean)
													})
						df <- do.call(rbind, meanV)
						VecIn<- c("i", indexTC)
						RepV<-rep(VecIn, length(indexTC))
						RepV1<-rep(RepV, length(ArrayDateVector))
						dfdef<-cbind(RepV1, df)
						return(dfdef)
						})
 res3Q<-apply(res,1,function(x){
						print(x)
						meanV<-lapply(x, function(y){	
													# print(y$Vec1)
													outMean<-y[y$Vec1 == "3rd Qu.",]
													return(outMean)
													})
						df <- do.call(rbind, meanV)
						VecIn<- c("i", indexTC)
						RepV<-rep(VecIn, length(indexTC))
						RepV1<-rep(RepV, length(ArrayDateVector))
						dfdef<-cbind(RepV1, df)
						return(dfdef)
						})
 resMAX<-apply(res,1,function(x){
						print(x)
						meanV<-lapply(x, function(y){	
													# print(y$Vec1)
													outMean<-y[y$Vec1 == "Max.   ",]
													return(outMean)
													})
						df <- do.call(rbind, meanV)
						VecIn<- c("i", indexTC)
						RepV<-rep(VecIn, length(indexTC))
						RepV1<-rep(RepV, length(ArrayDateVector))
						dfdef<-cbind(RepV1, df)
						return(dfdef)
						})

lenV<-length(names(resMean))
DfDef<-data.frame()
for (i in c(1:lenV)){
									print(names(resMean)[i])
									tab1<-rbind(resMean[[i]], resMedian[[i]], res1Q[[i]], res3Q[[i]], resMin[[i]], resMAX[[i]])
									tabDef<-cbind(names(resMean)[i], tab1)
									print(tabDef)
									DfDef<-rbind(DfDef,tabDef) 
									print("###DfDef")
									print(DfDef)
									# return(DfDef)
									}

colnames(DfDef)[1] <- "ID"
return(DfDef)
}
# /*/ResampleRaster Sentinel to 10m:
# # input:
# 1) Raster (a risoluzioni 10,20,60)
# # ouput:
# raster a 10 m.
# ###################################################
ResSample10Raster<-function(x){
											inRaster<-raster(x)
											print(paste0("print Disaggregate RasterIN : ", inRaster))
											resraster <- raster::res(inRaster)[[1]]
											print(resraster)
											if (resraster == '10'){ 
														print("If dentro 10")
														result <- inRaster}
											else if(resraster == '20'){
														print("If dentro 20")
														result <- disaggregate(inRaster, fact = 2)}
											else if(resraster == '60'){
														print("If dentro 60")
														result <-disaggregate(inRaster, fact = 6)}
											return (result)
											print(paste0("print Disaggregate RasterOUT : ", result))
											# str(rextract)
											# Vector<-unlist(data.frame((strsplit(summary(rextract[2], na.rm=TRUE), ":")))[2,], use.names=FALSE)
											# vector1<-as.numeric(as.character(Vector))
											# return(vector1)
										}
# ###################################################	
# /*/thresholding the model ML output exclunding PIXEL with thrshold VALUe for NDVI e GNDVI
# # input:
# 1) Raster (a risoluzioni 10,20,60)
# # ouput:
# raster a 10 m.
# ###################################################													
Filteringoutput_thre<- function(ypred, v2input, v3input, threNDVI, threGNDVI){

								vectPRED<-as.numeric(as.character(ypred))
								# print(paste0("###y_pred as numeric: ", length(vectPRED)))
								# print(vectPRED)
								tabFilter<-cbind(vectPRED, v2input, v3input)
								# print("###tabFilter")
								# print(tabFilter)
								outputVector<-apply(tabFilter, 1, function(x){	
														# print("###inputValues")
														# print(x[1])
														# print(x[2])
														# print(x[3])
														if (x[3] < threNDVI|| is.na(x[3]) || is.na(x[2]) || x[2] < threGNDVI )
														{output<- 0 
														# print("##INSIDE IF")
														}
														else{output<-x[1]}
														# if (is.na(x[1]) || x[1] == 0 )
														# {output<- 0 }
														# else{output<- x[1]}
														# print("###function output")
														# print(output)
														# print("###END print output")
														return(output)
														})

								# print("###output vector output")
								# print(outputVector)
								# print("###END outputVector")
								return(outputVector)
}
# ###################################################
# ###################################################
# ###################################################
# # # # # # MAIN
# #read polygon to calculate values
# shapeOUTArray<-fileListShape(shapefolder)
# shapeOUTArrayFileName<-fileListShapeFileName(shapefolder)
# shapeOUTArray
# # read raster
# listFolder<-list.files(folderpath, full = TRUE)
# GetFileN<-GetRASTERFilePathANDArray(listFolder) # creo la l'array con un vettore di path dei file per ogni SUBFOLDER
# ArrayRaster<-GetFileN$GetRasterArray()# creo un array di raster per ogni folder
# ArrayRasterDate<-GetFileN$GetRasterDATEvector() # creo una lista di date dei rilievi
# ArrayDateVector<-as.Date(unlist(ArrayRasterDate)) # converto la lista in vettore
# ######INDEX of Satellite INDEX to CALC as in RStoolCalcINDEX
# # IndexToCalc<-c("NDWI2", "SATVI") 
# # IndexToCalc<-c("NDVI", "SAVI", "CTVI", "GNDVI", "SR", "MSAVI", "NDWI", "MNDWI", "NDWI2") 
# # IndexToCalc<-c("SR", "NRVI", "MSAVI")
# IndexToCalc<-c("NDVI", "SAVI", "SR", "NDWI2")
# # IndexToCalc<-c("NDVI", "NDWI2")
# xLimMax<-c(2)
# xLimMin<-c(-0.5)
# # ###################################################
# # # # RUN PROCESS
# CalcIndexForEachFolder<-lapply(ArrayRaster, function(x){
									# print("##RUN Cycle##")
									# CalcIndexVector<-CalcIndex(x, shapeOUTArray, IndexToCalc )
									# VIindexes<-CalcIndexVector$RStoolCalcINDEX()
									# return(VIindexes)
									# })#calcolo gli indici vegetazionali e maschero con SCL raster
# # ##Elaborazione output Funzione						
# res<-mapply(funT, CalcIndexForEachFolder, ArrayDateVector) #formato output
# rownames(res)<-shapeOUTArrayFileName
# colnames(res)<-as.Date(ArrayDateVector)					
# TabDef<-ResElab(res)									
# TabDef<-TabDef[TabDef$RepV1 != "i",]
# # TabDef<-TabDef[TabDef$RepV1 == "SR",]
# TabDef$Vec2<-as.numeric(as.character(TabDef$Vec2))
# head(TabDef)
# # ##elaborazione per PLOT
# TabCast<-as.data.frame(dcast(TabDef,  ID+RepV1+Date~Vec1, mean))
# colnames(TabCast)[4:9]<-c("q1", "q3", "maxV", "meanV", "medianV", "minV")
# by(TabCast,TabCast$ID, function(x){ 
						# PlotDBtot <- (ggplot(x, aes(as.Date(Date), medianV)) +
													# geom_point(aes(color = factor(RepV1)),size=2.5) +
													# geom_line(aes(color = factor(RepV1)),size=1.5)+
													# theme(legend.text=element_text(size=8))+
													# theme(legend.position = "bottom")+
													# ylim(-1,xLimMax))+
													# geom_errorbar(aes(ymin=q1, ymax=q3), width=.9, position=position_dodge(0.05))+ 
													# labs(title=x$ID,
													# x ="Date", y = "index Value")+
												    # geom_smooth(se=FALSE, linetype="dashed", size=0.5)
						# tiff(paste(pathOUTplot, x$ID, ".tif", sep=""), width = 2000, height = 1800, units = 'px', res=200)
						# print(PlotDBtot)
						# dev.off()
						# })
# ##elaborazione per writeCSV output
# # csvout<-write.csv(TabCast, "C:/Users/tonpo/source/repos/PythonSentinelDownload/plot/outTable1.csv")
# csvout<-write.csv(TabCast, pathOUTtable)
# pathOUTtable
# #clean TEMP folder
# do.call(file.remove, list(list.files(tempFolderClean, full.names = TRUE)))





