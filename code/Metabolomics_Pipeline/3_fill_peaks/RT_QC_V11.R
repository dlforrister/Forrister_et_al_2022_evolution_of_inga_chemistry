# load necessary packages
library(xcms)
library(stringr)
library(RMySQL)
library(plyr)
library(CAMERA)

# set working directory
setwd("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master")

# This code analyzes all retention time index (RTI) files in order to minimize machine variance in signal over time for actual samples

### load XCMS / CAMERA information ####
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

# load mzs for RTI standard
RTI_STD_NEG <- read.csv(file="./data/RTI_QC/RTI_STD_NEG.CSV",header=T)
names(RTI_STD_NEG) <- c("compound","mz","rt","ion_type")
RTI_STD_POS <- read.csv(file="./data/RTI_QC/RTI_STD_POS.CSV",header=T)
names(RTI_STD_POS) <- c("compound","mz","rt","ion_type")
RTI_STD_POLAR <- read.csv(file="./data/RTI_QC/RTI_POLAR_STD.CSV",header=T)
RTI_STD_NEG["mode"]="NEG"
RTI_STD_POS["mode"]="POS"
RTI_STD_POLAR["mode"]="POLAR"
RTI_STD_ALL <- rbind(RTI_STD_NEG, RTI_STD_POS, RTI_STD_POLAR)

# connect to DB
all_cons <- dbListConnections(MySQL())
for(con in all_cons){
  dbDisconnect(con)}

mydb <- dbConnect(MySQL(), user='u***', password='**********', dbname='RTI_DB', host='chpc.utah.edu')

# get RTI files from UPLC results that haven't been processed
all_standards <- dbGetQuery(mydbnew,"SELECT file_name, mode FROM `RTI_DB` WHERE Sample_Type = 'Standard' and ms_mode='MS' and file_name like 'R%'")

files_only <- all_standards[,"file_name"]

processed_qc <- unique(dbGetQuery(mydbold,"SELECT RTI FROM `RTI_DB`"))[,"RTI"]

to_run <- all_standards[!files_only %in% processed_qc,]
row.names(to_run) <- 1:nrow(to_run)


##### set parameters for XCMS functions (same as initial XCMS peak picking for samples) #####
dparam1 <- PeakDensityParam(sampleGroups = all_standards, bw=10, binSize=0.05, minSamples=1, minFraction = 0.01)
oparam <- ObiwarpParam(binSize=1)
dparam2 <- PeakDensityParam(sampleGroups = all_standards, bw = 3, binSize = 0.025, minSamples = 1, minFraction = 0.01)
fpparam <- FillChromPeaksParam(expandMz = 0.25, expandRt = 0.5)



for(n in 1:nrow(to_run)) { 

  if(startsWith(to_run[n, "mode"], "NEG")) {
  filepath <- paste("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/data/mzxml_paths/RTI_QC/C18_NEG/", to_run[n, "file_name"], ".mzXML", sep="")
  polarity <- "negative"
  RTI_STD <- RTI_STD_ALL[RTI_STD_ALL$mode == "NEG",]
  cwparam <- CentWaveParam(ppm=15, peakwidth=c(4,12), snthresh=5, prefilter=c(10,500))
}
if(startsWith(to_run[n, "mode"], "POS")) {
  filepath <- paste("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/data/mzxml_paths/RTI_QC/C18_POS/", to_run[n, "file_name"], ".mzXML", sep="")
  polarity <- "positive"
  RTI_STD <- RTI_STD_ALL[RTI_STD_ALL$mode == "POS",]
  cwparam <- CentWaveParam(ppm=15, peakwidth=c(4,12), snthresh=5, prefilter=c(10,500))
}
if(grepl("AA", to_run[n, "file_name"], fixed=TRUE)) {
    filepath <- paste("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/data/mzxml_paths/RTI_QC/polar/", to_run[n, "file_name"], ".mzXML", sep="")
    polarity <- "negative"
    RTI_STD <- RTI_STD_ALL[RTI_STD_ALL$mode == "POLAR",]
    cwparam <- CentWaveParam(ppm=15, peakwidth=c(2,24), snthresh=1.5, prefilter=c(5,250))
}

loop <- data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
loop$standard <-as.character("standard")
loop$ion_type <- as.character("ion_type")
loop$RTI <- as.character("RTI")
loop$file <- as.character("file.raw")
names(loop) <- c("mz","rt_in_min","TIC_into","TIC_intb","TIC_maxo","sn","pkwidth","ppm","ppm_abs","mda","mda_abs", "rerror", "rerror_abs","standard","ion_type","RTI","file")
loop <- loop[rep(row.names(loop),nrow(RTI_STD)),]
row.names(loop) <- (1:nrow(RTI_STD))
loop$standard <- RTI_STD$compound
loop$ion_type <- as.character(RTI_STD$ion_type)
loop_2 <- loop

if(!file.exists(filepath)) next

try(raw <- readMSData(filepath, msLevel. = 1, mode="onDisk")) # read each file and ignore errors
tryCatch(raw1 <- findChromPeaks(object = raw, param = cwparam), error = function(n) raw1 <- data.frame()) #find all peaks

if(length(raw1) != 0 & !is.null(raw1@msFeatureData$chromPeaks)){
raw2 <- data.frame(raw1@msFeatureData$chromPeaks) # save peak as dataframe
raw2$rt_in_min <- (raw2$rt)/60 #convert rt to rt_in_minutes
raw2$pkwdth <- (raw2$rtmax)-(raw2$rtmin)
raw3 <- raw2[ ,c("mz", "rt", "into", "sn", "rt_in_min", "pkwdth")] #extract col of interest

# for each compound in RTI_STD extract the mzfeature which minimizing ppm and rt error
for(i in 1:nrow(raw3)) {
  ppm <- 1000000*((raw3$mz[i])-(RTI_STD$mz))/(RTI_STD$mz)
  ppm_abs <- abs(ppm)
  mda <- ((raw3$mz[i])-(RTI_STD$mz))*1000
  mda_abs<- abs(mda)
  rerror <- RTI_STD$rt - raw3$rt_in_min[i]
  rerror_abs <- abs(rerror)
  
  candidates <- ppm_abs <= 100 & rerror_abs <= 1
  if(sum(candidates) == 0) next 
  
  terror <- ppm_abs[candidates] + 100*rerror_abs[candidates]# array  with sum of ppm + rt_shift for each feature in raw3 compared with standard of interest.  
  
  compound_match <- RTI_STD$compound[candidates][which.min(terror)]
  compound_idx <- which(loop_2$standard == compound_match)
  
  if(is.na(loop_2[compound_idx, "mz"])) {
    loop_2$mz[compound_idx] <- raw3$mz[i]
    loop_2$rt_in_min[compound_idx] <- raw3$rt_in_min[i]
    loop_2$TIC_into[compound_idx] <- raw3$into[i]
    loop_2$sn[compound_idx] <- raw3$sn[i]
    loop_2$pkwidth[compound_idx] <- raw3$pkwdth[i]
    loop_2$ppm[compound_idx] <- ppm[compound_idx]
    loop_2$ppm_abs[compound_idx] <- ppm_abs[compound_idx]
    loop_2$mda[compound_idx] <- mda[compound_idx]
    loop_2$mda_abs[compound_idx] <- mda_abs[compound_idx]
    loop_2$rerror[compound_idx] <- rerror[compound_idx]
    loop_2$rerror_abs[compound_idx] <- rerror_abs[compound_idx]
    loop_2$RTI[compound_idx] <- to_run[n, "file_name"]
  } else {
    terror1 <- loop_2$ppm_abs[compound_idx] + 100*loop_2$rerror_abs[compound_idx]
    if(min(terror) < terror1) {
      loop_2$mz[compound_idx] <- raw3$mz[i]
      loop_2$rt_in_min[compound_idx] <- raw3$rt_in_min[i]
      loop_2$TIC_into[compound_idx] <- raw3$into[i]
      loop_2$sn[compound_idx] <- raw3$sn[i]
      loop_2$pkwidth[compound_idx] <- raw3$pkwdth[i]
      loop_2$ppm[compound_idx] <- ppm[compound_idx]
      loop_2$ppm_abs[compound_idx] <- ppm_abs[compound_idx]
      loop_2$mda[compound_idx] <- mda[compound_idx]
      loop_2$mda_abs[compound_idx] <- mda_abs[compound_idx]
      loop_2$rerror[compound_idx] <- rerror[compound_idx]
      loop_2$rerror_abs[compound_idx] <- rerror_abs[compound_idx]
      loop_2$RTI[compound_idx] <- to_run[n, "file_name"]
    }
  }
}

if(all(is.na(loop_2$mz))){
raw6 <- loop  #if no peaks are found add one row with a "no peaks found note"
raw6$RTI <- to_run[n, "file_name"]
raw6$ion_type <-NA
raw6$standard <- "NO PEAKS MATCHED"
raw6[1,1] <- 1.0
batch <- as.data.frame(strsplit(as.character(to_run[n,"file_name"]), "_", fixed = T))

if(startsWith(to_run[n, "file_name"], "RN")) {
  raw6$Batch <- batch[1,] 
  raw6$injection <- batch[2,] 
  } else {
    raw6$Batch <- NA
    raw6$injection <- NA
  }

if(startsWith(to_run[n, "file_name"], "RN_AA")) {
  raw6$Batch <- batch[2,] 
  raw6$injection <- batch[3,] 
} else {
  raw6$Batch <- NA
  raw6$injection <- NA
}

if(startsWith(to_run[n, "file_name"], "RP")) {
  raw6$Batch <- batch[1,] 
  raw6$injection <- batch[2,] 
} else {
  raw6$Batch <- NA
  raw6$injection <- NA
}

raw6 <- raw6[,c(16,18,19,14,15,1,9,8,11,10,2,13,12,3,4,5,6,7,17)]
raw7<-raw6[!is.na(raw6$mz),]  
} else {
batch <- as.data.frame(strsplit(as.character(to_run[n,"file_name"]), "_", fixed = T))

if(startsWith(to_run[n, "file_name"], "RN")){
  loop_2$Batch <- batch[1,] 
  loop_2$injection <- batch[2,] 
  } else {
  loop_2$Batch <- NA
  loop_2$injection <- NA
  }

if(startsWith(to_run[n, "file_name"], "RN_AA")) {
  loop_2$Batch <- batch[2,] 
  loop_2$injection <- batch[3,] 
} else {
  loop_2$Batch <- NA
  loop_2$injection <- NA
}

if(startsWith(to_run[n, "file_name"], "RP")){
  loop_2$Batch <- batch[1,] 
  loop_2$injection <- batch[2,] 
} else {
  loop_2$Batch <- NA
  loop_2$injection <- NA
}

raw6 <- loop_2[c(16,18,19,14,15,1,9,8,11,10,2,13,12,3,4,5,6,7,17)] # raw6 needs to match database....
raw7 <- raw6[!is.na(raw6$mz),]
}} else {
  
  raw6<- loop  #if no peaks are found add one row with a "no peaks found note"
raw6$RTI <- to_run[n, "file_name"]
raw6$ion_type <-NA
raw6$standard <- "NO PEAKS MATCHED"
raw6[1,1] <- 1.0
batch<- as.data.frame(strsplit(as.character(to_run[n,"file_name"]), "_", fixed = T))

if(startsWith(to_run[n, "file_name"], "RN")) {
  raw6$Batch <- batch[1,] 
  raw6$injection <- batch[2,] 
} else {
  raw6$Batch <- NA
  raw6$injection <- NA
}

if(startsWith(to_run[n, "file_name"], "RN_AA")) {
  raw6$Batch <- batch[2,] 
  raw6$injection <- batch[3,] 
} else {
  raw6$Batch <- NA
  raw6$injection <- NA
}

if(startsWith(to_run[n, "file_name"], "RP")) {
  raw6$Batch <- batch[1,] 
  raw6$injection <- batch[2,] 
} else {
  raw6$Batch <- NA
  raw6$injection <- NA
}

raw6<- raw6[,c(16,18,19,14,15,1,9,8,11,10,2,13,12,3,4,5,6,7,17)]
raw7 <- raw6[!is.na(raw6$mz),] 
raw7$mz <- NA
}

RTI <- to_run[n, "file_name"]
File <- raw7$file[1]

if(startsWith(to_run[n, "file_name"], "RN_AA")){
  Batch <- batch[2,]
  injection <- batch[3,]
} else {
  Batch <- batch[1,]
  injection <- batch[2,]
}

npeaks <- nrow(raw7)
avg_ppm_abs <- mean(raw7$ppm_abs)
avg_ppm <- mean(raw7$ppm)
avg_rerror_abs <- mean(raw7$rerror_abs)
avg_rerror <- mean(raw7$rerror)
avg_TIC_int <- mean(raw7$TIC_into)
avg_SN <- mean(raw7$sn)
avg_pkwidth <- mean(raw7$pkwidth)

Cat_CO2 <- if(length(which(raw7$standard == "epicatechin[M-H-C02]"))>0 & length(which(raw7$standard ==  "epicatichan[M-H]"))>0) {((raw7$TIC_into[which(raw7$standard == "epicatechin[M-H-C02]")]) / (raw7$TIC_into[which(raw7$standard == "epicatichan[M-H]")]))*100} else {Cat_CO2 <- NA}

Cat_1 <- if(length(which(raw7$standard == "epicatechin[insource-1]"))>0 & length(which(raw7$standard == "epicatichan[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "epicatechin[insource-1]")] / (raw7$TIC_into[which(raw7$standard == "epicatichan[M-H]")]))*100 } else { Cat_1 <- NA}

Cat_2 <- if(length(which(raw7$standard == "epicatechin[insource-2]"))>0 & length(which(raw7$standard == "epicatichan[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "epicatechin[insource-2]")] / (raw7$TIC_into[which(raw7$standard == "epicatichan[M-H]")]))*100 } else {Cat_2 <- NA}

Cat_3 <- if(length(which(raw7$standard == "epicatechin[insource-3]"))>0 & length(which(raw7$standard == "epicatichan[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "epicatechin[insource-3]")] / (raw7$TIC_into[which(raw7$standard == "epicatichan[M-H]")]))*100 } else {Cat_3 <- NA}

Mor_1 <- if(length(which(raw7$standard == "Morin[insource-1]"))>0 & length(which(raw7$standard == "Morin[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "Morin[insource-1]")] / (raw7$TIC_into[which(raw7$standard == "Morin[M-H]")]))*100 } else {Mor_1 <- NA}

Mor_2 <- if(length(which(raw7$standard == "Morin[insource-2]"))>0 & length(which(raw7$standard == "Morin[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "Morin[insource-2]")] / (raw7$TIC_into[which(raw7$standard == "Morin[M-H]")]))*100 } else {Mor_2 <- NA}

Mor_3 <- if(length(which(raw7$standard == "Morin[insource-3]"))>0 & length(which(raw7$standard == "Morin[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "Morin[insource-3]")] / (raw7$TIC_into[which(raw7$standard == "Morin[M-H]")]))*100 } else {Mor_3 <- NA}


Mor_4 <- if(length(which(raw7$standard == "Morin[insource-4]"))>0 & length(which(raw7$standard == "Morin[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "Morin[insource-4]")] / (raw7$TIC_into[which(raw7$standard == "Morin[M-H]")]))*100} else{Mor_4 <- NA}


Trypt_1 <- if(length(which(raw7$standard == "Tryptophan[insource]"))>0 & length(which(raw7$standard == "Tryptophan[M-H]"))>0){(raw7$TIC_into[which(raw7$standard == "Tryptophan[insource]")] / (raw7$TIC_into[which(raw7$standard == "Tryptophan[M-H]")]))*100 } else{ Trypt_1<- NA}


Current_RTI <-data.frame("RTI","File","Batch","injection","npeaks","avg_ppm_abs","avg_ppm","avg_rerror_abs","avg_rerror","avg_TIC_int","avg_SN","avg_pkwidth","Cat_CO2","Cat_1","Cat_2","Cat_3","Mor_1","Mor_2","Mor_3","Mor_4","Trypt_1")[1,]

if(raw7$standard[1] == "NO PEAKS MATCHED") {Current_RTI$NOTES <- "No data collected from RTI, bad sample or corrupt file"}




#NEXT STEP, get averages for all RTIs in processed RTI that match batch!! use said averages to generate rmarkdown graphs with line being the average of the batch and box plots for each M-H peak and(?) M-2H
#try to figure out RTI internal fragmentation? #ask anthony to help me with this!




fields <- "`RTI`, `batch`, `injection`, `standard`, `ion_type`, `mz`, `ppm_abs`, `ppm`, `mda_abs`, `mda`, `rt_in_min`, `rerror_abs`, `rerror`, `TIC_into`, `TIC_intb`, `TIC_maxo`, `sn`,`pkwidth`, 'file'"

dbWriteTable(mydb, 'RTI_DB', raw7 , field.types = fields, row.names = F, overwrite = FALSE, append = T)

}
