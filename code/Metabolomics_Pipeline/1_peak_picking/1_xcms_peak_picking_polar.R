####### Load required packages #######
library(multtest)
library(xcms)
library(vegan)
library(splitstackshape)
library(CAMERA)
library(reshape2)
library(RMySQL)

### set working directory ###
setwd("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master")

### load XCMS / CAMERA information ####
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

####### Load CAMERA file with list of possible adducts #######
rule_mod<-read.csv(file=paste("./data/extended_adducts_POS_NEG_modes",
                              list.files("./data/extended_adducts_POS_NEG_modes")[grep("current_neg",
                                                                                       list.files("./data/extended_adducts_POS_NEG_modes"))], sep="/"),
                   header= TRUE)

# Define absolute path to mzxml files, ending with folder containing folders named by site
mzxml_location <- "K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/mzxml_paths/polar/"

# set retention time cutoffs (in minutes) for peaks that should be kept in analysis (after 1 minute and before 10.5 minutes
# for polar column standard run)
min_rt <- 1.1
max_rt <- 10.5

# This code assumes converted .mzXML files are arranged into folders by species, with species being arranged in folders by site

# Within each species folder, there should be two folders, one named "Sample" containing all samples for that species
# And the other named "Blank" containing all blanks for that species
# ex 'K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/mzxml_paths/BCI/IngA/Sample/..mzXML' for extracted leaf samples
# and 'K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/mzxml_paths/BCI/IngA/Blank/..mzXML' for organic solvent blanks

# Create vector with all sites you wish to analyze
# For our data: BCI = BCI, Panama; Manaus = Manaus, Brazil; LA = Los Amigos, Peru; FG = Nouragues, French Guiana; Tiputini = Tiputini, Ecuador
sites <- c("BCI","Manaus","LA","FG","Tiputini")


# If you only want to analyze one or a few species, make a vector with the list of species and a vector of length 1 with the associated site. 
# If analyzing species from different sites, run them separately.
# Once setting these vectors, set i=1 and run all but the outermost loop.
# ex.
# sites <- "BCI"
# species <- c("IngA","IngC")



####### Loop through directories containing .mzXML files. Samples and their associated blanks should be in separate folders named "Sample" and "Blank". Depending on how files are arranged, you may or may not need to use the outer loop for site #######
for(i in 1:length(sites)) {
  # create list of species based on the folders present for each site in the mzxml directory
  species <- list.dirs(paste(mzxml_location,sites[i],sep=""), recursive=FALSE, full.names=FALSE)
for(j in 1:length(species)) {
  # keep track of which species is being processed in the R output
  print(species[j])
  
  # connect to MySQL database 
  all_cons <- dbListConnections(MySQL())
  for(con in all_cons){
    dbDisconnect(con)}
  mydb <- dbConnect(MySQL(), user='u***', password='**********', dbname='polar_XCMS_DB', host='chpc.utah.edu')
  
  # If some species have already been processed (i.e. j > 1 | i > 1), check database for processed species
  # This assumes that the database contains a "Species_code_sample" column containing "species" _ "sample_number"
  # ex. "IngA_1619"
  processed_sample <- unique(dbGetQuery(mydb, "SELECT Species_code_sample FROM `polar_XCMS_DB`"))[,"Species_code_sample"]
  processed <- sapply(processed_sample, function(x) {
    temp = unlist(strsplit(x, split = "_"))
    return(temp[1])
  })
  
  # Skip any species that have already been processed and are in the database
  if(species[j] %in% processed) next
  
  # make vector with relative filepaths (used by XCMS to find files)
  files_1 <- list.files(paste(mzxml_location, sites[i], "/",species[j],sep=""),
                        recursive = TRUE, full.names = TRUE)
  
  # make sure you are only using mzxml files, in case there are files of other types in your directory
  files <- files_1[endsWith(files_1, "mzXML")]
  
  # remove blank files from polar runs because no signal for this analysis
  files_polar <- files[!grepl("Blank",files)]
  
  # make vector with sample groups
  s_groups <- sapply(files_polar, function(x) {
    temp = unlist(strsplit(x, split = "/"))
    return(temp[length(temp)-1])
    })
  is.na(s_groups) <- 0
  
  # make vector with sample names (used by XCMS in peak list output)
  s_names <- sapply(files_polar, function(x) {
    temp = unlist(strsplit(x, split = "/"))
    return(sub(temp[length(temp)], pattern = ".mzXML", replacement = "", fixed = TRUE))
    })


  ##### set parameters for XCMS functions #####
  ## centwave is the method used to perform initial peak picking/finding -- optimized for amide column
  # ppm = minimum peak ppm for peak detection
  # peakwidth = upper and lower bound of peak width in seconds
  # snthresh = signal to noise  threshold
  # prefilter = (k,I) ; k = minimum number of centroids and I = minimum signal for peak detection and keeping
  
  cwparam <- CentWaveParam(ppm=15, peakwidth=c(0.2,5), snthresh=5, prefilter=c(1,500))
  
  
  ## peakdensity is the method used to group peaks across samples based on retention time. For the first peak grouping step, use more lenient parameters because retention time has not yet been corrected
  # sampleGroups = blank / sample group assignments from above
  # bw = bandwidth; standard deviation of the peak smoothing kernel
  # binSize = overlap of mz slices
  # minSamples = minimum number of samples for peak to be grouped / kept
  # minFraction = minimum fraction ofsamples in which peaks have to be present
  
  dparam1 <- PeakDensityParam(sampleGroups = s_groups, bw=5, binSize=0.05, minSamples=1, minFraction = 0.5)
  
  
  ## for polar files, use peak density is used for retention time correction
  # minFraction = minimum fraction of total samples a peak is present in for it to be aligned and RT to be adjusted
  # extraPeaks = maximum number of additional peaks that can be assigned to a peak group
  # smooth = function used to interpolate corrected RT
  # span = degree of smoothing for 'loess'
  
  # change smoothing pattern if only 1 file
  if (length(files_polar) < 2) {
  pgp <- PeakGroupsParam(minFraction= 0.5, extraPeaks = 100, smooth="linear", span = 0.6)
  } else {
  pgp <- PeakGroupsParam(minFraction= 0.5, extraPeaks = 100, smooth="loess", span = 0.6)
  }
  
  # second round of peak grouping uses more stringent parameters
  dparam2 <- PeakDensityParam(sampleGroups = s_groups, bw = 3, binSize = 0.025, minSamples = 1, minFraction = 0.01)

  ## FillChromPeaks look for all peaks in all samples to find peaks that were below threshold set during peak picking
  # expandMz = mz used to look on either side of peak to group / fill
  # expandRt = retention time (in min) used to look on either side of peak to group / fill
  
  fpparam <- FillChromPeaksParam(expandMz = 0.25, expandRt = 0.5)
  
  
  ## reads mzXML files into OnDiskMSnExp object
  # MS level one for initial peak picking
  raw_data <- readMSData(files, msLevel. = 1, mode="onDisk")
  
  ## Run XCMS functions with parameters set above
  xcmsexp <- findChromPeaks(object = raw_data, param = cwparam)
  xcmsexp <- groupChromPeaks(object = xcmsexp, param = dparam1)
  xcmsexp <- adjustRtime(object = xcmsexp, param = pgp)
  xcmsexp <- groupChromPeaks(object = xcmsexp, param = dparam2)
  xcmsexp <- fillChromPeaks(xcmsexp, param = fpparam)
  
  ## convert from OnDiskMSnExp object to xcmsset object so that CAMERA can be used
  ## to group adducts and isotopes
  xset <- as(xcmsexp, "xcmsSet")
  
  # set sample class to match sample grouping
  sampclass(xset) <- s_groups
  
  ## peak grouping and annotation
  # define polarity of MS run -- we use negative mode for our processing
  xset1 <- xsAnnotate(xs=xset, polarity="negative")
  
  # perfwhm = percentage of the full width at half maximum (FWHM)
  xset2 <- groupFWHM(xset1, perfwhm=0.7)
  
  # ppm = absolute error for ppm
  # mzabs = absolute error for mz
  # intval = general intensity (TIC) value used
  xset3 <- findIsotopes(xset2, ppm=20, mzabs=0.015,intval="intb")
  
  # cor_eic_th = correlation threshold (0-1)
  # pval = significant correlation threshold
  # graphMethod = method for grouping peaks after correlation into pseudospectra ; lpc = label propagation community algorithm
  # calcIso = use isotopic relationship for peak grouping
  # calcCis = use correlation inside samples for peak grouping
  # calcCas = use correlation across samples for peak grouping
  xset4 <- groupCorr(xset3, cor_eic_th=0.5, pval=0.5, graphMethod="lpc", calcIso = TRUE, calcCiS = TRUE, calcCaS = ifelse(
    sum(sampclass(xset) == "Sample") > 4, TRUE, FALSE))
  
  # find adducts using file loaded above, specify polarity
  xsetFA <- findAdducts(xset4, polarity="negative", rules = rule_mod)
  
  # find and group peaks based on isotope and adducts from above
  xset5 <- getPeaklist(xsetFA)
  
  #### END OF XCMS AND CAMERA CODE ####
  
  
  # Cleans up peak list by deleting internal standard and other known contaminants
  xset5[is.na(xset5)] <- 0
  
  # calculate average value of each peak in the samples
  xset5$TIC_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample")), drop=FALSE],na.rm=T)
  xset6 <- xset5
  
  ## remove raffinose, mz=91.008
  # raffinose
  if(nrow(xset6[which(abs(xset6$rt - 495) <= 20 & abs(xset6$mz-503.162) <=0.1),]) >0) {xset6 <- xset6[-which(abs(xset6$rt - 495) <= 20 & abs(xset6$mz-503.162) <=0.1),]}
  # mz=91.008
  if(nrow(xset6[which(abs(xset6$mz-91.008) <=0.01),]) > 0) {xset6 <- xset6[-which(abs(xset6$mz-91.008) <=0.01),]}
  
  
  # add columns with rounded MZ and rounded RT in minutes
  xset6$rt_in_min<- (xset6$rt)/60  
  xset6$mz_round <- round((xset6$mz),4)
  xset6$rt_round <- round((xset6$rt_in_min),4)
  
  # remove features before and after retention time cutoffs set above
  xset6 <- xset6[xset6$rt_in_min <= max_rt & xset6$rt_in_min >= min_rt, ]
  
  # change dataframe from wide to long format and rename columns
  xset6_long <- melt(data=xset6, id.vars=c("mz_round","rt_round","pcgroup"),measure.vars=c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample")))
  xset6_long$PC_ID <- paste(species[j], xset6_long$pcgroup,sep="_")
  xset6_long$Species_code_sample <- paste(species[j], gsub(".mzXML","",gsub("DN","", xset6_long$variable)), sep="_")
  feature_table_long <- xset6_long[,c(4,7,2,1,6,5)]
  names(feature_table_long) <- c("sample","Species_code_sample","RT","MZ","PC_ID","TIC")
  feature_table_long$TIC <- sapply(1:nrow(feature_table_long), function(x) round(feature_table_long$TIC[x], 4))
  fields <- names(feature_table_long)
  
  # reconnect to database to upload data for each species accession
  all_cons <- dbListConnections(MySQL())
  for(con in all_cons){
    dbDisconnect(con)}
  
  mydb <- dbConnect(MySQL(), user='u***', password='**********', dbname='polar_XCMS_DB', host='chpc.utah.edu')
  
  dbWriteTable(mydb, 'polar_xcms_feature_table_long', feature_table_long , field.types = fields, row.names = F, overwrite = FALSE, append = T)
}
}
