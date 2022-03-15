### This code is set up to run each sample as a separate batch 
### using the slurm and shell scripts included in this folder

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
cat("  START TIME:", date(),"\n")

if (length(args)!=1){
  stop("Wrong arguments -> need 1!")
} else {
  SPECIES_ID = as.integer(args[1])
  cat("  Sample_ID:", SPECIES_ID,"\n",sep="")
}

# load necessary packages (suppressMessages to simplify reading slurm output files)
suppressMessages(library(xcms))
suppressMessages(library(mzR))
suppressMessages(library(RMySQL))
suppressMessages(library(plyr))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(splitstackshape))
suppressMessages(library(dplyr))

# This code searches for all features identified in "2_group_features_compounds.py" in all samples
# Files should be in same folders used for "1_xcms_peak_picking.R"

# set working directory
setwd("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master")

####### load feature table built in '2_group_features_compounds.py' #######
features <- read.csv("./data/polar/polar_compound_feature_table.csv",stringsAsFactors=FALSE)

# make sure there are no duplicate features in the table
features <- unique(features[,c("feature_number","mz","rt")])

###### load in pkl table and get feature min and max info to set tolerance ###########
pkl <- read.csv("./data/polar/polar_comp_feat_pcid_mzrt_rtwindow.csv",row.names = 1)
# turn dataframe into data table
pkl_table <- data.table(pkl)
# summarize each feature by the min, max, and mean mz / rt, and the cluster mz and rt
feature_mzrt <- ddply(pkl_table,.(feature_number),summarize,cluster_mz=median(cluster_mz),mean_mz= mean(mz),min_mz=min(mz),max_mz=max(mz),cluster_rt=median(cluster_rt),mean_rt = mean(rt),min_rt=min(rt),max_rt=max(rt))


###### load in file with sample RTI and blank associations (from same day UPLC runs) ######
sb <- read.csv("./data/polar/amide_column_samples_quality_reviewed.csv")

# keep only sample files that have been deemed good quality
sb_row <- sb[sb$quality != "drop",]
# subset and rename desired columns
sample_blank_RTI <- sb_row[, c("file_name","blank","RTI","mzxml","species","site")]
colnames(sample_blank_RTI) <- c("file_name","asoc_blank","asoc_RTI","mzXML","species","site")

####### set TIC cutoff to use for finding features ######
tic_cutoff <- 1000

# Define absolute path to mzxml files, ending with folder containing folders named by site
# (the format of these folders are further explained in step 1 of this code - xcms peak picking)
mzxml_location <- "K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/mzxml_paths/polar/"

###### Create vector with all sites you wish to analyze
sites <- c("BCI","Manaus","LA","FG","Tiputini")

####### Load .mzXML file for each sample and search for all features grouped in step 2
####### using fillPeaks function from package xcms
# name FULL OR RELATIVE filepath for output file, including csv extention
filled_features_output_path <- paste("./data/polar/polar_fill_feat/",remaining_samples$species[n],"_",gsub(".mzXML","",basename(remaining_samples$mzxml[n])),"_filled_features_polar.csv",sep="")

# make sure sample hasn't already been run
species_ran <- gsub("_filled_features_polar.csv","",list.files("./data/polar/polar_fill_feat/"))
sample_blank_RTI$sample_name <- paste(sample_blank_RTI$species,gsub(".mzXML","",basename(sample_blank_RTI$mzxml)),sep = "_")
sample_blank_RTI$ran <- 0
sample_blank_RTI$ran[sample_blank_RTI$sample_name %in% species_ran] <- 1
remaining_samples <- sample_blank_RTI[sample_blank_RTI$ran == 0,]

# begin XCMS analysis
n = SPECIES_ID

print(remaining_samples$sample_name[n])

# make sure sample is labeled "keep" in quality review file
just_samp <- sapply(remaining_samples$mzxml[n], function(x) {
  temp = unlist(strsplit(x, split = "/"))
  return(sub(temp[length(temp)], pattern = ".mzXML", replacement = "", fixed = TRUE))
  })

if (just_samp %in% sample_blank_RTI$file_name){
  #get associated RTI (retention time index) data from DB
  asoc_rti <- as.character(sample_blank_RTI[sample_blank_RTI$file_name == just_samp, "asoc_RTI"])

  # connect to MySQL database 
  all_cons <- dbListConnections(MySQL())
  for(con in all_cons){
    dbDisconnect(con)}
  mydb <- dbConnect(MySQL(), user='u***', password='**********', dbname='XCMS_DB', host='chpc.utah.edu')
  
  # pull out RTI TIC to normalize sample TIC
  rti_tic <- dbGetQuery(mydbold, paste("SELECT RTI, Standard, TIC_into from `RTI_DB` where RTI like '%polar%'", sep = ""))
  rti_tic <-   rti_tic[rti_tic$RTI == asoc_rti,]         
  avg_rti_tic <- mean(rti_tic$TIC_into)
  
  # fill peaks for current sample
  sample_xcms <- xcmsRaw(paste(mzxml_location, remaining_samples$site[n], "/", remaining_samples$species[n], "/Sample/", remaining_samples$mzxml[n], sep = ""), profstep = 1, profmethod = "bin", includeMSn = FALSE, mslevel = 1)
  
  # create empty dataframe to append filled peaks to
  sample_peaks <- data.frame("feature_number"=numeric(0),
                                "cluster_mz"=numeric(0),
                                "cluster_rt"=numeric(0),
                                "into"=numeric(0),
                                "filled_tic"=numeric(0),
                                "norm_tic"=numeric(0),
                                "real_mz"=numeric(0),
                                "real_rt"=numeric(0))
  
  # search through sample for ALL features in feature table
  for (i in 1:nrow(features)){
    feature_number <- features$feature_number[i]
    # search mass window +/- 5 mDa from the max and min mz of given feature
    mzmin <- feature_mzrt$min_mz[feature_mzrt$feature_number==feature_number]-0.005
    mzmax <- feature_mzrt$max_mz[feature_mzrt$feature_number==feature_number]+0.005
    # search RT window +/- 30 seconds from the max and min RT of given feature
    rtmin <- 60*(feature_mzrt$min_rt[feature_mzrt$feature_number==feature_number]-0.5)
    rtmax <- 60*(feature_mzrt$max_rt[feature_mzrt$feature_number==feature_number]+0.5)
    # find the scan times in the actual sample corresponding to rtmin and rtmax
    scanmin <- min(which(sample_xcms@scantime >= rtmin))
    scanmax <- max(which(sample_xcms@scantime <= rtmax))
    # fill peak (integrate under chromatographic curve) in the specified mz, RT, and scan ranges
    int_scan <- rawEIC(sample_xcms, mzrange = c(mzmin,mzmax), rtrange = c(rtmin,rtmax),scanrange=c(scanmin,scanmax))
    # set the midpoint of the peak at the point of highest signal / intensity
    midpoint <- int_scan$scan[which.max(int_scan$intensity)]
    
    # evaluate the signal to determine if this peak exists in the sample
    # is the signal > 0 for the entirety of the mz, RT, and scan ranges?
    if (sum(int_scan$intensity > 0) > 0) {
      # is the midpoint well away from the end of the mz, RT, and scan ranges?
      if (!(midpoint %in% int_scan$scan[1:5]) && !(midpoint %in% int_scan$scan[(length(int_scan$scan)-5):length(int_scan$scan)])){
        # is the sum of the signal over the scan range > 1000 (our TIC cutoff)?
        if (sum(int_scan$intensity[(which(int_scan$scan==midpoint)-5):(which(int_scan$scan==midpoint)+5)]) > tic_cutoff){
          # if all of these conditions are true, proceed to fill the peak
          # calculate the mz and rt of the peak midpoint
          scan_mz <- getScan(sample_xcms,midpoint,mzrange=c(mzmin,mzmax))
          scan_rt <- sample_xcms@scantime[midpoint]
          # get percentage of each scan compared to the max scan (midpoint) within the mz/rt range
          int_pct <- int_scan$intensity / int_scan$intensity[which.max(int_scan$intensity)]
          # Smooth peak by averaging each scan with the 5 scans on either side
          int_5avg <- sapply(6:(length(int_pct)-6), function(x) mean(int_pct[(x-5):(x+5)]))
          # if any averaged scan values are greater than 20% of maximum peak TIC
          if (any(int_5avg > 0.2)){
            # find first and last scans of the smoothed peak
            first_scan <- int_scan$scan[min(which(int_5avg > 0.2))+5]
            last_scan <- int_scan$scan[max(which(int_5avg > 0.2))+5]
            # append all desired data to the sample_peaks dataframe
            temp_peaks <- data.frame("feature_number"=numeric(1),
                                     "cluster_mz"=numeric(1),
                                     "cluster_rt"=numeric(1),
                                     "into"=numeric(1),
                                     "filled_tic"=numeric(1),
                                     "norm_tic"=numeric(1),
                                     "real_mz"=numeric(1),
                                     "real_rt"=numeric(1))
            temp_peaks$feature_number <- feature_number
            temp_peaks$cluster_mz <- unique(features[features$feature_number == feature_number, "mz"])
            temp_peaks$cluster_rt <- unique(features[features$feature_number == feature_number, "rt"])
            temp_peaks$into <- int_scan$intensity[int_scan$scan==midpoint]
            temp_peaks$filled_tic <- sum(int_scan$intensity[which(int_scan$scan==first_scan):which(int_scan$scan==last_scan)])
            temp_peaks$norm_tic <- (temp_peaks$filled_tic / avg_rti_tic) * 1000000
            temp_peaks$real_mz <- scan_mz[1,1]
            temp_peaks$real_rt <- scan_rt/60
            sample_peaks <- rbind(temp_peaks,sample_peaks)
          }
        }
      }
    }
  }
  
  # rename columns to match output
  names(sample_peaks) <- c("feature_number","feature_mz","feature_rt", "into", "TIC","norm_TIC","actual_mz", "actual_rt")
  # remove any duplicate peaks
  sample_peaks_2 <- distinct(sample_peaks)
      
      
  # get organic solvent blank associated with the sample (run on the same day)
  sample_name <- remaining_samples$file_name[n]
  blank_name <- as.character(sample_blank_RTI[which(sample_blank_RTI[,"file_name"]==sample_name), "asoc_blank"])

  # use same process to check associated blank for all features
  blank_xcms <- xcmsRaw(paste(mzxml_location, remaining_samples$site[n], "/", remaining_samples$species[n], "/Blank/", blank_name, ".mzXML", sep = ""), profstep = 1, profmethod = "bin", includeMSn = FALSE, mslevel = 1)
  
  # find median TIC in blank run and set this as a second cutoff for peak filling
  # any peaks found in initial peak picking that fall below this cutoff are considered 'in the noise' of the UPLC-MS run
  median_tic <- median(blank_xcms@tic)
  tic_cutoff_2 <- median_tic
  
  # create empty dataframe to append filled peaks to
  blank_findpeaks <- data.frame("feature_number"=numeric(0),
                                 "cluster_mz"=numeric(0),
                                 "cluster_rt"=numeric(0),
                                 "into"=numeric(0),
                                 "filled_tic"=numeric(0),
                                 "norm_tic"=numeric(0),
                                 "real_mz"=numeric(0),
                                 "real_rt"=numeric(0))
  
  # search through blank for ALL features in feature table  
  for (i in 1:nrow(features)){
    feature_number <- features$feature_number[i]
    # search mass window +/- 5 mDa from the max and min mz of given feature
    mzmin <- feature_mzrt$min_mz[feature_mzrt$feature_number==feature_number]-0.005
    mzmax <- feature_mzrt$max_mz[feature_mzrt$feature_number==feature_number]+0.005
    # search RT window +/- 30 seconds from the max and min RT of given feature
    rtmin <- 60*(feature_mzrt$min_rt[feature_mzrt$feature_number==feature_number]-0.5)
    rtmax <- 60*(feature_mzrt$max_rt[feature_mzrt$feature_number==feature_number]+0.5)
    # find the scan times in the actual sample corresponding to rtmin and rtmax
    scanmin <- min(which(sample_xcms@scantime >= rtmin))
    scanmax <- max(which(sample_xcms@scantime <= rtmax))
    # fill peak (integrate under chromatographic curve) in the specified mz, RT, and scan ranges
    int_scan <- rawEIC(sample_xcms, mzrange = c(mzmin,mzmax), rtrange = c(rtmin,rtmax),scanrange=c(scanmin,scanmax))
    # set the midpoint of the peak at the point of highest signal / intensity
    midpoint <- int_scan$scan[which.max(int_scan$intensity)]
    
    # evaluate the signal to determine if this peak exists in the sample
    # is the signal > 0 for the entirety of the mz, RT, and scan ranges?
    if (sum(int_scan$intensity > 0) > 0) {
      # is the midpoint well away from the end of the mz, RT, and scan ranges?
      if (!(midpoint %in% int_scan$scan[1:5]) && !(midpoint %in% int_scan$scan[(length(int_scan$scan)-5):length(int_scan$scan)])){
        # is the sum of the signal over the scan range > 1000 (our TIC cutoff)?
        if (sum(int_scan$intensity[(which(int_scan$scan==midpoint)-5):(which(int_scan$scan==midpoint)+5)]) > tic_cutoff){
          # if all of these conditions are true, proceed to fill the peak
          # calculate the mz and rt of the peak midpoint
          scan_mz <- getScan(sample_xcms,midpoint,mzrange=c(mzmin,mzmax))
          scan_rt <- sample_xcms@scantime[midpoint]
          # get percentage of each scan compared to the max scan (midpoint) within the mz/rt range
          int_pct <- int_scan$intensity / int_scan$intensity[which.max(int_scan$intensity)]
          # Smooth peak by averaging each scan with the 5 scans on either side
          int_5avg <- sapply(6:(length(int_pct)-6), function(x) mean(int_pct[(x-5):(x+5)]))
          # if any averaged scan values are greater than 20% of maximum peak TIC
          if (any(int_5avg > 0.2)){
            # find first and last scans of the smoothed peak
            first_scan <- int_scan$scan[min(which(int_5avg > 0.2))+5]
            last_scan <- int_scan$scan[max(which(int_5avg > 0.2))+5]
            # append all desired data to the blank_findpeaks dataframe
            temp_peaks <- data.frame("feature_number"=numeric(1),
                                     "cluster_mz"=numeric(1),
                                     "cluster_rt"=numeric(1),
                                     "into"=numeric(1),
                                     "filled_tic"=numeric(1),
                                     "norm_tic"=numeric(1),
                                     "real_mz"=numeric(1),
                                     "real_rt"=numeric(1))
            temp_peaks$feature_number <- feature_number
            temp_peaks$cluster_mz <- unique(features[features$feature_number == feature_number, "mz"])
            temp_peaks$cluster_rt <- unique(features[features$feature_number == feature_number, "rt"])
            temp_peaks$into <- int_scan$intensity[int_scan$scan==midpoint]
            temp_peaks$filled_tic <- sum(int_scan$intensity[which(int_scan$scan==first_scan):which(int_scan$scan==last_scan)])
            temp_peaks$norm_tic <- (temp_peaks$filled_tic / avg_rti_tic) * 1000000
            temp_peaks$real_mz <- scan_mz[1,1]
            temp_peaks$real_rt <- scan_rt/60
            sample_peaks <- rbind(temp_peaks,blank_findpeaks)
          }
        }
      }
    }
  }
  
  # rename columns to match output
  names(blank_findpeaks) <- c("feature_number","feature_mz","feature_rt", "into", "TIC","norm_TIC","actual_mz", "actual_rt")
  # remove any duplicate features
  blank_findpeaks_2 <- distinct(blank_findpeaks)

  # remove all peaks that are not at least 5x as abundant in sample as in blank.
  sample_blank_peaks <- merge(sample_peaks_2, blank_findpeaks_2[,c("feature_number", "TIC")], by = "feature_number", all.x = TRUE, all.y = FALSE)
  sample_blank_peaks$TIC.y[is.na(sample_blank_peaks$TIC.y)] <- 0
  sample_peaks <- sample_blank_peaks[(sample_blank_peaks$TIC.x / sample_blank_peaks$TIC.y) > 5, ]
  sample_peaks_3 <- sample_peaks
      
  # in addition, remove all peaks that are not at least as abundant as the median blank intensity
  sample_peaks_3 <- sample_peaks_3[sample_peaks_3$TIC.x > tic_cutoff_2,]
      

  # remove duplicate mz_rt instances (different feature identified as same peak) and keep one with lower rt error compared to feature_rt
  sample_peaks$mzrt <- paste(sample_peaks$actual_mz, sample_peaks$actual_rt, sep="_")
  for(i in unique(sample_peaks$mzrt)) {
    curr_mzrt <- sample_peaks[sample_peaks$mzrt == i, ]
    if(nrow(unique(curr_mzrt)) > 1) {
      to_remove <- curr_mzrt$feature_number[-which.min(abs(curr_mzrt$actual_rt - curr_mzrt$feature_rt))]
      sample_peaks_3 <- sample_peaks_3[!sample_peaks_3$feature_number %in% to_remove, ]
    }
  }
  
  # put sample_peaks dataframe in the desired format
  sample_peaks_4 <- data.frame("feature_id" = sample_peaks_3$feature_number, "TIC" = sample_peaks_3$TIC.x, "norm_TIC" = sample_peaks_3$norm_TIC, "actual_mz" = sample_peaks_3$actual_mz, "actual_rt" = sample_peaks_3$actual_rt, "sample_name" = paste(remaining_samples_to_run$species[n], "_",unlist(strsplit(sample_name, split = "[.]"))[1],sep=""))
  # remove any NA feature_id
  sample_peaks_4 <- sample_peaks_4[!is.na(sample_peaks_4$feature_id),]
  
  ####### write results to CSV #######
  write.table(sample_peaks_4, filled_features_output_path, sep = ",", append = FALSE, row.names = FALSE, col.names = TRUE)
  
  }

# After each sample has been run and an output CSV created, you can append them into one long file using the code below:
## List filenames to be merged. 
filenames <- list.files(path="./data/polar/polar_fill_feat/",pattern="*.csv")

## Full path to csv filenames
fullpath=file.path("./data/polar/polar_fill_feat/",filenames)

## Merge listed files from the path above
dataset <- do.call("rbind",lapply(fullpath,FUN=function(files){ read.csv(files)}))

## write new CSV 
write.csv(dataset,"./data/polar/polar_fill_feat/filled_features_polar_allsamples.csv",row.names=F)