library(mzR)
library(MSnbase)

# function that extracts mz/rt/TIC data for peaks from mzXML files. Requires vector containing full or relative filepaths for .mzXML files of interest and MSlevel to extract data for
extract.peakdata <- function(files, MSlevel) {
  curr_header <- header(openMSfile(files[1]))
  msms_scans <- curr_header[curr_header$peaksCount >= 5 & curr_header$totIonCurrent >= 2000 & curr_header$basePeakIntensity >=1000 & curr_header$msLevel %in% MSlevel,]
  if(nrow(msms_scans) > 0) msms_scans$file_idx <- 1
  else msms_scans$file_idx <- numeric(0)
  if(length(files) > 1) {
  for(j in 2:length(files)) {
    print(files[j])
    curr_header <- header(openMSfile(files[j]))
    print("1")
    curr_header$file_idx <- j
    print("2")
    msms_scans <- rbind(msms_scans, curr_header[curr_header$peaksCount >= 5 & curr_header$totIonCurrent >= 2000 & curr_header$basePeakIntensity >=1000 & curr_header$msLevel %in% MSlevel,])
    print("3")
    }
  }
  return(msms_scans)
}

# Function searches a retention time and mz window around the MSMS data (peakdata) for all of the MS features associated with the species (features_spec)
# returns a dataframe with any matched features, their compound number, and the peakdata ID
match_features <- function(x){
  msms_scan <- peakdata[peakdata$id == x,]
  mz <- msms_scan$precursorMZ
  mzmin <- mz - 0.01
  mzmax <- mz + 0.01
  rt <- msms_scan$retentionTime
  rtmin <- rt-60
  rtmax <- rt+60
  pot_feature <- features_spec[features_spec$mz >= mzmin & features_spec$mz <= mzmax & features_spec$rt*60 >= rtmin & features_spec$rt*60 <= rtmax,]
  if(length(pot_feature$feature_number) == 0){return(data.frame(id = peakdata$id[x],feature_number = "no_match",compound_number="no_match"))}
  if(length(unique(pot_feature$feature_number)) == 1 ){return(data.frame(id =peakdata$id[x],feature_number= as.character(unique(pot_feature$feature_number)),compound_number=(as.character(unique(pot_feature$compound_number)))))}
  if(nrow(pot_feature) >1){
    unique_features <- ddply(pot_feature,.(feature_number),summarize,mz=mean(mz),rt=mean(rt))
    closest_feature <- unique_features$feature_number[which.min(abs(unique_features$mz - mz)+abs(unique_features$rt*60 - rt))]
    pot_feature <- pot_feature[pot_feature$feature_number == as.character(closest_feature),]
    return(data.frame(id =peakdata$id[x],feature_number= as.character(unique(pot_feature$feature_number)),compound_number=(as.character(unique(pot_feature$compound_number)))))}
}

# Function finds weighted average of multiple msms scans. Requires msms_scans dataframe produced by extract.peakdata, vector containing file paths for MSMS .mzXML files (IN SAME ORDER AS USED FOR extract.peakdata function) as well as mz and rt of peak of interest.
avg.msms.spec <- function(files,feature,comp) {
  scans_to_merge <- peakdata[peakdata$id %in% ms_ms_match_results_df_matched[ms_ms_match_results_df_matched$feature_number == as.character(feature) & ms_ms_match_results_df_matched$compound_number == as.numeric(comp),]$id,c("acquisitionNum", "file_idx")]
  if(nrow(scans_to_merge) < 1) return(NULL)
  spec1 <- lapply(1:nrow(scans_to_merge), function(j) mzR::peaks(openMSfile(file.paths[scans_to_merge$file_idx[j]]),scans_to_merge$acquisitionNum[j]))
  spec2 <- do.call(rbind, spec1)
  spec3 <- spec2[order(spec2[,2], decreasing = T),]
  row.names(spec3) <- 1:nrow(spec3)
  rows_done = numeric()
  spec4 <- array(numeric(), dim=c(0,2))
  while(length(setdiff(1:nrow(spec3), rows_done)) > 0) {
    i = min(setdiff(1:nrow(spec3), rows_done))
    to_merge = setdiff(which(abs(spec3[,1] - spec3[i,1]) <= 0.25), rows_done)
    mz = weighted.mean(spec3[to_merge, 1], spec3[to_merge, 2])
    TIC = sum(spec3[to_merge, 2])
    spec4 <- rbind(spec4, c(mz, TIC))
    rows_done <- c(rows_done, to_merge)
  }
  spec5 <- spec4[order(spec4[,1]),]
  if((length(spec5[spec5[,2]>1000,])/2)<5) return(NULL)
  TICmax <- max(spec5[,2])
  spec6<- spec5[spec5[,2] >= TICmax/100 & spec5[,2] >= 500,,drop=FALSE]
  return(spec6)
}


