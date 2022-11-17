### This code is set up to run each species accession as a separate batch 
### using the slurm and shell scripts included in this folder

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
cat("  START TIME:", date(),"\n")
if (length(args)!=1){
  stop("Wrong arguments -> need 1!")
} else {
  SPECIES_ID = as.integer(args[1])
  cat("  SPECIES_ID:", SPECIES_ID,"\n",sep="")
}

# load necessary packages (suppressMessages to simplify reading slurm output files)
suppressMessages(library(xcms))
suppressMessages(library(splitstackshape))
suppressMessages(library(RMySQL))
suppressMessages(library(splitstackshape))
suppressMessages(library(data.table))
suppressMessages(library(plyr))

# load average_multiple_msms_scans function from separate file (included in this folder)
source("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/code/average_multiple_msms_scans.R")

# set working directory
setwd("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master")


####### Load in mzrt from pickle object #######
compound_feature_table <- read.csv("./data/C18/c18_comp_feat_pcid_mzrt_rtwindow.csv")

# remove any features with signal (TIC) below 1000
compound_feature_table_1 <- compound_feature_table[compound_feature_table$TIC > 1000,]
# use PCID column to extract species for each row
compound_feature_table_1$species <- vapply(strsplit(as.character(compound_feature_table_1$pcid),"_"),'[',1,FUN.VALUE = character(1))
comp_table_final <- compound_feature_table_1


####### Search DDA and MSMS files for msms spectra of each compound #######
# make vector of all species and skip any already matched to MSMS 
species <- unique(comp_table_final$species)
species_ran <- gsub(".mgf","",gsub("all_msms_spectra_","",list.files("./data/C18/spectral_hits/matched_compounds/",pattern = ".mgf")))
species <- species[!species %in% species_ran]

# load file with names and pathways of all MSMS files
all_files <- read.csv("./data/all_MSMS_files.csv")

i <- SPECIES_ID
print(as.character(species[i]))

# subset the compound feature table for the species of interest
features_spec <- comp_table_final[comp_table_final$species == species[i],]
# find all MSMS files associated with that species
files <- all_files[all_files$species_code == species[i], c("file_name","project_name")]
if(nrow(files) < 1) next

# make a vector of the filepaths for associated files
file.paths <- sapply(1:nrow(files), function(x) paste("./data/DDA_projects/",files$project_name[x],"/mzXML/Sample/",files$file_name[x],".mzXML", sep=""))

# use function from 'average_multiple_msms_scans.R' to get precursor ion information from all relevant files
peakdata <- extract.peakdata(file.paths, MSlevel = 2) 

# enumerate each unique peak 
peakdata$id <- 1:nrow(peakdata)

# use function from 'average_multiple_msms_scans.R' to identify potential MSMS ion matches to MS features
ms_ms_match_results <- lapply(peakdata$id,FUN = match_features)
ms_ms_match_results_df <- do.call("rbind", ms_ms_match_results)

# print the number of compounds matched to MSMS data out of the original number of compounds identified in MS
print(paste("matched ",length(unique(ms_ms_match_results_df$feature_number[ms_ms_match_results_df$compound_number != "no_match"])), " out of ",length(unique(features_spec$feature_number)), " features and ", length(unique(ms_ms_match_results_df$compound_number[ms_ms_match_results_df$compound_number != "no_match"])), " out of ", length(unique(features_spec$compound_number))," compounds",sep=""))

# subset the match result dataframe so only potential matches are included
ms_ms_match_results_df_matched <- ms_ms_match_results_df[ms_ms_match_results_df$feature_number != "no_match",]

# create unique data table of MS features for the species of interest
species_ms_ms_table <- data.table(unique(features_spec[,c("feature_number","compound_number","cluster_mz","cluster_rt")]))
# add column for average TIC per feature number
species_ms_ms_table <- merge(species_ms_ms_table,aggregate(features_spec$TIC,by = list(features_spec$feature_number,features_spec$compound_number),FUN = mean),by.x = c("feature_number","compound_number"), by.y = c("Group.1","Group.2"))
names(species_ms_ms_table)[5]<- "TIC"
# add columns for the UPLC column used to collect data, the current species, and whether or not the MS feature was matched to MSMS
species_ms_ms_table$column <- "C18"
species_ms_ms_table$species  <- species[i]
species_ms_ms_table$matched <- sapply(1:nrow(species_ms_ms_table), function(x) as.numeric(species_ms_ms_table$feature_number[x] %in% ms_ms_match_results_df_matched$feature_number))

# assess the number of matches per each compound (may be more than one based on low cosine threshold)  
match_count <- species_ms_ms_table[, sum(matched), by = compound_number]
names(match_count)[2]  <- "n_match_per_comp"
# add match count information to species table
species_ms_ms_table <- merge(species_ms_ms_table,match_count,by= c("compound_number"))

# make separate table of MS compounds NOT matched to MSMS spectra
species_ms_ms_table_no_match <- species_ms_ms_table[species_ms_ms_table$n_match_per_comp == 0,]
# find MS compound withOUT a match that has the highest TIC
species_ms_ms_table_no_match_top_1 <- setorder(setDT(species_ms_ms_table_no_match), -TIC)[, head(.SD, 1), keyby = compound_number]

# make separate table of MS compounds matched to MSMS spectra
species_ms_ms_table_matched <- species_ms_ms_table[species_ms_ms_table$n_match_per_comp > 0,]
# find the 5 MS compounds with a match that have the highest TIC
species_ms_ms_table_matched_top_5 <- setorder(setDT(species_ms_ms_table_matched), -TIC)[, head(.SD, 5), keyby = compound_number]
# combine the top 1 no match and top 5 with MSMS matches into one dataframe and assign each a scan number
species_ms_ms_table_to_extract_msms <- rbind(species_ms_ms_table_no_match_top_1,species_ms_ms_table_matched_top_5)
species_ms_ms_table_to_extract_msms$scan_number <- 1:nrow(species_ms_ms_table_to_extract_msms)
# convert to a list of vectors
msms_spec <- vector("list", length(species_ms_ms_table_to_extract_msms$scan_number))
 
print(paste("averaging spectra for ", species[i],sep=""))

# second round of MSMS spectra searching for compounds and features
for(y in 1:nrow(species_ms_ms_table_to_extract_msms)){
  # don't average spectra if none have MSMS matches  
  if(species_ms_ms_table_to_extract_msms$matched[y] == 0) next
  # pull out feature and compound number of interest
  feature <- species_ms_ms_table_to_extract_msms$feature_number[y]
  comp <- species_ms_ms_table_to_extract_msms$compound_number[y]
  # use average MSMS spectra function from 'average_multiple_msms_scans.R' to combine the top 5 MSMS and top 1 MS features found for a compound where spectra are available
  current_spec <- avg.msms.spec(file.paths,feature,comp)
  if(is.null(current_spec)) next # if current_spec is NULL, a spectrum was not collected for this compound
  if(nrow(current_spec) < 5) next # filter out if spectrum contains fewer than 5 peaks
  msms_spec[[species_ms_ms_table_to_extract_msms$scan_number[y]]] <- current_spec
  }

# label each compound feature as matched or no_spec if no spectra was found in the files
species_ms_ms_table_to_extract_msms$matched_actual <- "matched"
species_ms_ms_table_to_extract_msms$matched_actual[which(sapply(1:length(msms_spec), function(x) is.null(msms_spec[[x]])))] <- "no_spec"
comps_with_spec <- species_ms_ms_table_to_extract_msms$scan_number[species_ms_ms_table_to_extract_msms$matched_actual == "matched"]

####### Write .mgf containing MSMS spectrum for each compound #######
print(paste("writing mgf and scan info for ", species[i]," with ",table(species_ms_ms_table_to_extract_msms$matched_actual)[1]," spectra",sep = ""))
text_to_write <- c()
for(k in species_ms_ms_table_to_extract_msms$scan_number) {
  if(k %in% comps_with_spec) {
    text_to_write = c(text_to_write, "BEGIN IONS",paste("PEPMASS=",species_ms_ms_table_to_extract_msms$cluster_mz[k],sep=""),"CHARGE=1-",paste("SCANS=",species_ms_ms_table_to_extract_msms$scan_number[k],sep=""),sapply(1:nrow(msms_spec[[species_ms_ms_table_to_extract_msms$scan_number[k]]]), function(x) paste(msms_spec[[species_ms_ms_table_to_extract_msms$scan_number[k]]][x,1],msms_spec[[species_ms_ms_table_to_extract_msms$scan_number[k]]][x,2],sep="\t")),"END IONS")
  }
  else {
    text_to_write = c(text_to_write, "BEGIN IONS", "END IONS")
  }
}
write(text_to_write, file=paste("./data/C18/spectral_hits/matched_compounds/all_msms_spectra_", species[i], ".mgf", sep = ""))

# add .mgf filename information to species MSMS dataframe
species_ms_ms_table_to_extract_msms$mgf_filename <- paste("all_msms_spectra_", species[i], ".mgf", sep = "")

# write species MSMS dataframe to a CSV
write.csv(species_ms_ms_table_to_extract_msms[,c("species","scan_number","feature_number","compound_number","TIC","cluster_mz","cluster_rt","column","matched","matched","n_match_per_comp","matched_actual")],file=paste("./data/C18/spectral_hits/matched_compounds/csv/all_msms_spectra_info", species[i], ".csv", sep = ""),row.names = F)
