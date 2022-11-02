########### This code is written to combine the data from both polar and C18 columns into one chemical similarity metric at the accession level ############
####### Load required packages #######
library(foreach)
library(doParallel)
library(plyr)
library(data.table)
library(RMySQL)

### set working directory ###
setwd("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/")

# load functions from separate create_pairwisecomps..R and chem_similarity_function.R files (included in this folder)
source("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/code/6_chemical_similarity/create_pairwiseComps_sampsByComps.R")
source("K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/code/6_chemical_similarity/chem_similarity_function.R")

####### Download 'view raw spectra' from GNPS results page and use make_pairwisecomps function from 'create_pairwiseComps_sampsByComps.R' to creating pairwise compound similarity matrix #######
all.pairwise.comps <- make_pairwisecomps("./data/GNPS_output/", method = "combined", through_network=TRUE)

####### Load filled compound tables created in 'fill_compounds.py' #######
polar_fill_comp <- read.csv("./data/polar/filled_compounds_polar.csv", stringsAsFactors = F)

c18_fill_comp <- read.csv("./data/C18/filled_compounds_c18.csv",stringsAsFactors = F)

## load in all compound info (cluster index / unique idenifiers) made from GNPS output and convert compound number to numeric
compound_info <- read.csv("./data/all_comps_info.csv")
compound_info$compound_number <- as.numeric(as.character(compound_info$compound_number))


#### Transform fill compounds dataframes to data tables and merge with compound info from GNPS
# C18
c18_fill_comps_table <- data.table(c18_fill_comp)
c18_fill_comps_table$column <- "C18"
## convert compound number to numeric
c18_fill_comps_table$compound_number <- as.numeric(as.character(c18_fill_comps_table$compound_number))
c18_fill_comps_table <- merge(c18_fill_comps_table,compound_info[,c("cluster.index","compound_number","column")],by=c("compound_number","column"))
## pull out species accession from the sample label
c18_fill_comps_table$species <- vapply(strsplit(as.character(c18_fill_comps_table$sample),"_"),'[',1,FUN.VALUE = character(1))

# polar
polar_fill_comps_table <- data.table(polar_fill_comp)
polar_fill_comps_table$column <- "polar"
## convert compound number to numeric
polar_fill_comps_table$compound_number <- as.numeric(as.character(polar_fill_comps_table$compound_number))
polar_fill_comps_table <- merge(polar_fill_comps_table,compound_info[,c("cluster.index","compound_number","column")],by=c("compound_number","column"))
## pull out species accession from the sample label
polar_fill_comps_table$species <- vapply(strsplit(as.character(polar_fill_comps_table$sample),"_"),'[',1,FUN.VALUE = character(1))


#### Run the code below only if you would like to exclude compounds not found in at least 40% of samples per species accession (excludes rare or low accumulation compounds) ####
#cutoff <- 0.4

# C18
## use ddply to count how many total samples are included in the dataset per species accession
#c18_species_count <- unique(ddply(c18_fill_comps_table, .(species), mutate, count = length(unique(sample)))[,c("species","count")])

## count how many samples of the total per species accession each compound appears in
#c18_compound_counts <- merge(c18_fill_comps_table[, list(Freq =.N), by=list(cluster.index,species)],c18_species_count,by = "species",all.x=T)

## calculate the percentage of total samples each compound appears in per species accession
#c18_compound_counts$perc <- round(c18_compound_counts$Freq/c18_compound_counts$count,2)

## label those above your predetermined cutoff as compounds to keep for that accession
#c18_compound_counts$keep <- c18_compound_counts$perc >= cutoff

## merge the count table with the fill compounds table and remove all that were not labeled to keep (not above the cutoff)
#c18_fill_comps_table_1 <- merge(c18_fill_comps_table,c18_compound_counts[,c("species","cluster.index","keep")], by = c("species","cluster.index"))
#c18_fill_comps_table_2 <- c18_fill_comps_table_1[c18_fill_comps_table_1$keep,-c("keep")]
#c18_fill_comp <- c18_fill_comps_table_2

# polar
## use ddply to count how many total samples are included in the dataset per species accession
#polar_species_count <- unique(ddply(polar_fill_comps_table, .(species), mutate, count = length(unique(sample)))[,c("species","count")])

## count how many samples of the total per species accession each compound appears in
#polar_compound_counts <- merge(polar_fill_comps_table[, list(Freq =.N), by=list(cluster.index,species)],polar_species_count,by = "species",all.x=T)

## calculate the percentage of total samples each compound appears in per species accession
#polar_compound_counts$perc <- round(polar_compound_counts$Freq/polar_compound_counts$count,2)

## label those above your predetermined cutoff as compounds to keep for that accession
#polar_compound_counts$keep <- polar_compound_counts$perc >= cutoff

## merge the count table with the fill compounds table and remove all that were not labeled to keep (not above the cutoff)
#polar_fill_comps_table_1 <- merge(polar_fill_comps_table,polar_compound_counts[,c("species","cluster.index","keep")], by = c("species","cluster.index"))
#polar_fill_comps_table_2 <- polar_fill_comps_table_1[polar_fill_comps_table_1$keep,-c("keep")]
#polar_fill_comp <- polar_fill_comps_table_2


## load in list of matches across both C18 and polar columns created in 'match_compounds_between_columns.py', remove anything that is not labeled as a true match
match_replacement <- read.csv("./data/matches_across_columns.csv",row.names=1)
match_replacement <- match_replacement[match_replacement$true_match!="N",]

###### Remove compounds from C18 / polar lists based on compounds matching ##########
c18.remove <- list()
polar.remove <- list()
for (i in 1:nrow(match_replacement)){
  # label each column based on structure of match replacement dataframe
  column1 <- "amide"
  column2 <- "C18"
  # determine the column of the compound chosen to keep of the true matches
  column_keep <- unlist(strsplit(as.character(match_replacement$comp2keep[i]),split="_"))[1]
  if (column1!=column_keep){
    # add compound to polar (amide) removal list if column_keep is not 'amide'
    column1.remove <- unlist(strsplit(as.character(match_replacement$compound_1[i]),split="_"))[2]
    polar.remove <- c(polar.remove,column1.remove)
  }
  if (column2!=column_keep){
    # add compound to C18 removal list if column_keep is not 'C18'
    column2.remove <- unlist(strsplit(as.character(match_replacement$compound_2[i]),split="_"))[2]
    c18.remove <- c(c18.remove,column2.remove)
  }
  
}
## reformat removal lists so they are a single list and not a list of lists
c18.remove <- unlist(c18.remove)
polar.remove <- unlist(polar.remove)

## remove compounds from the fill compounds tables
polar_fill_comp <- polar_fill_comp[!as.character(polar_fill_comp$compound_number) %in% polar.remove,]
c18_fill_comp <- c18_fill_comp[!as.character(c18_fill_comp$compound_number) %in% c18.remove,]

## combine the fill compounds tables so that it includes all samples from both polar and C18 columns
all_fill_comp <- rbind(c18_fill_comp,polar_fill_comp)

###### Alter format of combined filled compound table to create a species accession by compound matrix using function from 'create_pairwiseComps_sampsByComps.R' #######
all_accession_bycomps <- make_sampsByCompounds(all_fill_comp,by_species=T,method="combined")

####### Detect available cores and set up environment for parallelization #######
cores = detectCores()
cl <- makeCluster(cores[1] / 2) #not to overload your computer
registerDoParallel(cl)


####### Create list of phenolics and saponins based on mz (mass-to-charge) and RMD (residual mass defect) #######
# Load feature info associated with each compound from compound info loaded earlier
allcomps_feature_info <- compound_info[compound_info$column %in% c("C18","polar"),]
# saponins approximated as greater than 580 Da and RMD > 425 based on machine learning analysis
all_sap_compounds <- allcomps_feature_info[allcomps_feature_info$mz >=580 & allcomps_feature_info$parent_RMD >= 425,"cluster.index"]
# everything else with phenolics
all_phen_compounds <- allcomps_feature_info[!(allcomps_feature_info$compound_number %in% all_sap_compounds), "cluster.index"]


###### Create species by compound matrix for each compound class #######
# with raw TIC values
sap_accession_bycomps <- all_accession_bycomps[,names(all_accession_bycomps) %in% as.character(all_sap_compounds), drop=FALSE]
phen_accession_bycomps <- all_accession_bycomps[,names(all_accession_bycomps) %in% as.character(all_phen_compounds), drop=FALSE]

# standardize by species so that total compound abundance in each sample sums to 1.0
sap_accession_bycomps_stand <- standardizeByRow(sap_accession_bycomps)
phen_accession_bycomps_stand <- standardizeByRow(phen_accession_bycomps)


###### Run chemical similarity calculation from 'chem_similarity_function.R' for each of the saponin and phenolic matrices separately ######
accession_similarity_matrix_sap <- foreach(i = 1:nrow(sap_accession_bycomps_stand), .combine = rbind) %:% foreach(j = 1:nrow(sap_accession_bycomps_stand)) %dopar% {
  if(i <= j) {
    chemical_similarity_single(row.names(sap_accession_bycomps_stand)[i], row.names(sap_accession_bycomps_stand)[j], sap_accession_bycomps_stand, all.pairwise.comps)
  }
  else NA
}
# change format and fill in the rest of the pairwise sample similarity matrix
accession_similarity_matrix_sap <- as.data.frame(accession_similarity_matrix_sap)
# make row and column labels match the species accession names
names(accession_similarity_matrix_sap) <- row.names(sap_accession_bycomps_stand)
row.names(accession_similarity_matrix_sap) <- row.names(sap_accession_bycomps_stand)
# make the matrix square by filling in the other half with the corresponding values
for(i in 1:ncol(accession_similarity_matrix_sap)) accession_similarity_matrix_sap[,i] <- unlist(accession_similarity_matrix_sap[,i])
for(i in 1:nrow(accession_similarity_matrix_sap)) {
  for(j in i:nrow(accession_similarity_matrix_sap)) {
    accession_similarity_matrix_sap[j,i] <- accession_similarity_matrix_sap[i,j]
  }
}

# write an output CSV in case you have to reload this matrix and don't want to compute the calculation again
write.csv(accession_similarity_matrix_sap,"./results/accession_saponin_similarity_matrix.csv",row.names=T)

####### Repeat similarity calculation for phenolics #######
accession_similarity_matrix_phen <- foreach(i = 1:nrow(accession_similarity_matrix_phen), .combine = rbind) %:% foreach(j = 1:nrow(accession_similarity_matrix_phen)) %dopar% {
  if(i <= j) {
    chemical_similarity_single(row.names(accession_similarity_matrix_phen)[i], row.names(accession_similarity_matrix_phen)[j], accession_similarity_matrix_phen, all.pairwise.comps)
  }
  else NA
}
# change format and fill in the rest of the pairwise sample similarity matrix
accession_similarity_matrix_phen <- as.data.frame(accession_similarity_matrix_phen)
# make row and column labels match the species accession names
names(accession_similarity_matrix_phen) <- row.names(accession_similarity_matrix_phen)
row.names(accession_similarity_matrix_phen) <- row.names(accession_similarity_matrix_phen)
# make the matrix square by filling in the other half with the corresponding values
for(i in 1:ncol(accession_similarity_matrix_phen)) accession_similarity_matrix_phen[,i] <- unlist(accession_similarity_matrix_phen[,i])
for(i in 1:nrow(accession_similarity_matrix_phen)) {
  for(j in i:nrow(accession_similarity_matrix_phen)) {
    accession_similarity_matrix_phen[j,i] <- accession_similarity_matrix_phen[i,j]
  }
}

similarity_matrix_phen <- accession_similarity_matrix_phen

# write an output CSV in case you have to reload this matrix and don't want to compute the calculation again
write.csv(similarity_matrix_phen,"./results/accession_phenolic_similarity_matrix.csv",row.names=T)

# save each accession by compound and similarity matrix as a generic object
accession_comp <- all_accession_bycomps
accession_comp_sap <- sap_accession_bycomps_stand
accession_comp_phen <- phen_accession_bycomps_stand
similarity_matrix_sap <- accession_similarity_matrix_sap
similarity_matrix_phen <- accession_similarity_matrix_phen

###### Connect to new database
mydb <- dbConnect(MySQL(), user='u***', password='**********', dbname='XCMS_DB', host='chpc.utah.edu')

####### Load data on tyrosine content and extraction weights #######
tyr.pct <- dbGetQuery(mydb,"SELECT AVG(percent_tyrosine) AS percent_tyrosine, species_code FROM (SELECT * FROM Tyrosine WHERE keep_or_drop = 'keep') GROUP BY species_code")
extr.pct <- dbGetQuery(mydb,"SELECT AVG(a.percent_extracted) as percent_extracted, a.species_code FROM (SELECT extraction_weight.percent_extracted, chemistry.species_code FROM extraction_weight LEFT JOIN chemistry ON extraction_weight.`chem_number` = `chemistry`.`chem_number`) a GROUP BY a.species_code")


####### Combine phenolic and saponin similarity with tyrosine data and overall chemical investment (extraction weight) #######
# calculate total investment in phenolics and saponins based on sum of TIC for all compounds in each class 
comp.class.pcts <- data.frame("sample" = row.names(accession_comp), 
                              "phenSumTIC" = sapply(1:nrow(similarity_matrix_phen), function(x) sum(similarity_matrix_phen[x,])), 
                              "sapSumTIC" = sapply(1:nrow(accession_comp_sap), function(x) sum(accession_comp_sap[x,])), 
                              "species_code" = sapply(1:nrow(accession_comp),function(x) unlist(strsplit(row.names(accession_comp)[x], split = "_"))[1]),
                              stringsAsFactors=FALSE)
# take log of TIC for each class
comp.class.pcts$logPhen <- log(comp.class.pcts$phenSumTIC)
comp.class.pcts$logSap <- log(comp.class.pcts$sapSumTIC)
comp.class.pcts$logSap[is.infinite(comp.class.pcts$logSap)] <- 0

# calculate percent investment in saponins/phenolics
comp.class.pcts$phenTICpct <- comp.class.pcts$phenSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts$sapTICpct <- comp.class.pcts$sapSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)

# merge extraction weight and tyrosine % dry weight per accession onto the compound class investment dataframe
comp.class.pcts <- join(comp.class.pcts, extr.pct, by="species_code", type="left", match="all")
comp.class.pcts <- join(comp.class.pcts, tyr.pct, by="species_code", type="left", match="all")
# label any species accessions with no tyrosine data as 0 % tyrosine
comp.class.pcts$percent_tyrosine[is.na(comp.class.pcts$percent_tyrosine)] <- 0
# convert the percentages to proportions between 0 and 1
comp.class.pcts$percent_tyrosine <- comp.class.pcts$percent_tyrosine / 100
comp.class.pcts$percent_extracted <- comp.class.pcts$percent_extracted / 100
# calculate tyrosine as % of total secondary metabolite extractions
comp.class.pcts$tyr.final.pct <- comp.class.pcts$percent_tyrosine / (comp.class.pcts$percent_extracted + comp.class.pcts$percent_tyrosine)
# calculate phenolics and saponins as % of total secondary metabolite extractions
comp.class.pcts$phensap.final.pct <- comp.class.pcts$percent_extracted / (comp.class.pcts$percent_extracted + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phen.final.pct <- comp.class.pcts$phenTICpct * comp.class.pcts$phensap.final.pct
comp.class.pcts$sap.final.pct <- comp.class.pcts$sapTICpct * comp.class.pcts$phensap.final.pct

####### Create matrices of pairwise minimum phenolic, saponin, and tyrosine investment #######
pairwise.phen.min <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.sap.min <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.tyr.min <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.all.min <- outer(comp.class.pcts$phensap.final.pct, comp.class.pcts$phensap.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))

####### Combine matrices, using minimum investment in each compound class as weights
pairwise.spp <- similarity_matrix_phen*pairwise.phen.min + similarity_matrix_sap*pairwise.sap.min + pairwise.tyr.min

# make sure matrix is still between 0 and 1, with self comparisons being the highest values (1)
for (i in 1:nrow(pairwise.spp)){
  for (j in 1:nrow(pairwise.spp)){
    pairwise.spp[i,j] <- pairwise.spp[i,j] / max(pairwise.spp[i,i],pairwise.spp[j,j])
    
  }
}

# write final output CSV
write.csv(pairwise.spp, "./results/allcomps_accession_similarity_matrix.csv", row.names = TRUE)
