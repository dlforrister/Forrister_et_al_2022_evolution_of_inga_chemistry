#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
cat("  START TIME:", date(),"\n")
if (length(args)!=1){
  stop("Wrong arguments -> need 1!")
} else {
  SPECIES_ID = as.integer(args[1])
  cat("  SPECIES_ID:", SPECIES_ID,"\n",sep="")
}
itr <- SPECIES_ID

library(foreach)
library(doParallel)
library(plyr)
library(data.table)
library(RMySQL)
source("/uufs/chpc.utah.edu/common/home/inga-group1/4_directories_chem_evolution_2019_03_26/code/6_chemsim/create_pairwiseComps_sampsByComps_scankey_normTIC_compid.R")
source("/uufs/chpc.utah.edu/common/home/inga-group1/4_directories_chem_evolution_2019_03_26/code/6_chemsim/chem_similarity_function.R")
setwd("/uufs/chpc.utah.edu/common/home/inga-group1/4_directories_chem_evolution_2019_03_26/")

############### combined ###################
####### Download 'view raw spectra' from GNPS results page and use make_pairwisecomps function from 'create_pairwiseComps_sampsByComps.R' to creating pairwise compound similarity matrix #######
all.pairwise.comps <- read.csv("./data/0_combined/allcomps_pairwise_thru_network.csv",row.names=1)
names(all.pairwise.comps) <- gsub("X","",names(all.pairwise.comps))

## read in all null matricies
files <- list.files("/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/null_model",pattern = ".csv",full.names = T)

all_data <- data.frame()
for(file in files){
  file_ind <- read.csv(file,stringsAsFactors = F)
  file_ind$iteration <- gsub(".csv","",gsub("/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/null_model/Evolved_Random_no_cutoff_","",file))
  all_data <- rbind(all_data,file_ind)
}
#all_data[1:10,1:10]

names(all_data)[-c(1,which(names(all_data) == "iteration"))] <- gsub("X","",names(all_data)[-c(1,which(names(all_data) == "iteration"))])

all_sampsbycomps_species_actual <- read.csv("./data/0_combined/final_sampcomps/6_sampsbycomps_combined_nocutoff_species_2020_10_7.csv",row.names=1)
names(all_sampsbycomps_species_actual) <- gsub("X","",names(all_sampsbycomps_species_actual))

all_data$node_label[all_data$node_label == "M27"] <- "M27b"
all_data$node_label[all_data$node_label == "Zygia_mediana"] <- "Zygia mediana"

all_data_species <- all_data[all_data$node_label %in% row.names(all_sampsbycomps_species_actual),]


all_sampsbycomps_species_null <- all_data_species[all_data_species$iteration == itr,-which(names(all_data_species) == "iteration")]
row.names(all_sampsbycomps_species_null) <- all_sampsbycomps_species_null$node_label
all_sampsbycomps_species_null <- all_sampsbycomps_species_null[order(all_sampsbycomps_species_null$node_label),]
all_sampsbycomps_species_null <- all_sampsbycomps_species_null[,-1]

all_sampsbycomps_species <- all_sampsbycomps_species_null

####### Detect available cores and set up environment for parallelization #######
cores = detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

####### Create list of phenolics and saponins based on mz and rmd #######
# Load table listing all features associated with each compound created in 'run_compound_table_code.py'
feature_info <- read.table(file = './data/00_C18/all_comps_info_V3_2020_5_4_knowns_anchors.csv', sep = ',', header = TRUE)
allcomps_feature_info <- feature_info[feature_info$column %in% c("C18","polar"),]

saponin_comps <- read.csv("/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/data_general/2_Cluster_index_by_saponans_V4_11_17.csv")

# saponins approximated as greater than 580 Da and RMD > 425
all_sap_compounds <- allcomps_feature_info$cluster.index[allcomps_feature_info$cluster.index %in% saponin_comps$cluster.index[saponin_comps$Is_Saponin == 1]]

# everything else with phenolics
all_phen_compounds <- allcomps_feature_info$cluster.index[allcomps_feature_info$cluster.index %in% saponin_comps$cluster.index[saponin_comps$Is_Saponin == 0]]


###### Create species by compound matrix for each compound class #######
# with raw TIC values
all_sampsByCompoundsSap <- all_sampsbycomps_species[,names(all_sampsbycomps_species) %in% as.character(all_sap_compounds), drop=FALSE]
all_sampsByCompoundsPhen <- all_sampsbycomps_species[,names(all_sampsbycomps_species) %in% as.character(all_phen_compounds), drop=FALSE]

# standardize by species so that total compound abundance in each sample sums to 1.0
all_sampsCompsStandSap <- standardizeByRow(all_sampsByCompoundsSap)
all_sampsCompsStandPhen <- standardizeByRow(all_sampsByCompoundsPhen)

## species ##
print("Started species saponin similarity calculation")

all_species_similarity_matrix_sap <- foreach(i = 1:nrow(all_sampsCompsStandSap), .combine = rbind) %:% foreach(j = 1:nrow(all_sampsCompsStandSap)) %dopar% {
  if(i <= j) {
    chemical_similarity_single(row.names(all_sampsCompsStandSap)[i], row.names(all_sampsCompsStandSap)[j], all_sampsCompsStandSap, all.pairwise.comps)
  }
  else NA
}
# change format and fill in rest of pairwise sample similarity matrix
all_species_similarity_matrix_sap <- as.data.frame(all_species_similarity_matrix_sap)
names(all_species_similarity_matrix_sap) <- row.names(all_sampsCompsStandSap)
row.names(all_species_similarity_matrix_sap) <- row.names(all_sampsCompsStandSap)
for(i in 1:ncol(all_species_similarity_matrix_sap)) all_species_similarity_matrix_sap[,i] <- unlist(all_species_similarity_matrix_sap[,i])
for(i in 1:nrow(all_species_similarity_matrix_sap)) {
  for(j in i:nrow(all_species_similarity_matrix_sap)) {
    all_species_similarity_matrix_sap[j,i] <- all_species_similarity_matrix_sap[i,j]
  }
}

similarity_matrix_sap <- all_species_similarity_matrix_sap

for (i in 1:nrow(similarity_matrix_sap)){
  for (j in 1:nrow(similarity_matrix_sap)){
    if (i==j){
    similarity_matrix_sap[i,j] <- 1
    }
  }
}

write.csv(similarity_matrix_sap,paste("/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/null_model_chem_simulations/",itr,"_allcomps_sample_saponin_similarity.csv",sep=""),row.names=T)

####### Repeat similarity calculation for phenolics #######
print("Started species phenolic similarity calculation")
all_species_similarity_matrix_phen <- foreach(i = 1:nrow(all_sampsCompsStandPhen), .combine = rbind) %:% foreach(j = 1:nrow(all_sampsCompsStandPhen)) %dopar% {
  if(i <= j) {
    chemical_similarity_single(row.names(all_sampsCompsStandPhen)[i], row.names(all_sampsCompsStandPhen)[j], all_sampsCompsStandPhen, all.pairwise.comps)
  }
  else NA
}
all_species_similarity_matrix_phen <- as.data.frame(all_species_similarity_matrix_phen)
names(all_species_similarity_matrix_phen) <- row.names(all_sampsCompsStandPhen)
row.names(all_species_similarity_matrix_phen) <- row.names(all_sampsCompsStandPhen)
for(i in 1:ncol(all_species_similarity_matrix_phen)) all_species_similarity_matrix_phen[,i] <- unlist(all_species_similarity_matrix_phen[,i])
for(i in 1:nrow(all_species_similarity_matrix_phen)) {
  for(j in i:nrow(all_species_similarity_matrix_phen)) {
    all_species_similarity_matrix_phen[j,i] <- all_species_similarity_matrix_phen[i,j]
  }
}

similarity_matrix_phen <- all_species_similarity_matrix_phen
write.csv(similarity_matrix_phen,paste("/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/null_model_chem_simulations/",itr,"_allcomps_sample_phenolic_similarity.csv",sep=""),row.names=T)

sampsByCompounds <- all_sampsbycomps_species
sampsByCompoundsSap <- all_sampsByCompoundsSap
sampsByCompoundsPhen <- all_sampsByCompoundsPhen

###### Connect to new database
mydb = dbConnect(MySQL(), user='u6019685', password='xnukasq9gvenvk8NDglt', dbname='newinga', host='mysql.chpc.utah.edu')
myoldb = dbConnect(MySQL(), user='u6019685', password='xnukasq9gvenvk8NDglt', dbname='inga_2015_06_01', host='mysql.chpc.utah.edu')

####### Load data on tyrosine content and extraction weights #######
tyr.pct <- dbGetQuery(mydb,"SELECT AVG(a.percent_tyrosine) as percent_tyrosine, a.species_code FROM (
                      SELECT * FROM Tyrosine WHERE keep_or_drop = 'keep') a
                      group by a.species_code")
extr.pct <- dbGetQuery(myoldb,"SELECT AVG(a.percent_extracted) as percent_extracted, a.species_code FROM (SELECT extraction_weight.percent_extracted, chemistry.species_code FROM extraction_weight 
                       LEFT JOIN chemistry ON extraction_weight.`Chem#` = `chemistry`.`Chem#`) a
                       group by a.species_code")

# change recent species codes because extraction mostly from old DB
extr.pct$species_code[extr.pct$species_code=="T78"] <- "T67b"
extr.pct$species_code[extr.pct$species_code=="M27a"] <- "M57"
extr.pct$species_code[extr.pct$species_code=="M30"] <- "M8"

# for N30 because it was only run recently
extr.pct2 <- dbGetQuery(mydb,"SELECT AVG(real_percent_extracted) as percent_extracted, species_code FROM `Extraction_Percent` WHERE species_code = 'N30' and method in (19,20)")
extr.pct <- rbind(extr.pct,extr.pct2)

###
phy_code <- read.csv("../Figures_Evol_Inga_Chemistry/data_general/species_code_name_species_2020_8_25.csv")

extr.pct <- merge(extr.pct,phy_code,by.x="species_code",by.y="Species_Code")[,c("percent_extracted","Species")]

names(extr.pct) <- c("percent_extracted","species_code")

tyr.pct <- merge(tyr.pct,phy_code,by.x="species_code",by.y="Species_Code")[,c("percent_tyrosine","Species")]

names(tyr.pct) <- c("percent_tyrosine","species_code")

extr.pct <- setDT(extr.pct)[,.(percent_extracted = mean(percent_extracted)),by="species_code"][,c(2,1)]

tyr.pct <- setDT(tyr.pct)[,.(percent_tyrosine = mean(percent_tyrosine)),by="species_code"][,c(2,1)]

extr.pct$species_code <- as.character(extr.pct$species_code)
extr.pct$species_code[extr.pct$species_code == "Zygia_mediana"] <- "Zygia mediana"
extr.pct$species_code[extr.pct$species_code == "M27"] <- "M27b"

tyr.pct$species_code <- as.character(tyr.pct$species_code)
tyr.pct$species_code[tyr.pct$species_code == "M27"] <- "M27b"


####### Combine phenolic and saponin similarity with tyrosine data #######
# calculate total investment in phenolics and saponins based on sum of TIC for all compounds in each class 
comp.class.pcts <- data.frame("sample" = row.names(sampsByCompounds), 
                              "phenSumTIC" = sapply(1:nrow(sampsByCompoundsPhen), function(x) sum(sampsByCompoundsPhen[x,])), 
                              "sapSumTIC" = sapply(1:nrow(sampsByCompoundsSap), function(x) sum(sampsByCompoundsSap[x,])), 
                              "species_code" = sapply(1:nrow(sampsByCompoundsPhen), 
                                                      function(x) unlist(strsplit(row.names(sampsByCompoundsPhen)[x], split = "_"))[1]),
                              stringsAsFactors=FALSE)
# take log of TIC for each class
comp.class.pcts$logPhen <- log(comp.class.pcts$phenSumTIC)
comp.class.pcts$logSap <- log(comp.class.pcts$sapSumTIC)
comp.class.pcts$logSap[is.infinite(comp.class.pcts$logSap)] <- 0

# calculate percent investment in saponins/phenolics
comp.class.pcts$phenTICpct <- comp.class.pcts$phenSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)
comp.class.pcts$sapTICpct <- comp.class.pcts$sapSumTIC / (comp.class.pcts$phenSumTIC + comp.class.pcts$sapSumTIC)

# get intermediate polarity (phenolics+saponins) extraction weight and tyrosine % dry weight
comp.class.pcts <- join(comp.class.pcts, extr.pct, by="species_code", type="left", match="all")
comp.class.pcts <- join(comp.class.pcts, tyr.pct, by="species_code", type="left", match="all")
comp.class.pcts$percent_tyrosine[is.na(comp.class.pcts$percent_tyrosine)] <- 0
comp.class.pcts$percent_tyrosine <- comp.class.pcts$percent_tyrosine / 100
comp.class.pcts$percent_extracted <- comp.class.pcts$percent_extracted / 100
# calculate tyrosine as % of total secondary metabolite extractions
comp.class.pcts$tyr.final.pct <- comp.class.pcts$percent_tyrosine / (comp.class.pcts$percent_extracted + comp.class.pcts$percent_tyrosine)
# calculate phenolics+saponins as % of total secondary metabolite extractions
comp.class.pcts$phensap.final.pct <- comp.class.pcts$percent_extracted / (comp.class.pcts$percent_extracted + comp.class.pcts$percent_tyrosine)
comp.class.pcts$phen.final.pct <- comp.class.pcts$phenTICpct * comp.class.pcts$phensap.final.pct
comp.class.pcts$sap.final.pct <- comp.class.pcts$sapTICpct * comp.class.pcts$phensap.final.pct

####### Create matrices of pairwise minimum phenolic, saponin, and tyrosine investment #######
pairwise.phen.min <- outer(comp.class.pcts$phen.final.pct, comp.class.pcts$phen.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.sap.min <- outer(comp.class.pcts$sap.final.pct, comp.class.pcts$sap.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.tyr.min <- outer(comp.class.pcts$tyr.final.pct, comp.class.pcts$tyr.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))
pairwise.all.min <- outer(comp.class.pcts$phensap.final.pct, comp.class.pcts$phensap.final.pct, FUN = function(X,Y) sapply(1:length(X), function(i) min(X[i],Y[i])))

####### If using separate calculations for phenolics and saponins: Combine matrices, using minimum investment in each compound class as weights
pairwise.spp <- similarity_matrix_phen*pairwise.phen.min + similarity_matrix_sap*pairwise.sap.min + pairwise.tyr.min

for (i in 1:nrow(pairwise.spp)){
  for (j in 1:nrow(pairwise.spp)){
    pairwise.spp[i,j] <- pairwise.spp[i,j] / max(pairwise.spp[i,i],pairwise.spp[j,j])
    
  }
}

write.csv(pairwise.spp, paste("/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/null_model_chem_simulations/",itr,"_allcomps_sample_similarity.csv",sep=""), row.names = TRUE)
