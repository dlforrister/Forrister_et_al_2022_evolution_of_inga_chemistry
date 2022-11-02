pairwise.comps.all <- read.csv(here("./data/allcomps_pairwise_thru_network_2020_10_7.csv"),row.names=1)

names(pairwise.comps.all) <- gsub("X","",names(pairwise.comps.all))
pairwise.comps.all[1:10,1:10]
pairwise.comps.all <- pairwise.comps.all[order(as.numeric(row.names(pairwise.comps.all))),order(as.numeric(row.names(pairwise.comps.all)))]

dim(pairwise.comps.all)

#Read in Raw sampsByCompounds amtrix and fix names
sampsByCompounds <- read.csv(here("./data/6_sampsbycomps_combined_nocutoff_species_2020_10_7.csv"),row.names=1)
names(sampsByCompounds) <- gsub("X","",names(sampsByCompounds))
sampsByCompounds[1:10,1:10]
sampsByCompounds <- round(sampsByCompounds,0)

dim(sampsByCompounds)

#use the phylogeny code to convert species codes into sample names matching the phylogeny
phy_code <- read.csv(here("./data/species_code_name_species_2020_8_25.csv"))

#samples in samps by cmps missing from phylogeny
row.names(sampsByCompounds)[!row.names(sampsByCompounds) %in% phy_code$Species[phy_code$Overal_Keep == "Keep"]]
phy_code$Species[phy_code$Overal_Keep == "Keep"][!phy_code$Species[phy_code$Overal_Keep == "Keep"] %in% row.names(sampsByCompounds)]

row.names(sampsByCompounds)[row.names(sampsByCompounds) == "M27b"] <- "M27"
row.names(sampsByCompounds)[row.names(sampsByCompounds) == "Zygia mediana"] <- "Zygia_mediana"

sampsByCompounds <- sampsByCompounds[row.names(sampsByCompounds) %in% phy_code$Species[phy_code$Overal_Keep == "Keep"],]

row.names(sampsByCompounds)[!row.names(sampsByCompounds) %in% phy_code$Species[phy_code$Overal_Keep == "Keep"]]
phy_code$Species[phy_code$Overal_Keep == "Keep"][!phy_code$Species[phy_code$Overal_Keep == "Keep"] %in% row.names(sampsByCompounds)]

#Ok, samps by comps names match the phylogeny completely.

#remove compounds that are not found in any sample since we have removed some samples
sampsByCompounds <- sampsByCompounds[,colSums(sampsByCompounds>0) > 0]

#ok, now have pairwise comps match the compounds found in samps x comps

dim(pairwise.comps.all)
dim(sampsByCompounds)

#now remove the samples in pairwiss comps that do not occur in sampsByCompounds
pairwise.comps.all <- pairwise.comps.all[row.names(pairwise.comps.all) %in% names(sampsByCompounds),names(pairwise.comps.all) %in% names(sampsByCompounds)]

dim(pairwise.comps.all)
dim(sampsByCompounds)

#Convert this to a distance so that phytochemical diversity can be calculated using Hill function
pairwise.comps.all_dist <- 1-pairwise.comps.all

sampsByCompounds[1:10,1:10]

