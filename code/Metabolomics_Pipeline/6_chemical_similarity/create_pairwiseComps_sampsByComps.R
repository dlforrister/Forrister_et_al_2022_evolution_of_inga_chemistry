library(reshape2)
library(igraph)
library(splitstackshape)

# function takes network from GNPS and makes compound x compound similarity matrix using cosine scores and calculating through-network similarity (default), can prevent this by using through_network = FALSE
# requires path to folder that come from 'View Raw Spectra' in GNPS
make_pairwisecomps <- function(network_filepath, method, through_network = TRUE) {
# load .pairsinfo file from networkedges_selfloop folder and raw spectra file from main folder
  if(is.character(network_filepath)) {
  network_edges_path <- list.files(paste(network_filepath, "/networkedges_selfloop", sep = ""), full.names = TRUE)
  #raw_spectra_path <- list.files(network_filepath, full.names = TRUE)[grep("main", list.files(network_filepath))]
  network_clean <- read.delim(file = network_edges_path, header=TRUE, sep="\t")
  #compound_cluster_match <- read.delim(file = raw_spectra_path,header=TRUE, sep="\t")[,c("cluster.index","ScanNumber")]
  compound_cluster_match <- read.table(file = '/uufs/chpc.utah.edu/common/home/inga-group1/4_directories_chem_evolution_2019_03_26/data/00_C18/all_comps_info_V3_2020_5_4_knowns_anchors.csv', sep = ',', header = TRUE)[,c("cluster.index","compound_number","column")]
  compound_cluster_match$compound_number <- as.numeric(as.character(compound_cluster_match$compound_number))
  compound_cluster_match$cluster.index <- as.numeric(as.character(compound_cluster_match$cluster.index))
  names(network_clean) <- c("CLUSTERID1","CLUSTERID2","DeltaMZ","V4","Cosine","OtherScore","V7")
  network_clean <- network_clean[,!is.na(names(network_clean))]

  # convert network from clusterid to scan number(== compound number)
  names(compound_cluster_match) <- c("CLUSTERID1","NODE1","column")
  compound_network <- merge(network_clean, compound_cluster_match[,c("CLUSTERID1","NODE1")], by="CLUSTERID1",all.x=TRUE)
  names(compound_cluster_match) <- c("CLUSTERID2","NODE2","column")
  compound_network <- merge(compound_network, compound_cluster_match[,c("CLUSTERID2","NODE2")], by="CLUSTERID2",all.x=TRUE)
  # remove anchors and knowns
  compound_network <- compound_network[!is.na(compound_network$NODE1) & !is.na(compound_network$NODE2),]

  compound_network$Cosine <- as.numeric(as.character(compound_network$Cosine))
  # subset for c18 column
  if (method =="c18"){
    c18.remove <- as.numeric(unique(compound_cluster_match[compound_cluster_match$column !="C18","CLUSTERID2"]))
    compound_network <- compound_network[!(compound_network$CLUSTERID2 %in% c18.remove),]
    compound_network <- compound_network[,c("NODE1","NODE2","DeltaMZ","Cosine")]
  }
  if (method=="polar"){
    polar.remove <- as.numeric(unique(compound_cluster_match[compound_cluster_match$column !="polar","CLUSTERID2"]))
    compound_network <- compound_network[!(compound_network$CLUSTERID2 %in% polar.remove),]
    compound_network <- compound_network[,c("NODE1","NODE2","DeltaMZ","Cosine")]
  }
  if (method=="combined"){
    remove <- as.numeric(unique(compound_cluster_match[!as.character(compound_cluster_match$column) %in% c("polar","C18"),"CLUSTERID2"]))
    compound_network <- compound_network[!(compound_network$CLUSTERID2 %in% remove),]
    compound_network <- compound_network[,c("CLUSTERID1","CLUSTERID2","DeltaMZ","Cosine")]
    names(compound_network) <- c("NODE1","NODE2","DeltaMZ","Cosine")
  }
  
  
  } else {
   compound_network <- network_filepath 
  }
  
  # make changes to structure of network, calculate inverse cosine score for through-network calculations

  nodes <- unique(c(compound_network$NODE1, compound_network$NODE2))
  
  miss_n1 <- sort(nodes[which(!nodes %in% compound_network$NODE1)])
  miss_n2 <- sort(nodes[which(!nodes %in% compound_network$NODE2)])
  
  miss_n1_row <- data.frame("NODE1"=miss_n1,"NODE2"=miss_n1,"DeltaMZ"=rep(0,length(miss_n1)),"Cosine"=rep(1.0,length(miss_n1)))
  miss_n2_row <- data.frame("NODE1"=miss_n2,"NODE2"=miss_n2,"DeltaMZ"=rep(0,length(miss_n2)),"Cosine"=rep(1.0,length(miss_n2)))
  
  compound_network <- rbind(compound_network,miss_n1_row)
  compound_network <- rbind(compound_network,miss_n2_row)
  compound_network$inv_cos_score <- 1/compound_network$Cosine
  
  # compute shortest paths using inverse of cosine, then taking the inverse again
  g2 <- graph_from_data_frame(d=compound_network, vertices=data.frame(nodes), directed=FALSE)
  if(through_network == TRUE) {
  shortest_paths_list_inverse <- lapply(1:length(nodes), function(x) shortest.paths(g2,v=as.character(nodes[x]),weights=E(g2)$inv_cos_score))
  shortest_paths_df_inverse <- do.call(rbind.data.frame, shortest_paths_list_inverse)
  short_paths_inv_inv <- 1/shortest_paths_df_inverse
  short_paths_inv_inv[short_paths_inv_inv == Inf] <- 1
  pairwise.comps <- short_paths_inv_inv
  }
  else{
    pairwise.comps <- as.data.frame(matrix(0, ncol = length(nodes), nrow = length(nodes)))
    names(pairwise.comps) <- as.character(sort(nodes))
    row.names(pairwise.comps) <- as.character(sort(nodes))
    for(i in 1:nrow(compound_network)) {
      pairwise.comps[row.names(pairwise.comps) == as.character(compound_network$NODE1[i]),
                     names(pairwise.comps) == as.character(compound_network$NODE2[i])] <- compound_network$Cosine[i]
      pairwise.comps[row.names(pairwise.comps) == as.character(compound_network$NODE2[i]),
                     names(pairwise.comps) == as.character(compound_network$NODE1[i])] <- compound_network$Cosine[i]
    }
    for(i in 1:nrow(pairwise.comps)) {
      pairwise.comps[i,i] <- 1
    }
    }
  
  return(pairwise.comps) }

# function to create sampsByCompounds table
# requires compound tic table, which was created in python
# argument by_species = TRUE specifies that samples should be grouped by species. If FALSE (default), each sample will remain separate.
# if you want to remove some samples or species, specify them in samps_to_remove

make_sampsByCompounds <- function(compound_tic_table_path, samps_to_remove = character(0), by_species = FALSE,method=character(0)) {
  if(is.character(compound_tic_table_path)){ compound_tic_table <- read.csv(compound_tic_table_path)}
  else {compound_tic_table <- compound_tic_table_path}
  if (method=="combined"){
    sampsByCompounds <- dcast(compound_tic_table, cluster.index ~ sample, value.var = "norm_TIC")
     
    if(by_species == TRUE) {
      compound_tic_table_species <- data.frame(compound_tic_table,cSplit(compound_tic_table, 'sample', sep="_", type.convert=FALSE))
      compound_tic_table_species_1 <- aggregate(compound_tic_table_species$TIC, by=list(Category=compound_tic_table_species$cluster.index,compound_tic_table_species$sample_1),FUN=mean)
      names(compound_tic_table_species_1) <- c("cluster.index","species","TIC")
      sampsBySpecies <- dcast(compound_tic_table_species_1, cluster.index  ~ species, value.var = "TIC", fun.aggregate = sum)
      sampsByCompounds <-sampsBySpecies }
    
    if(by_species== FALSE){
      match_samples <- read.csv("/uufs/chpc.utah.edu/common/home/inga-group1/4_directories_chem_evolution_2019_03_26/data/0_combined/combine_polar_and_C18_sample_level_V1_2020_8_2.csv",row.names=1)
      match_samples$key <- paste(match_samples$spec_code,match_samples$sample_index,sep="_")  
      compound_tic_table_combined <- merge(compound_tic_table,match_samples[,c("sample","key")], by = "sample", all.x =T)
      compound_tic_table_combined <- compound_tic_table_combined[!is.na(compound_tic_table_combined$key)]
      compound_tic_table_combined_1 <- aggregate(compound_tic_table_combined$TIC, by=list(Category=compound_tic_table_combined$cluster.index,compound_tic_table_combined$key),FUN=mean)
      names(compound_tic_table_combined_1) <- c("cluster.index","sample","TIC")
      sampsByCompounds <- dcast(compound_tic_table_combined_1, cluster.index  ~ sample, value.var = "TIC", fun.aggregate = sum)}
    
    sampsByCompounds <- as.data.frame(sampsByCompounds)
    sampsByCompounds <- t(sampsByCompounds)
    sampsByCompounds[is.na(sampsByCompounds)] <- 0
    sampsByCompounds <- as.data.frame(sampsByCompounds)
    row.names(sampsByCompounds)
    names(sampsByCompounds) <- sampsByCompounds["cluster.index",]
    sampsByCompounds <- sampsByCompounds[!row.names(sampsByCompounds)%in%c("cluster.index",samps_to_remove),]
    
    return(sampsByCompounds) }   
  else {
    sampsByCompounds <- dcast(compound_tic_table, compound_number ~ sample, value.var = "norm_TIC")
    
    if(by_species == TRUE) {
      compound_tic_table_species <- data.frame(compound_tic_table,cSplit(compound_tic_table, 'sample', sep="_", type.convert=FALSE))
      compound_tic_table_species_1 <- aggregate(compound_tic_table_species$TIC, by=list(Category=compound_tic_table_species$compound_number,compound_tic_table_species$sample_1),FUN=mean)
      names(compound_tic_table_species_1) <- c( "compound_number","species","TIC")
      sampsBySpecies <- dcast(compound_tic_table_species_1, compound_number  ~ species, value.var = "TIC", fun.aggregate = sum)
      sampsByCompounds <-sampsBySpecies }
    
    sampsByCompounds <- as.data.frame(sampsByCompounds)
    sampsByCompounds <- t(sampsByCompounds)
    sampsByCompounds[is.na(sampsByCompounds)] <- 0
    sampsByCompounds <- as.data.frame(sampsByCompounds)
    row.names(sampsByCompounds)
    names(sampsByCompounds) <- sampsByCompounds["compound_number",]
    sampsByCompounds <- sampsByCompounds[!row.names(sampsByCompounds)%in%c("compound_number",samps_to_remove),]
    
    return(sampsByCompounds)
  }
   }

standardizeByRow <- function(dataframe) {
  for(i in 1:nrow(dataframe)) {
    dataframe[i, ] <- dataframe[i, ] / sum(dataframe[i, ])
  }
  return(dataframe)
}
