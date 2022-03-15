library(gdata)

# calculates chemical similarity between two samples. Requires sample names, sample by compound matrix standardized by sample (sampsCompsStand0) and a pairwise compound matrix with compound similarity scores (pairwise.comps).
chemical_similarity_single <- function(sample_1, sample_2, sampsCompsStand, pairwise.comps) {
  spp1 = sample_1
  i_compounds = which(sampsCompsStand[spp1,] != 0) # indices of compounds present in sample1
    spp2 = sample_2
    j_compounds = which(sampsCompsStand[spp2,] != 0) # indices of compounds present in sample2
    if(length(i_compounds) == 0 | length(j_compounds) == 0) {
      return(0) # if either sample contains no compounds, return score of 0
    }
    else {
      ij_compounds = sort(unique(c(i_compounds, j_compounds))) # indices of compounds present in sample1 or sample2
      i_tics = sampsCompsStand[spp1, ij_compounds] # TICs of ij compounds from sample1
      j_tics = sampsCompsStand[spp2, ij_compounds] # TICs of ij compounds from sample2
      shared_tics = sapply(1:length(ij_compounds), function(x) min(i_tics[x], j_tics[x])) # minimum TIC value (shared TIC) for each compound
      i_tics = i_tics - shared_tics # remaining TIC in sample1 after shared TIC is subtracted out
      j_tics = j_tics - shared_tics # remaining TIC in sample2 after shared TIC is subtracted out
      similarity = sum(shared_tics) # sum of TIC invested in the same compounds (fist component of similarity score)
      current_compounds = data.matrix(pairwise.comps[names(i_tics),names(j_tics)]) # compound similarity matrix for compounds of interest
      while(sum(i_tics) > 0.0001 & sum(j_tics) > 0.0001) { # as long as there is still TIC unaccounted for....
        i_tics = i_tics[, i_tics > 0.0001, drop = FALSE] # TICs for compounds with significant TIC unaccounted for in sample1
        j_tics = j_tics[, j_tics > 0.0001, drop = FALSE] # TICs for compounds with significant TIC unaccounted for in sample2
        if(length(i_tics) < 1 | length(j_tics) < 1) break
        current_compounds = data.matrix(pairwise.comps[names(i_tics),names(j_tics)]) # compound similarity matrix for compounds of interest
        current_cos = max(current_compounds) # highest available cosine score from compound similarity matrix of remaining compounds
        if(current_cos == min(current_compounds)) { # if current cosine score is also the lowest remaining cosine score
          similarity = similarity + min(sum(i_tics), sum(j_tics)) * current_cos # multiply remaining TIC by this cosine and add to similarity
          break
        }
        current_comp_locations = which(current_compounds == current_cos) # get indices for compound pairings that have current cosine
        i.numbers = current_comp_locations %% nrow(current_compounds) # index for these compounds in sample1
        i.numbers[i.numbers==0] = nrow(current_compounds)
        j.numbers = ceiling(current_comp_locations/nrow(current_compounds)) # index for these compounds in sample2
        
        i.repeats = sapply(1:length(i.numbers), function(x) sum(i.numbers == i.numbers[x]))
        j.repeats = sapply(1:length(j.numbers), function(x) sum(j.numbers == j.numbers[x]))
        ij.min = sapply(1:length(i.numbers), function(x) min(i_tics[i.numbers[x]]/i.repeats[x], j_tics[j.numbers[x]]/j.repeats[x]))
        
        similarity = similarity + sum(ij.min * current_cos)
        
        for(k in 1:length(i.numbers)) {
          i_tics[i.numbers[k]] = i_tics[i.numbers[k]] - ij.min[k]
          j_tics[j.numbers[k]] = j_tics[j.numbers[k]] - ij.min[k]
        }
        current_compounds[current_comp_locations] <- 0
      }
  return(similarity)
    }}


