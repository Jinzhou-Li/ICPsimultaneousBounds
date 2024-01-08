get_simul_FDbound <- function(R, pvalues_allsets, allsets){
  
  if(length(R)==0){return(0)} # empty set case
  
  subsets_R <- rev(get_allsubsets(R))
  
  FDbound <- length(R)
  for (index_subset in 1:length(subsets_R)) {
    I_temp <- subsets_R[[index_subset]]
    if(length(I_temp)==0){FDbound <- 0; break} # empty set case
    
    p_temp <- get_pvalue_set(I_temp, pvalues_allsets, allsets)
    if(p_temp > alpha) {FDbound <- length(I_temp); break}
  }
  
  return(FDbound)
}

####### sub-functions
set_index_func <- function(i_subset, subsets, allsets){
  set_compare <- sapply(allsets, function(set, targetset) {identical(set,targetset)}, subsets[[i_subset]])
  return(which(set_compare))
}

get_pvalue_set <- function(R, pvalues_allsets, allsets){
  p <- max(lengths(allsets))
  
  if(length(R)==0) {return(max(pvalues_allsets))}
  if(length(R)==p) {return(pvalues_allsets[1])}  # the p-value corresponds to empty set
  
  usevar <- (1:p)[-R]
  subsets <- get_allsubsets(usevar)

  set_index_vec <- sapply(1:length(subsets), set_index_func, subsets, allsets)
  p_R <- max(pvalues_allsets[set_index_vec])
  
  return(p_R)
}

#####
get_allsubsets <- function(usevar){
  subsets <- list()
  
  if(length(usevar)>0){
    for (ic in ((1:2^length(usevar))-1)){
      subsets[[ic+1]] <- usevar[which( ((ic %/% 2^(0:(length(usevar)-1))) %% 2 )==1)]
    }
  }
  subsets <- unique(subsets)
  le <- sapply(subsets,length)
  subsets <- subsets[order(le)]
  
  return(subsets)
}