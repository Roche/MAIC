#' Preprocess data
#'
#' Preprocess data before estimating weights for matching-adjusted indirect comparison (MAIC).
#'
#' @param intervention_input A data frame containing individual patient data from
#'   the intervention study.
#' @param target_pop A data frame containing aggregate dataset for the target population. 
#'   Variables are followed by one of the following suffixes to denote the type of summary: 
#'   varname_mean, varname_sd, varname_median, varname_prop. 
#'   After preprocessing these summary suffixes, intervention_input is 
#'   centered using the aggregate data averages. 
#' @return A list containing 2 objects. First, a data frame named analysis_data
#'   containing intervention_data with additional columns named wt (weights) and
#'   wt_rs (rescaled weights). Second, a vector called matching_vars of the
#'   names of the centered matching variables used.
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @export

preprocess_data <- function(intervention_input, target_pop){
  
  # Check intervention_data and target_pop are data frame and match_cov is a character vector
  if(!is.data.frame(intervention_input)){stop("intervention_input is expected to be a data frame")}
  if(!is.data.frame(target_pop)){stop("target_pop is expected to be a data frame")}
  
  # Check if target_pop is 1 row of aggregate data
  if(nrow(target_pop)!=1){stop("target_pop should have exactly 1 row")}
  
  # Strip off naming convention in the aggregate data
  varnames <- gsub("_([^_]+)$","", names(target_pop))
  vartype <- gsub("^.*_","", names(target_pop))
  
  #Preprocess standard deviation
  for(i in 1:dim(target_pop)[2]){
    if(vartype[i] == "sd"){
      
      # retrieve which variable sd was specified
      varwithsd <- varnames[i]
      
      if(!paste0(varwithsd, "_mean") %in% names(target_pop)){
        stop(paste0("Also need to provide mean for ", varwithsd, " when specifying sd"))
      }
      
      # derive squared mean term
      target_pop[,paste0(varwithsd, "_squared_mean")] <- target_pop[,paste0(varwithsd, "_mean")]^2 + target_pop[,paste0(varwithsd, "_sd")]^2
      
      # remove standard deviation from the data frame
      target_pop <- target_pop[,-which(colnames(target_pop) == paste0(varwithsd, "_sd"))]
    }
  }
  
  # Preprocess median
  for(i in 1:dim(target_pop)[2]){
    if(vartype[i] == "median"){
      
      # retrieve which variable median was specified
      varwithmedian <- varnames[i]
      
      # make median into binary category
      intervention_input[,varwithmedian] <- ifelse(intervention_input[,varwithmedian] > target_pop[,paste0(varwithmedian, "_median")], 1, 0)
      target_pop[,paste0(varwithmedian, "_median")] <- 0.5
    }
  }
  
  vartype <- gsub("^.*_","", names(target_pop))
  target_pop <- target_pop[,which(vartype %in% c("mean", "median", "prop"))]
  
  varnames <- gsub("_([^_]+)$","", names(target_pop))
  if(any(duplicated(varnames))){stop("Cannot have more than 1 summary stat for each variable")}
  names(target_pop) <- varnames
  
  # intervention_input is centered using the aggregate data averages.
  intervention_data <- intervention_input
  for(i in varnames){
    intervention_data[,paste0(i, "_centered")] <- intervention_input[,i] - target_pop[,i]
  }
  
  return(list(intervention_data = intervention_data, target_pop = target_pop))
}

