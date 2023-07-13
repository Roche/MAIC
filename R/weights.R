
#' Estimate MAIC propensity weights
#'
#' Estimate propensity weights for matching-adjusted indirect comparison (MAIC).
#'
#' @param intervention_data A data frame containing individual patient data from the intervention study. Intervention_data is assumed to have been preprocessed using preprocess_data (i.e. centered using aggregate data means)
#' @param match_cov A character vector giving the names of the covariates to use in matching. These names must match the column names in intervention_data.
#' @param method The method used for optimisation - The default is method = "BFGS". Refer to \code{\link[stats]{optim}} for options.
#' @param startVal a scalar, the starting value for all coefficients of the propensity score regression
#' @param ... Additional arguments to be passed to optimisation functions such as the method for maximum likelihood optimisation. Refer to \code{\link[stats]{optim}} for options.
#' @return A list with 4 elements
#' \itemize{
#'  \item{wt} a numeric vector of unscaled individual weights.
#'  \item{wt.rs} a numerical vector of rescaled individual weights, with summation equaling to sample size (# rows of input \code{EM})
#'  \item{ess} effective sample size, square of sum divided by sum of squares
#'  \item{opt} R object returned by \code{base::optim()}, for assess convergence and other details
#' }
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#' POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#' THE DECISION SUPPORT UNIT, December 2016
#' @seealso \code{\link{optim}}
#' @export

estimate_weights <- function(intervention_data, match_cov, startVal = 0, 
                             method = "BFGS", ...){
  
  # Check intervention_data is a data frame and match_cov is a character vector
  if(!is.data.frame(intervention_data)){stop("intervention_data is expected to be a data frame")}
  if(!is.character(match_cov)){stop("match_cov is expected to be a character vector")}
  
  # Check if there is any missingness in intervention_data
  missing <- apply(intervention_data[,match_cov], 1, function(x) any(is.na(x)))
  if(any(missing)){
    stop(paste0("Following rows have missing values: ", paste(which(missing), collapse = ",")))
  } 
  
  for(i in match_cov){
    # Check that match_vars is in one of the columns of intervention_data
    if(!paste0(i, "_centered") %in% colnames(intervention_data)){
      stop(paste0("Variable ", i, " is not one of intervention_data column names"))
    }
    
    # Check whether intervention_data has not been centered by the aggregate data means
    # by looking at whether binary variables have only values of 0 and 1 
    if(all(unique(intervention_data[,i]) == 2 & unique(intervention_data[,i]) %in% c(0,1))){
      stop("intervention_data does not seem to be centered by the aggregate data means")
    }
  }
  
  # Objective function
  objfn <- function(a1, X){
    sum(exp(X %*% a1))
  }
  
  # Gradient function
  gradfn <- function(a1, X){
    colSums(sweep(X, 1, exp(X %*% a1), "*"))
  }
  
  # Optimise Q(b) using Newton-Raphson techniques
  opt1 <- stats::optim(par = rep(startVal,dim(intervention_data[,match_cov])[2]),
                       fn = objfn,
                       gr = gradfn,
                       X = as.matrix(intervention_data[,paste0(match_cov,"_centered")]),
                       method = method,
                       control = list(maxit = 300, trace = 2),
                       ...)
  
  alpha <- opt1$par
  wt <- as.vector(exp(as.matrix(intervention_data[,paste0(match_cov,"_centered")]) %*% alpha))
  wt_rs <- (wt / sum(wt)) * nrow(intervention_data)
  
  output <- list(
    wt = wt,
    wt_rs = wt_rs,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1
  )
  return(output)
}

#' Summarize the weight values
#'
#' Produce a summary of the weights (minimum, maximum, median, mean, standard
#' deviation). Mean and standard deviation are provided for completeness.
#' In practice the distribution of weights may be skewed in which case mean and
#' SD should be interpreted with caution.
#'
#' @param weights output created from running \code{\link{estimate_weights}}#'
#' @return A data frame that includes a summary (minimum, maximum, median, mean,
#'   standard deviation) of the weights and rescaled weights.
#' @seealso \code{\link{estimate_weights}}
#' @export

summarize_wts <- function(weights){
  
  with(weights, {
    summary <- data.frame(
      type = c("Weights", "Rescaled weights"),
      mean = c(mean(wt), mean(wt_rs)),
      sd = c(stats::sd(wt), stats::sd(wt_rs)),
      median = c(stats::median(wt), stats::median(wt_rs)),
      min = c(min(wt), min(wt_rs)),
      max = c(max(wt), max(wt_rs))
    )
    return(summary)
  })
}

#' Summarize the weight values
#'
#' Produce a summary of the weights (minimum, maximum, median, mean, standard
#' deviation). Mean and standard deviation are provided for completeness.
#' In practice the distribution of weights may be skewed in which case mean and
#' SD should be interpreted with caution.
#'
#' @param intervention_data intervention data (IPD)
#' @param target_pop preprocessed target_pop data without prefixes
#' @param target_pop_N target population sample size
#' @param weights Weights estimated via \code{\link{estimate_weights}
#' @param covariate_list List of covariates to summarize. List can be larger than what was used in matching.
#' @return Weighted average of covariates for intervention, intervention weighted, and the comparator trial
#' @seealso \code{\link{estimate_weights}}
#' @export

check_weights <- function(intervention_data, target_pop, 
                          target_pop_N = NULL, weights, covariate_list){
  
  if(is.null(target_pop_N){
    stop("Please provide target population N for completeness")
  }
  
  ARM <- c("Intervention", "Intervention_weighted", "Comparator")
  ESS <- round(c(nrow(intervention_data), weights$ess,
                 target_pop$N))
  
  unweighted_cov <- intervention_data %>% summarise_at(covariate_list, list(~ mean(.)))
  weighted_cov <- intervention_data %>% summarise_at(covariate_list, list(~ weighted.mean(., weights$wt)))
  comparator_cov <- select(target_pop, all_of(covariate_list))
  
  cov <- rbind(unweighted_cov, weighted_cov, comparator_cov)
  baseline_summary <- cbind(ARM, ESS, cov)
  
  return(baseline_summary)
}


