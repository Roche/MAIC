
#' Estimate MAIC propensity weights
#'
#' Estimate propensity weights for matching-adjusted indirect comparison (MAIC).
#'
#' @param intervention_data A data frame containing individual patient data from 
#' the intervention study. Intervention_data is assumed to have been preprocessed using 
#' preprocess_data (i.e. centered using aggregate data means)
#' @param match_cov A character vector giving the names of the covariates to
#' use in matching. These names must match the column names in intervention_data.
#' @param method The method used for optimisation - The default is method =
#' "BFGS". Refer to \code{\link[stats]{optim}} for options.
#' @param startVal a scalar, the starting value for all coefficients of the propensity score 
#' regression
#' @param ... Additional arguments to be passed to optimisation functions such
#' as the method for maximum likelihood optimisation. Refer to \code{\link[stats]{optim}} 
#' for options.
#' @return A list with 4 elements
#'  \item{wt} a numeric vector of unscaled individual weights.
#'  \item{wt.rs} a numerical vector of rescaled individual weights, with summation equaling to sample size (# rows of input \code{EM})
#'  \item{ess} effective sample size, square of sum divided by sum of squares
#'  \item{opt} R object returned by \code{base::optim()}, for assess convergence and other details
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
  
  # Check if match_cov name is included in the IPD data
  
  
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



# Functions for summarizing the weights ----------------------------------------

#' Estimate effective sample size
#'
#' Estimate the effective sample size (ESS).
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#'
#' @details For a weighted estimate, the effective sample size (ESS) is the
#'   number of independent non-weighted individuals that would be required to
#'   give an estimate with the same precision as the weighted sample estimate. A
#'   small ESS, relative to the original sample size, is an indication that the
#'   weights are highly variable and that the estimate may be unstable. This
#'   often occurs if there is very limited overlap in the distribution of the
#'   matching variables between the populations being compared. If there is
#'   insufficient overlap between populations it may not be possible to obtain
#'   reliable estimates of the weights
#'
#' @return The effective sample size (ESS) as a numeric value.
#'
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @seealso \code{\link{estimate_weights}}
#'
#' @example inst/examples/MAIC_example_weight_diagnostics.R
#'
#' @export
estimate_ess <- function(data, wt_col="wt"){
  ess <- sum(data[,wt_col])^2/sum(data[,wt_col]^2)
  return(ess)
}

#' Summarize the weight values
#'
#' Produce a summary of the weights (minimum, maximum, median, mean, standard
#' deviation). Mean and standard deviation are provided for completeness.
#' In practice the distribution of weights may be skewed in which case mean and
#' SD should be interpreted with caution.
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#'
#' @return A data frame that includes a summary (minimum, maximum, median, mean,
#'   standard deviation) of the weights and rescaled weights.
#'
#' @seealso \code{\link{estimate_weights}}
#'
#' @example inst/examples/MAIC_example_weight_diagnostics.R
#'
#' @export
summarize_wts <- function(data, wt_col="wt", rs_wt_col="wt_rs"){
  summary <- data.frame(
    type = c("Weights", "Rescaled weights"),
    mean = c(mean(data[,wt_col]), mean(data[,rs_wt_col])),
    sd = c(stats::sd(data[,wt_col]), stats::sd(data[,rs_wt_col])),
    median = c(stats::median(data[,wt_col]), stats::median(data[,rs_wt_col])),
    min = c(min(data[,wt_col]), min(data[,rs_wt_col])),
    max = c(max(data[,wt_col]), max(data[,rs_wt_col]))
  )
  return(summary)
}


#' Produce histograms of weights and rescaled weights
#'
#' Produce a plot containing two histograms (one of the weights and one of the
#' rescaled weights).
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param rs_wt_col The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#' @param bin Number of bins to plot histogram. The default is 30.
#'
#' @return A histogram plot of the weights and rescaled weights.
#'
#' @seealso \code{\link{estimate_weights}}
#'
#' @example inst/examples/MAIC_example_weight_diagnostics.R
#'
#' @export
hist_wts <- function(data, wt_col="wt", rs_wt_col="wt_rs", bin = 30) {

  # create visible local bindings for R CMD check
  value <- `Rescaled weights` <- Weights <- NULL
  
  wt_data1 <- data %>%
    dplyr::select(tidyselect::all_of(c(wt_col, rs_wt_col))) %>% # select only the weights and rescaled weights
    dplyr::rename("Weights" = tidyselect::all_of(wt_col), "Rescaled weights" = tidyselect::all_of(rs_wt_col))  # rename so for plots
  
  # tidyr::gather() # weights and rescaled weights in one column for plotting
  
  wt_data <- dplyr::bind_rows(
    dplyr::transmute(wt_data1, key = "Weights", value = Weights),
    dplyr::transmute(wt_data1, key = "Rescaled weights", value = `Rescaled weights`)
  )

  hist_plot <- ggplot2::ggplot(wt_data) + ggplot2::geom_histogram(ggplot2::aes(value), bins = bin) +
    ggplot2::facet_wrap(~key,  ncol=1) + # gives the two plots (one on top of the other)
    ggplot2::theme_bw()+
    ggplot2::theme(axis.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 16)) +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Weight")

  return(hist_plot)
}


#' Produce a data frame of the weights assigned to patient profiles
#'
#' Select the patient characteristics used in the matching and the MAIC weights
#' and output a data frame of unique propensity weight values with the
#' associated summary baseline characteristics. This data frame helps to
#' understand how different patient profiles are contributing to the analyses by
#' illustrating the patient characteristics associated with different weight
#' values. For example, min, max and median weights. This function is most
#' useful when only matching on binary variables as there are fewer unique
#' values.
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using \code{\link{estimate_weights}}).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param wt_rs The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#' @param vars A character vector giving the variable names of the baseline
#'   characteristics (not centered). These names must match the column names in
#'   the data.
#'
#' @return A data frame that includes a summary of patient characteristics
#'   associated with each weight value.
#'
#' @seealso \code{\link{estimate_weights}}
#' @example inst/examples/MAIC_example_weight_diagnostics.R
#' @export

profile_wts <- function(data, wt_col="wt", wt_rs="wt_rs", vars){
  profile_data <-  data %>%
    dplyr::select(tidyselect::all_of(vars), tidyselect::all_of(wt_col), tidyselect::all_of(wt_rs))

  profile_wts <- profile_data %>%
    dplyr::distinct()

  return(profile_wts)
}

#' Weight diagnostics
#'
#' Produce a set of useful diagnostic metrics to summarize propensity weights
#' \itemize{
#'   \item ESS (\code{\link{estimate_ess}})
#'   \item Summary statistics of the weights: minimum, maximum, median, mean, SD (\code{\link{summarize_wts}})
#'   \item Patient profile associated with weight values (\code{\link{profile_wts}})
#' }
#'
#' @param data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param wt_col The name of the weights column in the data frame containing the
#'   intervention individual patient data and the MAIC propensity weights. The
#'   default is wt.
#' @param wt_rs The name of the rescaled weights column in the data frame
#'   containing the intervention individual patient data and the MAIC propensity
#'   weights. The default is wt_rs.
#' @param vars A character vector giving the variable names of the baseline
#'   characteristics (not centered). These names must match the column names in
#'   the data.
#'
#' @return List of the following:
#' \itemize{
#'   \item The effective sample size (ESS) as a numeric value.
#'   \item A data frame that includes a summary (minimum, maximum, median, mean,
#'   standard deviation) of the weights and rescaled weights.
#'   \item A data frame that includes a summary of patient characteristics
#'   associated with each weight value.
#' }
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{estimate_ess}}, \code{\link{summarize_wts}}, \code{\link{profile_wts}}
#' @example inst/examples/MAIC_example_weight_diagnostics.R
#' @export

wt_diagnostics <- function(data, wt_col="wt", wt_rs="wt_rs", vars){

  # ESS
  ESS <- estimate_ess(data, wt_col)

  # Summary
  summ_wts <- summarize_wts(data, wt_col, wt_rs)

  # Weight profiles
  profile <- profile_wts(data, wt_col, wt_rs, vars)

  output <- list("ESS" = ESS,
                 "Summary_of_weights" = summ_wts,
                 "Weight_profiles" = profile
  )
  return(output)
}

#' Checking whether optimization has worked
#'
#' Convenient function to check whether the re-weighted baseline characteristics for
#' the intervention-treated patients match those aggregate characteristics from the
#' comparator trial and outputs a summary that can be used for reporting
#' 
#' @param analysis_data A data frame containing individual patient data from
#'   the intervention study, including a column containing the weights (derived
#'   using estimate_weights).
#' @param match_covs A character vector giving the names of the covariates that were used to estimate weights 
#' @param target_pop_standard aggregate characteristics of the comparator trial with the same naming as the analysis_data
#' @example inst/examples/MAIC_example_weight_diagnostics.R
#' @seealso \code{\link{estimate_weights}}
#' @return Summary of patient characteristics before and after matching, including ESS and comparator trial aggregate summary
#' @export

check_weights <- function(analysis_data = NULL, matching_vars = NULL, 
                          target_pop_standard = NULL){
  
  ARM <- c("Intervention", "Intervention_weighted", "Comparator")
  ESS <- round(c(nrow(analysis_data), estimate_ess(analysis_data),
                 target_pop_standard$N))
  
  weighted_cov <- analysis_data %>% summarise_at(match_cov, list(~ weighted.mean(., wt)))
  unweighted_cov <-  analysis_data %>% summarise_at(match_cov, list(~ mean(.)))
  comparator_cov <- select(target_pop_standard, all_of(match_cov))
  
  cov <- rbind(unweighted_cov, weighted_cov, comparator_cov)
  baseline_summary <- cbind(ARM, ESS, cov)
  
  return(baseline_summary)
}