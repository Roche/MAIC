#' Bootstrapping for MAIC weighted hazard ratios
#'
#' @param intervention_data  A data frame containing individual patient data
#'   from the intervention study.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'Surv(Time, Event==1) ~ ARM'.
#'   Variable names need to match the corresponding columns in
#'   intervention_data.
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted hazard ratio (HR) from a
#'   Cox proportional hazards model. This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The HR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAIC_example_analysis.R
#'
#' @export
bootstrap_HR <- function(intervention_data, matching, i, model, comparator_data, min_weight = 0.0001){

  # create a visible binding for R CMD check
  wt <- NULL
  
  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching)

  # Give comparator data weights of 1
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")

  # Add the comparator data
  combined_data <- dplyr::bind_rows(perform_wt$analysis_data, comparator_data_wts)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")

  # set weights that are below eta to eta to avoid issues with 0 values
  combined_data$wt <- ifelse(combined_data$wt < min_weight, min_weight, combined_data$wt)
  
  # survival data stat
  cox_model <- survival::coxph(model, data = combined_data, weights = wt)
  HR <- exp(cox_model$coefficients)
}


#' Bootstrapping for MAIC weighted odds ratios
#'
#' @param intervention_data  A data frame containing individual patient data
#'   from the intervention study.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'endpoint ~ treatment_var'. Variable
#'   names need to match the corresponding columns in intervention_data.
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted odds ratio (OR) from a
#'   logistic regression model. This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The OR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAIC_example_analysis.R
#'
#' @export
bootstrap_OR <- function(intervention_data, matching, i, model, comparator_data, min_weight = 0.0001){

  # create a visible binding for R CMD check
  wt <- NULL
  
  # Samples the data
  bootstrap_data <- intervention_data[i,]

  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching)

  # Give comparator data weights of 1
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")

  # Add the comparator data
  combined_data <- dplyr::bind_rows(perform_wt$analysis_data, comparator_data_wts)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")

  # set weights that are below eta to eta to avoid issues with 0 values
  combined_data$wt <- ifelse(combined_data$wt < min_weight, min_weight, combined_data$wt)

  # Perform logistic regression and extract the OR estimate
  logistic.regr <- suppressWarnings(stats::glm(formula = model, family=stats::binomial(link="logit"), data = combined_data, weight = wt))
  OR <- exp(as.numeric(stats::coef(logistic.regr)[2]))
}

#' Bootstrapping for MAIC weighted  difference in means
#'
#' @param intervention_data  A data frame containing individual patient data
#'   from the intervention study.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param contvarname a character value that defines the variable name of the continuous variable to be analyzed
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted difference in means.
#'   This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The weighted mean as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAIC_example_analysis.R
#'
#' @export
bootstrap_means <- function(intervention_data, matching, i, comparator_data, contvarname,
                           min_weight = 0.0001){
  # create a visible binding for R CMD check
  wt <- NULL
  
  
  bootstrap_data <- intervention_data[i, ]
  perform_wt <- estimate_weights(intervention_data = bootstrap_data, 
                                 matching_vars = matching)
  
  # Give comparator data weights of 1
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")
  
  # Add the comparator data
  intervention_data_wts <- perform_wt$analysis_data
  
  # Exclude NA values for contvarname
  intervention_data_wts <- intervention_data_wts[!is.na(intervention_data_wts[,contvarname]),]
  
  combined_data <- dplyr::bind_rows(intervention_data_wts, comparator_data_wts)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")
  combined_data$wt <- ifelse(combined_data$wt < min_weight, 
                             min_weight, combined_data$wt)
  
  mn <- sum(combined_data[,"wt"] * combined_data[,contvarname])/sum(combined_data[,"wt"])
}

#' Bootstrapping for MAIC weighted rate ratios (RR)
#'
#' @param intervention_data  A data frame containing individual patient data
#'   from the intervention study.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param eventsvar Character variable name for number of events
#' @param expvar Character variable name for exposure
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted rates ratio.
#'   This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The rate as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAIC_example_analysis.R
#'
#' @export
bootstrap_rate <- function(intervention_data, matching, i, comparator_data, eventsvar, expvar,
                           min_weight = 0.0001){
  
  bootstrap_data <- intervention_data[i, ]
  perform_wt <- estimate_weights(intervention_data = bootstrap_data, 
                                 matching_vars = matching)
  
  intervention_data_wts <- perform_wt$analysis_data
  
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")
  
  combined_data <- dplyr::bind_rows(intervention_data_wts, comparator_data_wts)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")
  combined_data$wt <- ifelse(combined_data$wt < min_weight, 
                             min_weight, combined_data$wt)
    
  events <- sum(combined_data[,eventsvar]*combined_data[,"wt"])
  exposure <- sum(combined_data[,expvar]*combined_data[,"wt"])
  
  rate <- events/exposure 
}


#' Bootstrapping for MAIC weighted incidence rate from neg binom model to feed into rate ratios
#'
#' @param intervention_data  A data frame containing individual patient data
#'   from the intervention study.
#' @param matching A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.   
#' @param formula formula for the negative binomial regression
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted rates ratio.
#'   This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The IR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @example inst/examples/MAIC_example_analysis.R
#'
#' @export
bootstrap_ir <- function(intervention_data, matching, i, comparator_data, formula,
                           min_weight = 0.0001){
  
  bootstrap_data <- intervention_data[i,]
  # Estimates weights
  perform_wt <- estimate_weights(intervention_data=bootstrap_data, matching_vars=matching)
  
  intervention_data_wts <- perform_wt$analysis_data
  
  comparator_data_wts <- comparator_data %>% dplyr::mutate(wt=1, wt_rs=1, ARM="Comparator")
  
  combined_data <- dplyr::bind_rows(intervention_data_wts, comparator_data_wts)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")
  
  combined_data$wt <- ifelse(combined_data$wt < min_weight, 
                                           min_weight, combined_data$wt)
  
  nbmod_wt_boot <- MASS::glm.nb(as.formula(formula), data=combined_data, weights=wt)
  
  lsm_wt_boot <- emmeans::emmeans(nbmod_wt_boot, "ARM", offset=0, options=list(tran="log"), type="response", data=combined_data)
  
  IR <- summary(lsm_wt_boot)$response[summary(lsm_wt_boot)$ARM!="Comparator"]
  
}
