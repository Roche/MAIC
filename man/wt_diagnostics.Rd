% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{wt_diagnostics}
\alias{wt_diagnostics}
\title{Weight diagnostics}
\usage{
wt_diagnostics(data, wt_col = "wt", wt_rs = "wt_rs", vars)
}
\arguments{
\item{data}{A data frame containing individual patient data from
the intervention study, including a column containing the weights (derived
using estimate_weights).}

\item{wt_col}{The name of the weights column in the data frame containing the
intervention individual patient data and the MAIC propensity weights. The
default is wt.}

\item{wt_rs}{The name of the rescaled weights column in the data frame
containing the intervention individual patient data and the MAIC propensity
weights. The default is wt_rs.}

\item{vars}{A character vector giving the variable names of the baseline
characteristics (not centered). These names must match the column names in
the data.}
}
\value{
List of the following:
\itemize{
  \item The effective sample size (ESS) as a numeric value.
  \item A data frame that includes a summary (minimum, maximum, median, mean,
  standard deviation) of the weights and rescaled weights.
  \item A data frame that includes a summary of patient characteristics
  associated with each weight value.
}
}
\description{
Produce a set of useful diagnostic metrics to summarize propensity weights
\itemize{
  \item ESS (\code{\link{estimate_ess}})
  \item Summary statistics of the weights: minimum, maximum, median, mean, SD (\code{\link{summarize_wts}})
  \item Patient profile associated with weight values (\code{\link{profile_wts}})
}
}
\examples{

# This example code uses the weighted individual patient data, outputted from
# the estimate_weights function to perform weight diagnostics. The weighted data
# is saved within est_weights. To check the weighted aggregate baseline
# characteristics for 'intervention' match those in the comparator data,
# standardized data "target_pop_standard" is used. Please see the package
# vignette for more information on how to use the estimate_weights function and
# derive the "target_pop_standard" data.

library(dplyr)
library(MAIC)

# load est_weights
data(est_weights, package = "MAIC")

# load target_pop_standard
data(target_pop_standard, package = "MAIC")

# List out the uncentered variables used in the matching
match_cov <- c("AGE",
               "SEX",
               "SMOKE",
               "ECOG0")

# Are the weights sensible? ----------------------------------------------------

# The wt_diagnostics function requires the output from the estimate_weights
# function and will output:
# - the effective sample size (ESS)
# - a summary of the weights and rescaled weights (mean, standard deviation,
#   median, minimum and maximum)
# - a unique set of weights with the corresponding patient profile based on the
#   matching variables

diagnostics <- wt_diagnostics(est_weights$analysis_data,
                              vars = match_cov)

diagnostics$ESS
diagnostics$Summary_of_weights
diagnostics$Weight_profiles

# Each of the wt_diagnostics outputs can also be estimated individually
ESS <- estimate_ess(est_weights$analysis_data)
weight_summ <- summarize_wts(est_weights$analysis_data)
wts_profile <- profile_wts(est_weights$analysis_data, vars = match_cov)

# Plot histograms of unscaled and rescaled weights
# bin_width needs to be adapted depending on the sample size in the data set
histogram <- hist_wts(est_weights$analysis_data, bin = 50)
histogram


# Has the optimization worked? -------------------------------------------------

# The following code produces a summary table of the intervention baseline
# characteristics before and after matching compared with the comparator
# baseline characteristics:

# Create an object to hold the output
baseline_summary <- list('Intervention' = NA,
                         'Intervention_weighted' = NA,
                         'Comparator' = NA)

# Summarise matching variables for weighted intervention data
baseline_summary$Intervention_weighted <- est_weights$analysis_data \%>\%
  transmute(AGE, SEX, SMOKE, ECOG0, wt) \%>\%
  summarise_at(match_cov, list(~ weighted.mean(., wt)))

# Summarise matching variables for unweighted intervention data
baseline_summary$Intervention <- est_weights$analysis_data \%>\%
  transmute(AGE, SEX, SMOKE, ECOG0, wt) \%>\%
  summarise_at(match_cov, list(~ mean(.)))

# baseline data for the comparator study
baseline_summary$Comparator <- transmute(target_pop_standard,
                                         AGE,
                                         SEX,
                                         SMOKE,
                                         ECOG0)

# Combine the three summaries
# Takes a list of data frames and binds these together
trt <- names(baseline_summary)
baseline_summary <-  bind_rows(baseline_summary) \%>\%
  transmute_all(sprintf, fmt = "\%.2f") \%>\% #apply rounding for presentation
  transmute(ARM = as.character(trt), AGE, SEX, SMOKE, ECOG0)

# Insert N of intervention  as number of patients
baseline_summary$`N/ESS`[baseline_summary$ARM == "Intervention"] <- nrow(est_weights$analysis_data)

# Insert N for comparator from target_pop_standard
baseline_summary$`N/ESS`[baseline_summary$ARM == "Comparator"] <- target_pop_standard$N

# Insert the ESS as the sample size for the weighted data
# This is calculated above but can also be obtained using the estimate_ess function as shown below
baseline_summary$`N/ESS`[baseline_summary$ARM == "Intervention_weighted"] <- est_weights$analysis_data \%>\%
  estimate_ess(wt_col = 'wt')

baseline_summary <- baseline_summary \%>\%
  transmute(ARM, `N/ESS`=round(`N/ESS`,1), AGE, SEX, SMOKE, ECOG0)


}
\seealso{
\code{\link{estimate_weights}}, \code{\link{estimate_ess}}, \code{\link{summarize_wts}}, \code{\link{profile_wts}}
}
