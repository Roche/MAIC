
# This example code combines weighted individual patient data for 'intervention'
# and pseudo individual patient data for 'comparator' and performs analyses for
# two endpoints: overall survival (a time to event outcome) and objective
# response (a binary outcome)

library(dplyr)
library(boot)
library(survival)
library(MAIC)
library(ggplot2)
library(survminer)
library(flextable)
library(officer)

# load intervention data with weights saved in est_weights
data(est_weights, package = "MAIC")

intervention_mean <- 11
intervention_sd <- 4
est_weights$analysis_data$diff_from_baseline <- rnorm(nrow(est_weights$analysis_data),intervention_mean,intervention_sd)

est_weights$analysis_data$exp <- est_weights$analysis_data$Time


# Combine data -----------------------------------------------------------------

# Combine the the comparator pseudo data with the analysis data, outputted from
# the estimate_weights function

# Read in digitised pseudo survival data for the comparator
comparator_surv <- read.csv(system.file("extdata", "psuedo_IPD.csv",
                                        package = "MAIC", mustWork = TRUE)) %>%
                            rename(Time=Time, Event=Event)


# Simulate response data based on the known proportion of responders
comparator_n <- nrow(comparator_surv) # total number of patients in the comparator data
comparator_prop_events <- 0.4 # proportion of responders
# Calculate number with event
# Use round() to ensure we end up with a whole number of people
# number without an event = Total N - number with event to ensure we keep the
# same number of patients
n_with_event <- round(comparator_n*comparator_prop_events, digits = 0)
comparator_binary <- data.frame("response"= c(rep(1, n_with_event), rep(0, comparator_n - n_with_event)))

comparator_mean <- 10
comparator_sd <- 4
comparator_cont <- data.frame("diff_from_baseline"= rnorm(comparator_n,comparator_mean,comparator_sd))

comparator_exp <- rexp(comparator_n,0.000018)
comparator_events <- n_with_event
comparator_exposure <- data.frame("exp"= comparator_exp)



# Join survival and response comparator data
# Note the rows do not represent observations from a particular patient
# i.e. there is no relationship between the survival time and response status
# for a given row since this is simulated data
# In a real data set this relationship would be present
comparator_input <- cbind(comparator_surv, comparator_binary, comparator_cont,comparator_exposure) %>%
  mutate(wt=1, wt_rs=1, ARM="Comparator") # All patients have weight = 1
head(comparator_input)

# Join comparator data with the intervention data (naming of variables should be
# consistent between both datasets)
# Set factor levels to ensure "Comparator" is the reference treatment
combined_data <-  bind_rows(est_weights$analysis_data, comparator_input)
combined_data$ARM <- relevel(as.factor(combined_data$ARM), ref="Comparator")


#### Estimating the relative effect --------------------------------------------

### Example for survival data --------------------------------------------------

## Kaplan-Meier plot

# Unweighted intervention data
KM_int <- survfit(formula = Surv(Time, Event==1) ~ 1 ,
                  data = est_weights$analysis_data,
                  type="kaplan-meier")

# Weighted intervention data
KM_int_wtd <- survfit(formula = Surv(Time, Event==1) ~ 1 ,
                      data = est_weights$analysis_data,
                      weights = wt,
                      type="kaplan-meier")

# Comparator data
KM_comp <- survfit(formula = Surv(Time, Event==1) ~ 1 ,
                   data = comparator_input,
                   type="kaplan-meier")

# Combine the survfit objects ready for ggsurvplot
KM_list <- list(Intervention = KM_int,
                Intervention_weighted = KM_int_wtd,
                Comparator = KM_comp)

#Produce the Kaplan-Meier plot
KM_plot <- ggsurvplot(KM_list,
                      combine = TRUE,
                      risk.table=TRUE, # numbers at risk displayed on the plot
                      break.x.by=50,
                      xlab="Time (days)",
                      censor=FALSE,
                      legend.title = "Treatment",
                      title = "Kaplan-Meier plot of overall survival",
                      legend.labs=c("Intervention",
                                    "Intervention weighted",
                                    "Comparator"),
                      font.legend = list(size = 10)) +
  guides(colour=guide_legend(nrow = 2))


## Estimating the hazard ratio (HR)

# Fit a Cox model without weights to estimate the unweighted HR
unweighted_cox <- coxph(Surv(Time, Event==1) ~ ARM, data = combined_data)

HR_CI_cox <- summary(unweighted_cox)$conf.int %>%
  as.data.frame() %>%
  transmute("HR" = `exp(coef)`,
            "HR_low_CI" = `lower .95`,
            "HR_upp_CI" = `upper .95`)

# Fit a Cox model with weights to estimate the weighted HR
weighted_cox <- coxph(Surv(Time, Event==1) ~ ARM,
                      data = combined_data,
                      weights = wt)

HR_CI_cox_wtd <- summary(weighted_cox)$conf.int %>%
  as.data.frame() %>%
  transmute("HR" = `exp(coef)`,
            "HR_low_CI" = `lower .95`,
            "HR_upp_CI" = `upper .95`)

## Bootstrap the confidence interval of the weighted HR

HR_bootstraps <- boot(data = est_weights$analysis_data, # intervention data
                      statistic = bootstrap_HR, # bootstrap the HR (defined in the MAIC package)
                      R=1000, # number of bootstrap samples
                      comparator_data = comparator_input, # comparator pseudo data
                      matching = est_weights$matching_vars, # matching variables
                      model = Surv(Time, Event==1) ~ ARM # model to fit
  )

## Bootstrapping diagnostics
# Summarize bootstrap estimates in a histogram
# Vertical lines indicate the median and upper and lower CIs
hist(HR_bootstraps$t, main = "", xlab = "Boostrapped HR")
abline(v= quantile(HR_bootstraps$t, probs = c(0.025, 0.5, 0.975)), lty=2)

# Median of the bootstrap samples
HR_median <- median(HR_bootstraps$t)

# Bootstrap CI - Percentile CI
boot_ci_HR <- boot.ci(boot.out = HR_bootstraps, index=1, type="perc")

# Bootstrap CI - BCa CI
boot_ci_HR_BCA <- boot.ci(boot.out = HR_bootstraps, index=1, type="bca")

## Summary

# Produce a summary of HRs and CIs
HR_summ <-  rbind(HR_CI_cox, HR_CI_cox_wtd) %>% # Unweighted and weighted HRs and CIs from Cox models
  mutate(Method = c("HR (95% CI) from unadjusted Cox model",
                    "HR (95% CI) from weighted Cox model")) %>%

  # Median bootstrapped HR and 95% percentile CI
  rbind(data.frame("HR" = HR_median,
                   "HR_low_CI" = boot_ci_HR$percent[4],
                   "HR_upp_CI" = boot_ci_HR$percent[5],
                   "Method"="Bootstrap median HR (95% percentile CI)")) %>%

  # Median bootstrapped HR and 95% bias-corrected and accelerated bootstrap CI
  rbind(data.frame("HR" = HR_median,
                   "HR_low_CI" = boot_ci_HR_BCA$bca[4],
                   "HR_upp_CI" = boot_ci_HR_BCA$bca[5],
                   "Method"="Bootstrap median HR (95% BCa CI)")) %>%
  #apply rounding for numeric columns
  mutate_if(.predicate = is.numeric, sprintf, fmt = "%.3f") %>%
  #format for output
  transmute(Method, HR_95_CI = paste0(HR, " (", HR_low_CI, " to ", HR_upp_CI, ")"))

# Summarize the results in a table suitable for word/ powerpoint
HR_table <- HR_summ %>%
    regulartable() %>% #make it a flextable object
    set_header_labels(Method = "Method",  HR_95_CI = "Hazard ratio (95% CI)") %>%
    font(font = 'Arial', part = 'all') %>%
    fontsize(size = 14, part = 'all') %>%
    bold(part = 'header') %>%
    align(align = 'center', part = 'all') %>%
    align(j = 1, align = 'left', part = 'all') %>%
    border_outer(border = fp_border()) %>%
    border_inner_h(border = fp_border()) %>%
    border_inner_v(border = fp_border()) %>%
    autofit(add_w = 0.2, add_h = 2)


### Example for response data --------------------------------------------------

## Estimating the odds ratio (OR)

# Fit a logistic regression model without weights to estimate the unweighted OR
unweighted_OR <- glm(formula = response~ARM,
                     family = binomial(link="logit"),
                     data = combined_data)

# Log odds ratio
log_OR_CI_logit <- cbind(coef(unweighted_OR), confint.default(unweighted_OR, level = 0.95))[2,]

# Exponentiate to get Odds ratio
OR_CI_logit <- exp(log_OR_CI_logit)
#tidy up naming
names(OR_CI_logit) <- c("OR", "OR_low_CI", "OR_upp_CI")

# Fit a logistic regression model with weights to estimate the weighted OR
weighted_OR <- suppressWarnings(glm(formula = response~ARM,
                                    family = binomial(link="logit"),
                                    data = combined_data,
                                    weight = wt))

# Weighted log odds ratio
log_OR_CI_logit_wtd <- cbind(coef(weighted_OR), confint.default(weighted_OR, level = 0.95))[2,]

# Exponentiate to get weighted odds ratio
OR_CI_logit_wtd <- exp(log_OR_CI_logit_wtd)
#tidy up naming
names(OR_CI_logit_wtd) <- c("OR", "OR_low_CI", "OR_upp_CI")

## Bootstrap the confidence interval of the weighted OR
OR_bootstraps <- boot(data = est_weights$analysis_data, # intervention data
                      statistic = bootstrap_OR, # bootstrap the OR
                      R = 1000, # number of bootstrap samples
                      comparator_data = comparator_input, # comparator pseudo data
                      matching = est_weights$matching_vars, # matching variables
                      model = 'response ~ ARM' # model to fit
                      )

## Bootstrapping diagnostics
# Summarize bootstrap estimates in a histogram
# Vertical lines indicate the median and upper and lower CIs
hist(OR_bootstraps$t, main = "", xlab = "Boostrapped OR")
abline(v= quantile(OR_bootstraps$t, probs = c(0.025,0.5,0.975)), lty=2)

# Median of the bootstrap samples
OR_median <- median(OR_bootstraps$t)

# Bootstrap CI - Percentile CI
boot_ci_OR <- boot.ci(boot.out = OR_bootstraps, index=1, type="perc")

# Bootstrap CI - BCa CI
boot_ci_OR_BCA <- boot.ci(boot.out = OR_bootstraps, index=1, type="bca")

## Summary
# Produce a summary of ORs and CIs
OR_summ <- rbind(OR_CI_logit, OR_CI_logit_wtd) %>% # Unweighted and weighted ORs and CIs
  as.data.frame() %>%
  mutate(Method = c("OR (95% CI) from unadjusted logistic regression model",
                    "OR (95% CI) from weighted logistic regression model")) %>%

  # Median bootstrapped HR and 95% percentile CI
  rbind(data.frame("OR" = OR_median,
                   "OR_low_CI" = boot_ci_OR$percent[4],
                   "OR_upp_CI" = boot_ci_OR$percent[5],
                   "Method"="Bootstrap median HR (95% percentile CI)")) %>%

  # Median bootstrapped HR and 95% bias-corrected and accelerated bootstrap CI
  rbind(data.frame("OR" = OR_median,
                   "OR_low_CI" = boot_ci_OR_BCA$bca[4],
                   "OR_upp_CI" = boot_ci_OR_BCA$bca[5],
                   "Method"="Bootstrap median HR (95% BCa CI)")) %>%
  #apply rounding for numeric columns
  mutate_if(.predicate = is.numeric, sprintf, fmt = "%.3f") %>%
  #format for output
  transmute(Method, OR_95_CI = paste0(OR, " (", OR_low_CI, " to ", OR_upp_CI, ")"))

# turns the results to a table suitable for word/ powerpoint
OR_table <- OR_summ %>%
            regulartable() %>% #make it a flextable object
            set_header_labels(Method = "Method",  OR_95_CI = "Odds ratio (95% CI)")  %>%
            font(font = 'Arial', part = 'all') %>%
            fontsize(size = 14, part = 'all') %>%
            bold(part = 'header') %>%
            align(align = 'center', part = 'all') %>%
            align(j = 1, align = 'left', part = 'all') %>%
            border_outer(border = fp_border()) %>%
            border_inner_h(border = fp_border()) %>%
            border_inner_v(border = fp_border()) %>%
            autofit(add_w = 0.2)

### Example for continuous data --------------------------------------------------

comparator_n <- comparator_n
comparator_mean <-comparator_mean
comparator_sd <- comparator_sd
comparator_se <- comparator_sd/sqrt(comparator_n)
comparator_lcl <- comparator_mean - qnorm(0.975)*comparator_se
comparator_ucl <- comparator_mean + qnorm(0.975)*comparator_se
comparator_summary <- paste0(round(comparator_mean, 2), " (", round(comparator_lcl, 2), ", ", round(comparator_ucl, 2), ")")


intervention_n <- nrow(est_weights$analysis_data)
intervention_mean <- mean(est_weights$analysis_data[,"diff_from_baseline"])
intervention_sd <- sd(est_weights$analysis_data[,"diff_from_baseline"])
intervention_se <- intervention_sd/sqrt(intervention_n)
intervention_lcl <- intervention_mean - qnorm(0.975)*intervention_se
intervention_ucl <- intervention_mean + qnorm(0.975)*intervention_se
intervention_summary <- paste0(round(intervention_mean, 2), " (", round(intervention_lcl, 2), ", ", round(intervention_ucl, 2), ")")


intervention_n_w <- sum(est_weights$analysis_data[,"wt"])
intervention_mean_w <-     sum(est_weights$analysis_data[,"wt"] * est_weights$analysis_data[,"diff_from_baseline"])/sum(est_weights$analysis_data[,"wt"])
intervention_sd_w <- sqrt(Hmisc::wtd.var(est_weights$analysis_data[,"diff_from_baseline"], w=est_weights$analysis_data[,"wt"])) 
intervention_se_w <- intervention_sd_w/sqrt(intervention_n_w)
intervention_lcl_w <- intervention_mean_w - qnorm(0.975)*intervention_se_w
intervention_ucl_w <- intervention_mean_w + qnorm(0.975)*intervention_se_w
intervention_summary_w <- paste0(round(intervention_mean_w, 2), " (", round(intervention_lcl_w, 2), ", ", round(intervention_ucl_w, 2), ")")


## Compare arms

diff_u <- intervention_mean - comparator_mean
se_diff_u <- sqrt(intervention_se^2 + comparator_se^2)
lcl_diff_u <- diff_u-qnorm(0.975)*se_diff_u
ucl_diff_u <- diff_u+qnorm(0.975)*se_diff_u
pval_u <- 2*min(pnorm(as.numeric(diff_u/se_diff_u)), 1-pnorm(as.numeric(diff_u/se_diff_u)))

diff_w <- intervention_mean_w - comparator_mean
se_diff_w <- sqrt(intervention_se_w^2 + comparator_se^2)
lcl_diff_w <- diff_w-qnorm(0.975)*se_diff_w
ucl_diff_w <- diff_w+qnorm(0.975)*se_diff_w
pval_w <- 2*min(pnorm(as.numeric(diff_w/se_diff_w)), 1-pnorm(as.numeric(diff_w/se_diff_w)))


## Bootstrapping mean in intervention arm
mean_bootstraps <- boot(data = est_weights$analysis_data,
                        comparator_data = combined_data%>% filter(ARM=="Comparator"), 
                        statistic = bootstrap_means, # bootstrap the exp arm mean
                        R = 1000, # number of bootstrap samples
                        matching = est_weights$matching_vars, # matching variables
                        contvarname="diff_from_baseline"
)

# Bootstrap standard error of mean in intervention arm and use it to calculate bootstrap CI and p-value for difference

mean_bootvar <- var(mean_bootstraps$t)

se_diff_w_boot <- sqrt(mean_bootvar + comparator_se^2)
lcl_diff_w_boot <- diff_w-qnorm(0.975)*se_diff_w_boot
ucl_diff_w_boot <- diff_w+qnorm(0.975)*se_diff_w_boot
pval_w_boot <- 2*min(pnorm(as.numeric(diff_w/se_diff_w_boot)), 1-pnorm(as.numeric(diff_w/se_diff_w_boot)))


## Summarise and tabulate

diff_summary_standard   <- paste0(round(diff_u,2), " (", round(lcl_diff_u,2), ", ", round(ucl_diff_u, 2), ")")
diff_summary_standard_w <- paste0(round(diff_w,2), " (", round(lcl_diff_w,2), ", ", round(ucl_diff_w, 2), ")")

diff_summary_boot_w <- paste0(round(diff_w,2), " (", round(lcl_diff_w_boot,2), ", ", round(ucl_diff_w_boot, 2), ")")



## Example for rates data ----------------------------------------------------------------------------------------------------------------------------------


comparator_events <- comparator_events
comparator_exp <- flexsurv::mean_exp(0.00005)
comparator_n <- comparator_n
comparator_rate <- comparator_events/comparator_exp
comparator_summary <- paste0(round(comparator_rate, 4), " (", comparator_events, "/", round(comparator_exp,1), ")")

intervention_events <- sum(est_weights$analysis_data[,"response"])
intervention_exp <- sum(est_weights$analysis_data[,"exp"])
intervention_n <- nrow(est_weights$analysis_data)
intervention_rate <- intervention_events/intervention_exp
intervention_summary <- paste0(round(intervention_rate, 4), " (", intervention_events, "/", round(intervention_exp,1), ")")

intervention_events_w <- sum(est_weights$analysis_data[,"response"]*est_weights$analysis_data[,"wt"])
intervention_exp_w <- sum(est_weights$analysis_data[,"exp"]*est_weights$analysis_data[,"wt"])
intervention_n_w <- sum(est_weights$analysis_data[,"wt"])
intervention_rate_w <- intervention_events_w/intervention_exp_w
intervention_summary_w <- paste0(round(intervention_rate_w, 4), " (", round(intervention_events_w,1), "/", round(intervention_exp_w,1), ")")


## Compare arms

rr_u <- intervention_rate / comparator_rate
lrr_u <- log(rr_u)
se_lrr_u <- sqrt(1/intervention_events + 1/comparator_events)
lcl_rr_u <- exp(lrr_u-qnorm(0.975)*se_lrr_u)
ucl_rr_u <- exp(lrr_u+qnorm(0.975)*se_lrr_u)
pval_u <- 2*min(pnorm(as.numeric(lrr_u/se_lrr_u)), 1-pnorm(as.numeric(lrr_u/se_lrr_u)))

rr_w <- intervention_rate_w / comparator_rate
lrr_w <- log(rr_w)
se_lrr_w <- sqrt(1/intervention_events_w + 1/comparator_events)
lcl_rr_w <- exp(lrr_w-qnorm(0.975)*se_lrr_w)
ucl_rr_w <- exp(lrr_w+qnorm(0.975)*se_lrr_w)
pval_w <- 2*min(pnorm(as.numeric(lrr_w/se_lrr_w)), 1-pnorm(as.numeric(lrr_w/se_lrr_w)))


## Bootstrapping rate in intervention arm


rate_bootstraps <- boot(data = est_weights$analysis_data, 
                        statistic = bootstrap_rate, # bootstrap the exp arm rate
                        R = 1000, # number of bootstrap samples
                        comparator_data = combined_data%>% filter(ARM=="Comparator"), 
                        matching = est_weights$matching_vars, # matching variables
                        eventsvar="response",
                        expvar="exp"
)

# Bootstrap standard error of log rate in intervention arm and use it to calculate bootstrap CI and p-value for RR

lograte_bootvar <- var(log(rate_bootstraps$t))

se_lrr_w_boot <- sqrt(lograte_bootvar + 1/comparator_events)
lcl_rr_w_boot <- exp(lrr_w-qnorm(0.975)*se_lrr_w_boot)
ucl_rr_w_boot <- exp(lrr_w+qnorm(0.975)*se_lrr_w_boot)
pval_w_boot <- 2*min(pnorm(as.numeric(lrr_w/se_lrr_w_boot)), 1-pnorm(as.numeric(lrr_w/se_lrr_w_boot)))

## Summarise and tabulate

rr_summary_standard   <- paste0(round(rr_u,3), " (", round(lcl_rr_u,3), ", ", round(ucl_rr_u, 3), ")")
rr_summary_standard_w <- paste0(round(rr_w,3), " (", round(lcl_rr_w,3), ", ", round(ucl_rr_w, 3), ")")
rr_summary_boot_w <- paste0(round(rr_w,3), " (", round(lcl_rr_w_boot,3), ", ", round(ucl_rr_w_boot, 3), ")")




# Example for incidence rate to rr ------------------------------------------------------------------------------------------------------------------------

comparator_arr <- comparator_rate
comparator_n <- comparator_n
comparator_dur_yrs <- comparator_exp/comparator_n/365.25
comparator_durtot <- comparator_dur_yrs*comparator_n
comparator_etot <- comparator_n * (comparator_exp/365.25)/(comparator_n*0.9)
comparator_rtot_from_dur <- comparator_arr*comparator_durtot
comparator_rtot_from_etot <- comparator_arr*comparator_etot
comparator_summary_dur <- paste0(round(comparator_arr, 3), " (", round(comparator_rtot_from_dur, 1), "/", round(comparator_durtot, 1), ")")
comparator_summary_etot <- paste0(round(comparator_arr, 3), " (", round(comparator_rtot_from_etot, 1), "/", round(comparator_etot, 1), ")")

intervention_arr <- intervention_rate
intervention_n <- nrow(est_weights$analysis_data)
intervention_dur_yrs <- max(est_weights$analysis_data$Time)
intervention_durtot <- intervention_dur_yrs*intervention_n
intervention_etot <- sum(est_weights$analysis_data$exp)
intervention_rtot_from_dur <- intervention_arr*intervention_durtot
intervention_rtot_from_etot <- intervention_arr*intervention_etot
intervention_summary_dur <- paste0(round(intervention_arr, 3), " (", round(intervention_rtot_from_dur, 1), "/", round(intervention_durtot, 1), ")")
intervention_summary_etot <- paste0(round(intervention_arr, 3), " (", round(intervention_rtot_from_etot, 1), "/", round(intervention_etot, 1), ")")

intervention_arr_w <- intervention_rate_w
intervention_n_w <- intervention_n_w
intervention_dur_yrs_w <- max(est_weights$analysis_data$Time*est_weights$analysis_data$wt)
intervention_durtot_w <- intervention_dur_yrs_w*intervention_n_w
intervention_etot_w <- sum(est_weights$analysis_data$exp * est_weights$analysis_data$wt)
intervention_rtot_from_dur_w <- intervention_arr_w*intervention_durtot_w
intervention_rtot_from_etot_w <- intervention_arr_w*intervention_etot_w
intervention_summary_dur_w <- paste0(round(intervention_arr_w, 3), " (", round(intervention_rtot_from_dur_w, 1), "/", round(intervention_durtot_w, 1), ")")
intervention_summary_etot_w <- paste0(round(intervention_arr_w, 3), " (", round(intervention_rtot_from_etot_w, 1), "/", round(intervention_etot_w, 1), ")")

## Compare arms

rr_u <- intervention_arr / comparator_arr
lrr_u <- log(rr_u)
se_lrr_dur_u <- sqrt(1/intervention_rtot_from_dur + 1/comparator_rtot_from_dur)
lcl_rr_dur_u <- exp(lrr_u-qnorm(0.975)*se_lrr_dur_u)
ucl_rr_dur_u <- exp(lrr_u+qnorm(0.975)*se_lrr_dur_u)
se_lrr_etot_u <- sqrt(1/intervention_rtot_from_etot + 1/comparator_rtot_from_etot)
lcl_rr_etot_u <- exp(lrr_u-qnorm(0.975)*se_lrr_etot_u)
ucl_rr_etot_u <- exp(lrr_u+qnorm(0.975)*se_lrr_etot_u)
pval_dur_u <- 2*min(pnorm(as.numeric(lrr_u/se_lrr_dur_u)), 1-pnorm(as.numeric(lrr_u/se_lrr_dur_u)))
pval_etot_u <- 2*min(pnorm(as.numeric(lrr_u/se_lrr_etot_u)), 1-pnorm(as.numeric(lrr_u/se_lrr_etot_u)))

rr_w <- intervention_arr_w / comparator_arr
lrr_w <- log(rr_w)
se_lrr_dur_w <- sqrt(1/intervention_rtot_from_dur_w + 1/comparator_rtot_from_dur)
lcl_rr_dur_w <- exp(lrr_w-qnorm(0.975)*se_lrr_dur_w)
ucl_rr_dur_w <- exp(lrr_w+qnorm(0.975)*se_lrr_dur_w)
se_lrr_etot_w <- sqrt(1/intervention_rtot_from_etot_w + 1/comparator_rtot_from_etot)
lcl_rr_etot_w <- exp(lrr_w-qnorm(0.975)*se_lrr_etot_w)
ucl_rr_etot_w <- exp(lrr_w+qnorm(0.975)*se_lrr_etot_w)
pval_dur_w <- 2*min(pnorm(as.numeric(lrr_w/se_lrr_dur_w)), 1-pnorm(as.numeric(lrr_w/se_lrr_dur_w)))
pval_etot_w <- 2*min(pnorm(as.numeric(lrr_w/se_lrr_etot_w)), 1-pnorm(as.numeric(lrr_w/se_lrr_etot_w)))


## Bootstrapping incidence rate in intervention arm

IR_bootstraps <- boot(data = est_weights$analysis_data, 
                      statistic = bootstrap_ir, # bootstrap the exp arm incidence rate
                      R = 1000, # number of bootstrap samples
                      comparator_data = combined_data%>% filter(ARM=="Comparator"), 
                      formula = "response~ARM+offset(log(exp))",
                      matching = est_weights$matching_vars # matching variables
)

# Bootstrap standard error of log IR and calculate RR CI and p-value
IR_bootvar <- var(log(IR_bootstraps$t))

se_lrr_dur_w_boot <- sqrt(IR_bootvar + 1/comparator_rtot_from_dur)
lcl_rr_dur_w_boot <- exp(lrr_w-qnorm(0.975)*se_lrr_dur_w_boot)
ucl_rr_dur_w_boot <- exp(lrr_w+qnorm(0.975)*se_lrr_dur_w_boot)
se_lrr_etot_w_boot <- sqrt(IR_bootvar + 1/comparator_rtot_from_etot)
lcl_rr_etot_w_boot <- exp(lrr_w-qnorm(0.975)*se_lrr_etot_w_boot)
ucl_rr_etot_w_boot <- exp(lrr_w+qnorm(0.975)*se_lrr_etot_w_boot)
pval_dur_w_boot <- 2*min(pnorm(as.numeric(lrr_w/se_lrr_dur_w_boot)), 1-pnorm(as.numeric(lrr_w/se_lrr_dur_w_boot)))
pval_etot_w_boot <- 2*min(pnorm(as.numeric(lrr_w/se_lrr_etot_w_boot)), 1-pnorm(as.numeric(lrr_w/se_lrr_etot_w_boot)))


# Summarise and tabulate

rr_summary_standard   <- paste0(round(rr_u,3), " (", round(lcl_rr_dur_u,3), ", ", round(ucl_rr_dur_u, 3), ")")
rr_summary_standard_w <- paste0(round(rr_w,3), " (", round(lcl_rr_dur_w,3), ", ", round(ucl_rr_dur_w, 3), ")")

rr_summary_standard_etot   <- paste0(round(rr_u,3), " (", round(lcl_rr_etot_u,3), ", ", round(ucl_rr_etot_u, 3), ")")
rr_summary_standard_etot_w <- paste0(round(rr_w,3), " (", round(lcl_rr_etot_w,3), ", ", round(ucl_rr_etot_w, 3), ")")

rr_summary_boot_w <- paste0(round(rr_w,3), " (", round(lcl_rr_dur_w_boot,3), ", ", round(ucl_rr_dur_w_boot, 3), ")")
rr_summary_boot_etot_w <- paste0(round(rr_w,3), " (", round(lcl_rr_etot_w_boot,3), ", ", round(ucl_rr_etot_w_boot, 3), ")")

