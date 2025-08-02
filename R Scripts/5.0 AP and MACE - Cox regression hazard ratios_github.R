# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 13/11/2024

# Cox regression - Hazard ratios

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "data.table", "purrr", "ggplot2", "gtsummary", "gt", 
                   "mice", "tidylog", "survminer", "survival", "janitor", "WeightIt", "MatchThem", "sandwich", "nnet"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Set output directory

# Load data
load("cv_events_imputed_long.Rdata")

# Convert to mids
cv_events_imputed <- as.mids(cv_events_imputed_long)
remove(cv_events_imputed_long)

# OVERLAP WEIGHTING ###

# calculate overlap weights
imputed_weighted <- MatchThem::weightthem(trt_group ~ 
                                            age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                                            ethnicity_cat_cprdhes_4 + patprac_2019imd_quintile + gender +
                                            baselinedose_olanz  + I(baselinedose_olanz^2) + I(baselinedose_olanz^3) +
                                            priorhosp_psych_any + apuse_prior2years + diag_prev + 
                                            smidiag_to_cohortentry + I(smidiag_to_cohortentry^2) + 
                                            prioralcoholabuse + priorangina + priorarrhythmia + priordiabetes + priordyslipidaemia + 
                                            priorhypertension + priorliverdisease + priorrenaldisease + priorsubstanceabuse +
                                            anticoagulants_prior2years + antidepressant_prior2years + antidiabetics_prior2years + hypertensiondrugs_prior2years + antiplatelets_prior2years +
                                            anxiolytics_prior2years + benzodiazepines_prior2years + lipiddrugs_prior2years + moodstab_prior2years + zdrugs_prior2years +
                                            baseline_bmi_cat + prevcohortentry_year + I(prevcohortentry_year^2) +
                                            gpconsults_last6m +  I(gpconsults_last6m^2) + 
                                            priorhosp_physical_any + smoking_status_cat +
                                            
                                            # Interactions
                                            baselinedose_olanz*diag_prev + I(baselinedose_olanz^2)*diag_prev + I(baselinedose_olanz^3)*diag_prev,
                                          
                                          datasets = cv_events_imputed, 
                                          approach = "within",
                                          
                                          method = "glm", 
                                          include.obj = TRUE,
                                          estimand = "ATO") 

# Cox model, overlap weighted, missing covariates imputed
cox_model_table <- function(model) {
  table <- as.data.frame(tbl_regression(model, exp = FALSE,
                                        pvalue_fun = function(x) style_pvalue(x, digits = 3),
                                        tidy_fun = purrr::partial(tidy_robust, vcov = sandwich::vcovHC, vcov_args = list(type = "HC2")))) %>%
    janitor::clean_names() %>%
    separate(x95_percent_ci, into = c("ci_high", "ci_low"), 
             sep = ", ") %>%
    mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
    mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
    dplyr::select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
    filter(row_number() > 2) %>%
    mutate(result = paste0(format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
                           " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
           result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
           label = case_when(group == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                             group == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                             group == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                             TRUE ~ "Unknown")) %>%
    dplyr::select(label, result, p_value) %>%
    gt()
  
  return(table)}

# MACE
mace_itt_cox_model <- with(imputed_weighted, coxph(Surv(mace_itt_fuptime_cen, mace_itt_event) ~ trt_group))
mace_itt_cox_tab <- cox_model_table(mace_itt_cox_model)

mace_pp_cox_model <- with(imputed_weighted, coxph(Surv(mace_pp_fuptime_cen, mace_pp_event) ~ trt_group))
mace_pp_cox_tab <- cox_model_table(mace_pp_cox_model)

# Cox model - overlap weighted, complete case
load("cv_events_cohort.Rdata")

# Remove people with any missing covariates
covariate_names <- c(
  "age_atprevcohortentry", "ethnicity_cat_cprdhes_4", "patprac_2019imd_quintile", "gender","baselinedose_olanz", "priorhosp_psych_any", 
  "apuse_prior2years", "diag_prev","smidiag_to_cohortentry", "prioralcoholabuse", "priorangina", "priorarrhythmia", "priordiabetes", "priordyslipidaemia",
  "priorhypertension", "priorliverdisease", "priorrenaldisease", "priorsubstanceabuse","anticoagulants_prior2years", "antidepressant_prior2years", "antidiabetics_prior2years", "hypertensiondrugs_prior2years", "antiplatelets_prior2years",
  "anxiolytics_prior2years", "benzodiazepines_prior2years", "lipiddrugs_prior2years", "moodstab_prior2years", "zdrugs_prior2years", "baseline_bmi_cat", "prevcohortentry_year", "gpconsults_last6m", "priorhosp_physical_any", "smoking_status_cat")

cv_events_cohort <- cv_events_cohort %>%
  mutate(missing_covariates = if_any(all_of(covariate_names), is.na)) %>%
  filter(missing_covariates == FALSE) %>%
  dplyr::select(-missing_covariates)

# Calculate Overlap weights
ids <- cv_events_cohort$patid
trt <- cv_events_cohort$trt_group

multinom_model_manual <- multinom(trt_group ~
                                    age_atprevcohortentry + I(age_atprevcohortentry^2) +
                                    ethnicity_cat_cprdhes_4 + patprac_2019imd_quintile + gender +
                                    baselinedose_olanz  + I(baselinedose_olanz^2) + I(baselinedose_olanz^3) +
                                    priorhosp_psych_any + apuse_prior2years + diag_prev +
                                    smidiag_to_cohortentry + I(smidiag_to_cohortentry^2) +
                                    prioralcoholabuse + priorangina + priorarrhythmia + priordiabetes + priordyslipidaemia +
                                    priorhypertension + priorliverdisease + priorrenaldisease + priorsubstanceabuse +
                                    anticoagulants_prior2years + antidepressant_prior2years + antidiabetics_prior2years + hypertensiondrugs_prior2years + antiplatelets_prior2years +
                                    anxiolytics_prior2years + benzodiazepines_prior2years + lipiddrugs_prior2years + moodstab_prior2years + zdrugs_prior2years +
                                    baseline_bmi_cat + prevcohortentry_year + I(prevcohortentry_year^2) +
                                    gpconsults_last6m +  I(gpconsults_last6m^2) +
                                    priorhosp_physical_any + smoking_status_cat +
                                    # Interactions
                                    baselinedose_olanz*diag_prev + I(baselinedose_olanz^2)*diag_prev + I(baselinedose_olanz^3)*diag_prev,
                                  data = cv_events_cohort,
                                  trace = FALSE)

gps_matrix <- predict(multinom_model_manual, type = "probs")

# Convert the data to a data.table object
manual <- data.table(patid = ids, 
                     trt_group = trt, 
                     ps_a = gps_matrix[,1],
                     ps_o = gps_matrix[,2],
                     ps_q = gps_matrix[,3],
                     ps_r = gps_matrix[,4])

# Add the inverse probability columns
manual[, `:=`(ip_a = 1/ps_a,
              ip_o = 1/ps_o,
              ip_q = 1/ps_q,
              ip_r = 1/ps_r)]

# Calculate the sum of inverse probabilities for each row
manual[, ip_sum := ip_a + ip_o + ip_q + ip_r]

# Calculate the overlap weights based on treatment group
manual[, ow_manual := fifelse(trt_group == "Aripiprazole", ip_a / ip_sum,
                              fifelse(trt_group == "Olanzapine", ip_o / ip_sum,
                                      fifelse(trt_group == "Quetiapine", ip_q / ip_sum,
                                              fifelse(trt_group == "Risperidone", ip_r / ip_sum, NA_real_))))]

# Calculate the sum of original weights for each treatment group
ari_sum <- manual[trt_group == "Aripiprazole", sum(ow_manual)]
olanz_sum <- manual[trt_group == "Olanzapine", sum(ow_manual)]
quet_sum <- manual[trt_group == "Quetiapine", sum(ow_manual)]
risp_sum <- manual[trt_group == "Risperidone", sum(ow_manual)]

# Normalize the weights based on the treatment group
manual[, manual_weights_norm := fifelse(trt_group == "Aripiprazole", ow_manual / ari_sum,
                                        fifelse(trt_group == "Olanzapine", ow_manual / olanz_sum,
                                                fifelse(trt_group == "Quetiapine", ow_manual / quet_sum,
                                                        fifelse(trt_group == "Risperidone", ow_manual / risp_sum, NA_real_))))]

cv_events_cohort$w1 <- manual$manual_weights_norm

rm(manual, gps_matrix)

# Cox models

cox_model_table <- function(model) {
  table <- as.data.frame(tbl_regression(model, exp = FALSE,
                                        pvalue_fun = function(x) style_pvalue(x, digits = 3))) %>%
    janitor::clean_names() %>%
    separate(x95_percent_ci, into = c("ci_high", "ci_low"), 
             sep = ", ") %>%
    mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
    mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
    dplyr::select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
    filter(row_number() > 2) %>%
    mutate(result = paste0(format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
                           " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
           result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
           label = case_when(group == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                             group == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                             group == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                             TRUE ~ "Unknown")) %>%
    dplyr::select(label, result, p_value) %>%
    gt()
  
  return(table)}

output_dir <- 

# MACE 5y
mace_itt_cox_model_cc <- coxph(Surv(mace_itt_fuptime_cen, mace_itt_event) ~ trt_group, weights = w1, data = cv_events_cohort)
mace_itt_cox_model_cc_tab <- cox_model_table(mace_itt_cox_model_cc)

mace_pp_cox_model_cc <- coxph(Surv(mace_pp_fuptime_cen, mace_pp_event) ~ trt_group, weights = w1, data = cv_events_cohort)
mace_pp_cox_model_cc_tab <- cox_model_table(mace_pp_cox_model_cc)

# COVARIATE ADJUSTMENT 

# Cox model - with covariate adjustment instead of weighting
library(dplyr)
library(survival)

load("cv_events_imputed_long.Rdata")

cv_events_imputed_long <- cv_events_imputed_long %>%
  as.mids()

cox_itt_adjusted <- with(cv_events_imputed_long, 
                     coxph(Surv(mace_itt_fuptime_cen, mace_itt_event) ~ trt_group +
                             age_atprevcohortentry + I(age_atprevcohortentry^2) +
                             ethnicity_cat_cprdhes_4 + patprac_2019imd_quintile + gender +
                             baselinedose_olanz  + I(baselinedose_olanz^2) + I(baselinedose_olanz^3) +
                             priorhosp_psych_any + apuse_prior2years + diag_prev +
                             smidiag_to_cohortentry + I(smidiag_to_cohortentry^2) +
                             prioralcoholabuse + priorangina + priorarrhythmia + priordiabetes + priordyslipidaemia +
                             priorhypertension + priorliverdisease + priorrenaldisease + priorsubstanceabuse +
                             anticoagulants_prior2years + antidepressant_prior2years + antidiabetics_prior2years + hypertensiondrugs_prior2years + antiplatelets_prior2years +
                             anxiolytics_prior2years + benzodiazepines_prior2years + lipiddrugs_prior2years + moodstab_prior2years + zdrugs_prior2years +
                             baseline_bmi_cat + prevcohortentry_year + I(prevcohortentry_year^2) +
                             gpconsults_last6m +  I(gpconsults_last6m^2) +
                             priorhosp_physical_any + smoking_status_cat +
                             # Interactions
                             baselinedose_olanz*diag_prev + I(baselinedose_olanz^2)*diag_prev + I(baselinedose_olanz^3)*diag_prev))

itt_cox_adjusted_tab <- as.data.frame(tbl_regression(cox_itt_adjusted,
               estimate_fun = purrr::partial(style_ratio, digits = 8),
               pvalue_fun = purrr::partial(style_sigfig, digits = 3))) %>%
  clean_names() %>%
  filter(row_number() > 2) %>%
  filter(row_number() < 4) %>%
  separate(x95_percent_ci, into = c("ci_high", "ci_low"), 
           sep = ", ") %>%
  mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
  mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
  select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
  mutate(result = paste0(format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
                         " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
         result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
         label = case_when(grepl("Olanzapine", group) ~ "Aripiprazole vs. Olanzapine",
                           grepl("Quetiapine", group) ~ "Aripiprazole vs. Quetiapine",
                           grepl("Risperidone", group) ~ "Aripiprazole vs. Risperidone",
                           TRUE ~ "Unknown")) %>%
  select(label, result)

cox_pp_adjusted <- with(cv_events_imputed_long, 
                         coxph(Surv(mace_pp_fuptime_cen, mace_pp_event) ~ trt_group +
                                 age_atprevcohortentry + I(age_atprevcohortentry^2) +
                                 ethnicity_cat_cprdhes_4 + patprac_2019imd_quintile + gender +
                                 baselinedose_olanz  + I(baselinedose_olanz^2) + I(baselinedose_olanz^3) +
                                 priorhosp_psych_any + apuse_prior2years + diag_prev +
                                 smidiag_to_cohortentry + I(smidiag_to_cohortentry^2) +
                                 prioralcoholabuse + priorangina + priorarrhythmia + priordiabetes + priordyslipidaemia +
                                 priorhypertension + priorliverdisease + priorrenaldisease + priorsubstanceabuse +
                                 anticoagulants_prior2years + antidepressant_prior2years + antidiabetics_prior2years + hypertensiondrugs_prior2years + antiplatelets_prior2years +
                                 anxiolytics_prior2years + benzodiazepines_prior2years + lipiddrugs_prior2years + moodstab_prior2years + zdrugs_prior2years +
                                 baseline_bmi_cat + prevcohortentry_year + I(prevcohortentry_year^2) +
                                 gpconsults_last6m +  I(gpconsults_last6m^2) +
                                 priorhosp_physical_any + smoking_status_cat +
                                 # Interactions
                                 baselinedose_olanz*diag_prev + I(baselinedose_olanz^2)*diag_prev + I(baselinedose_olanz^3)*diag_prev))

pp_cox_adjusted_tab <- as.data.frame(tbl_regression(cox_pp_adjusted,
                                      estimate_fun = purrr::partial(style_ratio, digits = 8),
                                      pvalue_fun = purrr::partial(style_sigfig, digits = 3))) %>%
  clean_names() %>%
  filter(row_number() > 2) %>%
  filter(row_number() < 4) %>%
  separate(x95_percent_ci, into = c("ci_high", "ci_low"), 
           sep = ", ") %>%
  mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
  mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
  select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
  mutate(result = paste0(format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
                         " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
         result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
         label = case_when(grepl("Olanzapine", group) ~ "Aripiprazole vs. Olanzapine",
                           grepl("Quetiapine", group) ~ "Aripiprazole vs. Quetiapine",
                           grepl("Risperidone", group) ~ "Aripiprazole vs. Risperidone",
                           TRUE ~ "Unknown")) %>%
  select(label, result)

cox_covariate_adjusted_tab <- itt_cox_adjusted_tab %>%
  mutate(Type = "ITT") %>%
  bind_rows(pp_cox_adjusted_tab) %>%
  mutate(Type = case_when(is.na(Type) ~ "PP",
                          TRUE ~ Type)) %>%
  select(Type, Comparison = label, Estimate = result) %>%
  gt()


# IP WEIGHTING ####

# calculate IP weights
imputed_ipw <- MatchThem::weightthem(trt_group ~ 
                                       age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                                       ethnicity_cat_cprdhes_4 + patprac_2019imd_quintile + gender +
                                       baselinedose_olanz  + I(baselinedose_olanz^2) + I(baselinedose_olanz^3) +
                                       priorhosp_psych_any + apuse_prior2years + diag_prev + 
                                       smidiag_to_cohortentry + I(smidiag_to_cohortentry^2) + 
                                       prioralcoholabuse + priorangina + priorarrhythmia + priordiabetes + priordyslipidaemia + 
                                       priorhypertension + priorliverdisease + priorrenaldisease + priorsubstanceabuse +
                                       anticoagulants_prior2years + antidepressant_prior2years + antidiabetics_prior2years + hypertensiondrugs_prior2years + antiplatelets_prior2years +
                                       anxiolytics_prior2years + benzodiazepines_prior2years + lipiddrugs_prior2years + moodstab_prior2years + zdrugs_prior2years +
                                       baseline_bmi_cat + prevcohortentry_year + I(prevcohortentry_year^2) +
                                       gpconsults_last6m +  I(gpconsults_last6m^2) + 
                                       priorhosp_physical_any + smoking_status_cat +
                                       
                                       # Interactions
                                       baselinedose_olanz*diag_prev + I(baselinedose_olanz^2)*diag_prev + I(baselinedose_olanz^3)*diag_prev,
                                     
                                     datasets = cv_events_imputed, 
                                     approach = "within", 
                                     
                                     method = "glm", 
                                     include.obj = TRUE,
                                     estimand = "ATE") 

summary(imputed_ipw)

imputed_ipw_trimmed <- trim(imputed_ipw, at = 0.995, lower = TRUE)

summary(imputed_ipw_trimmed)

# Cox model
cox_model_table <- function(model) {
  table <- as.data.frame(tbl_regression(model, exp = FALSE,
                                        pvalue_fun = function(x) style_pvalue(x, digits = 3),
                                        tidy_fun = purrr::partial(tidy_robust, vcov = sandwich::vcovHC, vcov_args = list(type = "HC2")))) %>%
    janitor::clean_names() %>%
    separate(x95_percent_ci, into = c("ci_high", "ci_low"), 
             sep = ", ") %>%
    mutate_at(vars(2:5), as.numeric) %>% # convert to numeric
    mutate_at(vars(2:4), ~round(exp(. * -1), 2)) %>% # inverse the log hazards, exponentiate it, round to 2dp
    dplyr::select(group = characteristic, hr = log_hr, ci_low, ci_high, p_value) %>%
    filter(row_number() > 2) %>%
    mutate(result = paste0(format(hr, nsmall = 2), # nsmall ensures 2 digits are used in all cells (otherwise 1.00 might show as just 1)
                           " (", paste0(format(ci_low, nsmall = 2), ", ", format(ci_high, nsmall = 2), ")")),
           result = gsub("(\\()\\s+|(\\[)\\s+|\\s+(\\))|\\s+(\\])", "\\1\\2\\3\\4", result),
           label = case_when(group == "Olanzapine" ~ "Aripiprazole vs. Olanzapine",
                             group == "Quetiapine" ~ "Aripiprazole vs. Quetiapine",
                             group == "Risperidone" ~ "Aripiprazole vs. Risperidone",
                             TRUE ~ "Unknown")) %>%
    dplyr::select(label, result, p_value) %>%
    gt()
  
  return(table)}

# MACE
mace_itt_cox_ipw_trim_model <- with(imputed_ipw_trimmed, coxph(Surv(mace_itt_fuptime_cen, mace_itt_event) ~ trt_group))
mace_itt_cox_ipw_trim_tab <- cox_model_table(mace_itt_cox_ipw_trim_model)

mace_pp_cox_ipw_trim_model <- with(imputed_ipw_trimmed, coxph(Surv(mace_pp_fuptime_cen, mace_pp_event) ~ trt_group))
mace_pp_cox_ipw_trim_tab <- cox_model_table(mace_pp_cox_ipw_trim_model)

# Package versions
# packageVersion("dplyr") # 1.1.4
# packageVersion("tidyr") # 1.3.1
# packageVersion("stringr") # 1.5.1
# packageVersion("forcats") # 1.0.0
# packageVersion("lubridate") # 1.9.4
# packageVersion("data.table") # 1.16.4
# packageVersion("purrr") # 1.0.2
# packageVersion("ggplot2") # 3.5.1
# packageVersion("gtsummary") # 2.0.4
# packageVersion("gt") # 0.11.1
# packageVersion("mice") # 3.17.0
# packageVersion("tidylog") # 1.1.0
# packageVersion("survminer") # 0.5.0
# packageVersion("survival") # 3.7.0
# packageVersion("janitor") # 2.2.0
# packageVersion("WeightIt") # 1.3.2
# packageVersion("MatchThem") # 1.2.1
# packageVersion("sandwich") # 3.1.1
# packageVersion("nnet") # 7.3.19
