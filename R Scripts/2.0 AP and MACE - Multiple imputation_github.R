# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 13/11/2024

# MULTIPLE IMPUTATION

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "data.table", "purrr", "ggplot2", "gtsummary", "gt", "mice", "finalfit", "tidylog"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Load data
load(file = "cv_events_cohort.rdata")

#Select variables for inclusion in imputation model
adjustment_set <- c("age_atprevcohortentry", "gender", "ethnicity_cat_cprdhes_4", "diag_prev", "patprac_2019imd_quintile", "prevcohortentry_year", "gpconsults_last6m", "smoking_status_cat", 
                    "prioralcoholabuse", "priorsubstanceabuse", "priordyslipidaemia", "priordiabetes", "priorhypertension", "priorrenaldisease", "priorliverdisease", "priorarrhythmia", "priorangina",
                    "apuse_prior2years", "lipiddrugs_prior2years", "hypertensiondrugs_prior2years", "antidiabetics_prior2years", "antidepressant_prior2years", "moodstab_prior2years", "anticoagulants_prior2years", "antiplatelets_prior2years", 
                    "zdrugs_prior2years", "anxiolytics_prior2years", "benzodiazepines_prior2years", "baseline_bmi_cat", "baselinedose_olanz", "priorhosp_physical_any", "priorhosp_psych_any", "smidiag_to_cohortentry")

aux_variables <- c("patid", "trt_group", "baseline_weightkg", "baseline_bmi", "perprotocol")

# Variables required for MACE outcomes
outcome_vars <- c("mace_itt_fuptime_cen", "mace_itt_event", "mace_pp_fuptime_cen", "mace_pp_event", "mace_pp_fupend_reason")

#Select variables for inclusion in imputation model
imputationvariables <- cv_events_cohort %>%
  select(patid, trt_group, all_of(c(adjustment_set, aux_variables, outcome_vars))) %>%
  mutate_if(is.character, as.factor) # convert character columns to factors

# Check structure of data
str(imputationvariables)

# Exclude patid from imputation
predM <- finalfit::missing_predictorMatrix(imputationvariables,
                                           drop_from_imputed = c("patid"),
                                           drop_from_imputer = c("patid")) 

#Check updated predictor matrix
matrix <- as.data.frame(predM)
remove(matrix)

# Impute data using mice, in parralel (futuremice)
cv_events_imputed <- futuremice(imputationvariables,
                                m = 25, 
                                maxit = 5, 
                                predictorMatrix = predM, 
                                n.core = 5, 
                                future.plan = "multisession", 
                                parallelseed = 123)

# View logged events
cv_events_imputed$loggedEvents

# View summary of imputation
summary(cv_events_imputed)

# Convert to long format for saving
cv_events_imputed_long <- mice::complete(cv_events_imputed, action = "long", include = TRUE)

# Add some extra variables so that they are available in the imputed dataset
cv_events_imputed_long <- cv_events_imputed_long %>%
  left_join(select(cv_events_cohort, patid, cohortentrydate, notadherent_date, notadherent_event, 
                   mace_itt_fupend_reason, age_atfirstdiag, monthsstononadherence,
                   stroke_itt_fuptime_cen, stroke_itt_event, stroke_pp_fuptime_cen, stroke_pp_event,
                   mi_itt_fuptime_cen, mi_itt_event, mi_pp_fuptime_cen, mi_pp_event,
                   cvdeath_itt_fuptime_cen, cvdeath_itt_event, cvdeath_pp_fuptime_cen, cvdeath_pp_event), by = "patid")

# Save to file
n_distinct(cv_events_imputed_long$patid) #20,404
save(cv_events_imputed_long, file = "cv_events_imputed_long.rdata")

remove(cv_events_imputed_long, cv_events_imputed, cv_events_cohort, imputationvariables, predM)

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
# packageVersion("finalfit") # 1.0.8
# packageVersion("tidylog") # 1.1.0
