# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 13/11/24

# BASELINE ANALYSIS

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "readr", "data.table", "rlang", "reshape2", "ggplot2", "gtsummary", "gt", 
                   "mice", "finalfit", "tidylog", "survminer", "survival", "ggsurvfit", "RColorBrewer", "ggforce", "grid", "gridExtra", "survey"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Set output directory

# Load data
load("cv_events_cohort.rdata")
load("cv_events_svy.rds")

# ANALYSIS ###

# Baseline characteristics
baseline_chars <- cv_events_cohort %>%
  select(trt_group, age_atprevcohortentry, gender, ethnicity_cat_cprdhes_4, diag_prev, age_atfirstdiag, smidiag_to_cohortentry, patprac_2019imd_quintile, prevcohortentry_year,
         prioralcoholabuse, priorangina, priorarrhythmia, priordiabetes, priordyslipidaemia, priorhypertension, priorliverdisease, priorrenaldisease, priorsubstanceabuse,
         apuse_prior2years, anticoagulants_prior2years, antidepressant_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years, antiplatelets_prior2years, anxiolytics_prior2years, benzodiazepines_prior2years, lipiddrugs_prior2years, moodstab_prior2years, zdrugs_prior2years,
         gpconsults_last6m, smoking_status_cat,
         baseline_bmi_cat, baseline_bmi, baseline_weightkg, 
         priorhosp_psych_any, priorhosp_physical_any,
         baselinedose_olanz) %>%
  mutate(baseline_bmi_cat = if_else(baseline_bmi_cat == "Underweight", "<18.5",
                                  if_else(baseline_bmi_cat == "Healthy", "≥18.5 to <25",
                                          if_else(baseline_bmi_cat == "Overweight", "≥25 to <30",
                                                  if_else(baseline_bmi_cat == "Obese", "≥30", baseline_bmi_cat)))))

# Specify continuous variable statistics for baseline tables
  statistic_list <- list(
    age_atprevcohortentry ~ "{median} ({p25}, {p75})",
    age_atfirstdiag ~ "{median} ({p25}, {p75})",
    baseline_bmi ~ "{mean} ({sd})",
    gpconsults_last6m ~ "{median} ({p25}, {p75})",
    baselinedose_olanz ~ "{mean} ({sd})",
    prevcohortentry_year ~ "{median} ({p25}, {p75})")

# Specify continuous variable digits for baseline tables
digits_list <- list(
  all_categorical() ~ c(0, 1),
  age_atprevcohortentry ~ c(0, 0),
  age_atfirstdiag ~ c(0, 0),
  baseline_bmi ~ c(1, 1),
  gpconsults_last6m ~ c(0, 0),
  baselinedose_olanz ~ c(1, 1),
  prevcohortentry_year ~ c(0, 0))

# Specify labels variable digits for baseline tables
label_list <- list(
  age_atprevcohortentry ~ "Age (years), median (IQR)",
  gender ~ "Sex, n (%)",
  ethnicity_cat_cprdhes_4 ~ "Ethnicity, n (%)",
  patprac_2019imd_quintile = "Index of Multiple Deprivation (quintile), n (%)",
  diag_prev ~ "SMI diagnosis, n (%)",
  smidiag_to_cohortentry ~ "Years from SMI diagnosis to index date, median (IQR)",
  priorhosp_psych_any ~ "Psychiatric hospitalisation in prior 2y, n (%)",
  priorhosp_physical_any ~ "Physical health hospitalisation in prior 2y, n (%)",
  age_atfirstdiag ~ "Age at SMI diagnosis, median (IQR)",
  priorsubstanceabuse ~ "Substance misuse", 
  prioralcoholabuse ~ "Alcohol misuse",
  priorhypertension ~ "Hypertension", 
  priordiabetes ~ "Diabetes", 
  priorangina ~ "Angina",
  priorarrhythmia ~ "Arrhythmia",
  priorliverdisease ~ "Liver disease", 
  priorrenaldisease ~ "Renal disease", 
  priordyslipidaemia ~ "Dyslipidaemia",
  apuse_prior2years ~ "Prior antipsychotic use", 
  antidepressant_prior2years ~ "Antidepressants", 
  moodstab_prior2years ~ "Mood stabilisers", 
  lipiddrugs_prior2years ~ "Lipid regulating medications", 
  hypertensiondrugs_prior2years ~ "Antihypertensives", 
  antidiabetics_prior2years ~ "Antidiabetics", 
  anticoagulants_prior2years ~ "Anticoagulants", 
  antiplatelets_prior2years ~ "Antiplatelets",
  anxiolytics_prior2years ~ "Anxiolytics",
  zdrugs_prior2years ~ "Z drugs",
  benzodiazepines_prior2years ~ "Benzodiazepines",
  gpconsults_last6m ~ "Primary care consultations in last 6m, median (IQR)",
  smoking_status_cat ~ "Smoking status, n (%)",
  baseline_bmi_cat ~ "BMI category, n (%)",
  baseline_bmi ~ "BMI (kg/m2), mean (SD)",
  baseline_weightkg ~ "Body weight (kg), median (IQR)",
  baselinedose_olanz ~ "Starting daily dose (Olanzapine equivalent) (mg), mean (SD)",
  prevcohortentry_year ~ "Index year, median (IQR)")

statistic_list <- list(
  age_atprevcohortentry ~ "{median} ({p25}, {p75})",
  age_atfirstdiag ~ "{median} ({p25}, {p75})",
  baseline_bmi ~ "{mean} ({sd})",
  gpconsults_last6m ~ "{median} ({p25}, {p75})",
  baselinedose_olanz ~ "{mean} ({sd})",
  prevcohortentry_year ~ "{median} ({p25}, {p75})")

# Footnote text
missing_pat_imd <- sum(cv_events_cohort$pat_2019imd_quintile %in% NA & !is.na(cv_events_cohort$prac_2019imd_quintile))

unweighted_baseline_chararacteristics_table <- as.data.frame(baseline_chars %>%
  tbl_summary(
    by = trt_group,
    statistic = statistic_list,
    digits = digits_list,
    label = label_list,
    sort = everything() ~ "alphanumeric")) %>%
  rename_all(~gsub("\\*", "", .)) %>%
  mutate(Characteristic = str_replace_all(Characteristic, "_", ""),
         row_num = row_number()) %>%
  add_row(row_num = 23.5, Characteristic = "Comorbidities, n (%)") %>%
  add_row(row_num = 32.5, Characteristic = "Concomitant medications, n (%)") %>%
  arrange(row_num) %>%
  mutate_at(vars(-row_num), ~ifelse(row_num == 23, gsub("2,", "2", .), .)) %>% # remove commas in calendar years
  select(-row_num) %>%
  gt() %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_labels(everything())) %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_body(columns = 1, rows = c(1, 2, 5, 11, 15, 16, 17, 23, 24, 34, 46, 47, 52, 58, 60, 62, 63, 64))) %>%
  tab_footnote("BMI, body mass index; mg, miligrams") %>%
  tab_footnote(paste0("Quintile of the 2019 English Index of Multiple Deprivation - defined according to the patient's postcode or, where this was missing (n=", missing_pat_imd, "), the primary care practice postcode."), locations = cells_body(columns = 1, rows = 17)) %>%
  tab_footnote("Comorbidities determined at point in the patient's medical history up to and including the index date.", locations = cells_body(columns = 1, rows = 24)) %>%
  tab_footnote("Concomitant medications defined according to prescriptions on, or in the two years prior to, the index date (prior antipsychotic use considered only prescriptions prior).", locations = cells_body(columns = 1, rows = 34)) %>%
  tab_footnote("Calculated according to the Defined Daily Dose method and expressed as an olanzapine equivalent dose.", locations = cells_body(columns = 1, rows = 64)) %>%
  
  sub_missing(missing_text = "") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = 2:5)

# Weighted characteristics following imputation
weighted_baseline_chararacteristics_table <- as.data.frame(cv_events_svy %>%
                                                  tbl_svysummary(include = c(trt_group, age_atprevcohortentry, gender, ethnicity_cat_cprdhes_4, diag_prev, age_atfirstdiag, smidiag_to_cohortentry, patprac_2019imd_quintile, prevcohortentry_year,
                                                                             prioralcoholabuse, priorangina, priorarrhythmia, priordiabetes, priordyslipidaemia, priorhypertension, priorliverdisease, priorrenaldisease, priorsubstanceabuse,
                                                                             apuse_prior2years, anticoagulants_prior2years, antidepressant_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years, antiplatelets_prior2years, anxiolytics_prior2years, benzodiazepines_prior2years, lipiddrugs_prior2years, moodstab_prior2years, zdrugs_prior2years,
                                                                             gpconsults_last6m, smoking_status_cat,
                                                                             baseline_bmi_cat, baseline_bmi, baseline_weightkg, 
                                                                             priorhosp_psych_any, priorhosp_physical_any,
                                                                             baselinedose_olanz),
                                                by = trt_group,
                                                statistic = statistic_list,
                                                digits = digits_list,
                                                label = label_list,
                                                sort = everything() ~ "alphanumeric")) %>%
  rename_all(~gsub("\\*", "", .)) %>%
  mutate(Characteristic = str_replace_all(Characteristic, "_", ""),
         row_num = row_number()) %>%
  add_row(row_num = 22.5, Characteristic = "Comorbidities, n (%)") %>%
  add_row(row_num = 31.5, Characteristic = "Concomitant medications, n (%)") %>%
  arrange(row_num) %>%
  mutate_at(vars(-row_num), ~ifelse(row_num == 22, gsub("2,", "2", .), .)) %>% # remove commas in calendar years
  select(-row_num) %>%
  gt() %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_labels(everything())) %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_body(columns = 1, rows = c(1, 2, 5, 10, 14, 15, 16, 22, 23, 33, 45, 46, 50, 55, 56, 57, 58, 59))) %>%
  tab_footnote("BMI, body mass index; mg, miligrams. Characteristics calculated based on values from one imputed dataset.") %>%
  tab_footnote(paste0("Quintile of the 2019 English Index of Multiple Deprivation - defined according to the patient's postcode or, where this was missing (unweighted n=", missing_pat_imd, "), the primary care practice postcode."), locations = cells_body(columns = 1, rows = 16)) %>%
  tab_footnote("Comorbidities determined at point in the patient's medical history up to and including the index date.", locations = cells_body(columns = 1, rows = 23)) %>%
  tab_footnote("Concomitant medications defined according to prescriptions on, or in the two years prior to, the index date (prior antipsychotic use considered only prescriptions prior).", locations = cells_body(columns = 1, rows = 33)) %>%
  tab_footnote("Calculated according to the Defined Daily Dose method and expressed as an olanzapine equivalent dose.", locations = cells_body(columns = 1, rows = 59)) %>%
  
  sub_missing(missing_text = "") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = 2:5)

# Discontinuation
discont_unweighted <- cv_events_cohort %>%
  select(trt_group, cohortentrydate, perprotocol, notadherent_date, notadherent_event) %>%
  mutate(monthsstononadherence = as.numeric(notadherent_date - cohortentrydate)/30.4,
         notadherent_event = str_to_title(notadherent_event)) %>%
  select(-cohortentrydate, -notadherent_date) %>%
  tbl_summary(by = trt_group,
              digits = list(
                perprotocol ~ c(0, 1),
                notadherent_event ~ c(0, 1)),
              label = list(
                perprotocol ~ "Remained adherent",
                notadherent_event ~ "Non-adherence reason",
                monthsstononadherence ~ "Time to non-adherence (months)")) %>%
  add_overall() %>%
  as_gt()

discont_weighted <- cv_events_svy %>%
  tbl_svysummary(
    include = c(
      trt_group,
      perprotocol,
      notadherent_event,
      monthsstononadherence
    ),
    by = trt_group,
    digits = list(perprotocol ~ c(0, 1), notadherent_event ~ c(0, 1)),
    label = list(
      perprotocol ~ "Remained adherent",
      notadherent_event ~ "Non-adherence reason",
      monthsstononadherence ~ "Time to non-adherence (months)"
    )
  ) %>%
  add_overall() %>%
  as_gt()

# Overall outcome incidence

# unweighted
followupstatus <- cv_events_cohort %>%
  select(mace_itt_event, mace_itt_fupend_reason) %>%
  tbl_summary()

mace_timetoevent_unweighted <- cv_events_cohort %>%
  select(mace_itt_event, mace_itt_fuptime_cen) %>%
  mutate(mace_itt_fuptime_cen_months = mace_itt_fuptime_cen / 30.25) %>%
  select(-mace_itt_fuptime_cen) %>%
  tbl_summary(by = mace_itt_event) %>%
  add_overall() %>%
  as_gt()

mace_timetoevent_weighted <- cv_events_svy %>%
  tbl_svysummary(include = c(mace_itt_event, mace_itt_fuptime_cen_months),
                 by = mace_itt_event) %>%
  add_overall() %>%
  as_gt()

# Create a survival object using time in years
surv_obj <- Surv(time = cv_events_cohort$mace_itt_fuptime_cen / 365, 
                 event = cv_events_cohort$mace_itt_event)

# Fit a Kaplan-Meier survival curve
km_fit <- survfit(surv_obj ~ 1)

# Create a Kaplan-Meier plot
km_plot <- survminer::ggsurvplot(km_fit, data = cv_events_cohort, conf.int = TRUE, pval = TRUE, xlim = c(0, 5),
                                 xlab = "Time (years)", ylab = "Survival Probability",
                                 main = "Kaplan-Meier Survival Curve")

surv_probs <- tbl_survfit(km_fit, times = c(1,2,3,4,5),
                          estimate_fun = function(x) style_percent(x, symbol = TRUE, digits = 1, columns = "survival"), 
                          reverse = TRUE)
surv_probs %>%
  as_gt()

# Completeness of baseline covariates
bl_covariates <- cv_events_cohort %>%
  select(trt_group,
         age_atprevcohortentry, ethnicity_cat_cprdhes_4, patprac_2019imd_quintile, gender,
         baselinedose_olanz, priorhosp_psych_any, apuse_prior2years, diag_prev, smidiag_to_cohortentry,
         prioralcoholabuse, priorangina, priorarrhythmia, priordiabetes, priordyslipidaemia, priorhypertension, priorliverdisease, priorrenaldisease, priorsubstanceabuse,
         anticoagulants_prior2years, antidepressant_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years, antiplatelets_prior2years, anxiolytics_prior2years,
         benzodiazepines_prior2years, lipiddrugs_prior2years, moodstab_prior2years, zdrugs_prior2years,
         baseline_bmi_cat, prevcohortentry_year, gpconsults_last6m, priorhosp_physical_any, smoking_status_cat)
  
# Convert columns to binary indicators and calculate completeness
bl_covariates_binary <- bl_covariates
bl_covariates_binary[, 2:ncol(bl_covariates)] <- lapply(bl_covariates[, 2:ncol(bl_covariates)], function(x) ifelse(!is.na(x), 1, 0))

# Define variable labels
var_labels <- c(
  "age_atprevcohortentry" = "Age", "ethnicity_cat_cprdhes_4" = "Ethnicity",
  "patprac_2019imd_quintile" = "Index of Multiple Deprivation", "gender" = "Sex",
  "baselinedose_olanz" = "Starting daily dose", "priorhosp_psych_any" = "Prior psychiatric hospitalisation",
  "priorhosp_physical_any" = "Prior physical health hospitalisation", "apuse_prior2years" = "Prior antipsychotic use", 
  "diag_prev" = "SMI diagnosis", "smidiag_to_cohortentry" = "Years from SMI diagnosis to index date",
  "prioralcoholabuse" = "Alcohol misuse", "priorhypertension" = "Hypertension", 
  "priordiabetes" = "Diabetes", "priorangina" = "Angina", "priorarrhythmia" = "Arrhythmia",
  "priorliverdisease" = "Liver disease", "priorrenaldisease" = "Renal disease", 
  "priorsubstanceabuse" = "Substance misuse", "priordyslipidaemia" = "Dyslipidaemia",
  "antidepressant_prior2years" = "Antidepressant prescription", 
  "moodstab_prior2years" = "Mood stabiliser prescription", 
  "lipiddrugs_prior2years" = "Lipid regulating medication prescription", 
  "hypertensiondrugs_prior2years" = "Antihypertensive prescription", 
  "antidiabetics_prior2years" = "Antidiabetic prescription", 
  "anticoagulants_prior2years" = "Anticoagulant prescription", 
  "antiplatelets_prior2years" = "Antiplatelet prescription",
  "anxiolytics_prior2years" = "Anxiolytic prescription",
  "zdrugs_prior2years" = "Z-drug prescription",
  "benzodiazepines_prior2years" = "Benzodiazepine prescription",
  "gpconsults_last6m" = "Primary care consultations in last 6m",
  "smoking_status_cat" = "Smoking status",
  "baseline_bmi_cat" = "BMI category",
  "prevcohortentry_year" = "Index year")

# Create and format the completeness table
covariate_completeness <- data.frame(
  Variable = names(bl_covariates)[2:ncol(bl_covariates)],
  Completeness = sapply(bl_covariates[, 2:ncol(bl_covariates)], function(x) {
    round(100 * sum(!is.na(x)) / length(x), 2)})) %>%
  mutate(Variable = ifelse(Variable %in% names(var_labels), var_labels[Variable], Variable)) %>% # Map variable names to labels
  arrange(desc(Completeness), Variable) %>% # Sort by completeness (descending) and then by variable name
  gt() %>%
  fmt_number(columns = "Completeness", decimals = 2) %>%
  cols_label(
    Variable = "Variable",
    Completeness = "Completeness (%)") %>%
  tab_header(title = "Data Completeness Summary")

# Package versions
# packageVersion("dplyr") # 1.1.4
# packageVersion("tidyr") # 1.3.1
# packageVersion("stringr") # 1.5.1
# packageVersion("forcats") # 1.0.0
# packageVersion("lubridate") # 1.9.4
# packageVersion("readr") # 2.1.5
# packageVersion("data.table") # 1.16.4
# packageVersion("rlang") # 1.1.4
# packageVersion("reshape2") # 1.4.4
# packageVersion("ggplot2") # 3.5.1
# packageVersion("gtsummary") # 2.0.4
# packageVersion("gt") # 0.11.1
# packageVersion("mice") # 3.17.0
# packageVersion("finalfit") # 1.0.8
# packageVersion("tidylog") # 1.1.0
# packageVersion("survminer") # 0.5.0
# packageVersion("survival") # 3.7.0
# packageVersion("ggsurvfit") # 1.1.0
# packageVersion("RColorBrewer") # 1.1.3
# packageVersion("ggforce") # 0.4.2
# packageVersion("grid") # 4.4.2
# packageVersion("gridExtra") # 2.3
# packageVersion("survey") # 4.4.2