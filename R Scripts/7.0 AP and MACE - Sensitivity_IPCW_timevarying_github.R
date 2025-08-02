# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 16/12/2024

# PP - overlap weights, risk plots, bootstrapping (with IPCW - time-varying cumulative comorbidity count)

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("tidyr", "stringr", "lubridate", "data.table", "reshape2", "ggplot2", "gtsummary", "gt", 
                   "mice", "nnet", "dplyr", "scales", "future", "furrr", "vroom"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Set output directory

# Load data
load("cv_events_imputed_long.Rdata")
load("cv_events_cohort.Rdata")

cv_events_imputed_long <- cv_events_imputed_long %>%
  filter(.imp > 0) %>%
  dplyr::select(.imp, patid, trt_group,
                age_atprevcohortentry, ethnicity_cat_cprdhes_4, patprac_2019imd_quintile, gender, baselinedose_olanz, priorhosp_psych_any, apuse_prior2years, diag_prev, smidiag_to_cohortentry,
                prioralcoholabuse, priorangina, priorarrhythmia, priordiabetes, priordyslipidaemia, priorhypertension, priorliverdisease, priorrenaldisease, priorsubstanceabuse,
                anticoagulants_prior2years, antidepressant_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years, antiplatelets_prior2years,
                anxiolytics_prior2years, benzodiazepines_prior2years, lipiddrugs_prior2years, moodstab_prior2years, zdrugs_prior2years, baseline_bmi_cat, prevcohortentry_year, 
                gpconsults_last6m, priorhosp_physical_any, smoking_status_cat,
                mace_pp_event, mace_pp_fuptime_cen,
                stroke_pp_event, stroke_pp_fuptime_cen,
                mi_pp_event, mi_pp_fuptime_cen,
                cvdeath_pp_event, cvdeath_pp_fuptime_cen) %>%
  left_join(dplyr::select(cv_events_cohort, patid, mace_pp_fupend_reason), by = "patid") %>%
  mutate(censor = case_when(mace_pp_fupend_reason == "End of follow-up" ~ 0,
                            mace_pp_fupend_reason == "Had MACE" ~ 0,
                            TRUE ~ 1),
         mace_pp_fuptime_cen = pmax(round(mace_pp_fuptime_cen / 30.25, 0), 1)) # convert from days to months, convert 0 values to 1 (as the minimum ammount of f-up)

# Create time varying diagnoses
# Specify col types across HES files
hes_config <- list(
  patient = list(
    patid = col_character()),
  diagnosis_epi = list(
    spno = col_number(), epikey = col_number(), ICDx = col_skip(), epistart = col_date(format = "%d/%m/%Y"), 
    epiend = col_date(format = "%d/%m/%Y"), ICD = col_character(), d_order = col_number()),
  hosp = list(
    ICDx = col_skip(), admidate = col_date(format = "%d/%m/%Y"), discharged = col_date(format = "%d/%m/%Y")),
  episodes = list(
    spno = col_number(), admidate = col_date(format = "%d/%m/%Y"), epistart = col_date(format = "%d/%m/%Y"), 
    epiend = col_date(format = "%d/%m/%Y"), discharged = col_date(format = "%d/%m/%Y")),
  deaths = list(
    dod = col_date(format = "%d/%m/%Y")))

# Function to load HES data with dynamic column type handling
load_hes <- function(file_path, config) {
  # Determine suffix based on file_path contents
  suffix <- ifelse(grepl("gold", file_path, ignore.case = TRUE), "-G",
                   ifelse(grepl("aurum", file_path, ignore.case = TRUE), "-A",
                          stop("Database type (Gold or Aurum) not found in file_path.")))
  
  # Import data using vroom with specified column types from config
  data <- vroom(file.path(path, file_path), delim = "\t", escape_double = FALSE, 
                col_types = cols(
                  !!!config,  # Unquote-splice to directly use hes_config
                  patid = col_character()),  # Always include patid column type
                trim_ws = TRUE)
  
  # Update 'patid' based on database suffix
  data <- data %>% mutate(patid = paste0(patid, suffix))
  
  return(data)}

hes_diagnosis_epi_gold_1 <- load_hes("hes_diagnosis_epi_21_000729.txt", hes_config$diagnosis_epi)
hes_diagnosis_epi_gold_2 <- load_hes("hes_diagnosis_epi_21_000729_request2.txt", hes_config$diagnosis_epi)
hes_diagnosis_epi_aurum <- load_hes("hes_diagnosis_epi_21_000729.txt", hes_config$diagnosis_epi)

hes_primarydiag <- rbindlist(list(hes_diagnosis_epi_gold_1, hes_diagnosis_epi_gold_2, hes_diagnosis_epi_aurum)) %>%
  filter(epistart < '2020-01-01')

remove(hes_diagnosis_epi_gold_1, hes_diagnosis_epi_gold_2, hes_diagnosis_epi_aurum, hes_config, load_hes)
gc()

# Arrhythmia

arrhythmia_at_bl <- cv_events_cohort %>%
  filter(priorarrhythmia == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "arrhythmia", diff_months = 0)

arrhythmia_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priorarrhythmia) %>%
  filter(priorarrhythmia == "No")

load("SMIArrhythmia_C.Rdata")

arrhythmia_post <- arrhythmia_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIArrhythmia_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "arrhythmia") %>%
  select(patid, cohortentrydate, eventdate, condition)

arrhythmia_post_hes <- arrhythmia_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("I47", "I48", "I49"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate, epistart <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "arrhythmia") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

arrhythmia_post_combined <- arrhythmia_post %>%
  bind_rows(arrhythmia_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_arrhythmia <- arrhythmia_not %>%
  anti_join(arrhythmia_post_combined, by = "patid") %>%
  mutate(condition = "arrhythmia") %>%
  mutate(diff_months = 999)

arrhythmia_post_combined <- arrhythmia_post_combined %>%
  bind_rows(no_incident_arrhythmia) %>%
  bind_rows(arrhythmia_at_bl)

rm(arrhythmia_post, arrhythmia_post_hes, no_incident_arrhythmia, arrhythmia_at_bl, arrhythmia_not, SMIArrhythmia_C)

# Renal disease

renaldisease_at_bl <- cv_events_cohort %>%
  filter(priorrenaldisease == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "renaldisease", diff_months = 0)

renaldisease_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priorrenaldisease) %>%
  filter(priorrenaldisease == "No")

load("SMIRenalDisease_C.Rdata")

renaldisease_post <- renaldisease_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIRenalDisease_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "renaldisease") %>%
  select(patid, cohortentrydate, eventdate, condition)

renaldisease_post_hes <- renaldisease_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("K7"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate, epistart <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "renaldisease") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

renaldisease_post_combined <- renaldisease_post %>%
  bind_rows(renaldisease_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_renaldisease <- renaldisease_not %>%
  anti_join(renaldisease_post_combined, by = "patid") %>%
  mutate(condition = "renaldisease") %>%
  mutate(diff_months = 999)

renaldisease_post_combined <- renaldisease_post_combined %>%
  bind_rows(no_incident_renaldisease) %>%
  bind_rows(renaldisease_at_bl)

rm(renaldisease_post, renaldisease_post_hes, no_incident_renaldisease, renaldisease_at_bl, renaldisease_not, SMIRenalDisease_C)

# Liver
liverdisease_at_bl <- cv_events_cohort %>%
  filter(priorliverdisease == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "liverdisease", diff_months = 0)

liverdisease_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priorliverdisease) %>%
  filter(priorliverdisease == "No")

load("SMILiverDisease_C.Rdata")

liverdisease_post <- liverdisease_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMILiverDisease_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "liverdisease") %>%
  select(patid, cohortentrydate, eventdate, condition)

liverdisease_post_hes <- liverdisease_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("N17", "N18", "N19"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate, epistart <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "liverdisease") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

liverdisease_post_combined <- liverdisease_post %>%
  bind_rows(liverdisease_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_liverdisease <- liverdisease_not %>%
  anti_join(liverdisease_post_combined, by = "patid") %>%
  mutate(condition = "liverdisease") %>%
  mutate(diff_months = 999)

liverdisease_post_combined <- liverdisease_post_combined %>%
  bind_rows(no_incident_liverdisease) %>%
  bind_rows(liverdisease_at_bl)

rm(liverdisease_post, liverdisease_post_hes, no_incident_liverdisease, liverdisease_at_bl, liverdisease_not, SMILiverDisease_C)

# Angina
angina_at_bl <- cv_events_cohort %>%
  filter(priorangina == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "angina", diff_months = 0)

angina_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priorangina) %>%
  filter(priorangina == "No")

load("SMIAngina_C.Rdata")

angina_post <- angina_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIAngina_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "angina") %>%
  select(patid, cohortentrydate, eventdate, condition)

angina_post_hes <- angina_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("I20"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate, epistart <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "angina") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

angina_post_combined <- angina_post %>%
  bind_rows(angina_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_angina <- angina_not %>%
  anti_join(angina_post_combined, by = "patid") %>%
  mutate(condition = "angina") %>%
  mutate(diff_months = 999)

angina_post_combined <- angina_post_combined %>%
  bind_rows(no_incident_angina) %>%
  bind_rows(angina_at_bl)

rm(angina_post, angina_post_hes, no_incident_angina, angina_at_bl, angina_not, SMIAngina_C)

# Hypertension
hypertension_at_bl <- cv_events_cohort %>%
  filter(priorhypertension == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "hypertension", diff_months = 0)

hypertension_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priorhypertension) %>%
  filter(priorhypertension == "No")

load("SMIHypertension_C.Rdata")

hypertension_post <- hypertension_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIHypertension_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "hypertension") %>%
  select(patid, cohortentrydate, eventdate, condition)

hypertension_post_hes <- hypertension_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("I10", "I11", "I12", "I13", "I14", "I15", "K76.6"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate, epistart <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "hypertension") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

hypertension_post_combined <- hypertension_post %>%
  bind_rows(hypertension_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_hypertension <- hypertension_not %>%
  anti_join(hypertension_post_combined, by = "patid") %>%
  mutate(condition = "hypertension") %>%
  mutate(diff_months = 999)

hypertension_post_combined <- hypertension_post_combined %>%
  bind_rows(no_incident_hypertension) %>%
  bind_rows(hypertension_at_bl)

rm(hypertension_post, hypertension_post_hes, no_incident_hypertension, hypertension_at_bl, hypertension_not, SMIHypertension_C)

# Diabetes
diabetes_at_bl <- cv_events_cohort %>%
  filter(priordiabetes == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "diabetes", diff_months = 0)

diabetes_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priordiabetes) %>%
  filter(priordiabetes == "No")

load("SMIDiabetes_C.Rdata")

diabetes_post <- diabetes_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIDiabetes_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "diabetes") %>%
  select(patid, cohortentrydate, eventdate, condition)

diabetes_post_hes <- diabetes_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("E10", "E11", "E12", "E13", "E14", "E23.2", "H36.0", "H28.0"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate, epistart <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "diabetes") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

# Combine and filter for the minimum eventdate at the final step
diabetes_post_combined <- diabetes_post %>%
  bind_rows(diabetes_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_diabetes <- diabetes_not %>%
  anti_join(diabetes_post_combined, by = "patid") %>%
  mutate(condition = "diabetes") %>%
  mutate(diff_months = 999)

diabetes_post_combined <- diabetes_post_combined %>%
  bind_rows(no_incident_diabetes) %>%
  bind_rows(diabetes_at_bl)

rm(diabetes_post, diabetes_post_hes, no_incident_diabetes, diabetes_at_bl, diabetes_not, SMIDiabetes_C)

# Dyslipidaemia
dyslipid_at_bl <- cv_events_cohort %>%
  filter(priordyslipidaemia == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "dyslipidaemia", diff_months = 0)

dyslipid_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priordyslipidaemia) %>%
  filter(priordyslipidaemia == "No")

load("SMIDyslipidaemia_C.Rdata")

dyslipid_post <- dyslipid_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIDyslipidaemia_C, by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "dyslipidaemia") %>%
  select(patid, cohortentrydate, eventdate, condition)

dyslipid_post_hes <- dyslipid_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("E78.0", "E78.1", "E78.2", "E78.3", "E78.4", "E78.5"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate & epistart <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  mutate(condition = "dyslipidaemia") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

dyslipid_post_combined <- dyslipid_post %>%
  bind_rows(dyslipid_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_dyslipid <- dyslipid_not %>%
  anti_join(dyslipid_post_combined, by = "patid") %>%
  mutate(condition = "dyslipidaemia") %>%
  mutate(diff_months = 999)

dyslipid_post_combined <- dyslipid_post_combined %>%
  bind_rows(no_incident_dyslipid) %>%
  bind_rows(dyslipid_at_bl)

rm(dyslipid_post, dyslipid_post_hes, no_incident_dyslipid, dyslipid_at_bl, dyslipid_not, SMIDyslipidaemia_C)

# Substance misuse
load("SMIElixhauser_C.Rdata")

substanceabuse_at_bl <- cv_events_cohort %>%
  filter(priorsubstanceabuse == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "substanceabuse", diff_months = 0)

substanceabuse_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, priorsubstanceabuse) %>%
  filter(priorsubstanceabuse == "No")

substanceabuse_post <- substanceabuse_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIElixhauser_C %>%
               filter(condition == "Drug abuse" & !grepl("(?i)alcohol misuse", readterm)), by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "substanceabuse") %>%
  select(patid, cohortentrydate, eventdate, condition)

substanceabuse_post_hes <- substanceabuse_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("F11|F12|F13|F14|F15|F16|F18|F19|Z50.3|Z71.5|Z72.2"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate & epistart <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  mutate(condition = "substanceabuse") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

# Combine and filter for the minimum eventdate at the final step
substanceabuse_post_combined <- substanceabuse_post %>%
  bind_rows(substanceabuse_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_substanceabuse <- substanceabuse_not %>%
  anti_join(substanceabuse_post_combined, by = "patid") %>%
  mutate(condition = "substanceabuse") %>%
  mutate(diff_months = 999)

substanceabuse_post_combined <- substanceabuse_post_combined %>%
  bind_rows(no_incident_substanceabuse) %>%
  bind_rows(substanceabuse_at_bl)

rm(substanceabuse_post, substanceabuse_post_hes, no_incident_substanceabuse, substanceabuse_at_bl, substanceabuse_not)

# Alcohol
alcoholabuse_at_bl <- cv_events_cohort %>%
  filter(prioralcoholabuse == "Yes") %>%
  transmute(patid, cohortentrydate, condition = "alcoholabuse", diff_months = 0)

alcoholabuse_not <- cv_events_cohort %>%
  select(patid, cohortentrydate, prioralcoholabuse) %>%
  filter(prioralcoholabuse == "No")

alcoholabuse_post <- alcoholabuse_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIElixhauser_C %>%
               filter(condition == "Alcohol abuse"),
             by = "patid") %>%
  filter(eventdate > cohortentrydate, eventdate <= cohortentrydate + days(1826)) %>% # within 5 years
  mutate(condition = "alcoholabuse") %>%
  select(patid, cohortentrydate, eventdate, condition)

alcoholabuse_post_hes <- alcoholabuse_not %>%
  select(patid, cohortentrydate) %>%
  inner_join(
    hes_primarydiag %>%
      filter(grepl(paste0("(?i)", paste(c("F10|Z50.2|Z71.4|Z72.1"), collapse = "|")), ICD)),
    by = "patid") %>%
  filter(epistart > cohortentrydate & epistart <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  mutate(condition = "alcoholabuse") %>%
  select(patid, cohortentrydate, eventdate = epistart, condition)

# Combine and filter for the minimum eventdate at the final step
alcoholabuse_post_combined <- alcoholabuse_post %>%
  bind_rows(alcoholabuse_post_hes) %>%
  group_by(patid) %>%
  filter(eventdate == min(eventdate)) %>% # keep only the earliest event per patid
  ungroup() %>%
  distinct() %>%
  mutate(diff = eventdate - cohortentrydate) %>%
  mutate(diff_months = pmax(round(as.numeric(diff) / 30.25, 0), 1)) %>%
  select(-diff)

no_incident_alcoholabuse <- alcoholabuse_not %>%
  anti_join(alcoholabuse_post_combined, by = "patid") %>%
  mutate(condition = "alcoholabuse") %>%
  mutate(diff_months = 999)

alcoholabuse_post_combined <- alcoholabuse_post_combined %>%
  bind_rows(no_incident_alcoholabuse) %>%
  bind_rows(alcoholabuse_at_bl)

rm(alcoholabuse_post, alcoholabuse_post_hes, no_incident_alcoholabuse, alcoholabuse_at_bl, alcoholabuse_not, SMIElixhauser_C)

rm(cv_events_cohort, hes_primarydiag)

create_long_diagnosis_df <- function(df, diagnosis_name) {
  result <- df %>%
    select(patid, diff_months) %>%
    # Create a sequence of months for each patient
    rowwise() %>%
    mutate(month_data = list(seq(0, 59))) %>%
    unnest(month_data) %>%
    rename(time = month_data) %>%
    # Create diagnosis indicator based on diff_months
    mutate(
      !!paste0(diagnosis_name, "_indicator") := ifelse(time >= diff_months, 1, 0)) %>%
    # Select only relevant columns
    select(patid, time, !!paste0(diagnosis_name, "_indicator"))

  # Return the resulting data frame
  return(result)}

# Convert both dataframes to long format
diabetes_long <- create_long_diagnosis_df(diabetes_post_combined, "diabetes")
dyslipid_long <- create_long_diagnosis_df(dyslipid_post_combined, "dyslipidemia")
angina_long <- create_long_diagnosis_df(angina_post_combined, "angina")
hypertension_long <- create_long_diagnosis_df(hypertension_post_combined, "hypertension")
arrhythmia_long <- create_long_diagnosis_df(arrhythmia_post_combined, "arrhythmia")
liverdisease_long <- create_long_diagnosis_df(liverdisease_post_combined, "liverdisease")
renaldisease_long <- create_long_diagnosis_df(renaldisease_post_combined, "renaldisease")
substanceabuse_long <- create_long_diagnosis_df(substanceabuse_post_combined, "substanceabuse")
alcoholabuse_long <- create_long_diagnosis_df(alcoholabuse_post_combined, "alcoholabuse")

rm(diabetes_post_combined, dyslipid_post_combined, angina_post_combined, hypertension_post_combined,
   arrhythmia_post_combined, liverdisease_post_combined, renaldisease_post_combined, substanceabuse_post_combined, alcoholabuse_post_combined)

# Combine all long format dataframes
combined_long_df <- list(
  diabetes_long,
  dyslipid_long,
  angina_long,
  hypertension_long,
  arrhythmia_long,
  liverdisease_long,
  renaldisease_long,
  substanceabuse_long,
  alcoholabuse_long) %>%
  purrr::reduce(full_join, by = c("patid", "time")) %>%
  mutate(across(ends_with("_indicator"), ~ replace_na(., 0))) %>%
  mutate(sum_indicators = rowSums(select(., ends_with("_indicator"))))

rm(diabetes_long, dyslipid_long, angina_long, hypertension_long, arrhythmia_long, liverdisease_long,
   renaldisease_long, substanceabuse_long, alcoholabuse_long)

start <- Sys.time()

log_progress <- function(message) {
  futureCall(function(msg, dir) {
    log_file <- file.path(dir, "progress.log")
    cat(format(Sys.time()), msg, "\n", file = log_file, append = TRUE)
    cat(msg, "\n")
  }, list(message, output_dir))}

# Use function previously defined for bootstrap procedure (script 6)

# Package versions
packageVersion("tidyr") # 1.3.1
packageVersion("stringr") # 1.5.1
packageVersion("lubridate") # 1.9.4
packageVersion("data.table") # 1.16.4
packageVersion("reshape2") # 1.4.4
packageVersion("ggplot2") # 3.5.1
packageVersion("gtsummary") # 2.0.4
packageVersion("gt") # 0.11.1
packageVersion("mice") # 3.17.0
packageVersion("nnet") # 7.3.19
packageVersion("dplyr") # 1.1.4
packageVersion("scales") # 1.3.0
packageVersion("future") # 1.34.0
packageVersion("furrr") # 0.3.1
packageVersion("vroom") # 1.6.5