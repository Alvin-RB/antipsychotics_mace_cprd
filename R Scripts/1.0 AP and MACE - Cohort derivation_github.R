# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 13/11/2024

# DERIVE COHORT

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "readr", "readxl", "data.table", 
                   "purrr", "rlang", "reshape2", "doseminer", "drugprepr", "chlorpromazineR", "ggplot2", 
                   "gtsummary", "gt", "mice", "finalfit", "tidylog", "survminer", "survival", "vroom"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Load data

# CPRD
load("SMIcohort.rdata")

SMIcohort <- SMIcohort %>%
  select(patid, pracid, gender, yob, regstartdate, regenddate, lcd, ethnicity_cat_cprdhes, region, patprac_2019imd_quintile, pat_2019imd_quintile, prac_2019imd_quintile, hes_apc_e, ons_death_e, deathdate, age_atfirstdiag, first_diagnosis_date, first_smi_diag, age_atfirstdiag, oralap_rec_bin, injection_firstdate)

load("antipsychotics_combined.rdata")

antipsychotics_combined <- antipsychotics_combined %>%
  select(1:17) %>%
  select(-issueseq, -enterdate, -numpacks) %>%
  filter(AP != "Prochlorperazine") %>% 
  filter(issuedate < '2020-01-01')

# HES
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
                  !!!config,  
                  patid = col_character()),  
                trim_ws = TRUE)
  
  # Update patid based on database
  data <- data %>% mutate(patid = paste0(patid, suffix))
  
  return(data)}

hes_patient_gold_1 <- load_hes("hes_patient_21_000729.txt", hes_config$patient)
hes_patient_gold_2 <- load_hes("hes_patient_21_000729_request2.txt", hes_config$patient)
hes_patient_aurum <- load_hes("hes_patient_21_000729.txt", hes_config$patient)

hes_patient <- rbindlist(list(hes_patient_gold_1, hes_patient_gold_2, hes_patient_aurum))
remove(hes_patient_gold_1, hes_patient_gold_2, hes_patient_aurum)
gc()

hes_diagnosis_epi_gold_1 <- load_hes("hes_diagnosis_epi_21_000729.txt", hes_config$diagnosis_epi)
hes_diagnosis_epi_gold_2 <- load_hes("hes_diagnosis_epi_21_000729_request2.txt", hes_config$diagnosis_epi)
hes_diagnosis_epi_aurum <- load_hes("hes_diagnosis_epi_21_000729.txt", hes_config$diagnosis_epi)

hes_primarydiag <- rbindlist(list(hes_diagnosis_epi_gold_1, hes_diagnosis_epi_gold_2, hes_diagnosis_epi_aurum))
remove(hes_diagnosis_epi_gold_1, hes_diagnosis_epi_gold_2, hes_diagnosis_epi_aurum)
gc()

hes_hosp_gold_1 <- load_hes("hes_diagnosis_hosp_21_000729.txt", hes_config$hosp)
hes_hosp_gold_2 <- load_hes("hes_diagnosis_hosp_21_000729_request2.txt", hes_config$hosp)
hes_hosp_aurum <- load_hes("hes_diagnosis_hosp_21_000729.txt", hes_config$hosp)

hes_hosp <- rbindlist(list(hes_hosp_gold_1, hes_hosp_gold_2, hes_hosp_aurum))
remove(hes_hosp_gold_1, hes_hosp_gold_2, hes_hosp_aurum)
gc()

hes_episodes_gold_1 <- load_hes("hes_episodes_21_000729.txt", hes_config$episodes)
hes_episodes_gold_2 <- load_hes("hes_episodes_21_000729_request2.txt", hes_config$episodes)
hes_episodes_aurum <- load_hes("hes_episodes_21_000729.txt", hes_config$episodes)

hes_episodes <- rbindlist(list(hes_episodes_gold_1, hes_episodes_gold_2, hes_episodes_aurum)) 
remove(hes_episodes_gold_1, hes_episodes_gold_2, hes_episodes_aurum)
gc()

# DERIVE EXPOSURE COHORT ####

# Incident antipsychotic use

#Restrict to APs of interest
apsofinterest <- antipsychotics_combined %>%
  select(patid, AP, issuedate, deliverymethod) %>%
  distinct() %>%
  filter(AP %in% c("Aripiprazole", "Olanzapine", "Risperidone", "Quetiapine")) %>%
  group_by(patid, AP, deliverymethod) %>%
  mutate(min = min(issuedate),
         deliverymethod = str_to_lower(deliverymethod)) %>%
  ungroup()

# Ever received each AP of interest
extract_ap <- function(ap_name, ap_prefix) {
  apsofinterest %>%
    filter(AP == ap_name & min == issuedate) %>%
    pivot_wider(
      id_cols = patid,
      names_from = deliverymethod,
      values_from = issuedate,
      names_prefix = paste0(ap_prefix, "_first")
    ) %>%
    mutate(!!paste0(ap_prefix, "_ever") := 1) %>%
    select(patid, everything())}

ari <- extract_ap("Aripiprazole", "ari")
olanz <- extract_ap("Olanzapine", "olanz")
risp <- extract_ap("Risperidone", "risp")
quet <- extract_ap("Quetiapine", "quet")

#First AP of interest
first_wide <- apsofinterest %>%
  group_by(patid) %>%
  mutate(min = min(issuedate)) %>% 
  ungroup() %>%
  filter(min == issuedate) %>% # first prescription of an AP of interest
  distinct() %>%
  select(patid, AP) %>%
  reshape2::dcast(patid ~ AP) %>%
  mutate(across(2:5, ~ ifelse(.x > 1, 1, .x))) %>% # set values >1 to 1
  mutate(sum = rowSums(select(., Aripiprazole, Olanzapine, Risperidone, Quetiapine), na.rm = FALSE))

# Merge and code variables
cv_events_cohort <- list(ari, olanz, risp, quet, first_wide) %>% 
  purrr::reduce(full_join, by = "patid") %>%
  mutate(ari_firsttype = case_when(ari_ever == 1 & ari_firstoral == ari_firstinjection ~ "Oral & Injection",
                                   ari_ever == 1 & is.na(ari_firstoral) ~ "Injection",
                                   ari_ever == 1 & is.na(ari_firstinjection) ~ "Oral"),
         olanz_firsttype = case_when(olanz_ever == 1 & olanz_firstoral == olanz_firstinjection ~ "Oral & Injection",
                                     olanz_ever == 1 & is.na(olanz_firstoral) ~ "Injection",
                                     olanz_ever == 1 & is.na(olanz_firstinjection) ~ "Oral"),
         risp_firsttype = case_when(risp_ever == 1 & risp_firstoral == risp_firstinjection ~ "Oral & Injection",
                                    risp_ever == 1 & is.na(risp_firstoral) ~ "Injection",
                                    risp_ever == 1 & is.na(risp_firstinjection) ~ "Oral"),
         quet_firsttype = case_when(quet_ever == 1 ~ "Oral"),
         cohortentrydate = pmin(ari_firstoral, ari_firstinjection, risp_firstoral, risp_firstinjection, 
                                olanz_firstoral, olanz_firstinjection, quet_firstoral, na.rm = TRUE),
         trt_group = if_else(sum > 1, 0,
                             if_else(Aripiprazole == 1 & sum == 1, 1,
                                     if_else(Olanzapine == 1 & sum == 1, 2,
                                             if_else(Quetiapine == 1 & sum == 1, 3,
                                                     if_else(Risperidone == 1 & sum == 1, 4, 
                                                             NA_real_))))),
         trt_group = factor(trt_group, 
                            levels = c(0, 1, 2, 3, 4),
                            labels = c("Multiple", "Aripiprazole", "Olanzapine", "Quetiapine", "Risperidone"))) %>%
  select(-contains("firstoral"), -contains("firstinjection"))

remove(ari, olanz, risp, quet, first_wide, apsofinterest)

# ELIGIBLITY CRITERIA ####

# Inclusion criteria

# Eligible for linkage to HES/ONS
cv_events_cohort <- cv_events_cohort %>%
  full_join(SMIcohort, by = "patid") %>%
  filter(hes_apc_e == 1 & ons_death_e == 1) %>% #
  select(-hes_apc_e, -ons_death_e)

remove(SMIcohort)

#Received oral APs ever
cv_events_cohort <- cv_events_cohort %>%
  filter(oralap_rec_bin == 1)

receivedaps <- n_distinct(cv_events_cohort$patid)

#Received AP of interest
cv_events_cohort <- cv_events_cohort %>%
  filter(sum > 0) %>%
  select(-sum, -Aripiprazole, -Quetiapine, -Risperidone, -Olanzapine, -oralap_rec_bin)

receivedapsofinterest <- n_distinct(cv_events_cohort$patid)

# Received AP of interest 2005-2015
cv_events_cohort <- cv_events_cohort %>%
  filter(cohortentrydate >= '2005-01-01' & cohortentrydate < '2015-01-01')

receivedapsofinterestintimeperiod <- n_distinct(cv_events_cohort$patid)

# Restrict to age range
cv_events_cohort <- cv_events_cohort %>%
  mutate(prevcohortentry_year = lubridate::year(cohortentrydate),
         age_atprevcohortentry = prevcohortentry_year - yob) %>%
  filter(age_atprevcohortentry >= 40 & age_atprevcohortentry < 100)

agecriteria <- n_distinct(cv_events_cohort$patid)

#SMI diagnosis at time of cohort entry
load("PatSMI_C.Rdata")

diag_prev <- cv_events_cohort %>%
  select(patid, trt_group, cohortentrydate) %>%
  left_join(PatSMI_C %>%
              filter(eventdate > '1900-01-01'), # # remove invalid diagnosis dates
            by = "patid") %>%
  distinct() %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") <= 30) %>% # create date difference, Allows 30 days for potential recording delays
  group_by(patid) %>%
  mutate(max = max(eventdate)) %>%
  ungroup() %>%
  filter(max == eventdate) %>%
  group_by(patid) %>%
  mutate(schz_any = ifelse(group == "schizophrenia", 1, 0), # If more than one diagnosis on closest day, use schizophrenia, then bipolar, then other psychosis
         bipolar_any = ifelse(group == "bipolar", 1, 0),
         otherpsychosis_any = ifelse(group == "other psychosis", 1, 0)) %>%
  ungroup() %>%
  group_by(patid) %>%
  mutate(schz_sum = sum(schz_any),
         bipolar_sum = sum(bipolar_any),
         otherpsychosis_sum = sum(otherpsychosis_any),
         row = row_number(),
         n = n()) %>%
  ungroup() %>%
  filter(row == 1) %>%
  mutate(diag_prev = case_when(n > 1 & schz_sum > 0 ~ "schizophrenia",
                               n > 1 & bipolar_sum > 0 & schz_sum == 0 ~ "bipolar",
                               n > 1 & otherpsychosis_sum > 0 & schz_sum == 0 & bipolar_sum == 0 ~ "other psychosis",
                               TRUE ~ group)) %>%
  select(patid, diag_prev) %>%
  distinct() %>%
  mutate(diagosisatcohortentry = 1,
         diag_prev = factor(diag_prev, levels = c("schizophrenia", "bipolar", "other psychosis")))

# Restrict to those with SMI diagnosis at baseline (allowing 30 days for recording days)
cv_events_cohort <- cv_events_cohort %>%
  left_join(diag_prev, by = "patid") %>%
  filter(diagosisatcohortentry == 1) %>%
  select(-diagosisatcohortentry)

diagnosisatentry <- n_distinct(cv_events_cohort$patid)
remove(diag_prev, PatSMI_C)

# GP registration date is not after cohort entry, minimum 6 months registration required
cv_events_cohort <- cv_events_cohort %>%
  filter(difftime(cohortentrydate, regstartdate, units = "days") >= 180)

medhistory <- n_distinct(cv_events_cohort$patid)

# Multiple APs of interest
cv_events_cohort <- cv_events_cohort %>%
  filter(trt_group != "Multiple")

levels(cv_events_cohort$trt_group)
cv_events_cohort$trt_group <- droplevels(cv_events_cohort$trt_group, exclude = "Mulitple")

multiplepasofinterest <- n_distinct(cv_events_cohort$patid)

# Mulitple APs (one of interest with one not of interest)
napscohortentry <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  left_join(antipsychotics_combined, by = "patid") %>%
  filter(issuedate == cohortentrydate) %>%
  select(patid, AP) %>%
  distinct() %>%
  group_by(patid) %>%
  mutate(naps = n_distinct(AP)) %>%
  ungroup() %>%
  select(patid, naps) %>%
  distinct()

cv_events_cohort <- cv_events_cohort %>%
  left_join(napscohortentry) %>%
  filter(naps == 1) %>%
  select(-naps)

multipleapsnotofinterest <- n_distinct(cv_events_cohort$patid)

# First AP of interest was oral (not injection)
cv_events_cohort <- cv_events_cohort %>%
  mutate(notinjection = case_when(trt_group == "Olanzapine" & olanz_firsttype == "Oral" ~ 1,
                                  trt_group == "Aripiprazole" & ari_firsttype == "Oral" ~ 1,
                                  trt_group == "Quetiapine" & quet_firsttype == "Oral" ~ 1,
                                  trt_group == "Risperidone" & risp_firsttype == "Oral" ~ 1,
                                  TRUE ~ 0)) %>%
  filter(notinjection == 1) %>%
  select(-contains("firsttype"), -contains("_ever"), -notinjection)

firstapofinterestoral <- n_distinct(cv_events_cohort$patid)

# Most recent LAI
injection_mostrecent <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  left_join(antipsychotics_combined %>% filter(deliverymethod == "Injection"), by = "patid") %>%
  mutate(lastinjectionpriortobaseline = as.numeric(difftime(issuedate, cohortentrydate, units = "days"))) %>%
  filter(lastinjectionpriortobaseline <= 0) %>%
  group_by(patid) %>%
  filter(issuedate == max(issuedate)) %>%
  ungroup() %>%
  distinct(patid, lastinjectionpriortobaseline)

cv_events_cohort <- cv_events_cohort %>%
  left_join(injection_mostrecent, by = "patid") %>%
  filter(is.na(lastinjectionpriortobaseline) | lastinjectionpriortobaseline <= -90)

injectioninlast90days <- n_distinct(cv_events_cohort$patid)

remove(injection_mostrecent, napscohortentry)

# Incident prescription was not PRN
prncheck <- cv_events_cohort %>%
  inner_join(antipsychotics_combined, by = c("patid" = "patid", "cohortentrydate" = "issuedate")) %>%
  select(patid, trt_group, prn, dosage_text) %>%
  filter(!is.na(prn) | grepl("(?i)if needed|as needed|when needed|required|as req|prn|pro re nata", dosage_text)) %>%
  select(patid) %>%
  distinct() %>%
  mutate(prn = 1)

cv_events_cohort <- cv_events_cohort %>%
  left_join(prncheck, by = "patid") %>%
  filter(is.na(prn)) %>%
  select(-prn)

notprn <- n_distinct(cv_events_cohort$patid)
remove(prncheck)

# Remove patients with documented medical codes suggesting exposure to any one of the 4 APs prior to cohort entry
load(file = "SMIantipsychoticsmedical_C.Rdata")

tte_adversereactions <- cv_events_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(SMIAntipsychoticsMedical_C, by = "patid") %>%
  select(patid, cohortentrydate, trt_group, eventdate, readterm) %>%
  filter(eventdate < cohortentrydate) %>%
  filter(grepl("(?i)aripiprazole|quetiapine|olanzapine|risperidone", readterm)) %>%
  mutate(diff = cohortentrydate - eventdate) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorexposure = 1)

cv_events_cohort <- cv_events_cohort %>%
  left_join(tte_adversereactions, by = "patid") %>%
  filter(is.na(priorexposure)) %>%
  select(-priorexposure)

priorexposuresuggested <- n_distinct(cv_events_cohort$patid)
remove(SMIAntipsychoticsMedical_C, tte_adversereactions)

# Exclusion criteria

# Remained in cohort as of index date - had not died
cv_events_cohort <- cv_events_cohort %>%
  filter(cohortentrydate < deathdate | is.na(deathdate))

didnotdie <- n_distinct(cv_events_cohort$patid)

# Remained in cohort as of index date - had not deregistered
cv_events_cohort <- cv_events_cohort %>%
  filter(cohortentrydate < regenddate | is.na(regenddate))

remainedregistered <- n_distinct(cv_events_cohort$patid)

# Remained in cohort as of index date - before LCD
cv_events_cohort <- cv_events_cohort %>%
  filter(cohortentrydate < lcd)

notbeforelcd <- n_distinct(cv_events_cohort$patid)

# Dementia at baseline
load("SMIDementia_C.Rdata")
dementia_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIDementia_C, by = "patid") %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
  distinct(patid) %>%
  mutate(priordementia = 1)

dementia_hes <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag %>%
               filter(grepl("(?i)F00|F01|F02|F03|F09", ICD)),
             by = "patid") %>%
  filter(difftime(epistart, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priordementia_hes = 1)

dementia <- dementia_prev %>%
  full_join(dementia_hes, by = "patid") %>%
  select(patid) %>%
  distinct() %>%
  mutate(priordementia = 1)

# No dementia diagnosis
cv_events_cohort <- cv_events_cohort %>%
  left_join(dementia, by = "patid") %>%
  filter(is.na(priordementia)) %>%
  select(-priordementia)

priordementia <- n_distinct(cv_events_cohort$patid)

remove(dementia, dementia_hes, dementia_prev, SMIDementia_C)

#Identify cerebrovascular disease at baseline
load("SMIStroke_C.Rdata")

cerebrovascular_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIStroke_C, by = "patid")  %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorcerebrovasculardisease = 1)

cerebrovascular_hes <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag %>%
               filter(grepl("(?i)I60|I61|I62|I63|I64", ICD)),
             by = "patid") %>%
  filter(difftime(epistart, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorcerebrovasculardisease_hes = 1)

cerebrovascular <- cerebrovascular_prev %>%
  full_join(cerebrovascular_hes, by = "patid") %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorcerebrovasculardisease = 1)

# No prior stroke
cv_events_cohort <- cv_events_cohort %>%
  left_join(cerebrovascular, by = "patid") %>%
  filter(is.na(priorcerebrovasculardisease)) %>%
  select(-priorcerebrovasculardisease)

priorstroke <- n_distinct(cv_events_cohort$patid)

remove(cerebrovascular, cerebrovascular_hes, cerebrovascular_prev, SMIStroke_C)

#Identify MI at baseline
load("SMIMyocardialInfarction_C.Rdata")

myocardial_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIMyocardialInfarction_C, by = "patid")  %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priormyocardialinfarction = 1)

myocardial_hes <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag %>%
               filter(grepl("(?i)I21|I22|I23|I25.2", ICD)),
             by = "patid") %>%
  filter(difftime(epistart, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorcerebrovasculardisease_hes = 1)

myocardial <- myocardial_prev %>%
  full_join(myocardial_hes, by = "patid") %>%
  select(patid) %>%
  distinct() %>%
  mutate(priormyocardialinfarction = 1)

# No prior MI
cv_events_cohort <- cv_events_cohort %>%
  left_join(myocardial, by = "patid") %>%
  filter(is.na(priormyocardialinfarction)) %>%
  select(-priormyocardialinfarction)

remove(myocardial, myocardial_hes, myocardial_prev, SMIMyocardialInfarction_C)

priormi <- n_distinct(cv_events_cohort$patid)

# Identify duplicates by linking with HES
tte_hes_patient <- cv_events_cohort %>%
  select(patid, cohortentrydate, regstartdate, regenddate, lcd) %>%
  mutate(database = if_else(grepl("(?i)A", patid), "Aurum", "Gold")) %>%
  inner_join(hes_patient, by = "patid") %>%
  group_by(gen_hesid) %>%
  mutate(duplicate = n_distinct(patid)) %>% #
  ungroup() %>%
  filter(duplicate > 1) %>%
  mutate(is_longest_fup = if_else(!is.na(regenddate), regenddate - cohortentrydate, lcd - cohortentrydate)) %>%
  group_by(gen_hesid) %>%
  mutate(earliest_date = min(cohortentrydate),
         is_earliest = if_else(cohortentrydate == earliest_date, 1, 0),
         all_earliest = all(is_earliest == 1),
         is_aurum = if_else(database == "Aurum", 1, 0),
         all_aurum = all(is_aurum == 1),
         is_longest = ifelse(is_longest_fup == max(is_longest_fup), 1, 0),
         all_longest = all(is_longest == 1)) %>%  
  ungroup() %>%
  mutate(
    keep_earliest = case_when(all_earliest == FALSE & is_earliest == 1 ~ 1, # where the cohort entry dates are different, keep the earliest
                              TRUE ~ 0),
    keep_aurum = case_when(all_earliest == TRUE & all_aurum == FALSE & is_aurum == 1 ~ 1, # where the dates are the same, keep the aurum record
                           TRUE ~ 0),
    keep_longest = case_when(all_earliest == TRUE & all_aurum == TRUE & is_longest == 1 ~ 1, # where both earliest, both aurum, use the one with longest follow-up (verified that this ensures that just 1 of each duplicate pair is kept)
                             TRUE ~ 0)) %>%
  mutate(keep = case_when(keep_earliest == 1 ~ 1,
                          keep_aurum == 1 ~ 1,
                          keep_longest == 1 ~ 1,
                          TRUE ~ 0))

tte_hes_patient <- tte_hes_patient %>%
  select(patid, gen_hesid, keep)

cv_events_cohort <- cv_events_cohort %>%
  left_join(tte_hes_patient, by = "patid") %>%
  filter(is.na(keep) | keep == 1) %>%
  select(-keep)

remove(tte_hes_patient, hes_patient)

# Check if there are any remaining duplicates (there are none)
dupe_check <- cv_events_cohort %>%
  select(patid, gen_hesid) %>%
  group_by(gen_hesid) %>%
  mutate(duplicate = ifelse(!is.na(gen_hesid), n_distinct(patid), NA),
         n = ifelse(!is.na(gen_hesid), row_number(), NA)) %>%
  ungroup() %>%
  filter(n > 1)

if (nrow(dupe_check) > 0) {
  message("Warning: duplicates found")
} else {
  message("No duplicates found")}

remove(dupe_check)

deduplicated <- n_distinct(cv_events_cohort$patid)

# CONSORT table
consort1 <- data.frame(level = 1, Included = receivedaps, Excluded = 0, Reason = "Ever diagnosed with an SMI & prescribed an oral AP")
consort2 <- data.frame(level = 2, Included = receivedapsofinterest, Excluded = receivedaps - receivedapsofinterest, Reason = "Ever prescribed an AP of interest")
consort3 <- data.frame(level = 3, Included = receivedapsofinterestintimeperiod, Excluded = receivedapsofinterest - receivedapsofinterestintimeperiod, Reason = "First prescribed an AP of interest in study period")
consort4 <- data.frame(level = 4, Included = agecriteria, Excluded = receivedapsofinterestintimeperiod - agecriteria, Reason = "Aged 40-99")
consort5 <- data.frame(level = 5, Included = diagnosisatentry, Excluded = agecriteria - diagnosisatentry, Reason = "Had an SMI diagnosis")
consort6 <- data.frame(level = 6, Included = medhistory, Excluded = diagnosisatentry - medhistory, Reason = "Registered at primary care practice for at least 6m")
consort7 <- data.frame(level = 7, Included = multiplepasofinterest, Excluded = medhistory - multiplepasofinterest, Reason = "Not prescribed >1 AP of interest")
consort8 <- data.frame(level = 9, Included = multipleapsnotofinterest, Excluded = multiplepasofinterest - multipleapsnotofinterest, Reason = "Not prescribed >1 AP")
consort9 <- data.frame(level = 8, Included = firstapofinterestoral, Excluded = multipleapsnotofinterest - firstapofinterestoral, Reason = "Not prescribed an AP of interest as an LAI")
consort10 <- data.frame(level = 10, Included = injectioninlast90days, Excluded = firstapofinterestoral - injectioninlast90days, Reason = "Not prescribed an LAI AP in last 90 days")
consort11 <- data.frame(level = 11, Included = notprn, Excluded = injectioninlast90days - notprn, Reason = "First prescription was not PRN")
consort12 <- data.frame(level = 12, Included = priorexposuresuggested, Excluded = notprn - priorexposuresuggested, Reason = "No medical codes indicating receipt of an AP of interest prior to first identified prescription")
consort13 <- data.frame(level = 13, Included = didnotdie, Excluded = priorexposuresuggested - didnotdie, Reason = "Did not die")
consort14 <- data.frame(level = 14, Included = remainedregistered, Excluded = didnotdie - remainedregistered, Reason = "Did not end primary care registration")
consort15 <- data.frame(level = 15, Included = notbeforelcd, Excluded = remainedregistered - notbeforelcd, Reason = "Index date is not after practice last collection date")
consort16 <- data.frame(level = 16, Included = priordementia, Excluded = notbeforelcd - priordementia, Reason = "No dementia diagnosis")
consort17 <- data.frame(level = 17, Included = priorstroke, Excluded = priordementia - priorstroke, Reason = "No prior stroke")
consort18 <- data.frame(level = 18, Included = priormi, Excluded = priorstroke - priormi, Reason = "No prior MI")
consort19 <- data.frame(level = 19, Included = deduplicated, Excluded = priormi - deduplicated, Reason = "De-duplicated following HES linkage")
consort20 <- cv_events_cohort %>%
  count(trt_group) %>%
  rename(Included = n) %>%
  mutate(Reason = trt_group, Excluded = 0, level = 19 + row_number()) %>%
  select(level, Included, Excluded, Reason) %>%
  filter(Reason != "Multiple")

consort <- rbind(consort1, consort2, consort3, consort4, consort5, consort6, consort7, consort8, consort9, consort10, 
                 consort11, consort12, consort13, consort14, consort15, consort16, consort17, consort18, consort19, consort20) %>%
  gt()

rm(list = ls(pattern = "consort"))

# BASELINE CHARACTERISTICS ####

cv_events_cohort <- cv_events_cohort %>%
  mutate(smidiag_to_cohortentry = round(as.numeric(cohortentrydate - first_diagnosis_date) / 365.25, 2), # Time from SMI diagnosis to index date
         ethnicity_cat_cprdhes_4 = case_when(ethnicity_cat_cprdhes == "mixed" ~ "Mixed/Other", # Ethnicity with four levels (combined mixed/other due to low frequences)
                                             ethnicity_cat_cprdhes == "other" ~ "Mixed/Other",
                                             TRUE ~ ethnicity_cat_cprdhes)) %>%
  mutate(ethnicity_cat_cprdhes_4 = str_to_title(ethnicity_cat_cprdhes_4))

#Baseline antipsychotic dose
#Clean dosage texts with lookup tables
lookup <- read_excel("Data files/Misc/lookup.xlsx")
conditional_lookup <- read_excel("Data files/Misc/conditional_lookup.xlsx")

apdose_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(antipsychotics_combined, by = c("patid" = "patid", "cohortentrydate" = "issuedate")) %>%
  filter(AP != "Prochlorperazine") %>% 
  rename(prodcode = prodcodeid) %>%
  mutate(dosage_text_raw = dosage_text) %>%
  mutate(dosage_text = gsub("\\*", "", dosage_text)) %>%
  left_join(conditional_lookup, by = c("dosage_text", "productname")) %>%
  mutate(success = ifelse(is.na(code), 0, 1),
         dosage_text = ifelse(success == 1, code, dosage_text)) %>%
  select(-success, -code) %>%
  left_join(lookup, by = "dosage_text") %>%
  mutate(success = ifelse(is.na(code), 0, 1),
         dosage_text = ifelse(success == 1, code, dosage_text)) %>%
  select(-success, -code, -`Dupe check`) %>%
  distinct() # consider unique prescriptions only

#Extract free text prescription instructions and convert to numeric values, using doseminer, then merge with original data
free_text <- with(apdose_prev, dosage_text[!duplicated(dosage_text) & nchar(dosage_text) > 0])
extracted <- doseminer::extract_from_prescription(free_text)
apdose_prev <- merge(extracted, apdose_prev, by.x = 'raw', by.y = 'dosage_text', all = TRUE)

# Create new variables as per doseminer instructions
apdose_prev <- apdose_prev %>%
  separate(dose, c('min_dose', 'max_dose'), sep = '-',
           convert = TRUE, fill = 'right') %>%
  separate(itvl, c('min_itvl', 'max_itvl'), sep = '-',
           convert = TRUE, fill = 'right') %>%
  separate(freq, c('min_freq', 'max_freq'), sep = '-',
           convert = TRUE, fill = 'right') %>%
  mutate(dose = coalesce((min_dose + max_dose) / 2, min_dose),
         itvl = coalesce((min_itvl + max_itvl) / 2, min_itvl, 1),
         freq = coalesce((min_freq + max_freq) / 2, min_freq),
         ndd = freq * dose / itvl,
         dose_num = readr::parse_number(strength), # create numeric value of dose
         dose_num = case_when(grepl("micro", productname) | grepl("micro", strength) ~ dose_num/1000, #Convert mcg to mg
                              TRUE ~ dose_num)) %>%
  select(-min_dose, -max_dose, -min_itvl, -max_itvl, -min_freq, -max_freq, -unit) %>%
  select(patid, cohortentrydate, AP, productname, formulation, route, dosage_text_raw, dosage_text_processed = raw, strength, dose_num, everything())

apdose_prev <- apdose_prev %>%
  mutate(ndd_calc = case_when(!grepl("ml", strength) & quantity > 1 & quantity == duration & is.na(ndd) ~ round(quantity / duration, 1), # if qty/dur are the same & > 1 & AP is not prescribed in ml, calculate NDD
                              !grepl("ml", strength) & duration > 6 & duration <= 122 & quantity > 1 & quantity < 500 & is.na(ndd) ~ round(quantity / duration, 1), # if qty/dur are not the same, apply some limits
                              TRUE ~ NA_real_), # if criteria not met, set to NA
         ndd_calc = round(ndd_calc/.25)*.25, # round to nearest 0.25
         ndd_calc = case_when(ndd_calc == 0.0 ~ NA_real_, # if rounded value is 0, set to missing
                              TRUE ~ ndd_calc),
         ndd = coalesce(ndd, ndd_calc)) %>% # use newly calculated ndd if original ndd is missing
  select(-ndd_calc)

# Times dose in mg by number of tablets taken per day
apdose_prev <- apdose_prev %>%
  mutate(total_daily_dose = dose_num * ndd)

#Convert to chlorpromazine and olanzapine equivalent doses, using chlorpromazineR 
apdose_prev <- chlorpromazineR::to_ap(apdose_prev, convert_to_ap = "olanzapine", 
                                      convert_to_route = "oral", ap_label = "AP", 
                                      dose_label = "total_daily_dose", key = leucht2016)

# Where mulitple prescriptions on the same day, consider up to three prescriptions of each AP per patient, 
# then sum the total dose per prescription date, then filter to one row per patient
# where there are multiple prescriptions on the date, provide dose only if there is at least one known NDD
apdose_prev <- apdose_prev %>%
  group_by(patid, cohortentrydate, AP) %>%
  arrange(., desc(ndd)) %>%
  mutate(num = row_number()) %>%
  ungroup() %>%
  filter(num < 3) %>%
  group_by(patid, cohortentrydate) %>%
  mutate(doseknownforall = !any(is.na(ndd)),
         dose_total_rawunit = ifelse(doseknownforall != FALSE, sum(total_daily_dose), NA),
         dose_total_olanz = ifelse(doseknownforall != FALSE, sum(ap_eq), NA),
         dose_total_cpz = ifelse(doseknownforall != FALSE, sum(cpz_eq), NA),
         partial_dose_any = ifelse(any(ndd > 0), 1, 0),
         partial_dose_olanz = ifelse(partial_dose_any == 1, sum(ap_eq, na.rm = TRUE), NA),
         partial_dose_rawunit = ifelse(partial_dose_any == 1, sum(total_daily_dose, na.rm = TRUE), NA),
         partial_dose_cpz = ifelse(partial_dose_any == 1, sum(cpz_eq, na.rm = TRUE), NA),
         num = row_number()) %>%
  ungroup() %>%
  filter(num == 1) %>%
  select(patid, baselinedose_raw = partial_dose_rawunit, baselinedose_cpz = partial_dose_cpz, baselinedose_olanz = partial_dose_olanz)

# Prior AP in last 2 years (issuedate < cohortentrydate used as only want to count the non-study APs)
priorapuse_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(antipsychotics_combined, by = "patid") %>%
  select(patid, cohortentrydate, issuedate) %>%
  filter(issuedate >= cohortentrydate - days(731) & issuedate < cohortentrydate) %>%
  mutate(apuse_prior2years = 1) %>%
  select(patid, apuse_prior2years) %>%
  distinct()

# Switch and discontinuation
# Identify switch events
switch <- cv_events_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(antipsychotics_combined, by = "patid") %>%
  filter(issuedate >= cohortentrydate & issuedate < '2020-01-01') %>% # filter to records within study period
  filter(issuedate <= cohortentrydate + days(1826)) %>% # filter to records within 5 years for each patient
  filter(AP != "Levomepromazine") %>% # do not count this AP as treatment failure (Prochlorperazine has already been excluded)
  select(patid, cohortentrydate, trt_group, issuedate, AP) %>%
  arrange(patid, issuedate)

# Initialize data frame to store switch sequences
switch_sequence <- data.frame(patid = character(),
                              AP = character(),
                              seq_id = numeric(),
                              start_date = as.Date(character()),
                              stringsAsFactors = FALSE)

# Loop through each unique patid
for (participant in unique(switch$patid)) {
  
  # Subset data for the current participant
  participant_data <- switch[switch$patid == participant, ]
  
  # Initialize variables to store sequence information
  seq_id <- 1
  start_date <- participant_data$issuedate[1]
  current_drug <- participant_data$AP[1]
  
  # Loop through each prescription for the current participant
  for (i in 2:nrow(participant_data)) {
    if (!is.na(participant_data$AP[i]) && !is.na(current_drug) && participant_data$AP[i] != current_drug) {
      # End of a sequence detected
      switch_sequence <- rbind(switch_sequence, 
                               data.frame(patid = participant,
                                          AP = current_drug,
                                          seq_id = seq_id,
                                          start_date = start_date))
      
      # Update variables for the next sequence
      seq_id <- seq_id + 1
      start_date <- participant_data$issuedate[i]
      current_drug <- participant_data$AP[i]
    } else if (is.na(current_drug)) {
      # Initialize current_drug with the first valid value
      current_drug <- participant_data$AP[i]}}
  
  # Add the last sequence for the current participant
  switch_sequence_full <- rbind(switch_sequence, 
                                data.frame(patid = participant,
                                           AP = current_drug,
                                           seq_id = seq_id,
                                           start_date = start_date))}

# Identify discontuation events
discontinuation <- cv_events_cohort %>%
  select(patid, cohortentrydate, trt_group, regenddate, lcd, deathdate) %>%
  inner_join(antipsychotics_combined, by = "patid") %>%
  filter(issuedate >= cohortentrydate & issuedate < '2020-01-01') %>% # filter to records within study period
  filter(issuedate <= cohortentrydate + days(1826)) %>% # filter to records within 5 years for each patient
  select(patid, cohortentrydate, trt_group, issuedate, productname, AP, deliverymethod, strength, regenddate, lcd, deathdate) %>%
  arrange(patid, issuedate)

extract_discontinuation <- function(drug_name) {
  
  discontinuation %>%
    filter(trt_group == drug_name & AP == drug_name) %>% # process each medication separately
    group_by(patid) %>% 
    arrange(issuedate) %>%
    mutate(diff = as.numeric(issuedate - lag(issuedate)),
           first180 = if_else(diff >= 180, 1, 0),
           last = if_else(row_number() == n(), 1, 0),
           discontinue_date = if_else(last == 1, issuedate + days(180), NA_Date_)) %>%
    ungroup() %>%
    filter(first180 == 1 & last == 0 | last == 1 & first180 == 1 | last == 1 & first180 == 0 | last == 1 & is.na(first180)) %>% # identify first 180 day gap or last record
    group_by(patid) %>%
    mutate(min = min(issuedate)) %>%
    ungroup() %>%
    filter(min == issuedate) %>%
    select(-min) %>%
    group_by(patid) %>%
    mutate(n = row_number()) %>%
    ungroup() %>%
    filter(n == 1) %>% # filter one row per patient as assumed if they have more than one prescription on same day that it is on the same duration
    select(-n) %>%
    mutate(discontinue_date = issuedate + days(180), # generate discontinuation date
           afterdeath = discontinue_date > deathdate, # flag if after death
           afterregenddate = discontinue_date > regenddate, # flag if after regenddate
           afterlcd = discontinue_date > lcd, # flag if after lcd
           afterstudyend = discontinue_date > '2020-01-01', # flag if after study end date
           discontinue_date = case_when(afterdeath == TRUE | afterregenddate == TRUE | afterlcd == TRUE | afterstudyend == TRUE ~ NA, # if after any of the above dates, do not consider as discontinued
                                        TRUE ~ discontinue_date), # otherwise use discontinuation date
           discontinued = case_when(afterdeath == TRUE | afterregenddate == TRUE | afterlcd == TRUE ~ 0, # if after any of the above dates, do not consider as discontinued
                                    last == 1 & diff >= 90 ~ 1,
                                    last == 1 & is.na(diff) ~ 1,
                                    first180 == 1 ~ 1,
                                    TRUE ~ 0), # only had one prescription
           daystodiscontinuation = case_when(discontinued == 1 ~ as.numeric(discontinue_date - cohortentrydate),
                                             TRUE ~ NA),
           discontinue_date = case_when(discontinued == 1 ~ discontinue_date,
                                        TRUE ~ NA)) %>%
    select(patid, discontinued, discontinue_date, daystodiscontinuation)}

ari_discontinue <- extract_discontinuation("Aripiprazole")
olanz_discontinue <- extract_discontinuation("Olanzapine")
risp_discontinue <- extract_discontinuation("Risperidone")
quet_discontinue <- extract_discontinuation("Quetiapine")

discontinuation <- rbind(ari_discontinue, olanz_discontinue, quet_discontinue, risp_discontinue) %>%
  rename(start_date = discontinue_date) %>%
  select(-daystodiscontinuation) %>%
  filter(discontinued != 0)

switchanddiscontinuation <- switch_sequence_full[-1, ] %>% # Remove the first row (initialized empty row) from switch_sequence
  filter(seq_id == 2) %>% # seq_id 2 is the first switch
  bind_rows(discontinuation) %>%
  mutate(event = case_when(discontinued == 1 ~ "discontinued",
                           seq_id == 2 ~ "switched",
                           TRUE ~ NA)) %>%
  group_by(patid, start_date) %>%
  mutate(notadherent_event = paste(unique(event), collapse = "/")) %>%
  ungroup() %>%
  group_by(patid) %>%
  mutate(min = min(start_date)) %>%
  ungroup() %>%
  filter(min == start_date) %>%
  select(patid, notadherent_event, notadherent_date = start_date) %>%
  distinct()

load(file = "SMIantipsychoticsmedical_C.Rdata")

tte_adversereactions <- cv_events_cohort %>%
  select(patid, cohortentrydate, trt_group) %>%
  inner_join(SMIAntipsychoticsMedical_C, by = "patid") %>%
  select(patid, cohortentrydate, trt_group, eventdate, readterm) %>%
  filter(grepl("(?i)aripiprazole|quetiapine|olanzapine|risperidone", readterm)) %>%
  filter(eventdate >= cohortentrydate & eventdate < '2020-01-01') %>% # filter to records within study period
  filter(eventdate <= cohortentrydate + days(1826)) # filter to records within 5 years for each patient

switchdisc_adversereaction <- switchanddiscontinuation %>%
  bind_rows(tte_adversereactions) %>%
  group_by(patid) %>%
  # Fill down notadherent_event and notadherent_date within each patid
  fill(notadherent_event, notadherent_date, .direction = "downup") %>%
  mutate(
    had_rel_ae = as.integer(
      !is.na(eventdate) &
        (grepl("(?i)aripiprazole", readterm) & trt_group == "Aripiprazole" |
           grepl("(?i)quetiapine", readterm) & trt_group == "Quetiapine" |
           grepl("(?i)olanzapine", readterm) & trt_group == "Olanzapine" |
           grepl("(?i)risperidone", readterm) & trt_group == "Risperidone"))) %>%
  mutate(
    within_switch = as.integer(
      any(had_rel_ae == 1) & 
        notadherent_event == "switched" &
        !is.na(eventdate) &
        as.numeric(difftime(eventdate, notadherent_date, units = "days")) <= 30),
    within_disc = as.integer(
      any(had_rel_ae == 1) & 
        notadherent_event == "discontinued" &
        !is.na(eventdate) &
        as.numeric(difftime(eventdate, notadherent_date, units = "days")) <= 30)) %>%
  mutate(
    within_switch = coalesce(within_switch, 0),
    within_disc = coalesce(within_disc, 0),
    recode = as.integer(
      any(within_switch == 1 | within_disc == 1))) %>%
  ungroup() %>%
  select(patid, notadherent_event, notadherent_date, recode) %>%
  distinct() %>%
  filter(recode != 1) %>%
  select(-recode)

rm(list = ls(pattern = "_switch|_discontinue"), switch, discontinuation, switchanddiscontinuation, switch_sequence, switch_sequence_full, antipsychotics_combined, 
   extracted, lookup, conditional_lookup, apdose_fup, dosevars, DDD_consensus, participant_data, SMIAntipsychoticsMedical_C, tte_adversereactions)

# Concomitant medications in past two years
extract_baseline_medications <- function(medication_data, var_name) {
  
  # Construct file name to load  medication data
  load(paste0("", medication_data, ".Rdata"))
  
  # Create a new data frame by merging some fields from cv_events_cohort with the loaded medication data
  medication_prev <- cv_events_cohort %>%
    select(patid, cohortentrydate) %>%
    inner_join(get(medication_data) %>% select(patid, issuedate), by = "patid") %>%
    select(patid, cohortentrydate, issuedate) %>%
    filter(issuedate >= cohortentrydate - days(731) & issuedate <= cohortentrydate) %>%
    mutate(!!var_name := 1) %>%
    select(patid, var_name) %>%
    distinct()
  
  remove(medication_data)
  gc()
  
  # Return processed data frame
  return(medication_prev)}

antidepressants_prev <- extract_baseline_medications("SMIantidepressants", "antidepressant_prior2years")
moodstab_prev <- extract_baseline_medications("SMImoodstabilisers_C", "moodstab_prior2years")
lipiddrugs_prev <- extract_baseline_medications("SMIlipiddrugs_C", "lipiddrugs_prior2years")
hypertensiondrugs_prev <- extract_baseline_medications("SMIhypertensionmeds_C", "hypertensiondrugs_prior2years")
antidiabetics_prev <- extract_baseline_medications("SMIantidiabetics_C", "antidiabetics_prior2years")
anticoagulants_prev <- extract_baseline_medications("SMIanticoagulants_C", "anticoagulants_prior2years")
antiplatelets_prev <- extract_baseline_medications("SMIantiplatelets_C", "antiplatelets_prior2years")
zdrugs_prev <- extract_baseline_medications("SMIzdrugs_C", "zdrugs_prior2years")
anxiolytics_prev <- extract_baseline_medications("SMIanxiolytics_C", "anxiolytics_prior2years")
benzodiazepines_prev <- extract_baseline_medications("SMIbenzodiazepines_C", "benzodiazepines_prior2years")

# Comorbidities

extract_baseline_comorbidities <- function(diagnosis_data, ICD_codes, var_name) {
  
  # Construct file name to load diagnosis data
  load(paste0("", diagnosis_data, ".Rdata"))
  
  # CPRD data
  diagnosis_prev <- cv_events_cohort %>%
    select(patid, cohortentrydate) %>%
    inner_join(get(diagnosis_data), by = "patid") %>%
    filter(difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
    distinct(patid) %>%
    mutate(priordiagnosis = 1)
  
  # HES data
  diagnosis_hes <- cv_events_cohort %>%
    select(patid, cohortentrydate) %>%
    inner_join(hes_primarydiag %>%
                 filter(grepl(paste0("(?i)", paste(ICD_codes, collapse = "|")), ICD)), 
               by = "patid") %>%
    filter(difftime(epistart, cohortentrydate, units = "days") <= 0) %>%
    distinct(patid) %>%
    mutate(priordiagnosis_hes = 1)
  
  diagnosis <- diagnosis_prev %>%
    full_join(diagnosis_hes, by = "patid") %>%
    distinct(patid) %>%
    mutate(!!var_name := 1)
  
  remove(diagnosis_data)
  gc()
  
  return(diagnosis)}

# Call the function for baseline diagnoses
diabetes_prev <- extract_baseline_comorbidities("SMIDiabetes_C", c("E10", "E11", "E12", "E13", "E14", "E23.2", "H36.0", "H28.0"), "priordiabetes")
angina_prev <- extract_baseline_comorbidities("SMIAngina_C", "I20", "priorangina")
htn_prev <- extract_baseline_comorbidities("SMIHypertension_C", c("I10", "I11", "I12", "I13", "I14", "I15", "K76.6"), "priorhypertension")
arrhythmia_prev <- extract_baseline_comorbidities("SMIArrhythmia_C", c("I47", "I48", "I49"), "priorarrhythmia")
liver_prev <- extract_baseline_comorbidities("SMILiverDisease_C", c("N17", "N18", "N19"), "priorliverdisease")
renal_prev <- extract_baseline_comorbidities("SMIRenalDisease_C", "K7", "priorrenaldisease")
dyslipidaemia_prev <- extract_baseline_comorbidities("SMIDyslipidaemia_C", c("E78.0", "E78.1", "E78.2", "E78.3", "E78.4", "E78.5"), "priordyslipidaemia")

load("SMICharlson_C.Rdata")
load("SMIElixhauser_C.Rdata")

# Substance abuse at baseline
substance_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIElixhauser_C %>%
               filter(condition == "Drug abuse" & !grepl("(?i)alcohol misuse", readterm)), by = "patid") %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorsubstanceabuse = 1)

substance_hes <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag %>%
               filter(grepl("(?i)F11|F12|F13|F14|F15|F16|F18|F19|Z50.3|Z71.5|Z72.2", ICD)), by = "patid") %>%
  filter(difftime(epistart, cohortentrydate, units = "days") <= 0) %>%
  distinct(patid) %>%
  mutate(priorsubstanceabuse_hes = 1)

substance_prev <- substance_prev %>%
  full_join(substance_hes, by = "patid") %>%
  distinct(patid) %>%
  mutate(priorsubstanceabuse = 1)

#Alcohol abuse at baseline
alcohol_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIElixhauser_C %>%
               filter(condition == "Alcohol abuse"),
             by = "patid") %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(prioralcoholabuse = 1)

alcohol_hes <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag %>%
               filter(grepl("(?i)F10|Z50.2|Z71.4|Z72.1", ICD)), 
             by = "patid") %>%
  filter(difftime(epistart, cohortentrydate, units = "days") <= 0) %>%
  distinct(patid) %>%
  mutate(prioralcoholabuse_hes = 1)

alcohol_prev <- alcohol_prev %>%
  full_join(alcohol_hes, by = "patid") %>%
  distinct(patid) %>%
  mutate(prioralcoholabuse = 1)

remove(SMICharlson_C, SMIElixhauser_C)
remove(alcohol_hes, substance_hes)
gc()

# Number of GP consults in the last six months
load("SMIGPConsults_C.Rdata")

gpconsults_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIGPConsults_C, by = "patid") %>%
  filter(difftime(eventdate, cohortentrydate, units = "days") >= -180 & difftime(eventdate, cohortentrydate, units = "days") <= 0) %>%
  select(patid, eventdate) %>%
  distinct() %>%
  group_by(patid) %>%
  mutate(gpconsults_last6m = n()) %>%
  ungroup() %>%
  select(patid, gpconsults_last6m) %>%
  distinct()

remove(SMIGPConsults_C)

# Psychiatric hospitalisation
# Primary psychiatric diagnosis, at the episode level
tte_hes_primarydiag <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag, by = "patid") %>%
  filter(d_order == 1) %>% # filter to primary diagnosis 
  filter(grepl("(?i)F", ICD)) # select diagnoses of interest (any F diagnosis)

#Self-harm at any diagnosis level (as not recorded as primary diagnosis in HES)
tte_hes_selfharm <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag, by = "patid") %>%
  filter(grepl("(?i)X6|X7|X80|X81|X82|X83|X84|Y87.0", ICD))

# Data frame with psych hospital admissions within 2 years before cohort entry date
priorpsychhosp2y_prev <- tte_hes_primarydiag %>%
  bind_rows(tte_hes_selfharm) %>% # Bind primary psych and self-harm admissions
  select(patid, cohortentrydate, spno, ICD) %>%
  group_by(patid, spno) %>%
  mutate(eligible_ICD = paste(unique(ICD), collapse = "/"), # paste all primary diagnoses per patient and spellno into one column
         eligible_ICD = gsub("^/|/$", "", eligible_ICD)) %>% 
  ungroup() %>%
  select(-ICD) %>%
  distinct() %>%
  left_join(select(hes_hosp, patid, spno, admidate), by = c("patid", "spno")) %>% # link to hospital level info using spno to retrieve the admission date
  left_join(select(hes_episodes, patid, admidate, spno, classpat), by = c("patid", "spno", "admidate")) %>% # link to hospital level info using spno to retrieve the admission date
  filter(classpat == 1) %>% # remove regular day/night attendances and baby delivery admissions only
  select(-classpat, -spno) %>%
  distinct() %>%
  filter(admidate < cohortentrydate & admidate >= cohortentrydate - days(730)) %>% # filter to records within 2 years
  group_by(patid) %>%
  mutate(priorhosp_psych_any = 1,
         priorpsychhosp_num = n_distinct(admidate)) %>%
  ungroup() %>%
  select(patid, priorhosp_psych_any, priorpsychhosp_num) %>%
  distinct()

# Physical health admissions
priorphyshosp2y_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag, by = "patid") %>%
  filter(d_order == 1) %>% # filter to primary diagnosis 
  filter(!grepl("(?i)F|X6|X7|X80|X81|X82|X83|X84|Y87.0", ICD)) # remove psych diagnoses
  
self_harm_seplls <- tte_hes_selfharm %>%
  select(patid, spno) %>%
  mutate(selfharm = 1)

# Physical health hospital admissions within 2 years before cohort entry date
priorphyshosp2y_prev <- priorphyshosp2y_prev %>%
  left_join(self_harm_seplls, by = c("patid", "spno")) %>%
  filter(is.na(selfharm)) %>% # remove self-harm episodes
  select(patid, cohortentrydate, spno, ICD) %>%
  group_by(patid, spno) %>%
  mutate(eligible_ICD = paste(unique(ICD), collapse = "/"), # paste all primary diagnoses per patient and spellno into one column
         eligible_ICD = gsub("^/|/$", "", eligible_ICD)) %>% # remove leading/trailing slashes
  ungroup() %>%
  select(-ICD) %>%
  distinct() %>%
  left_join(select(hes_episodes, patid, admidate, spno, classpat), by = c("patid", "spno")) %>% # link to hospital level info using spno to retrieve the admission date
  filter(classpat == 1) %>% # remove regular day/night attendances and baby delivery admissions only
  select(-classpat, -spno) %>%
  distinct() %>%
  filter(admidate < cohortentrydate & admidate >= cohortentrydate - days(730)) %>% # filter to records within 2 years
  group_by(patid) %>%
  mutate(priorhosp_physical_any = 1,
         priorhosp_physical_num = n_distinct(admidate)) %>%
  ungroup() %>%
  select(patid, priorhosp_physical_any, priorhosp_physical_num) %>%
  distinct()

# Merge
df_list <- list(cv_events_cohort, apdose_prev, switchdisc_adversereaction, priorapuse_prev,
                antidepressants_prev, moodstab_prev, lipiddrugs_prev, hypertensiondrugs_prev, antidiabetics_prev, anticoagulants_prev,  antiplatelets_prev, zdrugs_prev, anxiolytics_prev, benzodiazepines_prev,
                substance_prev, alcohol_prev, htn_prev, diabetes_prev, arrhythmia_prev, angina_prev, liver_prev, renal_prev, dyslipidaemia_prev, 
                gpconsults_prev, priorphyshosp2y_prev, priorpsychhosp2y_prev)

cv_events_cohort <- df_list %>% purrr::reduce(full_join, by = "patid") %>%
  mutate(across(apuse_prior2years:priordyslipidaemia, ~ if_else(is.na(.), "No", if_else(. == 1, "Yes", "No"))), # set any missing concomitant medication values to 'No' and convert to factor
         across(apuse_prior2years:priordyslipidaemia, ~ factor(., levels = c("No", "Yes"))),
         across(gpconsults_last6m:priorpsychhosp_num, ~ if_else(is.na(.), 0, .)),
         perprotocol = case_when(is.na(notadherent_event) ~ 1,
                                 TRUE ~ 0),  # set any missing numeric values to 0
         monthsstononadherence = as.numeric(notadherent_date - cohortentrydate)/30.25)

rm(list = ls(pattern = "_prev"), df_list, hes_episodes, switchdisc_adversereaction, tte_hes_primarydiag, tte_hes_selfharm, self_harm_seplls)
gc()

# Smoking status
load("SMIsmokingstatus_C.Rdata")

smi_smoking <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  left_join(SMIsmokingstatus_C, by = "patid") %>%
  distinct()

#Create hierarchy (current, ex, never, missing) to handle ties
score_category <- function(status) {
  if (sum(!is.na(status)) == 0) {
    return(NA)
  } else if (sum(status == 2 & !is.na(status)) > 0) {
    return("Current")
  } else if (sum(status == 1 & !is.na(status)) > 0) {
    return("Ex")
  } else {
    return("Never")}}

# Look back over whole medical history
# Group by patid and score smoking categories for each patient
smoking_prev <- smi_smoking %>%
  filter(eventdate <= cohortentrydate) %>%
  group_by(patid) %>%
  mutate(ever = ifelse(smoking_status %in% c(1, 2), 1, 0),
         max = max(eventdate)) %>%
  ungroup() %>%
  filter(max == eventdate) %>%
  group_by(patid) %>%
  mutate(smoking_status_cat = score_category(smoking_status)) %>%
  ungroup() %>%
  select(patid, smoking_status_cat) %>%
  distinct()

cv_events_cohort <- merge(cv_events_cohort, smoking_prev, by = "patid", all = TRUE)

# For patients with missing smoking data, check future records to see if they only ever have never smoked codes
missing_smoking <- cv_events_cohort %>%
  filter(is.na(smoking_status_cat)) %>%
  select(patid, smoking_status_cat) %>%
  left_join(smi_smoking, by = "patid") %>%
  select(patid, medcode, readterm, eventdate, smoking_status) %>%
  group_by(patid) %>%
  mutate(smoking_status_cat = score_category(smoking_status)) %>% 
  ungroup() %>%
  group_by(patid) %>%
  summarize(only_never = all(smoking_status_cat == "Never")) %>%
  filter(only_never == TRUE)

#Merge into TTE data frame
cv_events_cohort <- merge(cv_events_cohort, missing_smoking, by = "patid", all = TRUE)

cv_events_cohort$smoking_status_cat[cv_events_cohort$only_never == TRUE] <- "Never"

cv_events_cohort <- cv_events_cohort %>%
  select(-only_never) %>%
  mutate(smoking_status_cat <- factor(smoking_status_cat, 
                                      levels = c(0, 1, 2), labels = c("Never", "Ex", "Current")))

rm(list = ls(pattern = "smoking"))

# Weight
load("SMIweight_C.Rdata")
weight_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIweight_C, by = "patid") %>%
  filter(eventdate >= cohortentrydate - days(731) & eventdate <= cohortentrydate) %>%
  group_by(patid) %>%
  mutate(baseline = max(eventdate)) %>%
  ungroup() %>%
  filter(baseline == eventdate) %>%
  distinct() %>%
  group_by(patid) %>%
  mutate(meanweightkg = mean(weightkg)) %>%
  ungroup() %>%
  select(patid, cohortentrydate, baseline_weightkg = meanweightkg, baseline_weightkg_date = eventdate) %>%
  distinct()

# BMI values recorded directly in CPRD
load("SMIBMI_combined.Rdata")
BMI_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIBMI_C, by = "patid") %>%
  filter(eventdate >= cohortentrydate - days(731) & eventdate <= cohortentrydate) %>%
  group_by(patid) %>%
  mutate(baseline = max(eventdate)) %>%
  ungroup() %>%
  filter(baseline == eventdate) %>%
  distinct() %>%
  group_by(patid) %>%
  mutate(meanbmi = mean(bmi)) %>%
  ungroup() %>%
  select(patid, baseline_bmi = meanbmi, baseline_bmi_date = eventdate) %>%
  distinct()

# Height
load("SMIheight_C.Rdata")
# Choose height from most recent recording (in past or future)
height_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  left_join(SMIheight_C, by = "patid") %>%
  select(patid, cohortentrydate, eventdate, height) %>%
  mutate(diff = eventdate - cohortentrydate,
         diff = abs(diff)) %>%
  group_by(patid) %>%
  filter(abs(diff-5)==min(abs(diff-5))) %>%
  ungroup() %>%
  select(patid, height) %>%
  distinct() %>%
  group_by(patid) %>%
  mutate(meanheight = mean(height)) %>%
  ungroup() %>%
  select(patid, height = meanheight) %>%
  distinct()

# BMI categories recorded directly in CPRD
load("SMIBMIcat_C.Rdata")
BMIcat_prev <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(SMIBMIcat_C, by = "patid") %>%
  filter(eventdate >= cohortentrydate - days(731) & eventdate <= cohortentrydate) %>%
  group_by(patid) %>%
  mutate(baseline = max(eventdate)) %>%
  ungroup() %>%
  filter(baseline == eventdate) %>%
  distinct() %>%
  group_by(patid, eventdate) %>%
  mutate(num = row_number(),
         anyobese = ifelse(BMIcat == "Obese", 1, 0), # handle multiple on the same day
         anyoverweight = ifelse(BMIcat == "Overweight", 1, 0),
         anyhealthy = ifelse(BMIcat == "Healthy", 1, 0),
         anyunderweight = ifelse(BMIcat == "Underweight", 1, 0)) %>%
  ungroup() %>%
  group_by(patid) %>%
  mutate(num_sum = sum(num),
         obese_sum = sum(anyobese),
         overweight_sum = sum(anyoverweight),
         healthy_sum = sum(anyhealthy),
         underweightweight_sum = sum(anyunderweight)) %>%
  ungroup() %>%
  mutate(BMIcat = case_when(num_sum > 1 & obese_sum > 0 ~ "Obese",
                            num_sum > 1 & overweight_sum > 0 & obese_sum == 0 ~ "Overweight",
                            num_sum > 1 & healthy_sum > 0 & obese_sum == 0 & overweight_sum == 0 ~ "Healthy",
                            num_sum > 1 & underweightweight_sum > 0 & obese_sum == 0 & overweight_sum == 0 & healthy_sum == 0 ~ "Underweight",
                            TRUE ~ BMIcat)) %>%
  filter(num == 1) %>%
  select(patid, baseline_bmi_cat = BMIcat, baseline_bmicatdate = eventdate)

# Merge all together
weightandBMI_prev <- weight_prev %>%
  full_join(BMI_prev, by = "patid") %>%
  full_join(height_prev, by = "patid") %>%
  full_join(BMIcat_prev, by = "patid") %>%
  mutate(bmi_calc = round(baseline_weightkg/height^2, 2),
         baseline_bmi = coalesce(baseline_bmi, bmi_calc),
         baseline_bmi_date = coalesce(baseline_bmi_date, baseline_weightkg_date),
         baseline_bmi_date = ifelse(is.na(baseline_bmi), NA, baseline_bmi_date),
         BMIcat_calc = case_when(baseline_bmi < 18.5 ~ "Underweight",
                                 baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Healthy",
                                 baseline_bmi >= 25 & baseline_bmi < 30 ~ "Overweight",
                                 baseline_bmi >= 30 ~ "Obese"),
         baseline_bmi_cat = coalesce(BMIcat_calc, baseline_bmi_cat),
         baseline_bmi_cat = factor(baseline_bmi_cat), # prioritise the category that matches the BMI value being used
         weight_calc = round(baseline_bmi*height^2, 2), # calculate weight from BMI and height
         baseline_weightkg_date = ifelse(!is.na(weight_calc) & is.na(baseline_weightkg), baseline_bmi_date, baseline_weightkg_date),
         baseline_weightkg = coalesce(baseline_weightkg, weight_calc)) %>% # use calculated weight for those that are missing %>%
  select(-weight_calc, -bmi_calc, -BMIcat_calc, -cohortentrydate, -baseline_bmicatdate)

cv_events_cohort <- merge(cv_events_cohort, weightandBMI_prev, by = "patid", all = TRUE)

rm(list = ls(pattern = "weight|BMI|height"))

# Mortality
cv_events_cohort <- cv_events_cohort %>%
  mutate(timetodeath = as.numeric(deathdate - cohortentrydate),
         diedby6m = case_when(timetodeath <= 180 ~ 1,
                              TRUE ~ 0),
         diedby5y = case_when(timetodeath <= 1826 ~ 1,
                              TRUE ~ 0))

# OUTCOMES ####

# Primary stroke or aMI diagnosis, at the episode level
tte_hes_primarydiag <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(hes_primarydiag, by = "patid") %>%
  filter(d_order == 1) %>% # filter to primary diagnosis 
  filter(grepl("(?i)I60|I61|I62|I63|I64|I21", ICD)) %>% # select diagnoses of interest
  group_by(patid, spno) %>%
  mutate(primary_ICD = paste(unique(ICD), collapse = "/")) %>% # paste all primary diagnoses per patient and spellno into one column
  ungroup() %>%
  select(-ICD, -d_order, -epiend, -epikey) %>%
  mutate(event = case_when(grepl("(?i)I21", primary_ICD) & grepl("(?i)I60|I61|I62|I63|I64", primary_ICD) ~ "Stroke/MI",
                           grepl("(?i)I21", primary_ICD) ~ "MI",
                           grepl("(?i)I60|I61|I62|I63|I64", primary_ICD) ~ "Stroke",
                           TRUE ~ "Unknown"))

remove(hes_primarydiag)

# Identify spells that are associated with relevant admissions, used in subsequent merges
spells <- tte_hes_primarydiag %>%
  select(patid, spno) %>%
  distinct()
n_distinct(spells$patid) # 1,047

# Hospital level information
tte_hes_hosp <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  inner_join(spells, by = "patid") %>% # inner join to keep only patients eligible for linkage and those that had relevant spells identified
  left_join(hes_hosp, by = c("patid", "spno")) %>% # left join to retrieve data at the hospital admission level
  select(-discharged, -ICD) %>%
  distinct()

#Full join of all HES data, at hospital spell level
hes_dataset <- tte_hes_primarydiag %>% 
  full_join(tte_hes_hosp, by = c("patid", "spno", "cohortentrydate")) %>% #  link to hospital level information (other diagnoses and admission date)
  distinct() %>%
  rename(event_date = admidate)

# CVD death
death_patient_gold_1 <- load_hes("death_patient_21_000729.txt", hes_config$deaths)
death_patient_gold_2 <- load_hes("death_patient_21_000729_request2.txt", hes_config$deaths) 
death_patient_aurum <- load_hes("death_patient_21_000729.txt", hes_config$deaths)

deaths_ons <- death_patient_aurum %>%
  bind_rows(death_patient_gold_1) %>%
  bind_rows(death_patient_gold_2) %>%
  select(patid, ons_dod = dod, cause, ons_match_rank = match_rank)

# Filter to relevant patients 
tte_deaths <- cv_events_cohort %>%
  select(patid, cohortentrydate) %>%
  left_join(deaths_ons, by = "patid") %>%
  filter(grepl("(?i)I", cause)) %>% # select diagnoses of interest (I= circulatory causes),
  mutate(event = "CV death") %>%
  rename(event_date = ons_dod, primary_ICD = cause) %>%
  select(-ons_match_rank)

remove(death_patient_aurum, death_patient_gold_1, death_patient_gold_2, hes_hosp, spells, tte_hes_hosp, hes_config, deaths_ons)

# MACE
mace_events <- tte_deaths %>%
  bind_rows(hes_dataset) %>%
  filter(event_date < "2020-01-01") %>%
  select(-spno, -epistart) %>%
  distinct()

remove(hes_dataset)

# Recode to CV death if died within 28 days of MI/stroke
mace_events <- mace_events %>%
  left_join(select(cv_events_cohort, patid, deathdate), by = "patid") %>% # merge in CPRD death date (which incorporates CPRD and ONS deaths)
  mutate(event = case_when(event != "CV death" & !is.na(deathdate) & abs(difftime(deathdate, event_date, units = "days")) <= 28 ~ "CV death",
                           TRUE ~ event)) %>% # if died within 28 days of MI/stroke - code to CV death instead
  distinct() %>%
  select(-deathdate)

mace_5y <- mace_events %>%
  filter(event_date >= cohortentrydate & event_date <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  group_by(patid) %>%
  mutate(mace_event = 1,
         mace_event_n = n_distinct(event_date),
         datediff = as.numeric(event_date - cohortentrydate),
         min = min(event_date)) %>%
  ungroup() %>%
  filter(event_date == min) %>% # filter to first only
  rename(mace_daystofirst = datediff, mace_first_date = event_date, mace_icd = primary_ICD) %>%
  group_by(patid, mace_first_date) %>%
  mutate(mace_icd = paste(unique(mace_icd), collapse = "/")) %>%
  ungroup() %>%
  distinct() %>%
  select(-cohortentrydate, -min)

# Some pts died and had stroke/MI on same date, keep only the death record
cv_events_cohort <- cv_events_cohort %>%
  full_join(mace_5y, by = "patid") %>%
  mutate(mace_event = case_when(is.na(mace_event) ~ 0,
                                TRUE ~ mace_event),
         ons_enddate = as.Date("2021-03-29"), # The end date for ONS linkage is 29Mar2021
         hes_enddate = as.Date("2021-03-31"), # HES APC data collection end date was 31Mar21
         study_enddate = as.Date("2019-12-31"), # End of study follow-up date
         
         fivey_date = cohortentrydate + 1826, # 5y follow-up date
         
         fupend_itt_date = pmin(hes_enddate, ons_enddate, study_enddate, deathdate, fivey_date, na.rm = TRUE), # ITT end of follow-up, earliest of death/hes/ons/study end date
         
         mace_itt_event = case_when(mace_first_date <= fupend_itt_date ~ 1, # if mace happened <= fupend, count it as an event
                                    TRUE ~ 0),
         
         mace_itt_fuptime_cen = case_when(mace_itt_event == 1 ~ as.numeric(mace_first_date - cohortentrydate), # if mace happened before fupend, count follow-up time to mace
                                          TRUE ~ as.numeric(fupend_itt_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]
         
         mace_itt_fupend_reason = case_when(mace_itt_event == 1 ~ "Had MACE",
                                            fupend_itt_date == deathdate ~ "Non-cardiovascular death",
                                            fupend_itt_date == study_enddate |  fupend_itt_date == fivey_date ~ "End of follow-up",
                                            fupend_itt_date == hes_enddate | fupend_itt_date == ons_enddate ~ "Administrative censoring",
                                            TRUE ~ "Unknown"),
         
         fupend_pp_date = pmin(hes_enddate, ons_enddate, study_enddate, deathdate, fivey_date, notadherent_date, regenddate, lcd, na.rm = TRUE),
         # PP end of follow-up, earliest of death, regend, lcd, 5y date, discontinue, switch, hes/ons/study end date (this requires continued primary care registration to ascertain if remained PP)
         
         mace_pp_event = case_when(mace_first_date <= fupend_pp_date ~ 1, # if mace happened <= fupend, count it as an event
                                   TRUE ~ 0),
         
         mace_pp_fuptime_cen = case_when(mace_pp_event == 1 ~ as.numeric(mace_first_date - cohortentrydate), # if mace happened before fupend, count follow-up time to mace
                                         TRUE ~ as.numeric(fupend_pp_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]
         
         mace_pp_fupend_reason = case_when(mace_pp_event == 1 ~ "Had MACE",
                                           fupend_pp_date == deathdate ~ "Non-cardiovascular death",
                                           fupend_pp_date == notadherent_date ~ "Switched/Discontinued",
                                           fupend_pp_date == regenddate ~ "Ended registration",
                                           fupend_pp_date == study_enddate | fupend_pp_date == fivey_date ~ "End of follow-up",
                                           fupend_pp_date == lcd | fupend_pp_date == hes_enddate | fupend_pp_date == ons_enddate ~ "Administrative censoring",
                                           TRUE ~ "Unknown"),
         
         firstmace_year = year(mace_first_date),
         firstmace_age = firstmace_year - yob,
         firstmace_age_cat = factor(case_when(firstmace_age < 40 ~ 1, # Age at cohort entry categories
                                              firstmace_age >= 40 & firstmace_age < 50 ~ 2,
                                              firstmace_age >= 50 & firstmace_age < 60 ~ 3,
                                              firstmace_age >= 60 & firstmace_age < 70 ~ 4,
                                              firstmace_age >= 70 & firstmace_age < 80 ~ 5,
                                              firstmace_age >= 80 ~ 6), 
                                    levels = c(1, 2, 3, 4, 5, 6), labels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+")))

stroke_5y <- mace_events %>%
  filter(grepl("(?i)Stroke", event)) %>%
  filter(event_date >= cohortentrydate & event_date <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  group_by(patid) %>%
  mutate(stroke_event = 1,
         stroke_event_n = n_distinct(event_date),
         datediff = as.numeric(event_date - cohortentrydate),
         min = min(event_date)) %>%
  ungroup() %>%
  filter(event_date == min) %>% # filter to first only
  rename(stroke_daystofirst = datediff, stroke_first_date = event_date, stroke_icd = primary_ICD) %>%
  distinct() %>%
  select(-cohortentrydate, -min, -event)

cv_events_cohort <- cv_events_cohort %>%
  full_join(stroke_5y, by = "patid") %>%
  mutate(stroke_event = case_when(is.na(stroke_event) ~ 0,
                                  TRUE ~ stroke_event),

         stroke_itt_event = case_when(stroke_first_date <= fupend_itt_date ~ 1, # if stroke happened <= fupend, count it as an event
                                      TRUE ~ 0),
         
         stroke_itt_fuptime_cen = case_when(stroke_itt_event == 1 ~ as.numeric(stroke_first_date - cohortentrydate), # if stroke happened before fupend, count follow-up time to stroke
                                        TRUE ~ as.numeric(fupend_itt_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]
         
         stroke_itt_fupend_reason = case_when(stroke_itt_event == 1 ~ "Had stroke",
                                              fupend_itt_date == deathdate ~ "Died",
                                              fupend_itt_date == study_enddate | fupend_itt_date == fivey_date ~ "End of follow-up",
                                              fupend_itt_date == hes_enddate | fupend_itt_date == ons_enddate ~ "Administrative censoring",
                                              TRUE ~ "Unknown"),
         
         stroke_pp_event = case_when(stroke_first_date <= fupend_pp_date ~ 1, # if stroke happened <= fupend, count it as an event
                                     TRUE ~ 0),
         
         stroke_pp_fuptime_cen = case_when(stroke_pp_event == 1 ~ as.numeric(stroke_first_date - cohortentrydate), # if stroke happened before fupend, count follow-up time to stroke
                                       TRUE ~ as.numeric(fupend_pp_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]

         stroke_pp_fupend_reason = case_when(stroke_pp_event == 1 ~ "Had stroke",
                                             fupend_pp_date == deathdate ~ "Died",
                                             fupend_pp_date == notadherent_date ~ "Switched/Discontinued",
                                             fupend_pp_date == regenddate ~ "Ended registration",
                                             fupend_pp_date == study_enddate | fupend_pp_date == fivey_date ~ "End of follow-up",
                                             fupend_pp_date == lcd | fupend_pp_date == hes_enddate | fupend_pp_date == ons_enddate ~ "Administrative censoring",
                                             TRUE ~ "Unknown"),
         
         firststroke_year = year(stroke_first_date),
         firststroke_age = firststroke_year - yob,
         firststroke_age_cat = factor(case_when(firststroke_age < 40 ~ 1, # Age at cohort entry categories
                                                firststroke_age >= 40 & firststroke_age < 50 ~ 2,
                                                firststroke_age >= 50 & firststroke_age < 60 ~ 3,
                                                firststroke_age >= 60 & firststroke_age < 70 ~ 4,
                                                firststroke_age >= 70 & firststroke_age < 80 ~ 5,
                                                firststroke_age >= 80 ~ 6), 
                                      levels = c(1, 2, 3, 4, 5, 6), labels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+")))

mi_5y <- mace_events %>%
  filter(grepl("(?i)mi", event)) %>%
  filter(event_date >= cohortentrydate & event_date <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  group_by(patid) %>%
  mutate(mi_event = 1,
         mi_event_n = n_distinct(event_date),
         datediff = as.numeric(event_date - cohortentrydate),
         min = min(event_date)) %>%
  ungroup() %>%
  filter(event_date == min) %>% # filter to first only
  rename(mi_daystofirst = datediff, mi_first_date = event_date, mi_icd = primary_ICD) %>%
  distinct() %>%
  select(-cohortentrydate, -min, -event)

cv_events_cohort <- cv_events_cohort %>%
  full_join(mi_5y, by = "patid") %>%
  mutate(mi_event = case_when(is.na(mi_event) ~ 0,
                              TRUE ~ mi_event),

         mi_itt_event = case_when(mi_first_date <= fupend_itt_date ~ 1, # if mi happened <= fupend, count it as an event
                                  TRUE ~ 0),
         
         mi_itt_fuptime_cen = case_when(mi_itt_event == 1 ~ as.numeric(mi_first_date - cohortentrydate), # if mi happened before fupend, count follow-up time to mi
                                    TRUE ~ as.numeric(fupend_itt_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]
         
         mi_itt_fupend_reason = case_when(mi_itt_event == 1 ~ "Had MI",
                                          fupend_itt_date == deathdate ~ "Died",
                                          fupend_itt_date == study_enddate | fupend_itt_date == fivey_date ~ "End of follow-up",
                                          fupend_itt_date == hes_enddate | fupend_itt_date == ons_enddate ~ "Administrative censoring",
                                          TRUE ~ "Unknown"),
         
         fupend_pp_date = pmin(hes_enddate, ons_enddate, study_enddate, deathdate, notadherent_date, regenddate, lcd, na.rm = TRUE), 
         # PP end of follow-up, earliest of death, regend, lcd, discontinue, switch, hes/ons/study end date (this requires continued primary care registration to ascertain if remained PP)
         
         mi_pp_event = case_when(mi_first_date <= fupend_pp_date ~ 1, # if mi happened <= fupend, count it as an event
                                 TRUE ~ 0),
         
         mi_pp_fuptime_cen = case_when(mi_pp_event == 1 ~ as.numeric(mi_first_date - cohortentrydate), # if mi happened before fupend, count follow-up time to mi
                                   TRUE ~ as.numeric(fupend_pp_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]
         
         mi_pp_fupend_reason = case_when(mi_pp_event == 1 ~ "Had MI",
                                         fupend_pp_date == deathdate ~ "Died",
                                         fupend_pp_date == notadherent_date ~ "Switched/Discontinued",
                                         fupend_pp_date == regenddate ~ "Ended registration",
                                         fupend_pp_date == study_enddate | fupend_pp_date == fivey_date ~ "End of follow-up",
                                         fupend_pp_date == lcd | fupend_pp_date == hes_enddate | fupend_pp_date == ons_enddate ~ "Administrative censoring",
                                         TRUE ~ "Unknown"),
         
         firstmi_year = year(mi_first_date),
         firstmi_age = firstmi_year - yob,
         firstmi_age_cat = factor(case_when(firstmi_age < 40 ~ 1, # Age at cohort entry categories
                                            firstmi_age >= 40 & firstmi_age < 50 ~ 2,
                                            firstmi_age >= 50 & firstmi_age < 60 ~ 3,
                                            firstmi_age >= 60 & firstmi_age < 70 ~ 4,
                                            firstmi_age >= 70 & firstmi_age < 80 ~ 5,
                                            firstmi_age >= 80 ~ 6), 
                                  levels = c(1, 2, 3, 4, 5, 6), labels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+")))

cvdeath_5y <- mace_events %>%
  filter(grepl("(?i)cv", event)) %>%
  filter(event_date >= cohortentrydate & event_date <= cohortentrydate + days(1826)) %>% # filter to records within 5 years of cohort entry
  group_by(patid) %>%
  mutate(cvdeath_event = 1,
         cvdeath_event_n = n_distinct(event_date),
         datediff = as.numeric(event_date - cohortentrydate),
         min = min(event_date)) %>%
  ungroup() %>%
  filter(event_date == min) %>% # filter to first only
  rename(cvdeath_daystofirst = datediff, cvdeath_first_date = event_date, cvdeath_icd = primary_ICD) %>%
  group_by(patid) %>%
  mutate(cvdeath_icd = paste(unique(cvdeath_icd), collapse = "/")) %>%
  ungroup() %>%
  distinct() %>%
  select(-cohortentrydate, -min, -event)

cv_events_cohort <- cv_events_cohort %>%
  full_join(cvdeath_5y, by = "patid") %>%
  mutate(cvdeath_event = case_when(is.na(cvdeath_event) ~ 0,
                                   TRUE ~ cvdeath_event),

         cvdeath_itt_event = case_when(cvdeath_first_date <= fupend_itt_date ~ 1, # if cvdeath happened <= fupend, count it as an event
                                       TRUE ~ 0),
         
         cvdeath_itt_fuptime_cen = case_when(cvdeath_itt_event == 1 ~ as.numeric(cvdeath_first_date - cohortentrydate), # if cvdeath happened before fupend, count follow-up time to cvdeath
                                         TRUE ~ as.numeric(fupend_itt_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]

         cvdeath_itt_fupend_reason = case_when(cvdeath_itt_event == 1 ~ "Had CV death",
                                               fupend_itt_date == deathdate ~ "Non-CV death",
                                               fupend_itt_date == study_enddate | fupend_itt_date == fivey_date ~ "End of follow-up",
                                               fupend_itt_date == hes_enddate | fupend_itt_date == ons_enddate ~ "Administrative censoring",
                                               TRUE ~ "Unknown"),
         
         cvdeath_pp_event = case_when(cvdeath_first_date <= fupend_pp_date ~ 1, # if cvdeath happened <= fupend, count it as an event
                                      TRUE ~ 0),
         
         cvdeath_pp_fuptime = case_when(cvdeath_pp_event == 1 ~ as.numeric(cvdeath_first_date - cohortentrydate), # if cvdeath happened before fupend, count follow-up time to cvdeath
                                        TRUE ~ as.numeric(fupend_pp_date - cohortentrydate)), # otherwise calculate follow-up time [this censors at max 5y, due to inclusion of 5y date]
         
         cvdeath_pp_fuptime_cen = case_when(cvdeath_pp_fuptime > 1825 ~ 1826, # censor at 5y
                                            TRUE ~ cvdeath_pp_fuptime),
         
         cvdeath_pp_fupend_reason = case_when(cvdeath_pp_event == 1 ~ "Had CV death",
                                              fupend_pp_date == deathdate ~ "Non-CV death",
                                              fupend_pp_date == notadherent_date ~ "Switched/Discontinued",
                                              fupend_pp_date == regenddate ~ "Ended registration",
                                              fupend_pp_date == study_enddate | fupend_pp_date == fivey_date ~ "End of follow-up",
                                              fupend_pp_date == lcd | fupend_pp_date == hes_enddate | fupend_pp_date == ons_enddate ~ "Administrative censoring",
                                              TRUE ~ "Unknown"),
         
         firstcvdeath_year = year(cvdeath_first_date),
         firstcvdeath_age = firstcvdeath_year - yob,
         firstcvdeath_age_cat = factor(case_when(firstcvdeath_age < 40 ~ 1, # Age at cohort entry categories
                                                 firstcvdeath_age >= 40 & firstcvdeath_age < 50 ~ 2,
                                                 firstcvdeath_age >= 50 & firstcvdeath_age < 60 ~ 3,
                                                 firstcvdeath_age >= 60 & firstcvdeath_age < 70 ~ 4,
                                                 firstcvdeath_age >= 70 & firstcvdeath_age < 80 ~ 5,
                                                 firstcvdeath_age >= 80 ~ 6), 
                                       levels = c(1, 2, 3, 4, 5, 6), labels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+")))

cv_events_cohort <- cv_events_cohort %>%
  mutate(diag_prev = case_when(diag_prev == "bipolar" ~ "Bipolar disorder",
                               diag_prev == "other psychosis" ~ "Other non-organic psychoses",
                               diag_prev == "schizophrenia" ~ "Schizophrenia",
                               TRUE ~ diag_prev),
         smoking_status_cat = case_when(smoking_status_cat == "Never" ~ "Never smoked",
                                        smoking_status_cat == "Ex" ~ "Ex-smoker",
                                        smoking_status_cat == "Current" ~ "Current smoker",
                                        TRUE ~ smoking_status_cat))

cv_events_cohort$gender <- droplevels(cv_events_cohort$gender, exclude = c("Indeterminate", "Unknown"))

#Save to file
n_distinct(cv_events_cohort$patid) #20,404

# Check for duplicates
dupes <-cv_events_cohort %>%
  group_by(patid) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1)

save(cv_events_cohort, file = "cv_events_cohort.rdata")

remove(cvdeath_5y, mace_events, mace_5y, stroke_5y, mi_5y, tte_deaths, dupes)

# Package versions
# packageVersion("dplyr") # 1.1.4
# packageVersion("tidyr") # 1.3.1
# packageVersion("stringr") # 1.5.1
# packageVersion("forcats") # 1.0.0
# packageVersion("lubridate") # 1.9.4
# packageVersion("readr") # 2.1.5
# packageVersion("readxl") # 1.4.3
# packageVersion("data.table") # 1.16.4
# packageVersion("purrr") # 1.0.2
# packageVersion("rlang") # 1.1.4
# packageVersion("reshape2") # 1.4.4
# packageVersion("doseminer") # 0.1.2
# packageVersion("drugprepr") # 0.0.4
# packageVersion("chlorpromazineR") # 0.2.0
# packageVersion("ggplot2") # 3.5.1
# packageVersion("gtsummary") # 2.0.4
# packageVersion("gt") # 0.11.1
# packageVersion("mice") # 3.17.0
# packageVersion("finalfit") # 1.0.8
# packageVersion("tidylog") # 1.1.0
# packageVersion("survminer") # 0.5.0
# packageVersion("survival") # 3.7.0
# packageVersion("vroom") # 1.6.5
