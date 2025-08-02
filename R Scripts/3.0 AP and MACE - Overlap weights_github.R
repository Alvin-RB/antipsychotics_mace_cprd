# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 13/11/2024

# Overlap weights

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "readr", "data.table", "purrr", "rlang", "reshape2", 
                   "ggplot2", "gtsummary", "gt", "mice", "finalfit", "tidylog", "survminer", "survival", "vroom", "webshot",
                   "WeightIt", "MatchThem", "grid", "gridExtra", "cobalt", "survey"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Set output directory

# Load data
load("cv_events_cohort.rdata")
load("cv_events_imputed_long.Rdata")

# Using all MI datasets 

# Select vars
cv_events_imputed <- cv_events_imputed_long %>%
  select(.imp, .id,
         patid, trt_group,
         
         # Covariates
         age_atprevcohortentry, ethnicity_cat_cprdhes_4, patprac_2019imd_quintile, gender,
         baselinedose_olanz, priorhosp_psych_any, apuse_prior2years, diag_prev, smidiag_to_cohortentry,
         prioralcoholabuse, priorangina, priorarrhythmia, priordiabetes, priordyslipidaemia, priorhypertension, priorliverdisease, priorrenaldisease, priorsubstanceabuse,
         anticoagulants_prior2years, antidepressant_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years, antiplatelets_prior2years, anxiolytics_prior2years,
         benzodiazepines_prior2years, lipiddrugs_prior2years, moodstab_prior2years, zdrugs_prior2years,
         baseline_bmi_cat, prevcohortentry_year, gpconsults_last6m, priorhosp_physical_any, smoking_status_cat) %>%
  as.mids() # convert to mids

# Generate weights
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

# Assess balance
covariate_balance_love_plot_imputed <- love.plot(imputed_weighted,
                                         binary = "std",
                                         abs = TRUE,
                                         threshold = c(m = .1),
                                         colors = c("darkorange", "darkgreen"),
                                         shapes = c("square", "circle"),
                                         sample.names = c("Unweighted", "Overlap weighted"),
                                         var.names = c(age_atprevcohortentry = "Age",
                                                       'I(age_atprevcohortentry^2)' = "Age (sq)",
                                                       gender = "Sex",
                                                       ethnicity_cat_cprdhes_4 = "Ethnicity",
                                                       diag_prev = "SMI diagnosis",
                                                       baselinedose_olanz = "Starting daily dose",
                                                       'I(baselinedose_olanz^2)' = "Starting daily dose (sq)",
                                                       'I(baselinedose_olanz^3)' = "Starting daily dose (c)",
                                                       smidiag_to_cohortentry = "Time from diagnosis to index",
                                                       'I(smidiag_to_cohortentry^2)' = "Time from diagnosis to index (sq)",
                                                       patprac_2019imd_quintile = "Deprivation quintile",
                                                       prevcohortentry_year = "Index year",
                                                       'I(prevcohortentry_year^2)' = "Index year (sq)",
                                                       gpconsults_last6m = "Number of primary care consults",
                                                       'I(gpconsults_last6m^2)' = "Number of primary care consults (sq)",
                                                       smoking_status_cat = "Smoking status",
                                                       prioralcoholabuse = "Alcohol misuse",
                                                       priorsubstanceabuse = "Substance misuse",
                                                       priordyslipidaemia = "Dyslipidaemia",
                                                       priordiabetes = "Diabetes",
                                                       priorhypertension = "Hypertension",
                                                       priorrenaldisease = "Renal disease",
                                                       priorliverdisease = "Liver disease",
                                                       priorarrhythmia = "Arrhythmia",
                                                       priorangina = "Angina",
                                                       apuse_prior2years = "Prior antipsychotic use",
                                                       lipiddrugs_prior2years = "Lipid-regulating medications",
                                                       hypertensiondrugs_prior2years = "Antihypertensives",
                                                       antidiabetics_prior2years = "Antidiabetics",
                                                       antidepressant_prior2years = "Antidepressants",
                                                       moodstab_prior2years = "Mood stabilisers",
                                                       anticoagulants_prior2years = "Anticoagulants",
                                                       antiplatelets_prior2years = "Antiplatelets",
                                                       zdrugs_prior2years = "Z-drugs",
                                                       anxiolytics_prior2years = "Anxiolytics",
                                                       benzodiazepines_prior2years = "Benzodiazepines",
                                                       baseline_bmi_cat = "Body mass index category",
                                                       priorhosp_physical_any = "Physical health hospitalisation",
                                                       priorhosp_psych_any = "Psychiatric hospitalisation"),
                                         limits = c(0, .70),
                                         position = c(.91, .10),
                                         which.treat = "Aripiprazole") +
  theme(legend.box.background = element_rect(), 
        legend.box.margin = margin(1, 1, 1, 1))

# Standardised differences before/after weighting for variables included in the PS model
absolutediffs_olanz <- as.data.frame(covariate_balance_love_plot_imputed$data) %>%
  select(var, stat, Sample, treat) %>%
  filter(treat == "Olanzapine vs. Aripiprazole") %>%
  select(-treat) %>%
  pivot_wider(names_from = Sample, values_from = stat)

absolutediffs_quet <- as.data.frame(covariate_balance_love_plot_imputed$data) %>%
  select(var, stat, Sample, treat) %>%
  filter(treat == "Quetiapine vs. Aripiprazole") %>%
  select(-treat) %>%
  pivot_wider(names_from = Sample, values_from = stat)

absolutediffs_risp <- as.data.frame(covariate_balance_love_plot_imputed$data) %>%
  select(var, stat, Sample, treat) %>%
  filter(treat == "Risperidone vs. Aripiprazole") %>%
  select(-treat) %>%
  pivot_wider(names_from = Sample, values_from = stat)

absolutediffs <- absolutediffs_olanz %>%
  left_join(absolutediffs_quet, by = "var") %>%
  left_join(absolutediffs_risp, by = "var") %>%
  mutate(var = gsub("_", " - ", var)) %>%
  mutate(var = gsub("- Yes", "", var)) %>%
  arrange(desc(row_number())) %>%  # Reverse the order of rows
  add_row(var = "Comorbidities") %>%
  add_row(var = "Concomitant medications") %>%
  mutate(row_num = row_number(),
         row_num = case_when(
           var == "Age" ~ 1,
           var == "Age (sq)" ~ 2,
             var == "Sex - Male" ~ 3,
             var == "Ethnicity - Asian" ~ 4,
             var == "Ethnicity - Black" ~ 5,
             var == "Ethnicity - Mixed/Other" ~ 6,
             var == "Ethnicity - White" ~ 7,
             var == "SMI diagnosis - Bipolar disorder" ~ 8,
             var == "SMI diagnosis - Other non-organic psychoses" ~ 9,
             var == "SMI diagnosis - Schizophrenia" ~ 10,
             var == "Time from diagnosis to index" ~ 11,
             var == "Time from diagnosis to index (sq)" ~ 12,
             var == "Deprivation quintile - 1 (Least deprived)" ~ 13,
             var == "Deprivation quintile - 2" ~ 14,
             var == "Deprivation quintile - 3" ~ 15,
             var == "Deprivation quintile - 4" ~ 16,
             var == "Deprivation quintile - 5 (Most deprived)" ~ 17,
             var == "Index year" ~ 18,
             var == "Index year (sq)" ~ 19,
             var == "Comorbidities" ~ 20,
             var == "Alcohol misuse " ~ 21,
             var == "Angina " ~ 22,
             var == "Arrhythmia " ~ 23,
             var == "Diabetes " ~ 24,
             var == "Dyslipidaemia " ~ 25,
             var == "Hypertension " ~ 26,
             var == "Liver disease " ~ 27,
             var == "Renal disease " ~ 28,
             var == "Substance misuse " ~ 29,
             var == "Concomitant medications" ~ 30,
             var == "Prior antipsychotic use " ~ 31,
             var == "Anticoagulants " ~ 32,
             var == "Antidepressants " ~ 33,
             var == "Antidiabetics " ~ 34,
             var == "Antihypertensives " ~ 35,
             var == "Antiplatelets " ~ 36,
             var == "Anxiolytics " ~ 37,
             var == "Benzodiazepines " ~ 38,
             var == "Lipid-regulating medications " ~ 39,
             var == "Mood stabilisers " ~ 40,
           var == "Z-drugs " ~ 41,
             var == "Number of primary care consults" ~ 42,
             var == "Number of primary care consults (sq)" ~ 43,
             var == "Smoking status - Current smoker" ~ 44,
             var == "Smoking status - Ex-smoker" ~ 45,
             var == "Smoking status - Never smoked" ~ 46,
             var == "Body mass index category - Underweight" ~ 47,
             var == "Body mass index category - Healthy" ~ 48,
             var == "Body mass index category - Overweight" ~ 49,
             var == "Body mass index category - Obese" ~ 50,
             var == "Psychiatric hospitalisation" ~ 51,
             var == "Physical health hospitalisation" ~ 52,
             var == "Starting daily dose" ~ 53,
             var == "Starting daily dose (sq)" ~ 54,
             var == "Starting daily dose (c)" ~ 55,
             TRUE ~ 999)) %>%
  arrange(row_num) %>%
  mutate(var = case_when(var == "Sex - Male" ~ "Sex", 
                           var == "Body mass index category - Underweight" ~ "BMI <18.5",
                             var == "Body mass index category - Healthy" ~ "BMI ≥18.5 to <25",
                             var == "Body mass index category - Overweight" ~  "BMI ≥25 to <30",
                             var == "Body mass index category - Obese" ~  "BMI ≥30",
                           TRUE ~ var)) %>%
  arrange(row_num) %>%
  select(-row_num) %>%
  gt() %>%
  tab_spanner(label = "Olanzapine vs. Aripiprazole",  # add spanner for N column
              columns = c(Unweighted.x, `Overlap weighted.x`), id = "olanz_ari") %>%
  tab_spanner(label = "Quetiapine vs. Aripiprazole",  # add spanner for N column
              columns = c(Unweighted.y, `Overlap weighted.y`), id = "quet_ari") %>%
  tab_spanner(label = "Risperidone vs. Aripiprazole",  # add spanner for N column
              columns = c(Unweighted, `Overlap weighted`), id = "quet_risp") %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_spanners(spanners = c("olanz_ari", "quet_ari", "quet_risp"))) %>%
  cols_label( # tidy column labels
    var = "Variable",
    Unweighted.x = "Unweighted",
    Unweighted.y = "Unweighted",
    Unweighted = "Unweighted",
    `Overlap weighted.x` = "Overlap weighted",
    `Overlap weighted.y` = "Overlap weighted",
    `Overlap weighted` = "Overlap weighted"
  ) %>%
  fmt_number(columns = c(2:7), decimals = 2) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  cols_align(
    align = "center", # center allign columns
    columns = everything()) %>% # add footnotes
  tab_footnote(paste0("BMI, body mass index. Calculated as the average across all ", (n_distinct(cv_events_imputed_long$.imp)-1), " imputations."))

# Add weights to imputed dataset
# Checked all weights are bounded between 0 and 1
mace_itt_imputed_long <- MatchThem::complete(imputed_weighted, action = "long", include = FALSE, all = TRUE)

# compute effective sample size
compute_ess <- function(weights) {
  round((sum(weights)^2) / sum(weights^2), 0)}

# Calculate effective sample size across imputed datasets
ess_by_group_df <- list()

for (i in unique(mace_itt_imputed_long$.imp[mace_itt_imputed_long$.imp > 0])) {
  imputed_data <- mace_itt_imputed_long %>% filter(.imp == i)
  ess_by_group <- imputed_data %>%
    group_by(trt_group) %>%
    summarise(ESS = compute_ess(weights), .groups = 'drop') %>%
    mutate(.imp = i)
  
  ess_by_group_df[[i]] <- ess_by_group}

# Combine all results into one data frame
ess_by_group_df <- bind_rows(ess_by_group_df) %>%
  spread(key = trt_group, value = ESS)

# Calculate mean across imputations
mean_ess_by_group <- ess_by_group_df %>%
  summarise(across(-.imp, mean, na.rm = TRUE)) %>%
  mutate(across(everything(), round, 0))

# Calculate the total weight for each treatment group
total_weight_by_group_df <- list()

for (i in unique(mace_itt_imputed_long$.imp[mace_itt_imputed_long$.imp > 0])) {
  imputed_data <- mace_itt_imputed_long %>% filter(.imp == i)
  total_weight_by_group <- imputed_data %>%
    group_by(trt_group) %>%
    summarise(total_weight = sum(weights), .groups = 'drop') %>%
    mutate(.imp = i)
  
  total_weight_by_group_df[[i]] <- total_weight_by_group}

# Combine all results into one data frame
total_weight_by_group_df <- bind_rows(total_weight_by_group_df) %>%
  spread(key = trt_group, value = total_weight)

# Calculate mean across imputations
mean_total_weight_by_group <- total_weight_by_group_df %>%
  summarise(across(-.imp, mean, na.rm = TRUE))

# Scale weights to effective sample size
mean_total_weight_long <- mean_total_weight_by_group %>%
  pivot_longer(
    cols = everything(),  
    names_to = "trt_group", 
    values_to = "mean_total_weight")  

mean_ess_by_group_long <- mean_ess_by_group %>%
  pivot_longer(
    cols = everything(), 
    names_to = "trt_group", 
    values_to = "mean_ess") 

extra_from_imputed <- cv_events_imputed_long %>%
  select(patid, .imp, cohortentrydate, baseline_weightkg, baseline_bmi, age_atfirstdiag, 
         mace_itt_event, mace_itt_fuptime_cen, 
         mace_pp_event, mace_pp_fuptime_cen,
         perprotocol, notadherent_date, notadherent_event, monthsstononadherence, mace_pp_fupend_reason) %>%
  mutate(mace_itt_fuptime_cen_months = round(mace_itt_fuptime_cen / 30.25, 0))
  
mace_itt_imputed_long <- mace_itt_imputed_long %>%
  left_join(mean_total_weight_long, by = "trt_group") %>%
  left_join(mean_ess_by_group_long, by = "trt_group") %>%
  left_join(extra_from_imputed, by = c("patid", ".imp")) %>%
  mutate(weights_scale = weights / mean_total_weight * mean_ess)

# Create svy object, one imputed dataset
mace_itt_imputed_long_no_pt <- mace_itt_imputed_long %>%
  filter(.imp == 5) %>%
  select(-patid, -.imp, -.id, -weights, -mean_total_weight, -mean_ess)

# Create the survey design object with the scaled weights 
# (this ensures that the characteristics are scaled by the effective sample size by group)
cv_events_svy <- survey::svydesign(
  ids = ~1,
  strata = ~trt_group,
  weights = ~weights_scale,
  data = mace_itt_imputed_long_no_pt)

# Plot the distribution of w_scale by treatment group
scaled_weight_distribution <- mace_itt_imputed_long %>%
  select(patid, trt_group, weights_scale, .imp) %>%
  group_by(patid) %>%
  mutate(mean = mean(weights_scale)) %>%
  ungroup() %>%
  select(-.imp, -weights_scale) %>%
  distinct() %>%
  ggplot(aes(x = mean, fill = trt_group, color = trt_group)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(NA, 3)) +  # Trim x-axis at 3
  labs(
       x = "Mean",  # X-axis label
       y = "Density",
       fill = "Treatment strategy",
       color = "Treatment strategy") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

# Save to file
n_distinct(mace_itt_imputed_long$patid) # 20,404
save(mace_itt_imputed_long, file = "mace_itt_imputed_long.rdata")

# Package versions
# packageVersion("dplyr") # 1.1.4
# packageVersion("tidyr") # 1.3.1
# packageVersion("stringr") # 1.5.1
# packageVersion("forcats") # 1.0.0
# packageVersion("lubridate") # 1.9.4
# packageVersion("readr") # 2.1.5
# packageVersion("data.table") # 1.16.4
# packageVersion("purrr") # 1.0.2
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
# packageVersion("vroom") # 1.6.5
# packageVersion("webshot") # 0.5.5
# packageVersion("WeightIt") # 1.3.2
# packageVersion("MatchThem") # 1.2.1
# packageVersion("grid") # 4.4.2
# packageVersion("gridExtra") # 2.3
# packageVersion("cobalt") # 4.5.5
# packageVersion("survey") # 4.4.2
