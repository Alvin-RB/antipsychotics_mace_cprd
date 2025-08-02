# -----------------------
# ANTIPSYCHOTICS AND MAJOR ADVERSE CARDIOVASCULAR EVENTS
# -----------------------
# CPRD 2023 data
# Last run: 02/12/2024

# Overlap weights and pooled logistic regression - Intention-to-treat (ITT) analysis

## For the per-protocol analysis, the ITT variables names are replaced with the PP versions. 
## For subgroup and sensitivity analyses, the primary function is adapted (e.g., to iterate over subgroup dataframes, or to use different weights) 

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("tidyr", "data.table", "reshape2", "nnet", "dplyr", "furrr", "future"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Set output directory

# Load data
load("cv_events_imputed_long.Rdata")

cv_events_imputed_long <- cv_events_imputed_long %>%
  filter(.imp > 0) %>%
  dplyr::select(.imp, patid, trt_group,
                age_atprevcohortentry, ethnicity_cat_cprdhes_4, patprac_2019imd_quintile, gender, baselinedose_olanz, priorhosp_psych_any, apuse_prior2years, diag_prev, smidiag_to_cohortentry,
                prioralcoholabuse, priorangina, priorarrhythmia, priordiabetes, priordyslipidaemia, priorhypertension, priorliverdisease, priorrenaldisease, priorsubstanceabuse,
                anticoagulants_prior2years, antidepressant_prior2years, antidiabetics_prior2years, hypertensiondrugs_prior2years, antiplatelets_prior2years,
                anxiolytics_prior2years, benzodiazepines_prior2years, lipiddrugs_prior2years, moodstab_prior2years, zdrugs_prior2years, baseline_bmi_cat, prevcohortentry_year, 
                gpconsults_last6m, priorhosp_physical_any, smoking_status_cat,
                mace_itt_event, mace_itt_fuptime_cen,
                stroke_itt_event, stroke_itt_fuptime_cen,
                mi_itt_event, mi_itt_fuptime_cen,
                cvdeath_itt_event, cvdeath_itt_fuptime_cen)

# convert from days to months, convert 0 values to 1 (as the minimum amount of f-up)
cv_events_imputed_long <- cv_events_imputed_long %>%
  mutate(mace_itt_fuptime_cen = pmax(round(mace_itt_fuptime_cen / 30.25, 0), 1),
         stroke_itt_fuptime_cen = pmax(round(stroke_itt_fuptime_cen / 30.25, 0), 1),
         mi_itt_fuptime_cen = pmax(round(mi_itt_fuptime_cen / 30.25, 0), 1),
         cvdeath_itt_fuptime_cen = pmax(round(cvdeath_itt_fuptime_cen / 30.25, 0), 1))

start <- Sys.time()

log_progress <- function(message) {
  futureCall(function(msg, dir) {
    log_file <- file.path(dir, "progress.log")
    cat(format(Sys.time()), msg, "\n", file = log_file, append = TRUE)
    cat(msg, "\n")
  }, list(message, output_dir))}

# Define function for bootstrap procedure
bootstrap_ow <- function(data, idvar, n.sessions, num.boot, min_imp, max_imp) {
  
  log_progress("MACE ITT primary: Starting analysis")
  
  # Get unique imp values
  imp_range <- seq(min_imp, max_imp)
  imp_values <- unique(imp_range)
  
  # Set up parallel processing
  plan(multisession, workers = n.sessions) 
  
  # Function to process each imp_value
  process_imp_value <- function(imp_value, data, idvar, num.boot) {
    log_progress(paste("MACE ITT primary: Loading imputed dataset #", imp_value))
    
    data0 <- data %>%
      filter(.imp == imp_value)
    
    set.seed(1231 + imp_value)
    seed <- floor(runif(num.boot) * 10^8)

    # Use furrr to parallelize the bootstrap iterations
    bootstrap_results <- future_map(0:num.boot, function(j) {
      log_progress(paste("MACE ITT primary: Starting imputed dataset #", imp_value, "bootstrap #", j))
      
      if (j == 0) {
        
        # calculate estimate in original sample
        boot.data <- copy(data0)
        boot.data$boot.id <- boot.data[[idvar]]
      } else {
        set.seed(seed[j])
        
        # Perform stratified bootstrapping by treatment group
        boot.data <- rbindlist(lapply(unique(data0$trt_group), function(trt) {
          data_sub <- data0[data0$trt_group == trt, ]
          sampled_ids <- sample(unique(data_sub[[idvar]]), length(unique(data_sub[[idvar]])), replace = TRUE)
          boot_sub <- data.table(id = sampled_ids)[, boot.id := .I]
          setnames(boot_sub, "id", idvar)
          merge(boot_sub, data_sub, by = idvar, all.x = TRUE)
        }), use.names = TRUE)}
      
      multinom_model <- multinom(trt_group ~
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
                                        data = boot.data,
                                        trace = FALSE)
      
      gps_matrix <- predict(multinom_model, type = "probs")
      
      # Calculate and normalise OWs
      manual <- data.table(
        patid = boot.data$patid, trt_group = boot.data$trt_group, 
        # Convert gps_matrix to named columns in the data.table
        setNames(as.data.frame(gps_matrix), c("ps_a", "ps_o", "ps_q", "ps_r")))[
          # Calculate inverse probabilities (ip) for each propensity score (ps)
          , c(paste0("ip_", c("a", "o", "q", "r"))) := lapply(.SD, function(x) 1/x), 
          .SDcols = patterns("^ps_")][
            # Calculate the sum of inverse probabilities for each row
            , ip_sum := rowSums(.SD), 
            .SDcols = patterns("^ip_")][
              # Calculate overlap weights (ow) based on treatment group
              , ow_manual := fcase(
                trt_group == "Aripiprazole", ip_a,
                trt_group == "Olanzapine", ip_o,
                trt_group == "Quetiapine", ip_q,
                trt_group == "Risperidone", ip_r) / ip_sum][
                  # Normalize weights within each treatment group
                  , manual_weights_norm := ow_manual / sum(ow_manual), 
                  by = trt_group]
      
      boot.data$w1 <- manual$manual_weights_norm
      
      rm(manual, gps_matrix, multinom_model)
      
      mace_boot.data <- boot.data %>%
        dplyr::select(patid, trt_group, mace_itt_event, mace_itt_fuptime_cen, w1) %>%
        uncount(weights = mace_itt_fuptime_cen, .remove = FALSE) %>%
        group_by(patid) %>% 
        mutate(time = row_number() - 1) %>%  # add variable for time-point
        ungroup() %>%
        mutate(timesq = time^2,  # add time squared
               event = ifelse(time == mace_itt_fuptime_cen - 1 & mace_itt_event == 1, 1, 0)) %>% # create new event indicator, and update its value so that it relates to the specific time-point
        dplyr::select(-mace_itt_event, -mace_itt_fuptime_cen)
      
      mace_itt_model_product <- glm(event ~ trt_group + time + timesq +
                                           trt_group * time +
                                           trt_group * timesq,
                                         family = binomial(link = 'logit'),
                                         weights = w1,
                                         data = mace_boot.data,
                                    control = glm.control(maxit = 50))
      
      mace_convergence_info <- data.table(
        imp = imp_value, iteration = j,
        converged = mace_itt_model_product$converged,
        iterations = mace_itt_model_product$iter,
        deviance = mace_itt_model_product$deviance)
      
      rm(mace_boot.data)
      gc()
      
      stroke_boot.data <- boot.data %>%
        dplyr::select(boot.id, patid, trt_group, stroke_itt_event, stroke_itt_fuptime_cen, w1) %>%
        uncount(weights = stroke_itt_fuptime_cen, .remove = FALSE) %>%
        group_by(patid) %>%
        mutate(time = row_number() - 1) %>%  # add variable for time-point
        ungroup() %>%
        mutate(timesq = time^2,  # add time squared
               event = ifelse(time == stroke_itt_fuptime_cen - 1 & stroke_itt_event == 1, 1, 0)) # create new event indicator, and update its value so that it relates to the specific time-point
      
      stroke_itt_model_product <- glm(event ~ trt_group + time + timesq +
                                             trt_group * time +
                                             trt_group * timesq,
                                           family = binomial(link = 'logit'),
                                           weights = w1,
                                           data = stroke_boot.data,
                                      control = glm.control(maxit = 50))

      stroke_convergence_info <- data.table(
        imp = imp_value, iteration = j,
        converged = stroke_itt_model_product$converged,
        iterations = stroke_itt_model_product$iter,
        deviance = stroke_itt_model_product$deviance)
      
      rm(stroke_boot.data)
      gc()
      
      mi_boot.data <- boot.data %>%
        dplyr::select(boot.id, patid, trt_group, mi_itt_event, mi_itt_fuptime_cen, w1) %>%
        uncount(weights = mi_itt_fuptime_cen, .remove = FALSE) %>%
        group_by(patid) %>%
        mutate(time = row_number() - 1) %>%  # add variable for time-point
        ungroup() %>%
        mutate(timesq = time^2,  # add time squared
               event = ifelse(time == mi_itt_fuptime_cen - 1 & mi_itt_event == 1, 1, 0)) # create new event indicator, and update its value so that it relates to the specific time-point
      
      mi_itt_model_product <- glm(event ~ trt_group + time + timesq +
                                         trt_group * time +
                                         trt_group * timesq,
                                       family = binomial(link = 'logit'),
                                       weights = w1,
                                       data = mi_boot.data,
                                  control = glm.control(maxit = 50))

      mi_convergence_info <- data.table(
        imp = imp_value, iteration = j,
        converged = mi_itt_model_product$converged,
        iterations = mi_itt_model_product$iter,
        deviance = mi_itt_model_product$deviance)
      
      rm(mi_boot.data)
      gc()
      
      cvdeath_boot.data <- boot.data %>%
        dplyr::select(boot.id, patid, trt_group, cvdeath_itt_event, cvdeath_itt_fuptime_cen, w1) %>%
        uncount(weights = cvdeath_itt_fuptime_cen, .remove = FALSE) %>%
        group_by(patid) %>%
        mutate(time = row_number() - 1) %>%  # add variable for time-point
        ungroup() %>%
        mutate(timesq = time^2,  # add time squared
               event = ifelse(time == cvdeath_itt_fuptime_cen - 1 & cvdeath_itt_event == 1, 1, 0)) # create new event indicator, and update its value so that it relates to the specific time-point
      
      cvdeath_itt_model_product <- glm(event ~ trt_group + time + timesq +
                                              trt_group * time +
                                              trt_group * timesq,
                                            family = binomial(link = 'logit'),
                                            weights = w1,
                                            data = cvdeath_boot.data,
                                       control = glm.control(maxit = 50))

      cvdeath_convergence_info <- data.table(
        imp = imp_value, iteration = j,
        converged = cvdeath_itt_model_product$converged,
        iterations = cvdeath_itt_model_product$iter,
        deviance = cvdeath_itt_model_product$deviance)
      
      rm(cvdeath_boot.data)
      gc()
      
      # Create data frame with treatment group, time, and timesq - for each outcome
      for (i in 0:3) {
        df <- data.frame(trt_group = c("Aripiprazole", "Olanzapine", "Quetiapine", "Risperidone")[i + 1], 
                         time = seq(0, 60), 
                         timesq = (seq(0, 60))^2)
        
        hazard_col_name <- paste0("hazard", i)  # Create hazard column name
        df[[hazard_col_name]] <- predict(mace_itt_model_product, df, type = "response")
        
        risk_col_name <- paste0("risk", i)  # Create risk column name
        df[[risk_col_name]] <- 1 - cumprod(1 - df[[hazard_col_name]])
        
        assign(paste0("mace_results", i), df)}
      
      for (i in 0:3) {
        df <- data.frame(trt_group = c("Aripiprazole", "Olanzapine", "Quetiapine", "Risperidone")[i + 1], 
                         time = seq(0, 60), 
                         timesq = (seq(0, 60))^2)
        
        hazard_col_name <- paste0("hazard", i)  # Create hazard column name
        df[[hazard_col_name]] <- predict(stroke_itt_model_product, df, type = "response")
        
        risk_col_name <- paste0("risk", i)  # Create risk column name
        df[[risk_col_name]] <- 1 - cumprod(1 - df[[hazard_col_name]])
        
        assign(paste0("stroke_results", i), df)}
      
      for (i in 0:3) {
        df <- data.frame(trt_group = c("Aripiprazole", "Olanzapine", "Quetiapine", "Risperidone")[i + 1], 
                         time = seq(0, 60), 
                         timesq = (seq(0, 60))^2)
        
        hazard_col_name <- paste0("hazard", i)  # Create hazard column name
        df[[hazard_col_name]] <- predict(mi_itt_model_product, df, type = "response")
        
        risk_col_name <- paste0("risk", i)  # Create risk column name
        df[[risk_col_name]] <- 1 - cumprod(1 - df[[hazard_col_name]])
        
        assign(paste0("mi_results", i), df)}
      
      for (i in 0:3) {
        df <- data.frame(trt_group = c("Aripiprazole", "Olanzapine", "Quetiapine", "Risperidone")[i + 1], 
                         time = seq(0, 60), 
                         timesq = (seq(0, 60))^2)
        
        hazard_col_name <- paste0("hazard", i)  # Create hazard column name
        df[[hazard_col_name]] <- predict(cvdeath_itt_model_product, df, type = "response")
        
        risk_col_name <- paste0("risk", i)  # Create risk column name
        df[[risk_col_name]] <- 1 - cumprod(1 - df[[hazard_col_name]])
        
        assign(paste0("cvdeath_results", i), df)}

      # Create df of MACE results
      mace_results <- data.frame(imp = imp_value, iteration = j)
      
      for (e in 0:60) {
        for (r in 0:3) {
          col_name <- sprintf("risk%s_%sm", r, e)
          mace_results[, col_name] <- get(sprintf("mace_results%s", r))[get(sprintf("mace_results%s", r))$time == e, sprintf("risk%s", r)]
        }}
      
      # Create df of Stroke results
      stroke_results <- data.frame(imp = imp_value, iteration = j)
      
      for (e in 0:60) {
        for (r in 0:3) {
          col_name <- sprintf("risk%s_%sm", r, e)
          stroke_results[, col_name] <- get(sprintf("stroke_results%s", r))[get(sprintf("stroke_results%s", r))$time == e, sprintf("risk%s", r)]
        }}
      
      # Create df of MI results
      mi_results <- data.frame(imp = imp_value, iteration = j)
      
      for (e in 0:60) {
        for (r in 0:3) {
          col_name <- sprintf("risk%s_%sm", r, e)
          mi_results[, col_name] <- get(sprintf("mi_results%s", r))[get(sprintf("mi_results%s", r))$time == e, sprintf("risk%s", r)]
        }}
      
      # Create df of cvdeath results
      cvdeath_results <- data.frame(imp = imp_value, iteration = j)
      
      for (e in 0:60) {
        for (r in 0:3) {
          col_name <- sprintf("risk%s_%sm", r, e)
          cvdeath_results[, col_name] <- get(sprintf("cvdeath_results%s", r))[get(sprintf("cvdeath_results%s", r))$time == e, sprintf("risk%s", r)]
        }}

      list(mace_results = mace_results, 
           stroke_results = stroke_results, 
           mi_results = mi_results,
           cvdeath_results = cvdeath_results,
           mace_convergence_info = mace_convergence_info, 
           stroke_convergence_info = stroke_convergence_info, 
           mi_convergence_info = mi_convergence_info, 
           cvdeath_convergence_info = cvdeath_convergence_info)}, .options = furrr_options(seed = TRUE))
    
    # Combine results
    bootstrap.results <- do.call(rbind, lapply(bootstrap_results, function(x) x$mace_results))
    stroke_bootstrap.results <- do.call(rbind, lapply(bootstrap_results, function(x) x$stroke_results))
    mi_bootstrap.results <- do.call(rbind, lapply(bootstrap_results, function(x) x$mi_results))
    cvdeath_bootstrap.results <- do.call(rbind, lapply(bootstrap_results, function(x) x$cvdeath_results))
    
    mace_convergence_log <- do.call(rbind, lapply(bootstrap_results, function(x) x$mace_convergence_info))
    stroke_convergence_log <- do.call(rbind, lapply(bootstrap_results, function(x) x$stroke_convergence_info))
    mi_convergence_log <- do.call(rbind, lapply(bootstrap_results, function(x) x$mi_convergence_info))
    cvdeath_convergence_log <- do.call(rbind, lapply(bootstrap_results, function(x) x$cvdeath_convergence_info))
    
    list(bootstrap.results = bootstrap.results,
         stroke_bootstrap.results = stroke_bootstrap.results,
         mi_bootstrap.results = mi_bootstrap.results,
         cvdeath_bootstrap.results = cvdeath_bootstrap.results,
         
         mace_convergence_log = mace_convergence_log,
         stroke_convergence_log = stroke_convergence_log,
         mi_convergence_log = mi_convergence_log,
         cvdeath_convergence_log = cvdeath_convergence_log)}
  
  # Use furrr to process each imp_value in parallel
  results_list <- future_map(imp_values, function(imp_value) {
    result <- process_imp_value(imp_value, data, idvar, num.boot)
    
    log_progress(paste("MACE ITT primary: Completed processing imputed dataset #", imp_value))  # New line
    
    result}, .options = furrr_options(seed = TRUE))
  
  # Combine results from all imp_values
  mace_results_list_all <- lapply(results_list, function(x) x$bootstrap.results)
  stroke_results_list_all <- lapply(results_list, function(x) x$stroke_bootstrap.results)
  mi_results_list_all <- lapply(results_list, function(x) x$mi_bootstrap.results)
  cvdeath_results_list_all <- lapply(results_list, function(x) x$cvdeath_bootstrap.results)
  
  mace_convergence_log_all <- lapply(results_list, function(x) x$mace_convergence_log)
  stroke_convergence_log_all <- lapply(results_list, function(x) x$stroke_convergence_log)
  mi_convergence_log_all <- lapply(results_list, function(x) x$mi_convergence_log)
  cvdeath_convergence_log_all <- lapply(results_list, function(x) x$cvdeath_convergence_log)
  
  bootstrap_results_mi <- list(
    mace = mace_results_list_all,
    stroke = stroke_results_list_all,
    mi = mi_results_list_all,
    cvdeath = cvdeath_results_list_all,
    
    convergence_mace = mace_convergence_log_all,
    convergence_stroke = stroke_convergence_log_all,
    convergence_mi = mi_convergence_log_all,
    convergence_cvdeath = cvdeath_convergence_log_all)
  
  # Define the primary and backup save locations
  primary_save_location <- 
  backup_save_location <- 

  log_progress("MACE ITT primary: Completed bootstrap analysis")
  
  return(bootstrap_results_mi)}

bootstrap_results_mi <- bootstrap_ow(data = cv_events_imputed_long, idvar = "patid", n.sessions = 25, num.boot = 500, min_imp = 1, max_imp = 25)

end <- Sys.time()
end - start
rm(bootstrap_results_mi)

# ANALYSIS ####

# Clear memory
rm(list = ls())
gc()

# Packages
invisible(sapply(c("tidyr", "data.table", "rlang", "reshape2", "gtsummary", "gt", "lubridate", "nnet", 
                   "dplyr", "forestplot", "ggplot2", "ggrepel"), 
                 library, character.only = TRUE))

# Set file path

# Set working directory

# Set output directory

# Load data
load("bootstrap_results_itt_mi_all.Rdata")


# Function to process bootstrap results with optional time-point
process_bootstrap_results <- function(bootstrap_results, time_point = 59) {
  combined_results <- do.call(rbind, bootstrap_results)
  
  # Create dynamic column names based on the time_point
  risk_cols <- paste0("risk", 0:3, "_", time_point, "m")
  
  combined_results <- combined_results %>%
    dplyr::select(imp, iteration, all_of(risk_cols)) %>%
    rename_with(~ paste0("risk", 0:3), all_of(risk_cols)) %>%
    mutate(r0r1_ratio = risk0 / risk1, r0r2_ratio = risk0 / risk2, r0r3_ratio = risk0 / risk3,
           r0r1_dif = risk0 - risk1, r0r2_dif = risk0 - risk2, r0r3_dif = risk0 - risk3)
  
  return(combined_results)}

# Process results
mace <- process_bootstrap_results(bootstrap_results_mi$mace)
mace_6m <- process_bootstrap_results(bootstrap_results_mi$mace, time_point = 6)
stroke <- process_bootstrap_results(bootstrap_results_mi$stroke)
mi <- process_bootstrap_results(bootstrap_results_mi$mi)
cvdeath <- process_bootstrap_results(bootstrap_results_mi$cvdeath)

# Compute point estimates and confidence intervals, with optional standardisation
process_risks <- function(data, risk_vars, standardise = FALSE) {
  original_data <- data[data$iteration == 0, ]
  resampled_data <- data[data$iteration > 0, ]
  
  # Apply standardisation, if required
  if (standardise) {
    original_data[risk_vars] <- lapply(original_data[risk_vars], function(x) x * 1000)
    resampled_data[risk_vars] <- lapply(resampled_data[risk_vars], function(x) x * 1000)
  }
  
  # Pre-calculate common values
  alpha <- 0.05
  z_alpha <- qnorm(c(alpha / 2, 1 - alpha / 2))
  
  results_df <- do.call(rbind, lapply(risk_vars, function(var) {
    
    # Calculate point estimate as the mean of the original sample across all imputations
    point_estimate <- mean(original_data[[var]], na.rm = TRUE)
    
    # Extract the bootstrap estimates for the variable
    boot_estimates <- resampled_data[[var]]
    mean_boot <- mean(boot_estimates, na.rm = TRUE)
    sd_boot <- sd(boot_estimates, na.rm = TRUE)
    
    # Calculate the absolute distances between bootstrapped estimates and point estimate
    distance <- abs(boot_estimates - point_estimate)
    
    # Initialize a list to store CIs
    ci_list <- list()
    
    ### METHODS USING THE FULL DISTRIBUTION OF BOOTSTRAP ESTIMATES ACROSS ALL IMPUTATIONS
    
    numerator <- sum((boot_estimates - mean_boot) ^ 3)
    denominator <- sum((boot_estimates - mean_boot) ^ 2)
    a <- if (denominator > 0) numerator / (6 * (denominator ^ 1.5)) else NA
  
    # BCa method - empirical distribution-based
    # Bias is determined by comparing each bootstrap estimate to the mean of bootstrap estimates (i.e. the empirical distribution)
    z0_empirical <- qnorm(mean(boot_estimates < mean_boot))
    p_adj_empirical <- pnorm(z0_empirical + (z0_empirical + z_alpha) / (1 - a * (z0_empirical + z_alpha)))
    ci_list$empirical <- quantile(boot_estimates, probs = p_adj_empirical, na.rm = TRUE)
    
    # Combine results into a data frame
    data.frame(
      variable = var,
      drug = variable_map[var],
      point_estimate = point_estimate,  # Using the pooled estimate from Rubin's Rules
      lower_ci_empirical = ci_list$empirical[1],
      upper_ci_empirical = ci_list$empirical[2])
  }))
  
  return(results_df)
}

# Define risk variables for both timeframes
risk_vars <- c("risk0", "risk1", "risk2", "risk3")
risk_vars_ratio <- c("r0r1_ratio", "r0r2_ratio", "r0r3_ratio")
risk_vars_dif <- c("r0r1_dif", "r0r2_dif", "r0r3_dif")

# Mapping of variables to drug names
variable_map <- c("risk0" = "Aripiprazole", "risk1" = "Olanzapine", "risk2" = "Quetiapine", "risk3" = "Risperidone",
                  "r0r1_ratio" = "Aripiprazole vs. Olanzapine", "r0r2_ratio" = "Aripiprazole vs. Quetiapine", "r0r3_ratio" = "Aripiprazole vs. Risperidone",
                  "r0r1_dif" = "Aripiprazole vs. Olanzapine", "r0r2_dif" = "Aripiprazole vs. Quetiapine", "r0r3_dif" = "Aripiprazole vs. Risperidone")

# Process the risk variables using the defined function

# select CIs (percentile CIs, combined using rubins rules)
select_cis <- function(data) {
 data %>%
    dplyr::select(1:3, lower_ci = lower_ci_empirical, upper_ci = upper_ci_empirical)}

results_df_5y <- process_risks(mace, risk_vars, standardise = TRUE) %>%
  select_cis()

results_df_ratio_5y <- process_risks(mace, risk_vars_ratio) %>%
  select_cis()

results_df_dif_5y <- process_risks(mace, risk_vars_dif, standardise = TRUE) %>%
  select_cis()

stroke_results_df_5y <- process_risks(stroke, risk_vars, standardise = TRUE) %>%
  select_cis()

stroke_results_df_ratio_5y <- process_risks(stroke, risk_vars_ratio) %>%
  select_cis()

stroke_results_df_dif_5y <- process_risks(stroke, risk_vars_dif, standardise = TRUE) %>%
  select_cis()

mi_results_df_5y <- process_risks(mi, risk_vars, standardise = TRUE) %>%
  select_cis()

mi_results_df_ratio_5y <- process_risks(mi, risk_vars_ratio) %>%
  select_cis()

mi_results_df_dif_5y <- process_risks(mi, risk_vars_dif, standardise = TRUE) %>%
  select_cis()

cvdeath_results_df_5y <- process_risks(cvdeath, risk_vars, standardise = TRUE) %>%
  select_cis()

cvdeath_results_df_ratio_5y <- process_risks(cvdeath, risk_vars_ratio) %>%
  select_cis()

cvdeath_results_df_dif_5y <- process_risks(cvdeath, risk_vars_dif, standardise = TRUE) %>%
  select_cis()

results_df_6m <- process_risks(mace_6m, risk_vars, standardise = TRUE) %>%
  select_cis()

results_df_ratio_6m <- process_risks(mace_6m, risk_vars_ratio) %>%
  select_cis()

results_df_dif_6m <- process_risks(mace_6m, risk_vars_dif, standardise = TRUE) %>%
  select_cis()

# Combine all results into one data frame
combined_results <- rbind(
  cbind(results_df_5y, Timepoint = "5y"),
  cbind(results_df_ratio_5y, Timepoint = "5y"),
  cbind(results_df_dif_5y, Timepoint = "5y"),
  cbind(results_df_6m, Timepoint = "6m"),
  cbind(results_df_ratio_6m, Timepoint = "6m"),
  cbind(results_df_dif_6m, Timepoint = "6m"))

stroke_combined_results <- rbind(
  cbind(stroke_results_df_5y, Timepoint = "5y"),
  cbind(stroke_results_df_ratio_5y, Timepoint = "5y"),
  cbind(stroke_results_df_dif_5y, Timepoint = "5y"))

mi_combined_results <- rbind(
  cbind(mi_results_df_5y, Timepoint = "5y"),
  cbind(mi_results_df_ratio_5y, Timepoint = "5y"),
  cbind(mi_results_df_dif_5y, Timepoint = "5y"))

cvdeath_combined_results <- rbind(
  cbind(cvdeath_results_df_5y, Timepoint = "5y"),
  cbind(cvdeath_results_df_ratio_5y, Timepoint = "5y"),
  cbind(cvdeath_results_df_dif_5y, Timepoint = "5y"))

# Reshape the data for presentation
reshape_data <- function(data, risk_vars, type = "RR") {
  
  if (type == "AR") {
    # If type is AR, use this logic
    result <- data %>%
      filter(variable %in% c(risk_vars)) %>%
      mutate(drug = factor(drug, levels = c("Aripiprazole", "Olanzapine", "Quetiapine", "Risperidone")),
             estimate = sprintf("%.1f (%.1f, %.1f)", round(point_estimate, 1), round(lower_ci, 1), round(upper_ci, 1))) %>%
      dplyr::select(Timepoint, drug, estimate) %>%
      spread(key = drug, value = estimate)
    
  } else if (type == "RD") {
    # Round to 1 decimal place for RD
    result <- data %>%
      filter(variable %in% c(risk_vars)) %>%
      mutate(estimate = sprintf("%.1f (%.1f, %.1f)", round(point_estimate, 1), round(lower_ci, 1), round(upper_ci, 1))) %>%
      dplyr::select(Timepoint, drug, estimate) %>%
      spread(key = drug, value = estimate)
    
  } else if (type == "RR") {
    # Round to 2 decimal places for RR
    result <- data %>%
      filter(variable %in% c(risk_vars)) %>%
      mutate(estimate = sprintf("%.2f (%.2f, %.2f)", round(point_estimate, 2), round(lower_ci, 2), round(upper_ci, 2))) %>%
      dplyr::select(Timepoint, drug, estimate) %>%
      spread(key = drug, value = estimate)
  }
  
  return(result)}

mace_absolute_risks <- reshape_data(combined_results, risk_vars, type = "AR") %>%
  mutate(Outcome = "MACE")

stroke_absolute_risks <- reshape_data(stroke_combined_results, risk_vars, type = "AR") %>%
  mutate(Outcome = "Non-fatal stroke")

mi_absolute_risks <- reshape_data(mi_combined_results, risk_vars, type = "AR") %>%
  mutate(Outcome = "Non-fatal myocardial infarction")

cvdeath_absolute_risks <- reshape_data(cvdeath_combined_results, risk_vars, type = "AR") %>%
  mutate(Outcome = "Cardiovascular death")

mace_risk_ratios <- reshape_data(combined_results, risk_vars_ratio, type = "RR")
stroke_risk_ratios <- reshape_data(stroke_combined_results, risk_vars_ratio, type = "RR")
mi_risk_ratios <- reshape_data(mi_combined_results, risk_vars_ratio, type = "RR")
cvdeath_risk_ratios <- reshape_data(cvdeath_combined_results, risk_vars_ratio, type = "RR")

mace_risk_differences <- reshape_data(combined_results, risk_vars_dif, type = "RD")
stroke_risk_differences <- reshape_data(stroke_combined_results, risk_vars_dif, type = "RD")
mi_risk_differences <- reshape_data(mi_combined_results, risk_vars_dif, type = "RD")
cvdeath_risk_differences <- reshape_data(cvdeath_combined_results, risk_vars_dif, type = "RD")

# Define a function to create gt tables
create_gt_table <- function(data, title, subtitle) {
  gt_table <- gt(data) %>%
    tab_header(
      title = title,
      subtitle = subtitle
    ) %>%
    tab_footnote(footnote = paste0("Analysis using ", max_imp_value, " imputed datasets with bias-corrected and accelerated percentile-based 95% confidence intervals generated using a stratified bootstrap procedure with ", 
                                   max(mace$iteration), " samples with replacement in each imputed dataset."))

  return(gt_table)}

# Convert reshaped data into gt tables
absolute_risks <- mace_absolute_risks %>%
  bind_rows(stroke_absolute_risks) %>%
  bind_rows(mi_absolute_risks) %>%
  bind_rows(cvdeath_absolute_risks) %>%
  dplyr::select(Outcome, everything())

absolute_risks_table <- create_gt_table(absolute_risks, "Absolute Risks (per 1,000)", "Bootstrap Results")


# PLOTS ####

# Set working directory

load("bootstrap_results_itt_mi_all.Rdata")

invisible(sapply(c("tidyr", "data.table", "reshape2", "patchwork", "grid", "dplyr", "ggplot2", "scales", "lubridate", "svglite"), 
                 library, character.only = TRUE))

# Function to prepare and merge results (used for both ITT and PP)
# Combined function to calculate point estimates and BCa confidence intervals
generate_cumulative_plots_data <- function(results_data, alpha = 0.05) {
  # Helper function to calculate BCa interval for a single variable
  calculate_bca <- function(x, alpha) {
    n <- length(x)
    boot_mean <- mean(x, na.rm = TRUE)
    
    # Calculate bias correction factor
    z0 <- qnorm(sum(x < boot_mean, na.rm = TRUE) / n)
    
    # Calculate acceleration factor
    numerator <- sum((x - boot_mean) ^ 3)
    denominator <- sum((x - boot_mean) ^ 2)
    a <- if (denominator > 0) numerator / (6 * (denominator ^ 1.5)) else NA
    
    # Calculate adjusted alpha values
    z_alpha <- qnorm(c(alpha / 2, 1 - alpha / 2))
    
    # Calculate adjusted percentiles
    p_adj <- pnorm(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
    
    # Return lower and upper CI as a vector
    quantile(x, probs = p_adj, na.rm = TRUE)}
  
  # Combine results data into a single data frame
  data <- do.call(rbind, results_data)
  
  # Calculate point estimates
  point_estimates <- data %>%
    filter(iteration == 0) %>%
    summarise(across(-c(imp, iteration), mean, na.rm = TRUE)) %>%
    pivot_longer(
      cols = starts_with("risk"),
      names_to = c("treatment_group", "time_point"),
      names_pattern = "risk(\\d{1,2})_(\\d+)m",
      values_to = "point_estimate") %>%
    mutate(time_point = as.numeric(time_point),
           treatment_group = as.numeric(treatment_group))
  
  # Calculate BCa confidence intervals
  ci_results <- data %>%
    select(matches("^risk\\d+_\\d+m$")) %>%
    summarise(across(everything(),
                     ~ {
                       ci <- calculate_bca(.x, alpha)
                       data.frame(lower = ci[1], upper = ci[2])
                     },
                     .names = "{col}_ci")) %>%
    unnest_wider(everything(), names_sep = "_", names_repair = "minimal") %>%
    pivot_longer(cols = everything(),
                 names_to = c("variable", ".value"),
                 names_sep = "_ci") %>%
    rename(ci_lower = `_lower`, ci_upper = `_upper`) %>%
    mutate(variable = gsub("risk", "", variable)) %>%
    separate(variable, into = c("treatment_group", "time_point"), 
             sep = "_", convert = TRUE) %>%
    mutate(time_point = gsub("m", "", time_point),
           time_point = as.numeric(time_point),
           treatment_group = as.numeric(treatment_group))
  
  # Combine point estimates and confidence intervals
  results <- point_estimates %>%
    left_join(ci_results, by = c("treatment_group", "time_point")) %>%
    mutate(time_updated = time_point + 1) %>%
    add_row(time_updated = 0, point_estimate = 0, ci_lower = 0, ci_upper = 0, treatment_group = 0) %>% 
    add_row(time_updated = 0, point_estimate = 0, ci_lower = 0, ci_upper = 0, treatment_group = 1) %>% 
    add_row(time_updated = 0, point_estimate = 0, ci_lower = 0, ci_upper = 0, treatment_group = 2) %>% 
    add_row(time_updated = 0, point_estimate = 0, ci_lower = 0, ci_upper = 0, treatment_group = 3) %>% 
    arrange(time_updated) %>%
    mutate(trt_group = case_when(treatment_group == 0 ~ "Aripiprazole",
                                 treatment_group == 1 ~ "Olanzapine",
                                 treatment_group == 2 ~ "Quetiapine",
                                 treatment_group == 3 ~ "Risperidone",
                                 TRUE ~ "Unknown"))
  
  return(results)}

# Process ITT data
itt_all <- generate_cumulative_plots_data(bootstrap_results_mi$mace)

# Function to create a plot
create_plot <- function(treatment_group_name, data) {
  ggplot(data %>% filter(trt_group %in% c(treatment_group_name, "Aripiprazole")),
         aes(x = time_updated, y = point_estimate, color = trt_group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = trt_group), alpha = 0.2) +
    labs(title = paste("Aripiprazole vs. ", treatment_group_name),
         x = "Time from index date (years)",
         y = "Cumulative risk",
         color = "Treatment Group",
         fill = "Treatment Group") +
    scale_x_continuous(
      breaks = seq(0, 61, 12),
      labels = seq(0, 5, 1),
      limits = c(0, 61)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, .065)) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())}

# Create combined plots
plot_group <- function(data, suffix) {
  plot1 <- create_plot("Olanzapine", data)
  plot2 <- create_plot("Quetiapine", data) + theme(axis.title.y = element_blank())
  plot3 <- create_plot("Risperidone", data) + theme(axis.title.y = element_blank())
  
  combined_plot <- (plot1 | plot2 | plot3) +
    plot_layout(ncol = 3, widths = c(1, 1.1, 1)) +
    plot_annotation(theme = theme(plot.title = element_text(size = 18, hjust = 0.5)))
  
  ggsave(plot = combined_plot,
         path = output_dir, 
         filename = paste0("cumulative_risks_mace_", suffix, "_", today(), ".png"),
         dpi = 300, width = 20, height = 8.5, bg = "white")
  
  ggsave(plot = combined_plot,
         path = output_dir, 
         filename = paste0("cumulative_risks_mace_", suffix, "_", today(), ".pdf"),
         dpi = 300, width = 20, height = 8.5, bg = "white")

  ggsave(plot = combined_plot,
         path = output_dir,
         filename = paste0("cumulative_risks_mace_", suffix, "_", today(), ".svg"),
         width = 20, height = 8.5, bg = "white")}

# Generate and save plots for ITT
plot_group(itt_all, "itt")

# Package versions
packageVersion("tidyr") # 1.3.1
packageVersion("data.table") # 1.16.4
packageVersion("reshape2") # 1.4.4
packageVersion("nnet") # 7.3.19
packageVersion("dplyr") # 1.1.4
packageVersion("furrr") # 0.3.1
packageVersion("future") # 1.34.0
packageVersion("gt") # 0.11.1
packageVersion("gtsummary") # 2.0.4
packageVersion("lubridate") # 1.9.4
packageVersion("ggplot2") # 3.5.1