# ============================================================
# 03_endpoint_analysis.R
# Endpoint/statistical analysis for simulated trials
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(emmeans)
  library(openxlsx)
  library(readr)
  library(pracma)
})

# ------------------------------------------------------------
# Generic helpers
# ------------------------------------------------------------
pick_time <- function(time, val, target, tol = 1e-8) {
  idx <- which(abs(time - target) < tol)
  if (length(idx) == 0) NA_real_ else val[idx[1]]
}

safe_read_trial_file <- function(file_path) {
  df_ctrl <- read_excel(file_path, sheet = "Control_ViralLoad") %>%
    mutate(arm = "Control")
  
  df_treat <- read_excel(file_path, sheet = "Treatment_ViralLoad") %>%
    mutate(arm = "Treatment")
  
  bind_rows(df_ctrl, df_treat) %>%
    mutate(arm = factor(arm, levels = c("Control", "Treatment")))
}

list_simulation_files <- function(sim_dir) {
  list.files(
    sim_dir,
    pattern = "Simulation_run_\\d+\\.xlsx$",
    full.names = TRUE
  )
}

# ============================================================
# RECOVERY-style ANCOVA
# ============================================================

process_one_recovery <- function(file_path, lod = 1.3) {
  
  df <- safe_read_trial_file(file_path)
  
  df_wide <- df %>%
    group_by(ID, sim, arm, sigma, age) %>%
    summarise(
      baseline_logVL = pick_time(time, log10V, sigma),
      day3_logVL     = pick_time(time, log10V, sigma + 2),
      day5_logVL     = pick_time(time, log10V, sigma + 4),
      .groups = "drop"
    ) %>%
    mutate(
      baseline_logVL = ifelse(baseline_logVL < lod, lod, baseline_logVL),
      day3_logVL     = ifelse(day3_logVL < lod, lod, day3_logVL),
      day5_logVL     = ifelse(day5_logVL < lod, lod, day5_logVL)
    )
  
  fit3 <- lm(day3_logVL ~ arm + baseline_logVL + age, data = df_wide)
  em3 <- emmeans(fit3, ~ arm)
  means3 <- summary(em3, type = "response")
  contrast3 <- summary(contrast(em3, "revpairwise"), type = "response")
  
  out3 <- data.frame(
    file = basename(file_path),
    control_D3 = means3$emmean[means3$arm == "Control"],
    treat_D3 = means3$emmean[means3$arm == "Treatment"],
    pval_D3 = contrast3$p.value,
    signif_D3 = contrast3$p.value < 0.05
  )
  
  fit5 <- lm(day5_logVL ~ arm + baseline_logVL + age, data = df_wide)
  em5 <- emmeans(fit5, ~ arm)
  means5 <- summary(em5, type = "response")
  contrast5 <- summary(contrast(em5, "revpairwise"), type = "response")
  
  out5 <- data.frame(
    file = basename(file_path),
    control_D5 = means5$emmean[means5$arm == "Control"],
    treat_D5 = means5$emmean[means5$arm == "Treatment"],
    pval_D5 = contrast5$p.value,
    signif_D5 = contrast5$p.value < 0.05
  )
  
  list(day3 = out3, day5 = out5)
}

run_recovery_analysis <- function(sim_dir,
                                  out_file = file.path(dirname(sim_dir), "ANCOVA_summary_results_LoD.xlsx"),
                                  lod = 1.3) {
  
  files <- list_simulation_files(sim_dir)
  res_list <- lapply(files, process_one_recovery, lod = lod)
  
  all_day3 <- bind_rows(lapply(res_list, function(x) x$day3))
  all_day5 <- bind_rows(lapply(res_list, function(x) x$day5))
  
  wb <- createWorkbook()
  addWorksheet(wb, "Day3")
  addWorksheet(wb, "Day5")
  
  writeData(wb, "Day3", all_day3)
  writeData(wb, "Day5", all_day5)
  
  saveWorkbook(wb, out_file, overwrite = TRUE)
  
  cat("RECOVERY results saved to:", out_file, "\n")
  invisible(list(day3 = all_day3, day5 = all_day5))
}

# ============================================================
# EPIC-HR ANCOVA
# Only onset_to_randomisation <= 3 days
# ============================================================

process_one_epichr <- function(file_path) {
  
  df <- safe_read_trial_file(file_path) %>%
    mutate(
      log10V_imp = ifelse(
        is.na(log10V),
        NA_real_,
        ifelse(log10V < 2, 1.7, log10V)
      )
    )
  
  df_wide <- df %>%
    group_by(ID, sim, arm, sigma, onset_to_randomisation) %>%
    summarise(
      baseline_logVL = pick_time(time, log10V_imp, sigma),
      day5_logVL = pick_time(time, log10V_imp, sigma + 5),
      .groups = "drop"
    ) %>%
    filter(
      !is.na(baseline_logVL),
      !is.na(day5_logVL),
      baseline_logVL > 1.7,
      onset_to_randomisation <= 3
    ) %>%
    mutate(
      change_D5 = day5_logVL - baseline_logVL
    )
  
  if (nrow(df_wide) == 0) return(NULL)
  
  fit <- lm(change_D5 ~ arm + baseline_logVL, data = df_wide)
  
  em <- emmeans(fit, ~ arm)
  em_sum <- summary(em)
  diff <- summary(contrast(em, "revpairwise"))
  
  data.frame(
    file = basename(file_path),
    population = "≤3 days",
    n = nrow(df_wide),
    mean_day5_ctrl = em_sum$emmean[em_sum$arm == "Control"],
    mean_day5_treat = em_sum$emmean[em_sum$arm == "Treatment"],
    diff_treat_vs_ctrl = diff$estimate,
    p_value = diff$p.value,
    significant = diff$p.value < 0.05
  )
}

run_epichr_analysis <- function(sim_dir,
                                out_file = file.path(dirname(sim_dir), "EPIC-HR_summary_results_LOD.xlsx")) {
  
  files <- list_simulation_files(sim_dir)
  
  all_results <- bind_rows(
    lapply(files, process_one_epichr)
  )
  
  write.xlsx(all_results, out_file, overwrite = TRUE)
  
  cat("EPIC-HR results saved to:", out_file, "\n")
  invisible(all_results)
}

# ============================================================
# PLATCOV summary from Stan outputs
# Use q2.5 and q97.5, same as original
# ============================================================

process_platcov_summary_file <- function(file_path) {
  
  summ <- read_csv(file_path, show_col_types = FALSE)
  
  beta_row <- summ %>% filter(variable == "beta_0")
  trt_row  <- summ %>% filter(grepl("trt_effect", variable))
  
  if (nrow(beta_row) == 0 || nrow(trt_row) == 0) {
    return(NULL)
  }
  
  slope_control <- beta_row$mean[1]
  slope_treatment <- beta_row$mean[1] * exp(trt_row$mean[1])
  
  faster_mean <- (exp(trt_row$mean[1]) - 1) * 100
  faster_q2.5 <- (exp(trt_row$q2.5[1]) - 1) * 100
  faster_q97.5 <- (exp(trt_row$q97.5[1]) - 1) * 100
  
  significant <- faster_q2.5 > 0
  
  data.frame(
    file = basename(file_path),
    slope_control = slope_control,
    slope_treatment = slope_treatment,
    faster_mean = faster_mean,
    faster_q2.5 = faster_q2.5,
    faster_q97.5 = faster_q97.5,
    significant = significant
  )
}

run_platcov_summary <- function(stan_out_dir,
                                out_file = file.path(stan_out_dir, "viral_clearance_summary.xlsx")) {
  
  files <- list.files(
    stan_out_dir,
    pattern = "sim[0-9]+_summary\\.csv",
    full.names = TRUE
  )
  
  results <- bind_rows(
    lapply(files, process_platcov_summary_file)
  )
  
  write.xlsx(results, out_file, rowNames = FALSE)
  
  cat("PLATCOV summary saved to:", out_file, "\n")
  invisible(results)
}

# ============================================================
# AUC summary and empirical power
# Power from t-test only
# ============================================================

compute_auc <- function(time, vl) {
  if (length(time) < 2) return(NA_real_)
  
  ord <- order(time)
  time <- time[ord]
  vl <- vl[ord]
  
  trapz(time, vl)
}

compute_patient_auc <- function(df, lod = 1.3) {
  
  df <- df %>%
    mutate(log10V = pmax(log10V, lod))
  
  df %>%
    group_by(ID, sigma) %>%
    group_modify(~ {
      sigma_i <- .y$sigma[[1]]
      
      sub_df <- .x %>%
        filter(time >= sigma_i, time <= sigma_i + 5)
      
      tibble(
        auc = compute_auc(sub_df$time, sub_df$log10V)
      )
    }) %>%
    ungroup() %>%
    filter(!is.na(auc))
}

process_one_auc <- function(file_path, lod = 1.3) {
  
  ctrl <- read_excel(file_path, sheet = "Control_ViralLoad")
  treat <- read_excel(file_path, sheet = "Treatment_ViralLoad")
  
  auc_ctrl <- compute_patient_auc(ctrl, lod = lod)
  auc_treat <- compute_patient_auc(treat, lod = lod)
  
  mean_ctrl <- mean(auc_ctrl$auc, na.rm = TRUE)
  mean_treat <- mean(auc_treat$auc, na.rm = TRUE)
  diff_val <- mean_treat - mean_ctrl
  
  if (nrow(auc_ctrl) > 1 && nrow(auc_treat) > 1) {
    t_res <- t.test(auc_ctrl$auc, auc_treat$auc)
    p_val <- t_res$p.value
  } else {
    p_val <- NA_real_
  }
  
  data.frame(
    file = basename(file_path),
    mean_auc_control = mean_ctrl,
    mean_auc_treatment = mean_treat,
    diff_auc = diff_val,
    p_value = p_val
  )
}

run_auc_power_analysis <- function(sim_dir,
                                   out_file = file.path(sim_dir, "AUC_power_results.xlsx"),
                                   lod = 1.3) {
  
  files <- list_simulation_files(sim_dir)
  
  results <- bind_rows(
    lapply(files, process_one_auc, lod = lod)
  )
  
  power <- mean(results$p_value < 0.05, na.rm = TRUE)
  
  write.xlsx(results, out_file, rowNames = FALSE)
  
  cat("AUC results saved to:", out_file, "\n")
  cat("Estimated power =", power, "\n")
  
  invisible(list(
    results = results,
    power = power
  ))
}
