# ============================================================
# 02_simulation.R
# Trial simulation pipeline
#
# Includes:
# - Covariate generation
# - Treatment/control assignment
# - Individual simulation under PK/PD + viral dynamics
# - Saving outputs into Excel files
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(doParallel)
  library(foreach)
  library(msm)        # rtnorm
  library(truncdist)
})

source(file.path("R", "01_functions.R"))

# ------------------------------------------------------------
# Utility helpers
# ------------------------------------------------------------
sample_age_profile <- function(n, mean_age, sd_age) {
  age <- msm::rtnorm(n = n, mean = mean_age, sd = sd_age)
  age <- round(age, 0)
  age_cat <- ifelse(age < 65, 0, 1)

  age <- ifelse(
    age_cat == 0,
    sample(18:64, n, replace = TRUE),
    sample(65:100, n, replace = TRUE)
  )

  data.frame(age = age, age_cat = age_cat)
}

sample_vaccination_profile <- function(n, vac_prop) {
  n_vac <- round(n * vac_prop)
  sample(c(rep(1, n_vac), rep(0, n - n_vac)))
}

sample_sigma_profile <- function(n,
                                 incubation_shape = 4.09,
                                 incubation_rate = 1.14,
                                 onset_shape = 7.2410,
                                 onset_rate = 2.4463) {
  incubation <- rgamma(n, shape = incubation_shape, rate = incubation_rate)
  onset_to_randomisation <- rgamma(n, shape = onset_shape, rate = onset_rate)
  sigma <- incubation + onset_to_randomisation

  data.frame(
    onset_to_randomisation = onset_to_randomisation,
    sigma = sigma
  )
}

simulate_arm <- function(df, parms, sim_id, antiviral, sd_err = 1.1673) {
  out_list <- lapply(seq_len(nrow(df)), function(i) {
    tryCatch({
      sigma_i <- round(df$sigma[i], 0)

      sim_parms <- simulate_pop_parms2(
        parms = parms,
        age   = df$age_cat[i],
        vac1  = df$vac[i],
        vac2  = 0,
        n     = 1,
        trunc = TRUE
      )

      fit <- fit_model(
        parms        = sim_parms,
        sigma        = sigma_i,
        antiviral    = antiviral,
        avg_efficacy = TRUE,
        verbose      = FALSE
      )

      Vfit <- fit$output[, c("time", "log10V")]
      Vfit$log10V_true <- Vfit$log10V
      Vfit$log10V <- Vfit$log10V + rnorm(nrow(Vfit), 0, sd_err)

      Vfit <- dplyr::mutate(
        Vfit,
        ID = i,
        sim = sim_id,
        age = df$age[i],
        age_cat = df$age_cat[i],
        vac_cat = df$vac[i],
        onset_to_randomisation = df$onset_to_randomisation[i],
        sigma = sigma_i,
        .before = 1
      )

      parm_df <- cbind(
        ID = i,
        sim = sim_id,
        age = df$age[i],
        age_cat = df$age_cat[i],
        vac_cat = df$vac[i],
        onset_to_randomisation = df$onset_to_randomisation[i],
        sigma = sigma_i,
        antiviral = antiviral,
        avg_efficacy = fit$avg_efficacy,
        sim_parms
      )

      list(
        Vfit = Vfit,
        parms = parm_df
      )
    }, error = function(e) {
      NULL
    })
  })

  out_list <- Filter(Negate(is.null), out_list)

  list(
    Vfit  = bind_rows(lapply(out_list, `[[`, "Vfit")),
    parms = bind_rows(lapply(out_list, `[[`, "parms"))
  )
}

save_simulation_workbook <- function(result, sim_id, data_dir) {
  wb <- createWorkbook()

  addWorksheet(wb, "Control_ViralLoad")
  addWorksheet(wb, "Treatment_ViralLoad")
  addWorksheet(wb, "Control_Parms")
  addWorksheet(wb, "Treatment_Parms")

  writeData(wb, "Control_ViralLoad", result$control$Vfit)
  writeData(wb, "Treatment_ViralLoad", result$treat$Vfit)
  writeData(wb, "Control_Parms", result$control$parms)
  writeData(wb, "Treatment_Parms", result$treat$parms)

  saveWorkbook(
    wb,
    file = file.path(data_dir, paste0("Simulation_run_", sprintf("%03d", sim_id), ".xlsx")),
    overwrite = TRUE
  )
}

# ------------------------------------------------------------
# Main trial simulation
# ------------------------------------------------------------
run_simulation <- function(
    n_control = 1126,
    n_treat = 1120,
    vac_prop = 0.2,
    mean_age = 72.5,
    sd_age = 13.9,
    n_sim = 100,
    sd_err = 1.1673,
    parm_file = file.path("data", "populationParameters.txt"),
    output_dir = file.path("results", "simulation_output"),
    seed = 123
) {
  data_dir <- file.path(output_dir, "data")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(parm_file)) {
    stop("Parameter file not found: ", parm_file, call. = FALSE)
  }

  parms <- read.csv(parm_file, row.names = 1)
  n <- n_control + n_treat

  set.seed(seed)

  ncores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr)
      library(openxlsx)
      library(msm)
      library(truncdist)
      library(tidyverse)
      library(deSolve)
      library(MCMCglmm)
      library(truncnorm)
    })
    source(file.path("R", "01_functions.R"))
    NULL
  })

  clusterExport(
    cl,
    varlist = c(
      "sample_age_profile",
      "sample_vaccination_profile",
      "sample_sigma_profile",
      "simulate_arm"
    ),
    envir = environment()
  )

  cat("Using", ncores, "cores for parallel execution.\n")

  results <- foreach(sim = 1:n_sim, .packages = c("dplyr")) %dopar% {
    age_df <- sample_age_profile(n = n, mean_age = mean_age, sd_age = sd_age)
    sigma_df <- sample_sigma_profile(n = n)
    vac <- sample_vaccination_profile(n = n, vac_prop = vac_prop)

    covariates <- data.frame(
      age = age_df$age,
      age_cat = age_df$age_cat,
      vac = vac,
      onset_to_randomisation = sigma_df$onset_to_randomisation,
      sigma = sigma_df$sigma
    )

    treat_idx <- sample(1:n, size = n_treat, replace = FALSE)
    treatment <- covariates[treat_idx, , drop = FALSE]
    control   <- covariates[-treat_idx, , drop = FALSE]

    list(
      control = simulate_arm(control, parms = parms, sim_id = sim, antiviral = FALSE, sd_err = sd_err),
      treat   = simulate_arm(treatment, parms = parms, sim_id = sim, antiviral = TRUE,  sd_err = sd_err)
    )
  }

  stopCluster(cl)

  for (sim in seq_along(results)) {
    save_simulation_workbook(results[[sim]], sim_id = sim, data_dir = data_dir)
  }

  cat("All", n_sim, "simulations finished.\n")
  invisible(results)
}
