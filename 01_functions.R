# ============================================================
# 01_functions.R
# Core functions for Digital Twin simulation
#
# Includes:
# - Population parameter extraction from Monolix output
# - Random sampling helper functions
# - PK/PD + viral dynamics model
# - Average efficacy calculation
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(deSolve)
  library(MCMCglmm)   # rtnorm
  library(truncnorm)
})

# ------------------------------------------------------------
# Extract population parameters
# ------------------------------------------------------------
get_pop_parms2 <- function(parms, age = 1, vac1 = 0, vac2 = 0, log = FALSE) {
  parm_names <- c(
    "log10beta_pop", "log10pi_pop", "log10phi_pop",
    "log10rho_pop", "delta_pop", "tau_pop", "IC50_pop"
  )

  if ((vac1 + vac2) >= 2) {
    stop("vac1 and vac2 cannot both be 1.", call. = FALSE)
  }

  if (!all(rownames(parms)[1:7] == parm_names)) {
    stop("Parameters have missing names or incorrect format.", call. = FALSE)
  }

  log10beta <- as.numeric(parms["log10beta_pop", "value"])
  log10pi   <- as.numeric(parms["log10pi_pop", "value"])
  log10phi  <- as.numeric(parms["log10phi_pop", "value"])
  log10rho  <- as.numeric(parms["log10rho_pop", "value"])
  delta     <- as.numeric(parms["delta_pop", "value"])
  tau       <- as.numeric(parms["tau_pop", "value"])
  IC50      <- as.numeric(parms["IC50_pop", "value"])

  if (log) {
    pop_parms <- data.frame(
      log10beta = log10beta,
      log10pi   = log10pi,
      log10phi  = log10phi,
      log10rho  = log10rho,
      delta     = 1.87 + age * (-0.08) + vac1 * 0.10 + vac2 * 0.10,
      tau       = tau,
      IC50      = 3.85 + age * 0.65 + vac1 * (-0.79) + vac2 * (-1.00)
    )
  } else {
    pop_parms <- data.frame(
      beta  = 10^log10beta,
      pi    = 10^log10pi,
      phi   = 10^log10phi,
      rho   = 10^log10rho,
      delta = 1.87 + age * (-0.08) + vac1 * 0.10 + vac2 * 0.10,
      tau   = tau,
      IC50  = 3.85 + age * 0.65 + vac1 * (-0.79) + vac2 * (-1.00)
    )
  }

  return(pop_parms)
}

# ------------------------------------------------------------
# Sampling helper functions
# ------------------------------------------------------------
normal_sample <- function(n = 1000, mean = 0, sd = 1, trunc = FALSE) {
  if (trunc) {
    return(
      rtnorm(
        n = n,
        mean = mean,
        sd = sd,
        lower = mean - (1.96 * sd),
        upper = mean + (1.96 * sd)
      )
    )
  } else {
    return(rnorm(n = n, mean = mean, sd = sd))
  }
}

lognormal_sample <- function(n = 1000, mean = 0, sd = 1, trunc = FALSE) {
  if (trunc) {
    return(
      exp(
        rtnorm(
          n = n,
          mean = log(mean),
          sd = sd,
          lower = log(mean) - (1.96 * sd),
          upper = log(mean) + (1.96 * sd)
        )
      )
    )
  } else {
    return(exp(rnorm(n = n, mean = log(mean), sd = sd)))
  }
}

logitnormal_sample <- function(n = 1000, mean = 0, sd = 1, min = 0, max = 1) {
  logit <- function(x, min = 0, max = 1) {
    log((x - min) / (max - x))
  }

  inverse_logit <- function(x, min = 0, max = 1) {
    min + ((max - min) / (1 + exp(-x)))
  }

  logit_mean <- logit(x = mean, min = min, max = max)
  norm_sample <- rnorm(n = n, mean = logit_mean, sd = sd)
  logit_sample <- inverse_logit(x = norm_sample, min = min, max = max)

  return(logit_sample)
}

# ------------------------------------------------------------
# Simulate individual-level parameters
# ------------------------------------------------------------
simulate_pop_parms2 <- function(parms, age = 1, vac1 = 0, vac2 = 0, n = 1000, trunc = TRUE) {
  pop_parms <- get_pop_parms2(parms, age = age, vac1 = vac1, vac2 = vac2, log = TRUE)

  epsilon_log10beta <- 0.2484
  epsilon_log10pi   <- 0.1175
  epsilon_log10phi  <- 0.01038
  epsilon_log10rho  <- 0.1399
  epsilon_delta     <- 0.2177
  epsilon_tau       <- 1.567
  epsilon_IC50      <- 1.627

  simulated_parms <- data.frame(
    beta  = 10^(normal_sample(n = n, mean = pop_parms$log10beta, sd = epsilon_log10beta, trunc = trunc)),
    pi    = 10^(normal_sample(n = n, mean = pop_parms$log10pi,   sd = epsilon_log10pi,   trunc = trunc)),
    phi   = 10^(normal_sample(n = n, mean = pop_parms$log10phi,  sd = epsilon_log10phi,  trunc = trunc)),
    rho   = 10^(normal_sample(n = n, mean = pop_parms$log10rho,  sd = epsilon_log10rho,  trunc = trunc)),
    delta = normal_sample(n = n, mean = pop_parms$delta, sd = epsilon_delta, trunc = trunc),
    tau   = logitnormal_sample(n = n, mean = pop_parms$tau, sd = epsilon_tau, min = 0, max = 20),
    IC50  = normal_sample(n = n, mean = pop_parms$IC50, sd = epsilon_IC50, trunc = trunc)
  )

  return(simulated_parms)
}

# ------------------------------------------------------------
# Compute average efficacy over treatment window
# ------------------------------------------------------------
average_efficacy <- function(out, t_start = 0, t_end = 5) {
  treatment_window <- out %>%
    as.data.frame() %>%
    dplyr::filter(time >= t_start, time <= t_end)

  if (nrow(treatment_window) < 2) {
    stop("Insufficient data points in the treatment window.", call. = FALSE)
  }

  times <- round(as.numeric(treatment_window$time), digits = 1)
  epsilon_t <- as.numeric(treatment_window$epsilon_t)

  integral <- 0
  for (i in 1:(length(times) - 1)) {
    dt <- times[i + 1] - times[i]
    integral <- integral + 0.5 * (epsilon_t[i] + epsilon_t[i + 1]) * dt
  }

  E_ave <- (1 / (t_end - t_start)) * integral
  return(E_ave)
}

# ------------------------------------------------------------
# PK/PD + viral dynamics model
# ------------------------------------------------------------
fit_model <- function(parms,
                      sigma = 0,
                      antiviral = TRUE,
                      avg_efficacy = TRUE,
                      Tmax = 30,
                      Tstep = 0.1,
                      verbose = FALSE) {

  beta  <- as.numeric(parms$beta)
  pi    <- as.numeric(parms$pi)
  phi   <- as.numeric(parms$phi)
  rho   <- as.numeric(parms$rho)
  delta <- as.numeric(parms$delta)
  tau   <- as.numeric(parms$tau)
  IC50  <- as.numeric(parms$IC50)

  # Fixed parameters
  k      <- 4
  gamma  <- 15
  k_a    <- 9.98
  k_PL   <- 1.58
  k_LP   <- 1.22
  k_CL   <- 4.96
  dose   <- 300
  Vol    <- 41743
  M      <- 499.5
  E_max  <- 0.999
  hill_n <- 3.16

  # Initial conditions
  y <- c(
    T   = 10^7,
    R   = 0,
    E   = 0,
    I   = 0,
    V   = 1,
    A_G = 0,
    A_P = 0,
    A_L = 0
  )

  times <- seq(from = 0, to = Tmax, by = Tstep)
  sigma <- as.numeric(sigma)

  dosing_times <- seq(from = sigma, to = sigma + 4.5, by = 0.5)

  dosing_events <- data.frame(
    var    = "A_G",
    time   = dosing_times,
    value  = dose * (antiviral == TRUE),
    method = "add"
  )

  treatment <- ifelse((times >= sigma) & (times < sigma + 5), as.logical(antiviral), 0)
  H_t <- approxfun(x = times, y = treatment, rule = 2)

  VL_PKPD <- function(times, y, parms) {
    with(as.list(c(parms, y)), {

      # PK model
      dA_G <- -(k_a * A_G)
      dA_P <- (k_a * A_G) + (k_LP * A_L) - (k_PL + k_CL) * A_P
      dA_L <- (k_PL * A_P) - (k_LP * A_L)

      C_P_mg_per_mL <- A_P / Vol
      C_P_uM <- (C_P_mg_per_mL / M) * 1e6

      epsilon <- ((C_P_uM^hill_n) / ((IC50^hill_n) + C_P_uM^hill_n)) * E_max
      epsilon_t <- ifelse(H_t(times) == 1, epsilon, 0)

      # Viral dynamics model
      dT <- -(beta * T * V) - (phi * I * T) + (rho * R)
      dR <- (phi * I * T) - (rho * R)
      dE <- (beta * T * V) - (k * E)
      dI <- (k * E) - (delta * I)
      dV <- ((1 - epsilon_t) * pi * I) - (gamma * V)

      return(list(c(dT, dR, dE, dI, dV, dA_G, dA_P, dA_L), epsilon_t = epsilon_t))
    })
  }

  out <- as.data.frame(
    ode(
      y = y,
      times = times,
      func = VL_PKPD,
      parms = parms,
      method = "bdf",
      events = list(data = dosing_events)
    )
  )

  out <- out %>%
    mutate(
      log10T = ifelse(T > 0, log10(T), NA_real_),
      log10R = ifelse(R > 0, log10(R), NA_real_),
      log10E = ifelse(E > 0, log10(E), NA_real_),
      log10I = ifelse(I > 0, log10(I), NA_real_),
      log10V = ifelse(V > 0, log10(V), NA_real_),
      .after = "time"
    ) %>%
    select(
      time,
      T, R, E, I, V,
      log10T, log10R, log10E, log10I, log10V,
      A_G, A_P, A_L, epsilon_t
    )

  if (verbose) {
    cat("---- PK/PD DEBUG ----\n")
    cat("IC50 =", IC50, "\n")
    cat("Max plasma conc (uM) =", max(out$A_P / Vol / M * 1e6, na.rm = TRUE), "\n")
    cat("Max epsilon_t =", max(out$epsilon_t, na.rm = TRUE), "\n")
    cat("----------------------\n\n")
  }

  if (avg_efficacy) {
    E_avg <- average_efficacy(out = out, t_start = sigma, t_end = sigma + 5)

    results <- list(
      output = out,
      t_dose = if (antiviral) sigma else NA_real_,
      avg_efficacy = if (antiviral) E_avg else NA_real_
    )
    return(results)
  } else {
    return(out)
  }
}
