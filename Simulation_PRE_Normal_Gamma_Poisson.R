########################################simulation with normal, gamma and Poisson bivariate distributions With PRE  

rm(list = ls())
library(MASS)
library(openxlsx)
library(copula)   
set.seed(123)

# --- Simulation Function ---
simulate_distribution <- function(dist = "normal", N = 100, n = 25, r = 0.887) {
  if (dist == "normal") {
    mu_X <- 44.45
    mu_Y <- 56.47
    var_X <- 3872.573
    var_Y <- 6430.019
    cov_XY <- r * sqrt(var_X * var_Y)
    mu <- c(mu_X, mu_Y)
    Sigma <- matrix(c(var_X, cov_XY,
                      cov_XY, var_Y), nrow = 2, byrow = TRUE)
    pop <- as.data.frame(mvrnorm(n = N, mu = mu, Sigma = Sigma))
    colnames(pop) <- c("X", "Y")
  }
  
  if (dist == "gamma") {
    # set arbitrary shapes/scales
    shape_X <- 5; scale_X <- 10
    shape_Y <- 7; scale_Y <- 8
    cop <- normalCopula(param = r, dim = 2)
    u <- rCopula(N, cop)
    pop <- data.frame(
      X = qgamma(u[,1], shape = shape_X, scale = scale_X),
      Y = qgamma(u[,2], shape = shape_Y, scale = scale_Y)
    )
  }
  
  if (dist == "poisson") {
    lambda_X <- 40
    lambda_Y <- 55
    cop <- normalCopula(param = r, dim = 2)
    u <- rCopula(N, cop)
    pop <- data.frame(
      X = qpois(u[,1], lambda = lambda_X),
      Y = qpois(u[,2], lambda = lambda_Y)
    )
  }
  
  # Add S and Z
  pop$S <- rnorm(N, mean = 0, sd = 2)
  pop$Z <- pop$Y + pop$S
  
  # Sample
  sample_indices <- sample(1:N, n)
  samp <- pop[sample_indices, ]
  
  # --- Parameters ---
  rho_yx <- cor(samp$Y, samp$X)
  rho_zx <- cor(samp$Z, samp$X)
  Y_bar  <- mean(samp$Y)
  X_bar  <- mean(samp$X)
  S_y2   <- var(samp$Y)
  S_x2   <- var(samp$X)
  sigma_s2 <- var(samp$S)
  C_y    <- sd(samp$Y) / Y_bar
  C_x    <- sd(samp$X) / X_bar
  C_z    <- sd(samp$Z) / mean(samp$Z)
  theta  <- X_bar / (X_bar + rho_zx)
  phi    <- (N - n) / (N * n)
  
  param_res <- data.frame(
    rho_yx = rho_yx,
    rho_zx = rho_zx,
    Y_bar = Y_bar,
    X_bar = X_bar,
    S_y2 = S_y2,
    S_x2 = S_x2,
    sigma_s2 = sigma_s2,
    C_y = C_y,
    C_x = C_x,
    C_z = C_z,
    theta = theta,
    phi = phi
  )
  
  
  # --- Estimators ---
  Z_bar   <- mean(samp$Z)
  X_pop   <- mean(pop$X)
  
  Y_hat_mean <- Z_bar
  Y_hat_R1   <- (Z_bar / X_bar) * X_pop
  Y_hat_R2   <- Z_bar * exp((X_pop - X_bar) / (X_pop + X_bar))
  Y_hat_R3   <- Z_bar + rho_zx * (C_y / C_z) * (X_pop - X_bar)
  
  # --- Bias & MSE ---
  bias_mean <- 0
  mse_mean  <- phi * (S_y2 + sigma_s2)
  
  bias_R1 <- phi * Y_bar * (C_x^2 - rho_zx * C_x * C_z)
  mse_R1  <- phi * (Y_bar^2) * (C_z^2 + C_x^2 - 2 * rho_zx * C_z * C_x)
  
  bias_R2 <- phi * Y_bar * ((3/8) * C_x^2 - (1/2))
  mse_R2  <- phi * (Y_bar^2) * (1/4) * (4 * C_z^2 - 4 * rho_zx * C_z * C_x + C_x^2)
  
  theta <- X_bar / (X_bar + rho_zx)
  bias_R3 <- phi * Y_bar * (theta^2 * C_x^2 - 2 * theta * rho_zx * C_z * C_x)
  mse_R3  <- phi * (Y_bar^2) * (theta^2 * C_x^2 + C_z^2 - 2 * theta * rho_zx * C_z * C_x)
  
  estimator_res <- data.frame(
    Estimator = c("Mean", "Ratio_I", "Exp_Ratio_II", "Estimator_III"),
    Bias = c(bias_mean, bias_R1, bias_R2, bias_R3),
    Estimate = c(Y_hat_mean, Y_hat_R1, Y_hat_R2, Y_hat_R3),
    MSE  = c(mse_mean, mse_R1, mse_R2, mse_R3)
  )
  
  # --- Percent Relative Efficiency (PRE) ---
  mse_mean_ref <- estimator_res$MSE[estimator_res$Estimator == "Mean"]
  estimator_res$PRE <- (mse_mean_ref / estimator_res$MSE) * 100
  
  return(list(Population = pop,
              Sample = samp,
              Parameters = param_res,
              Estimators = estimator_res))
}

# --- Run for all 3 distributions ---
results_normal  <- simulate_distribution("normal")
results_gamma   <- simulate_distribution("gamma")
results_poisson <- simulate_distribution("poisson")

# --- Save to Excel ---
wb <- createWorkbook()

# Normal
addWorksheet(wb, "Normal_Population")
addWorksheet(wb, "Normal_Sample")
addWorksheet(wb, "Normal_Parameters")
addWorksheet(wb, "Normal_Estimators")
writeData(wb, "Normal_Population", results_normal$Population)
writeData(wb, "Normal_Sample", results_normal$Sample)
writeData(wb, "Normal_Parameters", results_normal$Parameters)
writeData(wb, "Normal_Estimators", results_normal$Estimators)

# Gamma
addWorksheet(wb, "Gamma_Population")
addWorksheet(wb, "Gamma_Sample")
addWorksheet(wb, "Gamma_Parameters")
addWorksheet(wb, "Gamma_Estimators")
writeData(wb, "Gamma_Population", results_gamma$Population)
writeData(wb, "Gamma_Sample", results_gamma$Sample)
writeData(wb, "Gamma_Parameters", results_gamma$Parameters)
writeData(wb, "Gamma_Estimators", results_gamma$Estimators)

# Poisson
addWorksheet(wb, "Poisson_Population")
addWorksheet(wb, "Poisson_Sample")
addWorksheet(wb, "Poisson_Parameters")
addWorksheet(wb, "Poisson_Estimators")
writeData(wb, "Poisson_Population", results_poisson$Population)
writeData(wb, "Poisson_Sample", results_poisson$Sample)
writeData(wb, "Poisson_Parameters", results_poisson$Parameters)
writeData(wb, "Poisson_Estimators", results_poisson$Estimators)

# Save workbook
saveWorkbook(wb, "Simulation_PRE_Normal_Gamma_Poisson.xlsx", overwrite = TRUE)
