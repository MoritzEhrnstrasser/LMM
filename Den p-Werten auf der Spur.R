#### Replikation für p Werte ####

library(readxl)
library(lme4)
library(Matrix)
library(tmvtnorm)
library(lmerTest)

# Set working directory
setwd("~/Documents/Uni Kram/R")

# Load the data
meta <- read_xlsx("Metascript Multilevel HCCM Repli.xlsx")

# Number of simulations
simulations <- 500

# Simulation function
function_sim <- function(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2) {
  # Level-1 predictors
  predictors <- rbinom(N * Clstr, 1, 0.5)
  
  e <- rnorm(N * Clstr, mean = 0, sd = 1)
  
  level_2_predictors <- rbinom(Clstr, 1, 0.5)
  
  # Combine Level-1 and Level-2 data
  cluster_id <- rep(1:Clstr, each = N)
  level_2_expanded <- rep(level_2_predictors, each = N)
  x <- cbind(1, predictors, level_2_expanded)
  colnames(x) <- c("Intercept", "x1", "x2")
  
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = predictors, x2 = level_2_expanded, cluster_id = cluster_id)
  
  # Design matrix Z
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
  z <- do.call(bdiag, z_blocks)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u_base)
  u_vector <- as.vector(t(u))
  
  # Generate response variable y
  y <- x %*% all_pred + z %*% u_vector + e
  sim_data$y <- as.vector(y)
  
  # Formula for model
  m <- y ~ x1 + x2 + (1 + x1 | cluster_id)
  
  # Fit model using lmer
  fit <- lmer(m, data = sim_data, REML = TRUE)
  
  
  
  # Extract coefficients
  coef_ <- summary(fit)$coefficients
  p_value <- coef_["x1", "Pr(>|t|)"]
  
  return(p_value)
}

# Initialize list to store results
results_list <- list()

# Run the simulations for each row in the metadata
for (i in 1:nrow(meta)) {
  N <- as.numeric(meta[i, "N"]) # Number of observations per cluster
  Clstr <- as.numeric(meta[i, "Clstr"])  # Number of clusters
  lvl1_pred <- as.numeric(meta[i, grep("b", names(meta))])  # Fixed effects coefficients
  sigma <- as.numeric(meta[i, "sigma"])  # Level-1 error variance
  sigma_u_base <- matrix(as.numeric(c(meta[i, "var_u1"], meta[i, "cov_u1_u2"], meta[i, "cov_u1_u2"], meta[i, "var_u2"])), nrow = 2, ncol = 2)  # Random effects covariance matrix
  n_u <- ncol(sigma_u_base)
  lvl2_pred <- as.numeric(meta[i, grep("p", names(meta))])
  N_lvl2 <- length(lvl2_pred)
  all_pred <- c(lvl1_pred, lvl2_pred)
  sigma_lvl2 <- as.numeric(meta[i, "sigma_lvl2"])
  
  # Run simulations
  p_values <- sapply(1:simulations, function(j) function_sim(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2), simplify = "array")
  
  # Compute power
  power <- mean(p_values < 0.05) 
  
  # Store results
  results_list[[i]] <- list(
    condition = i,
    power = power,
    p_values = p_values
  )
}

# Summarize results
power_results <- sapply(results_list, function(res) res$power)
print(power_results)






#### wald KI #### 

library(readxl)
library(lme4)
library(Matrix)
library(tmvtnorm)
library(lmerTest)

# Set working directory
setwd("~/Documents/Uni Kram/R")

# Load the data
meta <- read_xlsx("Metascript Multilevel HCCM Repli.xlsx")

# Number of simulations
simulations <- 1500

# Simulation function
function_sim <- function(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2) {
  # Level-1 predictors
  predictors <- rbinom(N * Clstr, 1, 0.5)
  
  e <- rnorm(N * Clstr, mean = 0, sd = sigma)
  
  level_2_predictors <- rbinom(Clstr, 1, 0.5)
  
  # Combine Level-1 and Level-2 data
  cluster_id <- rep(1:Clstr, each = N)
  level_2_expanded <- rep(level_2_predictors, each = N)
  x <- cbind(1, predictors, level_2_expanded)
  colnames(x) <- c("Intercept", "x1", "x2")
  
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = predictors, x2 = level_2_expanded, cluster_id = cluster_id)
  
  # Design matrix Z
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
  z <- do.call(bdiag, z_blocks)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u_base)
  u_vector <- as.vector(t(u))
  
  # Generate response variable y
  y <- x %*% all_pred + z %*% u_vector + e
  sim_data$y <- as.vector(y)
  
  # Formula for model
  m <- y ~ x1 + x2 + (1 + x1 | cluster_id)
  
  # Fit model using lmer
  fit <- lmer(m, data = sim_data, REML = TRUE)
  
  # Calculate confidence intervals
  wald <- confint(fit, parm = "x1", method = "Wald" , level = 0.95)
  
  # Check if 0 is in the CI and assign P_value
  P_value <- ifelse(wald[1, 1] * wald[1, 2] < 0, 0, 1)
  
  return(P_value)
}

# Initialize list to store results
results_list1 <- list()

# Run the simulations for each row in the metadata
for (i in 1:nrow(meta)) {
  N <- as.numeric(meta[i, "N"]) # Number of observations per cluster
  Clstr <- as.numeric(meta[i, "Clstr"])  # Number of clusters
  lvl1_pred <- as.numeric(meta[i, grep("b", names(meta))])  # Fixed effects coefficients
  sigma <- as.numeric(meta[i, "sigma"])  # Level-1 error variance
  sigma_u_base <- matrix(as.numeric(c(meta[i, "var_u1"], meta[i, "cov_u1_u2"], meta[i, "cov_u1_u2"], meta[i, "var_u2"])), nrow = 2, ncol = 2)  # Random effects covariance matrix
  n_u <- ncol(sigma_u_base)
  lvl2_pred <- as.numeric(meta[i, grep("p", names(meta))])
  N_lvl2 <- length(lvl2_pred)
  all_pred <- c(lvl1_pred, lvl2_pred)
  sigma_lvl2 <- as.numeric(meta[i, "sigma_lvl2"])
  
  # Run simulations
  p_values <- sapply(1:simulations, function(j) {
    function_sim(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2)
  }, simplify = "array")
  
  # Compute power (proportion of simulations with P_value == 1)
  power <- mean(p_values == 1)
  
  # Store results
  results_list1[[i]] <- list(
    condition = i,
    power = power,
    p_values = p_values
  )
}

# Summarize results
power_results1 <- sapply(results_list1, function(res) res$power)
print(power_results1)








#### profile KI #### 

library(readxl)
library(lme4)
library(Matrix)
library(tmvtnorm)
library(lmerTest)

# Set working directory
setwd("~/Documents/Uni Kram/R")

# Load the data
meta <- read_xlsx("Metascript Multilevel HCCM Repli.xlsx")

# Number of simulations
simulations <- 1500

# Simulation function
function_sim <- function(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2) {
  tryCatch({
    # Level-1 predictors
    predictors <- rbinom(N * Clstr, 1, 0.5)
    
    e <- rnorm(N * Clstr, mean = 0, sd = sigma)
    
    level_2_predictors <- rbinom(Clstr, 1, 0.5)
    
    # Combine Level-1 and Level-2 data
    cluster_id <- rep(1:Clstr, each = N)
    level_2_expanded <- rep(level_2_predictors, each = N)
    x <- cbind(1, predictors, level_2_expanded)
    colnames(x) <- c("Intercept", "x1", "x2")
    
    sim_data <- data.frame(y = rep(0, N * Clstr), x1 = predictors, x2 = level_2_expanded, cluster_id = cluster_id)
    
    # Design matrix Z
    z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
    z <- do.call(bdiag, z_blocks)
    
    u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u_base)
    u_vector <- as.vector(t(u))
    
    # Generate response variable y
    y <- x %*% all_pred + z %*% u_vector + e
    sim_data$y <- as.vector(y)
    
    # Formula for model
    m <- y ~ x1 + x2 + (1 + x1 | cluster_id)
    
    # Fit model using lmer
    fit <- lmer(m, data = sim_data, REML = TRUE)
    
    # Calculate confidence intervals
    profile <- confint(fit, parm = "x1", method = "profile", devmatchtol = 0.05)
    
    # Check if 0 is in the CI and assign P_value
    P_value <- ifelse(profile[1, 1] * profile[1, 2] < 0, 0, 1)
    
    return(P_value)
  }, error = function(e) {
    # Return NA if any error occurs
    return(NA)
  })
}

# Initialize list to store results
results_list2 <- list()

# Run the simulations for each row in the metadata
for (i in 1:nrow(meta)) {
  N <- as.numeric(meta[i, "N"]) # Number of observations per cluster
  Clstr <- as.numeric(meta[i, "Clstr"])  # Number of clusters
  lvl1_pred <- as.numeric(meta[i, grep("b", names(meta))])  # Fixed effects coefficients
  sigma <- as.numeric(meta[i, "sigma"])  # Level-1 error variance
  sigma_u_base <- matrix(as.numeric(c(meta[i, "var_u1"], meta[i, "cov_u1_u2"], meta[i, "cov_u1_u2"], meta[i, "var_u2"])), nrow = 2, ncol = 2)  # Random effects covariance matrix
  n_u <- ncol(sigma_u_base)
  lvl2_pred <- as.numeric(meta[i, grep("p", names(meta))])
  N_lvl2 <- length(lvl2_pred)
  all_pred <- c(lvl1_pred, lvl2_pred)
  sigma_lvl2 <- as.numeric(meta[i, "sigma_lvl2"])
  
  # Run simulations
  p_values <- sapply(1:simulations, function(j) {
    function_sim(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2)
  }, simplify = "array")
  
  # Compute power (proportion of simulations with P_value == 1)
  power <- mean(p_values == 1)
  
  # Store results
  results_list2[[i]] <- list(
    condition = i,
    power = power,
    p_values = p_values
  )
}

# Summarize results
power_results2 <- sapply(results_list2, function(res) res$power)
print(power_results2)
