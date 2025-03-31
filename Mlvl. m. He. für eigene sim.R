library(MASS)
library(readxl)
library(Matrix)
library(sandwich)
library(lmtest)
library(tmvtnorm)
library(truncnorm)
library(ggplot2)
library(abind)
library(tidyverse)
library(metRology)
library(brms)
library(posterior)
library(lme4)
library(lmerTest)
library(nlme)
library(CR2)
library(microbenchmark)
library(parallel)
# Set working directory
setwd("~/Documents/Uni Kram/R")

# Load the data
meta <- read_xlsx("Metascript Multilevel HCCM eigene.xlsx")

# Number of simulations
simulations <- 8000

# Simulation function
function_sim <- function(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2) {
  # Level-1 predictors
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, length(lvl1_pred)-1), sigma = sigma)
  
  
  # calculate scaling factor
  # scaling_factor <- 3 * (predictors[, 1])^4
  
 # scaling_factor_normalized <- scaling_factor / sqrt(mean(scaling_factor^2))
  
  # heteroscedastic Residuals
 # e <- rnorm(N * Clstr, mean = 0, sd = scaling_factor_normalized)

  
   #e <- rnorm(N * Clstr, mean = 0, sd = 5 * (predictors[, 1])^8)
  
  
   
  e <- rnorm(N * Clstr, mean = 0, sd = 1 * exp(-0.5 * (predictors[, 1] - mean(predictors[,1]))^2))
  
  
  #e <- rnorm(N* Clstr, mean = 0, sd =  8 * (predictors[, 1] - min(predictors[,1])))
  
  e <- e*(sqrt(1)/sd(e))

  
  level_2_predictors <- matrix(rtmvnorm(Clstr, mean = rep(0, length(lvl2_pred)), sigma = sigma_lvl2), ncol = N_lvl2)
  
  # Combine Level-1 and Level-2 data
  cluster_id <- rep(1:Clstr, each = N)
  level_2_expanded <- do.call(rbind, lapply(1:Clstr, function(i) matrix(rep(level_2_predictors[i, ], each = N), nrow = N)))
  
  x <- cbind(1, predictors, level_2_expanded)
  colnames(x) <- c("Intercept", "x1", "x2","l1","l2")
  
  sim_data <- data.frame(y = rep(0, N * Clstr), cbind(x,cluster_id))
  sim_data <- sim_data[, !names(sim_data) %in% "Intercept"]
  # Design matrix Z
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
  z <- do.call(bdiag, z_blocks)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u_base)
  
  u_vector <- as.vector(t(u))
  
  # Generate response variable y
  y <- x %*% all_pred + z %*% u_vector + e
  sim_data$y <- as.vector(y)
  
  # Formula for model
  m <- y ~ x1 + x2 + l1 + l2 + (1 + x1 | cluster_id)
  
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
  sigma <- matrix(as.numeric(c(meta[i, "var_1"], meta[i, "cov_12"], meta[i, "cov_12"], meta[i, "var_2"])), nrow = 2, ncol = 2)  
  sigma_u_base <- matrix(as.numeric(c(meta[i, "var_u1"], meta[i, "cov_u12"], meta[i, "cov_u12"], meta[i, "var_u2"])), nrow = 2, ncol = 2)  # Random effects covariance matrix
  n_u <- ncol(sigma_u_base)
  lvl2_pred <- as.numeric(meta[i, grep("p", names(meta))])
  N_lvl2 <- length(lvl2_pred)
  all_pred <- c(lvl1_pred, lvl2_pred)
  sigma_lvl2 <- matrix(as.numeric(c(meta[i, "var_2_1"], meta[i, "cov_2_12"], meta[i, "cov_2_12"], meta[i, "var_2_2"])), nrow = 2, ncol = 2) 
 
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







#### Heteroskedastizität ####


e <- rnorm(N * Clstr, mean = 0, sd = 3*(predictors[, 1])^4) # Level 1 "öffnende Trompete"
e <- rnorm(N * Clstr, mean = 0, sd = 0.5 * exp(-0.5 * (predictors[, 1] - mean(predictors[,1]))^2)) # "bauchige" Form
e <- e * (sqrt(1) / sd(e))  # Standardisierung des Fehlers

# Skalierungsfaktor für die Varianz der zufälligen Effekte basierend auf Level-2-Prädiktoren
scaling_factor <- 3*(level_2_predictors[, 1])^4  # Beispiel: exponentielle Skalierung basierend auf Z1

# Generiere zufällige Effekte mit Heteroskedastizität
u <- matrix(0, nrow = Clstr, ncol = n_u)
for (i in 1:Clstr) {
  scaled_sigma_u <- sigma_u_base * scaling_factor[i]
  u[i, ] <- mvrnorm(1, mu = rep(0, n_u), Sigma = scaled_sigma_u)
}


# ohne hetero
e <- rnorm(N * Clstr, mean = 0, sd = 1)
 u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u_base)
 
 
 # mit konstanter Gesamtvarianz
 
 # calculate scaling factor
 scaling_factor <- 3 * (predictors[, 1])^4
 
 scaling_factor_normalized <- scaling_factor / sqrt(mean(scaling_factor^2))
 
 # heteroscedastic Residuals
 e <- rnorm(N * Clstr, mean = 0, sd = scaling_factor_normalized)
 
 # level 2
 
 # Skalierungsfaktor für die Varianz der zufälligen Effekte basierend auf Level-2-Prädiktoren
 scaling_factor <- 3 * (level_2_predictors[, 1])^4  # Beispiel: exponentielle Skalierung basierend auf Z1
 
 # Berechne die durchschnittliche Skalierung (Normierung)
 mean_scaling <- sqrt(mean(scaling_factor^2))  # Mittelwert der quadrierten Skalierungsfaktoren
 scaling_factor_normalized <- scaling_factor / mean_scaling  # Normierter Skalierungsfaktor
 
 # Generiere zufällige Effekte mit Heteroskedastizität
 u <- matrix(0, nrow = Clstr, ncol = n_u)
 for (i in 1:Clstr) {
   scaled_sigma_u <- sigma_u_base * scaling_factor_normalized[i]  # Nutze normalisierten Skalierungsfaktor
   u[i, ] <- mvrnorm(1, mu = rep(0, n_u), Sigma = scaled_sigma_u)
 }
 
 
 
 #### Ergebnisse ####
 
# ohne         0.04971429 0.05342857 0.05371429 0.05485714 0.04914286 0.05028571
 
# nur lvl 1    0.0678  0.0808  0.0894  0.0450  0.0478  0.0428
#              0.03114286 0.04328571 0.04442857 0.04157143 0.04385714 0.04785714
 #             0.051375 0.046750 0.050000 0.047250 0.047875 0.051375   mit konstanter Residualvarianz 
 
# nur lvl 2    0.0354  0.0474  0.0472  0.0336  0.0394  0.0434
 
# lvl 1 und 2  0.0690  0.0746  0.0778  0.0562  0.0476  0.0486

 # 0.04185714 0.04371429 0.04471429 0.03957143 0.05385714 0.04628571 # mit konstaner Residualvarianz
 
 #### Plots ####
 plot(predictors[, 1], e, 
      xlab = "Predictor", ylab = "Residuals")
 
 plot(level_2_predictors[, 1], u[, 1], 
      main = "Level-2 Residuals vs Level-2 Predictor",
      xlab = "Level-2 Predictor (Z1)", 
      ylab = "Level-2 Residuals (u)")
 
 
