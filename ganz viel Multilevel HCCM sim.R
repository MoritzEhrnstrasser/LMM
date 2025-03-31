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
library(MASS)
rm(list = ls())

### zweiter Versuch ####
N <- 50
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 80

sigma <- matrix(c(2,0,
                  0,2), nrow = 2, ncol = 2 )

sigma_u <- matrix(c(2,0,0,
                    0,2,0,
                    0,0,2), nrow = 3, ncol = 3 )



function_sim<- function() { 
  
  predictors <- rtmvnorm(N*Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  e<- rnorm(N*Clstr, mean = 0, sd = 1)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = 
                  sigma_u)
  
  x <- cbind(1, predictors)
  
  
  # Designmatrix Z als Blockmatrix erstellen
  
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  y <- x %*% fixed + z %*% as.vector(u) + e
  
  y <- as.vector(y)
  
  cluster_id <- rep(1:Clstr, each = N) #cluster ID hinzufügen
  
  # Dataframe erstellen
  sim_data <- data.frame(y = y, x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  
  # Fitte das Modell mit lmer
  fit <- lmer(y ~ x1 + x2 + (1 + x1 + x2 | cluster_id), data = sim_data)
  
  sum <- summary(fit)
  coef<- sum$coefficients
  varcorr<- VarCorr(fit)
  
  
  return(varcorr)
  
  }

results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# apply um Estimates zu mitteln
apply(results, c(1, 2), mean)

# funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)  






#### erster Versuch ####
# Settings

simulations <- 200
N <- 50
Clstr <- 80
c00 <- 1
c10 <- 1
n_pred_u <- 2

sigma <- matrix(c(2,0,
                  0,2), nrow = 2, ncol = 2 )
dat <- data.frame()
function_sim <- function() {
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_pred_u), sigma = 
                           sigma)
  
  u0 <- u[,1]
  u1 <- u[,2]
  
  x<- rnorm(N*Clstr, mean = 0, sd = 1)
  e<- rnorm(N*Clstr, mean = 0, sd = 1)
  
  B0 <- c00 + u0
  B1 <- c10 + u1
  
  B1 <- matrix(rep(B1, each = N), ncol = 1)
  B0 <- matrix(rep(B0, each = N), ncol = 1)
    
  
  # Berechnung von y 
  y <- B0 + B1 * x + e
  
  dat<- as.data.frame(cbind(y,x,e,B0,B1))
  
  cluster_id <- rep(1:Clstr, each = N)
  
  # Hinzufügen der Gruppierungsvariable zum DataFrame
  dat <- as.data.frame(cbind(y, x, e, B0, B1, cluster_id))
  
  m <- y ~  x + ( 1+ x | cluster_id)
  
  fit<-lmer(m, data=dat, REML = T)  
  summary<- summary(fit)
  varcorr<- VarCorr(fit)
 
  return(varcorr)
}

results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")


attr(results[[1]], "correlation")
attr(results[[1]], "stddev")
c(results[[1]])

allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)  

lapply(results [1:2], mean)



results$cluster_id




# nächste Aufgaben

# (2) Wir wollen flexibel eine beliebige Kombination von Level-1 und Level-2 Prädiktoren simulieren können 
# inkl. Interaktion (Tipp: model.matrix())
# (3) ist ChatGPTs Indizierungsvorschlag schneller als eine logische Abfrage nach dem Clusterindex
# (4) Wie baut man Heteroskedastizität auf Level 1 und Level 2 ein?
#


### benchmarken ####

z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
z <- bdiag(z_blocks)

y <- x %*% fixed + z %*% as.vector(u) + e


microbenchmark({
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z1 <- bdiag(z_blocks)
}, times = 1)

# noch schneller
microbenchmark({z2 <- sparseMatrix(i = rep(1:(N * Clstr), each = ncol(x)),
                                  j = rep(1:(Clstr * ncol(x)), times = N),
                                  x = as.vector(t(x)),
                                  dims = c(N * Clstr, Clstr * ncol(x)))}, times = 1)




all.equal(z1,z1)

#### ####







#### ####

library(MASS)
library(lme4)
library(Matrix)
library(tmvtnorm)

N <- 50
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 80
N_lvl2 <- 3

sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

sigma_u <- matrix(c(2, 0, 0,
                    0, 2, 0,
                    0, 0, 2), nrow = 3, ncol = 3)

# Funktion zur Simulation
function_sim <- function() { 
  
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[,1])^2)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u)
  
  x <- cbind(1, predictors)
  
  # Beispiel für Level-2 Prädiktoren
  level_2_predictors <- matrix(rnorm(Clstr * N_lvl2), ncol = N_lvl2)
  colnames(level_2_predictors) <- c("Z1", "Z2","Z3")
  
  # Cluster-ID hinzufügen und Level-2-Prädiktoren replizieren
  cluster_id <- rep(1:Clstr, each = N)
  level_2_data <- data.frame(cluster_id = 1:Clstr, level_2_predictors)
  level_2_expanded <- level_2_data[rep(1:Clstr, each = N), ]
  
  # Kombiniere Level-1- und Level-2-Daten
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  sim_data <- cbind(sim_data, level_2_expanded[, -1]) # Entferne die doppelte Spalte cluster_id
  
  # Formel inklusive Interaktionen
  m <- y ~ x1 + x2 + Z1 * Z2 + (1 + x1 + x2 | cluster_id)
  
  # Designmatrix Z als Blockmatrix erstellen
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  y <- x %*% fixed + z %*% as.vector(u) + e
  y <- as.vector(y)
  sim_data$y <- y
  
  # Fitte das Modell mit lmer
  
  fit <- lmer(m, data = sim_data)
  
  sum <- summary(fit)
  coef <- sum$coefficients
  varcorr <- VarCorr(fit)
  
  return(varcorr)
}

# Simulationen durchführen
results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# Funktion um Estimates zu mitteln
apply(results, c(1, 2), mean)

# Funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)




#### mit Heteroskedastizität ####
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
library(MASS)
library(Matrix)


N <- 90
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 120
N_lvl2 <- 3

sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

sigma_u_base <- matrix(c(2, 0, 0,
                         0, 2, 0,
                         0, 0, 2), nrow = 3, ncol = 3)

# Funktion zur Simulation
function_sim <- function() { 
  
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[,1])^2)
  
  # Beispiel für Level-2 Prädiktoren
  level_2_predictors <- matrix(rtmvnorm(Clstr , mean = rep(0, N_lvl2), sigma = 
                                          sigma_u_base), ncol = N_lvl2)
  
  colnames(level_2_predictors) <- c("Z1", "Z2", "Z3")
  
  # Skalierungsfaktor für die Varianz der zufälligen Effekte basierend auf Level-2-Prädiktoren
  scaling_factor <- 3*(level_2_predictors[, 1])^2  # Beispiel: exponentielle Skalierung basierend auf Z1
  
  # Generiere zufällige Effekte mit Heteroskedastizität
  u <- matrix(0, nrow = Clstr, ncol = n_u)
  for (i in 1:Clstr) {
    scaled_sigma_u <- sigma_u_base * scaling_factor[i]
    u[i, ] <- mvrnorm(1, mu = rep(0, n_u), Sigma = scaled_sigma_u)
  }
  
  x <- cbind(1, predictors)
  
  # Cluster-ID hinzufügen und Level-2-Prädiktoren replizieren
  cluster_id <- rep(1:Clstr, each = N)
  level_2_data <- data.frame(cluster_id = 1:Clstr, level_2_predictors)
  level_2_expanded <- level_2_data[rep(1:Clstr, each = N), ]
  
  # Kombiniere Level-1- und Level-2-Daten
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  sim_data <- cbind(sim_data, level_2_expanded[, -1]) # Entferne die doppelte Spalte cluster_id
  
  # Formel inklusive Interaktionen
  m <- y ~ x1 + x2 + Z1 * Z2 + (1 + x1 + x2 | cluster_id)
  
  # Designmatrix Z als Blockmatrix erstellen
  z <- sparseMatrix(i = rep(1:(N * Clstr), each = ncol(x)),
                    j = rep(1:(Clstr * ncol(x)), times = N),
                    x = as.vector(t(x)),
                    dims = c(N * Clstr, Clstr * ncol(x)))
  
  y <- x %*% fixed + z %*% as.vector(t(u)) + e
  y <- as.vector(y)
  sim_data$y <- y
  
  # Fitte das Modell mit lmer
  fit <- lmer(m, data = sim_data)
  
  sum <- summary(fit)
  coef <- sum$coefficients
  varcorr <- VarCorr(fit)
  
  return(coef)
}

# Simulationen durchführen
results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# Funktion um Estimates zu mitteln
apply(results, c(1, 2), mean)

# Funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)






#### sollte klappen #### 

library(MASS)
library(lme4)
library(Matrix)
library(tmvtnorm)

# Parameter
N <- 90 # Anzahl der Beobachtungen pro Cluster
simulations <- 200 # Anzahl der Simulationen
lvl1_pred <- c(1, 1, 1) # level 1 Prädiktoren
n_pred <- length(lvl1_pred) - 1 # Anzahl der Prädiktoren auf Level-1
n_u <- length(lvl1_pred) # Anzahl der zufälligen Effekte
Clstr <- 120 # Anzahl der Cluster
lvl2_pred <- c(1,1,1) #level 2 Prädiktoren
N_lvl2 <- length(lvl2_pred) # Anzahl der Prädiktoren auf Level-2
all_pred <- c(1,1,1,1,1,1)

# Kovarianzmatrix für die Level-1 Prädiktoren
sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

# Kovarianzmatrix für die random effects
sigma_u_base <- matrix(c(2, 0, 0,
                         0, 2, 0,
                         0, 0, 2), nrow = 3, ncol = 3)

# Kovarianzmatrix für die Level-2 Prädiktoren
sigma_lvl2 <-   matrix(c(1.5, 0, 0,
                         0, 1.5, 0,
                         0, 0, 1.5), nrow = 3, ncol = 3)

# Funktion zur Simulation
function_sim <- function() {
  
  # Generiere Level-1 Prädiktoren
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  # Generiere Level-1 Fehler
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[,1])^2)
  
  # Generiere Level-2 Prädiktoren
  level_2_predictors <- matrix(rtmvnorm(Clstr, mean = rep(0, N_lvl2), sigma = sigma_lvl2), ncol = N_lvl2)
  colnames(level_2_predictors) <- c("X3", "X4", "X5")
  
  # Skalierungsfaktor für die Varianz der zufälligen Effekte basierend auf Level-2-Prädiktoren
  scaling_factor <- 3 * (level_2_predictors[, 1])^2  
  
  # Generiere zufällige Effekte mit Heteroskedastizität
  u <- matrix(0, nrow = Clstr, ncol = n_u)
  for (i in 1:Clstr) {
    scaled_sigma_u <- sigma_u_base * scaling_factor[i]
    u[i, ] <- mvrnorm(1, mu = rep(0, n_u), Sigma = scaled_sigma_u)
  }
  
  # Cluster-ID hinzufügen und Level-2-Prädiktoren replizieren
  cluster_id <- rep(1:Clstr, each = N)
  level_2_data <- data.frame(cluster_id = 1:Clstr, level_2_predictors)
  level_2_expanded <- level_2_data[rep(1:Clstr, each = N), ]
  
  level_2_predictors <- level_2_expanded
  
  x <- cbind(1, predictors,level_2_predictors)
  x <- subset(x, select =-cluster_id)
  colnames(x) <- c("Intercept", "x1", "x2", "lx3", "lx4","lx5")
  x<- as.matrix(x)
  # Kombiniere Level-1- und Level-2-Daten
  sim_data <- data.frame(y = rep(0, N * Clstr), x[,-1], cluster_id = cluster_id)
  
  # Formel inklusive Interaktionen
  m <- y ~ x1 + x2 + lx3 + lx4 * lx5 + (1 + x1 + x2 | cluster_id)
  
  # Designmatrix Z als Blockmatrix erstellen
  z <- sparseMatrix(i = rep(1:(N * Clstr), each = ncol(x[,1:3])),
                    j = rep(1:(Clstr * ncol(x[,1:3])), times = N),
                    x = as.vector(t(x[,1:3])),
                    dims = c(N * Clstr, Clstr * ncol(x[,1:3])))
  
  # Berechne die abhängige Variable y
  y <- x %*% all_pred + z %*% as.vector(t(u)) + e
  y <- as.vector(y)
  sim_data$y <- y
  
  # Fitte das Modell mit lmer
  fit <- lmer(m, data = sim_data)
  
  sum <- summary(fit)
  coef <- sum$coefficients
  varcorr <- VarCorr(fit)
  
  return(coef)
}

# Simulationen durchführen
results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# Funktion um Estimates zu mitteln
mean_estimates <- apply(results, c(1, 2), mean)
print(mean_estimates)






#### aufgeräumt ####

# Libraries
library(tmvtnorm)
library(lme4)
library(MASS)
library(Matrix)

# Parameters
N <- 120
simulations <- 100
lvl1_pred <- c(1, 1, 1)
n_pred <- length(lvl1_pred) - 1
n_u <- length(lvl1_pred)
Clstr <- 150
lvl2_pred <- c(1, 1, 1)
N_lvl2 <- length(lvl2_pred)
all_pred <- c(1, 1, 1, 1, 1, 1)

# Covariance matrices
sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

sigma_u_base <- matrix(c(2, 0, 0,
                         0, 2, 0,
                         0, 0, 2), nrow = 3, ncol = 3)

sigma_lvl2 <- matrix(c(1.5, 0, 0,
                       0, 1.5, 0, 
                       0, 0, 1.5), nrow = 3, ncol = 3)

# Simulation function
function_sim <- function() {
  # Level-1 predictors
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[, 1])^2)
  
  # Level-2 predictors
  level_2_predictors <- matrix(rtmvnorm(Clstr, mean = rep(0, N_lvl2), sigma = sigma_lvl2), ncol = N_lvl2)
  colnames(level_2_predictors) <- c("X3", "X4", "X5")
  
  # Random effects variance scaling
  scaling_factor <- 3 * (level_2_predictors[, 1])^2
  u <- matrix(0, nrow = Clstr, ncol = n_u)
  for (i in 1:Clstr) {
    scaled_sigma_u <- sigma_u_base * scaling_factor[i]
    u[i, ] <- mvrnorm(1, mu = rep(0, n_u), Sigma = scaled_sigma_u)
  }
  
  # Combine Level-1 and Level-2 data
  cluster_id <- rep(1:Clstr, each = N)
  level_2_expanded <- level_2_predictors[rep(1:Clstr, each = N), ]
  x <- cbind(1, predictors, level_2_expanded)
  colnames(x) <- c("Intercept", "x1", "x2", "lx3", "lx4", "lx5")
  
  # Create data frame
  sim_data <- data.frame(y = rep(0, N * Clstr),  x[,-1] , cluster_id = cluster_id)
  
  # Formula for model
  m <- y ~ x1 + x2 + lx3 + lx4 + lx5 + (1  | cluster_id)
  
  # Design matrix Z
  z <- sparseMatrix(i = rep(1:(N * Clstr), each = ncol(x[, 1:3])), 
                    j = rep(1:(Clstr * ncol(x[, 1:3])), times = N), 
                    x = as.vector(t(x[, 1:3])), dims = c(N * Clstr, Clstr * ncol(x[, 1:3])))
  
  # Generate response variable y
  y <- x %*% all_pred + z %*% as.vector(t(u)) + e
  sim_data$y <- as.vector(y)
  
  # Fit model using lmer
  fit <- lmer(m, data = sim_data)
  
  # Extract coefficients
  varcorr <- VarCorr(fit)
  coef <- summary(fit)$coefficients
  return(coef)
}

# Run simulations
results <- sapply(1:simulations, function(i) function_sim(), simplify = "array")

# Compute mean estimates
mean_estimates <- apply(results, c(1, 2), mean)
print(mean_estimates)



# Funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)





#### Metaskript ####




# Libraries
library(tmvtnorm)
library(lme4)
library(MASS)
library(Matrix)
library(lmerTest)
library(readxl)


setwd("/Users/moritz/Desktop/Documents/R/Daten")
meta <-read_xlsx("Metascript Multilevel HCCM.xlsx")


simulations <- 500


#### Simulation function ####
function_sim <- function(N, Clstr, lvl1_pred, sigma, sigma_u_base,
                         n_u, lvl2_pred,N_lvl2, all_pred, sigma_lvl2  ) {
  # Level-1 predictors
  predictors <-  rbinom(N * Clstr, 1, 0.5)  #rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  e <- rnorm(N * Clstr, mean = 0, sd = 1)    # für hetero (predictors[, 1])^2
  
  level_2_predictors <- rbinom( Clstr, 1, 0.5)
  
  # Combine Level-1 and Level-2 data
  cluster_id <- rep(1:Clstr, each = N)
  level_2_expanded <- rep(level_2_predictors, each = N)
  x <- cbind(1, predictors, level_2_expanded)
  colnames(x) <- c("Intercept", "x1", "x2")
  sim_data <- data.frame(y = rep(0, N * Clstr),  x[,-1] , cluster_id = cluster_id)
  
  # Design matrix Z

  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u_base)
  
  
  # Generate response variable y
  
  x <- cbind(1, predictors, level_2_expanded)
  colnames(x) <- c("Intercept", "x1", "x2")
  sim_data <- data.frame(y = rep(0, N * Clstr),  x[,-1] , cluster_id = cluster_id)
  
  y <- x %*% all_pred + z %*% as.vector(t(u)) + e
  sim_data$y <- as.vector(y)
  
  # Formula for model
  m <- y ~ x1 + x2  + (1 + x1  | cluster_id) # + (1| Word)
  
  # Fit model using lmer
  fit <- lmer(m, data = sim_data, REML = T)
  
  # Extract coefficients
  varcorr <- VarCorr(fit)
  coef_ <- summary(fit)$coefficients
  p_value <- coef_["x1","Pr(>|t|)"]
  
  return( p_value)
}
#### ####
results_list <- list()


for (i in 1:nrow(meta)) {
  N <- as.numeric(meta[i, "N"]) # Number of observations per cluster
  Clstr <- as.numeric(meta[i, "Clstr"])  # Number of clusters
  lvl1_pred <- as.numeric(meta[i, grep("b", names(meta))])  # Fixed effects coefficients
  sigma <- as.numeric(meta[i, "sigma"])  # Level-1 error variance
  sigma_u_base <- matrix(as.numeric(c(meta[i, "var_u1"], meta[i, "cov_u12"], meta[i, "cov_u12"], meta[i, "var_u2"])), nrow = 2, ncol = 2)  # Random effects covariance matrix
  n_u <- length(lvl1_pred)
  lvl2_pred <- as.numeric(meta[i, grep("p", names(meta))])
  N_lvl2 <- length(lvl2_pred)
  all_pred <- c(lvl1_pred,lvl2_pred)
  sigma_lvl2 <- as.numeric(meta[i, "sigma_lvl2"])
  
  # Run simulations
  p_values <- sapply(1:simulations, function(j) function_sim(N, Clstr, lvl1_pred, sigma, sigma_u_base,
                                                             n_u, lvl2_pred,N_lvl2, all_pred, sigma_lvl2), simplify = "array")
  
  # Compute power
  power <- mean(p_values < 0.01)
  
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


# expand.grid









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
simulations <- 3000

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
  power <- mean(p_values < 0.05) # ich glaube, die nehmen 0.05
  
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








# sanity check

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
simulations <- 1000

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
  
  # Extract covariance matrix of random effects
  cov_matrix_u <- as.matrix(VarCorr(fit)$cluster_id)
  
  return(list(coef = coef_, cov_u = cov_matrix_u))
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
  
  # Run simulations and store coefficients and covariance matrices
  sim_results <- replicate(simulations, function_sim(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2), simplify = FALSE)
  
  # Separate coefficients and covariance matrices
  coef_list <- lapply(sim_results, function(res) res$coef)
  cov_u_list <- lapply(sim_results, function(res) res$cov_u)
  
  # Convert lists to arrays for averaging
  coef_array <- simplify2array(coef_list)
  cov_u_array <- simplify2array(cov_u_list)
  
  # Compute the mean of coefficients and covariance matrix of random effects across all simulations
  mean_coef <- apply(coef_array, 1:2, mean)
  mean_cov_u <- apply(cov_u_array, 1:2, mean)
  
  # Store results
  results_list[[i]] <- list(
    condition = i,
    mean_coef = mean_coef,
    mean_cov_u = mean_cov_u
  )
}

# Summarize results: extract mean coefficients and covariance matrices for each condition
mean_coef_results <- lapply(results_list, function(res) res$mean_coef)
mean_cov_u_results <- lapply(results_list, function(res) res$mean_cov_u)

print(mean_coef_results)
print(mean_cov_u_results)






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
simulations <- 1000

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
simulations <- 1000

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
  profile <- confint(fit, parm = "x1", method = "profile",devmatchtol = 0.05)
  
  # Check if 0 is in the CI and assign P_value
  P_value <- ifelse(profile[1, 1] * profile[1, 2] < 0, 0, 1)
  
  return(P_value)
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



#### wald und profile ####

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
simulations <- 1000

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
  wald <- confint(fit, parm = "x1", method = "Wald", level= 0.95)
  profile <- confint(fit, parm = "x1", method = "profile", devmatchtol = 0.05 , level = 0.95)
  
  # Check if 0 is in the CI and assign P_value
  P_wald <- ifelse(wald[1, 1] * wald[1, 2] < 0, 0, 1)
  P_profile <- ifelse(profile[1, 1] * profile[1, 2] < 0, 0, 1)
  
  return(list(P_wald = P_wald, P_profile = P_profile))
}

# Initialize list to store results
results_list3 <- list()

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
  simulation_results <- lapply(1:simulations, function(j) {
    tryCatch({
      function_sim(N, Clstr, lvl1_pred, sigma, sigma_u_base, n_u, lvl2_pred, N_lvl2, all_pred, sigma_lvl2)
    }, error = function(e) {
      return(list(P_wald = NA, P_profile = NA))
    })
  })
  
  # Extract P_values, Wald CIs, and Profile CIs
  p_wald <- sapply(simulation_results, function(res) res$P_wald)
  p_profile <- sapply(simulation_results, function(res) res$P_profile)
  
  # Compute power
  power_wald <- ifelse(is.null(p_wald), NA, mean(p_wald == 1, na.rm = TRUE))
  power_profile <- ifelse(is.null(p_profile), NA, mean(p_profile == 1, na.rm = TRUE))
  
  # Store results
  results_list3[[i]] <- list(
    condition = i,
    power_wald = power_wald,
    power_profile = power_profile
  )
}


# Summarize results
power_wald_results <- sapply(results_list3, function(res) res$power_wald)
print(power_wald_results)

power_profile_results <- sapply(results_list3, function(res) res$power_profile)
print(power_profile_results)






