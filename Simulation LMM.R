#### Simulation LMM ####

library(sandwich)
library(lmtest)
library(tmvtnorm)
library(truncnorm)
library(ggplot2)
library(abind)
library(tidyverse)
library(metRology)
library(lme4)
library(lmerTest)

rm(list = ls())

simulations <- 100

# Setzen der Parameter für die Simulation
n <- 100 # Anzahl der Beobachtungen pro Gruppe
m <- 50   # Anzahl der Gruppen
beta <- c(2,3) # Fixe Effekte
sigma_b <- 1    # Standardabweichung der Gruppeneffekte
sigma_e <- 2    # Standardabweichung der Residuen

function_sim <- function(data) {
# Erstellen der Gruppenindikatoren
group <- rep(1:m, each = n)

# Simulieren der Gruppeneffekte
b <- rnorm(m, 0, sigma_b)

# Simulieren der abhängigen Variablen
x <- rnorm(n * m)
y <- beta[1] + beta[2]*x + b[group] + rnorm(n * m, 0, sigma_e)

# Erstellen des Datenrahmens
data <- data.frame(y = y, x = x, group = factor(group))

# Anpassen des linearen gemischten Modells
fit <- lmer(y ~ x + (1 | group), data = data)

# Zusammenfassung des Modells
sum_fit <- summary(fit)

  
return(sum_fit$coefficients) 
}


results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

results_mean_df<-as.data.frame(apply(results, MARGIN = 1:2, mean))

