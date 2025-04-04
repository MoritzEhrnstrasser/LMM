
# ---- keep it maximal ----
plan(multisession, workers = 9)
setwd("~/Documents/Uni Kram/R")

library(lme4)
library(Matrix)
library(tmvtnorm)
library(lmerTest)
library(lmtest)
library(broom.mixed)
library(future)
library(future.apply)
library(CR2)
library(sloop)
library(MASS) 

# ---- AI ----
# Funktion zur Simulation eines Datensatzes
simulate_data <- function(n_subjects = 24, n_items = 12, b0 = 0, b1 = 0,
                          sigma_e = 1, 
                          mu_subject = c(0, 0), 
                          Sigma_subject = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                          sigma_item = 1) {
  # Erstelle Datenrahmen: jede Kombination von Subjekt und Item
  design <- expand.grid(Subject = factor(1:n_subjects),
                        Item = factor(1:n_items))
  
  # Simuliere zufällige Subjekt-Effekte (Intercept und Slope)
  subj_effects <- mvrnorm(n = n_subjects, mu = mu_subject, Sigma = Sigma_subject)
  colnames(subj_effects) <- c("S0", "S1")
  design$S0 <- subj_effects[as.numeric(design$Subject), "S0"]
  design$S1 <- subj_effects[as.numeric(design$Subject), "S1"]
  
  # Simuliere zufällige Item-Effekte (nur random intercepts)
  design$I0 <- rep(rnorm(n_items, mean = 0, sd = sigma_item), each = n_subjects)
  
  # Definiere den Treatment-Faktor X (innerhalb-Subjekte, z. B. deviation coding)
  design$X <- sample(c(-0.5, 0.5), nrow(design), replace = TRUE)
  
  # Simuliere den Fehlerterm
  design$epsilon <- rnorm(nrow(design), mean = 0, sd = sigma_e)
  
  # Generiere die abhängige Variable Y nach dem Modell:
  # Y = b0 + S0 + I0 + (b1 + S1) * X + epsilon
  design$Y <- b0 + design$S0 + design$I0 + (b1 + design$S1) * design$X + design$epsilon
  
  return(design)
}

# Funktion, die einen Simulationsdurchlauf durchführt, ein LMM fitet
# und wichtige Ergebnisse (z. B. Schätzer, t-Wert und p-Wert) zurückgibt
run_simulation <- function(n_subjects = 24, n_items = 12, 
                           b0 = 0, b1 = 0,
                           sigma_e = 1, 
                           mu_subject = c(0, 0), 
                           Sigma_subject = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                           sigma_item = 1) {
  
  # Simuliere Datensatz
  dat <- simulate_data(n_subjects, n_items, b0, b1, sigma_e, mu_subject, Sigma_subject, sigma_item)
  
  # Fit des LMM mit maximaler Random-Effects-Struktur:
  # Random intercept und random slope für Subject, random intercept für Item
  mod <- lmer(Y ~ X + (1 + X | Subject) + (1 | Item), data = dat, REML = FALSE)
  
  # Extrahiere Schätzer und t-Wert für den fixen Effekt X
  mod_summary <- summary(mod)
  est <- mod_summary$coefficients["X", "Estimate"]
  t_val <- mod_summary$coefficients["X", "t value"]
  # p-Wert
  p_val <- mod_summary$coefficients["X", "Pr(>|t|)"]
  
  return(c(Estimate = est, t_value = t_val, p_value = p_val))
}



# Simulation mehrmals mit lapply
n_sim <- 2000
sim_list <- future_lapply(1:n_sim, function(i) run_simulation())
# Ergebnisse in eine Matrix umwandeln
sim_results <- do.call(rbind, sim_list)

# Einen kleinen Ausschnitt der Ergebnisse anzeigen
head(sim_results)

alpha <- mean(sim_results[,"p_value"] < 0.05)

# Mittlere Werte für Estimate, t-Wert und p-Wert über alle Simulationen
sd_results <- apply(sim_results, 2, mean)
print(sd_results)






#---- eigener versuch ----

n_subjects = 24
n_items = 12
b0 = 0
b1 = 0
sigma_e = 1
mu_subject = c(0, 0)
mu_item = c(0, 0)
Sigma_subject = matrix(c(1, 0.5, 0.5, 1), nrow = 2)
Sigma_item = matrix(c(1, 0.5, 0.5, 1), nrow = 2)

simulate_data <- function(n_subjects, n_items , b0 , b1 ,
                          sigma_e, 
                          mu_subject, 
                          Sigma_subject, nrow),
                          sigma_item) {
                            
  # Erstelle Datenrahmen: jede Kombination von Subjekt und Item
  design <- expand.grid(Subject = factor(1:n_subjects),
                        Item = factor(1:n_items))
  
  # Simuliere zufällige Subjekt-Effekte (Intercept und Slope)
  subj_effects <- mvrnorm(n = n_subjects, mu = mu_subject, Sigma = Sigma_subject)
  colnames(subj_effects) <- c("S0", "S1")
  design$S0 <- subj_effects[as.numeric(design$Subject), "S0"]
  design$S1 <- subj_effects[as.numeric(design$Subject), "S1"]
  
  
  # Simuliere zufällige Itemeffekte (Intercept und Slope) aus MVN
  item_effects <- mvrnorm(n = n_items, mu = mu_item, Sigma = Sigma_item)
  colnames(item_effects) <- c("I0", "I1")
  # Ordne die Itemeffekte anhand der Item-ID zu
  design$I0 <- item_effects[as.numeric(design$Item), "I0"]
  design$I1 <- item_effects[as.numeric(design$Item), "I1"]
  
  
  
  
  
  
  # Definiere den Treatment-Faktor X (innerhalb-Subjekte, z. B. deviation coding)
  design$X <- sample(c(-0.5, 0.5), nrow(design), replace = TRUE)
  
  # Simuliere den Fehlerterm
  design$epsilon <- rnorm(nrow(design), mean = 0, sd = sigma_e)
  
  # Generiere die abhängige Variable Y nach dem Modell:
  # Y = b0 + S0 + I0 + (b1 + S1) * X + epsilon
  design$Y <- b0 + design$S0 + design$I0 + (b1 + design$S1 + design$I1) * design$X + design$epsilon
  
  mod <- lmer(Y ~ X + (1 + X | Subject) + (1 | Item), data = design, REML = T)
  summ <- summary(mod)
  coef <- summ$
  return()
}