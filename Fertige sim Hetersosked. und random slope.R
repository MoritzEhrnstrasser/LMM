setwd("~/Documents/Uni Kram/R")
plan(multicore, workers = 9)
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


 #### mit Parallelisierung ####
 
 plan(multisession, workers = 9)
 
 
 simulate_all_models <- function(isample, N, Clstr, all_pred) {
   
   # --- Simulation ohne zufällige Steigung (heteroskedastische Fehler, nur zufälliger Intercept) ---
   simulate_noSlope <- function(N, Clstr, all_pred) {
     lvl1_pred <- 1
     sigma <- diag(lvl1_pred)
     # Level-1 Prädiktor
     predictors <- rtmvnorm(N * Clstr, mean = rep(0, lvl1_pred), sigma = sigma)
     # Heteroskedastische Residuen
     e <- rnorm(N * Clstr, mean = 0, sd = 3 * (predictors[, 1])^4)
     e <- e * (sqrt(1) / sd(e))
     # Cluster-IDs und Designmatrix X
     cluster_id <- rep(1:Clstr, each = N)
     x <- cbind(1, predictors)
     colnames(x) <- c("Intercept", "x1")
     sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, "x1"], cluster_id = cluster_id)
     # Designmatrix Z nur für zufälligen Intercept
     z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), "Intercept", drop = FALSE])
     z <- do.call(bdiag, z_blocks)
     # Zufallseffekte (nur Intercept)
     sigma_u_base <- diag(1)
     u <- rtmvnorm(Clstr, mean = 0, sigma = sigma_u_base)
     u_vector <- as.vector(u)
     
     # Antwortvariable y generieren
     y <- x %*% all_pred + z %*% u_vector + e
     sim_data$y <- as.vector(y)
     
     # Modelle formulieren:
     m  <- y ~ x1 + (1 + x1 | cluster_id)  # ws-Modell
     m2 <- y ~ x1 + (1 | cluster_id)        # ns-Modell
     
     # Modelle schätzen
     fit_ws    <- lmer(m, data = sim_data, REML = TRUE)
     fit_ns    <- lmer(m2, data = sim_data, REML = TRUE)
     fit_ws_cr <- robust_mixed(lmer(m, data = sim_data, REML = TRUE), type = "CR2")
     fit_ns_cr <- robust_mixed(lmer(m2, data = sim_data, REML = TRUE), type = "CR2")
     
     # Ergebnisse extrahieren
     result_ws <- list(
       coefficients = summary(fit_ws)$coefficients,
       p_value      = summary(fit_ws)$coefficients["x1", "Pr(>|t|)"],
       varcorr      = VarCorr(fit_ws)
     )
     result_ns <- list(
       coefficients = summary(fit_ns)$coefficients,
       p_value      = summary(fit_ns)$coefficients["x1", "Pr(>|t|)"],
       varcorr      = VarCorr(fit_ns)
     )
     # Bei den robusten Fits wird das varcorr-Objekt nicht zurückgegeben:
     result_ws_cr <- list(
       coefficients = fit_ws_cr$results[, "Estimate"],
       p_value      = fit_ws_cr$results["x1", 6]
     )
     result_ns_cr <- list(
       coefficients = fit_ns_cr$results[, "Estimate"],
       p_value      = fit_ns_cr$results["x1", 6]
     )
     
     return(list(ws = result_ws, ns = result_ns, ws_cr = result_ws_cr, ns_cr = result_ns_cr))
   }
   
   # --- Simulation mit zufälliger Steigung (homoskedastische Fehler, zufälliger Intercept und zufälliger Slope) ---
   simulate_withSlope <- function(N, Clstr, all_pred) {
     lvl1_pred <- 1
     sigma <- diag(lvl1_pred)
     # Level-1 Prädiktor
     predictors <- rtmvnorm(N * Clstr, mean = rep(0, lvl1_pred), sigma = sigma)
     # Homoskedastische Residuen
     e <- rnorm(N * Clstr, mean = 0, sd = 1)
     e <- e * (sqrt(1) / sd(e))
     # Cluster-IDs und Designmatrix X
     cluster_id <- rep(1:Clstr, each = N)
     x <- cbind(1, predictors)
     colnames(x) <- c("Intercept", "x1")
     sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, "x1"], cluster_id = cluster_id)
     # Designmatrix Z für zufälligen Intercept und zufälligen Slope
     z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
     z <- do.call(bdiag, z_blocks)
     # Zufallseffekte (Intercept und Slope)
     sigma_u_base <- diag(c(1, 0.1))
     u <- rtmvnorm(Clstr, mean = c(0, 0), sigma = sigma_u_base)
     u_vector <- as.vector(t(u))
     
     # Antwortvariable y generieren
     y <- x %*% all_pred + z %*% u_vector + e
     sim_data$y <- as.vector(y)
     
     # Modelle formulieren:
     m  <- y ~ x1 + (1 + x1 | cluster_id)  # ws-Modell
     m2 <- y ~ x1 + (1 | cluster_id)        # ns-Modell
     
     # Modelle schätzen
     fit_ws    <- lmer(m, data = sim_data, REML = TRUE)
     fit_ns    <- lmer(m2, data = sim_data, REML = TRUE)
     fit_ws_cr <- robust_mixed(lmer(m, data = sim_data, REML = TRUE), type = "CR2")
     fit_ns_cr <- robust_mixed(lmer(m2, data = sim_data, REML = TRUE), type = "CR2")
     
     # Ergebnisse extrahieren
     result_ws <- list(
       coefficients = summary(fit_ws)$coefficients,
       p_value      = summary(fit_ws)$coefficients["x1", "Pr(>|t|)"],
       varcorr      = VarCorr(fit_ws)
     )
     result_ns <- list(
       coefficients = summary(fit_ns)$coefficients,
       p_value      = summary(fit_ns)$coefficients["x1", "Pr(>|t|)"],
       varcorr      = VarCorr(fit_ns)
     )
     result_ws_cr <- list(
       coefficients = fit_ws_cr$results[, "Estimate"],
       p_value      = fit_ws_cr$results["x1", 6]
     )
     result_ns_cr <- list(
       coefficients = fit_ns_cr$results[, "Estimate"],
       p_value      = fit_ns_cr$results["x1", 6]
     )
     
     return(list(ws = result_ws, ns = result_ns, ws_cr = result_ws_cr, ns_cr = result_ns_cr))
   }
   
   # --- Simulationen parallel durchführen (z. B. mittels future_lapply) ---
   results_noSlope  <- future_lapply(1:isample, function(i) simulate_noSlope(N, Clstr, all_pred))
   results_withSlope <- future_lapply(1:isample, function(i) simulate_withSlope(N, Clstr, all_pred))
   
   # --- Funktion zur Aggregation der Ergebnisse ---
   aggregate_results <- function(results, model_types = c("ws", "ns", "ws_cr", "ns_cr")) {
     agg_list <- list()
     for(mt in model_types) {
       # Berechnung der Power (Anteil signifikanter p-Werte)
       p_values <- sapply(results, function(res) res[[mt]]$p_value)
       power    <- mean(p_values < 0.05)
       # Aggregation der geschätzten Koeffizienten
       if (mt %in% c("ws", "ns")) {
         coef_matrix <- sapply(results, function(res) res[[mt]]$coefficients[, "Estimate"])
       } else {
         coef_matrix <- sapply(results, function(res) res[[mt]]$coefficients)
       }
       mean_estimates <- rowMeans(coef_matrix)
       
       # Nur für die lmer-basierten Modelle (ws, ns) wird varcorr aggregiert
       if (mt %in% c("ws", "ns")) {
         allVarcorr <- t(sapply(results, function(res) {
           if (!is.null(res[[mt]]$varcorr)) {
             var_df <- as.data.frame(res[[mt]]$varcorr)
             return(var_df$vcov)
           } else {
             return(rep(NA, length(as.data.frame(results[[1]][[mt]]$varcorr)$vcov)))
           }
         }, simplify = TRUE))
         mean_varcorr <- apply(allVarcorr, 2, mean, na.rm = TRUE)
         agg_list[[mt]] <- list(power = power, mean_estimates = mean_estimates, mean_varcorr = mean_varcorr)
       } else {
         agg_list[[mt]] <- list(power = power, mean_estimates = mean_estimates)
       }
     }
     return(agg_list)
   }
   
   # --- Ergebnisse für beide Simulationsdesigns aggregieren ---
   aggregated_noSlope   <- aggregate_results(results_noSlope)
   aggregated_withSlope <- aggregate_results(results_withSlope)
   
   return(list(noSlope = aggregated_noSlope, withSlope = aggregated_withSlope))
 }
 
 # Beispielaufruf (parallel vorausgesetzt, z. B. mit plan(multisession)):
 result_all <- simulate_all_models(isample = 2500, N = 50, Clstr = 50, all_pred = c(0, 0))
 
 save(result_all, file = "result_LMM_heterosked_50.RData")
 # Zum Laden später:
 load("result_all.RData")
 
 # für Git
 
 