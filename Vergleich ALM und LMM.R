#### ALM ####

isample = 10000
N = 2500
B0 = 0
B1 = 0
sd = 1

# Generate data
 generate_data<-function(N, B0, B1, sd) {
 x <- rnorm(N, mean=0, sd = 1)
 e <- rnorm(N , mean = 0, sd = 3*x^2)
 e <- e*(sqrt(1)/sd(e))
 y <- B0 + B1 * x + e

 fit <- lm(y ~ x)
sum_fit <-  summary(fit)$coefficients
 
 p_value <- sum_fit["x","Pr(>|t|)"]
    return(p_value)
}

 p_values <- sapply(1:isample, function(i) generate_data(N, B0, B1, sd))

 mean(p_values < 0.05) # 0.05 is the significance level

 # 0.3878 bei ^2
 
 
 #### LMM random intercept only daten generierend #### 
 
 library(lme4)
 library(Matrix)
 library(tmvtnorm)
 library(lmerTest)
 library(lmtest)
 library(CR2)
 
 # Anzahl der Simulationen
 isample <- 10000
 N <- 100       # Beobachtungen pro Cluster
 Clstr <- 50   # Anzahl der Cluster
 
 # Funktion zur Simulation
 simulate_random_intercept <- function(N, Clstr) {
   # Anzahl der Cluster und Level-1-Prädiktor
   lvl1_pred <- 1  # Ein Level-1-Prädiktor
   sigma <- diag(lvl1_pred)  # Kovarianzmatrix für den Prädiktor (Identitätsmatrix)
   
   # Level-1 Prädiktor
   predictors <- rtmvnorm(N * Clstr, mean = rep(0, lvl1_pred), sigma = sigma)
   
   # Heteroskedastische Residuen
   e <- rnorm(N * Clstr, mean = 0, sd = 3 * (predictors[, 1])^2) 
   e <- e * (sqrt(1) / sd(e))  # Standardisierung des Fehlers
   
   # Cluster-IDs
   cluster_id <- rep(1:Clstr, each = N)
   
   # Designmatrix X
   x <- cbind(1, predictors)
   colnames(x) <- c("Intercept", "x1")
   sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, "x1"], cluster_id = cluster_id)
   
   # Designmatrix Z für Zufallseffekte (Intercept)
   z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), "Intercept", drop = FALSE])
   z <- do.call(bdiag, z_blocks)
   
   # Zufallseffekte u
   sigma_u_base <- diag(1)  # Varianzkomponente für den zufälligen Intercept
   u <- rtmvnorm(Clstr, mean = 0, sigma = sigma_u_base)
   u_vector <- as.vector(u)
   
   # Antwortvariable y
   all_pred <- c(0, 0)  # Koeffizienten für Intercept und Prädiktor (x1)
   y <- x %*% all_pred + z %*% u_vector + e
   sim_data$y <- as.vector(y)
   
   # Modellformel (random intercept only)
   m <- y ~ x1 + (1 + x1 | cluster_id)
   
   # Modell anpassen
   fit <- lmer(m, data = sim_data, REML = TRUE)
  
   
   # Koeffizienten und p-Wert extrahieren
   coef_ <- summary(fit)$coefficients
   p_value <- coef_["x1", "Pr(>|t|)"]
   
   return(p_value)
 }
 
 # Simulation
 p_values <- sapply(1:isample, function(i) simulate_random_intercept(N, Clstr))
 
 # Ergebnis: Anteil der signifikanten p-Werte
 mean(p_values < 0.05)  # 0.05 ist das Signifikanzniveau
 
 
 
 # 0.3736 bei ^2
 
 
 #### LMM random Intercept random slope ####
 
 library(lme4)
 library(Matrix)
 library(tmvtnorm)
 library(lmerTest)
 library(lmtest)
 library(CR2)
 
 # Anzahl der Simulationen
 isample <- 3000
 N <- 100      # Beobachtungen pro Cluster
 Clstr <- 50   # Anzahl der Cluster
 
 # Funktion zur Simulation mit random slope
 simulate_random_slope <- function(N, Clstr) {
   # Anzahl der Cluster und Level-1-Prädiktor
   lvl1_pred <- 1  # Ein Level-1-Prädiktor
   sigma <- diag(lvl1_pred)  # Kovarianzmatrix für den Prädiktor (Identitätsmatrix)
   
   # Level-1 Prädiktor
   predictors <- rtmvnorm(N * Clstr, mean = rep(0, lvl1_pred), sigma = sigma)
   
   # Heteroskedastische Residuen
   e <- rnorm(N * Clstr, mean = 0, sd = 1)#3 * (predictors[, 1])^2)
   e <- e * (sqrt(1) / sd(e))  # Standardisierung des Fehlers
   
   # Cluster-IDs
   cluster_id <- rep(1:Clstr, each = N)
   
   # Designmatrix X
   x <- cbind(1, predictors)
   colnames(x) <- c("Intercept", "x1")
   sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, "x1"], cluster_id = cluster_id)
   
   # Designmatrix Z für Zufallseffekte (Intercept und Slope)
   z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
   z <- do.call(bdiag, z_blocks)
   
   # Zufallseffekte u (random intercept und random slope)
   sigma_u_base <- diag(c(1, 0.1))  # Varianzkomponenten für Intercept und Slope
   u <- rtmvnorm(Clstr, mean = c(0, 0), sigma = sigma_u_base)
   u_vector <- as.vector(t(u))
   
   # Antwortvariable y
   all_pred <- c(0, 0)  # Koeffizienten für Intercept und Prädiktor (x1)
   y <- x %*% all_pred + z %*% u_vector + e
   sim_data$y <- as.vector(y)
   
   # Modellformel (random intercept und slope)
   m <- y ~ x1 + (1  | cluster_id)
   
   # Modell anpassen
   
   
   fit <- robust_mixed(lmer(m, data = sim_data, REML = TRUE),type="CR2")
   p_value <- fit$results["x1","p.val"]
   

   
   # Koeffizienten und p-Wert extrahieren
   #coef_ <- summary(fit)$coefficients
   #p_value <- coef_["x1", "Pr(>|t|)"]
   
   return(p_value)
 }
 
 # Simulation
 p_values <- sapply(1:isample, function(i) simulate_random_slope(N, Clstr))
 
 # Ergebnis: Anteil der signifikanten p-Werte
 mean(p_values < 0.05)  # 0.05 ist das Signifikanzniveau
 
 
 # 0.0497 bei ^2 und N = 50 
 # 0.0525 bei N = 1000
 # 0.049 bei N = 2000
 
### ALM und random Intercep hetero selbes Problem, sobald random slope hinzukommt, 
 #ist das Problem gelöst. Warum?
 
 # wenn Random slope in daten generierendem aber nicht geschätzt und keine heterosked. gibt es 
 # maßgebliche Abweichungen vom p-Wert bei slope var von 0.1 = p-Wert 0.42 bei N=50 
 # problem wird größer mit höherer ANzahl an Level 1 Einheiten

 
 
 
 ####  fertige Simulation ####
 
 library(lme4)
 library(Matrix)
 library(tmvtnorm)
 library(lmerTest)
 library(lmtest)
 library(broom.mixed)
 
 # Anzahl der Simulationen
 isample <- 1000
 N <- 50       # Beobachtungen pro Cluster
 Clstr <- 50   # Anzahl der Cluster
 all_pred <- c(0, 0)  # Koeffizienten für Intercept und Prädiktor (x1)
 
 # Funktion zur Simulation
 simulate_random_intercept <- function(N, Clstr, all_pred) {
   # Anzahl der Cluster und Level-1-Prädiktor
   lvl1_pred <- 1  # Ein Level-1-Prädiktor
   sigma <- diag(lvl1_pred)  # Kovarianzmatrix für den Prädiktor (Identitätsmatrix)
   
   # Level-1 Prädiktor
   predictors <- rtmvnorm(N * Clstr, mean = rep(0, lvl1_pred), sigma = sigma)
   
   # Heteroskedastische Residuen
   
   e <- rnorm(N * Clstr, mean = 0, sd = 3 * (predictors[, 1])^4) 
   #e <- rnorm(N * Clstr, mean = 0, sd = 0.5 * exp(-0.5 * (predictors[, 1] - mean(predictors[,1]))^2))
   
   e <- e * (sqrt(1) / sd(e))  # Standardisierung des Fehlers
   
   # Cluster-IDs
   cluster_id <- rep(1:Clstr, each = N)
   
   # Designmatrix X
   x <- cbind(1, predictors)
   colnames(x) <- c("Intercept", "x1")
   sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, "x1"], cluster_id = cluster_id)
   
   # Designmatrix Z für Zufallseffekte (Intercept)
   z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), "Intercept", drop = FALSE])
   z <- do.call(bdiag, z_blocks)
   
   # Zufallseffekte u
   sigma_u_base <- diag(1)  # Varianzkomponente für den zufälligen Intercept
   u <- rtmvnorm(Clstr, mean = 0, sigma = sigma_u_base)
   u_vector <- as.vector(u)
   
   
   # Antwortvariable y
   
   y <- x %*% all_pred + z %*% u_vector + e
   sim_data$y <- as.vector(y)
   
   # Modellformel (random intercept only)
   m <- y ~ x1 + (1  + x1 | cluster_id) #ws
   m2 <- y ~ x1 + (1 | cluster_id) #ns
   
   # Modell anpassen
   fit_ws <- lmer(m, data = sim_data, REML = TRUE)
   
   fit_ns  <- lmer(m2, data = sim_data, REML = TRUE)
   
   fit_ws_cr <- robust_mixed(fit_ws,type="CR2")
   
   fit_ns_cr <- robust_mixed(fit_ns,type="CR2")
   
   
   
   
   # Koeffizienten und p-Wert extrahieren
   result_ws <- list(
     coefficients = summary(fit_ws)$coefficients,
     p_value = summary(fit_ws)$coefficients["x1", "Pr(>|t|)"],
     varcorr = VarCorr(fit_ws)
   )
   
   result_ns <- list(
     coefficients2 = summary(fit_ns)$coefficients,
     p_value2 = summary(fit_ns)$coefficients["x1", "Pr(>|t|)"],
     varcorr2 = VarCorr(fit_ns)
   )
  
   result_ws_cr <- list(
     coefficients = fit_ws_cr$results[,"Estimate"],
     p_value = fit_ws_cr$results["x1", 6],
     varcorr = fit_ws_cr[["vcov"]]
   )
   
   result_wn_cr <- list(
     coefficients = fit_ns_cr$results[,"Estimate"],
     p_value = fit_ns_cr$results["x1", 6],
     varcorr = fit_ns_cr[["vcov"]]
   )

   result_combined <- list(
     ws = result_ws,
     ns = result_ns,
     ws_cr = result_ws_cr,
     ns_cr = result_wn_cr)
     
   
   return(result_combined)
 }
 
 results <- lapply(1:isample, function(i) 
   simulate_random_intercept(N, Clstr, all_pred)
 )
 
 ### Für das "ws"-Modell ###
 
 # 1. p-Wert mitteln:
 ws_p_values <- sapply(results, function(res) res$ws$p_value)
 mean_ws_p <- mean(0.05 > ws_p_values)
 
 # 2. Koeffizienten mitteln (hier werden die "Estimate"-Werte pro Parameter gemittelt):
 ws_coef_matrix <- sapply(results, function(res) res$ws$coefficients[ , "Estimate"])
 # Das liefert eine Matrix, in der jede Zeile einen Parameter (z.B. "Intercept" und "x1") repräsentiert.
 mean_ws_coef <- rowMeans(ws_coef_matrix)
 
 allCov_ws <- t(sapply(results, FUN = function(iSample){
   if (!is.null(iSample$ws$varcorr)) {
     # Konvertiere VarCorr-Objekt in Dataframe und extrahiere Varianzen
     var_df <- as.data.frame(iSample$ws$varcorr)
     return(var_df$vcov)  # Nimmt nur die Varianz-Kovarianz-Werte
   } else {
     return(rep(NA, length(results[[1]]$ws$varcorr)))  # Falls NULL, ersetze mit NAs
   }
 }, simplify = TRUE))
 
 
 mean_varcov_ws <- apply(allCov_ws, MARGIN = 2, mean, na.rm = TRUE)
 
 ### Für das "ns"-Modell ###
 
 # 1. p-Wert mitteln:
 ns_p_values <- sapply(results, function(res) res$ns$p_value)
 mean_ns_p <- mean(0.05>ns_p_values)
 
 # 2. Koeffizienten mitteln:
 ns_coef_matrix <- sapply(results, function(res) res$ns$coefficients[ , "Estimate"])
 mean_ns_coef <- rowMeans(ns_coef_matrix)
 
 allCov_ns <- t(sapply(results, FUN = function(iSample){
   if (!is.null(iSample$ns$varcorr)) {
     # Konvertiere VarCorr-Objekt in Dataframe und extrahiere Varianzen
     var_df <- as.data.frame(iSample$ns$varcorr)
     return(var_df$vcov)  # Nimmt nur die Varianz-Kovarianz-Werte
   } else {
     return(rep(NA, length(results[[1]]$ns$varcorr)))  # Falls NULL, ersetze mit NAs
   }
 }, simplify = TRUE))
 
 
 mean_varcov_ns <- apply(allCov_ns, MARGIN = 2, mean, na.rm = TRUE)
 
 
 ### Berechnungen für das "ws_cr"-Modell (cluster robust) ###
 
 # 1. Anteil signifikanter p-Werte
 ws_cr_p_values <- sapply(results, function(res) res$ws_cr$p_value)
 mean_ws_cr_p <- mean(0.05 > ws_cr_p_values)
 
 # 2. Koeffizienten mitteln
 ws_cr_coef_matrix <- sapply(results, function(res) res$ws_cr$coefficients)
 mean_ws_cr_coef <- rowMeans(ws_cr_coef_matrix)
 
 # 3. Varianz-Kovarianzen mitteln
 allCov_ws_cr <- t(sapply(results, FUN = function(iSample) {
   if (!is.null(iSample$ws_cr$varcorr)) {
     return(as.vector(iSample$ws_cr$varcorr))
   } else {
     return(rep(NA, length(as.vector(results[[1]]$ws_cr$varcorr))))
   }
 }, simplify = TRUE))
 mean_varcov_ws_cr <- apply(allCov_ws_cr, MARGIN = 2, mean, na.rm = TRUE)
 
 
 ### Berechnungen für das "ns_cr"-Modell (cluster robust) ###
 
 # 1. Anteil signifikanter p-Werte
 ns_cr_p_values <- sapply(results, function(res) res$ns_cr$p_value)
 mean_ns_cr_p <- mean(0.05 > ns_cr_p_values)
 
 # 2. Koeffizienten mitteln
 ns_cr_coef_matrix <- sapply(results, function(res) res$ns_cr$coefficients)
 mean_ns_cr_coef <- rowMeans(ns_cr_coef_matrix)
 
 # 3. Varianz-Kovarianzen mitteln
 allCov_ns_cr <- t(sapply(results, FUN = function(iSample) {
   if (!is.null(iSample$ns_cr$varcorr)) {
     return(as.vector(iSample$ns_cr$varcorr))
   } else {
     return(rep(NA, length(as.vector(results[[1]]$ns_cr$varcorr))))
   }
 }, simplify = TRUE))
 mean_varcov_ns_cr <- apply(allCov_ns_cr, MARGIN = 2, mean, na.rm = TRUE)
 
 ### Zusammenfassen der gemittelten Ergebnisse ###
 result_mean <- list(
   ws = list(
     mean_p_value     = mean_ws_p,
     mean_coefficients = mean_ws_coef,
     mean_varcorr     = mean_varcov_ws
   ),
   ns = list(
     mean_p_value     = mean_ns_p,
     mean_coefficients = mean_ns_coef,
     mean_varcorr     = mean_varcov_ns
   ),
   ws_cr = list(
     mean_p_value     = mean_ws_cr_p,
     mean_coefficients = mean_ws_cr_coef,
     mean_varcorr     = mean_varcov_ws_cr
   ),
   ns_cr = list(
     mean_p_value     = mean_ns_cr_p,
     mean_coefficients = mean_ns_cr_coef,
     mean_varcorr     = mean_varcov_ns_cr
   )
 )
 
 # Ausgabe der gemittelten Ergebnisse:
 print(result_mean) 
 
 
 
 
 
 


#### mit slope in daten ####
 
 library(lme4)
 library(Matrix)
 library(tmvtnorm)
 library(lmerTest)
 library(lmtest)
 library(broom.mixed)
 
 # Anzahl der Simulationen
 isample <- 1000
 N <- 50       # Beobachtungen pro Cluster
 Clstr <- 50   # Anzahl der Cluster
 all_pred <- c(0, 0)  # Koeffizienten für Intercept und Prädiktor (x1)
 
 # Funktion zur Simulation
 simulate_random_intercept <- function(N, Clstr, all_pred) {
   # Anzahl der Cluster und Level-1-Prädiktor
   lvl1_pred <- 1  # Ein Level-1-Prädiktor
   sigma <- diag(lvl1_pred)  # Kovarianzmatrix für den Prädiktor (Identitätsmatrix)
   
   # Level-1 Prädiktor
   predictors <- rtmvnorm(N * Clstr, mean = rep(0, lvl1_pred), sigma = sigma)
   
   # Heteroskedastische Residuen
   
   e <- rnorm(N * Clstr, mean = 0, sd = 1) # 3 * (predictors[, 1])^4) 
   #e <- rnorm(N * Clstr, mean = 0, sd = 0.5 * exp(-0.5 * (predictors[, 1] - mean(predictors[,1]))^2))
   
   e <- e * (sqrt(1) / sd(e))  # Standardisierung des Fehlers
   
   # Cluster-IDs
   cluster_id <- rep(1:Clstr, each = N)
   
   # Designmatrix X
   x <- cbind(1, predictors)
   colnames(x) <- c("Intercept", "x1")
   sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, "x1"], cluster_id = cluster_id)
   
   # Designmatrix Z für Zufallseffekte (Intercept und Slope)
   z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), c("Intercept", "x1")])
   z <- do.call(bdiag, z_blocks)
   
   # Zufallseffekte u (random intercept und random slope)
   sigma_u_base <- diag(c(1, 0.1))  # Varianzkomponenten für Intercept und Slope
   u <- rtmvnorm(Clstr, mean = c(0, 0), sigma = sigma_u_base)
   u_vector <- as.vector(t(u))
   
   
   # Antwortvariable y
   
   y <- x %*% all_pred + z %*% u_vector + e
   sim_data$y <- as.vector(y)
   
   # Modellformel (random intercept only)
   m <- y ~ x1 + (1  + x1 | cluster_id) #ws
   m2 <- y ~ x1 + (1 | cluster_id) #ns
   
   # Modell anpassen
   fit_ws <- lmer(m, data = sim_data, REML = TRUE)
   
   fit_ns  <- lmer(m2, data = sim_data, REML = TRUE)
   
   fit_ws_cr <- robust_mixed(lmer(m, data = sim_data, REML = TRUE),type="CR2")
   
   fit_ns_cr <- robust_mixed(lmer(m2, data = sim_data, REML = TRUE),type="CR2")
   
   
   # Koeffizienten und p-Wert extrahieren
   result_ws <- list(
     coefficients = summary(fit_ws)$coefficients,
     p_value = summary(fit_ws)$coefficients["x1", "Pr(>|t|)"],
     varcorr = VarCorr(fit_ws)
   )
   
   result_ns <- list(
     coefficients2 = summary(fit_ns)$coefficients,
     p_value2 = summary(fit_ns)$coefficients["x1", "Pr(>|t|)"],
     varcorr2 = VarCorr(fit_ns)
   )
   
   result_ws_cr <- list(
     coefficients = fit_ws_cr$results[,"Estimate"],
     p_value = fit_ws_cr$results["x1", 6],
     varcorr = fit_ws_cr[["vcov"]]
   )
   
   result_wn_cr <- list(
     coefficients = fit_ns_cr$results[,"Estimate"],
     p_value = fit_ns_cr$results["x1", 6],
     varcorr = fit_ns_cr[["vcov"]]
   )
   
   result_combined <- list(
     ws = result_ws,
     ns = result_ns,
     ws_cr = result_ws_cr,
     ns_cr = result_wn_cr)
   
   
   return(result_combined)
 }
 
 results2 <- lapply(1:isample, function(i) 
   simulate_random_intercept(N, Clstr, all_pred)
 )
 
 ### Für das "ws"-Modell ###
 
 # 1. p-Wert mitteln:
 ws_p_values <- sapply(results2, function(res) res$ws$p_value)
 mean_ws_p <- mean(0.05 > ws_p_values)
 
 # 2. Koeffizienten mitteln (hier werden die "Estimate"-Werte pro Parameter gemittelt):
 ws_coef_matrix <- sapply(results2, function(res) res$ws$coefficients[ , "Estimate"])
 # Das liefert eine Matrix, in der jede Zeile einen Parameter (z.B. "Intercept" und "x1") repräsentiert.
 mean_ws_coef <- rowMeans(ws_coef_matrix)
 
 allCov_ws <- t(sapply(results2, FUN = function(iSample){
   if (!is.null(iSample$ws$varcorr)) {
     # Konvertiere VarCorr-Objekt in Dataframe und extrahiere Varianzen
     var_df <- as.data.frame(iSample$ws$varcorr)
     return(var_df$vcov)  # Nimmt nur die Varianz-Kovarianz-Werte
   } else {
     return(rep(NA, length(results2[[1]]$ws$varcorr)))  # Falls NULL, ersetze mit NAs
   }
 }, simplify = TRUE))
 
 
 mean_varcov_ws <- apply(allCov_ws, MARGIN = 2, mean, na.rm = TRUE)
 
 ### Für das "ns"-Modell ###
 
 # 1. p-Wert mitteln:
 ns_p_values <- sapply(results2, function(res) res$ns$p_value)
 mean_ns_p <- mean(0.05>ns_p_values)
 
 # 2. Koeffizienten mitteln:
 ns_coef_matrix <- sapply(results2, function(res) res$ns$coefficients[ , "Estimate"])
 mean_ns_coef <- rowMeans(ns_coef_matrix)
 
 allCov_ns <- t(sapply(results2, FUN = function(iSample){
   if (!is.null(iSample$ns$varcorr)) {
     # Konvertiere VarCorr-Objekt in Dataframe und extrahiere Varianzen
     var_df <- as.data.frame(iSample$ns$varcorr)
     return(var_df$vcov)  # Nimmt nur die Varianz-Kovarianz-Werte
   } else {
     return(rep(NA, length(results2[[1]]$ns$varcorr)))  # Falls NULL, ersetze mit NAs
   }
 }, simplify = TRUE))
 
 
 mean_varcov_ns <- apply(allCov_ns, MARGIN = 2, mean, na.rm = TRUE)
 
 
 ### Berechnungen für das "ws_cr"-Modell (cluster robust) ###
 
 # 1. Anteil signifikanter p-Werte
 ws_cr_p_values <- sapply(results2, function(res) res$ws_cr$p_value)
 mean_ws_cr_p <- mean(0.05 > ws_cr_p_values)
 
 # 2. Koeffizienten mitteln
 ws_cr_coef_matrix <- sapply(results2, function(res) res$ws_cr$coefficients)
 mean_ws_cr_coef <- rowMeans(ws_cr_coef_matrix)
 
 # 3. Varianz-Kovarianzen mitteln
 allCov_ws_cr <- t(sapply(results2, FUN = function(iSample) {
   if (!is.null(iSample$ws_cr$varcorr)) {
     return(as.vector(iSample$ws_cr$varcorr))
   } else {
     return(rep(NA, length(as.vector(results2[[1]]$ws_cr$varcorr))))
   }
 }, simplify = TRUE))
 mean_varcov_ws_cr <- apply(allCov_ws_cr, MARGIN = 2, mean, na.rm = TRUE)
 
 
 ### Berechnungen für das "ns_cr"-Modell (cluster robust) ###
 
 # 1. Anteil signifikanter p-Werte
 ns_cr_p_values <- sapply(results2, function(res) res$ns_cr$p_value)
 mean_ns_cr_p <- mean(0.05 > ns_cr_p_values)
 
 # 2. Koeffizienten mitteln
 ns_cr_coef_matrix <- sapply(results2, function(res) res$ns_cr$coefficients)
 mean_ns_cr_coef <- rowMeans(ns_cr_coef_matrix)
 
 # 3. Varianz-Kovarianzen mitteln
 allCov_ns_cr <- t(sapply(results2, FUN = function(iSample) {
   if (!is.null(iSample$ns_cr$varcorr)) {
     return(as.vector(iSample$ns_cr$varcorr))
   } else {
     return(rep(NA, length(as.vector(results2[[1]]$ns_cr$varcorr))))
   }
 }, simplify = TRUE))
 mean_varcov_ns_cr <- apply(allCov_ns_cr, MARGIN = 2, mean, na.rm = TRUE)
 
 ### Zusammenfassen der gemittelten Ergebnisse ###
 result_mean2 <- list(
   ws = list(
     mean_p_value     = mean_ws_p,
     mean_coefficients = mean_ws_coef,
     mean_varcorr     = mean_varcov_ws
   ),
   ns = list(
     mean_p_value     = mean_ns_p,
     mean_coefficients = mean_ns_coef,
     mean_varcorr     = mean_varcov_ns
   ),
   ws_cr = list(
     mean_p_value     = mean_ws_cr_p,
     mean_coefficients = mean_ws_cr_coef,
     mean_varcorr     = mean_varcov_ws_cr
   ),
   ns_cr = list(
     mean_p_value     = mean_ns_cr_p,
     mean_coefficients = mean_ns_cr_coef,
     mean_varcorr     = mean_varcov_ns_cr
   )
 )
 
 # Ausgabe der gemittelten Ergebnisse:
 print(result_mean2) 
 
 
 
 
 
 
 
 
 


 
 
#Bei "bauchiger" Heterosked wird der alpha fehler konservativer, slope fängt das nicht ab, 
#hat keinen Einfluss auf die geschätzte Slope Varianz
 
 
# random slope im daten generierenden aber nicht schätzen -> müsste zu erhöhtem alpha fehler führen
 
 
        #random Intercept                         random slope

 #[1]  0.9937493458 -0.0002335415 -0.0002335415  0.0018448485 ohne heterosked. , wenn random slope var = 0
 
 
 # [1] 1.0051579771 0.0001851261 0.0001851261 0.1187376644 mit heterosked. (^4), wenn random slope var = 0
 #und N = 50
 
 # [1]  1.0035063358 -0.0002438043 -0.0002438043  0.0666080637 mit heterosked. (^4), wenn random slope var = 0
 #und N = 100
 
 # [1] 0.997786140 0.003527678 0.003527678 0.617769948 wenn man random slope auf 0.5 setzt und
 # heterosked. anstellt überschätzt er den slope 
 
 
 # 1.0051765478 0.0003594919 0.0003594919 0.0701559717
 
 # p wert  0.0488 bei N = 100 slope var = 0 und heterosked. (^4)
 # p wert  0.042 bei N = 50 slope var = 0 ohne heterosked. (^4)
 # p wert  0.0456 bei N = 50 slope var = 0 und heterosked. (^2)
 
 
# 0.05366667
 
 
 
 
 