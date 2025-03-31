# Lade notwendige Pakete
library(rstan)
library(MASS)  # Falls nötig, um `rmvnorm` für die Datensimulation zu nutzen

# 1. Simuliere die Daten (aus den vorherigen Schritten)
set.seed(222)
n <- 50
beta <- c(10, 0, 0, 1)
sigma_e <- 10

# Korrelation der Prädiktoren X
r_x <- 0.3
sigma_x <- matrix(r_x, nrow = length(beta) - 1, ncol = length(beta) - 1)
diag(sigma_x) <- 1
X <- cbind(1, mvrnorm(n = n, mu = rep(0, length(beta) - 1), Sigma = sigma_x))
y <- X %*% beta + rnorm(n, sd = sigma_e)

# In ein Dataframe für den Überblick
daten <- data.frame(X[,-1], y)

# 2. Bereite die Datenliste für Stan vor
dat <- list(
  N = nrow(X),
  Y = as.vector(y),
  K = ncol(X),
  X = X,
  Kc = ncol(X) - 1,
  prior_only = 0
)

# 3. Speichere den Stan-Code als Datei (z.B. "model.stan")
stan_code <- "
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma2;  // variance parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  real sigma = sqrt(sigma2);  // standard deviation as square root of variance
  
  lprior += normal_lpdf(b | 0, 1e+05);
  lprior += normal_lpdf(Intercept | 0, 1e+05);
  lprior += inv_gamma_lpdf(sigma2 | 1e-04, 1e-04);  // Inverse-Gamma prior on variance
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
"
writeLines(stan_code, "model25.6.stan")

# 4. Führe das Modell aus
fit <- stan(file = "model25.6.stan", data = dat, iter = 20000, chains = 4)

# 5. Zusammenfassen und Ergebnisse anzeigen
print(fit)
plot(fit)


