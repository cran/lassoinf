## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      eval = requireNamespace("glmnet", quietly = TRUE))

## -----------------------------------------------------------------------------
library(lassoinf)
set.seed(0)

n <- 100; p <- 20
signal_strength <- 4
X <- matrix(rnorm(n * p), n, p)
beta <- numeric(p)
beta[1:3] <- c(2, 1, -1) * signal_strength / sqrt(n)
sigma2 <- 2
y <- drop(X %*% beta) + rnorm(n, sd = sqrt(sigma2))

## -----------------------------------------------------------------------------
pi_selection <- 0.9
Q <- crossprod(X)
Z_full <- drop(crossprod(X, y))
sigma2_hat <- sum(lm.fit(X, y)$residuals^2) / (n - p)
Sigma <- sigma2_hat * Q
sigma2_rand <- (1 - pi_selection) / pi_selection * sigma2_hat
Sigma_noise <- sigma2_rand * Q
y_noisy <- y + rnorm(n, sd = sqrt(sigma2_rand))
c(sigma2_hat = sigma2_hat, sigma2_rand = sigma2_rand)

## -----------------------------------------------------------------------------
lam <- 2 * sqrt(n)
fit <- glmnet::glmnet(X, y_noisy, lambda = lam / n, intercept = FALSE,
                      standardize = FALSE, control = list(thresh = 1e-14))
beta_hat <- as.numeric(coef(fit))[-1]
G_hat <- drop(crossprod(X, X %*% beta_hat - y_noisy))
which(beta_hat != 0)

## -----------------------------------------------------------------------------
li <- LassoInference$new(beta_hat, G_hat, Q, rep(lam, p), NULL, NULL,
                         Z_full, Sigma, Sigma_noise, level = 0.9)
summ <- li$summary()
summ

## -----------------------------------------------------------------------------
Z_mu <- drop(Q %*% beta)
summ$truth <- sapply(as.character(summ$index), function(j) sum(li$contrasts[[j]] * Z_mu))
summ$cover <- summ$lower_conf <= summ$truth & summ$truth <= summ$upper_conf
summ

## -----------------------------------------------------------------------------
li$`_naive`
li$`_splitting`

## -----------------------------------------------------------------------------
one_replication <- function() {
  X <- matrix(rnorm(n * p), n, p)
  y <- drop(X %*% beta) + rnorm(n, sd = sqrt(sigma2))
  Q <- crossprod(X)
  sigma2_hat <- sum(lm.fit(X, y)$residuals^2) / (n - p)
  sigma2_rand <- (1 - pi_selection) / pi_selection * sigma2_hat
  y_noisy <- y + rnorm(n, sd = sqrt(sigma2_rand))
  fit <- glmnet::glmnet(X, y_noisy, lambda = lam / n, intercept = FALSE,
                        standardize = FALSE, control = list(thresh = 1e-14))
  beta_hat <- as.numeric(coef(fit))[-1]
  if (all(beta_hat == 0)) return(NULL)
  G_hat <- drop(crossprod(X, X %*% beta_hat - y_noisy))
  li <- LassoInference$new(beta_hat, G_hat, Q, rep(lam, p), NULL, NULL,
                           drop(crossprod(X, y)), sigma2_hat * Q, sigma2_rand * Q, level = 0.9)
  summ <- li$summary()
  truth <- sapply(as.character(summ$index), function(j) sum(li$contrasts[[j]] * drop(Q %*% beta)))
  covers <- function(df) df$lower_conf <= truth & truth <= df$upper_conf
  data.frame(method = rep(c("randomized", "naive", "splitting"), each = nrow(summ)),
             cover = c(covers(summ), covers(li$`_naive`), covers(li$`_splitting`)),
             length = c(summ$upper_conf - summ$lower_conf,
                        li$`_naive`$upper_conf - li$`_naive`$lower_conf,
                        li$`_splitting`$upper_conf - li$`_splitting`$lower_conf))
}

set.seed(1)
results <- do.call(rbind, replicate(200, one_replication(), simplify = FALSE))
data.frame(coverage = tapply(results$cover, results$method, mean),
           median_length = tapply(results$length, results$method, median))

