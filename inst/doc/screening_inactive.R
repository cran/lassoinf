## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      eval = requireNamespace("glmnet", quietly = TRUE))

## -----------------------------------------------------------------------------
library(lassoinf)

n <- 100; p <- 300; k <- 15
pi_selection <- 0.8
scalar_noise <- (1 - pi_selection) / pi_selection   # Var(omega) = scalar_noise * Sigma

make_data <- function() {
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p)
  beta[1:k] <- 3 / sqrt(n) * sample(c(-1, 1), k, replace = TRUE)
  mu <- drop(X %*% beta)
  y <- mu + rnorm(n)
  y_noisy <- y + rnorm(n, sd = sqrt(scalar_noise))
  list(X = X, beta = beta, mu = mu, y = y, y_noisy = y_noisy)
}

set.seed(0)
d <- make_data()

## -----------------------------------------------------------------------------
lam <- 0.8 * sqrt(2 * n * log(p))

lasso_inference <- function(d) {
  fit <- glmnet::glmnet(d$X, d$y_noisy, lambda = lam / n, intercept = FALSE,
                        standardize = FALSE, control = list(thresh = 1e-14))
  beta_hat <- as.numeric(coef(fit))[-1]
  beta_hat[abs(beta_hat) < 1e-6] <- 0
  G_hat <- drop(crossprod(d$X, d$X %*% beta_hat - d$y_noisy))
  Q <- xtvx_operator(d$X)
  LassoInference$new(beta_hat, G_hat, Q, rep(lam, p), NULL, NULL,
                     drop(crossprod(d$X, d$y)), Q, NULL, scalar_noise, level = 0.9)
}

li <- lasso_inference(d)
E <- li$E + 1   # the package reports 0-based indices
cat(length(E), "active variables,", sum(d$beta[E] != 0), "of them true signals;",
    sum(d$beta[-E] != 0), "true signals were left out\n")

## -----------------------------------------------------------------------------
screen <- ScreenedSelection$new(li, top_k = 10)
screen$screened + 1
screen$screened_signs
c(first_dropped = screen$first_dropped + 1, sign = screen$first_dropped_sign)
c(lasso_rows = li$A$rows(), screening_rows = screen$A_screen$rows())

## -----------------------------------------------------------------------------
coords <- SelectionCoordinates$new(li)
Z_mu <- drop(crossprod(d$X, d$mu))

screened_df <- inactive_summary(screen)
screened_df$truth <- sapply(screened_df$index, function(j) estimand_value(coords$inactive_coef(j), Z_mu))
screened_df$signal <- d$beta[screened_df$index + 1] != 0
screened_df

## -----------------------------------------------------------------------------
naive_p_value <- function(li, est) {
  sd <- sqrt(li$si$compute_contrast(est$eta)$naive_variance)
  2 * pnorm(-abs(estimand_value(est, li$Z_full)) / sd)
}

data.frame(variable = screen$screened + 1,
           screened = inactive_summary(screen)$p_value,
           lasso_only = inactive_summary(li, variables = screen$screened)$p_value,
           naive = sapply(screen$screened, function(j) naive_p_value(li, coords$inactive_coef(j))),
           signal = d$beta[screen$screened + 1] != 0)

## -----------------------------------------------------------------------------
head(inactive_summary(screen, basis = "lasso"))

## -----------------------------------------------------------------------------
exact <- ScreenedSelection$new(li, top_k = 10, conditioning = "exact")
c(first_dropped = screen$A_screen$rows(), exact = exact$A_screen$rows())
cols <- c("lower_conf", "upper_conf", "p_value")
cbind(index = screen$screened,
      first_dropped = inactive_summary(screen)[, cols],
      exact = inactive_summary(exact)[, cols])

## -----------------------------------------------------------------------------
one_replication <- function() {
  d <- make_data()
  li <- lasso_inference(d)
  screen <- ScreenedSelection$new(li, top_k = 10)
  coords <- SelectionCoordinates$new(li)
  Z_mu <- drop(crossprod(d$X, d$mu))
  q <- qnorm(0.95)
  rows <- list()
  for (j in screen$screened) {
    est <- coords$inactive_coef(j)
    target <- estimand_value(est, Z_mu)
    sd <- sqrt(li$si$compute_contrast(est$eta)$naive_variance)
    theta <- estimand_value(est, li$Z_full)
    signal <- d$beta[j + 1] != 0
    for (method in c("screened", "lasso only")) {
      res <- contrast_inference(if (method == "screened") screen else li, est)
      rows[[length(rows) + 1]] <- data.frame(method = method,
        cover = res$lower_conf <= target && target <= res$upper_conf,
        length = res$upper_conf - res$lower_conf, reject = res$p_value < 0.1, signal = signal)
    }
    rows[[length(rows) + 1]] <- data.frame(method = "naive", cover = abs(theta - target) <= q * sd,
      length = 2 * q * sd, reject = abs(theta) > q * sd, signal = signal)
  }
  do.call(rbind, rows)
}

set.seed(1)
results <- do.call(rbind, replicate(50, one_replication(), simplify = FALSE))
data.frame(coverage = tapply(results$cover, results$method, mean),
           median_length = tapply(results$length, results$method, median))

## -----------------------------------------------------------------------------
with(results, tapply(reject, list(method, signal = signal), mean))

