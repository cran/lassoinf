## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      eval = requireNamespace("glmnet", quietly = TRUE))

## -----------------------------------------------------------------------------
library(lassoinf)
set.seed(0)

n <- 600; p <- 30
beta <- c(0.6, -0.6, 0.4, rep(0, p - 3))
make_data <- function(n) {
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, plogis(drop(X %*% beta)))
  list(X = X, y = y)
}
d <- make_data(n)

lam <- 0.03
lambda_path <- lam * exp(seq(log(10), 0, length.out = 30))
fit <- glmnet::glmnet(d$X, d$y, family = "binomial", lambda = lambda_path, control = list(thresh = 1e-14))

## -----------------------------------------------------------------------------
li <- glmnet_inference(fit, d$X, d$y, s = lam, level = 0.9)
summ <- li$summary()
summ$variable <- ifelse(summ$index == 0, "(Intercept)", paste0("x", summ$index))
summ

## -----------------------------------------------------------------------------
target <- function(li) drop(li$W %*% (li$Q_hat %*% c(0, beta))[li$E + 1])
summ$truth <- target(li)
summ$cover <- summ$lower_conf <= summ$truth & summ$truth <= summ$upper_conf
summ[, c("variable", "lower_conf", "upper_conf", "p_value", "truth", "cover")]

## -----------------------------------------------------------------------------
li$`_naive`

## -----------------------------------------------------------------------------
prob <- glmnet_problem_from_fit(fit, d$X, d$y, s = lam)
Z <- -prob$G_hat + drop(prob$Q_hat %*% prob$beta_hat)
li_manual <- LassoInference$new(prob$beta_hat, prob$G_hat, prob$Q_hat, prob$D, prob$L, prob$U,
                                Z_full = Z, Sigma = prob$Q_hat / n, scalar_noise = 0,
                                level = 0.9)
all.equal(li_manual$summary(), li$summary())

## -----------------------------------------------------------------------------
one_replication <- function() {
  d <- make_data(n)
  fit <- glmnet::glmnet(d$X, d$y, family = "binomial", lambda = lambda_path, control = list(thresh = 1e-14))
  li <- glmnet_inference(fit, d$X, d$y, s = lam, level = 0.9)
  summ <- li$summary()
  keep <- summ$index != 0   # report the selected variables, not the intercept
  truth <- target(li)[keep]
  tabs <- list(polyhedral = summ[keep, ], naive = li$`_naive`[keep, ])
  data.frame(method = rep(names(tabs), each = sum(keep)),
             cover = unlist(lapply(tabs, function(df) df$lower_conf <= truth & truth <= df$upper_conf)),
             length = unlist(lapply(tabs, function(df) df$upper_conf - df$lower_conf)),
             signal = rep(beta[summ$index[keep]] != 0, length(tabs)))
}

set.seed(1)
results <- do.call(rbind, replicate(100, one_replication(), simplify = FALSE))
with(results, tapply(cover, list(method, signal = signal), mean))
with(results, tapply(length, list(method, signal = signal), median))

