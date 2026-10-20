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

## -----------------------------------------------------------------------------
pi_selection <- 0.5
lam <- 0.06
lambda_path <- lam * exp(seq(log(10), 0, length.out = 30))

select <- function(d) {
  rows <- sample(nrow(d$X), round(pi_selection * nrow(d$X)))
  fit <- glmnet::glmnet(d$X[rows, ], d$y[rows], family = "binomial",
                        lambda = lambda_path, control = list(thresh = 1e-14))
  list(fit = fit, rows = rows)
}
sel <- select(d)
which(coef(sel$fit, s = lam)[-1] != 0)

## -----------------------------------------------------------------------------
li <- glmnet_inference(sel$fit, d$X, d$y, s = lam, selection_rows = sel$rows, level = 0.9)
li$scalar_noise
summ <- li$summary()
summ$variable <- ifelse(summ$index == 0, "(Intercept)", paste0("x", summ$index))
summ

## -----------------------------------------------------------------------------
target <- function(li) drop(li$W %*% (li$Q_hat %*% c(0, beta))[li$E + 1])
summ$truth <- target(li)
summ$cover <- summ$lower_conf <= summ$truth & summ$truth <= summ$upper_conf
summ[, c("variable", "lower_conf", "upper_conf", "p_value", "truth", "cover")]

## -----------------------------------------------------------------------------
splitting <- function(d, rows, li, level = 0.9) {
  held_out <- setdiff(seq_len(nrow(d$X)), rows)
  X1 <- cbind(1, d$X)[held_out, li$E + 1, drop = FALSE]
  fit <- glm.fit(X1, d$y[held_out], family = binomial())
  est <- unname(coef(fit))
  mu <- fit$fitted.values
  sd <- sqrt(diag(solve(crossprod(X1, X1 * (mu * (1 - mu))))))
  q <- qnorm(1 - (1 - level) / 2)
  data.frame(index = li$E, estimate = est, lower_conf = est - q * sd, upper_conf = est + q * sd)
}
li$`_naive`
splitting(d, sel$rows, li)

## -----------------------------------------------------------------------------
one_replication <- function() {
  d <- make_data(n)
  sel <- select(d)
  li <- glmnet_inference(sel$fit, d$X, d$y, s = lam, selection_rows = sel$rows, level = 0.9)
  summ <- li$summary()
  keep <- summ$index != 0   # report the selected variables, not the intercept
  truth <- target(li)[keep]
  covers <- function(df) (df$lower_conf <= truth & truth <= df$upper_conf)
  lengths <- function(df) df$upper_conf - df$lower_conf
  tabs <- list(carving = summ[keep, ], naive = li$`_naive`[keep, ],
               splitting = splitting(d, sel$rows, li)[keep, ])
  data.frame(method = rep(names(tabs), each = sum(keep)),
             cover = unlist(lapply(tabs, covers)), length = unlist(lapply(tabs, lengths)),
             signal = rep(beta[summ$index[keep]] != 0, length(tabs)))
}

set.seed(1)
results <- do.call(rbind, replicate(100, one_replication(), simplify = FALSE))
data.frame(coverage = tapply(results$cover, results$method, mean),
           median_length = tapply(results$length, results$method, median))

## -----------------------------------------------------------------------------
with(results, tapply(cover, list(method, signal = signal), mean))

