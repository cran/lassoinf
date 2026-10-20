## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      eval = requireNamespace("glmnet", quietly = TRUE))

## -----------------------------------------------------------------------------
library(lassoinf)

n <- 100; p <- 300; k <- 15
pi_selection <- 0.8
scalar_noise <- (1 - pi_selection) / pi_selection   # Var(omega) = scalar_noise * Sigma
signals <- seq(1, p, by = p / k)

make_data <- function(rho, amp) {
  # AR(1) columns: corr(x_i, x_l) = rho^|i - l|; rho = 0 gives independent columns
  X <- matrix(rnorm(n * p), n, p)
  for (i in 2:p) X[, i] <- rho * X[, i - 1] + sqrt(1 - rho^2) * X[, i]
  beta <- numeric(p)
  beta[signals] <- amp / sqrt(n) * sample(c(-1, 1), k, replace = TRUE)
  mu <- drop(X %*% beta)
  y <- mu + rnorm(n)
  y_noisy <- y + rnorm(n, sd = sqrt(scalar_noise))
  list(X = X, beta = beta, mu = mu, y = y, y_noisy = y_noisy)
}

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

## -----------------------------------------------------------------------------
swap_estimands <- function(li) {
  E <- li$E + 1
  Q_EE <- solve(li$W)
  Q_EE <- (Q_EE + t(Q_EE)) / 2
  cache <- new.env()
  column <- function(j) {
    key <- as.character(j)
    if (is.null(cache[[key]])) {
      e_j <- numeric(p); e_j[j] <- 1
      cache[[key]] <- drop(li$Q_hat$multiply(e_j))
    }
    cache[[key]]
  }
  # F = (E minus {s}) plus {j}, with j last
  model <- function(s, j) c(E[E != s], j)
  swap_coef <- function(s, j) {
    if (!(s %in% E) || j %in% E) stop("need s in E and j not in E")
    Q_j <- column(j)
    keep <- E != s
    F <- model(s, j)
    m <- length(F)
    Q_FF <- matrix(0, m, m)
    Q_FF[-m, -m] <- Q_EE[keep, keep]
    Q_FF[-m, m] <- Q_FF[m, -m] <- Q_j[F[-m]]
    Q_FF[m, m] <- Q_j[j]
    eta <- numeric(p)
    eta[F] <- solve(Q_FF, c(numeric(m - 1), 1))
    contrast_estimand(eta)
  }
  list(model = model, swap_coef = swap_coef)
}

## -----------------------------------------------------------------------------
swap_table <- function(selection, swaps, selected, screened) {
  pairs <- expand.grid(added = screened, dropped = selected)[, c("dropped", "added")]
  estimands <- Map(swaps$swap_coef, pairs$dropped, pairs$added)
  names(estimands) <- paste(pairs$dropped, pairs$added, sep = ":")
  df <- estimand_summary(selection, estimands)
  list(table = cbind(pairs, df[, -1]), estimands = estimands)
}

analyze <- function(d) {
  li <- lasso_inference(d)
  screen <- ScreenedSelection$new(li, top_k = 10)
  swaps <- swap_estimands(li)
  E <- li$E + 1; screened <- screen$screened + 1
  cat(length(E), "active variables,", sum(d$beta[E] != 0), "of them true signals;",
      sum(d$beta[-E] != 0), "true signals were left out\n")
  cat("screened:", screened, "\n")
  st <- swap_table(screen, swaps, E, screened)
  tab <- st$table
  tab$truth <- sapply(st$estimands, estimand_value, Z = drop(crossprod(d$X, d$mu)))
  tab$dropped_signal <- d$beta[tab$dropped] != 0
  tab$added_signal <- d$beta[tab$added] != 0
  tab$adjacent <- abs(tab$dropped - tab$added) <= 2
  rownames(tab) <- NULL
  list(li = li, screen = screen, swaps = swaps, E = E, screened = screened,
       table = tab, estimands = st$estimands)
}

p_value_grid <- function(tab) round(xtabs(p_value ~ dropped + added, tab), 3)

## -----------------------------------------------------------------------------
naive_p_value <- function(li, est) {
  sd <- sqrt(li$si$compute_contrast(est$eta)$naive_variance)
  2 * pnorm(-abs(estimand_value(est, li$Z_full)) / sd)
}

compare_p_values <- function(a) {
  lasso_only <- swap_table(a$li, a$swaps, a$E, a$screened)$table
  data.frame(a$table[, c("dropped", "added")],
             screened = a$table$p_value,
             lasso_only = lasso_only$p_value,
             naive = sapply(a$estimands, naive_p_value, li = a$li), row.names = NULL)
}

## -----------------------------------------------------------------------------
one_replication <- function(rho, amp) {
  d <- make_data(rho, amp)
  li <- lasso_inference(d)
  screen <- ScreenedSelection$new(li, top_k = 10)
  swaps <- swap_estimands(li)
  Z_mu <- drop(crossprod(d$X, d$mu))
  q <- qnorm(0.95)
  rows <- list()
  for (s in li$E + 1) for (j in screen$screened + 1) {
    est <- swaps$swap_coef(s, j)
    target <- estimand_value(est, Z_mu)
    sd <- sqrt(li$si$compute_contrast(est$eta)$naive_variance)
    theta <- estimand_value(est, li$Z_full)
    info <- data.frame(dropped_signal = d$beta[s] != 0, added_signal = d$beta[j] != 0,
                       adjacent = abs(s - j) <= 2)
    for (method in c("screened", "lasso only")) {
      res <- contrast_inference(if (method == "screened") screen else li, est)
      rows[[length(rows) + 1]] <- data.frame(method = method,
        cover = res$lower_conf <= target && target <= res$upper_conf,
        length = res$upper_conf - res$lower_conf, reject = res$p_value < 0.1, info)
    }
    rows[[length(rows) + 1]] <- data.frame(method = "naive", cover = abs(theta - target) <= q * sd,
      length = 2 * q * sd, reject = abs(theta) > q * sd, info)
  }
  do.call(rbind, rows)
}

simulate <- function(rho, amp, nsim, seed = 1) {
  set.seed(seed)
  do.call(rbind, replicate(nsim, one_replication(rho, amp), simplify = FALSE))
}

coverage_table <- function(results) {
  data.frame(coverage = tapply(results$cover, results$method, mean),
             median_length = tapply(results$length, results$method, median))
}

## ----include = FALSE----------------------------------------------------------
# The simulations take a few minutes, so their results are saved in
# surrogate_screening_sims.rds (simulate() sets its own seed, so they are the same).
# Set LASSOINF_RERUN_SIMS=true and rebuild the vignette to recompute them.
sims_file <- "surrogate_screening_sims.rds"
if (!file.exists(sims_file) || identical(Sys.getenv("LASSOINF_RERUN_SIMS"), "true")) {
  saveRDS(list(indep = simulate(rho = 0, amp = 3, nsim = 50),
               ar = simulate(rho = 0.8, amp = 5, nsim = 30)), sims_file)
}
sims <- readRDS(sims_file)
results_indep <- sims$indep
results_ar <- sims$ar

## -----------------------------------------------------------------------------
set.seed(0)
d <- make_data(rho = 0, amp = 3)
a <- analyze(d)

## -----------------------------------------------------------------------------
s <- a$E[1]; j <- a$screened[1]
F <- a$swaps$model(s, j)
c(ols = unname(tail(qr.solve(d$X[, F], d$y), 1)),
  swap = estimand_value(a$swaps$swap_coef(s, j), a$li$Z_full))

## -----------------------------------------------------------------------------
head(a$table, 10)

## -----------------------------------------------------------------------------
p_value_grid(a$table)

## -----------------------------------------------------------------------------
head(compare_p_values(a), 10)

## ----eval = -1----------------------------------------------------------------
results_indep <- simulate(rho = 0, amp = 3, nsim = 50)
coverage_table(results_indep)
with(results_indep, tapply(reject, list(method, added_signal = added_signal), mean))

## -----------------------------------------------------------------------------
set.seed(0)
d <- make_data(rho = 0.8, amp = 5)
a <- analyze(d)

## -----------------------------------------------------------------------------
a$table[a$table$adjacent, ]

## -----------------------------------------------------------------------------
p_value_grid(a$table)

## -----------------------------------------------------------------------------
surrogate_pairs <- a$table[a$table$adjacent & !a$table$dropped_signal & a$table$added_signal, ]
surrogate_pairs$median_other_rows <- sapply(seq_len(nrow(surrogate_pairs)), function(i) {
  col <- a$table[a$table$added == surrogate_pairs$added[i], ]
  median(col$p_value[col$dropped != surrogate_pairs$dropped[i]])
})
surrogate_pairs[, c("dropped", "added", "p_value", "median_other_rows")]

## -----------------------------------------------------------------------------
head(compare_p_values(a), 10)

## ----eval = -1----------------------------------------------------------------
results_ar <- simulate(rho = 0.8, amp = 5, nsim = 30)
coverage_table(results_ar)
with(results_ar, tapply(reject, list(method, dropped_signal = dropped_signal,
                                     added_signal = added_signal), mean))

## -----------------------------------------------------------------------------
surrogate <- subset(results_ar, adjacent & !dropped_signal & added_signal)
data.frame(pairs = tapply(surrogate$reject, surrogate$method, length),
           coverage = tapply(surrogate$cover, surrogate$method, mean),
           reject = tapply(surrogate$reject, surrogate$method, mean))

