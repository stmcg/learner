penalty_example <- function() {
  set.seed(421)
  list(source = matrix(rnorm(20), 5, 4),
       target = matrix(rnorm(20), 5, 4))
}

test_that("equal separate penalties preserve common and positional calls", {
  dat <- penalty_example()
  for (with_missing in c(FALSE, TRUE)) {
    target <- dat$target
    if (with_missing) target[c(2, 7)] <- NA_real_
    control <- list(max_iter = 12, threshold = 0)
    common <- learner(dat$source, target, 1, 3, 1, 0.01, control)
    separate <- learner(dat$source, target, r = 1, lambda_2 = 1,
                        step_size = 0.01, control = control,
                        lambda_1_row = 3, lambda_1_col = 3)
    expect_equal(separate, common, tolerance = 1e-12)
  }
})

test_that("transposing data swaps the two independent space penalties", {
  dat <- penalty_example()
  fit <- function(source, target, row, col) {
    learner(source, target, r = 1, lambda_1_row = row,
            lambda_1_col = col, lambda_2 = 1, step_size = 0.02,
            control = list(max_iter = 12, threshold = 0))
  }
  original <- fit(dat$source, dat$target, 3, 0.25)
  transposed <- fit(t(dat$source), t(dat$target), 0.25, 3)
  expect_equal(original$learner_estimate, t(transposed$learner_estimate),
               tolerance = 1e-9)
  expect_equal(original$objective_values, transposed$objective_values,
               tolerance = 1e-9)
})

# A finite-difference oracle checks the C++ updates against the objective,
# without copying the analytic gradient formulas used by the implementation.
numeric_penalty_fit <- function(source, target, row, col, balance,
                                step_size, iterations, weights = NULL) {
  sources <- if (is.matrix(source)) list(source) else source
  bases <- lapply(sources, svd, nu = 1, nv = 1)
  s <- bases[[1]]
  if (is.null(weights)) weights <- rep(1, length(bases))
  weights <- weights / sum(weights)
  p <- nrow(sources[[1]])
  q <- ncol(sources[[1]])
  pu <- Reduce(`+`, Map(function(b, w) w * tcrossprod(b$u), bases, weights))
  pv <- Reduce(`+`, Map(function(b, w) w * tcrossprod(b$v), bases, weights))
  u <- s$u * sqrt(s$d[1])
  v <- s$v * sqrt(s$d[1])
  theta <- c(u, v)
  observed <- !is.na(target)
  objective <- function(x) {
    U <- matrix(x[seq_len(p)], p, 1)
    V <- matrix(x[p + seq_len(q)], q, 1)
    residual <- (U %*% t(V) - target)[observed]
    sum(residual^2) / mean(observed) +
      row * sum(((diag(p) - pu) %*% U)^2) +
      col * sum(((diag(q) - pv) %*% V)^2) +
      balance * (sum(U^2) - sum(V^2))^2
  }
  best <- theta
  best_value <- objective(theta)
  values <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gradient <- vapply(seq_along(theta), function(j) {
      direction <- numeric(length(theta))
      direction[j] <- 1e-6
      (objective(theta + direction) - objective(theta - direction)) / 2e-6
    }, numeric(1))
    for (indices in list(seq_len(p), p + seq_len(q))) {
      theta[indices] <- theta[indices] - step_size *
        sqrt(sum(theta[indices]^2)) /
        (sqrt(sum(gradient[indices]^2)) + 1e-12) * gradient[indices]
    }
    values[i] <- objective(theta)
    if (values[i] < best_value) {
      best <- theta
      best_value <- values[i]
    }
  }
  list(estimate = outer(best[seq_len(p)], best[p + seq_len(q)]),
       objective = values)
}

test_that("unequal and zero penalties follow the intended objective", {
  dat <- penalty_example()
  for (with_missing in c(FALSE, TRUE)) {
    target <- dat$target
    if (with_missing) target[c(2, 7)] <- NA_real_
    for (penalties in list(c(3, 0.25), c(0.25, 3), c(0, 3), c(3, 0))) {
      reference <- numeric_penalty_fit(dat$source, target, penalties[1],
                                       penalties[2], 1, 0.02, 12)
      actual <- learner(dat$source, target, r = 1,
                        lambda_1_row = penalties[1],
                        lambda_1_col = penalties[2], lambda_2 = 1,
                        step_size = 0.02,
                        control = list(max_iter = 12, threshold = 0,
                                       max_value = 1e6))
      expect_equal(actual$objective_values, reference$objective,
                   tolerance = 1e-6)
      expect_equal(actual$learner_estimate, reference$estimate,
                   tolerance = 1e-6)
    }
  }
})

test_that("invalid or ambiguous space penalties fail clearly", {
  dat <- penalty_example()
  fit <- function(...) learner(dat$source, dat$target, r = 1,
                               lambda_2 = 1, step_size = 0.01, ...)
  expect_error(fit(), "Provide lambda_1")
  expect_error(fit(lambda_1_row = 1), "provided together")
  expect_error(fit(lambda_1_col = 1), "provided together")
  expect_error(fit(lambda_1 = 1, lambda_1_row = 1, lambda_1_col = 1),
               "not both modes")
  expect_error(fit(lambda_1 = 1, lambda_1_row = 1), "not both modes")
  invalid <- list(numeric(0), c(1, 2), NA_real_, NaN, Inf, -1, "1", TRUE)
  for (value in invalid) {
    expect_error(fit(lambda_1 = value), "lambda_1 must be")
    expect_error(fit(lambda_1_row = value, lambda_1_col = 1),
                 "lambda_1_row must be")
    expect_error(fit(lambda_1_row = 1, lambda_1_col = value),
                 "lambda_1_col must be")
  }
})


test_that("multiple-source updates match finite differences of the weighted objective", {
  dat <- penalty_example()
  set.seed(25)
  sources <- list(dat$source, matrix(rnorm(20), 5, 4), matrix(rnorm(20), 5, 4))
  for (masked in c(FALSE, TRUE)) {
    y <- dat$target
    if (masked) y[c(2, 7)] <- NA_real_
    reference <- numeric_penalty_fit(sources, y, 3, 0.25, 1, 0.02, 12,
                                     weights = c(0.2, 0.3, 0.5))
    actual <- learner(sources, y, r = 1, lambda_1_row = 3, lambda_1_col = 0.25,
                      lambda_2 = 1, step_size = 0.02, source_weights = c(2, 3, 5),
                      control = list(max_iter = 12, threshold = 0, max_value = 1e6))
    expect_equal(actual$objective_values, reference$objective, tolerance = 1e-6)
    expect_equal(actual$learner_estimate, reference$estimate, tolerance = 1e-6)
  }
})
