multi_example <- function() {
  set.seed(705)
  list(a = matrix(rnorm(35), 7, 5), b = matrix(rnorm(35), 7, 5),
       y = matrix(rnorm(35), 7, 5))
}

# Rebuild the held-out predictions explicitly to check grid axes and masking.
manual_cv_scores <- function(sources, target, row, col, balance, weights = NULL,
                              seed = 92, folds = 3) {
  set.seed(seed)
  observed <- which(!is.na(target))
  indices <- sample(observed, length(observed), replace = FALSE)
  groups <- lapply(seq_len(folds), function(fold) {
    indices[seq.int(floor((fold - 1) * length(indices) / folds) + 1,
                    floor(fold * length(indices) / folds))]
  })
  scores <- array(NA_real_, c(length(row), length(col), length(balance)))
  for (i in seq_along(row)) for (j in seq_along(col)) for (k in seq_along(balance)) {
    scores[i, j, k] <- sum(vapply(groups, function(held_out) {
      train <- target
      train[held_out] <- NA_real_
      fit <- learner(sources, train, r = 2, lambda_1_row = row[i],
                     lambda_1_col = col[j], lambda_2 = balance[k],
                     source_weights = weights, step_size = 0.01,
                     control = list(max_iter = 9, threshold = 0))
      sum((fit$learner_estimate[held_out] - target[held_out])^2)
    }, numeric(1)))
  }
  scores
}

test_that("independent CV axes and scores match explicit held-out predictions", {
  d <- multi_example()
  for (multi in c(FALSE, TRUE)) for (missing in c(FALSE, TRUE)) {
    sources <- if (multi) list(d$a, d$b) else d$a
    weights <- if (multi) c(2, 1) else NULL
    y <- d$y
    if (missing) y[c(3, 10, 14)] <- NA_real_
    row <- c(0, 3)
    col <- c(0.2, 2, 5)
    balance <- c(0, 1)
    reference <- manual_cv_scores(sources, y, row, col, balance, weights)
    set.seed(92)
    actual <- cv.learner(sources, y, r = 2, lambda_1_row_all = row,
                         lambda_1_col_all = col, lambda_2_all = balance,
                         source_weights = weights, step_size = 0.01, n_folds = 3,
                         control = list(max_iter = 9, threshold = 0))
    expect_equal(actual$mse_all, reference, tolerance = 1e-10)
    best <- arrayInd(which.min(reference), dim(reference))
    expect_equal(actual$lambda_1_row_min, row[best[1]])
    expect_equal(actual$lambda_1_col_min, col[best[2]])
    expect_equal(actual$lambda_2_min, balance[best[3]])
    expect_true(is.na(actual$lambda_1_min))
  }
})

test_that("tied CV is the diagonal of the independent grid", {
  d <- multi_example()
  fit <- function(...) {
    set.seed(18)
    cv.learner(list(d$a, d$b), d$y, r = 2, lambda_2_all = c(0, 1),
               step_size = 0.01, control = list(max_iter = 8), ...)
  }
  tied <- fit(lambda_1_all = c(1, 3))
  separate <- fit(lambda_1_row_all = c(1, 3), lambda_1_col_all = c(1, 3))
  expect_named(tied, c('lambda_1_min', 'lambda_2_min', 'mse_all', 'r'))
  expect_equal(dim(tied$mse_all), c(2L, 2L))
  for (i in 1:2) expect_equal(tied$mse_all[i, ], separate$mse_all[i, i, ], tolerance = 1e-10)
  single <- fit(lambda_1_row_all = 1, lambda_1_col_all = 3)
  expect_equal(dim(single$mse_all), c(1L, 1L, 2L))
  parallel <- fit(lambda_1_row_all = c(1, 3), lambda_1_col_all = c(1, 3), n_cores = 2)
  expect_equal(parallel, separate, tolerance = 1e-12)
})

test_that("source representations and normalized weights are consistent", {
  d <- multi_example()
  fit <- function(source, weights = NULL) {
    learner(source, d$y, r = 2, lambda_1_row = 3, lambda_1_col = 0.5,
            lambda_2 = 1, step_size = 0.01, source_weights = weights,
            control = list(max_iter = 12, threshold = 0))
  }
  expect_equal(fit(d$a), fit(list(d$a)), tolerance = 1e-12)
  expect_equal(fit(d$a), fit(list(d$a, d$a)), tolerance = 1e-10)
  expect_equal(fit(d$a), fit(list(d$a, d$b), c(1, 0)), tolerance = 1e-10)
  expect_equal(fit(list(d$a, d$b), c(2, 1)),
               fit(list(d$a, d$b), c(2/3, 1/3)), tolerance = 1e-12)
  expect_equal(fit(list(d$a, d$b)), fit(list(d$a, d$b), c(1, 1)), tolerance = 1e-12)
  expect_gt(max(abs(fit(list(d$a, d$b), c(1, 0))$learner_estimate -
                      fit(list(d$a, d$b), c(1, 3))$learner_estimate)), 1e-6)
})

test_that("D-LEARNER matches weighted projectors and preserves matrix calls", {
  d <- multi_example()
  ua <- svd(d$a, nu = 2, nv = 2)
  ub <- svd(d$b, nu = 2, nv = 2)
  pu <- 0.7 * tcrossprod(ua$u) + 0.3 * tcrossprod(ub$u)
  pv <- 0.7 * tcrossprod(ua$v) + 0.3 * tcrossprod(ub$v)
  result <- dlearner(list(d$a, d$b), d$y, 2, c(7, 3))
  expect_equal(result$dlearner_estimate, pu %*% d$y %*% pv, tolerance = 1e-12)
  expect_equal(dlearner(d$a, d$y, 2), dlearner(list(d$a), d$y, 2), tolerance = 1e-12)
  expect_equal(dlearner(d$a, d$y, 2), dlearner(list(d$a, d$b), d$y, 2, c(1, 0)), tolerance = 1e-12)
  expect_equal(result, dlearner(list(d$b, d$a), d$y, 2, c(3, 7)), tolerance = 1e-12)
})

test_that("invalid sources, weights, ranks, folds, and grids are rejected", {
  d <- multi_example()
  fit <- function(source = list(d$a, d$b), weights = NULL, rank = 2) {
    learner(source, d$y, r = rank, source_weights = weights,
            lambda_1 = 1, lambda_2 = 1, step_size = 0.01)
  }
  expect_error(fit(list()), 'nonempty list')
  expect_error(fit(list(d$a, 2)), 'numeric matrix')
  expect_error(fit(list(d$a, d$b[1:3, ])), 'same dimensions')
  bad <- d$b; bad[1] <- NA_real_
  expect_error(fit(list(d$a, bad)), 'finite values')
  for (w in list(c(0, 0), c(-1, 2), c(1, NA), c(1, Inf), 1, 'x')) {
    expect_error(fit(weights = w), 'source_weights')
  }
  for (r in list(0, 1.5, 6, NA, Inf, c(1, 2))) expect_error(fit(rank = r), 'r must')
  cv <- function(...) cv.learner(d$a, d$y, r = 2, step_size = 0.01, ...)
  expect_error(cv(lambda_1_row_all = 1, lambda_2_all = 1), 'provided together')
  expect_error(cv(lambda_1_all = 1, lambda_1_row_all = 1,
                  lambda_1_col_all = 1, lambda_2_all = 1), 'not both modes')
  for (grid in list(numeric(), NA, Inf, -1, '1')) {
    expect_error(cv(lambda_1_all = grid, lambda_2_all = 1), 'nonempty vector')
  }
  expect_error(cv(lambda_1_all = 1, lambda_2_all = numeric()), 'lambda_2_all')
  for (folds in c(1, 36, 2.5, NA)) {
    expect_error(cv(lambda_1_all = 1, lambda_2_all = 1, n_folds = folds), 'n_folds')
  }
  y <- matrix(NA_real_, 7, 5)
  expect_error(learner(d$a, y, 2, 1, 1, 0.01), 'observed values')
  y[1] <- 1
  expect_error(cv.learner(d$a, y, 2, 1, 1, 0.01), 'n_folds')
})
