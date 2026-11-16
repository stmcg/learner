auto_fixture <- function() {
  set.seed(191)
  list(sources = list(matrix(rnorm(48), 8, 6), matrix(rnorm(48), 8, 6)),
       target = matrix(rnorm(48), 8, 6))
}

# Independent reference: construct full projectors as in the paper/helper.
reference_weights <- function(sources, target, r) {
  target[is.na(target)] <- 0
  b0 <- svd(target, nu = r, nv = r)
  d <- vapply(sources, function(x) {
    b <- svd(x, nu = r, nv = r)
    sum((tcrossprod(b$u) - tcrossprod(b0$u))^2) +
      sum((tcrossprod(b$v) - tcrossprod(b0$v))^2)
  }, numeric(1))
  w <- 1 / (d + 1e-12)
  unname(w / sum(w))
}

test_that('automatic weights match Equation 3 with complete and missing targets', {
  d <- auto_fixture()
  for (missing in c(FALSE, TRUE)) {
    y <- d$target
    if (missing) y[c(1, 9, 35)] <- NA
    expected <- reference_weights(d$sources, y, 2)
    fit <- learner(d$sources, y, 2, 1, 1, 0.01, control = list(max_iter = 4))
    expect_equal(fit$source_weights, expected, tolerance = 1e-12)
    expect_equal(fit$source_weight_method, 'auto')
    fixed <- learner(d$sources, y, 2, 1, 1, 0.01, source_weights = expected,
                     control = list(max_iter = 4))
    expect_equal(fit$learner_estimate, fixed$learner_estimate, tolerance = 1e-12)
  }
  direct <- dlearner(d$sources, d$target, 2)
  expected <- reference_weights(d$sources, d$target, 2)
  expect_equal(direct$source_weights, expected, tolerance = 1e-12)
  expect_equal(direct$dlearner_estimate,
               dlearner(d$sources, d$target, 2, expected)$dlearner_estimate,
               tolerance = 1e-12)
})

test_that('weights handle coincident, orthogonal, and full-rank spaces', {
  a <- diag(c(4, 0, 0, 0)); b <- diag(c(0, 3, 0, 0))
  fit <- dlearner(list(a, b), a, r = 1)
  expect_gt(fit$source_weights[1], 1 - 1e-10)
  expect_true(all(is.finite(fit$source_weights)))
  expect_equal(sum(fit$source_weights), 1)
  expect_equal(dlearner(list(a, 2*a), a, 1)$source_weights, c(.5, .5))
  expect_equal(dlearner(list(a, b), a, 4)$source_weights, c(.5, .5))
  expect_equal(dlearner(a, a, 1)$source_weights, 1)
  expect_equal(dlearner(list(a, b), a, 1, 'equal')$source_weights, c(.5, .5))
  expect_equal(dlearner(list(a, b), a, 1, NULL)$source_weights, c(.5, .5))
  expect_equal(dlearner(list(a, b), a, 1, c(2, 1))$source_weights, c(2/3, 1/3))
  expect_equal(dlearner(list(b, a), a, 1)$source_weights, rev(fit$source_weights))
})

test_that('automatic CV scores and weights match independent fold-specific fits', {
  d <- auto_fixture(); d$target[c(1, 10)] <- NA
  set.seed(92)
  cv <- cv.learner(d$sources, d$target, r = 2, lambda_1_row_all = c(0, 2),
                  lambda_1_col_all = c(1, 3), lambda_2_all = c(0, .1),
                  step_size = .01, n_folds = 3,
                  control = list(max_iter = 9, threshold = 0))
  expected <- manual_cv_scores(d$sources, d$target, c(0, 2), c(1, 3),
                              c(0, .1), weights = 'auto')
  expect_equal(cv$mse_all, expected, tolerance = 1e-10)
  expect_equal(cv$source_weights, reference_weights(d$sources, d$target, 2),
               tolerance = 1e-12)
  set.seed(92)
  observed <- which(!is.na(d$target))
  indices <- sample(observed, length(observed), replace = FALSE)
  for (fold in 1:3) {
    held <- indices[seq.int(floor((fold-1)*length(indices)/3)+1,
                            floor(fold*length(indices)/3))]
    train <- d$target; train[held] <- NA
    expect_equal(cv$fold_source_weights[fold, ],
                 reference_weights(d$sources, train, 2), tolerance = 1e-12)
    # Changing held-out values cannot affect that fold's training weights.
    changed <- d$target; changed[held] <- changed[held] + 1e4
    set.seed(92)
    perturbed <- cv.learner(d$sources, changed, 2, 1, 1, .01, n_folds = 3,
                            control = list(max_iter = 2))
    expect_equal(perturbed$fold_source_weights[fold, ],
                 cv$fold_source_weights[fold, ], tolerance = 1e-12)
  }
})

test_that('automatic weighting does not change RNG and rejects invalid modes', {
  d <- auto_fixture(); seed <- .Random.seed
  dlearner(d$sources, d$target, 2)
  expect_identical(.Random.seed, seed)
  for (mode in list('automatic', NA_character_, c('auto','equal'), character())) {
    expect_error(dlearner(d$sources, d$target, 2, mode), 'source_weights')
  }
})

test_that('single-source weighting helper needs no target decomposition', {
  d <- auto_fixture()
  bases <- list(svd(d$sources[[1]], nu = 2, nv = 2))
  # A single source must receive weight one even without observed target data;
  # the helper should return before attempting a target SVD.
  expect_identical(source_subspace_weights(bases, d$target * NA_real_, 2), 1)
})

test_that('CV kernel rejects incompatible fold-weight matrix dimensions', {
  d <- auto_fixture()
  call_kernel <- function(weights) {
    cv_learner_cpp(
      Y_source = d$sources, Y_target = d$target,
      source_weights = c(.5, .5), row_grid = 1, col_grid = 1,
      balance_grid = 1, step_size = .01, max_iter = 2L,
      threshold = 0, n_cores = 1L, r = 2L, max_value = 10,
      index_set = list(0:23, 24:47), relative = FALSE, window = 1L,
      fold_source_weights = weights)
  }
  expect_error(call_kernel(matrix(.5, 1, 2)),
               'fold_source_weights has incompatible dimensions', fixed = TRUE)
  expect_error(call_kernel(matrix(1, 2, 1)),
               'fold_source_weights has incompatible dimensions', fixed = TRUE)
  expect_error(call_kernel(matrix(numeric(), 0, 2)),
               'fold_source_weights has incompatible dimensions', fixed = TRUE)
  expect_no_error(call_kernel(matrix(.5, 2, 2)))
})
