# Small, deterministic matrices keep error-path tests independent of fixtures.
failure_example <- function() {
  list(source = diag(c(4, 3, 2, 1)), target = diag(c(1, 2, 3, 4)))
}

test_that('all fitting interfaces reject malformed target matrices', {
  d <- failure_example()
  infinite <- d$target
  infinite[1, 1] <- Inf
  calls <- list(
    function(y) learner(d$source, y, r = 2, lambda_1 = 1,
                         lambda_2 = 1, step_size = 0.001),
    function(y) cv.learner(d$source, y, r = 2, lambda_1_all = 1,
                            lambda_2_all = 1, step_size = 0.001),
    function(y) dlearner(d$source, y, r = 2))
  for (target in list(as.vector(d$target), matrix('x', 4, 4),
                      matrix(numeric(), 0, 4), infinite)) {
    for (fit in calls) {
      expect_error(fit(target), 'Y_target must be a nonempty numeric matrix',
                   fixed = TRUE)
    }
  }
})

test_that('missing-target messages distinguish an empty target from D-LEARNER requirements', {
  d <- failure_example()
  missing <- d$target
  missing[] <- NA_real_
  expect_error(learner(d$source, missing, r = 2, lambda_1 = 1,
                        lambda_2 = 1, step_size = 0.001),
               '^Y_target must contain observed values\\.$')
  expect_error(cv.learner(d$source, missing, r = 2, lambda_1_all = 1,
                           lambda_2_all = 1, step_size = 0.001),
               '^Y_target must contain observed values\\.$')
  d$target[1, 1] <- NA_real_
  expect_error(dlearner(d$source, d$target, r = 2),
               'dlearner requires a complete Y_target without missing values', fixed = TRUE)
})

test_that('learner and CV reject invalid optimization settings', {
  d <- failure_example()
  calls <- list(
    function(...) learner(d$source, d$target, r = 2, lambda_1 = 1,
                           lambda_2 = 1, ...),
    function(...) cv.learner(d$source, d$target, r = 2, lambda_1_all = 1,
                              lambda_2_all = 1, ...))
  for (fit in calls) {
    expect_error(fit(step_size = 0.001, control = 3),
                 'control must be a list', fixed = TRUE)
    for (bad in c(0, -1, Inf, NA_real_)) {
      expect_error(fit(step_size = bad),
                   'step_size must be a finite positive number', fixed = TRUE)
      expect_error(fit(step_size = 0.001, control = list(max_value = bad)),
                   'control$max_value must be a finite positive number', fixed = TRUE)
    }
  }
})

test_that('a non-finite update stops immediately and retains the finite initialization', {
  d <- failure_example()
  # A finite but extreme step overflows the objective on the first update.
  fit <- learner(d$source, d$target, r = 2, lambda_1 = 1, lambda_2 = 1,
                  step_size = 1e200, control = list(max_iter = 3))
  expect_equal(fit$convergence_criterion, 3)
  expect_length(fit$objective_values, 1)
  expect_false(is.finite(fit$objective_values[1]))
  expect_true(all(is.finite(fit$learner_estimate)))
  expect_equal(fit$learner_estimate, diag(c(4, 3, 0, 0)), tolerance = 1e-12)
})

test_that('CV reports when no candidate has a finite validation score', {
  d <- failure_example()
  set.seed(42)
  expect_error(
    cv.learner(d$source, d$target, r = 2, lambda_1_all = c(0, 1),
                lambda_2_all = 1, step_size = 1e200, n_folds = 2,
                control = list(max_iter = 3)),
    'No candidate produced finite validation errors', fixed = TRUE)
})

test_that('CV excludes a failed candidate and preserves the valid candidate score', {
  d <- failure_example()
  cv <- function(grid) {
    set.seed(42)
    cv.learner(d$source, d$target, r = 2, lambda_1_all = grid,
                lambda_2_all = 1, step_size = 0.001, n_folds = 2,
                control = list(max_iter = 3))
  }
  reference <- cv(0)
  # An extreme finite penalty causes non-finite updates in just that candidate.
  expect_warning(result <- cv(c(1e308, 0)),
                 'Some candidates produced non-finite validation errors and were excluded',
                 fixed = TRUE)
  expect_equal(result$mse_all[1, 1], Inf)
  expect_true(is.finite(result$mse_all[2, 1]))
  expect_equal(result$mse_all[2, 1], reference$mse_all[1, 1], tolerance = 1e-12)
  expect_equal(result$lambda_1_min, reference$lambda_1_min)
  expect_equal(result$lambda_2_min, reference$lambda_2_min)
})
