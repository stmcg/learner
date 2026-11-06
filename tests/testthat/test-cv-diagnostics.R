test_that('saved diagnostics and scores agree with independent fold refits', {
  source <- dat_highsim$Y_source[1:8, 1:6]
  target <- dat_highsim$Y_target[1:8, 1:6]
  target[c(2, 9, 18)] <- NA_real_
  n_folds <- 2L
  for (separate in c(FALSE, TRUE)) {
    for (relative in c(FALSE, TRUE)) {
      control <- list(max_iter = 12, threshold = if (relative) .1 else 0,
                      stopping_rule = if (relative) 'relative' else 'absolute',
                      window = if (relative) 3 else 1)
      args <- if (separate) {
        list(lambda_1_row_all = c(2, 0), lambda_1_col_all = c(1, 0))
      } else list(lambda_1_all = c(2, 0))
      set.seed(91)
      cv <- do.call(cv.learner, c(list(Y_source = list(source, source * .9),
        Y_target = target, source_weights = c(.7, .3), r = 2,
        lambda_2_all = c(.1, .01), n_folds = n_folds,
        step_size = .001, control = control), args))
      d <- cv$cv_diagnostics
      expect_equal(nrow(d), length(cv$mse_all) * n_folds)
      expect_identical(d$candidate, rep(seq_along(cv$mse_all), each = n_folds))
      expect_identical(d$fold, rep(seq_len(n_folds), length(cv$mse_all)))
      expect_identical(d$selected, d$candidate == which.min(cv$mse_all))
      set.seed(91)
      idx <- sample(which(!is.na(target)))
      n <- length(idx)
      folds <- lapply(seq_len(n_folds), function(k) {
        idx[seq.int(floor((k - 1) * n / n_folds) + 1, floor(k * n / n_folds))]
      })
      expected_scores <- numeric(length(cv$mse_all))
      for (i in seq_len(nrow(d))) {
        held_out <- folds[[d$fold[i]]]
        training <- target
        training[held_out] <- NA_real_
        fit <- learner(list(source, source * .9), training,
          source_weights = c(.7, .3), r = 2,
          lambda_1_row = d$lambda_1_row[i], lambda_1_col = d$lambda_1_col[i],
          lambda_2 = d$lambda_2[i], step_size = .001, control = control)
        expect_equal(d$stop_code[i], fit$convergence_criterion)
        expect_equal(d$updates[i], length(fit$objective_values))
        g <- d$candidate[i]
        expected_scores[g] <- expected_scores[g] +
          sum((target[held_out] - fit$learner_estimate[held_out])^2)
      }
      expect_equal(as.vector(cv$mse_all), expected_scores, tolerance = 1e-12)
      if (!separate) expect_equal(d$lambda_1_row, d$lambda_1_col)
    }
  }
})

test_that('non-finite candidates record failure and mark skipped folds as missing', {
  source <- diag(c(4, 3, 2, 1))
  target <- diag(c(1, 2, 3, 4))
  set.seed(42)
  expect_warning(cv <- cv.learner(source, target, r = 2,
    lambda_1_all = c(1e308, 0), lambda_2_all = 1, step_size = .001,
    n_folds = 3, control = list(max_iter = 3)), 'non-finite validation errors')
  d <- cv$cv_diagnostics
  expect_equal(d$stop_code[1:3], c(3L, NA_integer_, NA_integer_))
  expect_equal(d$updates[1:3], c(1L, NA_integer_, NA_integer_))
  expect_false(any(d$selected[1:3]))
  expect_false(anyNA(d$stop_code[4:6]))
  expect_true(all(d$selected[4:6]))
})

test_that('CV records tolerance stopping including the singleton grid', {
  source <- matrix(0, 4, 3)
  set.seed(12)
  cv <- cv.learner(source, source, r = 1, lambda_1_all = 0,
    lambda_2_all = 0, step_size = .001, n_folds = 2,
    control = list(max_iter = 10))
  expect_equal(cv$cv_diagnostics$stop_code, c(1L, 1L))
  expect_equal(cv$cv_diagnostics$updates, c(2L, 2L))
  expect_true(all(cv$cv_diagnostics$selected))
})

test_that('diagnostics match refits when stopping differs across candidates and folds', {
  sources <- dat_multisource$Y_source
  target <- dat_multisource$Y_target
  weights <- c(0.7, 0.3)
  settings <- list(
    list(n_folds = 3L, args = list(lambda_1_all = c(0, 1, 10))),
    list(n_folds = 2L, args = list(lambda_1_row_all = c(0, 1, 10),
                                   lambda_1_col_all = c(0, 1, 10))))
  for (setting in settings) {
    n_folds <- setting$n_folds
    set.seed(815)
    cv <- do.call(cv.learner, c(list(Y_source = sources, Y_target = target,
      r = 3, source_weights = weights, lambda_2_all = c(0.1, 1),
      n_folds = n_folds, step_size = 0.0003, control = list(max_iter = 2000)),
      setting$args))
    d <- cv$cv_diagnostics
    # The rows must differ, or a row stored in the wrong place would go unnoticed.
    expect_gt(nrow(unique(d[c('stop_code', 'updates')])), 3L)
    expect_gt(length(unique(d$stop_code)), 1L)
    set.seed(815)
    idx <- sample(which(!is.na(target)))
    n <- length(idx)
    folds <- lapply(seq_len(n_folds), function(k) {
      idx[seq.int(floor((k - 1) * n / n_folds) + 1, floor(k * n / n_folds))]
    })
    for (i in seq_len(nrow(d))) {
      training <- target
      training[folds[[d$fold[i]]]] <- NA_real_
      fit <- learner(sources, training, r = 3, source_weights = weights,
        lambda_1_row = d$lambda_1_row[i], lambda_1_col = d$lambda_1_col[i],
        lambda_2 = d$lambda_2[i], step_size = 0.0003,
        control = list(max_iter = 2000))
      expect_equal(d$stop_code[i], fit$convergence_criterion)
      expect_equal(d$updates[i], length(fit$objective_values))
    }
  }
})
