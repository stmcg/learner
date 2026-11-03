test_that('CV stores the exact grids and observed count in both modes', {
  target <- dat_highsim$Y_target[1:8, 1:6]
  source <- dat_highsim$Y_source[1:8, 1:6]
  target[c(2, 9)] <- NA_real_
  for (separate in c(FALSE, TRUE)) {
    grids <- if (separate) {
      list(lambda_1_row_all = c(3, 0), lambda_1_col_all = c(2, 1))
    } else list(lambda_1_all = c(3, 0))
    set.seed(17)
    cv <- do.call(cv.learner, c(list(Y_source = source, Y_target = target,
      r = 2, lambda_2_all = c(.1, .01), step_size = .0003,
      n_folds = 2, control = list(max_iter = 5)), grids))
    expect_equal(cv$n_observed, 46)
    expect_equal(cv$lambda_2_all, c(.1, .01))
    for (name in names(grids)) expect_identical(cv[[name]], grids[[name]])
    auto <- prepare_cv_plot(cv)
    expect_equal(auto$scores, cv$mse_all / 46)
    expect_identical(auto$scale, 'Validation MSE')
    old <- cv[setdiff(names(cv), c(names(grids), 'lambda_2_all', 'n_observed'))]
    manual <- do.call(prepare_cv_plot, c(list(x = old, lambda_2_all = c(.1, .01),
                                           n_observed = 46), grids))
    expect_equal(auto, manual)
    raw <- do.call(prepare_cv_plot, c(list(x = old, lambda_2_all = c(.1, .01)), grids))
    expect_equal(raw$scores, cv$mse_all)
    expect_identical(raw$scale, 'Validation SSE')
  }
})

test_that('explicit plot inputs override metadata and stored inputs are validated', {
  x <- list(mse_all = matrix(c(4, 8, 12, 16), 2), lambda_1_min = 1,
            lambda_2_min = .1, lambda_1_all = c(1, 2),
            lambda_2_all = c(.1, 1), n_observed = 4)
  expect_equal(prepare_cv_plot(x, n_observed = 2)$minimum, 2)
  x$lambda_1_all <- c(1, 2, 3)
  expect_error(prepare_cv_plot(x), 'dimension')
  expect_equal(prepare_cv_plot(x, lambda_1_all = c(1, 2))$minimum, 1)
  x$lambda_1_all <- c(1, 2)
  x$n_observed <- 0
  expect_error(prepare_cv_plot(x), 'n_observed')
  expect_equal(prepare_cv_plot(x, n_observed = 4)$minimum, 1)
})

test_that('automatic MSE scaling is reflected in heatmap labels', {
  skip_if_not_installed('ggplot2', minimum_version = '3.4.0')
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  x <- list(mse_all = matrix(c(4, 4.1, 4.2, 4.3), 2), lambda_1_min = 1,
            lambda_2_min = .1, lambda_1_all = c(1, 2),
            lambda_2_all = c(.1, 1), n_observed = 4)
  z <- plot_cv(x)
  expect_equal(z$minimum, 1)
  expect_match(z$plot$labels$subtitle, 'Minimum MSE = +1')
  expect_match(z$plot$labels$caption, 'MSE range', fixed = TRUE)
  expect_match(plot_cv(x, n_observed = 2)$plot$labels$subtitle,
               'Minimum MSE = +2')
})
