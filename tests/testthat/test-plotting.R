test_that('objective plots retain the update indices and restore graphics state', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(3, 3, 2, 1))
  old <- graphics::par(c('mfrow', 'mar'))
  fit <- list(objective_values = c(9, 5, 6, 5.5), convergence_criterion = 2L)
  original <- fit
  result <- withVisible(plot_objective(fit, last = 2))
  expect_false(result$visible)
  expect_equal(result$value$iteration, 1:4)
  expect_equal(result$value$last_indices, 3:4)
  expect_equal(result$value$recorded_min_index, 2L)
  expect_identical(fit, original)
  expect_equal(graphics::par(c('mfrow', 'mar')), old)
  expect_equal(plot_objective(fit, last = 100)$last_indices, 1:4)
})

test_that('objective plots support singleton, flat, and failed histories', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  for (values in list(4, c(4, 4, 4), c(4, Inf, NaN))) {
    fit <- list(objective_values = values, convergence_criterion = 3L)
    result <- plot_objective(fit, last = 1)
    expect_equal(result$recorded_min_index, 1L)
    expect_equal(result$objective, values)
  }
  expect_error(plot_objective(list(objective_values = c(Inf, NA))), 'No finite')
  expect_error(plot_objective(list()), 'objective_values')
  fit <- list(objective_values = 4, convergence_criterion = 1L)
  expect_error(plot_objective(fit, last = 0), 'last')
  expect_error(plot_objective(fit, last = 1.5), 'last')
  expect_error(plot_objective(fit, main = NA_character_), 'main')
  fit$convergence_criterion <- 9L
  expect_error(plot_objective(fit), 'convergence_criterion')
})

test_that('CV plots retain candidate order, mark the global minimum, and scale SSE', {
  testthat::skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  x <- cv_plot_example()
  original <- x
  old <- graphics::par(c('mfrow','mar'))
  args <- list(x = x, lambda_1_row_all = c(10, 0),
               lambda_1_col_all = c(0.2, 3), lambda_2_all = c(0, 1))
  sse <- do.call(plot_cv, args)
  mse <- do.call(plot_cv, c(args, list(n_observed = 40)))
  expect_identical(sse$grids$row, c(10, 0))
  expect_identical(sse$best_index, c(2L, 2L, 2L))
  expect_equal(sse$scores, x$mse_all)
  expect_equal(mse$scores, x$mse_all / 40)
  expect_equal(sse$scale, 'Validation SSE')
  expect_equal(mse$scale, 'Validation MSE')
  expect_equal(mse$color_limits, range(x$mse_all) / 40)
  expect_identical(x, original)
  expect_equal(graphics::par(c('mfrow','mar')), old)
  subset <- do.call(plot_cv, c(args, list(lambda_2_index = 1)))
  expect_equal(subset$color_limits, sse$color_limits)
  expect_equal(subset$shown_lambda_2_index, 1)
  expect_equal(subset$best_index, sse$best_index)
  common <- plot_cv(cv_plot_example(FALSE), lambda_1_all = c(10, 0), lambda_2_all = c(0, 1))
  expect_identical(common$best_index, c(2L, 2L))
})

test_that('CV plots handle singleton axes, duplicate candidates, ties and unavailable scores', {
  testthat::skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  x <- list(mse_all = array(c(5, 5, Inf, NA), c(1, 2, 2)),
            lambda_1_row_min = 0, lambda_1_col_min = 2, lambda_2_min = 1)
  z <- plot_cv(x, lambda_1_row_all = 0, lambda_1_col_all = c(2, 2),
               lambda_2_all = c(1, 10))
  expect_identical(z$best_index, c(1L, 1L, 1L))
  expect_equal(z$unavailable_count, 2)
  expect_gt(diff(z$color_limits), 0)
  x <- list(mse_all = matrix(0, 1, 1), lambda_1_min = 0, lambda_2_min = 0)
  z <- plot_cv(x, lambda_1_all = 0, lambda_2_all = 0)
  expect_identical(z$best_index, c(1L, 1L))
  expect_gt(diff(z$color_limits), 0)
  # Seven balance slices exercise pagination without inventing new candidates.
  x <- list(mse_all = array(7:1, c(1, 1, 7)), lambda_1_row_min = 0,
            lambda_1_col_min = 2, lambda_2_min = 7)
  z <- plot_cv(x, lambda_1_row_all = 0, lambda_1_col_all = 2, lambda_2_all = 1:7)
  expect_equal(z$shown_lambda_2_index, 1:7)
  expect_identical(z$best_index, c(1L, 1L, 7L))
})

test_that('CV plots reject missing, mismatched, and ambiguous grid information', {
  testthat::skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  x <- cv_plot_example()
  args <- list(x = x, lambda_1_row_all = c(10, 0),
               lambda_1_col_all = c(0.2, 3), lambda_2_all = c(0, 1))
  expect_error(plot_cv(x), 'lambda_1_row_all')
  bad <- args; bad$lambda_1_row_all <- c(0,10)
  expect_error(do.call(plot_cv,bad), 'selected penalties')
  bad <- args; bad$lambda_1_row_all <- 0
  expect_error(do.call(plot_cv,bad), 'score-array dimension')
  expect_error(do.call(plot_cv,c(args,list(lambda_1_all=1))), 'Do not supply')
  expect_error(do.call(plot_cv,c(args,list(n_observed=0))), 'n_observed')
  expect_error(do.call(plot_cv,c(args,list(lambda_2_index=c(1,1)))), 'distinct')
  expect_error(do.call(plot_cv,c(args,list(lambda_2_index=3))), 'distinct')
  expect_error(do.call(plot_cv,c(args,list(lambda_2_index=NA_real_))), 'distinct')
  x$mse_all[] <- Inf
  expect_error(plot_cv(x), 'No finite')
  x$mse_all[] <- -1
  expect_error(plot_cv(x), 'cannot be negative')
  expect_error(plot_cv(list(mse_all=1:3)), 'two- or three-dimensional')
})

test_that('plot functions accept actual fits and both CV modes without changing outputs', {
  testthat::skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  set.seed(623)
  s <- matrix(rnorm(48), 8, 6)
  t <- matrix(rnorm(48), 8, 6); t[c(1,3)] <- NA_real_
  fit <- learner(s,t,r=2,lambda_1=1,lambda_2=1,step_size=0.001,
                 control=list(max_iter=6))
  expect_equal(plot_objective(fit)$objective,fit$objective_values)
  for (separate in c(FALSE,TRUE)) {
    args <- list(Y_source=list(s,s*0.9),Y_target=t,r=2,step_size=0.001,
                 lambda_2_all=c(0,1),n_folds=2,control=list(max_iter=6))
    grid <- if(separate) list(lambda_1_row_all=c(0,1),lambda_1_col_all=c(0.5,2)) else list(lambda_1_all=c(0,1))
    cv <- do.call(cv.learner,c(args,grid))
    original <- cv
    result <- do.call(plot_cv,c(list(x=cv,lambda_2_all=c(0,1),n_observed=sum(!is.na(t))),grid))
    expect_equal(result$scores,cv$mse_all/46)
    expect_identical(cv,original)
  }
})


test_that('paper-style log labels preserve scores and candidate selection', {
  testthat::skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  x <- list(mse_all = matrix(c(4, 3, 2, 1), 2, 2),
            lambda_1_min = 100, lambda_2_min = 1)
  raw <- plot_cv(x, lambda_1_all = c(10, 100), lambda_2_all = c(0.01, 1))
  log <- plot_cv(x, lambda_1_all = c(10, 100), lambda_2_all = c(0.01, 1), log_labels = TRUE)
  expect_equal(log$scores, raw$scores)
  expect_equal(log$best_index, raw$best_index)
  expect_equal(log$grids, raw$grids)
  expect_true(log$log_labels)
  expect_error(plot_cv(x, lambda_1_all = c(0, 100), lambda_2_all = c(0.01, 1),
                       log_labels = TRUE), 'strictly positive')
  expect_error(plot_cv(x, log_labels = NA), 'TRUE or FALSE')
  expect_error(plot_cv(x, log_labels = 'yes'), 'TRUE or FALSE')
})
