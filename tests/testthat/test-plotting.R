test_that('objective plots show the requested iterations', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  # Decrease to the minimum at iteration 31, then values that stay above it.
  fit <- list(objective_values = c(seq(100, 10, length.out = 30), rep(c(9, 9.5), 50)),
              convergence_criterion = 2L)
  original <- fit
  result <- withVisible(plot_objective(fit))
  expect_false(result$visible)
  expect_equal(result$value$iteration, 1:130)
  expect_equal(result$value$objective, fit$objective_values)
  expect_equal(plot_objective(fit, iterations = 2:4)$iteration, 2:4)
  expect_equal(plot_objective(fit, iterations = c(5, 2, 2))$iteration, c(2, 5))
  # From the smallest objective value (iteration 31) through the last iteration.
  expect_equal(plot_objective(fit, iterations = 'near_optimum')$iteration, 31:130)
  expect_identical(fit, original)
})

test_that('near-optimum plots share the horizontal axis of the full history', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  fit <- list(objective_values = c(seq(100, 10, length.out = 30), rep(c(9, 9.5), 50)),
              convergence_criterion = 2L)
  plot_objective(fit)
  usr_all <- graphics::par('usr')[1:2]
  plot_objective(fit, iterations = 'near_optimum')
  expect_equal(graphics::par('usr')[1:2], usr_all)
  # A history still decreasing at the end shows its final stretch with its own axis.
  fit <- list(objective_values = seq(100, 10, length.out = 100), convergence_criterion = 1L)
  plot_objective(fit)
  usr_all <- graphics::par('usr')[1:2]
  plot_objective(fit, iterations = 'near_optimum')
  expect_gt(graphics::par('usr')[1], usr_all[1] + 50)
})

test_that('near-optimum plots show the final stretch when the minimum is at the end', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  # A history that is still decreasing at the last iteration.
  fit <- list(objective_values = seq(100, 10, length.out = 100), convergence_criterion = 1L)
  expect_equal(plot_objective(fit, iterations = 'near_optimum')$iteration, 91:100)
  # Short histories are shown in full.
  fit <- list(objective_values = c(9, 5, 4), convergence_criterion = 1L)
  expect_equal(plot_objective(fit, iterations = 'near_optimum')$iteration, 1:3)
})

test_that('near-optimum plots handle alternating and missing values and keep graphics settings', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  # Alternating values after a decrease; near_optimum draws them as a band.
  values <- c(seq(50, 11, length.out = 40), 10 + rep(c(0, 0.5), 200))
  values[100] <- NA
  fit <- list(objective_values = values, convergence_criterion = 2L)
  old <- graphics::par('mar')
  result <- plot_objective(fit)
  expect_equal(result$iteration, seq_along(values))
  expect_equal(plot_objective(fit, iterations = 'near_optimum')$iteration, 41:440)
  expect_equal(graphics::par('mar'), old)
})

test_that('objective plots handle flat, singleton, and non-finite histories', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  for (values in list(4, c(4, 4, 4), c(4, Inf, NaN), c(Inf, 3, 1, NA))) {
    fit <- list(objective_values = values, convergence_criterion = 2L)
    expect_equal(plot_objective(fit)$objective, values)
    expect_true(length(plot_objective(fit, iterations = 'near_optimum')$iteration) >= 1L)
  }
  fit <- list(objective_values = c(2, 1, 3), convergence_criterion = 2L)
  for (bad in list(0, 4, 1.5, NA, numeric(), 'near', TRUE)) {
    expect_error(plot_objective(fit, iterations = bad), 'iterations')
  }
  expect_error(plot_objective(list(objective_values = c(Inf, NA))), 'learner')
  expect_error(plot_objective(list()), 'learner')
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

test_that('objective plots reject ranges with no finite values before drawing', {
  path <- plot_device()
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  fit <- list(objective_values = c(3, 1, NA_real_, Inf, NaN))
  old <- graphics::par(c('mar', 'usr'))
  expect_error(
    withCallingHandlers(plot_objective(fit, iterations = 3:5),
                        warning = function(w) stop('Unexpected warning: ', conditionMessage(w))),
    'selected iterations contain no finite objective values')
  expect_equal(graphics::par(c('mar', 'usr')), old)
  expect_equal(plot_objective(fit, iterations = 2:5)$objective,
               fit$objective_values[2:5])
})
