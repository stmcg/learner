# This heatmap demonstration requires the optional ggplot2 package.
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    utils::packageVersion("ggplot2") < "3.4.0") {
  stop('Install ggplot2 >= 3.4.0 to run this example: install.packages("ggplot2")')
}
# Synthetic data demonstrate the interface, not statistical performance.
# Use a small aligned subset to keep this demonstration quick.
sources <- lapply(learner::dat_multisource$Y_source,
                  function(x) x[1:20, 1:10])
target <- learner::dat_multisource$Y_target[1:20, 1:10]
weights <- c(0.7, 0.3)
row_grid <- c(0, 1, 10)
col_grid <- c(0, 1, 3)
balance_grid <- c(0.1, 1)

set.seed(815)
cv <- learner::cv.learner(
  sources, target, r = 3, source_weights = weights,
  lambda_1_row_all = row_grid, lambda_1_col_all = col_grid,
  lambda_2_all = balance_grid, step_size = 0.0003, n_folds = 3,
  control = list(max_iter = 500))
learner::plot_cv(
  cv, lambda_1_row_all = row_grid, lambda_1_col_all = col_grid,
  lambda_2_all = balance_grid, n_observed = sum(!is.na(target)))

fit <- learner::learner(
  sources, target, r = cv$r, source_weights = weights,
  lambda_1_row = cv$lambda_1_row_min, lambda_1_col = cv$lambda_1_col_min,
  lambda_2 = cv$lambda_2_min, step_size = 0.0003,
  control = list(max_iter = 2000))
learner::plot_objective(fit, iterations = "near_optimum")

# Common-penalty results are supported as well.
set.seed(815)
common_grid <- c(0, 1, 10)
common_cv <- learner::cv.learner(
  sources, target, r = 3, source_weights = weights,
  lambda_1_all = common_grid, lambda_2_all = balance_grid,
  step_size = 0.0003, n_folds = 3, control = list(max_iter = 500))
learner::plot_cv(common_cv, lambda_1_all = common_grid,
                 lambda_2_all = balance_grid, n_observed = sum(!is.na(target)))

# Inspect stopping codes as well as curves. Plotting does not prove convergence.
