# Run after installing the plotting development version of learner.
# Synthetic data demonstrate the interface, not statistical performance.
set.seed(803)
source_A <- learner::dat_highsim$Y_source[1:20, 1:10]
source_B <- source_A + matrix(rnorm(length(source_A), sd = 0.3), 20, 10)
target <- learner::dat_highsim$Y_target[1:20, 1:10]
sources <- list(source_A, source_B)
weights <- c(0.7, 0.3)
row_grid <- c(0, 1, 10)
col_grid <- c(0, 1, 3)
balance_grid <- c(0.1, 1)

set.seed(815)
cv <- learner::cv.learner(
  sources, target, r = 2, source_weights = weights,
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
learner::plot_objective(fit, last = 100)

# Common-penalty results are supported as well.
set.seed(815)
common_grid <- c(0, 1, 10)
common_cv <- learner::cv.learner(
  sources, target, r = 2, source_weights = weights,
  lambda_1_all = common_grid, lambda_2_all = balance_grid,
  step_size = 0.0003, n_folds = 3, control = list(max_iter = 500))
learner::plot_cv(common_cv, lambda_1_all = common_grid,
                 lambda_2_all = balance_grid, n_observed = sum(!is.na(target)))

# Inspect stopping codes as well as curves. Plotting does not prove convergence.
