# Demonstration only: these are simulated matrices, not actual population data.
# The second source is constructed here to illustrate the multiple-source API.
set.seed(803)
source_A <- learner::dat_highsim$Y_source
source_B <- source_A + matrix(rnorm(length(source_A), sd = 0.3),
                              nrow(source_A), ncol(source_A))
target <- learner::dat_highsim$Y_target
sources <- list(source_A, source_B)
weights <- c(0.7, 0.3)

# Tune the U-space, V-space, and balance penalties separately.
set.seed(815)
cv_result <- learner::cv.learner(
  Y_source = sources, Y_target = target, r = 3,
  source_weights = weights,
  lambda_1_row_all = c(1, 10),
  lambda_1_col_all = c(1, 3),
  lambda_2_all = c(0.1, 1),
  step_size = 0.003, n_folds = 3,
  control = list(max_iter = 100)
)
print(c(row = cv_result$lambda_1_row_min,
        col = cv_result$lambda_1_col_min,
        balance = cv_result$lambda_2_min))
print(dim(cv_result$mse_all)) # row grid x column grid x balance grid

# Refit with the selected penalties using all observed target entries.
fit <- learner::learner(
  Y_source = sources, Y_target = target, r = cv_result$r,
  source_weights = weights,
  lambda_1_row = cv_result$lambda_1_row_min,
  lambda_1_col = cv_result$lambda_1_col_min,
  lambda_2 = cv_result$lambda_2_min,
  step_size = 0.003, control = list(max_iter = 100)
)
print(fit$learner_estimate[1:5, 1:5])
print(fit$convergence_criterion)
# 1: objective-change threshold; 2: iteration limit; 3: excessive growth/breakdown.
# This short example tests the workflow, not scientific performance or convergence.

# D-LEARNER also accepts multiple sources with the same weights.
direct_fit <- learner::dlearner(sources, target, r = 3, source_weights = weights)
print(dim(direct_fit$dlearner_estimate))
