# Demonstration using the built-in simulated multiple-source data.
# Each population has independently generated noise; source A is more precise.
sources <- learner::dat_multisource$Y_source
target <- learner::dat_multisource$Y_target
# Illustrative weights, not weights estimated or tuned from these data.
weights <- c(0.7, 0.3)

# Tune the U-space, V-space, and balance penalties separately.
set.seed(815)
cv_result <- learner::cv.learner(
  Y_source = sources, Y_target = target, r = 3,
  source_weights = weights,
  lambda_1_row_all = c(1, 10),
  lambda_1_col_all = c(1, 3),
  lambda_2_all = c(0.1, 1),
  step_size = 0.0003, n_folds = 3,
  control = list(max_iter = 2000)
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
  step_size = 0.0003, control = list(max_iter = 5000)
)
print(fit$learner_estimate[1:5, 1:5])
print(fit$convergence_criterion)
# 1: objective-change threshold; 2: iteration limit; 3: excessive growth/breakdown.
# This example demonstrates the workflow, not scientific performance or convergence.

# D-LEARNER also accepts multiple sources with the same weights.
direct_fit <- learner::dlearner(sources, target, r = 3, source_weights = weights)
print(dim(direct_fit$dlearner_estimate))

# Compare source choices using D-LEARNER at the same fixed rank.
# D-LEARNER avoids optimization and penalty tuning as additional differences.
# Truth is available only because this is a simulated dataset.
source_choices <- list(
  A_only = sources$A,
  B_only = sources$B,
  Weighted_A_B = sources
)
comparison <- do.call(rbind, lapply(names(source_choices), function(label) {
  estimate <- learner::dlearner(
    source_choices[[label]], target, r = 3,
    source_weights = if (label == 'Weighted_A_B') weights else NULL
  )$dlearner_estimate
  data.frame(source = label,
             frobenius_error = norm(estimate - learner::dat_multisource$Theta_target, 'F'))
}))
print(comparison, row.names = FALSE)
# Lower Frobenius error is better. These are single-dataset results, not CV MSE.
# The weights 0.7/0.3 are illustrative; combining sources need not beat source A.

# A second setting: both sources have appreciable noise.
# Keep the built-in data unchanged; variances are target 0.10, A 0.25, B 0.50.
# Reuse the true signals, with new independent noise for every population.
set.seed(20261016)
add_noise <- function(signal, variance) {
  signal + matrix(rnorm(length(signal), sd = sqrt(variance)),
                  nrow(signal), ncol(signal))
}
target_noisy <- add_noise(learner::dat_multisource$Theta_target, variance = 0.10)
sources_noisy <- list(
  A = add_noise(learner::dat_multisource$Theta_source$A, variance = 0.25),
  B = add_noise(learner::dat_multisource$Theta_source$B, variance = 0.50)
)
estimates_noisy <- list(
  Target_only = learner::dlearner(target_noisy, target_noisy, r = 3)$dlearner_estimate,
  A_only = learner::dlearner(sources_noisy$A, target_noisy, r = 3)$dlearner_estimate,
  B_only = learner::dlearner(sources_noisy$B, target_noisy, r = 3)$dlearner_estimate,
  Weighted_A_B = learner::dlearner(sources_noisy, target_noisy, r = 3,
                         source_weights = c(0.7, 0.3))$dlearner_estimate
)
errors_noisy <- vapply(estimates_noisy, function(estimate) {
  norm(estimate - learner::dat_multisource$Theta_target, "F")
}, numeric(1))
print(round(errors_noisy, 3))
# Here the weighted combination improves on either source alone, but the
# target-only SVD has the smallest error. Fixed weights and one noise draw
# illustrate a possible outcome, not a general performance guarantee.
