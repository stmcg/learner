library(learner)
sources <- dat_multisource$Y_source
target <- dat_multisource$Y_target
set.seed(815)
cv <- cv.learner(sources, target, r = 3, source_weights = NULL,
                 lambda_1_all = c(0, 1, 10), lambda_2_all = c(.1, 1),
                 step_size = .0003, n_folds = 3,
                 control = list(max_iter = 2000))
cv$fold_source_weights # Training-only weights, one row per fold.
cv$source_weights      # All-observed-data weights for final refitting.
fit <- learner(sources, target, r = cv$r, source_weights = NULL,
               lambda_1 = cv$lambda_1_min, lambda_2 = cv$lambda_2_min,
               step_size = .0003, control = list(max_iter = 2000))
stopifnot(isTRUE(all.equal(fit$source_weights, cv$source_weights)))
fit$convergence_criterion
cv$cv_diagnostics
# Supply equal numeric weights to give each source the same influence.
direct <- dlearner(sources, target, r = 3, source_weights = c(1, 1))
direct$source_weights
