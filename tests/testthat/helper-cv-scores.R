# Rebuild the held-out predictions explicitly to check grid axes and masking.
manual_cv_scores <- function(sources, target, row, col, balance, weights = NULL,
                              seed = 92, folds = 3) {
  set.seed(seed)
  observed <- which(!is.na(target))
  indices <- sample(observed, length(observed), replace = FALSE)
  groups <- lapply(seq_len(folds), function(fold) {
    indices[seq.int(floor((fold - 1) * length(indices) / folds) + 1,
                    floor(fold * length(indices) / folds))]
  })
  scores <- array(NA_real_, c(length(row), length(col), length(balance)))
  for (i in seq_along(row)) for (j in seq_along(col)) for (k in seq_along(balance)) {
    scores[i, j, k] <- sum(vapply(groups, function(held_out) {
      train <- target
      train[held_out] <- NA_real_
      fit <- learner(sources, train, r = 2, lambda_1_row = row[i],
                     lambda_1_col = col[j], lambda_2 = balance[k],
                     source_weights = weights, step_size = 0.01,
                     control = list(max_iter = 9, threshold = 0))
      sum((fit$learner_estimate[held_out] - target[held_out])^2)
    }, numeric(1)))
  }
  scores
}
