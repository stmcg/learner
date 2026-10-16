# Normalize one source or a source list without changing the user's ordering.
# The first source determines the default rank and LEARNER initialization.
prepare_sources <- function(Y_source, Y_target, source_weights = NULL,
                            allow_missing_target = TRUE) {
  if (!is.matrix(Y_target) || !is.numeric(Y_target) ||
      any(dim(Y_target) == 0L) || any(is.infinite(Y_target))) {
    stop('Y_target must be a nonempty numeric matrix without infinite values.')
  }
  if ((!allow_missing_target && anyNA(Y_target)) || all(is.na(Y_target))) {
    stop('Y_target must contain observed values; dlearner does not allow NA values.')
  }
  sources <- if (is.matrix(Y_source)) list(Y_source) else Y_source
  if (!is.list(sources) || length(sources) == 0L) {
    stop('Y_source must be a matrix or a nonempty list of matrices.')
  }
  for (source in sources) {
    if (!is.matrix(source) || !is.numeric(source) || any(!is.finite(source))) {
      stop('Every source must be a numeric matrix with finite values and no NA values.')
    }
    if (!identical(dim(source), dim(Y_target))) {
      stop('All source matrices and Y_target must have the same dimensions.')
    }
  }
  if (is.null(source_weights)) source_weights <- rep(1, length(sources))
  if (!is.numeric(source_weights) || length(source_weights) != length(sources) ||
      any(!is.finite(source_weights)) || any(source_weights < 0) ||
      !any(source_weights > 0)) {
    stop('source_weights must have one finite nonnegative value per source and a positive total.')
  }
  weights <- source_weights / max(source_weights)
  list(matrices = sources, weights = weights / sum(weights))
}

resolve_rank <- function(r, source) {
  if (is.null(r)) {
    r <- max(ScreeNOT::adaptiveHardThresholding(
      Y = source, k = min(dim(source)) / 3)$r, 1)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) ||
      r < 1 || r != floor(r) || r > min(dim(source))) {
    stop('r must be an integer between 1 and the smaller matrix dimension.')
  }
  as.integer(r)
}

validate_penalty_grid <- function(values, name) {
  if (!is.numeric(values) || length(values) == 0L ||
      any(!is.finite(values)) || any(values < 0)) {
    stop(name, ' must be a nonempty vector of finite nonnegative numbers.', call. = FALSE)
  }
}

validate_count <- function(value, name, lower = 1, upper = .Machine$integer.max) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper) {
    stop(name, ' must be an integer between ', lower, ' and ', upper, '.', call. = FALSE)
  }
}

prepare_control <- function(control, step_size) {
  if (!is.list(control)) stop('control must be a list.')
  if (!is.numeric(step_size) || length(step_size) != 1L ||
      !is.finite(step_size) || step_size <= 0) {
    stop('step_size must be a finite positive number.')
  }
  if (is.null(control$max_iter)) control$max_iter <- 100L
  if (is.null(control$threshold)) control$threshold <- 0.001
  if (is.null(control$max_value)) control$max_value <- 10
  validate_count(control$max_iter, 'control$max_iter')
  validate_space_penalty(control$threshold, 'control$threshold')
  if (!is.numeric(control$max_value) || length(control$max_value) != 1L ||
      !is.finite(control$max_value) || control$max_value <= 0) {
    stop('control$max_value must be a finite positive number.')
  }
  control
}
