# Normalize one source or a source list without changing the user's ordering.
# The first source determines the default rank and LEARNER initialization.
prepare_sources <- function(Y_source, Y_target, source_weights = NULL,
                            allow_missing_target = TRUE) {
  if (!is.matrix(Y_target) || !is.numeric(Y_target) ||
      any(dim(Y_target) == 0L) || any(is.infinite(Y_target))) {
    stop('Y_target must be a nonempty numeric matrix without infinite values.')
  }
  if (all(is.na(Y_target))) {
    stop('Y_target must contain observed values.')
  }
  if (!allow_missing_target && anyNA(Y_target)) {
    stop('dlearner requires a complete Y_target without missing values.')
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
  mode <- "manual"
  if (is.character(source_weights)) {
    if (length(source_weights) != 1L || is.na(source_weights) ||
        !source_weights %in% c("auto", "equal")) {
      stop("source_weights must be 'auto', 'equal', NULL, or numeric relative weights.")
    }
    mode <- source_weights
    source_weights <- rep(1, length(sources))
  } else if (is.null(source_weights)) {
    # NULL uses the empirical subspace-distance rule.
    mode <- "auto"
    source_weights <- rep(1, length(sources))
  }
  if (!is.numeric(source_weights) || length(source_weights) != length(sources) ||
      any(!is.finite(source_weights)) || any(source_weights < 0) ||
      !any(source_weights > 0)) {
    stop('source_weights must have one finite nonnegative value per source and a positive total.')
  }
  weights <- source_weights / max(source_weights)
  list(matrices = sources, weights = weights / sum(weights), mode = mode)
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

# k and epsilon are accepted as other names for window and threshold. Only one
# name may be given for each setting. The value is stored under the main name, and
# the label tells the validators which name the user actually wrote.
use_control_alias <- function(control, name, alias) {
  label <- paste0('control$', name)
  alias_value <- control[[alias]]
  if (is.null(alias_value)) return(list(control = control, label = label))
  if (!is.null(control[[name, exact = FALSE]])) {
    stop('Give control$', name, ' or control$', alias, ', not both.', call. = FALSE)
  }
  control[[name]] <- alias_value
  control[[alias]] <- NULL
  list(control = control, label = paste0('control$', alias))
}

prepare_control <- function(control, step_size) {
  if (!is.list(control)) stop('control must be a list.')
  if (!is.numeric(step_size) || length(step_size) != 1L ||
      !is.finite(step_size) || step_size <= 0) {
    stop('step_size must be a finite positive number.')
  }
  window_alias <- use_control_alias(control, 'window', 'k')
  control <- window_alias$control
  threshold_alias <- use_control_alias(control, 'threshold', 'epsilon')
  control <- threshold_alias$control
  if (is.null(control$max_iter)) control$max_iter <- 100L
  if (is.null(control$stopping_rule)) control$stopping_rule <- 'absolute'
  if (!is.character(control$stopping_rule) || length(control$stopping_rule) != 1L ||
      is.na(control$stopping_rule) || !control$stopping_rule %in% c('absolute', 'relative')) {
    stop("control$stopping_rule must be 'absolute' or 'relative'.")
  }
  if (is.null(control$window)) {
    control$window <- if (control$stopping_rule == 'relative') 5L else 1L
  }
  validate_count(control$window, window_alias$label)
  if (is.null(control$threshold)) {
    control$threshold <- if (control$stopping_rule == 'relative') 1e-6 else 0.001
  }
  if (is.null(control$max_value)) control$max_value <- 10
  validate_count(control$max_iter, 'control$max_iter')
  validate_space_penalty(control$threshold, threshold_alias$label)
  if (!is.numeric(control$max_value) || length(control$max_value) != 1L ||
      !is.finite(control$max_value) || control$max_value <= 0) {
    stop('control$max_value must be a finite positive number.')
  }
  control
}

# Equation (3): inverse summed squared distances of row/column projectors.
# Use thin bases to avoid allocating p-by-p or q-by-q projection matrices.
source_subspace_weights <- function(bases, target, r) {
  if (length(bases) == 1L) return(1)
  # A spectral pilot for incomplete targets: only observed values enter the SVD.
  # The inverse observation-fraction scaling would not change singular vectors.
  target[is.na(target)] <- 0
  pilot <- svd(target, nu = r, nv = r)
  distances <- vapply(bases, function(b) {
    du <- pilot$u - b$u %*% crossprod(b$u, pilot$u)
    dv <- pilot$v - b$v %*% crossprod(b$v, pilot$v)
    2 * (sum(du^2) + sum(dv^2))
  }, numeric(1))
  # Same 1e-12 stabilizer as learner-paper's multisource simulation helper.
  inverse <- 1 / (distances + 1e-12)
  unname(inverse / sum(inverse))
}

resolve_source_weights <- function(sources, target, r, bases = NULL) {
  if (sources$mode != "auto" || length(sources$matrices) == 1L) {
    return(sources$weights)
  }
  if (is.null(bases)) bases <- lapply(sources$matrices, svd, nu = r, nv = r)
  source_subspace_weights(bases, target, r)
}
