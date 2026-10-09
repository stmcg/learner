#' Cross-validation for LEARNER
#'
#' Select common or separate source-space penalties and a balance penalty for
#' \code{\link{learner}} by k-fold cross-validation on observed target entries.
#' One source matrix or a list of source matrices can be supplied.
#'
#' @param lambda_1_all nonempty vector of finite nonnegative candidate common penalties. Use this or both separate grids, not both modes.
#' @param lambda_2_all nonempty vector of finite nonnegative candidate balance penalties.
#' @param lambda_1_row_all nonempty vector of finite nonnegative candidate penalties on U; requires \code{lambda_1_col_all}.
#' @param lambda_1_col_all nonempty vector of finite nonnegative candidate penalties on V; requires \code{lambda_1_row_all}.
#' @inheritParams learner
#' @param step_size numeric scalar specifying the step size for the scaled gradient steps in the numerical optimization algorithm, as in \code{\link{learner}}
#' @param n_folds integer from 2 to the number of observed target entries. The default is \code{4}.
#' @param n_cores positive integer specifying the requested number of OpenMP
#' threads. Candidate penalty combinations are evaluated in parallel when the
#' package is built with OpenMP support; otherwise they run sequentially.
#' The default is \code{1}.
#' @param control a list of parameters for controlling the stopping criteria for the numerical optimization algorithm, as in \code{\link{learner}}.
#'
#' @return A list with the following elements:
#' \item{lambda_1_min}{selected common penalty; \code{NA} in separate mode.}
#' \item{lambda_1_row_min, lambda_1_col_min}{selected separate penalties, returned only in separate mode.}
#' \item{lambda_2_min}{balance penalty in the combination with the smallest held-out squared-error sum.}
#' \item{mse_all}{held-out squared-error sums, retaining the historical scale despite the name. Common mode returns a matrix indexed by \code{lambda_1_all} and \code{lambda_2_all}. Separate mode returns a three-dimensional array indexed by \code{lambda_1_row_all}, \code{lambda_1_col_all}, and \code{lambda_2_all}, in that order.}
#' \item{r}{rank value used.}
#'
#' @details
#' Each observed target entry is held out once. The same folds are used for all
#' candidates. Dividing \code{mse_all} by the number of observed target entries
#' gives the mean squared error without changing the selected candidate. Ties
#' select the first minimum in R array order. Existing missing entries are
#' excluded from validation. Rank, source spaces, source weights, and the
#' first-source initialization are fixed across folds. The function tunes
#' penalties, not rank or source weights. Call \code{set.seed()} before
#' \code{cv.learner()} to reproduce the fold allocation.
#' For independent penalties, supply both separate grids and \code{lambda_2_all}.
#' Candidates with non-finite errors are excluded with a warning; if all fail,
#' the function stops. Finite fits that reach the iteration limit are retained,
#' so selecting a candidate does not certify convergence of its fold fits.
#' Assess the optimization settings as well as the final refit.
#'
#' The return value contains scores and selected parameters, not a final fitted
#' target matrix. Refit with \code{learner()} using all observed target entries
#' and the selected penalties, while retaining the same sources, weights, and rank.
#'
#' Observed target entries are randomly partitioned into approximately equal
#' folds. Each fold is set to missing in the training matrix, and the fitted
#' matrix is used to predict the held-out entries. Errors are summed across
#' folds for each candidate pair (common penalties) or triple (separate
#' penalties). Source entries are not held out. See McGrath et al. (2024).
#'
#' @references
#' McGrath, S., Zhu, C,. Guo, M. and Duan, R. (2024). \emph{LEARNER: A transfer learning method for low-rank matrix estimation}. arXiv preprint	arXiv:2412.20605.
#'
#' @examples
#' res <- cv.learner(Y_source = dat_highsim$Y_source,
#'                   Y_target = dat_highsim$Y_target,
#'                   lambda_1_all = c(1, 10, 100),
#'                   lambda_2_all = c(1, 10, 100),
#'                   step_size = 0.003)
#'
#' # A small simulated example: source B is a perturbation of source A,
#' # not an independent population. Weights are fixed demonstration values.
#' set.seed(803)
#' source_A <- dat_highsim$Y_source[1:20, 1:10]
#' source_B <- source_A + matrix(rnorm(length(source_A), sd = 0.3),
#'                               nrow(source_A), ncol(source_A))
#' sources <- list(A = source_A, B = source_B)
#' target <- dat_highsim$Y_target[1:20, 1:10]
#' weights <- c(0.7, 0.3)
#'
#' set.seed(815)
#' res_separate <- cv.learner(
#'   Y_source = sources, Y_target = target, r = 2,
#'   source_weights = weights,
#'   lambda_1_row_all = c(1, 10), lambda_1_col_all = c(1, 3),
#'   lambda_2_all = c(0.1, 1), step_size = 0.0003,
#'   n_folds = 2, control = list(max_iter = 100))
#' dim(res_separate$mse_all) # row grid x column grid x balance grid
#'
#' fit <- learner(
#'   Y_source = sources, Y_target = target, r = res_separate$r,
#'   source_weights = weights,
#'   lambda_1_row = res_separate$lambda_1_row_min,
#'   lambda_1_col = res_separate$lambda_1_col_min,
#'   lambda_2 = res_separate$lambda_2_min, step_size = 0.0003,
#'   control = list(max_iter = 100))
#' fit$convergence_criterion
#' # This short example demonstrates the workflow; check convergence before
#' # interpreting results or comparing statistical performance.
#'
#' @export
cv.learner <- function(Y_source, Y_target, r, lambda_1_all = NULL, lambda_2_all,
                       step_size, n_folds = 4, n_cores = 1, control = list(),
                       lambda_1_row_all = NULL, lambda_1_col_all = NULL,
                       source_weights = NULL) {
  sources <- prepare_sources(Y_source, Y_target, source_weights)
  r <- resolve_rank(if (missing(r)) NULL else r, sources$matrices[[1]])
  control <- prepare_control(control, step_size)
  validate_penalty_grid(lambda_2_all, 'lambda_2_all')

  # A common penalty uses a 2D grid; independent penalties use a 3D grid.
  separate <- !is.null(lambda_1_row_all) || !is.null(lambda_1_col_all)
  if (separate) {
    if (!is.null(lambda_1_all)) {
      stop('Provide either lambda_1_all or both separate grids, not both modes.')
    }
    if (is.null(lambda_1_row_all) || is.null(lambda_1_col_all)) {
      stop('Both lambda_1_row_all and lambda_1_col_all must be provided together.')
    }
    validate_penalty_grid(lambda_1_row_all, 'lambda_1_row_all')
    validate_penalty_grid(lambda_1_col_all, 'lambda_1_col_all')
    grid <- expand.grid(row = lambda_1_row_all, col = lambda_1_col_all,
                        balance = lambda_2_all)
    grid_dim <- c(length(lambda_1_row_all), length(lambda_1_col_all), length(lambda_2_all))
  } else {
    validate_penalty_grid(lambda_1_all, 'lambda_1_all')
    grid <- expand.grid(row = lambda_1_all, balance = lambda_2_all)
    grid$col <- grid$row
    grid_dim <- c(length(lambda_1_all), length(lambda_2_all))
  }

  # Only observed target entries can be validation entries.
  available_indices <- which(!is.na(Y_target))
  n_indices <- length(available_indices)
  validate_count(n_folds, 'n_folds', lower = 2, upper = n_indices)
  validate_count(n_cores, 'n_cores')
  indices <- sample(available_indices, size = n_indices, replace = FALSE) - 1L
  index_set <- vector('list', n_folds)
  for (fold in seq_len(n_folds)) {
    start <- floor((fold - 1) * n_indices / n_folds) + 1L
    end <- floor(fold * n_indices / n_folds)
    index_set[[fold]] <- indices[seq.int(start, end)]
  }

  scores <- cv_learner_cpp(sources$matrices, Y_target, sources$weights,
                           grid$row, grid$col, grid$balance, step_size,
                           control$max_iter, control$threshold, n_cores, r,
                           control$max_value, index_set)
  if (!any(is.finite(scores))) stop('No candidate produced finite validation errors.')
  if (any(!is.finite(scores))) warning('Some candidates produced non-finite validation errors and were excluded.')
  scores[!is.finite(scores)] <- Inf
  best <- which.min(scores)
  # Keep the historical mse_all scale: summed held-out squared errors.
  mse_all <- array(scores, dim = grid_dim)
  if (!separate) {
    return(list(lambda_1_min = grid$row[best], lambda_2_min = grid$balance[best],
                mse_all = mse_all, r = r))
  }
  list(lambda_1_min = NA_real_, lambda_2_min = grid$balance[best],
       lambda_1_row_min = grid$row[best], lambda_1_col_min = grid$col[best],
       mse_all = mse_all, r = r)
}

#' Latent space-based transfer learning
#'
#' Estimate a low-rank target matrix using information from one or more source
#' populations with the LatEnt spAce-based tRaNsfer lEaRning (LEARNER) method
#' (McGrath et al. 2024). The two source-space penalties can share a coefficient
#' or use separate coefficients.
#'
#' @param Y_target nonempty numeric matrix containing target population data.
#' Missing entries (\code{NA} or \code{NaN}) are allowed, but at least one entry
#' must be observed. Infinite values are not allowed.
#' @param Y_source numeric source matrix or nonempty list of numeric source matrices. Every source must have the same dimensions and corresponding row/column order as \code{Y_target}, and contain no missing or infinite values.
#' @param r optional common rank for the source spaces and target factors, between 1 and the smaller matrix dimension. If omitted, ScreeNOT selects it from the first source.
#' @param source_weights optional finite nonnegative numeric vector with one
#' weight per source and at least one positive value. Weights correspond to the
#' order of matrices in \code{Y_source}; names are not used for matching. They
#' are normalized to sum to one. \code{NULL} gives equal weights. For a single
#' source, a supplied weight must be a positive scalar and is normalized to one.
#' @param lambda_1 finite nonnegative numeric scalar specifying the common space penalty. Supply either this parameter or both \code{lambda_1_row} and \code{lambda_1_col}, not both modes.
#' @param lambda_2 finite nonnegative numeric scalar for the balance penalty (see Details)
#' @param lambda_1_row optional finite nonnegative numeric scalar for the source-space penalty on \eqn{U} (the SNP side when SNPs are rows). Must be supplied together with \code{lambda_1_col}.
#' @param lambda_1_col optional finite nonnegative numeric scalar for the source-space penalty on \eqn{V} (the phenotype side when phenotypes are columns). Must be supplied together with \code{lambda_1_row}.
#' @param step_size numeric scalar specifying the step size for the scaled gradient steps in the numerical optimization algorithm
#' @param control list controlling the numerical optimization:
#' \describe{
#' \item{\code{max_iter}}{Positive integer giving the maximum number of updates;
#' default \code{100}.}
#' \item{\code{threshold}}{Finite nonnegative absolute objective-change threshold;
#' default \code{0.001}. Starting with the second update, the algorithm stops
#' when \eqn{|L_t-L_{t-1}|} is strictly less than this value. This rule alone does
#' not establish factor stability or a global optimum.}
#' \item{\code{max_value}}{Finite positive multiplier; default \code{10}.
#' Starting with the second update, the algorithm stops if the current objective
#' exceeds this multiplier times the preceding objective (unless the
#' objective-change stopping rule has already been met).}}
#'
#' @return A list with the following elements:
#' \item{learner_estimate}{target signal estimate with the same dimensions as \code{Y_target}, formed from the factors with the lowest objective encountered, including the initialization.}
#' \item{objective_values}{objective values after each update, excluding the initialization. The last value need not correspond to the returned estimate.}
#' \item{convergence_criterion}{stopping code: \code{1}, the absolute objective-change threshold was met; \code{2}, the iteration limit was reached; \code{3}, excessive objective growth or non-finite objective/factors were detected.}
#' \item{r}{rank value used.}
#'
#' @details
#'
#' \strong{Data and source spaces:}
#'
#' Let \eqn{Y_0} be the \eqn{p \times q} target matrix, and let
#' \eqn{Y_1,\ldots,Y_K} be aligned source matrices of the same dimensions.
#' The target signal is estimated by \eqn{UV^\top}, with
#' \eqn{U \in \mathbb{R}^{p \times r}} and
#' \eqn{V \in \mathbb{R}^{q \times r}}. For each source, the first \eqn{r}
#' left and right singular vectors form bases \eqn{B_{U,k}} and \eqn{B_{V,k}}.
#' The source operators are
#' \deqn{\bar P_U = \sum_{k=1}^K w_k B_{U,k}B_{U,k}^\top, \qquad
#'       \bar P_V = \sum_{k=1}^K w_k B_{V,k}B_{V,k}^\top,}
#' where \eqn{w_k \geq 0} and \eqn{\sum_k w_k=1}.
#' These combine source spaces, not the raw source matrices. The averaged
#' operators are generally not orthogonal projectors themselves.
#' Rows and columns must already be aligned by the caller; names are not matched.
#'
#' \strong{Objective and penalties:}
#'
#' The numerical algorithm seeks to minimize
#' \deqn{L(U,V) = \rho^{-1}\|\mathcal{P}_{\Omega}(UV^\top-Y_0)\|_F^2
#'   + \lambda_{1,\mathrm{row}}\|(I-\bar P_U)U\|_F^2
#'   + \lambda_{1,\mathrm{col}}\|(I-\bar P_V)V\|_F^2
#'   + \lambda_2\|U^\top U-V^\top V\|_F^2.}
#' Here \eqn{\Omega} contains the observed target entries,
#' \eqn{\mathcal{P}_{\Omega}} keeps residuals at those entries and sets the
#' others to zero, and \eqn{\rho=|\Omega|/(pq)} is the observed fraction.
#' With no missing entries, \eqn{\rho=1}. Source matrices must be complete.
#'
#' Supply \code{lambda_1} for a common space penalty, or supply both
#' \code{lambda_1_row} and \code{lambda_1_col} and leave \code{lambda_1=NULL}.
#' Mixing the modes or providing only one separate coefficient is an error.
#' Equal separate coefficients reproduce the common-penalty calculation under
#' the same inputs and optimization settings. The balance penalty
#' \code{lambda_2} is required in either mode.
#' The space penalties square the residual from the weighted source operator;
#' they are not weighted sums of separately squared per-source residuals.
#' Use \code{\link{cv.learner}} to select common or separate penalty values.
#'
#' \strong{Initialization and stopping:}
#'
#' If \code{r} is omitted, ScreeNOT selects it from the first source, with a
#' minimum of one. A supplied \code{r} is used directly. The first source is
#' always used to initialize the target factors, including when its weight is
#' zero. Reordering sources with their weights therefore preserves the weighted
#' operators for fixed \code{r}, but can change initialization and fitted results.
#' The initial factors are the first source's retained left and right singular
#' vectors, each multiplied by the square roots of its retained singular values.
#'
#' This function uses scaled gradient descent steps. Both gradients are
#' evaluated at the current \eqn{U} and \eqn{V}, then both matrices are updated.
#' Updates are repeated until a stopping criterion is satisfied. The returned
#' estimate uses the best objective encountered, not necessarily the final
#' iterate. Inspect \code{objective_values} and \code{convergence_criterion};
#' code 1 only indicates that the stated objective-change rule was met.
#'
#' @references
#' McGrath, S., Zhu, C,. Guo, M. and Duan, R. (2024). \emph{LEARNER: A transfer learning method for low-rank matrix estimation}. arXiv preprint	arXiv:2412.20605.
#'
#' Donoho, D., Gavish, M. and Romanov, E. (2023). \emph{ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise}. The Annals of Statistics, 51(1), pp.122-148.
#'
#' @examples
#' res <- learner(Y_source = dat_highsim$Y_source,
#'                Y_target = dat_highsim$Y_target,
#'                lambda_1 = 1, lambda_2 = 1,
#'                step_size = 0.003)
#'
#' res_separate <- learner(Y_source = dat_highsim$Y_source,
#'                         Y_target = dat_highsim$Y_target,
#'                         lambda_1_row = 10, lambda_1_col = 1,
#'                         lambda_2 = 1, step_size = 0.003)
#'
#' # Multiple sources with fixed demonstration weights.
#' # Source B is a perturbed copy, not an independent population.
#' set.seed(803)
#' source_A <- dat_highsim$Y_source[1:20, 1:10]
#' source_B <- source_A + matrix(rnorm(length(source_A), sd = 0.3),
#'                               nrow(source_A), ncol(source_A))
#' res_multi <- learner(
#'   Y_source = list(A = source_A, B = source_B),
#'   Y_target = dat_highsim$Y_target[1:20, 1:10], r = 2,
#'   source_weights = c(0.7, 0.3),
#'   lambda_1_row = 10, lambda_1_col = 1, lambda_2 = 1,
#'   step_size = 0.0003, control = list(max_iter = 100))
#' dim(res_multi$learner_estimate)
#' res_multi$convergence_criterion
#' # Short runs demonstrate usage; they need not meet the stopping threshold.
#'
#' @export
learner <- function(Y_source, Y_target, r, lambda_1 = NULL, lambda_2, step_size,
                    control = list(), lambda_1_row = NULL, lambda_1_col = NULL,
                    source_weights = NULL) {
  sources <- prepare_sources(Y_source, Y_target, source_weights)
  r <- resolve_rank(if (missing(r)) NULL else r, sources$matrices[[1]])
  control <- prepare_control(control, step_size)
  # Resolve the original common penalty or the two independent penalties.
  separate <- !is.null(lambda_1_row) || !is.null(lambda_1_col)
  if (separate){
    if (!is.null(lambda_1)){
      stop('Provide either lambda_1 or both lambda_1_row and lambda_1_col, not both modes.')
    }
    if (is.null(lambda_1_row) || is.null(lambda_1_col)){
      stop('Both lambda_1_row and lambda_1_col must be provided together.')
    }
    validate_space_penalty(lambda_1_row, 'lambda_1_row')
    validate_space_penalty(lambda_1_col, 'lambda_1_col')
  } else {
    if (is.null(lambda_1)){
      stop('Provide lambda_1 or both lambda_1_row and lambda_1_col.')
    }
    validate_space_penalty(lambda_1, 'lambda_1')
    lambda_1_row <- lambda_1
    lambda_1_col <- lambda_1
  }
  validate_space_penalty(lambda_2, 'lambda_2')

  result <- learner_cpp(sources$matrices, Y_target, sources$weights, r,
                        lambda_1_row, lambda_1_col, lambda_2, step_size,
                        control$max_iter, control$threshold, control$max_value)
  return(result)
}

validate_space_penalty <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value < 0){
    stop(name, ' must be a finite nonnegative numeric scalar.', call. = FALSE)
  }
}
