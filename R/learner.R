#' Cross-validation for LEARNER
#'
#' @description
#' Performs k-fold cross-validation to select the source-space and balance
#' penalties for \code{\link{learner}}. Returns prediction errors and the
#' penalty combination with the smallest error.
#'
#' @param Y_source Source population data, as in \code{\link{learner}}.
#' Supply a matrix for one source or a list of matrices for multiple sources.
#' @param Y_target Target population data matrix, as in \code{\link{learner}}.
#' Only observed entries are used to evaluate predictions.
#' @param r Rank used in every fold, as in \code{\link{learner}}. It is not
#' selected by cross-validation.
#' @param lambda_1_all Candidate values for the common source-space penalty
#' \eqn{\lambda_1}. All combinations with \code{lambda_2_all} are evaluated.
#' @param lambda_2_all Candidate values for the balance penalty \eqn{\lambda_2}.
#' @param step_size Step size for each fit, as in \code{\link{learner}}.
#' @param n_folds Number of cross-validation folds. Default is \code{4}.
#' @param n_cores Number of CPU cores used to evaluate the candidate penalty
#' combinations in parallel. Default is \code{1}, i.e., no parallelization.
#' @param control Optimization settings for each fit; see \code{\link{learner}}.
#' @param lambda_1_row_all,lambda_1_col_all Candidate values for separate
#' penalties on the row and column factors. Supply both instead of
#' \code{lambda_1_all}; all combinations with \code{lambda_2_all} are evaluated.
#' @param source_weights Relative source weights, as in \code{\link{learner}}. The weights
#' are not selected by cross-validation. Default is equal weights.
#'
#' @details
#' Given sets of candidate values of \eqn{\lambda_1} and \eqn{\lambda_2}, this
#' function performs k-fold cross-validation to select the pair
#' \eqn{(\lambda_1, \lambda_2)} with the smallest held out error. This function
#' randomly partitions the entries of \code{Y_target} into \eqn{k}
#' (approximately) equally sized subsamples. The training data sets are obtained
#' by removing one of the \eqn{k} subsamples and the corresponding test data
#' sets are based on the held out subsamples. The \code{\link{learner}} function is
#' applied to each training data set. The held out error is computed by comparing
#' the entries in the test data sets with those imputed from the LEARNER
#' estimates. See McGrath et al. (2024) for further details.
#'
#' With separate penalties, the procedure selects the row, column, and balance
#' penalties together. Both modes can be used with multiple sources. Call
#' \code{set.seed()} first to reproduce the random partition.
#'
#' To obtain the final estimate, apply \code{\link{learner}} with the selected
#' penalties to the full data (see Examples). Fits that reach the iteration
#' limit still contribute to the held out error, so check the settings in
#' \code{control} before interpreting the selection.
#'
#' @return
#' A list containing:
#' \item{lambda_1_min}{Common source-space penalty in the selected combination.
#' \code{NA} when separate penalties are used.}
#' \item{lambda_2_min}{Balance penalty in the selected combination.}
#' \item{lambda_1_row_min, lambda_1_col_min}{Row and column penalties in the
#' selected combination. Returned only when separate penalties are used.}
#' \item{mse_all}{Sum of held-out squared errors for each candidate combination.
#' For common penalties, rows correspond to \code{lambda_1_all} and columns to
#' \code{lambda_2_all}. For separate penalties, the three dimensions correspond
#' to the row, column, and balance grids, in that order. Divide by the number
#' of observed target entries to obtain mean squared errors.}
#' \item{r}{Rank used in the fits.}
#'
#' @references
#' McGrath, S., Zhu, C,. Guo, M. and Duan, R. (2024). \emph{LEARNER: A transfer learning method for low-rank matrix estimation}. arXiv preprint	arXiv:2412.20605.
#'
#' @seealso
#' \code{\link{learner}} for fitting with the selected penalties.
#'
#' @examples
#' res <- cv.learner(Y_source = dat_highsim$Y_source,
#'                   Y_target = dat_highsim$Y_target,
#'                   lambda_1_all = c(1, 10, 100),
#'                   lambda_2_all = c(1, 10, 100),
#'                   step_size = 0.003)
#'
#' # Refit using the selected penalties and all observed target entries.
#' fit <- learner(dat_highsim$Y_source, dat_highsim$Y_target, r = res$r,
#'                lambda_1 = res$lambda_1_min, lambda_2 = res$lambda_2_min,
#'                step_size = 0.003)
#' fit$convergence_criterion
#'
#' # Multiple sources with fixed demonstration weights.
#' sources <- dat_multisource$Y_source
#' target <- dat_multisource$Y_target
#' weights <- c(0.7, 0.3)
#'
#' set.seed(815)
#' res_separate <- cv.learner(
#'   Y_source = sources, Y_target = target, r = 3,
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
#' @description
#' This function applies the LatEnt spAce-based tRaNsfer lEaRning (LEARNER)
#' method (McGrath et al. 2024) to leverage data from a source population to improve
#' estimation of a low-rank matrix in an underrepresented target population.
#'
#' @param Y_source Matrix containing the source population data, with the same rows
#' and columns as \code{Y_target}. For multiple sources, supply a list of
#' matrices.
#' @param Y_target Matrix containing the target population data. Missing entries
#' are allowed.
#' @param r Rank of the knowledge graphs. If omitted, ScreeNOT (Donoho et al. 2023)
#' is applied to the source population data to select the rank.
#' @param lambda_1 Common source-space penalty. Larger values borrow more
#' information from the source population; zero does not use the source spaces.
#' Use \code{\link{cv.learner}} to select penalty values.
#' @param lambda_2 Balance penalty. Larger values more strongly encourage
#' \eqn{U} and \eqn{V} to have comparable scale.
#' @param step_size Step size in the optimization algorithm. Smaller values may
#' reduce oscillation but require more iterations.
#' @param control List of settings controlling when optimization stops:
#' \describe{
#' \item{\code{max_iter}}{Maximum number of updates; default \code{100}.}
#' \item{\code{threshold}}{Stop when the absolute change between consecutive
#' objective values is smaller than this value; default \code{0.001}.}
#' \item{\code{max_value}}{Stop if the objective grows beyond this multiple of
#' its preceding value; default \code{10}.}}
#' @param lambda_1_row,lambda_1_col Separate source-space penalties for the row
#' factor \eqn{U} and column factor \eqn{V}, respectively. Use these when the two
#' sides should borrow different amounts of source information. Supply both
#' instead of \code{lambda_1}; equal values reproduce the common penalty.
#' @param source_weights Relative weights of the sources, in the order of
#' \code{Y_source}. Default is equal weights.
#'
#' @details
#' \strong{Data and notation:}
#'
#' The data consists of a matrix in the target population
#' \eqn{Y_0\in\mathbb{R}^{p\times q}} and the source population
#' \eqn{Y_1\in\mathbb{R}^{p\times q}}. Let
#' \eqn{\hat U_k\hat\Lambda_k\hat V_k^\top} denote the rank-\eqn{r} truncated
#' singular value decomposition (SVD) of \eqn{Y_k}, \eqn{k = 0, 1}.
#'
#' For \eqn{k = 0, 1}, one can view \eqn{Y_k} as a noisy version of
#' \eqn{\Theta_k}, referred to as the knowledge graph. The target of inference
#' is the target population knowledge graph, \eqn{\Theta_0}.
#'
#' \strong{Estimation:}
#'
#' This method estimates \eqn{\Theta_0} by \eqn{\tilde U\tilde V^\top}, where
#' \eqn{(\tilde U, \tilde V)} is the solution to the following optimization problem
#'
#' \deqn{\mathrm{arg\,min}_{U\in\mathbb{R}^{p\times r},\,V\in\mathbb{R}^{q\times r}}
#'   \big\{\|UV^\top-Y_0\|_F^2
#'   +\lambda_1\|\mathcal{P}_\perp(\hat U_1)U\|_F^2
#'   +\lambda_1\|\mathcal{P}_\perp(\hat V_1)V\|_F^2
#'   +\lambda_2\|U^\top U-V^\top V\|_F^2\big\}}
#'
#' where \eqn{\mathcal{P}_\perp(\hat U_1) = I - \hat U_1\hat U_1^\top} and
#' \eqn{\mathcal{P}_\perp(\hat V_1) = I - \hat V_1\hat V_1^\top}. The four terms measure
#' target fit, departures from the two source spaces, and imbalance between the
#' target factors, respectively. If \code{Y_target} has missing entries, the
#' first term is computed over the observed entries and divided by the observed
#' fraction (the number of observed entries divided by \eqn{pq}). The penalty
#' terms are unchanged. The same adjustment applies to the extensions below.
#'
#' At each iteration, \eqn{U} and \eqn{V} are updated by a gradient
#' descent step, starting from the truncated SVD of the source population data.
#' These updates are repeated until convergence.
#'
#' \strong{Extensions:}
#'
#' Separate penalties \code{lambda_1_row} and \code{lambda_1_col} replace
#' \eqn{\lambda_1} in the second and third terms. For multiple sources, the
#' source projectors are replaced by weighted combinations:
#'
#' \deqn{\bar P_U=\sum_{k=1}^K w_k\hat U_k\hat U_k^\top,\qquad
#'       \bar P_V=\sum_{k=1}^K w_k\hat V_k\hat V_k^\top.}
#'
#' where \eqn{\hat U_k} and \eqn{\hat V_k} come from the truncated SVD of
#' the \eqn{k}th source and the weights \eqn{w_k} are given by
#' \code{source_weights}. The full objective is
#'
#' \deqn{\begin{aligned}
#' L(U,V) ={}& \|UV^\top-Y_0\|_F^2 \\
#' &+\lambda_{1,\mathrm{row}}\|(I_p-\bar P_U)U\|_F^2 \\
#' &+\lambda_{1,\mathrm{col}}\|(I_q-\bar P_V)V\|_F^2 \\
#' &+\lambda_2\|U^\top U-V^\top V\|_F^2.
#' \end{aligned}}
#'
#' where \eqn{I_p} and \eqn{I_q} are identity matrices. Equal space penalties
#' give the common-penalty version. With multiple sources, the first source is
#' used to initialize the algorithm and, if \code{r} is not specified, to
#' select the rank.
#'
#' @return
#' A list with the following elements:
#' \item{learner_estimate}{matrix containing the LEARNER estimate of the target population knowledge graph.}
#' \item{objective_values}{numeric vector containing the values of the objective function at each iteration.}
#' \item{convergence_criterion}{stopping code: \code{1}, the absolute objective-change threshold was met; \code{2}, the iteration limit was reached; \code{3}, excessive objective growth or non-finite objective/factors were detected.}
#' \item{r}{rank value used.}
#'
#' @references
#' McGrath, S., Zhu, C,. Guo, M. and Duan, R. (2024). \emph{LEARNER: A transfer learning method for low-rank matrix estimation}. arXiv preprint	arXiv:2412.20605.
#'
#' Donoho, D., Gavish, M. and Romanov, E. (2023). \emph{ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise}. The Annals of Statistics, 51(1), pp.122-148.
#'
#' @seealso
#' \code{\link{cv.learner}} for penalty selection;
#' \code{\link{dlearner}} for direct projection onto source spaces.
#'
#' @examples
#' # Single source with a common source-space penalty.
#' res <- learner(Y_source = dat_highsim$Y_source,
#'                Y_target = dat_highsim$Y_target,
#'                lambda_1 = 1, lambda_2 = 1,
#'                step_size = 0.003)
#'
#' # Same source, with separate penalties for the two factors.
#' res_separate <- learner(Y_source = dat_highsim$Y_source,
#'                         Y_target = dat_highsim$Y_target,
#'                         lambda_1_row = 10, lambda_1_col = 1,
#'                         lambda_2 = 1, step_size = 0.003)
#'
#' # Multiple sources with fixed demonstration weights.
#' res_multi <- learner(
#'   Y_source = dat_multisource$Y_source,
#'   Y_target = dat_multisource$Y_target, r = 3,
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
