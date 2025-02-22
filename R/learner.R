#' Cross-validation for LEARNER
#'
#' This function performs k-fold cross-validation to select the nuisance parameters \eqn{(\lambda_1, \lambda_2)} for \code{\link{learner}}.
#'
#' @param Y_target matrix containing the target population data, as in \code{\link{learner}}
#' @param Y_source matrix containing the source population data, as in \code{\link{learner}}
#' @param r (optional) integer specifying the rank of the knowledge graphs, as in \code{\link{learner}}
#' @param lambda_1_all vector of numerics specifying the candidate values of \eqn{\lambda_1} (see Details)
#' @param lambda_2_all vector of numerics specifying the candidate values of \eqn{\lambda_2} (see Details)
#' @param step_size numeric scalar specifying the step size for the Newton steps in the numerical optimization algorithm, as in \code{\link{learner}}
#' @param n_folds an integer specify the number of cross-validation folds. The default is \code{4}.
#' @param n_cores an integer specifying the number of CPU cores in OpenMP parallelization. Parallelization is performed across the different candidate \eqn{(\lambda_1, \lambda_2)} pairs. The default is \code{1}, i.e., no parallelization.
#' @param threshold Convergence threshold.
#' @param max_iter Maximum number of iterations.
#'
#' @return A list with the following elements:
#' \item{lambda_1_min}{value of \eqn{\lambda_1} with the smallest MSE}
#' \item{lambda_2_min}{value of \eqn{\lambda_2} with the smallest MSE}
#' \item{mse_all}{matrix containing MSE value for each \eqn{(\lambda_1, \lambda_2)} pair. The rows correspond to the \eqn{\lambda_1} values, and the columns correspond to the \eqn{\lambda_2} values.}
#' \item{r}{rank value used.}
#'
#' @details
#' Given sets of candidate values of \eqn{\lambda_1} and \eqn{\lambda_2}, this function performs k-fold cross-validation to select the pair \eqn{(\lambda_1, \lambda_2)} with the smallest held out error. This function randomly partitions the entries of \code{Y_target} into \eqn{k} (approximately) equally sized subsamples. The training data sets are obtained by removing one of the \eqn{k} subsamples and the corresponding test data sets are based on the held out subsamples. The \code{\link{learner}} function is applied to each training data set. The held out error is computed by the mean squared error comparing the entries in the test data sets with those imputed from the LEARNER estimates. See McGrath et al. (2024) for further details.
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
#'
#'
#' @export
cv.learner <- function(Y_source, Y_target, lambda_1_all, lambda_2_all, step_size,
                       n_folds = 4, max_iter = 100, threshold = 1e-3, n_cores = 1, r) {
  ## Error catching
  if (!identical(dim(Y_source), dim(Y_target))){
    stop('Y_source and Y_target must have the same dimensions')
  }
  if (any(is.na(Y_source))){
    stop('Y_source cannot have NA values.')
  }

  if (missing(r)){
    p <- nrow(Y_source); q <- ncol(Y_source)
    max_rank <- min(p, q) / 3
    r <- max(ScreeNOT::adaptiveHardThresholding(Y = Y_source, k = max_rank)$r, 1)
  }

  result <- cv_learner_cpp(Y_source, Y_target, lambda_1_all, lambda_2_all, step_size,
                           n_folds, max_iter, threshold, n_cores, r)
  return(result)
}

#' Learner Optimization
#'
#' @param Y_source Matrix containing the source population data.
#' @param Y_target Matrix containing the target population data.
#' @param lambda_1 Regularization parameter lambda1.
#' @param lambda_2 Regularization parameter lambda2.
#' @param step_size Step size for optimization.
#' @param max_iter Maximum number of iterations.
#' @param threshold Convergence threshold.
#' @param r Rank used for the truncated SVD.
#' @return A list containing the learner estimate, objective values, convergence criterion, and rank.
#' @export
learner <- function(Y_source, Y_target, lambda_1, lambda_2, step_size,
                    max_iter = 100, threshold = 1e-3, r) {
  ## Error catching
  if (!identical(dim(Y_source), dim(Y_target))){
    stop('Y_source and Y_target must have the same dimensions')
  }
  if (any(is.na(Y_source))){
    stop('Y_source cannot have NA values.')
  }
  if (length(lambda_1) > 1){
    stop('lambda_1 must be a scalar. See the function cv.learner for selecting a suitable lambda_1 value.')
  }
  if (length(lambda_2) > 1){
    stop('lambda_2 must be a scalar. See the function cv.learner for selecting a suitable lambda_2 value.')
  }

  if (missing(r)){
    p <- nrow(Y_source); q <- ncol(Y_source)
    max_rank <- min(p, q) / 3
    r <- max(ScreeNOT::adaptiveHardThresholding(Y = Y_source, k = max_rank)$r, 1)
  }

  result <- learner_cpp(Y_source, Y_target, r, lambda_1, lambda_2, step_size, max_iter, threshold)
  return(result)
}


