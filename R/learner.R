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
#' @param control a list of parameters for controlling the stopping criteria for the numerical optimization algorithm, as in \code{\link{learner}}.
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
cv.learner <- function(Y_source, Y_target, r, lambda_1_all, lambda_2_all,
                       step_size, n_folds = 4, n_cores = 1, control = list()) {
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
  if (!('max_iter' %in% names(control))){
    control$max_iter <- 100
  }
  if (!('threshold' %in% names(control))){
    control$threshold <- 0.001
  }
  if (!('max_value' %in% names(control))){
    control$max_value <- 10
  }

  if (n_folds < 2){
    stop('n_folds must be 2 or greater.')
  }

  ## Creating training and testing data sets
  available_indices <- which(!is.na(Y_target))
  n_indices <- length(available_indices)
  indices <- sample(available_indices, size = n_indices, replace = FALSE) - 1

  index_set <- vector(mode = "list", length = n_folds)
  for (fold in 1:n_folds){
    index_set[[fold]] <- indices[floor((fold - 1) * n_indices / n_folds + 1):floor(fold * n_indices / n_folds)]
  }

  result <- cv_learner_cpp(Y_source, Y_target, lambda_1_all, lambda_2_all, step_size,
                           n_folds, control$max_iter, control$threshold, n_cores, r, control$max_value,
                           index_set)
  return(result)
}

#' Latent space-based transfer learning
#'
#' This function applies the LatEnt spAce-based tRaNsfer lEaRning (LEARNER) method (McGrath et al. 2024) to leverage data from a source population to improve
#' estimation of a low rank matrix in an underrepresented target population.
#'
#' @param Y_target matrix containing the target population data
#' @param Y_source matrix containing the source population data
#' @param r (optional) integer specifying the rank of the knowledge graphs. By default, ScreeNOT (Donoho et al. 2023) is applied to the source population knowledge graph to select the rank.
#' @param lambda_1 numeric scalar specifying the value of \eqn{\lambda_1} (see Details)
#' @param lambda_2 numeric scalar specifying the value of \eqn{\lambda_2} (see Details)
#' @param step_size numeric scalar specifying the step size for the Newton steps in the numerical optimization algorithm
#' @param control a list of parameters for controlling the stopping criteria for the numerical optimization algorithm. The list may include the following components:
#' \tabular{ll}{
#' \code{max_iter} \tab integer specifying the maximum number of iterations \cr
#' \code{threshold} \tab numeric scalar specifying a convergence threshold. The algorithm converges when \eqn{|\epsilon_t - \epsilon_{t-1}| < }\code{threshold}, where \eqn{\epsilon_t} denotes the value of the objective function at iteration \eqn{t}. \cr
#' \code{max_value} \tab numeric scalar used to specify the maximum value of the objective function allowed before terminating the algorithm. Specifically, the algorithm will terminate if the value of the objective function exceeds \code{max_value}\eqn{\times \epsilon_0}, where \eqn{\epsilon_0} denotes the value of the objective function at the initial point. This is used to prevent unnecessary computation time after the optimization algorithm diverges. \cr}
#'
#' @return A list with the following elements:
#' \item{learner_estimate}{matrix containing the LEARNER estimate of the target population knowledge graph}
#' \item{objective_values}{numeric vector containing the values of the objective function at each iteration}
#' \item{convergence_criterion}{integer specifying the criterion that was satisfied for terminating the numerical optimization algorithm. A value of 1 indicates the convergence threshold was satisfied; A value of 2 indicates that the maximum number of iterations was satisfied; A value of 3 indicates that the maximum value of the objective function was satisfied.}
#' \item{r}{rank value used.}
#'
#' @details
#'
#' \strong{Data and notation:}
#'
#' The data consists of a matrix in the target population \eqn{Y_0 \in \mathbb{R}^{p \times q}} and the source population \eqn{Y_1 \in \mathbb{R}^{p \times q}}.
#' Let \eqn{\hat{U}_{k} \hat{\Lambda}_{k} \hat{V}_{k}^{\top}} denote the truncated singular value decomposition (SVD) of \eqn{Y_k}, \eqn{k = 0, 1}.
#'
#' For \eqn{k = 0, 1}, one can view \eqn{Y_k} as a noisy version of \eqn{\Theta_k}, referred to as the knowledge graph. The target of inference is the target population knowledge graph, \eqn{\Theta_0}.
#'
#' \strong{Estimation:}
#'
#' This method estimates \eqn{\Theta_0} by \eqn{\tilde{U}\tilde{V}^{\top}}, where \eqn{(\tilde{U}, \tilde{V})} is the solution to the following optimization problem
#' \deqn{\mathrm{arg\,min}_{U \in \mathbb{R}^{p \times r}, V \in \mathbb{R}^{q \times r}} \big\{ \| U V^{\top} - Y_0 \|_F^2 + \lambda_1\| \mathcal{P}_{\perp}(\hat{U}_{1})U \|_F^2 + \lambda_1\|  \mathcal{P}_{\perp}(\hat{V}_{1})V \|_F^2  + \lambda_2 \| U^{\top} U - V^{\top} V \|_F^2 \big\}}
#' where \eqn{\mathcal{P}_{\perp}(\hat{U}_{1}) = I - \hat{U}_{1}^{\top}\hat{U}_{1}} and \eqn{\mathcal{P}_{\perp}(\hat{V}_{1}) = I - \hat{V}_{1}^{\top}\hat{V}_{1}}.
#'
#' This function uses an alternating minimization strategy to solve the optimization problem. That is, this approach updates \eqn{U} by minimizing the objective function (via a gradient descent step) treating \eqn{V} as fixed. Then, \eqn{V} is updated treating \eqn{U} as fixed. These updates of \eqn{U} and \eqn{V} are repeated until convergence.
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
#' @export
learner <- function(Y_source, Y_target, r, lambda_1, lambda_2, step_size,
                    control = list()) {
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
  if (!('max_iter' %in% names(control))){
    control$max_iter <- 100
  }
  if (!('threshold' %in% names(control))){
    control$threshold <- 0.001
  }
  if (!('max_value' %in% names(control))){
    control$max_value <- 10
  }

  result <- learner_cpp(Y_source, Y_target, r, lambda_1, lambda_2, step_size, control$max_iter, control$threshold, control$max_value)
  return(result)
}


