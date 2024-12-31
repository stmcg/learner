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
#' @param n_cores an integer specifying the number of CPU cores in parallelization. Parallelization is performed across the different candidate \eqn{(\lambda_1, \lambda_2)} pairs. The default is \code{1}, i.e., no parallelization.
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
#' @import foreach
#' @import doParallel
#'
#' @export



cv.learner <- function(Y_source, Y_target, r, lambda_1_all, lambda_2_all,
                       step_size, n_folds = 4, n_cores = 1, control = list()){

  ## Error catching
  if (!identical(dim(Y_source), dim(Y_target))){
    stop('Y_source and Y_target must have the same dimensions')
  }
  if (any(is.na(Y_source))){
    stop('Y_source cannot have NA values.')
  }

  ## Setting the rank and other parameters
  p <- nrow(Y_source)
  q <- ncol(Y_source)
  if (missing(r)){
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

  svd_source <- svd(Y_source, nu = r, nv = r)
  P_U <- diag(rep(1, p)) - svd_source$u %*% t(svd_source$u)
  P_V <- diag(rep(1, q)) - svd_source$v %*% t(svd_source$v)
  U_init <- svd_source$u %*% diag(sqrt(svd_source$d[1:r]), nrow = r, ncol = r)
  V_init <- svd_source$v %*% diag(sqrt(svd_source$d[1:r]), nrow = r, ncol = r)

  if (n_folds < 2){
    stop('n_folds must be 2 or greater.')
  }
  n_lambda_1 <- length(lambda_1_all)
  n_lambda_2 <- length(lambda_2_all)

  ## Creating training and testing data sets
  n_indices <- p * q
  indices <- sample(1:(p * q), size = n_indices, replace = FALSE)

  index_set <- vector(mode = "list", length = n_folds)
  for (fold in 1:n_folds){
    index_set[[fold]] <- indices[floor((fold - 1) * n_indices / n_folds + 1):floor(fold * n_indices / n_folds)]
  }

  ## Applying LEARNER for each candidate lambda_1 and lambda_2 values
  if (n_cores > 1){
    ## Parallelization
    registerDoParallel(cores = n_cores)
    mse_all <- foreach(lambda_1_ind = 1:n_lambda_1, .combine = 'cbind') %:%
      foreach(lambda_2_ind = 1:n_lambda_2, .combine = 'c') %dopar% {
        cv.learner_helper(Y_target = Y_target, index_set = index_set, n_folds = n_folds,
                          lambda1 = lambda_1_all[lambda_1_ind], lambda2 = lambda_2_all[lambda_2_ind],
                          U_init = U_init, V_init = V_init,
                          P_U = P_U, P_V = P_V, step_size = step_size,
                          control = control)
      }
  } else {
    ## Sequential
    mse_all <- foreach(lambda_1_ind = 1:n_lambda_1, .combine = 'cbind') %:%
      foreach(lambda_2_ind = 1:n_lambda_2, .combine = 'c') %do% {
        cv.learner_helper(Y_target = Y_target, index_set = index_set, n_folds = n_folds,
                          lambda1 = lambda_1_all[lambda_1_ind], lambda2 = lambda_2_all[lambda_2_ind],
                          U_init = U_init, V_init = V_init,
                          P_U = P_U, P_V = P_V, step_size = step_size,
                          control = control)
      }
  }
  mse_all <- t(mse_all)

  min_ind <- which(mse_all == min(mse_all, na.rm = TRUE), arr.ind = TRUE)
  lambda_1_min <- lambda_1_all[min_ind[1]]
  lambda_2_min <- lambda_2_all[min_ind[2]]

  return(list(
    lambda_1_min = lambda_1_min,
    lambda_2_min = lambda_2_min,
    mse_all = mse_all,
    r = r
  ))
}

cv.learner_helper <- function(Y_target, index_set, n_folds,
                              lambda1, lambda2,
                              U_init, V_init, P_U, P_V,
                              step_size, control){
  norm_temp <- 0
  p <- nrow(Y_target); q <- ncol(Y_target)
  for (fold in 1:n_folds){
    Y_training <- Y_target
    myset <- index_set[[fold]]
    Y_training[myset] <- NA
    perc_nonmissing <- sum(!is.na(Y_training)) / (p * q)

    fit <- try(learner_helper(theta_hat = Y_training,
                              lambda1 = lambda1,
                              lambda2 = lambda2,
                              U_init = U_init, V_init = V_init,
                              P_U = P_U, P_V = P_V, missing = TRUE,
                              perc_nonmissing = perc_nonmissing,
                              step_size = step_size, control = control))
    if ('try-error' %in% class(fit)){
      return(NA)
    } else {
      norm_temp <- norm_temp + sum((fit$learner_estimate[myset] - Y_target[myset])^2)
    }
  }
  return(norm_temp)
}
