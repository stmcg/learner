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
                    control = list()){
  ## Error catching
  if (!identical(dim(Y_source), dim(Y_target))){
    stop('Y_source and Y_target must have the same dimensions')
  }
  if (any(is.na(Y_source))){
    stop('Y_source cannot have NA values.')
  }

  ## Setting the rank and other parameters
  missing <- any(is.na(Y_target))
  p <- nrow(Y_source)
  q <- ncol(Y_source)
  perc_nonmissing <- sum(!is.na(Y_target)) / (p * q)
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

  out <- learner_helper(control = control, U_init = U_init, V_init = V_init, lambda1 = lambda_1,
                        lambda2 = lambda_2, theta_hat = Y_target, P_U = P_U, P_V = P_V, missing = missing,
                        perc_nonmissing = perc_nonmissing, step_size = step_size)

  return(out)
}

learner_helper <- function(control, U_init, V_init, lambda1, lambda2, theta_hat, P_U, P_V, missing,
                           perc_nonmissing, step_size){
  U <- U_init
  V <- V_init
  max_iter <- control$max_iter
  threshold <- control$threshold
  max_value <- control$max_value

  obj_func_init <- f(U = U, V = V, lambda1 = lambda1, lambda2 = lambda2, theta_hat = theta_hat, P_U = P_U, P_V = P_V, missing = missing,
                     perc_nonmissing = perc_nonmissing)
  obj_func_best <- obj_func_init; U_best <- U; V_best <- V
  obj_func <- rep(NA, times = max_iter)
  U_norm <- norm(U, type = 'F')
  V_norm <- norm(V, type = 'F')

  for (iter in 1:max_iter){
    f_delta_U <- f_prime_u(U = U, V = V, lambda1 = lambda1, lambda2 = lambda2, theta_hat = theta_hat, P_U = P_U, missing = missing, perc_nonmissing = perc_nonmissing)
    f_delta_V <- f_prime_v(U = U, V = V, lambda1 = lambda1, lambda2 = lambda2, theta_hat = theta_hat, P_V = P_V, missing = missing, perc_nonmissing = perc_nonmissing)

    f_delta_U_norm <- norm(f_delta_U, type = 'F')
    f_delta_V_norm <- norm(f_delta_V, type = 'F')

    U <- U - (step_size * U_norm / f_delta_U_norm) * f_delta_U
    V <- V - (step_size * V_norm / f_delta_V_norm) * f_delta_V

    obj_func[iter] <- f(U = U, V = V, lambda1 = lambda1, lambda2 = lambda2, theta_hat = theta_hat, P_U = P_U, P_V = P_V, missing = missing, perc_nonmissing = perc_nonmissing)
    if (obj_func[iter] < obj_func_best){
      obj_func_best <- obj_func[iter]; U_best <- U; V_best <- V
    }
    if (iter > 1){
      if (abs(obj_func[iter] - obj_func[iter - 1]) < threshold){
        convergence_criterion <- 1
        break()
      }
      if (obj_func[iter] > max_value * obj_func_init){
        convergence_criterion <- 3
        break()
      }
    }
    U_norm <- norm(U, type = 'F')
    V_norm <- norm(V, type = 'F')

    if (iter == max_iter){
      convergence_criterion <- 2
    }
  }
  theta_hat <- U_best %*% t(V_best)
  return(list(learner_estimate = theta_hat, objective_values = obj_func,
              convergence_criterion = convergence_criterion,
              r = ncol(U)))
}


norm_2_sq <- function(x){
  return(sum(x^2))
}

f <- function(U, V, lambda1, lambda2, theta_hat, P_U, P_V, missing, perc_nonmissing){
  if (missing){
    t1 <- sum((U %*% t(V) - theta_hat)^2, na.rm = TRUE) / perc_nonmissing
  } else {
    t1 <- norm(U %*% t(V) - theta_hat, type = 'F')^2 / perc_nonmissing
  }
  t2 <- lambda1 * norm(P_U %*% U, type = 'F')^2
  t3 <- lambda1 * norm(P_V %*% V, type = 'F')^2
  t4 <- lambda2 * norm(t(U) %*% U - t(V) %*% V, type = 'F')^2
  return(t1 + t2 + t3 + t4)
}

f_prime_u <- function(U, V, lambda1, lambda2, theta_hat, P_U, missing, perc_nonmissing){
  U_tilde <- t(U) %*% U
  V_tilde <- t(V) %*% V

  if (missing){
    temp <- U %*% t(V)
    missing_entries <- which(is.na(theta_hat), arr.ind = T)
    theta_hat[missing_entries] <- temp[missing_entries]
  }
  t1 <- (2 / perc_nonmissing) * (U %*% V_tilde - theta_hat %*% V)
  t2 <- lambda1 * 2 * P_U %*% U
  t4 <- lambda2 * 4 * U %*% (U_tilde - V_tilde)
  return(t1 + t2 + t4)
}

f_prime_v <- function(U, V, lambda1, lambda2, theta_hat, P_V, missing, perc_nonmissing){
  U_tilde <- t(U) %*% U
  V_tilde <- t(V) %*% V

  if (missing){
    temp <- U %*% t(V)
    missing_entries <- which(is.na(theta_hat), arr.ind = T)
    theta_hat[missing_entries] <- temp[missing_entries]
  }
  t1 <- (2 / perc_nonmissing) * (V %*% U_tilde - t(theta_hat) %*% U)
  t3 <- lambda1 * 2 * P_V %*% V
  t4 <- lambda2 * 4 * V %*% (V_tilde - U_tilde)
  return(t1 + t3 + t4)
}
