#' Latent space-based transfer learning
#'
#' This function applies the Direct project LatEnt spAce-based tRaNsfer lEaRning (D-LEARNER) method (McGrath et al. 2024) to leverage data from a source population to improve
#' estimation of a low rank matrix in an underrepresented target population.
#'
#' @param Y_target matrix containing the target population data
#' @param Y_source matrix containing the source population data
#' @param r (optional) integer specifying the rank of the knowledge graphs. By default, ScreeNOT (Donoho et al. 2023) is applied to the source population knowledge graph to select the rank.
#'
#' @return A list with the following components:
#' \item{dlearner_estimate}{matrix containing the D-LEARNER estimate of the target population knowledge graph.}
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
#' This method estimates \eqn{\Theta_0} by \eqn{\hat{U}_{1}^{\top}\hat{U}_{1} Y_0 \hat{V}_{1}^{\top}\hat{V}_{1}}.
#'
#' @references
#' Donoho, D., Gavish, M. and Romanov, E. (2023). \emph{ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise}. The Annals of Statistics, 51(1), pp.122-148.
#'
#' @examples
#' res <- dlearner(Y_source = dat_highsim$Y_source,
#'                 Y_target = dat_highsim$Y_target)
#'
#' @export

dlearner <- function(Y_source, Y_target, r){
  # Error catching
  if (!identical(dim(Y_source), dim(Y_target))){
    stop('Y_source and Y_target must have the same dimensions')
  }
  if (any(is.na(Y_source))){
    stop('Y_source cannot have NA values.')
  }
  if (any(is.na(Y_target))){
    stop('Y_target cannot have NA values')
  }

  p <- nrow(Y_source)
  q <- ncol(Y_source)

  if (missing(r)){
    max_rank <- min(p, q) / 3
    r <- max(ScreeNOT::adaptiveHardThresholding(Y = Y_source, k = max_rank)$r, 1)
  }

  svd_source <- svd(Y_source, nu = r, nv = r)
  dlearner_estimate <- (svd_source$u %*% t(svd_source$u)) %*%
    Y_target %*% (svd_source$v %*% t(svd_source$v))

  colnames(dlearner_estimate) <- colnames(Y_source)
  rownames(dlearner_estimate) <- rownames(Y_source)

  return(list(dlearner_estimate = dlearner_estimate, r = r))
}

