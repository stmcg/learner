#' Direct projection transfer learning
#'
#' Estimate a target signal matrix by applying source-space operators directly
#' to the target data, using the D-LEARNER method (McGrath et al. 2024).
#' One source matrix or a list of source matrices can be supplied.
#'
#' @param Y_target nonempty numeric target matrix with only finite values.
#' Unlike \code{learner()}, this function does not allow missing target entries.
#' @param r optional common rank of the retained source spaces. If omitted, ScreeNOT selects it from the first source.
#' @inheritParams learner
#'
#' @return A list with the following components:
#' \item{dlearner_estimate}{matrix containing the D-LEARNER estimate of the target population knowledge graph.}
#' \item{r}{rank value used.}
#'
#' @details
#'
#' Each source has the same dimensions as the target matrix \eqn{Y_0}.
#' Let \eqn{B_{U,k}} and \eqn{B_{V,k}} contain the first \eqn{r} left and
#' right singular vectors of source \eqn{k}. The estimate is
#' \deqn{\widehat\Theta_0 = \bar P_U Y_0 \bar P_V, \qquad
#' \bar P_U=\sum_k w_k B_{U,k}B_{U,k}^\top, \qquad
#' \bar P_V=\sum_k w_k B_{V,k}B_{V,k}^\top.}
#' The nonnegative weights are normalized to sum to one and default to equal
#' values. With one source this reduces to projection onto its retained left
#' and right singular spaces. With multiple sources, the averaged operators
#' need not be idempotent, and the estimate can have rank larger than \code{r}.
#' This is not an average of the raw source matrices or of separately computed
#' single-source D-LEARNER estimates.
#'
#' D-LEARNER does not iteratively fit target factors and has no space or balance
#' penalty parameters. All inputs must have aligned rows and columns and no
#' missing values. Names are not matched. If \code{r} is omitted, it is selected
#' from the first source, including when that source has zero weight.
#' Output dimension names are taken from the first source.
#'
#' @references
#' Donoho, D., Gavish, M. and Romanov, E. (2023). \emph{ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise}. The Annals of Statistics, 51(1), pp.122-148.
#'
#' @examples
#' res <- dlearner(Y_source = dat_highsim$Y_source,
#'                 Y_target = dat_highsim$Y_target)
#'
#' # A second simulated source, made by perturbing the first.
#' set.seed(803)
#' source_A <- dat_highsim$Y_source[1:20, 1:10]
#' source_B <- source_A + matrix(rnorm(length(source_A), sd = 0.3),
#'                               nrow(source_A), ncol(source_A))
#' res_multi <- dlearner(
#'   Y_source = list(A = source_A, B = source_B),
#'   Y_target = dat_highsim$Y_target[1:20, 1:10], r = 2,
#'   source_weights = c(0.7, 0.3))
#' dim(res_multi$dlearner_estimate)
#'
#' @export

dlearner <- function(Y_source, Y_target, r, source_weights = NULL) {
  sources <- prepare_sources(Y_source, Y_target, source_weights,
                             allow_missing_target = FALSE)
  r <- resolve_rank(if (missing(r)) NULL else r, sources$matrices[[1]])
  bases <- lapply(sources$matrices, svd, nu = r, nv = r)
  if (length(bases) == 1L) {
    # Preserve the original single-source multiplication order.
    s <- bases[[1]]
    estimate <- s$u %*% (t(s$u) %*% Y_target %*% s$v) %*% t(s$v)
  } else {
    # Average source projectors, not raw source matrices.
    middle <- matrix(0, nrow(Y_target), ncol(Y_target))
    for (k in seq_along(bases)) {
      u <- bases[[k]]$u
      middle <- middle + sources$weights[k] * (u %*% crossprod(u, Y_target))
    }
    estimate <- matrix(0, nrow(Y_target), ncol(Y_target))
    for (k in seq_along(bases)) {
      v <- bases[[k]]$v
      estimate <- estimate + sources$weights[k] * ((middle %*% v) %*% t(v))
    }
  }
  dimnames(estimate) <- dimnames(sources$matrices[[1]])
  list(dlearner_estimate = estimate, r = r)
}
