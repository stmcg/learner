#' Direct projection transfer learning
#'
#' @description
#' This function applies the Direct projection LatEnt spAce-based tRaNsfer
#' lEaRning (D-LEARNER) method (McGrath et al. 2024) to leverage data from a source
#' population to improve estimation of a low-rank matrix in an underrepresented
#' target population.
#'
#' @param Y_source Matrix containing the source population data, with the same rows
#' and columns as \code{Y_target}. For multiple sources, supply a list of
#' matrices.
#' @param Y_target Matrix containing the target population data. Missing entries
#' are not allowed.
#' @param r Rank of the knowledge graphs. If omitted, ScreeNOT (Donoho et al. 2023)
#' is applied to the source population data to select the rank.
#' @param source_weights Optional numeric vector of nonnegative relative weights,
#' with one value per source in list order and at least one positive value.
#' The default, \code{NULL}, estimates weights from the source and target spaces.
#' Numeric weights are normalized to sum to one; use equal values for equal
#' weighting. The strings \code{"auto"} and \code{"equal"} remain supported for
#' compatibility. See Details for the empirical weighting rule.
#'
#' @details
#' \strong{Data and notation:}
#'
#' The data consists of a matrix in the target population
#' \eqn{Y_0 \in \mathbb{R}^{p \times q}} and the source population
#' \eqn{Y_1 \in \mathbb{R}^{p \times q}}. Let
#' \eqn{\hat U_k\hat\Lambda_k\hat V_k^\top} denote the rank-\eqn{r} truncated
#' singular value decomposition (SVD) of \eqn{Y_k}, \eqn{k = 0, 1}.
#'
#' For \eqn{k = 0, 1}, one can view \eqn{Y_k} as a noisy version of
#' \eqn{\Theta_k}, referred to as the knowledge graph. The target of inference
#' is the target population knowledge graph, \eqn{\Theta_0}.
#'
#' \strong{Estimation:}
#'
#' This method estimates \eqn{\Theta_0} by
#'
#' \deqn{\hat U_1\hat U_1^\top Y_0 \hat V_1\hat V_1^\top.}
#'
#' That is, the target population data are projected onto the latent row and
#' column spaces of the source population, without iterative fitting or penalty
#' selection.
#'
#' \strong{Multiple sources:}
#'
#' Supply a list of source matrices to combine information from several
#' populations. Their projection operators are averaged using
#' \code{source_weights}, and the estimate becomes
#'
#' \deqn{\widehat\Theta_0=\bar P_U Y_0\bar P_V,\qquad
#' \bar P_U=\sum_{k=1}^K w_k\hat U_k\hat U_k^\top,\qquad
#' \bar P_V=\sum_{k=1}^K w_k\hat V_k\hat V_k^\top.}
#'
#' where \eqn{\hat U_k} and \eqn{\hat V_k} come from the truncated SVD of
#' the \eqn{k}th source and the weights \eqn{w_k} are given by
#' \code{source_weights}. If \code{r} is not specified, the rank is selected
#' using the first source.
#'
#' \strong{Empirical source weights:}
#'
#' With \code{source_weights = NULL}, weights follow Equation (3) of McGrath
#' et al. (2026): \eqn{w_k \propto 1/(d_k+10^{-12})}, where
#' \deqn{d_k=\|\hat U_k\hat U_k^\top-\hat U_0\hat U_0^\top\|_F^2+
#' \|\hat V_k\hat V_k^\top-\hat V_0\hat V_0^\top\|_F^2.}
#' All bases use rank \code{r}. The small constant stabilizes coincident spaces.
#' A single source receives weight one. The same weights apply to both spaces.
#' D-LEARNER requires a complete target matrix. Supply numeric weights to choose
#' a fixed combination.
#'
#' @return
#' A list with the following components:
#' \item{dlearner_estimate}{matrix containing the D-LEARNER estimate of the target population knowledge graph.}
#' \item{r}{rank value used.}
#' \item{source_weights}{Normalized source weights used for this fit.}
#' \item{source_weight_method}{Weighting rule: \code{"auto"} for empirical weights,
#' \code{"equal"} for the compatibility option, or \code{"manual"} for numeric weights.}
#'
#' @references
#' McGrath, S., Zhu, C., O'Dea, R., Guo, M. and Duan, R. (2026).
#' LEARNER: A transfer learning method for low-rank matrix estimation.
#' Statistics in Medicine, 45, e70669. doi:10.1002/sim.70669.
#'
#' Donoho, D., Gavish, M. and Romanov, E. (2023). \emph{ScreeNOT: Exact MSE-optimal singular value thresholding in correlated noise}. The Annals of Statistics, 51(1), pp.122-148.
#'
#' @seealso
#' \code{\link{learner}} for a target fit that can adapt the source
#' spaces through penalized estimation.
#'
#' @examples
#' res <- dlearner(Y_source = dat_highsim$Y_source,
#'                 Y_target = dat_highsim$Y_target)
#'
#' # Multiple sources with fixed demonstration weights.
#' res_multi <- dlearner(
#'   Y_source = dat_multisource$Y_source,
#'   Y_target = dat_multisource$Y_target, r = 3,
#'   source_weights = c(0.7, 0.3))
#' dim(res_multi$dlearner_estimate)
#'
#' @export

dlearner <- function(Y_source, Y_target, r, source_weights = NULL) {
  sources <- prepare_sources(Y_source, Y_target, source_weights,
                             allow_missing_target = FALSE)
  r <- resolve_rank(if (missing(r)) NULL else r, sources$matrices[[1]])
  bases <- lapply(sources$matrices, svd, nu = r, nv = r)
  sources$weights <- resolve_source_weights(sources, Y_target, r, bases)
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
  list(dlearner_estimate = estimate, r = r, source_weights = sources$weights,
       source_weight_method = sources$mode)
}
