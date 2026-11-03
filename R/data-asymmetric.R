#' Simulated data set: One shared latent space
#'
#' One source and one target with the same true U-space but different true
#' V-spaces. This example illustrates how to fit LEARNER
#' with separate penalties on the two spaces.
#'
#' @docType data
#' @format A list of four aligned matrices, each with 100 rows and 50 columns:
#' \describe{
#'   \item{Y_source}{Observed source data.}
#'   \item{Y_target}{Observed target data.}
#'   \item{Theta_source}{True rank-3 source signal.}
#'   \item{Theta_target}{True rank-3 target signal.}
#' }
#' @details The target signal is the rank-3 truncated SVD of a standard normal
#' matrix. The source keeps its U-space, reverses its three singular values,
#' and perturbs its V basis with independent uniform entries between
#' \eqn{-1/(2\sqrt{50})} and \eqn{1/(2\sqrt{50})}, followed by QR orthonormalization.
#' This adapts the U-shared, V-perturbed "low similarity" scenario in Section
#' 2.4.1 of the Supporting Information of McGrath et al. (2026) to a small
#' example (the paper uses 500 by 500 matrices of rank 4). It is not a
#' replication of that study. The Supporting Information calls this setting
#' "same column spaces but different row spaces": the span of U is the column
#' space of the signal matrix, and the span of V is its row space.
#'
#' Independent mean-zero normal noise is added to each signal, with variance
#' 0.10 for the target and 0.01 for the source. The exact shared-space statement
#' applies to the true signals; noise also affects estimated spaces.
#' True signals are supplied to illustrate estimation error, not for fitting
#' or penalty selection. Generation code is in \file{data-raw/gen-asymmetric.R}
#' in the source repository.
#'
#' @references
#' McGrath, S., Zhu, C., O'Dea, R., Guo, M. and Duan, R. (2026). \emph{LEARNER: A transfer learning method for low-rank matrix estimation}. Statistics in Medicine, 45(15-17), e70669. \doi{10.1002/sim.70669}.
#'
#' @examples
#' data(dat_asymmetric)
#' d <- dat_asymmetric
#' # Demonstration penalties; use cv.learner() to select them from data.
#' fit <- learner(d$Y_source, d$Y_target, r = 3, lambda_1_row = 100,
#'                lambda_1_col = 10, lambda_2 = 0.1, step_size = 0.003,
#'                control = list(max_iter = 2000))
#' dim(fit$learner_estimate)
#' fit$convergence_criterion
#' # See vignette("advanced-learner") for cross-validation and diagnostics.
#' @seealso \code{\link{dat_highsim}}, \code{\link{dat_multisource}},
#'   \code{\link{learner}}, \code{\link{cv.learner}}
#' @keywords datasets
"dat_asymmetric"
