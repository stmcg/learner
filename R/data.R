#' Simulated data set: High similarity in the latent spaces
#'
#' This data set contains simulated data in the source and target populations where there is a high degree of similarity in the underlying latent spaces between these populations.
#'
#' @docType data
#'
#' @format A list containing the observed and true matrices in the source and target populations. The list contains the following components:
#' \describe{
#'   \item{\code{Y_source}}{A matrix of size \eqn{100 \times 50} representing the observed source population matrix.}
#'   \item{\code{Y_target}}{A matrix of size \eqn{100 \times 50} representing the observed target population matrix.}
#'   \item{\code{Theta_source}}{A matrix of size \eqn{100 \times 50} (rank 3) representing the true source population matrix.}
#'   \item{\code{Theta_target}}{A matrix of size \eqn{100 \times 50} (rank 3) representing the true target population matrix.}
#' }
#' @details
#' In this data set, there is a high degree of similarity in the latent spaces between the source and target populations. Specifically, the true source population matrix was obtained by reversing the order of the singular values of the true target population matrix.
#' The observed target population matrix was obtained by adding independent and identically distributed noise to the entries of the true source population matrix. The noise was generated from a normal distribution with mean 0 and standard deviation of 1. The observed source population matrix was generated analogously, where the noise had a standard deviation of 0.5.
#'
#' @seealso \code{\link{dat_modsim}}
#'
#' @keywords datasets
"dat_highsim"

#' Simulated data set: Moderate similarity in the latent spaces
#'
#' This data set contains simulated data in the source and target populations where there is a moderate degree of similarity in the underlying latent spaces between these populations.
#'
#' @docType data
#'
#' @format A list containing the observed and true matrices in the source and target populations. The list contains the following components:
#' \describe{
#'   \item{\code{Y_source}}{A matrix of size \eqn{100 \times 50} representing the observed source population matrix.}
#'   \item{\code{Y_target}}{A matrix of size \eqn{100 \times 50} representing the observed target population matrix.}
#'   \item{\code{Theta_source}}{A matrix of size \eqn{100 \times 50} (rank 3) representing the true source population matrix.}
#'   \item{\code{Theta_target}}{A matrix of size \eqn{100 \times 50} (rank 3) representing the true target population matrix.}
#' }
#' @details
#' In this data set, there is a moderate degree of similarity in the latent spaces between the source and target populations. Specifically, the true source population matrix was obtained by (i) reversing the order of the singular values of the true target population matrix and (ii) adding perturbations to the left and right singular vectors of the true target population matrix.
#' The observed target population matrix was obtained by adding independent and identically distributed noise to the entries of the true source population matrix. The noise was generated from a normal distribution with mean 0 and standard deviation of 1. The observed source population matrix was generated analogously, where the noise had a standard deviation of 0.5.
#'
#' @seealso \code{\link{dat_modsim}}
#'
#' @keywords datasets
"dat_modsim"
