# Adapted from simulations/helper-uhigh-vlow.R, set_theta_Ushared_Vdiff(),
# in stmcg/learner-paper, commit cb092c7e239b9a7b8eb8fcea80bce26f8017be86.
# Uses dif_factor_v = 2 from runsim-uhigh-vlow.R. Smaller dimensions and fixed
# rank make this a package example, not a reproduction of the simulation study.
local({
  set.seed(20261022)
  p <- 100L; q <- 50L; r <- 3L
  basis <- svd(matrix(rnorm(p * q), p, q), nu = r, nv = r)
  strengths <- basis$d[seq_len(r)]
  Theta_target <- basis$u %*% diag(strengths) %*% t(basis$v)
  threshold_v <- 1 / (2 * sqrt(q))
  delta_v <- matrix(runif(q * r, -threshold_v, threshold_v), q, r)
  source_v <- qr.Q(qr(basis$v + delta_v))
  Theta_source <- basis$u %*% diag(rev(strengths)) %*% t(source_v)
  # Independent population noise; arguments to rnorm() are standard deviations.
  Y_target <- Theta_target + matrix(rnorm(p * q, sd = sqrt(0.10)), p, q)
  Y_source <- Theta_source + matrix(rnorm(p * q, sd = sqrt(0.01)), p, q)
  dat_asymmetric <- list(Y_source = Y_source, Y_target = Y_target,
                         Theta_source = Theta_source, Theta_target = Theta_target)
  usethis::use_data(dat_asymmetric, overwrite = TRUE, version = 2)
})
