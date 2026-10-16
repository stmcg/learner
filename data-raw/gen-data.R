################################################################################
## High similarity data
################################################################################

set.seed(1234)
p <- 100; q <- 50; r <- 3

# Generating Theta_target
Theta_target <- matrix(rnorm(p * r), p, r) %*% matrix(rnorm(r * q), r, q)

# Generating Theta_source
Theta_target_svd <- svd(Theta_target, nu = r, nv = r)
Theta_target_singular_values <- Theta_target_svd$d[1:r]
Theta_source <- Theta_target_svd$u %*% diag(rev(Theta_target_singular_values), nrow = r, ncol = r) %*% t(Theta_target_svd$v)

# Generating Y_source and Y_target
Y_source <- Theta_source + matrix(rnorm(p * q, sd = 0.5), p, q)
Y_target <- Theta_target + matrix(rnorm(p * q, sd = 1), p, q)

dat_highsim <- list(Y_source = Y_source,
                    Y_target = Y_target,
                    Theta_source = Theta_source,
                    Theta_target = Theta_target)

usethis::use_data(dat_highsim, overwrite = TRUE)

################################################################################
## Moderate similarity data
################################################################################

set.seed(12345)

# Generating Theta_target
Theta_target <- matrix(rnorm(p * r), p, r) %*% matrix(rnorm(r * q), r, q)

# Generating Theta_source
Theta_target_svd <- svd(Theta_target, nu = r, nv = r)
Theta_target_singular_values <- Theta_target_svd$d[1:r]

threshold_p <- 1/(8 * sqrt(p)); threshold_q <- 1/(8 * sqrt(q))
delta1 <- matrix(runif(p * r, min = -threshold_p, max = threshold_p), nrow = p, ncol = r)
delta2 <- matrix(runif(q * r, min = -threshold_q, max = threshold_q), nrow = q, ncol = r)
u_new <-  qr.Q(qr(Theta_target_svd$u + delta1))
v_new <- qr.Q(qr(Theta_target_svd$v + delta2))

Theta_source <- u_new %*% diag(rev(Theta_target_singular_values), nrow = r, ncol = r) %*% t(v_new)

# Generating Y_source and Y_target
Y_source <- Theta_source + matrix(rnorm(p * q, sd = 0.5), p, q)
Y_target <- Theta_target + matrix(rnorm(p * q, sd = 1), p, q)

dat_modsim <- list(Y_source = Y_source,
                   Y_target = Y_target,
                   Theta_source = Theta_source,
                   Theta_target = Theta_target)

usethis::use_data(dat_modsim, overwrite = TRUE)

################################################################################
## Multiple sources with different noise levels
################################################################################

# This block can be run on its own without regenerating the existing datasets.
local({
  set.seed(20261014)
  p <- 100; q <- 50; r <- 3

  Theta_target <- matrix(rnorm(p * r), p, r) %*%
    matrix(rnorm(r * q), r, q)
  target_svd <- svd(Theta_target, nu = r, nv = r)
  # Both sources share the target's latent spaces, with reversed signal strengths.
  Theta_source <- target_svd$u %*%
    diag(rev(target_svd$d[seq_len(r)]), nrow = r, ncol = r) %*%
    t(target_svd$v)

  # Variances, not standard deviations. Each call draws independent noise.
  noise_variance <- c(target = 0.10, A = 0.01, B = 0.50)
  Y_target <- Theta_target + matrix(
    rnorm(p * q, sd = sqrt(noise_variance['target'])), p, q)
  Y_A <- Theta_source + matrix(
    rnorm(p * q, sd = sqrt(noise_variance['A'])), p, q)
  Y_B <- Theta_source + matrix(
    rnorm(p * q, sd = sqrt(noise_variance['B'])), p, q)

  dat_multisource <- list(
    Y_source = list(A = Y_A, B = Y_B),
    Y_target = Y_target,
    Theta_source = list(A = Theta_source, B = Theta_source),
    Theta_target = Theta_target
  )
  usethis::use_data(dat_multisource, overwrite = TRUE, version = 2)
})
