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
