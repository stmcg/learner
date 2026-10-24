test_that('asymmetric example shares only the true U-space', {
  d <- dat_asymmetric
  expect_named(d, c('Y_source', 'Y_target', 'Theta_source', 'Theta_target'))
  for (m in d) {
    expect_identical(dim(m), c(100L, 50L))
    expect_true(all(is.finite(m)))
  }
  target <- svd(d$Theta_target)
  source <- svd(d$Theta_source)
  expect_equal(sum(target$d > 1e-8), 3)
  expect_equal(sum(source$d > 1e-8), 3)
  expect_equal(tcrossprod(target$u[, 1:3]), tcrossprod(source$u[, 1:3]),
               tolerance = 1e-10)
  expect_gt(norm(tcrossprod(target$v[, 1:3]) - tcrossprod(source$v[, 1:3]), 'F'), 0.1)
  noise <- cbind(target = as.vector(d$Y_target - d$Theta_target),
                 source = as.vector(d$Y_source - d$Theta_source))
  ratios <- apply(noise, 2, var) / c(0.10, 0.01)
  expect_true(all(ratios > 0.8 & ratios < 1.2))
  expect_lt(abs(cor(noise)[1, 2]), 0.05)
})

test_that('source reverses the pairing of target singular values at the low-similarity V scale', {
  d <- dat_asymmetric
  target <- svd(d$Theta_target, nu = 3, nv = 3)
  source <- svd(d$Theta_source, nu = 3, nv = 3)
  # Reversing the pairing leaves the set of singular values unchanged; a dropped
  # QR step or a rescaled signal would not.
  expect_equal(source$d[1:3], target$d[1:3], tolerance = 1e-8)
  # The top three singular values of a 100 x 50 standard normal matrix lie in
  # about 14.5 to 18 (2000 simulated draws).
  expect_true(all(target$d[1:3] > 14 & target$d[1:3] < 19))
  # The strength the source carries on each shared left direction is the
  # reversed target strength.
  strength_on_u <- sqrt(rowSums(crossprod(target$u[, 1:3], d$Theta_source)^2))
  expect_equal(strength_on_u, rev(target$d[1:3]), tolerance = 1e-8)
  # The V perturbation is at the low-similarity scale: 1/(4 sqrt(q)) gives about
  # 0.34 and 1/sqrt(q) about 1.1, against about 0.65 for 1/(2 sqrt(q)).
  d_v <- norm(tcrossprod(target$v[, 1:3]) - tcrossprod(source$v[, 1:3]), 'F')
  expect_gt(d_v, 0.5)
  expect_lt(d_v, 0.8)
})

test_that('noise is mean-zero with the documented scales', {
  d <- dat_asymmetric
  noise <- cbind(target = as.vector(d$Y_target - d$Theta_target),
                 source = as.vector(d$Y_source - d$Theta_source))
  se <- sqrt(c(0.10, 0.01) / nrow(noise))
  expect_true(all(abs(colMeans(noise)) < 5 * se))
  ratios <- apply(noise, 2, var) / c(0.10, 0.01)
  expect_true(all(abs(ratios - 1) < 5 * sqrt(2 / nrow(noise))))
})
