test_that('multiple-source example contains aligned observations and rank-three signals', {
  d <- dat_multisource
  expect_named(d, c('Y_source', 'Y_target', 'Theta_source', 'Theta_target'))
  expect_named(d$Y_source, c('A', 'B'))
  expect_named(d$Theta_source, c('A', 'B'))
  matrices <- c(d$Y_source, list(d$Y_target), d$Theta_source, list(d$Theta_target))
  for (m in matrices) {
    expect_identical(dim(m), c(100L, 50L))
    expect_true(all(is.finite(m)))
  }
  for (m in c(d$Theta_source, list(d$Theta_target))) {
    expect_equal(sum(svd(m, nu = 0, nv = 0)$d > 1e-8), 3)
  }
  expect_equal(d$Theta_source$A, d$Theta_source$B)
  target <- svd(d$Theta_target, nu = 3, nv = 3)
  source <- svd(d$Theta_source$A, nu = 3, nv = 3)
  expect_equal(tcrossprod(target$u), tcrossprod(source$u), tolerance = 1e-10)
  expect_equal(tcrossprod(target$v), tcrossprod(source$v), tolerance = 1e-10)
})

test_that('saved example has the intended noise scales without shared source noise', {
  d <- dat_multisource
  noise <- cbind(target = as.vector(d$Y_target - d$Theta_target),
                 A = as.vector(d$Y_source$A - d$Theta_source$A),
                 B = as.vector(d$Y_source$B - d$Theta_source$B))
  # Fixed saved data: broad checks catch variance/SD mistakes and reused noise.
  ratios <- apply(noise, 2, var) / c(target = 0.10, A = 0.01, B = 0.50)
  expect_true(all(ratios > 0.8 & ratios < 1.2))
  correlations <- cor(noise)
  expect_true(all(abs(correlations[upper.tri(correlations)]) < 0.05))
})
