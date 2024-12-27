test_that("cv.learner (sequential) does not fail", {
  expect_no_error(
    cv.learner(Y_source = dat_highsim$Y_source,
               Y_target = dat_highsim$Y_target,
               lambda_1_all = c(1, 10),
               lambda_2_all = c(1, 10),
               step_size = 0.003))
})

test_that("cv.learner (parallel) does not fail", {
  expect_no_error(
    cv.learner(Y_source = dat_highsim$Y_source,
               Y_target = dat_highsim$Y_target,
               lambda_1_all = c(1, 10),
               lambda_2_all = c(1, 10),
               step_size = 0.003,
               n_cores = 2))
})

test_that("learner does not fail", {
  expect_no_error(
    learner(Y_source = dat_highsim$Y_source,
            Y_target = dat_highsim$Y_target,
            lambda_1 = 1, lambda_2 = 1,
            step_size = 0.003))
})

test_that("dlearner does not fail", {
  expect_no_error(
    dlearner(Y_source = dat_highsim$Y_source,
             Y_target = dat_highsim$Y_target))
})

