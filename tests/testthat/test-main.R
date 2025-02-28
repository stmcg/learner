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

test_that("cv.learner output is correct", {
  expected_result <- t(matrix(c(5494.916, 5312.724,
                                5341.662, 5233.828,
                                5009.216, 5013.827), nrow = 2, ncol = 3))
  set.seed(1234)
  result <- cv.learner(Y_source = dat_highsim$Y_source,
                       Y_target = dat_highsim$Y_target,
                       lambda_1_all = c(1, 10, 100),
                       lambda_2_all = c(1, 10),
                       step_size = 0.003)
  expect_equal(result$mse, expected_result, tolerance = 1e-3)
})

test_that("learner does not fail", {
  expect_no_error(
    learner(Y_source = dat_highsim$Y_source,
            Y_target = dat_highsim$Y_target,
            lambda_1 = 1, lambda_2 = 1,
            step_size = 0.003))
})

test_that("learner output is correct", {
  expected_result <- t(matrix(c(0.1578405, -1.54883716,  1.2205437,  0.03377010,  0.0767873,
                                0.1189457, -0.12793612, -0.9755801,  0.01372606, -0.3902793,
                                -0.4384613,  1.21518704, -1.0718671,  0.60388502, -0.4104101,
                                1.4426139, -3.13260948,  0.9175863, -1.79041482,  0.5015381,
                                0.5602868, -0.05068081, -1.3878535, -0.78502749, -0.1315952), nrow = 5, ncol = 5))
  result <- learner(Y_source = dat_highsim$Y_source,
                    Y_target = dat_highsim$Y_target,
                    lambda_1 = 1, lambda_2 = 1,
                    step_size = 0.003)
  expect_equal(result$learner_estimate[1:5, 1:5], expected_result, tolerance = 1e-5)
})

test_that("dlearner does not fail", {
  expect_no_error(
    dlearner(Y_source = dat_highsim$Y_source,
             Y_target = dat_highsim$Y_target))
})

test_that("dlearner output is correct", {
  expected_result <- t(matrix(c(0.0959171, -1.72143637,  1.15500149,  0.005478273,  0.3258085,
                                0.1077137, -0.32243480, -1.03557177, -0.108195636, -0.4093049,
                                -0.7237275,  0.86634922, -0.43528206,  1.105802194, -0.5961351,
                                1.6347696, -2.91985925,  0.04876288, -2.357504089,  0.8127403,
                                0.4757450,  0.05665543, -1.50373470, -0.744076401, -0.2876033), nrow = 5, ncol = 5))
  result <- dlearner(Y_source = dat_highsim$Y_source,
                    Y_target = dat_highsim$Y_target)
  expect_equal(result$dlearner_estimate[1:5, 1:5], expected_result, tolerance = 1e-5)
})


# Missing data cases
set.seed(1234)
Y_target_highsim_missing <- dat_highsim$Y_target
Y_target_highsim_missing[round(runif(100, 1, 100 * 50))] <- NA

test_that("learner does not fail with missing data", {
  expect_no_error(
    learner(Y_source = dat_highsim$Y_source,
            Y_target = Y_target_highsim_missing,
            lambda_1 = 1, lambda_2 = 1,
            step_size = 0.003))
})

test_that("learner output is correct with missing data", {
  expected_result <- t(matrix(c(0.1402368, -1.53788299,  1.294196,  0.06241446,  0.1034029,
                            0.2234374, -0.16286860, -1.017956, -0.13353267, -0.3327750,
                            -0.4183411,  1.20040593, -1.115950,  0.59002710, -0.4299424,
                            1.4181535, -3.09632714,  1.011979, -1.75238013,  0.5554817,
                            0.5603679, -0.06539629, -1.381598, -0.75764518, -0.1373639), nrow = 5, ncol = 5))
  result <- learner(Y_source = dat_highsim$Y_source,
                    Y_target = Y_target_highsim_missing,
                    lambda_1 = 1, lambda_2 = 1,
                    step_size = 0.003)
  expect_equal(result$learner_estimate[1:5, 1:5], expected_result, tolerance = 1e-5)
})

test_that("cv.learner does not fail with missing data", {
  expect_no_error(
    cv.learner(Y_source = dat_highsim$Y_source,
               Y_target = Y_target_highsim_missing,
               lambda_1_all = c(1, 10),
               lambda_2_all = c(1, 10),
               step_size = 0.003))
})

test_that("cv.learner output is correct with missing data", {
  expected_result <- t(matrix(c(5457.513, 5254.740,
                                5291.389, 5175.707,
                                4933.222, 4947.356), nrow = 2, ncol = 3))
  set.seed(1234)
  result <- cv.learner(Y_source = dat_highsim$Y_source,
                       Y_target = Y_target_highsim_missing,
                       lambda_1_all = c(1, 10, 100),
                       lambda_2_all = c(1, 10),
                       step_size = 0.003)
  expect_equal(result$mse, expected_result, tolerance = 1e-3)
})
