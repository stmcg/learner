test_that('unchanged consecutive objectives meet a positive threshold', {
  d <- stopping_example()
  fit <- learner(d$source, d$target, r = 1, lambda_1 = 1, lambda_2 = 1,
                  step_size = 0.01, control = list(max_iter = 5, threshold = 0.001))
  expect_equal(fit$convergence_criterion, 1)
  expect_equal(fit$objective_values, c(0, 0), tolerance = 1e-12)
  expect_equal(fit$learner_estimate, d$target, tolerance = 1e-12)
})

test_that('a zero threshold uses a strict inequality and reaches the budget', {
  d <- stopping_example()
  fit <- learner(d$source, d$target, r = 1, lambda_1 = 1, lambda_2 = 1,
                  step_size = 0.01, control = list(max_iter = 5, threshold = 0))
  expect_equal(fit$convergence_criterion, 2)
  expect_length(fit$objective_values, 5)
  expect_equal(fit$learner_estimate, d$target, tolerance = 1e-12)
})

test_that('the returned estimate can retain the initialization after a bad step', {
  d <- stopping_example()
  # Opposite target and an intentionally huge step make the first update worse.
  fit <- learner(d$source, -d$target, r = 1, lambda_1 = 1, lambda_2 = 1,
                  step_size = 10, control = list(max_iter = 1, threshold = 0))
  initial_objective <- sum((d$target - (-d$target))^2)
  expect_equal(fit$convergence_criterion, 2)
  expect_length(fit$objective_values, 1)
  expect_gt(fit$objective_values[1], initial_objective)
  expect_equal(fit$learner_estimate, d$target, tolerance = 1e-12)
})

test_that('excessive objective growth returns code three and the best estimate', {
  d <- stopping_example()
  fit <- learner(d$source, -d$target, r = 1, lambda_1 = 1, lambda_2 = 1,
                  step_size = 10,
                  control = list(max_iter = 5, threshold = 0, max_value = 2))
  expect_equal(fit$convergence_criterion, 3)
  expect_length(fit$objective_values, 2)
  expect_gt(fit$objective_values[2], 2 * fit$objective_values[1])
  expect_equal(fit$learner_estimate, d$target, tolerance = 1e-12)
})
