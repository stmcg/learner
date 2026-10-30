test_that('relative stopping requires K consecutive changes, including at zero', {
  source <- matrix(0, 4, 3); source[cbind(1:3,1:3)] <- c(4,2,1)
  target <- matrix(0,4,3); target[1,1] <- 4
  fit <- function(...) learner(source,target,r=1,lambda_1=1,lambda_2=1,
    step_size=.01,control=c(list(stopping_rule='relative',max_iter=10),list(...)))
  expect_length(fit()$objective_values,6)
  expect_equal(fit()$convergence_criterion,1)
  expect_length(fit(window=1)$objective_values,2)
  expect_length(fit(threshold=0)$objective_values,10)
  expect_equal(fit(threshold=0)$convergence_criterion,2)
  expect_equal(fit(window=20)$convergence_criterion,2)
})

test_that('stopping matches rolling maxima of recorded changes', {
  source <- dat_highsim$Y_source[1:10,1:8]
  target <- dat_highsim$Y_target[1:10,1:8]
  run <- function(control) learner(source,target,r=2,lambda_1=1,lambda_2=.1,
                                  step_size=.01,control=control)
  full <- run(list(max_iter=80,threshold=0,max_value=1e10))$objective_values
  for (rule in c('absolute','relative')) {
    changes <- abs(diff(full))
    if (rule=='relative') changes <- changes/abs(head(full,-1))
    for (k in c(1,3,5)) {
      threshold <- as.numeric(quantile(changes,.6))
      ends <- seq.int(k,length(changes))
      pass <- vapply(ends,function(i) max(changes[(i-k+1):i]) < threshold,logical(1))
      expected <- if (any(pass)) ends[which(pass)[1]]+1L else length(full)
      fit <- run(list(max_iter=80,threshold=threshold,stopping_rule=rule,
                      window=k,max_value=1e10))
      expect_length(fit$objective_values,expected)
      expect_equal(fit$objective_values,head(full,expected))
      expect_equal(fit$convergence_criterion,if(any(pass)) 1 else 2)
    }
  }
})

test_that('relative settings propagate into each CV fold', {
  source <- dat_highsim$Y_source[1:8,1:6]
  target <- dat_highsim$Y_target[1:8,1:6]
  control <- list(stopping_rule='relative',window=3,threshold=.1,max_iter=20)
  set.seed(91)
  indices <- sample(seq_along(target))
  folds <- split(indices,rep(1:2,each=length(indices)/2))
  expected <- sum(vapply(folds,function(idx) {
    train <- target;train[idx] <- NA
    fit <- learner(source,train,r=2,lambda_1=1,lambda_2=.1,
                   step_size=.001,control=control)
    sum((target[idx]-fit$learner_estimate[idx])^2)
  },numeric(1)))
  set.seed(91)
  cv <- cv.learner(source,target,r=2,lambda_1_all=1,lambda_2_all=.1,
                   n_folds=2,step_size=.001,control=control)
  expect_equal(as.numeric(cv$mse_all),expected)
})

test_that('invalid stopping options fail early', {
  for (rule in list('other',NA_character_,c('relative','absolute'),1))
    expect_error(prepare_control(list(stopping_rule=rule),.01),'stopping_rule')
  for (window in list(0,-1,1.5,Inf,NA,TRUE))
    expect_error(prepare_control(list(window=window),.01),'window')
  expect_equal(prepare_control(list(stopping_rule='relative'),.01)$threshold,1e-6)
})

test_that('relative stopping retains failure protection even with a loose tolerance', {
  d <- stopping_example()
  fit <- learner(d$source,-d$target,r=1,lambda_1=1,lambda_2=1,step_size=10,
    control=list(stopping_rule='relative',window=1,threshold=1e100,
                 max_value=2,max_iter=10))
  expect_equal(fit$convergence_criterion,3)
})

test_that('k and epsilon are other names for window and threshold', {
  pc <- prepare_control(list(k=3,epsilon=1e-5),.01)
  expect_equal(pc$window,3)
  expect_equal(pc$threshold,1e-5)
  expect_null(pc$k)
  expect_null(pc$epsilon)
  # Defaults still depend on the rule when only one of the two is given.
  expect_equal(prepare_control(list(stopping_rule='relative',k=2),.01)$threshold,1e-6)
  expect_equal(prepare_control(list(epsilon=0),.01)$window,1)
  # A NULL alias counts as not given.
  expect_equal(prepare_control(list(k=NULL,window=4),.01)$window,4)
  source <- dat_highsim$Y_source[1:10,1:8]
  target <- dat_highsim$Y_target[1:10,1:8]
  fit <- function(control) learner(source,target,r=2,lambda_1=1,lambda_2=.1,
                                   step_size=.01,control=control)
  for (rule in c('absolute','relative')) {
    by_main <- fit(list(stopping_rule=rule,window=3,threshold=1e-3,max_iter=200))
    by_alias <- fit(list(stopping_rule=rule,k=3,epsilon=1e-3,max_iter=200))
    expect_identical(by_alias,by_main)
  }
  # k really changes when the fit stops: a longer run of small changes is needed.
  updates <- vapply(c(1,3,5),function(k) length(fit(list(
    stopping_rule='relative',k=k,epsilon=1e-3,max_iter=200))$objective_values),numeric(1))
  expect_true(all(diff(updates) > 0))
})

test_that('giving both names for one setting is an error', {
  expect_error(prepare_control(list(window=3,k=3),.01),'not both')
  expect_error(prepare_control(list(threshold=1e-3,epsilon=1e-3),.01),'not both')
  expect_error(prepare_control(list(window=3,k=5),.01),'control\\$window or control\\$k')
  expect_error(prepare_control(list(threshold=1e-3,epsilon=1e-5),.01),
               'control\\$threshold or control\\$epsilon')
})

test_that('errors name the setting the way the user wrote it', {
  for (k in list(0,-1,1.5,Inf,NA,TRUE,c(1,2)))
    expect_error(prepare_control(list(k=k),.01),'control\\$k')
  for (epsilon in list(-1,Inf,NA,'a',c(1e-3,1e-4)))
    expect_error(prepare_control(list(epsilon=epsilon),.01),'control\\$epsilon')
  expect_error(prepare_control(list(window=0),.01),'control\\$window')
  expect_error(prepare_control(list(threshold=-1),.01),'control\\$threshold')
})

test_that('k and epsilon reach every cross-validation fold', {
  source <- dat_highsim$Y_source[1:8,1:6]
  target <- dat_highsim$Y_target[1:8,1:6]
  cv <- function(control) {
    set.seed(91)
    cv.learner(source,target,r=2,lambda_1_all=1,lambda_2_all=.1,
               n_folds=2,step_size=.001,control=control)
  }
  by_main <- cv(list(stopping_rule='relative',window=3,threshold=.1,max_iter=20))
  by_alias <- cv(list(stopping_rule='relative',k=3,epsilon=.1,max_iter=20))
  expect_identical(by_alias$mse_all,by_main$mse_all)
  loose <- cv(list(stopping_rule='relative',k=1,epsilon=1e100,max_iter=20))
  expect_false(isTRUE(all.equal(loose$mse_all,by_main$mse_all)))
})
