test_that('CV annotations distinguish small differences from exact ties and numeric boundaries', {
  testthat::skip_if_not_installed("ggplot2", minimum_version = "3.4.0")
  path<-plot_device();on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  x<-list(mse_all=matrix(c(1.06,1.051,1.07,1.06,1.052,1.07),3),
          lambda_1_min=50,lambda_2_min=0.01)
  d<-plot_cv(x,lambda_1_all=c(10,50,100),lambda_2_all=c(0.01,1))
  expect_equal(d$exact_minima,1)
  expect_equal(d$score_range,c(1.051,1.07))
  expect_equal(d$lambda_2_slice_range,c(1.051,1.052))
  expect_identical(d$boundary,c(common=FALSE,balance=TRUE))
  x$mse_all[2,2]<-1.051
  expect_equal(plot_cv(x,lambda_1_all=c(10,50,100),lambda_2_all=c(0.01,1))$exact_minima,2)
  x$mse_all<-matrix(c(3,1,2),3);x$lambda_1_min<-1;x$lambda_2_min<-1
  expect_identical(plot_cv(x,lambda_1_all=c(3,1,2),lambda_2_all=1)$boundary,
                   c(common=TRUE,balance=FALSE))
})

test_that('CV profiles sort horizontal values without reordering stored scores', {
  path<-plot_device();on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  old<-graphics::par(c('mfrow','mar'))
  x<-list(mse_all=matrix(6:1,3),lambda_1_min=10,lambda_2_min=1)
  d<-plot_cv_profile(x,lambda_1_all=c(100,1,10),lambda_2_all=c(.1,1),log_labels=TRUE)
  expect_equal(d$x_order,c(2L,3L,1L));expect_equal(d$scores,x$mse_all)
  expect_equal(d$best_index,c(3L,2L))
  expect_equal(graphics::par(c('mfrow','mar')),old)
  expect_error(plot_cv_profile(x,lambda_1_all=c(10,1,10),lambda_2_all=c(.1,1)), 'distinct')
  d<-plot_cv_profile(cv_plot_example(),lambda_1_row_all=c(10,0),
                    lambda_1_col_all=c(.2,3),lambda_2_all=c(0,1),lambda_2_index=2)
  expect_equal(d$shown_lambda_2_index,2)
  expect_equal(d$best_index,c(2L,2L,2L))
})

test_that('space plots compute full-matrix spaces before taking common submatrices', {
  path<-plot_device();on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  # This rank-one matrix has analytically known normalized left/right vectors.
  u<-c(1,2,3);v<-c(2,1);y<-outer(u,v)
  set.seed(1);before<-.Random.seed
  old<-graphics::par(c('mfrow','mar'))
  d<-plot_source_spaces(y,list(A=y,B=2*y),r=1,row_index=c(3,1),column_index=2:1)
  expect_equal(d$row_projectors[[1]],outer(u,u)[c(3,1),c(3,1)]/sum(u^2))
  expect_equal(d$column_projectors[[1]],outer(v,v)[2:1,2:1]/sum(v^2))
  expect_equal(d$row_projectors[[2]],d$row_projectors[[3]])
  expect_identical(.Random.seed,before)
  expect_equal(graphics::par(c('mfrow','mar')),old)
  expect_equal(d$limits,c(row=9/14,column=4/5))
  clipped<-plot_source_spaces(y,y,r=1,limits=c(.1,.1))
  expect_gt(sum(unlist(clipped$clipped_count)),0)
  expect_equal(clipped$row_projectors[[1]],outer(u,u)/sum(u^2))
  expect_error(plot_source_spaces(y,y,r=2),'zero singular')
  expect_error(plot_source_spaces(y,y,r=1,row_index=c(1,1)),'distinct')
  expect_error(plot_source_spaces(y,y,r=1,limits=c(0,.1)),'positive')
  y[1]<-NA
  expect_error(plot_source_spaces(y,y,r=1),'complete finite')
})

test_that('contributions use squared normalized singular vectors, not raw factors', {
  path<-plot_device();on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  u<-c(1,2,3);v<-c(2,1);y<-outer(u,v)
  d<-plot_contributions(y,r=1,row_index=c(3,1),column_index=2)
  expect_equal(drop(d$row_scores),u^2/sum(u^2))
  expect_equal(drop(d$column_scores),v^2/sum(v^2))
  expect_equal(colSums(d$row_scores),1)
  expect_equal(colSums(d$column_scores),1)
  expect_equal(d$row_index,c(3L,1L))
  expect_equal(plot_contributions(-y,r=1)$row_scores,d$row_scores)
  expect_equal(plot_contributions(list(learner_estimate=y,r=1))$column_scores,d$column_scores)
  expect_equal(plot_contributions(list(dlearner_estimate=y,r=1))$row_scores,d$row_scores)
  expect_error(plot_contributions(y),'Supply r')
  expect_error(plot_contributions(y,r=2),'zero singular')
  expect_error(plot_contributions(y,r=1,column_index=3),'valid integer')
})

test_that('simulation plots summarize replicate errors and preserve graphics state', {
  path<-plot_device();on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  d<-data.frame(scenario='High',rank=2,variance_ratio=c(1,1,3),method='LEARNER',error=c(2,4,1))
  old<-graphics::par(c('mfrow','mar'));original<-d
  for(uncertainty in c('none','sd','se')) {
    out<-plot_simulation_results(d,uncertainty=uncertainty)
    expect_equal(out$mean,c(3,1));expect_equal(out$n,c(2L,1L))
    expect_equal(out$sd,c(sqrt(2),NA));expect_equal(out$se,c(1,NA))
    expect_equal(graphics::par(c('mfrow','mar')),old)
  }
  expect_identical(d,original)
  d$error[1]<-NA;expect_error(plot_simulation_results(d),'finite')
  d$error[1]<- -1;expect_error(plot_simulation_results(d),'nonnegative')
  d$error[1]<-1;d$variance_ratio[1]<-0;expect_error(plot_simulation_results(d),'positive')
  expect_error(plot_simulation_results(data.frame()),'nonempty')
  # Empty scenario/rank combinations and pages must not synthesize zero errors.
  d<-data.frame(scenario=rep(c('A','B','C','D'),each=3),rank=rep(1:3,4),
                variance_ratio=1,method='Other',error=1:12)
  expect_equal(nrow(plot_simulation_results(d)),12)
})
