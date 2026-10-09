# Paper-inspired displays using simulated data, not published study results.
# Run after loading the plotting development version of learner.
set.seed(803)
space_source_A <- learner::dat_highsim$Y_source
space_source_B <- space_source_A + matrix(rnorm(length(space_source_A), sd=0.3),
                                         nrow(space_source_A), ncol(space_source_A))
space_target <- learner::dat_highsim$Y_target
space_sources <- list(A=space_source_A, B=space_source_B)
learner::plot_source_spaces(space_target, space_sources, r=3,
                            limits=c(0.08,0.15), target_name='Target (simulated)')
space_fit <- learner::dlearner(space_sources,space_target,r=3,source_weights=c(0.7,0.3))
learner::plot_contributions(space_fit, main='D-LEARNER contributions (simulated data)')

# A small, reproducible simulation demonstrates the comparison-plot interface.
# Fixed penalties and only three replicates: this is not a performance study.
# The target matrix is reused across noise ratios within each replicate.
# The error is the Frobenius norm, not MSE or RMSE.
run_toy_simulation <- function() {
  set.seed(910)

  p <- 24L; q <- 16L
  ratios <- c(1,3,5,10)
  results <- list(); index <- 0L
  for (similarity in c('High similarity','Moderate similarity','Low similarity')) {
    alignment <- switch(similarity,'High similarity'=0.98,'Moderate similarity'=0.75,
                         'Low similarity'=0.35)
    for (rank in c(2L,3L)) for (replicate in 1:3) {
      u_all <- qr.Q(qr(matrix(rnorm(p*2*rank),p,2*rank)))
      v_all <- qr.Q(qr(matrix(rnorm(q*2*rank),q,2*rank)))
      u <- u_all[,seq_len(rank),drop=FALSE];v <- v_all[,seq_len(rank),drop=FALSE]
      source_u <- alignment*u+sqrt(1-alignment^2)*u_all[,rank+seq_len(rank),drop=FALSE]
      source_v <- alignment*v+sqrt(1-alignment^2)*v_all[,rank+seq_len(rank),drop=FALSE]
      truth <- u %*% diag(seq(12,8,length.out=rank)) %*% t(v)
      source_truth <- source_u %*% diag(seq(12,8,length.out=rank)) %*% t(source_v)
      target <- truth+matrix(rnorm(p*q),p,q)
      source_noise <- matrix(rnorm(p*q),p,q)
      s <- svd(target,nu=rank,nv=rank)
      baseline <- s$u %*% diag(s$d[seq_len(rank)]) %*% t(s$v)
      for (ratio in ratios) {
        source <- source_truth+source_noise/sqrt(ratio)
        direct <- learner::dlearner(source,target,r=rank)
        learned <- learner::learner(source,target,r=rank,lambda_1=2,lambda_2=0.1,
                                    step_size=0.001,control=list(max_iter=10000))
        estimates <- list('Target-only SVD'=baseline,'D-LEARNER'=direct$dlearner_estimate,
                           'LEARNER'=learned$learner_estimate)
        for (method in names(estimates)) {
          index <- index+1L
          results[[index]] <- data.frame(scenario=similarity,rank=rank,variance_ratio=ratio,
                  method=method,replicate=replicate,error=norm(estimates[[method]]-truth,'F'),
                  stop_code=if(method=='LEARNER') learned$convergence_criterion else NA_integer_)
        }
      }
    }
  }
  simulation_results <- do.call(rbind,results)
  simulation_results
}
simulation_results <- run_toy_simulation()
learner::plot_simulation_results(simulation_results,
          main='Toy simulation: fixed penalties, 3 replicates',ylab='Frobenius estimation error')
# Report optimizer stopping codes before interpreting method comparisons.
print(table(simulation_results$stop_code,useNA='no'))
# Error bars are optional; these are SEs, not confidence intervals.
# learner::plot_simulation_results(simulation_results, uncertainty='se')
