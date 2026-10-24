# Population spaces and contributions using simulated multiple-source data.
space_target <- learner::dat_multisource$Y_target
space_sources <- learner::dat_multisource$Y_source
learner::plot_source_spaces(space_target, space_sources, r=3,
                            target_name='Target (simulated)')
space_fit <- learner::dlearner(space_sources,space_target,r=3,source_weights=c(0.7,0.3))
learner::plot_contributions(space_fit, main='D-LEARNER contributions (simulated data)')
