# These expected values come from Sean's base-R learnerv2 0.3.0 implementation,
# not from the function under test. See fixtures/README.md for provenance.
reference_fixture <- function() {
  readRDS(test_path('fixtures', 'multisource-v030.rds'))
}

with_reference_rng <- function(seed, kind, code) {
  old_kind <- RNGkind()
  had_seed <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get('.Random.seed', envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign('.Random.seed', old_seed, envir = .GlobalEnv)
    } else if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)) {
      rm('.Random.seed', envir = .GlobalEnv)
    }
  })
  do.call(RNGkind, as.list(kind))
  set.seed(seed)
  force(code)
}

for (case_name in c('complete_common', 'complete_separate',
                    'missing_common', 'missing_separate')) {
  test_that(paste('multisource fit matches the saved reference:', case_name), {
    fixture <- reference_fixture()
    case <- fixture$fit[[case_name]]
    actual <- do.call(learner, case$args)
    # Includes every estimated entry, every objective value, rank, and stop code.
    expect_equal(actual, case$expected, tolerance = fixture$provenance$tolerance)
  })
}

for (case_name in c('complete', 'missing')) {
  test_that(paste('separate-penalty CV matches the saved reference:', case_name), {
    fixture <- reference_fixture()
    case <- fixture$cv[[case_name]]
    actual <- with_reference_rng(case$seed, fixture$provenance$rng_kind,
                                  do.call(cv.learner, case$args))
    # Includes all eight scores and the selected penalties, not only dimensions.
    expect_equal(actual, case$expected, tolerance = fixture$provenance$tolerance)
  })
}

test_that('multisource D-LEARNER matches the saved reference', {
  fixture <- reference_fixture()
  actual <- do.call(dlearner, fixture$direct$args)
  expect_equal(actual, fixture$direct$expected,
               tolerance = fixture$provenance$tolerance)
})
