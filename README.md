
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/stmcg/learner/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stmcg/learner/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/stmcg/learner/graph/badge.svg)](https://app.codecov.io/gh/stmcg/learner)
<!-- badges: end -->

# learner

The `learner` package implements transfer learning methods for low-rank
matrix estimation. These methods leverage similarity in the latent row
and column spaces between the source and target populations to improve
estimation in the target population. The methods include the LatEnt
spAce-based tRaNsfer lEaRning (LEARNER) method and the direct projection
LEARNER (D-LEARNER) method described by [McGrath et
al. (2024)](https://doi.org/10.48550/arXiv.2412.20605).

## Installation

You can install the development version of `learner` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("stmcg/learner")
```

## Example

We illustrate an example of how `learner` can be used. We first load the
package.

``` r
library(learner)
```

In this illustration, we will use one of the toy data sets in the
package (`dat_highsim`) that has a high degree of similarity between the
latent spaces of the source and target populations. The object
`dat_highsim` is a list which contains the observed source population
data matrix `Y_source` and the target population data matrix `Y_target`.
Since the data was simulated, the true values of the matrices are
included in `dat_highsim` as `Theta_source` and `Theta_target`.

#### LEARNER Method

We can apply the LEARNER method via the `learner` function.

This method allows for flexible patterns of heterogeneity across the
source and target populations. It consequently requires specifying the
tuning parameters `lambda_1` and `lambda_2` which control the degree of
transfer learning between the populations. These values can be selected
based on cross-validation via the `cv.learner` function. For example, we
can specify candidate values of 1, 10, and 100 for `lambda_1` and
`lambda_2` and select the optimal values based on cross-validation as
follows:

``` r
res_cv.learner <- cv.learner(Y_source = dat_highsim$Y_source, 
                             Y_target = dat_highsim$Y_target, 
                             lambda_1_all = c(1, 10, 100), 
                             lambda_2_all = c(1, 10, 100), 
                             step_size = 0.003)
res_cv.learner$lambda_1_min
#> [1] 100
res_cv.learner$lambda_2_min
#> [1] 1
```

Next, we apply the `learner` function with these values of `lambda_1`
and `lambda_2`:

``` r
res_learner <- learner(Y_source = dat_highsim$Y_source, 
                       Y_target = dat_highsim$Y_target,
                       lambda_1 = 100, lambda_2 = 1, 
                       step_size = 0.003)
```

The LEARNER estimate is given by the `learner_estimate` component in the
output of the `learner` function, e.g.,

``` r
res_learner$learner_estimate[1:5, 1:5]
#>            [,1]        [,2]       [,3]         [,4]       [,5]
#> [1,]  0.1468889 -1.63414688  1.1857308  0.008093962  0.2196070
#> [2,]  0.1087766 -0.24392770 -1.0061737 -0.058119721 -0.4012974
#> [3,] -0.6213028  1.04206753 -0.6907776  0.905709331 -0.5324437
#> [4,]  1.5811588 -3.00001490  0.3978080 -2.124285120  0.6909386
#> [5,]  0.5046133  0.01845138 -1.4798209 -0.765601712 -0.2230227
```

We can check the convergence of the numerical optimization algorithm by
plotting the values of the objective function at each iteration:

``` r
plot(res_learner$objective_values, 
     xlab = 'Iteration', ylab = 'Objective function',
     type = 'l', col = 'blue', lwd = 2)
```

<img src="README_files/figure-gfm/unnamed-chunk-6-1.png" width="75%" />

#### D-LEARNER Method

We can apply the D-LEARNER method via the `dlearner` function. This
method makes stronger assumptions on the heterogeneity across the source
and target populations. It consequently does not rely on choosing tuning
parameters. The `dlearner` function can be applied as follows:

``` r
res_dlearner <- dlearner(Y_source = dat_highsim$Y_source, 
                         Y_target = dat_highsim$Y_target)
res_dlearner$dlearner_estimate[1:5, 1:5]
#>            [,1]        [,2]        [,3]         [,4]       [,5]
#> [1,]  0.0959171 -1.72143637  1.15500149  0.005478273  0.3258085
#> [2,]  0.1077137 -0.32243480 -1.03557177 -0.108195636 -0.4093049
#> [3,] -0.7237275  0.86634922 -0.43528206  1.105802194 -0.5961351
#> [4,]  1.6347696 -2.91985925  0.04876288 -2.357504089  0.8127403
#> [5,]  0.4757450  0.05665543 -1.50373470 -0.744076401 -0.2876033
```
