if(getRversion() >= "2.15.1"){
  # To remove 'no visible binding for global variable ...' notes
  utils::globalVariables(c('lambda_1_ind', 'lambda_2_ind'))
}

.onAttach <- function(libname, pkgname) {
  max_threads <- tryCatch(omp_max_threads(), error = function(e) 1)
  if (max_threads < 2) {
    packageStartupMessage("*******\n
                          This installation of learner has not detected OpenMP support\n
                          It will still work but will not support multithreading via the `n_cores` argument
                          If you plan to use multithreading, please ensure you are using R>=3.4.0 and have OpenMP installed\n
                          *******")
  }
}
