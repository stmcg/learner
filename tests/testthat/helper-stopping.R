# Analytic examples distinguish a stopping code from optimization accuracy.
# The source has distinct singular values 4, 2, 1; its rank-one approximation
# and initial target estimate are exactly diag(4, 0, 0) in a 4-by-3 matrix.
stopping_example <- function() {
  source <- matrix(0, 4, 3)
  source[cbind(1:3, 1:3)] <- c(4, 2, 1)
  target <- matrix(0, 4, 3)
  target[1, 1] <- 4
  list(source = source, target = target)
}
