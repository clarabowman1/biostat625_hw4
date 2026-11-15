# Simulated data shared by the tests: three well-separated Gaussian blobs in
# two dimensions, with n_each points per blob. Rows 1 to n_each belong to the
# first blob, the next n_each to the second etc 
sim_blobs <- function(n_each = 50) {
  set.seed(625)
  blob_means <- rbind(c(0, 0), c(10, 10), c(0, 10))
  truth <- rep(1:3, each = n_each)
  x <- blob_means[truth, ] + matrix(rnorm(3 * n_each * 2), ncol = 2)
  return(x)
}
