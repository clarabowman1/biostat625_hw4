#' K-means clustering with Lloyd's algorithm
#'
#' Groups the rows of `x` into clusters by repeatedly assigning each point to
#' its nearest center and moving each center to the mean of its points, until
#' the assignments stop changing.
#'
#' @details
#' When `centers` is a number k, k distinct rows of `x` are chosen at random
#' as the starting centers, so results depend on the random seed; call
#' [set.seed()] first for reproducible results. Because the final clusters
#' depend on where the algorithm starts, `nstart` can be used to try several
#' random starts and keep the best one. When `centers` is a matrix, the result
#' is fully determined and `nstart` is ignored.
#'
#' Started from the same centers, the result matches
#' `stats::kmeans(x, centers, algorithm = "Lloyd")`. The two functions differ
#' in two cases:
#' * With `nstart` of 2 or more, they draw their random starts differently,
#'   so the same seed gives different starting centers.
#' * If a cluster ends up with no points, `kmeans625()` stops with an error,
#'   while [stats::kmeans()] gives a warning.
#'
#' The algorithm is written in plain R, so it is slower than
#' [stats::kmeans()], whose iteration is compiled C code.
#'
#' @param x A numeric matrix or data frame, one observation per row.
#' @param centers Either the number of clusters or a matrix of starting
#'   centers with one row per cluster.
#' @param iter.max The maximum number of iterations. Default = 10
#' @param nstart Number of random starts to try when `centers` is a number.
#'   The run with the lowest `tot.withinss` is kept. Ignored when `centers`
#'   is a matrix. Default = 1
#'
#' @return A list with the same components as [stats::kmeans()]:
#' \describe{
#'   \item{`cluster`}{Integer vector with one entry per row of `x`, giving
#'     the cluster (from 1 to k) that each point is assigned to.}
#'   \item{`centers`}{Matrix of final cluster centers, with one row per
#'     cluster and one column per column of `x`.}
#'   \item{`totss`}{Total sum of squares: the sum of squared distances from
#'     each point to the overall mean of `x`.}
#'   \item{`withinss`}{Vector of within-cluster sums of squares, one per
#'     cluster: the sum of squared distances from each point to the center
#'     of its own cluster.}
#'   \item{`tot.withinss`}{Total within-cluster sum of squares,
#'     `sum(withinss)`. Smaller values mean tighter clusters.}
#'   \item{`betweenss`}{Between-cluster sum of squares,
#'     `totss - tot.withinss`.}
#'   \item{`size`}{Integer vector of the number of points in each cluster.}
#'   \item{`iter`}{Number of iterations performed.}
#' }
#' When `nstart` is greater than 1, these describe the run with the smallest
#' `tot.withinss`.
#'
#' @seealso [stats::kmeans()], which this function re-implements.
#'
#' @examples
#' x <- as.matrix(iris[, 1:4])
#'
#' # Three clusters from one random start
#' set.seed(625)
#' fit <- kmeans625(x, centers = 3)
#' fit$size
#' fit$centers
#' fit$tot.withinss
#'
#' # A single start can settle on a poor split; keeping the best of ten
#' # random starts gives a much smaller tot.withinss
#' set.seed(625)
#' fit10 <- kmeans625(x, centers = 3, iter.max = 50, nstart = 10)
#' fit10$tot.withinss
#' table(cluster = fit10$cluster, species = iris$Species)
#'
#' # Supply the starting centers (one row from each species)
#' init <- x[c(1, 51, 101), ]
#' fit_init <- kmeans625(x, centers = init)
#' fit_init$size
#'
#' # Same result as stats::kmeans() from the same starting centers
#' ref <- stats::kmeans(x, centers = init, algorithm = "Lloyd")
#' all.equal(fit_init$cluster, ref$cluster, check.attributes = FALSE)
#'
#' @export
kmeans625 <- function(x,
                      centers,
                      iter.max = 10,
                      nstart = 1) {
  # Check inputs
  # Convert x to matrix if data frame
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  # Stop and return error if x is not correct form
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix or a data frame of numeric columns")
  }
  if (anyNA(x)) {
    stop("x must not contain missing values")
  }
  if (length(nstart) != 1 ||
      nstart < 1 || nstart != round(nstart)) {
    stop("nstart must be a whole number that is at least 1")
  }
  if (length(iter.max) != 1 || iter.max < 1 || iter.max != round(iter.max)) {
    stop("iter.max must be a whole number that is at least 1")
  }
  
  # Set starting centers
  n <- nrow(x)
  if (is.matrix(centers)) {
    # Use matrix of starting centers if supplied
    if (ncol(centers) != ncol(x)) {
      stop("centers must have the same number of columns as x")
    }
    k <- nrow(centers)
    random_start <- FALSE
  } else {
    # Number of clusters supplied --> pick k rows of x at random as centers
    k <- centers
    if (length(k) != 1 || k < 1 || k > n || k != round(k)) {
      stop("centers must be a whole number between 1 and nrow(x), or a matrix")
    }
    centers <- x[sample.int(n, k), , drop = FALSE]
    random_start <- TRUE
  }
  
  # Run the algorithm once from the first set of starting centers
  best <- lloyd_run(x, centers, iter.max)
  
  # With random starts, try nstart - 1 more and keep the best run
  if (random_start && nstart >= 2) {
    for (i in 2:nstart) {
      centers <- x[sample.int(n, k), , drop = FALSE]
      fit <- lloyd_run(x, centers, iter.max)
      if (fit$tot.withinss < best$tot.withinss) {
        best <- fit
      }
    }
  }
  return(best)
}

# Internal helper functions -- not exported as part of package

# Run Lloyd's algorithm once from a given set of starting centers.
#   x         numeric matrix of observations
#   centers   numeric matrix of starting centers (k rows)
#   iter.max  maximum number of iterations
# Returns the list of results that kmeans625() returns.
lloyd_run <- function(x, centers, iter.max) {
  n <- nrow(x)
  k <- nrow(centers)
  
  # Iterate until the assignments are constant up to iter.max times
  cluster <- rep(0L, n)  # no point is assigned yet
  converged <- FALSE
  for (iter in 1:iter.max) {
    new_cluster <- assign_clusters(x, centers)
    if (all(new_cluster == cluster)) {
      converged <- TRUE
      break
    }
    cluster <- new_cluster
    centers <- update_centers(x, cluster, k)
  }
  if (!converged) {
    warning("did not converge in ", iter.max, " iterations")
  }
  
  # Squared distance from each point to the center of its own cluster
  sq_dist <- rowSums((x - centers[cluster, , drop = FALSE])^2)
  withinss <- as.vector(rowsum(sq_dist, cluster))
  tot.withinss <- sum(withinss)
  # Squared distance from each point to the overall mean
  totss <- sum((x - rep(colMeans(x), each = n))^2)
  
  # Return a list with the same names stats::kmeans() uses
  return(
    list(
      cluster = cluster,
      centers = centers,
      totss = totss,
      withinss = withinss,
      tot.withinss = tot.withinss,
      betweenss = totss - tot.withinss,
      size = tabulate(cluster, nbins = k),
      iter = iter
    )
  )
}

# Assign each point to its nearest center.
#   x        numeric matrix of observations (n rows)
#   centers  numeric matrix of current centers (k rows)
# Returns an integer vector of length n with values 1-k. Ties go to the
# lowest-numbered cluster so results are reproducible.
assign_clusters <- function(x, centers) {
  k <- nrow(centers)
  # Calculate distances to centers and assign to closest ones
  # NOTE: may rewrite to speed up/avoid using stats::dist
  D <- as.matrix(stats::dist(rbind(centers, x)))[-(1:k), 1:k, drop = FALSE]
  return(max.col(-D, ties.method = "first"))
}

# Move each center to the mean of its assigned points.
#   x        numeric matrix of observations
#   cluster  integer vector of cluster labels for each row of x
#   k        number of clusters
# Returns a numeric matrix with k rows and ncol(x) columns.
update_centers <- function(x, cluster, k) {
  # Number of points in each cluster
  size <- tabulate(cluster, nbins = k)
  if (any(size == 0)) {
    stop("empty cluster: try a different set of initial centers")
  }
  # Column sums within each cluster (k x p), divided by cluster sizes
  return(rowsum(x, cluster) / size)
}
