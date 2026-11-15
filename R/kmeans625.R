# ---------------------------------------------------------------------------
# 1. Roxygen block for the main function
# ---------------------------------------------------------------------------
# Replace this outline with a roxygen block (lines starting with #'), in this
# order:
#
#   Title line   e.g. K-means clustering with Lloyd's algorithm
#   Description  what the function does: assign each point to its nearest
#                center, move each center to its cluster mean, repeat until
#                the assignments stop changing
#   @param x         numeric matrix or data frame, one observation per row
#   @param centers   a number of clusters k, or a matrix of starting centers
#                    with one row per cluster
#   @param iter.max  maximum number of iterations
#   @return          a list; name each component (see step e below)
#   @examples        a short run on iris[, 1:4] with set.seed() and 3 clusters
#   @export
#
# Then run devtools::document() to generate the help page.

# ---------------------------------------------------------------------------
# 2. Main function
# ---------------------------------------------------------------------------
kmeans625 <- function(x, centers, iter.max = 10) {
  # a. Check inputs
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix or a data frame of numeric columns")
  }
  if (anyNA(x)) {
    stop("x must not contain missing values")
  }
  n <- nrow(x)

  # b. Set starting centers
  if (is.matrix(centers)) {
    # A matrix of starting centers was supplied: use it as given
    if (ncol(centers) != ncol(x)) {
      stop("centers must have the same number of columns as x")
    }
    k <- nrow(centers)
  } else {
    # A number of clusters was supplied: pick k distinct rows of x at random
    k <- centers
    if (length(k) != 1 || k < 1 || k > n || k != round(k)) {
      stop("centers must be a whole number between 1 and nrow(x), or a matrix")
    }
    centers <- x[sample.int(n, k), , drop = FALSE]
  }

  # c. Iterate until the assignments stop changing, at most iter.max rounds
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

  # d. Summarise
  # Squared distance from each point to the center of its own cluster
  sq_dist <- rowSums((x - centers[cluster, , drop = FALSE])^2)
  withinss <- as.vector(rowsum(sq_dist, cluster))
  tot.withinss <- sum(withinss)
  # Squared distance from each point to the overall mean
  totss <- sum((x - rep(colMeans(x), each = n))^2)

  # e. Return a list with the same names stats::kmeans() uses
  return(list(
    cluster = cluster,
    centers = centers,
    totss = totss,
    withinss = withinss,
    tot.withinss = tot.withinss,
    betweenss = totss - tot.withinss,
    size = tabulate(cluster, nbins = k),
    iter = iter
  ))
}

# ---------------------------------------------------------------------------
# 3. Internal helpers (not exported)
# ---------------------------------------------------------------------------

# Assign each point to its nearest center.
#   x        numeric matrix of observations (n rows)
#   centers  numeric matrix of current centers (k rows)
# Returns an integer vector of length n with values in 1..k. Ties go to the
# lowest-numbered cluster so that results are reproducible.
assign_clusters <- function(x, centers) {
  k <- nrow(centers)
  # Calculate distances to centers and assign to closest ones
  # NOTE: may rewrite to speed up/avoid using stats::dist
  D <- as.matrix(stats::dist(rbind(centers, x)))[-(1:k), 1:k, drop = FALSE]  
  max.col(-D, ties.method = 'first')
}

# Move each center to the mean of its assigned points.
#   x        numeric matrix of observations
#   cluster  integer vector of cluster labels, one per row of x
#   k        number of clusters
# Returns a numeric matrix with k rows and ncol(x) columns.
update_centers <- function(x, cluster, k) {
  # Number of points in each cluster
  size <- tabulate(cluster, nbins = k)
  if (any(size == 0)) {
    stop("empty cluster: try a different set of initial centers")
  }
  # Column sums within each cluster (k x p), divided by cluster sizes
  rowsum(x, cluster) / size
}

# ---------------------------------------------------------------------------
# Suggested order to write this file
# ---------------------------------------------------------------------------
# 1. assign_clusters(), tested by hand on a tiny matrix with obvious groups
# 2. update_centers(), tested the same way
# 3. the loop in kmeans625(), starting from a matrix of centers only
# 4. compare with stats::kmeans(x, centers = init, algorithm = "Lloyd") using
#    the same `init`
# 5. add the k-as-a-number option, the input checks and the summaries
# 6. write the roxygen block, then devtools::document() and devtools::check()
#
# Use devtools::load_all() after each edit to try the functions.
# Delete this section and the outline comments once the file is finished.
