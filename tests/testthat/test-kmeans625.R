# Input checking (invalid inputs should give correct errors)

test_that("x must be a numeric matrix or data frame", {
  expect_error(kmeans625(letters, 2), "x must be a numeric matrix")
  expect_error(kmeans625(1:10, 2), "x must be a numeric matrix")
  # A data frame with a non-numeric column becomes a character matrix
  expect_error(kmeans625(iris, 3), "x must be a numeric matrix")
})

test_that("x must not contain missing values", {
  x <- matrix(c(1, 2, NA, 4), ncol = 2)
  expect_error(kmeans625(x, 1), "missing values")
})

test_that("a matrix of centers must have the same number of columns as x", {
  x <- as.matrix(iris[, 1:4])
  init <- matrix(0, nrow = 3, ncol = 2)
  expect_error(kmeans625(x, init), "same number of columns")
})

test_that("a number of centers must be a whole number between 1 and nrow(x)", {
  x <- as.matrix(iris[, 1:4])
  expect_error(kmeans625(x, 0), "centers must be a whole number")
  expect_error(kmeans625(x, 2.5), "centers must be a whole number")
  expect_error(kmeans625(x, nrow(x) + 1), "centers must be a whole number")
  # A plain vector is neither a single number nor a matrix
  expect_error(kmeans625(x, c(1, 2, 3)), "centers must be a whole number")
})

test_that("iter.max must be a whole number that is at least 1", {
  x <- as.matrix(iris[, 1:4])
  expect_error(kmeans625(x, 3, iter.max = 0), "iter.max must be")
  expect_error(kmeans625(x, 3, iter.max = 1.5), "iter.max must be")
  expect_error(kmeans625(x, 3, iter.max = c(5, 10)), "iter.max must be")
})

test_that("nstart must be a whole number that is at least 1", {
  x <- as.matrix(iris[, 1:4])
  expect_error(kmeans625(x, 3, nstart = 0), "nstart must be")
  expect_error(kmeans625(x, 3, nstart = 2.5), "nstart must be")
  expect_error(kmeans625(x, 3, nstart = c(1, 2)), "nstart must be")
})

# Error and warning while algorithm runs

test_that("an empty cluster gives a warning and matches stats::kmeans", {
  x <- matrix(c(0, 0, 10, 10, 0, 1, 10, 11), ncol = 2)
  # The third center is far from every point, so its cluster stays empty
  init <- rbind(c(0, 0), c(10, 10), c(1000, 1000))
  expect_warning(ours <- kmeans625(x, init), "empty cluster")
  ref <- suppressWarnings(stats::kmeans(x, init, algorithm = "Lloyd"))
  expect_equal(ours$cluster, ref$cluster, ignore_attr = TRUE)
  expect_equal(ours$centers, ref$centers, ignore_attr = TRUE)
  expect_equal(ours$withinss, ref$withinss)
  expect_equal(ours$size, ref$size)
  expect_equal(ours$iter, ref$iter)
})

test_that("reaching iter.max without converging gives a warning", {
  x <- as.matrix(iris[, 1:4])
  init <- x[c(1, 51, 101), ]
  expect_warning(kmeans625(x, init, iter.max = 1), "did not converge")
})

# Agreement with stats::kmeans(algorithm = "Lloyd")
# Both functions get the same starting centers (or the same seed with one
# start)

test_that("matches stats::kmeans on iris from the same starting centers", {
  x <- as.matrix(iris[, 1:4])
  init <- x[c(1, 51, 101), ]
  ours <- kmeans625(x, init, iter.max = 50)
  ref <- stats::kmeans(x, init, iter.max = 50, algorithm = "Lloyd")
  expect_equal(ours$cluster, ref$cluster, ignore_attr = TRUE)
  expect_equal(ours$centers, ref$centers, ignore_attr = TRUE)
  expect_equal(ours$tot.withinss, ref$tot.withinss)
})

test_that("every returned field matches stats::kmeans on simulated blobs", {
  x <- sim_blobs()
  init <- x[c(1, 51, 101), ]
  ours <- kmeans625(x, init, iter.max = 50)
  ref <- stats::kmeans(x, init, iter.max = 50, algorithm = "Lloyd")
  expect_equal(ours$cluster, ref$cluster, ignore_attr = TRUE)
  expect_equal(ours$centers, ref$centers, ignore_attr = TRUE)
  expect_equal(ours$totss, ref$totss)
  expect_equal(ours$withinss, ref$withinss)
  expect_equal(ours$tot.withinss, ref$tot.withinss)
  expect_equal(ours$betweenss, ref$betweenss)
  expect_equal(ours$size, ref$size)
  expect_equal(ours$iter, ref$iter)
})

test_that("a number of centers matches stats::kmeans after the same seed", {
  x <- sim_blobs()
  set.seed(1)
  ours <- kmeans625(x, centers = 3, iter.max = 50)
  set.seed(1)
  ref <- stats::kmeans(x, centers = 3, iter.max = 50, algorithm = "Lloyd")
  expect_equal(ours$cluster, ref$cluster, ignore_attr = TRUE)
  expect_equal(ours$tot.withinss, ref$tot.withinss)
})

test_that("a single-column x is handled", {
  x <- as.matrix(iris[, "Petal.Length", drop = FALSE])
  init <- x[c(1, 51), , drop = FALSE]
  ours <- kmeans625(x, init, iter.max = 50)
  ref <- stats::kmeans(x, init, iter.max = 50, algorithm = "Lloyd")
  expect_equal(dim(ours$centers), c(2, 1))
  expect_equal(ours$cluster, ref$cluster, ignore_attr = TRUE)
  expect_equal(ours$tot.withinss, ref$tot.withinss)
})

# Other valid inputs

test_that("a data frame gives the same result as the equivalent matrix", {
  x <- as.matrix(iris[, 1:4])
  init <- x[c(1, 51, 101), ]
  from_matrix <- kmeans625(x, init, iter.max = 50)
  from_df <- kmeans625(iris[, 1:4], init, iter.max = 50)
  expect_equal(from_df, from_matrix)
})

test_that("one cluster puts every point at the overall mean", {
  x <- sim_blobs()
  fit <- kmeans625(x, centers = 1)
  expect_true(all(fit$cluster == 1))
  expect_equal(as.vector(fit$centers), colMeans(x))
  expect_equal(fit$tot.withinss, fit$totss)
  expect_equal(fit$size, nrow(x))
})

# Random starts (nstart >= 2)
# Not compared to stats::kmeans: with nstart >= 2 it draws its random starts
# differently, so the same seed gives different starting centers. These tests
# check properties of the result instead

test_that("more random starts never give a worse fit", {
  x <- sim_blobs()
  set.seed(2)
  one_start <- kmeans625(x, centers = 3, iter.max = 50, nstart = 1)
  ten_starts <- kmeans625(x, centers = 3, iter.max = 50, nstart = 10)
  expect_lte(ten_starts$tot.withinss, one_start$tot.withinss)
})

test_that("several random starts recover well-separated blobs", {
  x <- sim_blobs(n_each = 50)
  set.seed(2)
  fit <- kmeans625(x, centers = 3, iter.max = 50, nstart = 10)
  expect_equal(sort(fit$size), c(50, 50, 50))
  # Points simulated from the same blob share one label
  expect_length(unique(fit$cluster[1:50]), 1)
  expect_length(unique(fit$cluster[51:100]), 1)
  expect_length(unique(fit$cluster[101:150]), 1)
})

test_that("nstart is ignored when centers is a matrix", {
  x <- sim_blobs()
  init <- x[c(1, 51, 101), ]
  one_start <- kmeans625(x, init, iter.max = 50, nstart = 1)
  five_starts <- kmeans625(x, init, iter.max = 50, nstart = 5)
  expect_equal(five_starts, one_start)
})

# Test internal helpers directly 

test_that("assign_clusters picks the nearest center", {
  x <- matrix(c(0, 0, 10, 10, 0, 1, 10, 11), ncol = 2)
  centers <- x[c(1, 3), ]
  expect_equal(assign_clusters(x, centers), c(1, 1, 2, 2))
  # With one center every point gets label 1
  expect_equal(assign_clusters(x, x[1, , drop = FALSE]), c(1, 1, 1, 1))
})

test_that("assign_clusters sends ties to the lowest-numbered cluster", {
  x <- matrix(c(5, 5), ncol = 2)
  centers <- rbind(c(0, 5), c(10, 5))
  expect_equal(assign_clusters(x, centers), 1)
})

test_that("update_centers returns the mean of each cluster", {
  x <- matrix(c(0, 0, 10, 10, 0, 1, 10, 11), ncol = 2)
  centers <- update_centers(x, cluster = c(1, 1, 2, 2), k = 2)
  expect_equal(centers, rbind(c(0, 0.5), c(10, 10.5)), ignore_attr = TRUE)
})
