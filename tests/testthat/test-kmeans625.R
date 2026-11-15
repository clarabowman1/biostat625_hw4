test_that("matches stats::kmeans from the same starting centers", {
  x <- as.matrix(iris[, 1:4])
  init <- x[c(1, 51, 101), ]
  ours <- kmeans625(x, init, iter.max = 50)
  ref <- stats::kmeans(x, init, iter.max = 50, algorithm = "Lloyd")
  expect_equal(ours$cluster, ref$cluster, ignore_attr = TRUE)
  expect_equal(ours$centers, ref$centers, ignore_attr = TRUE)
  expect_equal(ours$tot.withinss, ref$tot.withinss)
})