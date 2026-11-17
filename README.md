# kmeans625

<!-- badges: start -->
[![R-CMD-check](https://github.com/clarabowman1/biostat625_hw4/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/clarabowman1/biostat625_hw4/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/clarabowman1/biostat625_hw4/graph/badge.svg)](https://app.codecov.io/gh/clarabowman1/biostat625_hw4)
<!-- badges: end -->

`kmeans625` is an R package that implements k-means clustering from
scratch, using the standard algorithm (Lloyd's algorithm). It was written
for Biostat 625 Homework 4 as a re-implementation of `stats::kmeans()`.

The algorithm repeats two steps until the cluster assignments stop
changing:

1. assign each point to its nearest center;
2. move each center to the mean of the points assigned to it.

## Installation

```r
# install.packages("remotes")
remotes::install_github("clarabowman1/biostat625_hw4", build_vignettes = TRUE)
```

## Example

Cluster the four `iris` measurements into three groups, keeping the best
of ten random starts:

```r
library(kmeans625)

x <- as.matrix(iris[, 1:4])
set.seed(625)
fit <- kmeans625(x, centers = 3, iter.max = 50, nstart = 10)

fit$size
#> [1] 38 62 50

table(cluster = fit$cluster, species = iris$Species)
#>        species
#> cluster setosa versicolor virginica
#>       1      0          2        36
#>       2      0         48        14
#>       3     50          0         0
```

The result is a list with the same components as `stats::kmeans()`:
`cluster`, `centers`, `totss`, `withinss`, `tot.withinss`, `betweenss`,
`size` and `iter`.

## Comparison with `stats::kmeans()`

Given the same starting centers, `kmeans625()` returns the same result as
`stats::kmeans(algorithm = "Lloyd")`. It is written in plain R, so it is
slower than `stats::kmeans()`, whose iteration is compiled C code.

## Learn more

The tutorial vignette walks through the function and compares it with
`stats::kmeans()` for correctness and speed:

```r
browseVignettes("kmeans625")
```
