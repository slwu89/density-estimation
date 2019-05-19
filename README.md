# density-estimation

## tsde94

This is an R package that runs the algorithm described in Spear, R. C., Grieb, T. M., & Shang, N. (1994). Parameter uncertainty and interaction in complex environmental models. Water Resources Research. https://doi.org/10.1029/94WR01732

To install this package, you must have first installed the `Rcpp` and `igraph` packages, which it depends on.  Then you can install it from github directly with `devtools::install_github("slwu89/density-estimation", subdir="tsde94")`

To test it with some sample data provided with the package, you can do:

```R
library(tsde94)
data(sdata)
tree <- TSDE(sdata)
g <- TSDE_maketree(df = tree,pass = sdata,varnames = LETTERS[1:5],plot = TRUE)
```

For help, just do `?TSDE` or `?TSDE_maketree`

## det11

This is an R package that runs the algorithm described in Ram, P., & Gray, A. G. (2011). Density estimation trees. https://doi.org/10.1145/2020408.2020507

It depends on code from mlpack, see http://www.mlpack.org.

To install this package, you must have first installed the `Rcpp`, `RcppArmadillo`, and `igraph` packages, which it depends on. Then you can install it from github directly with `devtools::install_github("slwu89/density-estimation", subdir="det11")`

To test it with some sample data, you can do (remember to transpose the data):

```R
library(tsde94)
data(sdata)
library(det11)
tree <- ml_det(dataset = t(sdata))
processed_tree <- DET_maketree(x = tree)
```

For help, just do `?ml_det`
