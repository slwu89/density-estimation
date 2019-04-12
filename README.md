# density-estimation

## tsde94

This is an R package that runs the algorithm described in Spear, R. C., Grieb, T. M., & Shang, N. (1994). Parameter uncertainty and interaction in complex environmental models. Water Resources Research. https://doi.org/10.1029/94WR01732

To test it with some sample data provided with the package, you can do:

```R
library(tsde94)
data(sdata)
TSDE(sdata)
```
