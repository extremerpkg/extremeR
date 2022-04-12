# extremeR <img src="man/figures/hex.png" align="right" width="120"/>

Repository for package accompanying paper "Explaining Recruitment to Violent Extremism: A Bayesian Case-Control Approach."

## Overview

The goal of this package is to make it easy for the user to:

1.  **Prepare** data for estimation procedure (`data.prep()`)
2.  **Estimate** Bayesian case-control models (`fit()`)

You can begin using the package by loading into memory using:

``` r
devtools::install_github("extremerpkg/extremeR")
```

## Getting started

```{r}
library(extremeR)
```

There are two principal <tt>extremeR</tt> functions.

-   **Preparing** the data is achieved with `data.prep()`, which allows the user to bundle the main relevant data source (case-control data and shapefile data) into a list object for estimation.
-   **Estimating** models with different specifications is achieved with `fit()`, which allows the user to estimate models described in the associated article, and to specify hierarchical structure, priors, contamination layer etc.
