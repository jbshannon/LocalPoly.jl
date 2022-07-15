# LocalPoly

`LocalPoly.jl` is a Julia implementation of the local polynomial regression methods outlined in [Fan and Gijbels (1996)](https://doi.org/10.1201/9780203748725). This package is still experimental, and the API is subject to change.

## Overview

This package provides the function `lpreg`, which computes the local polynomial regression coefficients and standard errors at a vector of evaluation knots. This function also (optionally) implements the linear binning method to speed up the computations by reducing the dimensionality of the data. The number of bins is controlled by the keyword argument  `nbins` (set to 0 for no binning).


## References

1. Fan, J., & Gijbels, I. (1996). Local Polynomial Modelling and its Applications (1st ed.). Chapman & Hall.
