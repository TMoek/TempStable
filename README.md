
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TempStable

<!-- badges: start -->
<!-- badges: end -->

TODO: Describe what the package does… How to write readme:
<https://monashbioinformaticsplatform.github.io/2017-11-16-open-science-training/topics/rmarkdown.html>

<!-- Start of my description -->

A collection of methods to estimate parameters of different tempered
stable distributions. Currently, there are three different tempered
stable distributions to choose from: Tempered stable subordinator
distribution, classical tempered stable distribution, normal tempered
stable distribution. The package also provides functions to compute
characteristic functions and tools to run Monte Carlo simulations.

The main function of this package are briefly described below:

-   Main function: TemperedEstim() computes all the information about
    the estimator. It allows the user to choose the preferred method and
    several related options.
-   Characteristic function, density function, probability function and
    other functions for every tempered stable distribution mentioned
    above. E.g. charSTS(), dCTS(), …
-   Monte Carlo simulation: a tool to run a Monte Carlo simulation
    (TemperedEstim_Simulation()) is provided and can save output files
    or produce statistical summary.

The package was developed by Till Massing and Cedric Jüssen and is
structurally based on the “StableEstim” package by Tarak Kharrat and
Georgi N. Boshnakov.

<!-- End of my description -->

## Installation

You can install the development version of TempStable from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cedricjuessen/TempStable")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TempStable)
## basic example code
```
