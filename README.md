
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TempStable <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->
<!-- badges: end -->
<!-- Start of my description -->

A collection of methods to estimate parameters of different tempered
stable distributions. Currently, there are three different tempered
stable distributions to choose from: Tempered stable subordinator
distribution, classical tempered stable distribution, normal tempered
stable distribution. The package also provides functions to compute
density and probability functions and tools to run Monte Carlo
simulations.

The main function of this package are briefly described below:

- Main function: TemperedEstim() computes all the information about the
  estimator. It allows the user to choose the preferred method and
  several related options.
- Characteristic function, density function, probability function and
  other functions for every tempered stable distribution mentioned
  above. E.g. charTSS(), dCTS(), …
- Monte Carlo simulation: a tool to run a Monte Carlo simulation
  (TemperedEstim_Simulation()) is provided and can save output files or
  produce statistical summary.

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
# Such a simulation can take a very long time. Therefore, it can make sense to 
# parallelise after Monte Carlo runs. Parallelisation of the simulation is not 
# yet part of the package. 

# For testing purposes, the amount of runs and parameters is greatly reduced. 
# Therefore, the result is not meaningful. To start a meaningful simulation, the
# SampleSize could be, for example, 1000 and MCParam also 1000.
thetaT <- c(1.5,1,1,1,1,0)
res_CTS_ML_size10 <- TemperedEstim_Simulation(ParameterMatrix = rbind(thetaT),
                                               SampleSizes = c(10), MCparam = 10,
                                               TemperedType = "Classic", Estimfct = "ML",
                                               saveOutput = FALSE)
#> ---------------- Alpha=1.5 *** DeltaP=1 *** DeltaM=1 *** LambdaP=1 *** LambdaM=1 *** mu=0 --------------- 
#> *** Iter 1/10 *** Estimated Remaining Time: 0h2min6sec. *** 
#> *** Iter 2/10 *** Estimated Remaining Time: 0h2min4sec. *** 
#> *** Iter 3/10 *** Estimated Remaining Time: 0h4min8sec. *** 
#> *** Iter 4/10 *** Estimated Remaining Time: 0h2min59sec. ***
#> Warning in log(densis): NaNs wurden erzeugt

#> Warning in log(densis): NaNs wurden erzeugt
#> *** Iter 5/10 *** Estimated Remaining Time: 0h0min45sec. ***
#> Warning in log(densis): NaNs wurden erzeugt
#> *** Iter 6/10 *** Estimated Remaining Time: 0h0min30sec. *** 
#> *** Iter 7/10 *** Estimated Remaining Time: 0h0min50sec. ***
#> Warning in log(densis): NaNs wurden erzeugt
#> *** Iter 8/10 *** Estimated Remaining Time: 0h0min15sec. ***
#> Warning in log(densis): NaNs wurden erzeugt
#> *** Iter 9/10 *** Estimated Remaining Time: 0h0min8sec. ***
#> Warning in log(densis): NaNs wurden erzeugt
#> *** Iter 10/10 *** Estimated Remaining Time: 0h0min0sec. ***

colMeans(sweep(res_CTS_ML_size10$outputMat[,9:14],2,thetaT), na.rm = TRUE)
#>     alphaE    delta+E    delta-E   lambda+E   lambda-E        muE 
#> -1.3042931  2.6728708  3.5903234  1.2565026  1.5549668  0.4091464
```
