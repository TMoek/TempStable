#' TempStable: A collection of methods to estimate parameters of different
#' tempered stable distributions.
#'
#' A collection of methods to estimate parameters of different tempered stable
#' distributions. Currently, there are three different tempered stable
#' distributions to choose from: Tempered stable subordinator distribution,
#' classical tempered stable distribution, normal tempered stable distribution.
#' The package also provides functions to compute characteristic functions and
#' tools to run Monte Carlo simulations.
#'
#' The package was developed by Till Massing and Cedric Juessen and is
#' structurally based on the "StableEstim" package by Tarak Kharrat and Georgi
#' N. Boshnakov.
#'
#' @section Brief description of functions:
#'
#'  \strong{TemperedEstim()} [TemperedEstim()]computes all the information about the
#' estimator. It allows the user to choose the preferred method and several
#' related options.
#'
#' Characteristic function, density function, probability function and other
#' functions for every tempered stable distribution mentioned above.
#' E.g. [charTSS()], [dCTS()], ...
#'
#' \strong{Monte Carlo simulation:} a tool to run a Monte Carlo simulation
#' [TemperedEstim_Simulation()] is provided and can save output files or
#' produce statistical summary.To parallelize this function, you can use
#' [parallelizeMCsimulation()].
#'
#'
#' @examples
#' ## basic example code
#' # Such a simulation can take a very long time. Therefore, it can make sense
#' # to parallelize after Monte Carlo runs. Parallelization of the simulation is
#' # now possible with [parallelizeMCsimulation()].
#'
#' # For testing purposes, the amount of runs and parameters is greatly reduced.
#' # Therefore, the result is not meaningful. To start a meaningful simulation,
#' # the SampleSize could be, for example, 1000 and MCParam also 1000.
#' \donttest{
#' thetaT <- c(1.5,1,1,1,1,0)
#' res_CTS_ML_size4 <- TemperedEstim_Simulation(ParameterMatrix =
#'                                                 rbind(thetaT),
#'                                               SampleSizes = c(4),
#'                                               MCparam = 4,
#'                                               TemperedType = "CTS",
#'                                               Estimfct = "ML",
#'                                               saveOutput = FALSE)
#'
#' colMeans(sweep(res_CTS_ML_size4$outputMat[,9:14],2,thetaT), na.rm = TRUE)
#' }
#'
#' @docType package
#' @name TempStable
#' @aliases TempStable-package
NULL
