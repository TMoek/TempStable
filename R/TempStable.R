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
#' @section Brief describtion of functions:
#'
#'  \strong{TemperedEstim()} computes all the information about the
#' estimator. It allows the user to choose the preferred method and several
#' related options.
#'
#' Characteristic function, density function, probability function and other
#' functions for every tempered stable distribution mentioned above.
#' E.g. charSTS(), dCTS(), ...
#'
#' \strong{Monte Carlo simulation:} a tool to run a Monte Carlo simulation
#' (TemperedEstim_Simulation()) is provided and can save output files or
#' produce statistical summary.
#'
#'
#' @examples
#' TODO
#'
#'
#'
#'
#' @docType package
#' @name TempStable
NULL
