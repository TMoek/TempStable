# TempStable 0.2.1

## Breaking changes

* Package moved on GitHub from cedricjuessen to TMoek

## Minor improvements and fixes 

* @aliases TempStable-package added in documentation

* Fisher Information for TSS patched

# TempStable 0.2.0

## New features

* 4 new tempered stable distribution were implemented: Generalized classical TS, 
  modified TS, rapid decreasing TS, and Kim-Rachev TS.

* With each new distribution, among other things, function for determining the 
  characteristic function and density function were implemented.

* In addition, a random variates generation function was implemented for each 
  new distribution.

* All new distributions were also implemented in the Monte Carlo and the 
  estimation function.

## Changes to existing functions and code

* A new method to generate random variates of a tempered stable distribution
  has been added to the existing TS distributions, which is also the current
  default method

* For estimator CGMM for NTS there was a small error in calculation: There was
  a small error in the calculation of the variable Cmat in CGMM_NTS.R line 359 
  and 360


# TempStable 0.1.1

## Small patch due to CRAN request (no functional changes)

* Changes in DESCRIPTION file due to "imports" and "depends" 

* Some changes in test files


# TempStable 0.1.0

* First version, no downstream dependencies.
