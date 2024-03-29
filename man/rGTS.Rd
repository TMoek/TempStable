% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basicfcts_TS_extended.R
\name{rGTS}
\alias{rGTS}
\title{Function to generate random variates of GTS distribution.}
\usage{
rGTS(
  n,
  alphap = NULL,
  alpham = NULL,
  deltap = NULL,
  deltam = NULL,
  lambdap = NULL,
  lambdam = NULL,
  mu = NULL,
  theta = NULL,
  methodR = "AR",
  k = 10000,
  c = 1
)
}
\arguments{
\item{n}{sample size (integer).}

\item{alphap, alpham}{Stability parameter. A real number between 0 and 2.}

\item{deltap}{Scale parameter for the right tail. A real number > 0.}

\item{deltam}{Scale parameter for the left tail. A real number > 0.}

\item{lambdap}{Tempering parameter for the right tail. A real number > 0.}

\item{lambdam}{Tempering parameter for the left tail. A real number > 0.}

\item{mu}{A location parameter, any real number.}

\item{theta}{Parameters stacked as a vector.}

\item{methodR}{A String. Either "TM","AR" or "SR".}

\item{k}{integer: the level of truncation, if \code{methodR == "SR"}. 10000
by default.}

\item{c}{A real number. Only relevant for \code{methodR == "AR"}.
1 by default.}
}
\value{
Generates \code{n} random numbers of the CTS distribution.
}
\description{
Generates \code{n} random numbers distributed according to the generalized
classical tempered stable (GTS) distribution.
}
\details{
\code{theta} denotes the parameter vector \code{(alphap, alpham, deltap,
deltam, lambdap, lambdam, mu)}. Either provide the parameters individually OR
provide \code{theta}.
"AR" stands for the approximate Acceptance-Rejection Method and "SR" for a
truncated infinite shot noise series representation.

It is recommended to check the generated random numbers once for each
distribution using the density function. If the random numbers are shifted,
e.g. for the method "SR", it may be worthwhile to increase k.

For more details, see references.
}
\examples{
rGTS(2,1.5,0.5,1,1,1,1,0,NULL,"SR")
rGTS(2,1.5,0.5,1,1,1,1,1,NULL,"aAR")

}
\references{
Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'

Kawai, R & Masuda, H (2011), 'On simulation of tempered stable random
variates' \doi{10.1016/j.cam.2010.12.014}

Hofert, M (2011), 'Sampling Exponentially Tilted Stable Distributions'
\doi{10.1145/2043635.2043638}
}
\seealso{
\code{\link[copula:retstable]{copula::retstable()}} as "TM" uses this function and \code{\link[=rCTS]{rCTS()}}.
}
