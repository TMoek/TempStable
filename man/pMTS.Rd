% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basicfcts_TS_extended.R
\name{pMTS}
\alias{pMTS}
\title{Cumulative probability function of the  modified tempered stable (MTS)
distribution}
\usage{
pMTS(
  q,
  alpha = NULL,
  delta = NULL,
  lambdap = NULL,
  lambdam = NULL,
  mu = NULL,
  theta = NULL,
  dens_method = "FFT",
  a = -40,
  b = 40,
  nf = 2048,
  ...
)
}
\arguments{
\item{q}{A vector of real numbers where the CF is evaluated.}

\item{alpha}{Stability parameter. A real number between 0 and 2.}

\item{delta}{Scale parameter. A real number > 0.}

\item{lambdap, lambdam}{Tempering parameter. A real number > 0.}

\item{mu}{A location parameter, any real number.}

\item{theta}{Parameters stacked as a vector.}

\item{dens_method}{A method to get the density function. Here, only "FFT" is
available.}

\item{a}{Starting point of FFT, if \code{dens_method == "FFT"}. -20
by default.}

\item{b}{Ending point of FFT, if \code{dens_method == "FFT"}. 20
by default.}

\item{nf}{Pieces the transformation is divided in. Limited to power-of-two
size.}

\item{...}{Possibility to modify \code{stats::integrate()}.}
}
\value{
As \code{q} is a numeric vector, the return value is also a numeric
vector of probabilities.
}
\description{
The cumulative probability distribution function (CDF) of the  modified
tempered stable distribution.
}
\details{
\code{theta} denotes the parameter vector \code{(alpha, delta,
lambdap, lambdam, mu)}. Either provide the parameters individually OR
provide \code{theta}.
The function integrates the PDF numerically with \code{integrate()}.
}
