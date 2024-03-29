% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basicfcts_TS.R
\name{dCTS}
\alias{dCTS}
\title{Density function of the classical tempered stable (CTS) distribution}
\usage{
dCTS(
  x,
  alpha = NULL,
  deltap = NULL,
  deltam = NULL,
  lambdap = NULL,
  lambdam = NULL,
  mu = NULL,
  theta = NULL,
  dens_method = "FFT",
  a = -20,
  b = 20,
  nf = 2048,
  ...
)
}
\arguments{
\item{x}{A numeric vector of quantiles.}

\item{alpha}{Stability parameter. A real number between 0 and 2.}

\item{deltap}{Scale parameter for the right tail. A real number > 0.}

\item{deltam}{Scale parameter for the left tail. A real number > 0.}

\item{lambdap}{Tempering parameter for the right tail. A real number > 0.}

\item{lambdam}{Tempering parameter for the left tail. A real number > 0.}

\item{mu}{A location parameter, any real number.}

\item{theta}{Parameters stacked as a vector.}

\item{dens_method}{Algorithm for numerical evaluation. Choose between \code{
"FFT"} (default) and \code{"Conv"}.}

\item{a}{Starting point of FFT, if \code{dens_method == "FFT"}. -20
by default.}

\item{b}{Ending point of FFT, if \code{dens_method == "FFT"}. 20
by default.}

\item{nf}{Pieces the transformation is divided in. Limited to power-of-two
size. 2048 by default.}

\item{...}{Possibility to modify \code{\link[=charCTS]{charCTS()}}.}
}
\value{
As \code{x} is a numeric vector, the return value is also a numeric
vector of densities.
}
\description{
The probability density function (PDF) of the classical tempered stable
distributions is not available in closed form.
Relies on fast Fourier transform (FFT) applied to the characteristic
function.
}
\details{
\code{theta} denotes the parameter vector \code{(alpha, deltap, deltam,
lambdap, lambdam, mu)}. Either provide the parameters individually OR
provide \code{theta}. Methods include the FFT or alternatively by convolving
two totally positively skewed tempered stable distributions, see Massing
(2022).

The "FFT" method is automatically selected for Mac users, as the "Conv"
method causes problems.
}
\examples{
x <- seq(0,15,0.25)
y <- dCTS(x,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048)
plot(x,y)

}
\references{
Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
}
