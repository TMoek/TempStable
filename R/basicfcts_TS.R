imagN <- (0 + (0 + (0 + (0+1i))))

####Subordinator Tempered Stable (STS)####

#' Characteristic function of the tempered stable subordinator
#'
#' When \code{alpha} is a real number between 0 and 1, the characterisic
#' function of the tempered stable distribution (\code{alpha}, \code{delta},
#' \code{lambda}) can be simulated exactly through acceptance-rejection
#' sampling.
#' For whole derivation and meaning of single letters see also Kawai et. Masuda
#' (2011).
#' Basically, the derivation starts from a one-dimensional stable distribution
#' which is later added with a centered and totally positively skewed tempered
#' stable (Levy) process.
#'
#' \deqn{ \varphi_{TSS}(t;\theta):=\mathbb{E}_{\theta}\left[
#' \mathrm{e}^{\mathrm{i}tY}\right]= \exp\left(\delta\Gamma(-\alpha)
#' \left((\lambda-\mathrm{i}t)^{\alpha}-\lambda^{\alpha}\right)\right)}
#'
#' @param t A positive integer.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#'
#' @return The result of spectral positive tempered stable process
#'
#' @references
#' Massing, T. (2022), 'Parametric Estimation of Tempered Stbale Laws'
#'
#' Reiichiro K. & Hiroki M. (2011), 'On simulation of tempered stable random
#' variates' \url{https://doi.org/10.1016/j.cam.2010.12.014}
#'
#' @examples
#' charSTS(1000,0.5,1,0.3)
#' charSTS(500,0.9,1,0.3)
#'
#' @export
charSTS <- function(t, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL){
    if ((missing(alpha) | missing(delta) | missing(lambda)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, delta, lambda)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(delta) | missing(lambda)) {
      alpha <- theta[1]
      delta <- theta[2]
      lambda <- theta[3]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)
    return(exp(delta * gamma(-alpha) *
                 ((lambda - t * imagN)^alpha - lambda^alpha)))
}

#' Density function of the tempered stable subordinator (TSS) distribution
#'
#' The probability density function of tempered stable distributions is
#' generally not available in closed form. However, many software packages (like
#' the \code{stabledist} package) have fast computation routines based on series
#' or integral representations. Combining such series representations with the
#' density function of the TSS distribution, we obtain a series representation
#' for the TSS distribution.
#'
#' \deqn{f_{TSS}(y;\theta)=\mathrm{e}^{-\lambda y-\lambda^{\alpha}\delta
#' \Gamma(-\alpha)}\frac{-1}{\pi}\sum_{k=1}^{\infty}\frac{(-1)^k}{k!}
#' \Gamma(1+\alpha k)\Gamma(1-\alpha)^k\left(\frac{\delta}{\alpha}\right)^ky^
#' {-(1+\alpha k)}\sin(\alpha\pi k)}
#'
#' @param x A numeric vector of quantiles.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#'
#' @return As \code{x} is a numeric vector, the return value is also a numeric
#' vector.
#'
#' @references
#' Massing, Till (2022), 'Parametric Estimation of Tempered Stbale Laws'
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distributions and
#' processes' \url{https://doi.org/10.1016/j.spa.2013.06.012}
#'
#' @examples
#' dSTS(1000,0.5,1,0.3)
#' dSTS(1,0.5,1,0.3)
#'
#' @importFrom stabledist dstable
#' @export
dSTS <- function(x, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL) {
    if ((missing(alpha) | missing(delta) | missing(lambda)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, delta, lambda)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(delta) | missing(lambda)) {
      alpha <- theta[1]
      delta <- theta[2]
      lambda <- theta[3]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)
    sigma <- ((1/alpha * gamma(1 - alpha) * delta *
                 cos(alpha * pi/2))^(1/alpha))

    return(dstable(x, alpha, 1, sigma, 0, pm = 1) *
             exp(-lambda * x - lambda^alpha * delta * gamma(-alpha)))

}

#' Quantile function of the tempered stable subordinator distribution
#'
#' This function returns the quantile function of the tempered stable
#' suborinator distribution.
#'
#'
#' @param q A numeric vector of positive quantiles.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#' @param pmethod A string. If not "integrate", the function \code{chartocdf()}
#' will be triggered.
#' @param N integer: the number of replications, if
#' \code{pmethod != "integrate"}. 10000 by default.
#' @param ... Possibility to modify \code{stats::integrate()}.
#'
#' @return  As \code{x} is a numeric vector, the return value is also a numeric
#' vector.
#'
#' @seealso
#' See also the [dSTS()] density-function.
#'
#' @examples
#' pSTS(3,0.7,1.354,0.3)
#' pSTS(1,0.5,10,300,NULL,"")
#'
#' @export
pSTS <- function(q, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL,
                   pmethod = "integrate", N = 10000, ...) {
    if ((missing(alpha) | missing(delta) | missing(lambda)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, delta, lambda)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(delta) | missing(lambda)) {
      alpha <- theta[1]
      delta <- theta[2]
      lambda <- theta[3]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)

    if(pmethod =="integrate"){
    p <- sapply(q, function(z, ...) min(stats::integrate(dSTS, 0, z,
                                                         alpha = alpha,
                                                         delta = delta,
                                                         lambda = lambda,
                                                         ...)$value, 1 - 1e-07))

  } else {
    m <- gamma(1 - alpha)*delta/lambda^(1-alpha)
    s <- sqrt(gamma(2-alpha)*delta/lambda^(2-alpha))
    p <- sapply(q, chartocdf, N = N, m = m, s = s, char = charSTS,
                alpha = alpha, delta = delta, lambda = lambda)
  }

  return(p)
}

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param n A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#' @param method A String. Either "AR" or "SR".
#' @param k integer: the number of replications, if \code{method == "SR"}. 100
#' by default.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' rSTS(100,0.5,1,1)
#' rSTS(100,0.5,1,1,NULL,"SR",50)
#'
#' @export
rSTS <- function(n, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL,
                   method = "AR", k = 100, ...) {
    if ((missing(alpha) | missing(delta) | missing(lambda)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, delta, lambda)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(delta) | missing(lambda)) {
      alpha <- theta[1]
      delta <- theta[2]
      lambda <- theta[3]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)
    x <- switch(method,
                AR = rSTS_AR(n = n, alpha = alpha, delta = delta,
                             lambda = lambda),
                SR = rSTS_SR(n = n, alpha = alpha, delta = delta,
                             lambda = lambda, k = k))
    return(x)
}

# No export.
#' @importFrom stabledist rstable
rSTS_AR <- function(n, alpha, delta, lambda) {
    sigma <- ((1/alpha * gamma(1 - alpha) * delta *
                 cos(alpha * pi/2))^(1/alpha))

    returnVector <- c()

    i <- 0
    while (i < n) {
        U <- 2
        V <- 0

        while (U > exp(-min(c(lambda * V, 700)))) {
            U <- stats::runif(1, 0, 1)
            V <- rstable(1, alpha, 1, sigma, 0, pm = 1)
        }

        returnVector <- append(returnVector, V)
        i <- i + 1
    }

    return(returnVector)

}

# No export.
rSTS_SR <- function(n, alpha, delta, lambda, k) {
    base::replicate(n = n, rSTS_SR1(alpha = alpha, delta = delta,
                                    lambda = lambda, k = k))

}

# No export.
rSTS_SR1 <- function(alpha, delta, lambda, k) {
    parrivalslong <- cumsum(stats::rexp(k * 1.1))
    parrivals <- parrivalslong[parrivalslong <= k]
    E1 <- stats::rexp(length(parrivals))
    U <- stats::runif(length(parrivals))
    X <- cbind((alpha * parrivals/delta)^(-1/alpha), E1 * U^(1/alpha)/lambda)
    return(sum(apply(X, 1, FUN = min)))
}

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param p A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#' @param qmin,qmax Limits of the interval.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' qSTS(0.5,0.5,5,0.01)
#' qSTS(0.5,0.9,1,10,NULL)
#'
#' @importFrom stabledist qstable
#' @export
qSTS <- function(p, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL,
                  qmin = NULL, qmax = NULL, ...) {
    if (missing(alpha) & missing(delta) & missing(lambda) & is.null(theta))
      stop("No parameters supplied")
    theta0 <- c(alpha, delta, lambda)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(delta) | missing(lambda)) {
      alpha <- theta[1]
      delta <- theta[2]
      lambda <- theta[3]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)

    froot <- function(x, y, ...){
      pSTS(q = x, alpha = alpha, delta = delta, lambda = lambda, ...) - y
    }

    qroot <- function(y, qmin, qmax) {
      if(is.null(qmin)|is.null(qmax)){
        qmini <- 0
        qmaxi <-  abs(qstable(y, alpha = alpha, beta = 1,
                        gamma = min(((1/alpha*gamma(1-alpha)*(delta)*
                                        cos(alpha*pi/2))^(1/alpha)),100),
                        delta = 0, pm = 1))
      } else {
        qmini <- qmin
        qmaxi <- qmax
      }
      .unirootNA(froot,interval = c(qmini,qmaxi),
                              extendInt = "yes", y = y, ...)}

    q <- sapply(p, qroot, qmin = qmin, qmax = qmax)
    return(q)
}

##### Classical Tempered Stable (CTS)#####

#' Characteristic function of the classical tempered stable (CTS) distribution
#'
#' Let $X\sim CTS(\alpha,\delta_+,\delta_-,\lambda_+,\lambda_-,\mu)$.
#' The characteristic function is given by (Kuechler & Tappe 2013).
#'
#' \deqn{\varphi_{CTS}(t;\theta):=\mathbb{E}_{\theta}\left
#' [\mathrm{e}^{\mathrm{i}tX}\right]&=\exp\left(\mathrm{i}t\mu+\delta_+
#' \Gamma(-\alpha)\left((\lambda_+-\mathrm{i}t)^{\alpha}-\lambda_+^{\alpha}+
#' \mathrm{i}t\alpha\lambda_+^{\alpha-1}\right)\right.\\&\ \left. +\delta_-
#' \Gamma(-\alpha)\left((\lambda_-+\mathrm{i}t)^{\alpha}-\lambda_-^{\alpha}-
#' \mathrm{i}t\alpha\lambda_-^{\alpha-1}\right)\right)}
#'
#' @param t A positive integer.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#'
#' @return The result of spectral positive tempered stable process
#'
#' @references
#' Massing, T. (2022), 'Parametric Estimation of Tempered Stbale Laws';
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distributions and
#' processes' \url{https://doi.org/10.1016/j.spa.2013.06.012}
#'
#' @examples
#' charCTS(2,0.5,1,1.5,2,2.5,3)
#' charCTS(2,1.5,200,150,300,2500,1,NULL)
#'
#' @export
charCTS <- function(t, alpha = NULL, deltap = NULL, deltam = NULL,
                    lambdap = NULL, lambdam = NULL, mu = NULL, theta = NULL) {
    if ((missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
         missing(lambdam) | missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, deltap, deltam, lambdap, lambdam, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
        missing(lambdam) | missing(mu)) {
      alpha <- theta[1]
      deltap <- theta[2]
      deltam <- theta[3]
      lambdap <- theta[4]
      lambdam <- theta[5]
      mu <- theta[6]
    }
    stopifnot(0 < alpha, alpha < 2, 0 < deltap, 0 < deltam, 0 < lambdap, 0 <
                lambdam)
    return(exp(imagN * t * mu + deltap * gamma(-alpha) *
                 ((lambdap - imagN * t)^alpha - lambdap^alpha + imagN * t *
                    alpha * lambdap^(alpha - 1)) +
                 deltam * gamma(-alpha) *
                 ((lambdam + imagN * t)^alpha - lambdam^alpha - imagN * t *
                    alpha * lambdam^(alpha - 1))))

}

#' Density function of the classic tempered stable (CTS) distribution
#'
#' The probability density function of tempered stable distributions is
#' generally not available in closed form. Crucially, even a simple relationship
#' with a stable density as for the TSS distribution is not available.
#' For numerical evaluations it is therefore necessary to rely on algorithms
#' like the fast Fourier transform applied to the characteristic function.
#'
#' @param x A numeric vector of quantiles.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param dens_method Algorithm for numerical evaluation. Choose between \code{
#' "FFT"} and \code{"Conv"}.
#' @param a Starting point of FFT, if \code{dens_method == "FFT"}. -20
#' by default.
#' @param b Ending point of FFT, if \code{dens_method == "FFT"}. 20
#' by default.
#' @param nf Pieces the transformation is divided in. Limited to power-of-two
#' size.
#'
#' @return As \code{x} is a numeric vector, the return value is also a numeric
#' vector.
#'
#' @references
#' Massing, Till (2022), 'Parametric Estimation of Tempered Stbale Laws'
#'
#' @examples
#' dCTS(1,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048)
#' dCTS(1,0.6,1,1,1,1,1,NULL,"Conv")
#'
#' @export
dCTS <- function(x, alpha = NULL, deltap = NULL, deltam = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, dens_method = "FFT",
                 a = -20, b = 20, nf = 2048) {
    if ((missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
         missing(lambdam) | missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, deltap, deltam, lambdap, lambdam, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
        missing(lambdam) | missing(mu)) {
      alpha <- theta[1]
      deltap <- theta[2]
      deltam <- theta[3]
      lambdap <- theta[4]
      lambdam <- theta[5]
      mu <- theta[6]
    }
    stopifnot(0 < alpha, alpha < 2, 0 < deltap, 0 < deltam, 0 < lambdap, 0 <
                lambdam)

    if (dens_method == "FFT") {
        d <- sapply(x, dCTS_FFT, alpha = alpha, deltap = deltap,
                    deltam = deltam, lambdap = lambdap, lambdam = lambdam,
                    mu = mu, a = a, b = b, nf = nf)
    } else {
        d <- sapply(x, dCTS_Conv, alpha = alpha, deltap = deltap,
                    deltam = deltam, lambdap = lambdap, lambdam = lambdam,
                    mu = mu)
    }
    return(d)
}


# No export.
dCTS_FFT <- function(x, alpha, deltap, deltam, lambdap, lambdam,
                          mu, a, b, nf) {
    dx <- ((b - a)/nf)
    dt <- (2 * pi/(nf * dx))

    sq <- seq(from = 0, to = nf - 1, 1)
    t <- (-nf/2 * dt + sq * dt)
    xgrid <- (a + dx * sq)

    cft <- charCTS(t, alpha, deltap, deltam, lambdap, lambdam, mu)

    tX <- (exp(-imagN * sq * dt * a) * cft)

    # fast Fourier transform
    y <- stats::fft(tX, inverse = FALSE)

    densityW <- Re((dt/(2 * pi)) * exp(-imagN * (nf/2 * dt) * xgrid) * y)

    return (as.numeric(stats::approx(xgrid, densityW, xout=x,
                              yleft = 1e-18, yright = 1e-18)[2]))
}

# No export.
#' @importFrom stabledist dstable
dCTS_Conv <- function(x, alpha, deltap, deltam, lambdap, lambdam, mu) {
    integrandDTSclass <- function(y) {

        Sigmap <- ((1/alpha * gamma(1 - alpha) * deltap *
                      cos(alpha * pi/2))^(1/alpha))
        Sigmam <- ((1/alpha * gamma(1 - alpha) * deltam *
                      cos(alpha * pi/2))^(1/alpha))
        z <- x - mu

        exp((-lambdap * y - lambdap^(alpha) * deltap * gamma(-alpha) *
               (alpha + 1)) +
              (-lambdam * (y + z) - lambdam^(alpha) * deltam * gamma(-alpha) *
                 (alpha + 1)) +
              log(dstable((y - gamma(1 - alpha) * deltap * lambdap^(alpha - 1)),
                          alpha, 1, Sigmap, 0, pm = 1) *
                    dstable((y - gamma(1 - alpha) * deltam *
                               lambdam^(alpha - 1) + z),
                            alpha, 1, Sigmam, 0, pm = 1)))

    }

    return(stats::integrate(f = integrandDTSclass, lower = -Inf, upper = Inf,
                            abs.tol = 0)$value)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param q A gap holder.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param a A gap holder. -20 by default.
#' @param b A gap holder. 20 by default.
#' @param nf A gap holder.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' pCTS(0.5,0.5,1,1,1,1,1)
#' pCTS(0.5,0.9,1,2,3,4,5)
#'
#' @export
pCTS <- function(q, alpha = NULL, deltap = NULL, deltam = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL,
                 a = -40, b = 40, nf = 2^13, ...) {
    if ((missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
         missing(lambdam) | missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, deltap, deltam, lambdap, lambdam, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
        missing(lambdam) | missing(mu)) {
      alpha <- theta[1]
      deltap <- theta[2]
      deltam <- theta[3]
      lambdap <- theta[4]
      lambdam <- theta[5]
      mu <- theta[6]
    }
    stopifnot(0 < alpha, alpha < 2, 0 < deltap, 0 < deltam, 0 < lambdap, 0 <
                lambdam)

    p <- numeric(length(q))

    p <- sapply(q,
                function(z) {
                  if(z<a) 0
                  else if (z>b) 1
                  else min(stats::integrate(dCTS, lower = a, upper = z, alpha = alpha,
                                 deltap = deltap, deltam = deltam,
                                 lambdap = lambdap, lambdam = lambdam,
                                 mu = mu,...)$value, 1 - 1e-07)})

    return(p)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param n A gap holder.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param method A String. Either "aAR" or "SR".
#' @param k integer: the number of replications, if \code{method == "SR"}. 100
#' by default.
#' @param c A real number. Only relevant for \code{method == "aAR"}.
#' TODO(Till): Choose a suitable default value.
#'
#' @return Gap holder for return.
#'
#' @examples
#' rCTS(10,0.5,1,1,1,1,1,NULL,"SR",10)
#' rCTS(10,0.5,1,1,1,1,1,NULL,"aAR")
#'
#' @export
rCTS <- function(n, alpha = NULL, deltap = NULL, deltam = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, method = "aAR",
                 k = 100, c = 1) {
    if ((missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
         missing(lambdam) | missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, deltap, deltam, lambdap, lambdam, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
        missing(lambdam) | missing(mu)) {
      alpha <- theta[1]
      deltap <- theta[2]
      deltam <- theta[3]
      lambdap <- theta[4]
      lambdam <- theta[5]
      mu <- theta[6]
    }
    stopifnot(0 < alpha, alpha < 2, 0 < deltap, 0 < deltam, 0 < lambdap, 0 <
                lambdam)

    x <- switch(method,
                aAR = rCTS_aAR(n = n, alpha = alpha, deltap = deltap,
                                       deltam = deltam, lambdap = lambdap,
                                       lambdam = lambdam, mu = mu, c = c),
                SR = rCTS_SR(n = n, alpha = alpha, deltap = deltap,
                             deltam = deltam, lambdap = lambdap,
                             lambdam = lambdam, mu = mu, k = k))
    return(x)
}

# No export.
rCTS_aAR <- function(n, alpha, deltap, deltam, lambdap, lambdam, mu, c) {
    return(mu +
             rCTS_aARp(n, alpha = alpha, delta = deltap, lambda = lambdap,
                       c = c) -
             rCTS_aARp(n, alpha = alpha, delta = deltam, lambda = lambdam,
                       c = c))
}

# No export.
#' @importFrom stabledist rstable
rCTS_aARp <- function(n, alpha, delta, lambda, c) {
    sigma <- ((1/alpha * gamma(1 - alpha) * delta *
                 cos(alpha * pi/2))^(1/alpha))
    cc <- ifelse(alpha < 1, 0, c)
    returnVector <- c()

    i <- 0
    while (i < n) {
        U <- 2
        V <- 0
        Y <- 0

        while (U > exp(-min(c(lambda * (V + cc), 700)))) {
            U <- stats::runif(1, 0, 1)
            V <- rstable(1, alpha, 1, sigma, 0, pm = 1)
            Y <- V - gamma(1 - alpha) * delta * lambda^(alpha - 1)
        }

        returnVector <- append(returnVector, Y)
        i <- i + 1
    }

    return(returnVector)

}

# No export.
rCTS_SR <- function(n, alpha, deltap, deltam, lambdap, lambdam, mu, k) {
    replicate(n = n, rCTS_SR1(alpha = alpha, deltap = deltap, deltam = deltam,
                              lambdap = lambdap, lambdam = lambdam,
                              mu = mu, k = k))

}

# No export.
rCTS_SR1 <- function(alpha, deltap, deltam, lambdap, lambdam, mu, k) {
    return(rCTS_SRp(alpha = alpha, delta = deltap, lambda = lambdap, k = k)
           -rCTS_SRp(alpha = alpha, delta = deltam, lambda = lambdam, k = k)
           +mu)
}



# No export.
#' @importFrom VGAM zeta
rCTS_SRp <- function(alpha, delta, lambda, k) {
    parrivalslong <- cumsum(stats::rexp(k * 2))
    parrivals <- parrivalslong[parrivalslong <= k]
    n <- length(parrivals)
    E1 <- stats::rexp(n)
    U <- stats::runif(n)
    cntr <- ((1:n) * alpha/delta)^(-1/alpha)
    gam <- (delta/alpha)^(1/alpha) * VGAM::zeta(1/alpha) - gamma(1 - alpha) *
      delta * lambda^(alpha - 1)
    xBig <- base::cbind((alpha * parrivals/delta)^(-1/alpha),
                        E1 * U^(1/alpha)/lambda)
    jumps <- apply(xBig, 1, FUN = min) - cntr
    x <- sum(jumps) + gam
    return(x)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param p A gap holder.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param qmin,qmax Limits of the interval.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' \donttest{
#'   qCTS(0.5,1.5,10,10,10,10,10)
#'   qCTS(0.5,1.5,1,1,1,1,1)
#' }
#'
#' @export
qCTS <- function(p, alpha = NULL, deltap = NULL, deltam = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, qmin = NULL,
                 qmax = NULL, ...) {
    if ((missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
         missing(lambdam) | missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, deltap, deltam, lambdap, lambdam, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(deltap) | missing(deltam) | missing(lambdap) |
        missing(lambdam) | missing(mu)) {
      alpha <- theta[1]
      deltap <- theta[2]
      deltam <- theta[3]
      lambdap <- theta[4]
      lambdam <- theta[5]
      mu <- theta[6]
    }
    stopifnot(0 < alpha, alpha < 2, 0 < deltap, 0 < deltam, 0 < lambdap, 0 <
                lambdam)
    froot <- function(x, y, ...){
      pCTS(q = x, alpha = alpha, deltap = deltap, deltam = deltam,
           lambdap = lambdap, lambdam = lambdam, mu = mu, ...) - y
    }

    qroot <- function(y, qmin, qmax) {
      if(is.null(qmin)|is.null(qmax)){
        x1 <- stats::qcauchy(y, location = mu,
                      scale = min(((1/alpha*gamma(1-alpha)*(deltap+deltam)*
                                      cos(alpha*pi/2))^(1/alpha)),100))
        qmini <- min(c(x1,-x1, mu-1, mu+1))
        qmaxi <- max(c(x1,-x1, mu-1, mu+1))
      } else {
        qmini <- qmin
        qmaxi <- qmax
      }
      .unirootNA(froot, interval = c(qmini,qmaxi), extendInt = "yes",
                y = y, ...)}

    q <- sapply(p, qroot, qmin = qmin, qmax = qmax)
    return(q)
}

##### Normal Tempered Stable#####

#' Characteristic function of the normal tempered stable (NTS) distribution
#'
#' We can obtain the NTS distribution by tempering a stable distribution.
#' The corresponding tempering function can be found in Rachev et al. (2011).
#'
#' \deqn{\varphi_{NTS}(t;\theta)=E\left[\mathrm{e}^{\mathrm{i}tZ}\right]= \exp
#' \left(\mathrm{i}t\mu+\delta\Gamma(-\alpha)\left((\lambda-\mathrm{i}t
#' \beta+t^2/2)^{\alpha}-\lambda^{\alpha}\right)\right)
#' }
#'
#' @param t A positive integer.
#' @param alpha A real number between 0 and 1.
#' @param beta Any real number.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#'
#' @return The result of spectral positive tempered stable process
#'
#' @references
#' Massing, T. (2022), 'Parametric Estimation of Tempered Stbale Laws'
#'
#' Rachev, S., Kim, Y., Bianchi, M. & Fabozzi, F. (2011), 'Financial Models with
#' Levy Processes and Volatility Clustering'
#' \url{https://books.google.de/books?id=XKvUUrcS_twC}
#'
#' @examples
#' charNTS(0.1,0.9,10,20,30,40)
#' charNTS(0.1,0.9,1,2,3,4, NULL)
#'
#' @export
charNTS <- function(t, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                      mu = NULL, theta = NULL) {
    if ((missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
         missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, beta, delta, lambda, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
        missing(mu)) {
      alpha <- theta[1]
      beta <- theta[2]
      delta <- theta[3]
      lambda <- theta[4]
      mu <- theta[5]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)
    return(exp(imagN * t * mu + delta * gamma(-alpha) *
                 ((lambda - imagN * t * beta + t^2/2)^alpha - lambda^alpha)))

}

#' Density function of the normal tempered stable (NTS) distribution
#'
#' As for the CTS distribution, the density function is not available in closed
#' form and numerical computation relies on the fast Fourier transform (FFT).
#'
#' @param x A numeric vector of quantiles.
#' @param alpha A real number between 0 and 1.
#' @param beta Any real number.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param dens_method Currently, useless param, as it does nothing and FFT is
#' always used.
#' @param a Starting point of FFT, if \code{dens_method == "FFT"}. -20
#' by default.
#' @param b Ending point of FFT, if \code{dens_method == "FFT"}. 20
#' by default.
#' @param nf Pieces the transformation is divided in. Limited to power-of-two
#' size.
#'
#' @return As \code{x} is a numeric vector, the return value is also a numeric
#' vector.
#'
#' @references
#' Massing, Till (2022), 'Parametric Estimation of Tempered Stbale Laws'
#'
#' @examples
#' dNTS(1,0.8,1,1,1,1)
#' dNTS(0.5,0.5,20,20,20,20, a = -2000, b = 2000, nf = 8192)
#'
#' @export
dNTS <- function(x, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, dens_method = "FFT",
                 a = -20, b = 20, nf = 2048) {
    if ((missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
         missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, beta, delta, lambda, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
        missing(mu)) {
      alpha <- theta[1]
      beta <- theta[2]
      delta <- theta[3]
      lambda <- theta[4]
      mu <- theta[5]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)

    d <- sapply(x, dNTS_FFT, alpha = alpha , beta = beta, delta = delta,
                lambda = lambda, mu = mu, a = a, b = b , nf = nf)
    return(d)
}


# No export
dNTS_FFT <- function(x, alpha, beta, delta, lambda, mu, a, b, nf) {
    dx <- ((b - a)/nf)
    dt <- (2 * pi/(nf * dx))
    sq <- seq(from = 0, to = nf - 1, 1)
    t <- (-nf/2 * dt + sq * dt)
    xgrid <- (a + dx * sq)

    cft <- charNTS(t, alpha, beta, delta, lambda, mu)

    tX <- (exp(-(1i) * sq * dt * a) * cft)

    y <- stats::fft(tX, inverse = FALSE)

    densityW <- Re((dt/(2 * pi)) * exp(-(1i) * (nf/2 * dt) * xgrid) * y)

    # The return value of stats::approx (RStudio) and Interpolation (Wolfram
    # Mathematica) differ slightly. The larger nf, the smaller the difference.
    return (as.numeric(stats::approx(xgrid, densityW, xout=x,
                              yleft = 1e-18, yright = 1e-18)[2]))
}

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param q A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param a A gap holder. -20 by default.
#' @param b A gap holder. 20 by default.
#' @param nf A gap holder. 2048 by default. I don't know if this param is used
#' in this function.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' pNTS(0.1,0.5,1,1,1,1)
#' pNTS(0.1,0.5,1,1,1,1,NULL, -20, 20, 2^6)
#'
#' @export
pNTS <- function(q, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, a = -40, b = 40, nf = 2^13, ...) {
    if ((missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
         missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, beta, delta, lambda, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
        missing(mu)) {
      alpha <- theta[1]
      beta <- theta[2]
      delta <- theta[3]
      lambda <- theta[4]
      mu <- theta[5]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)

    p <- numeric(length(q))

    p <- sapply(q,
              function(z) {
                if(z<a) 0
                else if
                (z>b) 1
                else min(stats::integrate(dNTS, lower = a, upper = z, alpha = alpha,
                               beta = beta, delta = delta, lambda = lambda,
                               mu = mu, ...)$value, 1-1e-7)})

    return(p)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param n A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder.
#' @param delta  A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param method A String. Either "AR" or "SR". "AR" by default.
#' @param k integer: the number of replications, if \code{method == "SR"}. 100
#' by default.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' rNTS(100, 0.5, 1,1,1,1)
#' rNTS(10, 0.6, 0,1,1,0)
#' rNTS(10, 0.5, 1,1,1,1, NULL, "SR", 100)
#'
#' @export
rNTS <- function(n, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, method = "AR", k = 100, ...) {
    if ((missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
         missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, beta, delta, lambda, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
        missing(mu)) {
      alpha <- theta[1]
      beta <- theta[2]
      delta <- theta[3]
      lambda <- theta[4]
      mu <- theta[5]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)
    x <- switch(method,
                AR = rNTS_AR(n = n, alpha = alpha, beta = beta,
                             delta = delta, lambda = lambda, mu = mu, ...),
                SR = rNTS_SR(n = n, alpha = alpha, beta = beta, delta = delta,
                             lambda = lambda, mu = mu, k = k))
    return(x)
}

# No export.
rNTS_AR <- function(n, alpha, beta, delta, lambda, mu) {
    z <- stats::rnorm(n = n)
    y <- rSTS_AR(n = n, alpha = alpha, delta = delta, lambda = lambda)
    x <- sqrt(y) * z + beta * y + mu
    return(x)
}


# No export.
rNTS_SR <- function(n, alpha, beta, delta, lambda, mu, k) {
    z <- stats::rnorm(n = n)
    y <- rSTS_SR(n = n, alpha = alpha, delta = delta, lambda = lambda, k = k)
    x <- sqrt(y) * z + beta * y + mu
    return(x)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param p A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder.
#' @param delta  A real number > 0.
#' @param lambda A  real number >= 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param qmin,qmax Limits of the interval.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' qNTS(0.1,0.5,1,1,1,1)
#' qNTS(0.3,0.6,1,1,1,1,NULL)
#'
#' @export
qNTS <- function(p, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, qmin = NULL, qmax = NULL, ...) {
    if ((missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
         missing(mu)) & is.null(theta))
      stop("No or not enough parameters supplied")
    theta0 <- c(alpha, beta, delta, lambda, mu)
    if (!is.null(theta) & !is.null(theta0)) {
      if (!all(theta0 == theta))
        stop("Parameters do not match")
    }
    if (missing(alpha) | missing(beta) | missing(delta) | missing(lambda) |
        missing(mu)) {
      alpha <- theta[1]
      beta <- theta[2]
      delta <- theta[3]
      lambda <- theta[4]
      mu <- theta[5]
    }
    stopifnot(0 < alpha, alpha < 1, 0 < delta, 0 < lambda)

    froot <- function(x, y, ...){
      pNTS(q = x, alpha = alpha, beta = beta, delta = delta,
           lambda = lambda, mu = mu, ...) - y
    }

    qroot <- function(y, qmin, qmax) {
      if(is.null(qmin)|is.null(qmax)){
        x1 <- stats::qcauchy(y, location = mu,
                      scale = min(((1/alpha*gamma(1-alpha)*(delta)*
                                      cos(alpha*pi/2))^(1/alpha)),100))
        qmini <- min(c(x1,-x1, mu-1, mu+1))
        qmaxi <- max(c(x1,-x1, mu-1, mu+1))
      } else {
        qmini <- qmin
        qmaxi <- qmax
      }
      .unirootNA(froot,interval = c(qmini,qmaxi), extendInt = "yes",
                 y = y, ...)}

    q <- sapply(p, qroot, qmin = qmin, qmax = qmax)
    return(q)
}

##### CGMY#####

# param t A gap holder.
# param C  A real number > 0.
# param G,M A  real number > 0.
# param Y A real number between 0 and 2.
# examples
# charCGMY(1,1,1,1,0.5)
# charCGMY(10,1,1,1,0.5)
charCGMY <- function(t, C, G, M, Y) {
    charCTS(t = t, alpha = Y, deltap = C, deltam = C, lambdap = G, lambdam = M,
            mu = 0)
}

# examples
# dCGMY(1,1,1,1,0.5)
# dCGMY(1,1,1,1,0.5,"", -2, 2, 2^4)
dCGMY <- function(x, C, G, M, Y, dens_method = "FFT", a = -20, b = 20,
                  nf = 2048, ...) {
    dCTS(x = x, alpha = Y, deltap = C, deltam = C, lambdap = G, lambdam = M,
         mu = 0, dens_method = dens_method, a = a, b = b, nf = 2048, ...)
}



# examples
# rCGMY(100,1,1,1,0.5)
# rCGMY(100,1,1,1,0.5, "SR", k= 100)
rCGMY <- function(n, C, G, M, Y, method = "SR", k = 100, ...) {
    rCTS(n = n, alpha = Y, deltap = C, deltam = C, lambdap = G, lambdam = M,
         mu = 0, method = method, k = k, ...)
}


#####Bohman Char to CDF#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param q A gap holder.
#' @param N A gap holder.
#' @param m A gap holder.
#' @param s A gap holder.
#' @param char A gap holder.
#' @param ... A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' chartocdf(0.5,10,1,1,charNTS, alpha=0.5, beta = 1, delta = 1, lambda = 1,
#'  mu = 1)
#' chartocdf(0.5,10,1,1,charCTS, alpha=0.5, deltap = 1, deltam = 1, lambdap = 1,
#'  lambdam = 1, mu = 1)
#'
#' @export
chartocdf <- function(q, N, m , s, char, ...){
    lambda_bohman <- pi/(m+3*s+2*q)

    summand_fct <- function(n, N, q, char, ...){
      sin(lambda_bohman*n*q)/n*c_err(n/N)*Re(char(lambda_bohman*n, ...))}

    summands <- sapply(X = 1:N,
                       FUN = summand_fct, N = N, q = q, char = char, ...)

    cdf_q <- lambda_bohman*q/pi + 2/pi * sum(summands)

    return(cdf_q)
}


# No export.
c_err <- function(t) {if(abs(t)<1){(1-t)*cos(pi*t)+1/pi*abs(sin(pi*t))} else 0}


#' @importFrom stats uniroot
.unirootNA <-
  function(f,
           interval,
           ...,
           lower = min(interval),
           upper = max(interval),
           f.lower = f(lower, ...),
           f.upper = f(upper, ...),
           extendInt = c("no", "yes", "downX", "upX"),
           check.conv = FALSE,
           tol = .Machine$double.eps ^ 0.25,
           maxiter = 1000,
           trace = 0)
  {
    # Arguments:
    #   see 'uniroot'

    # Value:
    #   Returns the x value of f where the root is located. If
    #   no root exists,  NA will be returned instead. In that case,
    #   the function doesn't terminate with an error  as
    #   the standard function uniroot().

    # Example:
    #   .unirootNA(sin, c(1, 2)); .unirootNA(sin, c(-1, 1))

    # If there is no Root:
    if (is.na(f.lower) || is.na(f.upper) || f.lower * f.upper > 0)
      return(NA)
    ## else there is one :
    stats::uniroot(
      f,
      interval = interval,
      ...,
      lower = lower,
      upper = upper,
      f.lower = f.lower,
      f.upper = f.upper,
      extendInt = extendInt,
      check.conv = check.conv,
      tol = tol,
      maxiter = maxiter,
      trace = trace
    )$root
  }
