imagN <- (0 + (0 + (0 + (0+1i))))

####Subordinator Tempered Stable (STS)####

#' Function title.
#'
#' When \code{alpha} is a real number between 0 and 1, the tempered stable
#' distribution TS'(\code{alpha}, \code{delta}, \code{lambda}) can be simulated
#' exactly through acceptance-rejection sampling.
#'
#' For whole derivation and meaning of single letters see also Kawai et. Masuda
#' (2011).
#' Basically, the derivation starts from a one-dimensional stable distribution
#' which is later added with a centered and totally positively skewed tempered
#' stable (Levy) process.
#'
#' Reiichiro Kawai, & Hiroki Masuda (2011).
#' On simulation of tempered stable random variates.
#' Journal of Computational and Applied Mathematics, 235(8), 2873-2887.
#' https://doi.org/10.1016/j.cam.2010.12.014
#' @seealso(\url{https://www.sciencedirect.com/science/article/pii/S0377042710006643})
#'
#' @param t A positive integer.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#'
#' @return The result of spectral positive tempered stable process
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

#' PDF TS subordinator by Stable subordinator
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param x A positive integer.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#'
#' @return As \code{x} is a numeric vector, the return value is also a numeric
#' vector
#'
#' @examples
#' dSTS(1000,0.5,1,0.3)
#' dSTS(1,0.5,1,0.3) <-- Warnings will show up.
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

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param q A numeric vector of quantiles.
#' @param alpha A real number between 0 and 1.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param theta A vector of all other arguments.
#' @param pmethod A gap holder.
#' @param N integer: the number of replications, if \code{pmethod != "integrate"}. 10000
#' by default.
#'
#' @return A gap holder.
#'
#' @examples
#' pSTS(3,0.7,1.354,0.3) <-- Warnings will show up.
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
    p <- sapply(q, function(z, ...) min(integrate(dSTS, 0, z, alpha = alpha,
                                         delta = delta, lambda = lambda, ...)$value,
                                   1 - 1e-07))

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

#' No export.
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
            U <- runif(1, 0, 1)
            V <- rstable(1, alpha, 1, sigma, 0, pm = 1)
        }

        returnVector <- append(returnVector, V)
        i <- i + 1
    }

    return(returnVector)

}

#' No export.
rSTS_SR <- function(n, alpha, delta, lambda, k) {
    base::replicate(n = n, rSTS_SR1(alpha = alpha, delta = delta,
                                    lambda = lambda, k = k))

}

#' No export.
rSTS_SR1 <- function(alpha, delta, lambda, k) {
    parrivalslong <- cumsum(rexp(k * 1.1))
    parrivals <- parrivalslong[parrivalslong <= k]
    E1 <- rexp(length(parrivals))
    U <- runif(length(parrivals))
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
#'
#' @return Gap holder for return.
#'
#' @examples
#' pSTS(0.5,0.5,5,0.01) <-- Warnings will show up.
#' pSTS(1,0.9,1,10,NULL) <-- Warnings will show up.
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
                                        cos(alpha*pi/2))^(1/alpha)),100), delta = 0, pm = 1))
      } else {
        qmini <- qmin
        qmaxi <- qmax
      }
      stabledist:::.unirootNA(froot,interval = c(qmini,qmaxi),
                              extendInt = "yes", y = y)}

    q <- sapply(p, qroot, qmin = qmin, qmax = qmax)
    return(q)
}

##### Classical Tempered Stable (CTS)#####

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param t A gap holder.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#'
#' @return Gap holder for return.
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

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param x A gap holder.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdam A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param dens_method A gap holder.
#' @param a A gap holder, if \code{dens_method == "FFT"}. -20
#' by default.
#' @param b A gap holder, if \code{dens_method == "FFT"}. 20
#' by default.
#' @param nf A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' dCTS(1,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048)
#' dCTS(1,0.6,1,1,1,1,1,NULL,"Conv") <-- Warnings will show up.
#'
#' @export
dCTS <- function(x, alpha = NULL, deltap = NULL, deltam = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, dens_method = "FFT",
                 a = -20, b = 20, nf = 2048, ...) {
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


#' No export.
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
    de

    densityW <- Re((dt/(2 * pi)) * exp(-imagN * (nf/2 * dt) * xgrid) * y)

    return (as.numeric(approx(xgrid, densityW, xout=x,
                              yleft = 1e-18, yright = 1e-18)[2]))
}

#' No export.
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
#' @param lambdap,lambdap A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param a A gap holder. -20 by default.
#' @param b A gap holder. 20 by default.
#' @param nf A gap holder.
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
                  else min(integrate(dCTS, lower = a, upper = z, alpha = alpha,
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
#' @param lambdap,lambdap A  real number > 0.
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
                 lambdam = NULL, mu = NULL, theta = NULL, method = "SR",
                 k = 100, c = 0.1) {
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

#' No export.
rCTS_aAR <- function(n, alpha, deltap, deltam, lambdap, lambdam, mu, c) {
    mu + rCTS_aARp(n, alpha = alpha, delta = deltap, lambda = lambdap, c = c)
    -rCTS_aARp(n, alpha = alpha, delta = deltam, lambda = lambdam, c = c)
}

#' No export.
#' @importFrom stabledist rstable
rCTS_aARp <- function(n, alpha, delta, lambda, c) {
    sigma <- ((1/alpha * gamma(1 - alpha) * delta *
                 cos(alpha * pi/2))^(1/alpha))

    returnVector <- c()

    i <- 0
    while (i < n) {
        U <- 2
        V <- 0
        Y <- 0

        while (U > exp(-min(c(lambda * (V + c), 700)))) {
            U <- runif(1, 0, 1)
            V <- rstable(1, alpha, 1, sigma, 0, pm = 1)
            Y <- V - gamma(1 - alpha) * delta * lambda^(alpha - 1)
        }

        returnVector <- append(returnVector, Y)
        i <- i + 1
    }

    return(returnVector)

}

#' No export.
rCTS_SR <- function(n, alpha, deltap, deltam, lambdap, lambdam, mu, k) {
    replicate(n = n, rCTS_SR1(alpha = alpha, deltap = deltap, deltam = deltam,
                              lambdap = lambdap, lambdam = lambdam,
                              mu = mu, k = k))

}

#' No export.
rCTS_SR1 <- function(alpha, deltap, deltam, lambdap, lambdam, mu, k) {
    x <- rCTS_SRp(alpha = alpha, delta = deltap, lambda = lambdap, k = k)
    -rCTS_SRp(alpha = alpha, delta = deltam, lambda = lambdam, k = k)
    +mu
    return(x)
}


#' No export.
#' @importFrom VGAM zeta
rCTS_SRp <- function(alpha, delta, lambda, k) {
    parrivalslong <- cumsum(rexp(k * 1.5))
    parrivals <- parrivalslong[parrivalslong <= k]
    n <- length(parrivals)
    E1 <- rexp(n)
    U <- stats::runif(n)
    cntr <- ((1:n) * alpha/delta)^(-1/alpha)
    gam <- (delta/alpha)^(1/alpha) * VGAM::zeta(1/alpha)
    -gamma(1 - alpha) * delta * lambda^(alpha - 1)
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
#' @param q A gap holder.
#' @param alpha A real number between 0 and 2.
#' @param deltap,deltam  A real number > 0.
#' @param lambdap,lambdap A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param qmin,qmax Limits of the interval.
#'
#' @return Gap holder for return.
#'
#' @examples
#' qCTS(0.5,1.5,10,10,10,10,10)
#' qCTS(0.5,1.5,1,1,1,1,1)
#'
#' @importFrom stabledist qstable
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
        x1 <- qcauchy(y, location = mu,
                      scale = min(((1/alpha*gamma(1-alpha)*(deltap+deltam)*
                                      cos(alpha*pi/2))^(1/alpha)),100))
        qmini <- min(c(x1,-x1, mu-1, mu+1))
        qmaxi <- max(c(x1,-x1, mu-1, mu+1))
      } else {
        qmini <- qmin
        qmaxi <- qmax
      }
      stabledist:::.unirootNA(froot, interval = c(qmini,qmaxi), extendInt = "yes",
                              y = y, ...)}

    q <- sapply(p, qroot, qmin = qmin, qmax = qmax)
    return(q)
}

##### Normal Tempered Stable#####

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param t A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#'
#' @return Gap holder for return.
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

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param x A gap holder.
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param dens_method Currently useless param, as it does nothing. "FFT" by
#' default.
#' @param a A gap holder. -20 by default.
#' @param b A gap holder. 20 by default.
#' @param nf A gap holder. 2048 by default.
#'
#' @return Gap holder for return.
#'
#' @examples
#' dNTS(1,0.8,1,1,1,1)
#' dNTS(0.5,0.5,20,20,20,20, a = -2000, b = 2000, nf = 8192)
#'
#' @export
dNTS <- function(x, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, dens_method = "FFT",
                 a = -20, b = 20, nf = 2048,...) {
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

#' @seealso{dNTS}
#' No export
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

    # The return value of approx (RStudio) and Interpolation (Wolfram
    # Mathematica) differ slightly. The larger nf, the smaller the difference.
    return (as.numeric(approx(xgrid, densityW, xout=x,
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
                else min(integrate(dNTS, lower = a, upper = z, alpha = alpha,
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

#' No export.
rNTS_AR <- function(n, alpha, beta, delta, lambda, mu) {
    z <- rnorm(n = n)
    y <- rSTS_AR(n = n, alpha = alpha, delta = delta, lambda = lambda)
    x <- sqrt(y) * z + beta * y + mu
    return(x)
}


#' No export.
rNTS_SR <- function(n, alpha, beta, delta, lambda, mu, k) {
    z <- rnorm(n = n)
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
        x1 <- qcauchy(y, location = mu,
                      scale = min(((1/alpha*gamma(1-alpha)*(delta)*
                                      cos(alpha*pi/2))^(1/alpha)),100))
        qmini <- min(c(x1,-x1, mu-1, mu+1))
        qmaxi <- max(c(x1,-x1, mu-1, mu+1))
      } else {
        qmini <- qmin
        qmaxi <- qmax
      }
      stabledist:::.unirootNA(froot,interval = c(qmini,qmaxi),
                              extendInt = "yes", y = y, ...)}

    q <- sapply(p, qroot, qmin = qmin, qmax = qmax)
    return(q)
}

##### CGMY#####

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param t A gap holder.
#' @param C  A real number > 0.
#' @param G,M A  real number > 0.
#' @param Y A real number between 0 and 2.
#'
#' @return Gap holder for return.
#'
#' @examples
#' charCGMY(1,1,1,1,0.5)
#' charCGMY(10,1,1,1,0.5)
#'
#' @export
charCGMY <- function(t, C, G, M, Y) {
    charCTS(t = t, alpha = Y, deltap = C, deltam = C, lambdap = G, lambdam = M,
            mu = 0)
}

#' Title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param x A gap holder.
#' @param C  A real number > 0.
#' @param G,M A  real number > 0.
#' @param Y A real number between 0 and 2.
#' @param dens_method A gap holder. "FFT" by default.
#' @param a A gap holder, if \code{dens_method == "FFT"}. -20
#' by default.
#' @param b A gap holder, if \code{dens_method == "FFT"}. 20
#' by default.
#' @param nf A gap holder.
#'
#' @return Gap holder for return.
#'
#' @examples
#' dCGMY(1,1,1,1,0.5)
#' dCGMY(1,1,1,1,0.5,"", -2, 2, 2^4)
#'
#' @export
dCGMY <- function(x, C, G, M, Y, dens_method = "FFT", a = -20, b = 20,
                  nf = 2048, ...) {
    dCTS(x = x, alpha = Y, deltap = C, deltam = C, lambdap = G, lambdam = M,
         mu = 0, dens_method = dens_method, a = a, b = b, nf = 2048, ...)
}

#' Function title
#'
#' Gap holder for description.
#'
#' Gap holder for details.
#'
#' @param n A gap holder.
#' @param Y A real number between 0 and 2.
#' @param C  A real number > 0.
#' @param G,M A  real number > 0.
#' @param method A String. Either "aAR" or "SR". "SR" by default.
#' @param k integer: the number of replications, if \code{method == "SR"}. 100
#' by default.
#'
#' @return Gap holder for return.
#'
#' @examples
#' rCGMY(100,1,1,1,0.5)
#' rCGMY(100,1,1,1,0.5, "SR", k= 100)
#'
#' @export
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


#' No export.
c_err <- function(t) {if(abs(t)<1){(1-t)*cos(pi*t)+1/pi*abs(sin(pi*t))} else 0}
