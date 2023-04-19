imagN <- (0 + (0 + (0 + (0+1i))))

####Tempered stable subordinator (TSS)####

#' Characteristic function of the tempered stable subordinator
#'
#' Theoretical characteristic function (CF) of the distribution of the tempered
#' stable subordinator. See Kawai & Masuda (2011) for details.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, delta, lambda)}.
#' Either provide the parameters \code{alpha}, \code{delta}, \code{lambda}
#' individually OR provide \code{theta}.
#' \deqn{\varphi_{TSS}(t;\theta):=E_{\theta}\left[
#' \mathrm{e}^{\mathrm{i}tY}\right]= \exp\left(\delta\Gamma(-\alpha)
#' \left((\lambda-\mathrm{i}t)^{\alpha}-\lambda^{\alpha}\right)\right)}
#'
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambda Tempering parameter. A real number > 0.
#' @param theta Parameters stacked as a vector.
#'
#' @return The CF of the tempered stable subordinator distribution.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kawai, R. & Masuda, H. (2011), 'On simulation of tempered stable random
#' variates' \doi{10.1016/j.cam.2010.12.014}
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distributions and
#' processes' \doi{10.1016/j.spa.2013.06.012}
#'
#' @examples
#' x <- seq(-10,10,0.25)
#' y <- charTSS(x,0.5,1,1)
#'
#' @export
charTSS <- function(t, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL){
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
#' The probability density function (PDF) of tempered stable subordinator distribution.
#' It can be computed via the stable distribution (see details)
#' using the \code{stabledist} package.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, delta, lambda)}. Either provide the parameters
#' \code{alpha}, \code{delta}, \code{lambda} individually OR provide \code{theta}.
#' \deqn{f_{TSS}(y;\theta)=\mathrm{e}^{-\lambda y-\lambda^{\alpha}\delta\Gamma(-\alpha)}f_{S(\alpha,\delta)}(y),}
#' where \deqn{f_{S(\alpha,\delta)}} is the density of the stable subordinator.
#'
#' @param x A numeric vector of positive quantiles.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambda Tempering parameter. A real number > 0.
#' @param theta Parameters stacked as a vector.
#'
#' @return As \code{x} is a numeric vector, the return value is also a numeric
#' vector of probability densities.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kawai, R. & Masuda, H. (2011), 'On simulation of tempered stable random
#' variates' \doi{10.1016/j.cam.2010.12.014}
#'
#' @examples
#' x <- seq(0,15,0.25)
#' y <- dTSS(x,0.5,1,0.3)
#' plot(x,y)
#'
#' @importFrom stabledist dstable
#' @export
dTSS <- function(x, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL) {
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

#' Cumulative probability distribution function of the tempered stable subordinator
#' distribution
#'
#' The cumulative probability distribution function (CDF) of the tempered
#' stable subordinator distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, delta, lambda)}. Either provide the parameters
#' \code{alpha}, \code{delta}, \code{lambda} individually OR provide \code{theta}.
#' The function integrates the PDF numerically with \code{integrate()}.
#'
#' @param q A numeric vector of positive quantiles.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambda Tempering parameter. A real number > 0.
#' @param theta Parameters stacked as a vector.
#' @param pmethod A string. If not "integrate", the function \code{chartocdf()}
#' will be triggered.
#' @param N is a power of two & N >= 1024. if
#' \code{pmethod != "integrate"}. 8192 by default. Relevant for
#' @param ... Possibility to modify \code{stats::integrate()}.
#'
#' @return  As \code{q} is a numeric vector, the return value is also a numeric
#' vector of probabilities.
#'
#' @seealso
#' See also the [dTSS()] density-function.
#'
#' @examples
#' \donttest{
#' x <- seq(0,15,0.5)
#' y <- pTSS(x,0.7,1.354,0.3)
#' plot(x,y)
#' }
#'
#' @export
pTSS <- function(q, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL,
                   pmethod = "integrate", N = 8192, ...) {
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
    p <- sapply(q, function(z, ...) min(stats::integrate(dTSS, 0, z,
                                                         alpha = alpha,
                                                         delta = delta,
                                                         lambda = lambda,
                                                         ...)$value, 1 - 1e-07))

  } else {
    m <- gamma(1 - alpha)*delta/lambda^(1-alpha)
    s <- sqrt(gamma(2-alpha)*delta/lambda^(2-alpha))
    p <- sapply(q, chartocdf, N = N, m = m, s = s, char = charTSS,
                alpha = alpha, delta = delta, lambda = lambda)
  }

  return(p)
}


#' Function to generate random variates of the TSS distribution.
#'
#' Generates \code{n} random numbers distributed according
#' of the tempered stable subordinator distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, delta, lambda)}. Either provide the parameters
#' \code{alpha}, \code{delta}, \code{lambda} individually OR provide \code{theta}.
#' "AR" stands for the Acceptance-Rejection Method and "SR" for a truncated infinite shot
#' noise series representation. "AR" is the standard method used.
#' For more details, see references.
#'
#' @param n sample size (integer).
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambda Tempering parameter. A real number > 0.
#' @param theta Parameters stacked as a vector.
#' @param methodR A String. Either "AR" or "SR".
#' @param k integer: the level of truncation, if \code{methodR == "SR"}. 10000
#' by default.
#'
#' @return Generates \code{n} random numbers.
#'
#' @examples
#' rTSS(100,0.5,1,1)
#' rTSS(100,0.5,1,1,NULL,"SR",50)
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kawai, R & Masuda, H (2011), 'On simulation of tempered stable random
#' variates' \doi{10.1016/j.cam.2010.12.014}
#'
#' @export
rTSS <- function(n, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL,
                   methodR = "AR", k = 10000) {
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
    x <- switch(methodR,
                AR = rTSS_AR(n = n, alpha = alpha, delta = delta,
                             lambda = lambda),
                SR = rTSS_SR(n = n, alpha = alpha, delta = delta,
                             lambda = lambda, k = k))
    return(x)
}

# No export.
#' @importFrom stabledist rstable
rTSS_AR <- function(n, alpha, delta, lambda) {
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





##Tests
rTSS_MM <- function(n, alpha, delta, lambda){

  sigma <- ((1/alpha * gamma(1 - alpha) * delta *
               cos(alpha * pi/2))^(1/alpha))

  #copula::retstableR und copula::retstable(C) geben unterschiedliche Werte, obwohl beide das
  # gleiche ausgeben sollten und für rTSS_MM(n,0.5,1,1) hier sogar die gleiche Methode benutzen sollten ...

  # Für copula::rstable1 und stabledist::rstable kommen die gleichen Werte raus. Der Unterschied in den Methoden ist,
  # dass in rTSS_AR sigma(in der Funktion ist es delta) ein fester Wert ist, während in rTSS_MM delta == .gamm ist,
  # welcher zum Teil aus V0 besteht und ein random-Vektor ist


  #Ansatz MM: Alles direkt in der neune Funktion eingeben
  # --> Funktioniert noch nicht
  #returnVector <- copula::retstableR(alpha = alpha,
  #                                  V0 = (stats::runif(n, 0, 1)/cos(alpha*pi/2))^(alpha),
  #                                  h = sigma)*8



  #Ansatz retstable ohne doppelte Zufallszahlen:
  # Performance rausholen, indem man sigma zum n-langen Vektor macht
  # Funktioniert perfekt für c(0.5,1,3)
  returnVector <- c()

  for(i in 1:n){
    returnVector <- append(returnVector, copula::retstable(alpha,sigma) /
                             (gamma(1-alpha)*delta*lambda^(alpha-1))/lambda)
  }

  return(returnVector)

  #Ansatz Till:
  returnVector <- c()

  V0 <- -delta*gamma(-alpha)

  for(i in 1:n){
    returnVector <- append(returnVector, copula::retstable(alpha,V0,lambda))
  }

  return(returnVector)










  #Ansatz Test
  # returnVector <- c()
  #
  # for(i in 1:n){
  #   returnVector <- append(returnVector, copula::retstable(alpha,sigma) /
  #                            log(-lambda-delta*(alpha+1)*gamma(-alpha)*lambda^alpha))
  # }
  #
  # return(returnVector)


  # Austauschen von rstable in der ursprünglichen Funktion. Code zum Testen:
  # plot(x, dTSS(x,0.5,1,3), col = "red")
  # lines(density(rTSS(2000,0.5,1,3))) ~23Sekunden
  # lines(density(rTSS_MM(2000,0.5,1,3))) ~21Sekunden + copula::rstable1 liegt eher auf der Dichtefunktion
  # returnVector <- c()
  #
  # i <- 0
  # while (i < n) {
  #   U <- 2
  #   V <- 0
  #
  #   while (U > exp(-min(c(lambda * V, 700)))) {
  #     U <- stats::runif(1, 0, 1)
  #     V <- copula::rstable1(1, alpha, 1, sigma, 0, pm = 1)
  #   }
  #
  #   returnVector <- append(returnVector, V)
  #   i <- i + 1
  # }
  #
  # return(returnVector)

  #Ansatz gamm. rauskürzen
  returnVector <- c()

  V0 <- (sigma^(alpha)*m.opt.retst(sigma))/cospi2(alpha)

  returnVector <- c()

  for(i in 1:n){
    returnVector <- append(returnVector, retstableR(alpha,sigma,lambda = lambda))
  }

  return(returnVector)

}





# No export.
rTSS_SR <- function(n, alpha, delta, lambda, k) {
    base::replicate(n = n, rTSS_SR1(alpha = alpha, delta = delta,
                                    lambda = lambda, k = k))

}

# No export.
rTSS_SR1 <- function(alpha, delta, lambda, k) {
    parrivalslong <- cumsum(stats::rexp(k * 1.1))
    parrivals <- parrivalslong[parrivalslong <= k]
    E1 <- stats::rexp(length(parrivals))
    U <- stats::runif(length(parrivals))
    X <- cbind((alpha * parrivals/delta)^(-1/alpha), E1 * U^(1/alpha)/lambda)
    return(sum(apply(X, 1, FUN = min)))
}

#' Quantile function of the tempered stable subordinator distribution
#'
#' The quantile function of the tempered stable
#' subordinator distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, delta, lambda)}.
#' Either provide the parameters \code{alpha}, \code{delta}, \code{lambda}
#' individually OR provide \code{theta}. The function searches for a root
#' between \code{qmin} and \code{qmax} with \code{uniroot}. Boundaries can
#' either be supplied by the user or a built-in approach using the stable
#' distribution is used.
#'
#' @param p A numeric vector of probabilities. Each probability must be a real
#' number >0 and <1.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambda Tempering parameter. A real number > 0.
#' @param theta Parameters stacked as a vector.
#' @param qmin,qmax Limits of the interval. Will be computed if
#' \code{==NULL}.
#' @param ... Modify [pTSS()] and [stats::uniroot()].
#'
#' @return  As \code{p} is a numeric vector, the return value is also a numeric
#' vector of quantiles.
#'
#' @seealso
#' See also the [pTSS()] probability function.
#'
#' @examples
#' \donttest{
#' qTSS(0.5,0.5,5,0.01)
#' qTSS(0.5,0.9,1,10,NULL)
#' }
#'
#' @importFrom stabledist qstable
#' @export
qTSS <- function(p, alpha = NULL, delta = NULL, lambda = NULL, theta = NULL,
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
      pTSS(q = x, alpha = alpha, delta = delta, lambda = lambda, ...) - y
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
#' Theoretical characteristic function (CF) of the classical tempered
#' stable distribution. See Kuechler & Tappe (2013) for details.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, deltap, deltam,
#' lambdap, lambdam, mu)}. Either provide the parameters individually OR
#' provide \code{theta}.
#' \deqn{\varphi_{CTS}(t;\theta):=
#' E_{\theta}\left[
#' \mathrm{e}^{\mathrm{i}tX}\right]=
#' \exp\left(\mathrm{i}t\mu+\delta_+\Gamma(-\alpha)
#' \left((\lambda_+-\mathrm{i}t)^{\alpha}-\lambda_+^{\alpha}+
#' \mathrm{i}t\alpha\lambda_+^{\alpha-1}\right)\right.\\}
#' \deqn{\left. +\delta_-\Gamma(-\alpha)
#' \left((\lambda_-+\mathrm{i}t)^{\alpha}-\lambda_-^{\alpha}-\mathrm{i}t\alpha
#' \lambda_-^{\alpha-1}\right)
#' \right)}
#'
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alpha Stability parameter. A real number between 0 and 2.
#' @param deltap Scale parameter for the right tail. A real number > 0.
#' @param deltam  Scale parameter for the left tail. A real number > 0.
#' @param lambdap Tempering parameter for the right tail. A real number > 0.
#' @param lambdam Tempering parameter for the left tail. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#'
#' @return The CF of the tempered stable subordinator distribution.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kuechler, U. & Tappe, S. (2013), 'Tempered stable distributions and
#' processes' \doi{10.1016/j.spa.2013.06.012}
#'
#' @examples
#' x <- seq(-10,10,0.25)
#' y <- charCTS(x,1.5,1,1,1,1,0)
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
#' The probability density function (PDF) of the classical tempered stable
#' distributions is not available in closed form.
#' Relies on fast Fourier transform (FFT) applied to the characteristic
#' function.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, deltap, deltam,
#' lambdap, lambdam, mu)}. Either provide the parameters individually OR
#' provide \code{theta}. Methods include the FFT or alternatively by convolving
#' two totally positively skewed tempered stable distributions, see Massing
#' (2022).
#'
#' The "FFT" method is automatically selected for Mac users, as the "Conv"
#' method causes problems.
#'
#' @param x A numeric vector of quantiles.
#' @param alpha Stability parameter. A real number between 0 and 2.
#' @param deltap Scale parameter for the right tail. A real number > 0.
#' @param deltam  Scale parameter for the left tail. A real number > 0.
#' @param lambdap Tempering parameter for the right tail. A real number > 0.
#' @param lambdam Tempering parameter for the left tail. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#' @param dens_method Algorithm for numerical evaluation. Choose between \code{
#' "FFT"} (default) and \code{"Conv"}.
#' @param a Starting point of FFT, if \code{dens_method == "FFT"}. -20
#' by default.
#' @param b Ending point of FFT, if \code{dens_method == "FFT"}. 20
#' by default.
#' @param nf Pieces the transformation is divided in. Limited to power-of-two
#' size. 2048 by default.
#'
#' @return As \code{x} is a numeric vector, the return value is also a numeric
#' vector of densities.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' @examples
#' x <- seq(0,15,0.25)
#' y <- dCTS(x,0.6,1,1,1,1,1,NULL,"FFT",-20,20,2048)
#' plot(x,y)
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

    if (dens_method == "FFT" || .Platform$OS.type != "windows") {
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
              log(stabledist::dstable((y - gamma(1 - alpha) * deltap *
                                         lambdap^(alpha - 1)),
                          alpha, 1, Sigmap, 0, pm = 1) *
                    stabledist::dstable((y - gamma(1 - alpha) * deltam *
                               lambdam^(alpha - 1) + z),
                            alpha, 1, Sigmam, 0, pm = 1)))

    }

    return(stats::integrate(f = integrandDTSclass, lower = -Inf, upper = Inf,
                            abs.tol = 0)$value)

}

#' Cumulative probability function of the classic tempered stable (CTS)
#' distribution
#'
#' The cumulative probability distribution function (CDF) of the classic
#' tempered stable distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, deltap, deltam,
#' lambdap, lambdam, mu)}. Either provide the parameters individually OR
#' provide \code{theta}.
#' The function integrates the PDF numerically with \code{integrate()}.
#'
#' @param q A numeric vector of quantiles.
#' @param alpha Stability parameter. A real number between 0 and 2.
#' @param deltap Scale parameter for the right tail. A real number > 0.
#' @param deltam  Scale parameter for the left tail. A real number > 0.
#' @param lambdap Tempering parameter for the right tail. A real number > 0.
#' @param lambdam Tempering parameter for the left tail. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#' @param a Starting point of FFT, if \code{dens_method == "FFT"}. -20
#' by default.
#' @param b Ending point of FFT, if \code{dens_method == "FFT"}. 20
#' by default.
#' @param nf Pieces the transformation is divided in. Limited to power-of-two
#' size.
#' @param ... Possibility to modify \code{stats::integrate()}.
#'
#' @return As \code{q} is a numeric vector, the return value is also a numeric
#' vector of probabilities.
#'
#' @seealso
#' See also the [dCTS()] density-function.
#'
#' @examples
#' \donttest{
#' x <- seq(-5,5,0.25)
#' y <- pCTS(x,0.5,1,1,1,1,1)
#' plot(x,y)
#' }
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
                  else min(stats::integrate(dCTS, lower = a, upper = z,
                                            alpha = alpha, deltap = deltap,
                                            deltam = deltam, lambdap = lambdap,
                                            lambdam = lambdam, mu = mu,...)
                           $value, 1 - 1e-07)})

    return(p)
}

#' Function to generate random variates of CTS distribution.
#'
#' Generates \code{n} random numbers distributed according to the classic
#' tempered stable (CTS) distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, deltap, deltam,
#' lambdap, lambdam, mu)}. Either provide the parameters individually OR
#' provide \code{theta}.
#' "AR" stands for the approximate Acceptance-Rejection Method and "SR" for a
#' truncated infinite shot noise series representation. "AR" is the standard
#' method used.
#' For more details, see references.
#'
#' @param n sample size (integer).
#' @param alpha Stability parameter. A real number between 0 and 2.
#' @param deltap Scale parameter for the right tail. A real number > 0.
#' @param deltam  Scale parameter for the left tail. A real number > 0.
#' @param lambdap Tempering parameter for the right tail. A real number > 0.
#' @param lambdam Tempering parameter for the left tail. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#' @param methodR A String. Either "AR" or "SR".
#' @param k integer: the level of truncation, if \code{methodR == "SR"}. 10000
#' by default.
#' @param c A real number. Only relevant for \code{methodR == "AR"}.
#' 1 by default.
#'
#' @return Generates \code{n} random numbers.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kawai, R & Masuda, H (2011), 'On simulation of tempered stable random
#' variates' \doi{10.1016/j.cam.2010.12.014}
#'
#' @examples
#' rCTS(10,0.5,1,1,1,1,1,NULL,"SR",10)
#' rCTS(10,0.5,1,1,1,1,1,NULL,"aAR")
#'
#' @export
rCTS <- function(n, alpha = NULL, deltap = NULL, deltam = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, methodR = "AR",
                 k = 10000, c = 1) {
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

    x <- switch(methodR,
                AR = rCTS_aAR(n = n, alpha = alpha, deltap = deltap,
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

#' Quantile function of the classic tempered stable (CTS)
#'
#' The quantile function of the classic tempered stable (CTS).
#'
#' \code{theta} denotes the parameter vector \code{(alpha, deltap, deltam,
#' lambdap, lambdam, mu)}. Either provide the parameters individually OR
#' provide \code{theta}.
#' The function searches for a root between \code{qmin} and \code{qmax} with
#' \code{uniroot}. Boundaries can either be supplied by the user or a built-in
#' approach using the stable distribution is used.
#'
#' @param p A numeric vector of probabilities. Each probability must be a real
#' number >0 and <1.
#' @param alpha Stability parameter. A real number between 0 and 2.
#' @param deltap Scale parameter for the right tail. A real number > 0.
#' @param deltam  Scale parameter for the left tail. A real number > 0.
#' @param lambdap Tempering parameter for the right tail. A real number > 0.
#' @param lambdam Tempering parameter for the left tail. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#' @param qmin,qmax Limits of the interval. Will be computed if
#' \code{==NULL}.
#' @param ... Modify [pTSS()] and [stats::uniroot()].
#'
#' @return  As \code{p} is a numeric vector, the return value is also a numeric
#' vector of quantiles.
#'
#' @seealso
#' See also the [pCTS()] probability function.
#'
#' @examples
#' \donttest{
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
#' Theoretical characteristic function (CF) of the normal tempered
#' stable distribution.
#' See Rachev et al. (2011) for details.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, beta, delta, lambda,
#' mu)}. Either provide the parameters individually OR provide \code{theta}.
#' \deqn{\varphi_{NTS}(t;\theta)=E\left[\mathrm{e}^{\mathrm{i}tZ}\right]= \exp
#' \left(\mathrm{i}t\mu+\delta\Gamma(-\alpha)\left((\lambda-\mathrm{i}t
#' \beta+t^2/2)^{\alpha}-\lambda^{\alpha}\right)\right)
#' }
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param beta Skewness parameter. Any real number.
#' @param delta Scale parameter. A real number > 0.
#' @param lambda Tempering parameter. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#'
#' @return The CF of the normal tempered stable distribution.
#'
#' @references
#' Massing, T. (2022), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Rachev, S., Kim, Y., Bianchi, M. & Fabozzi, F. (2011), 'Financial Models with
#' Levy Processes and Volatility Clustering' \doi{10.1002/9781118268070}
#'
#' @examples
#' x <- seq(-10,10,0.25)
#' y <- charNTS(x,0.5,1,1,1,0)
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
#' The probability density function (PDF) of the normal tempered stable
#' distributions is not available in closed form.
#' Relies on fast Fourier transform (FFT) applied to the characteristic
#' function.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, beta, delta, lambda,
#' mu)}. Either provide the parameters individually OR provide \code{theta}.
#' Currently, the only method is FFT.
#'
#' @param x A numeric vector of quantile.
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
#' vector of densities.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' @examples
#' x <- seq(0,15,0.25)
#' y <- dNTS(x,0.8,1,1,1,1)
#' plot(x,y)
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
    return (as.numeric(stats::approx(xgrid, densityW, xout=x,
                              yleft = 1e-18, yright = 1e-18)[2]))
}

#' Cumulative probability function of the normal tempered stable (NTS)
#' distribution
#'
#' The cumulative probability distribution function (CDF) of the normal
#' tempered stable distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, beta, delta, lambda,
#' mu)}. Either provide the parameters individually OR provide \code{theta}.
#' The function integrates the PDF numerically with \code{integrate()}.
#'
#' @param q A numeric vector of quantile.
#' @param alpha A real number between 0 and 1.
#' @param beta Any real number.
#' @param delta A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param a Starting point integrate density function. -40 by default.
#' @param b Ending point of integrate density function. 40 by default.
#' @param nf Pieces the fast Fourier transformation is divided in. Limited to
#' power-of-two size. 2^11 by default.
#' @param ... Change parameters in [dNTS()]
#'
#' @return  As \code{q} is a numeric vector, the return value is also a numeric
#' vector of probabilities.
#'
#' @seealso
#' See also the [dNTS()] density-function.
#'
#' @examples
#' \donttest{
#' x <- seq(-5,5,0.25)
#' y <- pNTS(x,0.5,1,1,1,1)
#' plot(x,y)
#' }
#'
#' @export
pNTS <- function(q, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, a = -40, b = 40, nf = 2^11, ...) {
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
                else min(stats::integrate(dNTS, lower = a, upper = z,
                                          alpha = alpha, beta = beta,
                                          delta = delta, lambda = lambda,
                                          mu = mu, ...)$value, 1-1e-7)})

    return(p)
}

#' Function to generate random variates of NTS distribution.
#'
#' Generates \code{n} random numbers distributed according
#' of the normal tempered stable distribution.
#'
#' \code{theta} denotes the parameter vector \code{(alpha, beta, delta, lambda,
#' mu)}. Either provide the parameters individually OR provide \code{theta}.
#' Works by a normal variance-mean mixture with a TSS distribution. Method
#' parameter is for the method of simulating the TSS random variable, see the
#' [rTSS()] function.
#' "AR" stands for the Acceptance-Rejection Method and "SR" for a truncated
#' infinite shot noise series representation. "AR" is the standard method used.
#'
#' For more details, see references.
#'
#' @param n sample size (integer).
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder.
#' @param delta  A real number > 0.
#' @param lambda A  real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param methodR A String. Either "AR" or "SR". "AR" by default.
#' @param k integer: the number of replications, if \code{methodR == "SR"}. 10000
#' by default.
#'
#' @return Generates \code{n} random numbers.
#'
#' @seealso
#' See also the [rTSS()] function.
#'
#' @references
#' Massing, T. (2023), 'Parametric Estimation of Tempered Stable Laws'
#'
#' Kawai, R & Masuda, H (2011), 'On simulation of tempered stable random
#' variates' \doi{10.1016/j.cam.2010.12.014}
#'
#' @examples
#' rNTS(100, 0.5, 1,1,1,1)
#' rNTS(10, 0.6, 0,1,1,0)
#' rNTS(10, 0.5, 1,1,1,1, NULL, "SR", 100)
#'
#' @export
rNTS <- function(n, alpha = NULL, beta = NULL, delta = NULL, lambda = NULL,
                 mu = NULL, theta = NULL, methodR = "AR", k = 10000) {
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
    x <- switch(methodR,
                AR = rNTS_AR(n = n, alpha = alpha, beta = beta,
                             delta = delta, lambda = lambda, mu = mu),
                SR = rNTS_SR(n = n, alpha = alpha, beta = beta, delta = delta,
                             lambda = lambda, mu = mu, k = k))
    return(x)
}

# No export.
rNTS_AR <- function(n, alpha, beta, delta, lambda, mu) {
    z <- stats::rnorm(n = n)
    y <- rTSS_AR(n = n, alpha = alpha, delta = delta, lambda = lambda)
    x <- sqrt(y) * z + beta * y + mu
    return(x)
}


# No export.
rNTS_SR <- function(n, alpha, beta, delta, lambda, mu, k) {
    z <- stats::rnorm(n = n)
    y <- rTSS_SR(n = n, alpha = alpha, delta = delta, lambda = lambda, k = k)
    x <- sqrt(y) * z + beta * y + mu
    return(x)
}

#' Quantile function of the normal tempered stable (NTS)
#'
#' The quantile function of the normal tempered stable (CTS).
#'
#' \code{theta} denotes the parameter vector \code{(alpha, beta, delta, lambda,
#' mu)}. Either provide the parameters individually OR provide \code{theta}.
#' The function searches for a root between \code{qmin} and \code{qmax} with
#' \code{uniroot}.
#' Boundaries can either be supplied by the user or a built-in approach using
#' the stable distribution is used.
#'
#' @param p A numeric vector of probabilities. Each probability must be a real
#' number >0 and <1.
#' @param alpha A real number between 0 and 1.
#' @param beta A gap holder.
#' @param delta  A real number > 0.
#' @param lambda A  real number >= 0.
#' @param mu A location parameter, any real number.
#' @param theta A vector of all other arguments.
#' @param qmin,qmax Limits of the interval. Will be computed if
#' \code{==NULL}.
#' @param ... Modify [pNTS()] and [stats::uniroot()].
#'
#' @return As \code{p} is a numeric vector, the return value is also a numeric
#' vector of quantiles.
#'
#' @seealso
#' See also the [pNTS()] probability function.
#'
#' @examples
#' \donttest{
#' qNTS(0.1,0.5,1,1,1,1)
#' qNTS(0.3,0.6,1,1,1,1,NULL)
#' }
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
rCGMY <- function(n, C, G, M, Y, methodR = "SR", k = 100, ...) {
    rCTS(n = n, alpha = Y, deltap = C, deltam = C, lambdap = G, lambdam = M,
         mu = 0, methodR = methodR, k = k, ...)
}


#####Bohman Char to CDF#####

#Cumulative distribution function from characteristic function
# No export.
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

