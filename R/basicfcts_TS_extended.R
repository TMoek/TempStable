
####Modified Tempered Stable Distribution (MTS)####

#' Characteristic function of the modified tempered stable distribution
#'
#' Theoretical characteristic function (CF) of the modified tempered stable
#' distribution. Since the parameterisation can be different for this
#' characteristic function in different approaches, the respective approach can
#' be selected with [functionOrigin]. For the estimation function
#' [TemperedEstim] and therefore also the Monte Carlo function
#' [TemperedEstim_Simulation] and the calculation of the density function
#' [dMTS], however, only the approach of Kim et al. (2008) or rachev11 can be selected. If you want to use the approach of kim09 for these functions, you have to clone the package from GitHub and adapt the functions accordingly.
#'
#' \strong{Origin of functions}
#' \describe{
#'   \item{kim09}{Ansatz aus: Kim et al. 2009 The modified tempered stable
#'   distribution, GARCH-models and option pricing. Alpha darf hier von -Inf
#'   bis 1 gehen. Ausgenommen wird alpha == 1/2.}
#'   \item{kim08}{Ansatz aus: Kim et al. 2008? Financial market models with
#'   Levy processes and time-varying volatility.}
#'   \item{rachev11}{ Ansatz aus: Rachev et al. 2011 Financial Models with Levy
#'   Processes and time-varying volatility. Rechnerisch gleich wie kim08, sieht
#'   nur etwas anders aus.
#'   }
#' }
#'
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambdap,lambdam Tempering parameter. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#' @param functionOrigin A string. Either "kim09", "rachev11" or "kim08".
#' Default is "kim08".
#'
#' @return The CF of the the modified tempered stable distribution.
#'
#' @references
#' Kim, Y. s.; Rachev, S. T.; Bianchi, M. L. & Fabozzi, F. J. (2008), 'Financial
#' market models with lévy processes and time-varying volatility'
#' \doi{10.1016/j.jbankfin.2007.11.004}
#'
#' @examples
#' gapholder
#'
#' @export
charMTS <- function(t, alpha = NULL, delta = NULL, lambdap = NULL,
                    lambdam = NULL, mu = NULL, theta = NULL,
                    functionOrigin = "kim08") {
  if ((missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
       missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
      missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }

  if (functionOrigin == "kim08" || functionOrigin == "rachev11"){
    stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 < lambdam)
  }
  else {
    #Function from: Kim et al. 2009 The modified tempered stable distribution ..
    stopifnot(0.5 != alpha, alpha < 1, 0 < delta, 0 < lambdap, 0 < lambdam)
  }

  #Ansatz aus: Rachev et al. 2011 Financial Models with Levy Processes...
  if(functionOrigin == "rachev11"){
    subfunctionR <- function(x, alpha, lambda){
      2^(-(alpha+3)/2) * sqrt(pi) * gamma(-alpha/2) *
        ((lambda^2+x^2)^(alpha/2) - lambda^alpha)
    }

    subfunctionI <- function(x, alpha, lambda){
      2^(-(alpha+1)/2) * gamma((1-alpha)/2) * lambda^(alpha-1) *
        (modifiedHyperGeoXr(1, (1 - alpha)/2, 3/2, -(x^2)/(lambdam^2)) -1)
    }

    return(exp(imagN * mu * t +
                 delta*(subfunctionR(t,alpha,lambdap) +
                          subfunctionR(t,alpha,lambdam)) +
                 imagN*t*delta*(subfunctionI(t,alpha,lambdap) -
                                  subfunctionI(t,alpha,lambdam))))
  }

  else if(functionOrigin == "kim08"){
    subfunctionR <- function(t, alpha, delta, lambdap, lambdam){
      sqrt(pi) * delta * gamma(-alpha/2) / (2^((alpha + 3) / 2))*
        ((lambdap^2 + t^2)^(alpha/2) - lambdap^alpha +
           (lambdam^2 + t^2)^(alpha/2) - lambdam^alpha)
    }

    subfunctionI <- function(t, alpha, delta, lambdap, lambdam){
      imagN * t * delta * gamma((1 - alpha) / 2) / (2^((alpha + 1) / 2)) *
        (lambdap^(alpha - 1) * modifiedHyperGeoXr(1, (1 - alpha)/2, 3/2,
                                                -(t^2)/(lambdap^2)) -
           lambdam^(alpha - 1) * modifiedHyperGeoXr(1, (1 - alpha)/2, 3/2,
                                                  -(t^2)/(lambdam^2))
        )
    }

    return(exp(imagN * mu * t + subfunctionR(t, alpha, delta, lambdap, lambdam) +
                 subfunctionI(t, alpha, delta, lambdap, lambdam)))
  }

  else{
    #Function from: Kim et al. 2009 The modified tempered stable distribution ..
    subfunctionR <- function(t, alpha, delta, lambdap, lambdam){
      if(alpha == 0){
        sqrt(pi)*2^(3/2)*delta*
          (log((lambdap^2)/((lambdap^2)+(t^2)))+
             log((lambdam^2)/((lambdam^2)+(t^2))))
      }
      else{
        sqrt(pi) * 2^(-alpha-(3/2)) * delta * gamma(-alpha) *
          (((lambdap^2) + (t^2))^alpha - lambdap^(2*alpha) +
             ((lambdam^2) + (t^2))^alpha - lambdam^(2*alpha))
      }
    }

    subfunctionI <- function(t, alpha, delta, lambdap, lambdam){
      (imagN * t * delta * gamma(1/2 - alpha)) / (2^(alpha + (1/2))) *
        (lambdap^(2*alpha-1) *
           modifiedHyperGeoXr(1, (1/2)-alpha, 3/2, - (t^2)/(lambdap^2)) -
           lambdam^(2*alpha-1) *
           modifiedHyperGeoXr(1, (1/2)-alpha, 3/2, - (t^2)/(lambdam^2))
        )
    }

    return(exp(imagN * mu * t + subfunctionR(t, alpha, delta, lambdap, lambdam)
               + subfunctionI(t, alpha, delta, lambdap, lambdam)))
  }
}


# xMTS <- seq(-0.25,0.25,0.005)
# yMTS <- dMTS(x, 0.7,0.005,50,50,0)
dMTS <- function(x, alpha = NULL, delta = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, dens_method = "FFT",
                 a = -20, b = 20, nf = 256) {
  if ((missing(alpha) | missing(delta) | missing(lambdap) |
       missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) |
      missing(lambdam) | missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 <
              lambdam)

  if (dens_method == "FFT" || .Platform$OS.type != "windows") {
    d <- sapply(x, d_FFT, charFunc = charMTS,
                theta = c(alpha, delta, lambdap, lambdam, mu),
                a = a, b = b, nf = nf)
  } else {
    d <- NULL
  }
  return(d)
}


pMTS <- function(q, alpha = NULL, delta = NULL, lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, dens_method = "FFT",
                 a = -40, b = 40, nf = 2048, ...) {
  if ((missing(alpha) | missing(delta) | missing(lambdap) |
       missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) |
      missing(lambdam) | missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 <
              lambdam)

  p <- numeric(length(q))

  p <- sapply(q,
              function(z) {
                if(z<a) 0
                else if (z>b) 1
                else min(stats::integrate(dMTS, lower = a, upper = z,
                                          alpha = alpha, delta = delta,
                                          lambdap = lambdap, lambdam = lambdam,
                                          mu = mu,...)
                         $value, 1 - 1e-07)})

  return(p)
}

rMTS <- function(n, alpha = NULL, delta = NULL, lambdap = NULL, lambdam = NULL,
                 mu = NULL, theta = NULL, methodR = "TM", k = 10000) {
  if ((missing(alpha) | missing(delta) | missing(lambdap) |
       missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) |
      missing(lambdam) | missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 <
              lambdam)
  #TODO
  x <- switch(methodR,
              AR = 0,
              SR = 0,
              TM = 0)
  return(x)
}

rMTS_SR <- function(n, alpha, delta, lambdap, lambdam, mu, k) {
  replicate(n = n,rMTS_SR_Ro(alpha = alpha, delta = delta,
                             lambdap = lambdap, lambdam = lambdam, k = k) + mu)
}

rMTS_SR_Ro <- function (alpha, delta, lambdap, lambdam, k){
  parrivalslong <- cumsum(stats::rexp(k * 1.1))
  parrivals <- parrivalslong[parrivalslong <=  k]
  E1 <- stats::rexp(length(parrivals))
  U <- stats::runif(length(parrivals))

  #Im Vergleich zu Rachev11 sind in sigma und damit auch in V + und - vertauscht.
  # Der rest ist gleich
  sigma <- 2^((alpha+1)/2) * delta * gamma(alpha/2+1/2)
  V <- rMTS_SR_rVj(length(parrivals), sigma, alpha, delta, lambdap, lambdam, k)
  b <- -2^(-(alpha+1)/2) * delta * gamma((1-alpha)/2) *
    (lambdap^(alpha-1)-lambdam^(alpha-1))
  X <- cbind((alpha * parrivals / sigma)^(-1/alpha),
             sqrt(2) * E1^(1/2) * U^(1/alpha)/abs(V))
  Xreturn <- sum((apply(X, 1, FUN = min)*V/abs(V)))+b

  return(Xreturn)
}

rMTS_SR_dVj <- function(x, sigma, alpha, delta, lambdap, lambdam){
  returnVec <- NULL
  for(xi in x){

    Ip <- 0
    Im <- 0

    if (xi>0){Ip <- 1}
    else if (xi<0){Im <- 1}

    if(xi == 0){y <- 0}
    else {
      y <- delta / sigma *
        (lambdap^(alpha-1) * exp(-(lambdap^2*xi^2)/2) * Ip +
           lambdam^(alpha-1) * exp(-(lambdam^2*xi^2)/2) * Im )
    }
    returnVec <- append(returnVec,y)
  }
  return(returnVec)
}

rMTS_SR_rVj <- function(n, sigma, alpha, delta, lambdap, lambdam, k){
  dX <- (20*2/k)
  x <- seq(-5,5,dX)
  y <- rMTS_SR_dVj(x, sigma, alpha, delta, lambdap, lambdam)
  cumY <- cumsum(y)
  rV <- runif(n, min(cumY), max(cumY))

  returnVector <- NULL
  for(s in rV){
    pos <- which.min(abs(cumY - s))
    returnVector <- append(returnVector, x[pos])
  }

  return(returnVector)
}

#### Generalized Classical Tempered Stable Distribution ####

#' Characteristic function of the generalized classical tempered stable
#' distribution
#'
#' gapholder
#'
#' gapholder
#'
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alphap,alpham Stability parameter. A real number between 0 and 1.
#' @param deltap,deltam Scale parameter. A real number > 0.
#' @param lambdap,lambdam Tempering parameter. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#'
#' @return The CF of the the generalized classical tempered stable distribution.
#'
#' @references
#' gapholder
#'
#' @examples
#' gapholder
#'
#' @export
charGTS <- function(t, alphap = NULL, alpham = NULL, deltap = NULL,
                    deltam = NULL, lambdap = NULL, lambdam = NULL, mu = NULL,
                    theta = NULL) {
  if ((missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
       missing(lambdap) | missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alphap, alpham, deltap, deltam, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
      missing(lambdap) | missing(lambdam) | missing(mu)) {
    alphap <- theta[1]
    alpham <- theta[2]
    deltap <- theta[3]
    deltam <- theta[4]
    lambdap <- theta[5]
    lambdam <- theta[6]
    mu <- theta[7]
  }
  stopifnot(0 < alphap, alphap < 2, 0 < alpham, alpham < 2, 0 < deltap,
            0 < deltam, 0 < lambdap, 0 < lambdam)

  return(exp(imagN * t * mu
             - imagN * t * gamma(1 - alphap)*(deltap * lambdap^(alphap-1))
             + imagN * t * gamma(1 - alpham)*(deltam * lambdam^(alpham-1))
             + deltap * gamma(-alphap) * ((lambdap - imagN * t)^(alphap) -
                                            lambdap^(alphap))
             + deltam * gamma(-alpham) * ((lambdam + imagN * t)^(alpham) -
                                            lambdam^(alpham))))
}

dGTS <- function(x, alphap = NULL, alpham = NULL, deltap = NULL,
                 deltam = NULL, lambdap = NULL, lambdam = NULL, mu = NULL,
                 theta = NULL, dens_method = "FFT",
                 a = -20, b = 20, nf = 2048) {
  if ((missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
       missing(lambdap) | missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alphap, alpham, deltap, deltam, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
      missing(lambdap) | missing(lambdam) | missing(mu)) {
    alphap <- theta[1]
    alpham <- theta[2]
    deltap <- theta[3]
    deltam <- theta[4]
    lambdap <- theta[5]
    lambdam <- theta[6]
    mu <- theta[7]
  }
  stopifnot(0 < alphap, alphap < 2, 0 < alpham, alpham < 2, 0 < deltap,
            0 < deltam, 0 < lambdap, 0 < lambdam)

  if (dens_method == "FFT" || .Platform$OS.type != "windows") {
    d <- sapply(x, d_FFT, charFunc = charGTS,
                theta = c(alphap, alpham, deltap, deltam, lambdap, lambdam, mu),
                a = a, b = b, nf = nf)
  } else {
    d <- NULL
  }
  return(d)
}


pGTS <- function(q, alphap = NULL, alpham = NULL, deltap = NULL,
                 deltam = NULL, lambdap = NULL, lambdam = NULL, mu = NULL,
                 theta = NULL, dens_method = "FFT",
                 a = -40, b = 40, nf = 2048, ...) {
  if ((missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
       missing(lambdap) | missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alphap, alpham, deltap, deltam, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
      missing(lambdap) | missing(lambdam) | missing(mu)) {
    alphap <- theta[1]
    alpham <- theta[2]
    deltap <- theta[3]
    deltam <- theta[4]
    lambdap <- theta[5]
    lambdam <- theta[6]
    mu <- theta[7]
  }
  stopifnot(0 < alphap, alphap < 2, 0 < alpham, alpham < 2, 0 < deltap,
            0 < deltam, 0 < lambdap, 0 < lambdam)

  p <- numeric(length(q))

  p <- sapply(q,
              function(z) {
                if(z<a) 0
                else if (z>b) 1
                else min(stats::integrate(dGTS, lower = a, upper = z,
                                          alphap = alphap, alpham = alpham,
                                          deltap = deltap, deltam = deltam,
                                          lambdap = lambdap, lambdam = lambdam,
                                          mu = mu, ...)
                         $value, 1 - 1e-07)})

  return(p)
}


rGTS <- function(n, alphap = NULL, alpham = NULL, deltap = NULL, deltam = NULL,
                 lambdap = NULL,
                 lambdam = NULL, mu = NULL, theta = NULL, methodR = "AR",
                 k = 10000, c = 1) {
  if ((missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
       missing(lambdap) | missing(lambdam) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alphap, alpham, deltap, deltam, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alphap) | missing(alpham) | missing(deltap) | missing(deltam) |
      missing(lambdap) | missing(lambdam) | missing(mu)) {
    alphap <- theta[1]
    alpham <- theta[2]
    deltap <- theta[3]
    deltam <- theta[4]
    lambdap <- theta[5]
    lambdam <- theta[6]
    mu <- theta[7]
  }
  stopifnot(0 < alphap, alphap < 2, 0 < alpham, alpham < 2, 0 < deltap,
            0 < deltam, 0 < lambdap, 0 < lambdam)

  #TODO Insert TM
  if(methodR == "TM"){
    methodR <- "AR"
  }

  x <- switch(methodR,
              AR = rGTS_aAR(n = n, alphap = alphap, alpham = alpham,
                            deltap = deltap, deltam = deltam,
                            lambdap = lambdap, lambdam = lambdam,
                            mu = mu, c = c),
              SR = rGTS_SR(n = n, alphap = alphap, alpham = alpham,
                            deltap = deltap, deltam = deltam,
                            lambdap = lambdap, lambdam = lambdam,
                            mu = mu, k = k))
  return(x)
}

rGTS_aAR <- function(n, alphap, alpham, deltap, deltam, lambdap, lambdam, mu, c)
{
  reVal <- rCTS_aAR(n/2, alphap,deltap,deltam,lambdap,lambdam,mu,0)
  reVal <- append(reVal,
                  rCTS_aAR(n/2, alpham,deltap,deltam,lambdap,lambdam,mu,c))
  return(reVal)
}

rGTS_SR <- function(n, alphap, alpham, deltap, deltam, lambdap, lambdam, mu, k) {
  reVal <- replicate(n = n/2,
            rCTS_SRp(alpha = alphap, delta = deltap, lambda = lambdap, k = k)
            -rCTS_SRp(alpha = alphap, delta = deltam, lambda = lambdam, k = k)
            +mu)
  reVal <- append(reVal,
                  replicate(n = n/2,
                            rCTS_SRp(alpha = alpham, delta = deltap,
                                     lambda = lambdap, k = k)
                            -rCTS_SRp(alpha = alpham, delta = deltam,
                                      lambda = lambdam, k = k)
                            +mu))
  return(reVal)
}


#### Kim-Rachev Tempered Stable Distribution ####

#' Characteristic function of the Kim-Rachev tempered stable distribution
#'
#' gapholder
#'
#' Gapholder. Currently different function with "functionOrigin" is no option,
#' as the second is not working.
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param kp,km Gapholder.
#' @param rp,rm Gapholder.
#' @param pp,pm Gapholder.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#'
#' @return The CF of the the Kim-Rachev tempered stable distribution.
#'
#' @references
#' gapholder
#'
#' @examples
#' gapholder
#'
#' @export
charKRTS <- function(t, alpha = NULL, kp = NULL, km = NULL, rp = NULL,
                     rm = NULL, pp = NULL, pm = NULL, mu = NULL, theta = NULL){
  if ((missing(alpha) | missing(kp) | missing(km) | missing(rp) |
       missing(rm) | missing(pp) | missing(pm) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, kp, km, rp, rm, pp, pm, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(kp) | missing(km) | missing(rp) |
      missing(rm) | missing(pp) | missing(pm) | missing(mu)) {
    alpha <- theta[1]
    kp <- theta[2]
    km <- theta[3]
    rp <- theta[4]
    rm <- theta[5]
    pp <- theta[6]
    pm <- theta[7]
    mu <- theta[8]
  }

  #In GMM kam es vor, dass pp == 0 ist für eine Simulation.
  #Deswegen die Anpassung
  if(pp == 0) pp <- 0.0001
  if(pm == 0) pm <- 0.0001
  if(pp == -1) pp <- -0.9999
  if(pm == -1) pm <- -0.9999

  stopifnot(0 < alpha, alpha < 2, alpha != 1,  0 < kp, 0 < km, 0 < rp,
            0 < rm, pp > -alpha, pp != -1, pp != 0, pm > -alpha, pm != -1,
            pm != 0)

  # Ansatz 1: Rachev et al. 2011 S.74
  subfunctionH <- function(x, alpha, r, p){
    (gamma(-alpha)) / p * (modifiedHyperGeoXc(p, -alpha, 1 + p, r*x) - 1)
  }

  return(exp(imagN * t * mu - imagN * t * gamma(1 - alpha) *
               ((kp*rp) / (pp + 1) - (km*rm) / (pm + 1))
             + kp * subfunctionH(imagN*t, alpha, rp, pp)
             + km * subfunctionH(-imagN*t, alpha, rm, pm)))

  # Ansatz 2: Kim et al 2009: A New Tempered Stable Distribution
  # Der Ansatz lässt sich nicht mit Zahlen aus dem Paper überprüfen, bzw. legt
  # nahe, dass etwas falsch ist
  # subfunctionH <- function(alpha, t, a, h, p){
  #   (a*gamma(-alpha))/p * (modifiedHyperGeoXc(p, -alpha, 1 + p, imagN*h*t) -1)
  # }
  #
  # return(exp(subfunctionH(alpha, t, kp, rp, pp)
  #            + subfunctionH(alpha, -t, km, rm, pm)
  #            + imagN * t * (mu + alpha * gamma(- alpha) *
  #                             (((kp*rp) / (pp + 1)) - ((km*rm) / (pm + 1))))
  #            )
  #        )
}


# Bsp. aus Rachev et al. 2011 S.74
# xKRTS <- seq(-3,3,0.01)
# alpha <- 1.25; pp <- 10; pm <- 10; rp <- 1/10; rm <- 1/2;
# kp <- (alpha + pp)*rp^(-alpha); km <- (alpha + pm)*rm^(-alpha); mu <- 0;
# yKRTS <- dKRTS(xKRTS, alpha,kp,km,rp,rm,pp,pm,mu)
dKRTS <- function(x, alpha = NULL, kp = NULL, km = NULL, rp = NULL,
                  rm = NULL, pp = NULL, pm = NULL, mu = NULL, theta = NULL,
                  dens_method = "FFT", a = -20, b = 20, nf = 256){
  if ((missing(alpha) | missing(kp) | missing(km) | missing(rp) |
       missing(rm) | missing(pp) | missing(pm) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, kp, km, rp, rm, pp, pm, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(kp) | missing(km) | missing(rp) |
      missing(rm) | missing(pp) | missing(pm) | missing(mu)) {
    alpha <- theta[1]
    kp <- theta[2]
    km <- theta[3]
    rp <- theta[4]
    rm <- theta[5]
    pp <- theta[6]
    pm <- theta[7]
    mu <- theta[8]
  }
  stopifnot(0 < alpha, alpha < 2, alpha != 1,  0 < kp, 0 < km, 0 < rp,
            0 < rm, pp > -alpha, pp != -1, pp != 0, pm > -alpha, pm != -1,
            pm != 0)

  if (dens_method == "FFT" || .Platform$OS.type != "windows") {
    d <- sapply(x, d_FFT, charFunc = charKRTS,
                theta = c(alpha, kp, km, rp, rm, pp, pm, mu),
                a = a, b = b, nf = nf)
  } else {
    d <- NULL
  }
  return(d)
}


pKRTS <- function(q, alpha = NULL, kp = NULL, km = NULL, rp = NULL,
                     rm = NULL, pp = NULL, pm = NULL, mu = NULL, theta = NULL,
                     dens_method = "FFT", a = -40, b = 40, nf = 2048, ...){
  if ((missing(alpha) | missing(kp) | missing(km) | missing(rp) |
       missing(rm) | missing(pp) | missing(pm) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, kp, km, rp, rm, pp, pm, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(kp) | missing(km) | missing(rp) |
      missing(rm) | missing(pp) | missing(pm) | missing(mu)) {
    alpha <- theta[1]
    kp <- theta[2]
    km <- theta[3]
    rp <- theta[4]
    rm <- theta[5]
    pp <- theta[6]
    pm <- theta[7]
    mu <- theta[8]
  }
  stopifnot(0 < alpha, alpha < 2, alpha != 1,  0 < kp, 0 < km, 0 < rp,
            0 < rm, pp > -alpha, pp != -1, pp != 0, pm > -alpha, pm != -1,
            pm != 0)

  p <- numeric(length(q))

  p <- sapply(q,
              function(z) {
                if(z<a) 0
                else if (z>b) 1
                else min(stats::integrate(dKRTS, lower = a, upper = z,
                                          alpha = alpha,
                                          kp = kp, km = km,
                                          rp = rp, rm = rm,
                                          pp = pp, pm = pm,
                                          mu = mu, ...)
                         $value, 1 - 1e-07)})

  return(p)
}

rKRTS <- function(n, alpha = NULL, kp = NULL, km = NULL, rp = NULL, rm = NULL,
                  pp = NULL, pm = NULL, mu = NULL, theta = NULL, methodR = "SR",
                  k = 10000) {
  if ((missing(alpha) | missing(kp) | missing(km) | missing(rp) |
       missing(rm) | missing(pp) | missing(pm) | missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, kp, km, rp, rm, pp, pm, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(kp) | missing(km) | missing(rp) |
      missing(rm) | missing(pp) | missing(pm) | missing(mu)) {
    alpha <- theta[1]
    kp <- theta[2]
    km <- theta[3]
    rp <- theta[4]
    rm <- theta[5]
    pp <- theta[6]
    pm <- theta[7]
    mu <- theta[8]
  }
  stopifnot(0 < alpha, alpha < 2, alpha != 1,  0 < kp, 0 < km, 0 < rp,
            0 < rm, pp > -alpha, pp != -1, pp != 0, pm > -alpha, pm != -1,
            pm != 0)

  #TODO:Insert other methods
  if(methodR == "TM" || methodR == "AR"){
    methodR <- "SR"
  }

  x <- switch(methodR,
              AR = 0,
              SR = rKRTS_SR(n, alpha, kp, km, rp, rm, pp, pm, mu, k),
              TM = 0)
  return(x)
}

rKRTS_SR <- function(n, alpha, kp, km, rp, rm, pp, pm, mu, k){
  replicate (n = n, rKRTS_SR_Ro(alpha, kp, km, rp, rm, pp, pm, k) + mu)
}

rKRTS_SR_Ro <- function(alpha, kp, km, rp, rm, pp, pm, k){
  parrivalslong <- cumsum(stats::rexp(k * 1.1))
  parrivals <- parrivalslong[parrivalslong <=  k]
  E1 <- stats::rexp(length(parrivals))
  U <- stats::runif(length(parrivals))

  sigma <- (kp*(rp^alpha))/(alpha+pp) + (km*(rm^alpha))/(alpha+pm)
  V <- rKRTS_SR_rVj(length(parrivals), sigma, alpha, kp, km, rp, rm, pp, pm, k)

  x0 <- sigma^(-1) * ((kp*(rp^alpha))/(alpha+pp) - (km*(rm^alpha))/(alpha+pm))
  x1 <- (kp*rp)/(pp+1) - (km*rm)/(pm+1)
  b <- alpha^(-1/alpha) * VGAM::zeta(1/alpha) * (sigma)^(1/alpha) *
    x0 - gamma(1-alpha) * x1
  cntr <- sum((alpha*(1:length(parrivals))/(sigma))^(-1/alpha)*x0)

  X <- cbind((alpha * parrivals/(sigma ))^(-1/alpha), E1 * U^(1/alpha)/abs(V))
  Xreturn <- sum((apply(X, 1, FUN = min)*V/abs(V)))-cntr+b
  return(Xreturn)
}

rKRTS_SR_dVj <- function(x, sigma, alpha, kp, km, rp, rm, pp, pm){
  returnVec <- NULL

  for(xi in x){
    Irp <- 0
    if (xi>(1/rp)) Irp <- 1

    Irm <- 0
    if (xi<(-1/rm)) Irm <- 1

    if(xi == 0){
      y <- 0
    }
    else {
      y <- (1/sigma * (kp*rp^(-pp)*Irp*(abs(xi)^(-alpha-pp-1)) +
                         km*rm^(-pm)*Irm*(abs(xi)^(-alpha-pm-1))))
    }
    returnVec <- append(returnVec,y)
  }

  return(returnVec)
}

rKRTS_SR_rVj <- function(n, sigma, alpha, kp, km, rp, rm, pp, pm, k){
  dX <- (20*2/k)
  x <- seq(-20,20,dX)
  y <- rKRTS_SR_dVj(x, sigma, alpha, kp, km, rp, rm, pp, pm)
  cumY <- cumsum(y)
  rV <- runif(n, min(cumY), max(cumY))

  returnVector <- NULL
  for(s in rV){
    pos <- which.min(abs(cumY - s))
    returnVector <- append(returnVector, x[pos])
  }

  return(returnVector)
}


#### Rapidly Decreasing Tempered Stable Distribution ####

#' Characteristic function of the rapidly decreasing tempered stable
#' distribution
#'
#' gapholder
#'
#' Gapholder.
#'
#' @param t A vector of real numbers where the CF is evaluated.
#' @param alpha Stability parameter. A real number between 0 and 1.
#' @param delta Scale parameter. A real number > 0.
#' @param lambdap,lambdam Tempering parameter. A real number > 0.
#' @param mu A location parameter, any real number.
#' @param theta Parameters stacked as a vector.
#'
#' @return The CF of the the rapidly decreasing tempered stable distribution.
#'
#' @references
#' gapholder
#'
#' @examples
#' gapholder
#'
#' @export
charRDTS <- function(t, alpha = NULL, delta = NULL, lambdap = NULL,
                     lambdam = NULL, mu = NULL, theta = NULL) {
  if ((missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
       missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
      missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, alpha != 1, 0 < delta, 0 < lambdap,
            0 < lambdam)

  subfunctionG <- function(x, alpha, lambda){
    (2^(-(alpha/2) - 1) * lambda^(alpha) * gamma(-(alpha/2)) *
       (hypergeo::genhypergeo(-alpha/2, 1/2, (x^2)/(2*lambda^2)) - 1) +
       2^(-alpha/2 - 1/2) * lambda^(alpha-1) * x * gamma((1-alpha)/2) *
       (hypergeo::genhypergeo((1-alpha)/2, 3/2, (x^2)/(2*lambda^2)) - 1))
  }

  returnVec <- NULL
  nextG <- NULL
  counterL <- 1

  # Each value is compared with the two values before the examined value and
  # the one after, because there are errors in the function
  # hypergeo::genhypergeo().
  for(x in t){

    if(is.null(nextG)){
      g <- exp(imagN * x * mu +
                 delta * (subfunctionG(imagN * x, alpha, lambdap)
                          + subfunctionG(-imagN * x, alpha, lambdam)))
    }
    else g <- nextG

    if(counterL < length(t)){
      nextG <- exp(imagN * t[counterL+1] * mu +
                     delta * (subfunctionG(imagN * t[counterL+1], alpha,
                                           lambdap)
                              + subfunctionG(-imagN * t[counterL+1], alpha,
                                             lambdam)))
    }

    if (is.nan(nextG) || is.null(nextG) || is.infinite(Re(nextG)) ||
        is.infinite(Im(nextG))){
      nextG <- 0 + 0i
    }

    if (is.nan(g) || is.infinite(Re(g)) || is.infinite(Im(g))){
      returnVec <- append(returnVec,0 + 0i)
    }
    else if(!is.null(returnVec)){
      avLast2Re <- mean(c(Re(g),Re(returnVec[length(returnVec)]), Re(nextG)))
      #avLast2Im <- mean(c(Im(g),Im(returnVec[length(returnVec)]), Im(nextG)))

      # Only real numbers must be checked as all errors can be catched with it.
      if (Re(g) != 0 && (avLast2Re/Re(g) < 0.75 || avLast2Re/Re(g) > 1.5)){
        returnVec <- append(returnVec,0 + 0i)
      }
      else returnVec <- append(returnVec, g)
    }
    else returnVec <- append(returnVec, g)

    counterL <- counterL +1

  }

  return(returnVec)
}


dRDTS <- function(x, alpha = NULL, delta = NULL, lambdap = NULL,
                  lambdam = NULL, mu = NULL, theta = NULL, dens_method = "FFT",
                  a = -20, b = 20, nf = 256) {
  if ((missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
       missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
      missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 < lambdam)

  if (dens_method == "FFT" || .Platform$OS.type != "windows") {
    d <- sapply(x, d_FFT, charFunc = charRDTS,
                theta = c(alpha, delta, lambdap, lambdam, mu),
                a = a, b = b, nf = nf)
  } else {
    d <- NULL
  }
  return(d)

}

pRDTS <- function(p, alpha = NULL, delta = NULL, lambdap = NULL,
                  lambdam = NULL, mu = NULL, theta = NULL, dens_method = "FFT",
                  a = -130, b = 130, nf = 2048, ...) {
  if ((missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
       missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
      missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 < lambdam)

  p <- numeric(length(q))

  p <- sapply(q,
              function(z) {
                if(z<a) 0
                else if (z>b) 1
                else min(stats::integrate(dRDTS, lower = a, upper = z,
                                          alpha = alpha,
                                          delta = delta,
                                          lambdap = lambdap,
                                          lambdam = lambdam,
                                          mu = mu, ...)
                         $value, 1 - 1e-07)})

  return(p)

}

rRDTS <- function(n, alpha = NULL, delta = NULL, lambdap = NULL, lambdam = NULL,
                  mu = NULL, theta = NULL, methodR = "TM", k = 10000){
  if ((missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
       missing(mu)) & is.null(theta))
    stop("No or not enough parameters supplied")
  theta0 <- c(alpha, delta, lambdap, lambdam, mu)
  if (!is.null(theta) & !is.null(theta0)) {
    if (!all(theta0 == theta))
      stop("Parameters do not match")
  }
  if (missing(alpha) | missing(delta) | missing(lambdap) | missing(lambdam) |
      missing(mu)) {
    alpha <- theta[1]
    delta <- theta[2]
    lambdap <- theta[3]
    lambdam <- theta[4]
    mu <- theta[5]
  }
  stopifnot(0 < alpha, alpha < 2, 0 < delta, 0 < lambdap, 0 < lambdam)

  #TODO:Insert other methods
  if(methodR == "TM" || methodR == "AR"){
    methodR <- "SR"
  }

  x <- switch(methodR,
              AR = 0,
              SR = rRDTS_SR(n = n, alpha = alpha, delta = delta, lambdap =
                              lambdap, lambdam = lambdam, mu = mu, k = k),
              TM = 0)
  return(x)

}

# Function from here: Kim et al 2010 Tempered stable and tempered infinitely
# divisble GARCH models. BUT, its modified by factor: 1/(1+((1-alpha)/2)), as
# function gave wrong values before
rRDTS_SR <- function(n, alpha, delta, lambdap, lambdam, mu, k) {
  replicate(n = n, (rTSS_SR2(alpha = alpha, delta = delta,
                            lambda = lambdap, k = k) -
              rTSS_SR2(alpha = alpha, delta = delta,
                       lambda = lambdam, k = k))
            - (delta*gamma((1-alpha)/2)/(2^((alpha+1)/2)))*
              (lambdap^(alpha-1)-lambdam^(alpha-1))
            + mu)
}


#### Supportive functions ####

#Die FFT wird bei vielen temperierten Verteilungen angewendet, um an die
# Dichtefunktion zu gelangen. Der Weg dahin kann für alle Verteilungen
# verallgemeinert werden.
d_FFT <- function(x, charFunc, theta, a, b, nf){
  dx <- ((b - a)/nf)
  dt <- (2 * pi/(nf * dx))

  sq <- seq(from = 0, to = nf - 1, 1)
  t <- (-nf/2 * dt + sq * dt)
  xgrid <- (a + dx * sq)

  cft <- charFunc(t = t, theta = theta)

  tX <- (exp(-imagN * sq * dt * a) * cft)

  # fast Fourier transform
  y <- stats::fft(tX, inverse = FALSE)

  densityW <- Re((dt/(2 * pi)) * exp(-imagN * (nf/2 * dt) * xgrid) * y)

  return (as.numeric(stats::approx(xgrid, densityW, xout=x,
                                   yleft = 1e-18, yright = 1e-18)[2]))
}


# Die Hypergeometrische Funktion funktioniert nur für einen sehr kleinen Werte-
# bereich. Hiermit wird sie angepasst
#Modified hypergeometric function for real values.
modifiedHyperGeoXr <- function(a, b, c, x){
  returnVec <- NULL

  for (z in x){
    if(z == 0 | z == 0+0i){
      returnVec <- append(returnVec, 1)
    }

    else if(abs(z/(z-1)) < 1){
      returnVec <- append(returnVec,
                          (1 - z) ^ (-b) * gsl::hyperg_2F1(b, c - a, c, z /
                                                             (z - 1)))
    }

    else if (abs(1/z) < 1){
      returnVec <-
        append(returnVec,
               (-z)^(-a) * (gamma(c)*gamma(b-a))/(gamma(c-a)*gamma(b)) *
                 gsl::hyperg_2F1(a, a-c+1, a-b+1, 1/z) +
                 (-z)^(-b) * (gamma(c)*gamma(a-b))/(gamma(c-b)*gamma(a)) *
                 gsl::hyperg_2F1(b, b-c+1, b-a+1, 1/z)
               )
    }

    else{
      returnVec <- append(returnVec,
                          gsl::hyperg_2F1(a, b, c, z)
      )
    }
  }
  returnVec
}


#Modified hypergeometric function for complex values.
modifiedHyperGeoXc <- function(a, b, c, x){
  returnVec <- NULL

  for (z in x){
    if(z == 0 | z == 0+0i){
      returnVec <- append(returnVec, 1)
    }

    else if(abs(Re(z)/(Re(z)-1)) < 1){
      returnVec <- append(returnVec,
                          (1 - z) ^ (-b) * hypergeo::hypergeo(b, c - a, c, z /
                                                                (z - 1)))
    }

    else if (abs(1/Re(z)) < 1){
      returnVec <-
        append(returnVec,
               (-z)^(-a) * (gamma(c)*gamma(b-a))/(gamma(c-a)*gamma(b)) *
                 hypergeo::hypergeo(a, a-c+1, a-b+1, 1/z) +
                 (-z)^(-b) * (gamma(c)*gamma(a-b))/(gamma(c-b)*gamma(a)) *
                 hypergeo::hypergeo(b, b-c+1, b-a+1, 1/z)
               )
    }

    else{
      returnVec <- append(returnVec,
                          hypergeo::hypergeo(a, b, c, z)
      )
    }
  }
  returnVec
}

modifiedGenHyperGeo <- function(a, c, x){
  returnVec <- NULL

  for (z in x){
    g <- hypergeo::genhypergeo(a,c,z)
    if (is.nan(g)){
      returnVec <- append(returnVec,0 + 0i)
    }
    else if (Re(g) < 1e-8 && Im(g) < 1e-8 || Re(g) > 1e+8 && Im(g) > 1e+8){
      returnVec <- append(returnVec,0 + 0i)
    }
    else returnVec <- append(returnVec,g)
  }
  returnVec
}


rVariantesTS <- function(n, densityFunction, x, theta, ...){
  y <- densityFunction(x = x, theta = theta, ...)
  cumY <- cumsum(y)
  rV <- runif(n, min(cumY), max(cumY))

  returnVector <- NULL
  for(s in rV){
    pos <- which.min(abs(cumY - s))
    returnVector <- append(returnVector, x[pos])
  }

  return(returnVector)
}

















