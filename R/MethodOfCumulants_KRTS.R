#' @importFrom moments all.moments
CumFinder_KRTS <- function(x, jmax = 8) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}


MoCObjective_KRTS <- function(x, parms) {
    c(F1 = x[8] - parms[1],
      F2 = gamma(2-x[1]) * ((x[2]*(x[4]^2))/(x[6]+2) + (x[3]*(x[5]^2))/(x[7]+2))
        - parms[2],
      F3 = gamma(3-x[1]) * ((x[2]*(x[4]^3))/(x[6]+3) - (x[3]*(x[5]^3))/(x[7]+3))
        - parms[3],
      F4 = gamma(4-x[1]) * ((x[2]*(x[4]^4))/(x[6]+4) + (x[3]*(x[5]^4))/(x[7]+4))
        - parms[4],
      F5 = gamma(5-x[1]) * ((x[2]*(x[4]^5))/(x[6]+5) - (x[3]*(x[5]^5))/(x[7]+5))
        - parms[5],
      F6 = gamma(6-x[1]) * ((x[2]*(x[4]^6))/(x[6]+6) + (x[3]*(x[5]^6))/(x[7]+6))
        - parms[6],
      F7 = gamma(7-x[1]) * ((x[2]*(x[4]^7))/(x[6]+7) - (x[3]*(x[5]^7))/(x[7]+7))
        - parms[7],
      F8 = gamma(8-x[1]) * ((x[2]*(x[4]^8))/(x[6]+8) + (x[3]*(x[5]^8))/(x[7]+8))
        - parms[8])
}

#' @importFrom rootSolve multiroot
MoC_KRTS <- function(x, theta0 = c(1.5, 1, 1, 1, 1, 1, 1, 0), eps = 1e-06) {
    cumulants <- CumFinder_KRTS(x)
    utils::capture.output(
      parroot <-
        rootSolve::multiroot(MoCObjective_KRTS, theta0, parms = cumulants)
    )
    theta <- parroot$root
    if (theta[1] < 0) {
        theta <- theta0
    }
    if (theta[1] > 2) {
        theta <- theta0
    }
    return(theta)
}
