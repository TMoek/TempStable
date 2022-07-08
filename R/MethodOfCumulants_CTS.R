#' @importFrom moments all.moments
CumFinder_CTS <- function(x, jmax = 6) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}


MoCObjective_CTS <- function(x, parms) {
    c(F1 = x[6] - parms[1],
      F2 = gamma(2 - x[1]) * (x[2]/x[4]^(2 - x[1]) + x[3]/x[5]^(2 - x[1]))
        - parms[2],
      F3 = gamma(3 - x[1]) * (x[2]/x[4]^(3 - x[1]) - x[3]/x[5]^(3 - x[1]))
        - parms[3],
      F4 = gamma(4 - x[1]) * (x[2]/x[4]^(4 - x[1]) + x[3]/x[5]^(4 - x[1]))
        - parms[4],
      F5 = gamma(5 - x[1]) * (x[2]/x[4]^(5 - x[1]) - x[3]/x[5]^(5 - x[1]))
        - parms[5],
      F6 = gamma(6 - x[1]) * (x[2]/x[4]^(6 - x[1]) + x[3]/x[5]^(6 - x[1]))
        - parms[6])
}

#' @importFrom rootSolve multiroot
MoC_CTS <- function(x, theta0 = c(1.5, 1, 1, 1, 1, 0), eps = 1e-06) {
    cumulants <- CumFinder_CTS(x)
    parroot <- rootSolve::multiroot(MoCObjective_CTS, theta0, parms = cumulants)
    theta <- parroot$root
    if (theta[1] < 0) {
        theta <- theta0
    }
    if (theta[1] > 2) {
        theta <- theta0
    }
    return(theta)
}
