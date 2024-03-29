#' @importFrom moments all.moments
CumFinder_CGMY <- function(x, jmax = 5) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}


MoCObjective_CGMY <- function(x, parms) {
    c(F1 = gamma(2 - x[4]) * (x[1]/x[3]^(2 - x[4]) + x[1]/x[4]^(2 - x[4])) -
        parms[2],
      F2 = gamma(3 - x[4]) * (x[1]/x[3]^(3 - x[4]) - x[1]/x[4]^(3 - x[4])) -
        parms[3],
      F3 = gamma(4 - x[4]) * (x[1]/x[3]^(4 - x[4]) + x[1]/x[4]^(4 - x[4])) -
        parms[4],
      F4 = gamma(5 - x[4]) * (x[1]/x[3]^(5 - x[4]) - x[1]/x[4]^(5 - x[4])) -
        parms[5])
}

#' @importFrom rootSolve multiroot
MoC_CGMY <- function(x, theta0 = c(1, 1, 1, 1.5), eps = 1e-06) {
    cumulants <- CumFinder_CGMY(x)
    parroot <- rootSolve::multiroot(MoCObjective_CGMY, theta0, parms = cumulants)
    theta <- parroot$root
    if (theta[4] < 0) {
        theta <- theta0
    }
    if (theta[4] > 2) {
        theta <- theta0
    }
    return(theta)
}
