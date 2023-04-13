#' @importFrom moments all.moments
CumFinder_GTS <- function(x, jmax = 7) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}


MoCObjective_GTS <- function(x, parms) {
    c(F1 = x[7] - parms[1],
      F2 = gamma(2 - x[1]) * (x[3]/x[5]^(2 - x[1]) + x[4]/x[6]^(2 - x[2])) -
        parms[2],
      F3 = gamma(3 - x[1]) * (x[3]/x[5]^(3 - x[1]) - x[4]/x[6]^(3 - x[2])) -
        parms[3],
      F4 = gamma(4 - x[1]) * (x[3]/x[5]^(4 - x[1]) + x[4]/x[6]^(4 - x[2])) -
        parms[4],
      F5 = gamma(5 - x[1]) * (x[3]/x[5]^(5 - x[1]) - x[4]/x[6]^(5 - x[2])) -
        parms[5],
      F6 = gamma(6 - x[1]) * (x[3]/x[5]^(6 - x[1]) + x[4]/x[6]^(6 - x[2])) -
        parms[6],
      F7 = gamma(7 - x[1]) * (x[3]/x[5]^(7 - x[1]) - x[4]/x[6]^(7 - x[2])) -
        parms[7])
}

#' @importFrom rootSolve multiroot
MoC_GTS <- function(x, theta0 = c(1.5, 1.5, 1, 1, 1, 1, 0), eps = 1e-06) {
    cumulants <- CumFinder_GTS(x)
    utils::capture.output(
      parroot <-
        rootSolve::multiroot(MoCObjective_GTS, theta0, parms = cumulants)
    )
    theta <- parroot$root
    if (theta[1] < 0 && theta[2] < 0) {
        theta <- theta0
    }
    if (theta[1] > 2 && theta[2] > 2) {
        theta <- theta0
    }
    return(theta)
}
