#' @importFrom moments all.moments
CumFinder_NTS <- function(x, jmax = 5) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}


MoCObjective_NTS <- function(x, parms) {
    c(F1 = x[5] + x[2] * x[3]/x[4]^(1 - x[1]) * gamma(1 - x[1]) - parms[1], F2 = gamma(2 - x[1]) * x[2]^2 * x[3]/x[4]^(2 -
        x[1]) - x[3]/x[4]^(1 - x[1]) * gamma(1 - x[1]) - parms[2], F3 = gamma(3 - x[1]) * x[2]^3 * x[3]/x[4]^(3 - x[1]) -
        3 * x[2] * x[3]/x[4]^(2 - x[1]) * gamma(2 - x[1]) - parms[3], F4 = gamma(4 - x[1]) * x[2]^4 * x[3]/x[4]^(4 - x[1]) -
        6 * x[2]^2 * x[3]/x[4]^(3 - x[1]) * gamma(3 - x[1]) + 3 * x[3]/x[4]^(2 - x[1]) * gamma(2 - x[1]) - parms[4], F5 = gamma(5 -
        x[1]) * x[2]^5 * x[3]/x[4]^(5 - x[1]) - 10 * x[2]^3 * x[3]/x[4]^(4 - x[1]) * gamma(4 - x[1]) + 15 * x[2] * x[3]/x[4]^(3 -
        x[1]) * gamma(3 - x[1]) - parms[5])
}

MoC_NTS <- function(x, theta0 = c(0.5, 0, 1, 1, 0), eps = 1e-06) {
    cumulants <- CumFinder_NTS(x)
    parroot <- multiroot(MoCObjective_NTS, theta0, parms = cumulants)
    theta <- parroot$root
    if (theta[1] < 0) {
        theta <- theta0
    }
    if (theta[1] > 1) {
        theta <- theta0
    }
    return(theta)
}
