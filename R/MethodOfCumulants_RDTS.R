#' @importFrom moments all.moments
CumFinder_RDTS <- function(x, jmax = 5) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}


MoCObjective_RDTS <- function(x, parms) {
    c(F1 = x[5] - parms[1],
      F2 = 2^((2-x[1]-2)/2) * x[2] * gamma((2-x[1])/2) *
        (x[3]^(x[1]-2)+x[4]^(x[1]-2)),
      F3 = 2^((3-x[1]-2)/2) * x[2] * gamma((3-x[1])/2) *
        (x[3]^(x[1]-3)-x[4]^(x[1]-3)),
      F4 = 2^((4-x[1]-2)/2) * x[2] * gamma((4-x[1])/2) *
        (x[3]^(x[1]-4)+x[4]^(x[1]-4)),
      F5 = 2^((5-x[1]-2)/2) * x[2] * gamma((5-x[1])/2) *
        (x[3]^(x[1]-5)-x[4]^(x[1]-5)),
        )
}

#' @importFrom rootSolve multiroot
MoC_RDTS <- function(x, theta0 = c(0.5, 1, 1, 1, 0), eps = 1e-06) {
    cumulants <- CumFinder_RDTS(x)
    utils::capture.output(parroot <-
        rootSolve::multiroot(MoCObjective_RDTS, theta0, parms = cumulants))
    theta <- parroot$root
    if (theta[1] < 0) {
        theta <- theta0
    }
    if (theta[1] > 2) {
        theta <- theta0
    }
    return(theta)
}
