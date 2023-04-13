#' @importFrom moments all.moments
CumFinder_MTS <- function(x, jmax = 5) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}

# Function origin is kim09. See charMTS().
MoCObjective_MTS <- function(x, parms) {

  c(F1 = x[5] + 2^(-x[1]-1/2)* x[2]*gamma((1/2)-x[1]) *
      (x[3]^(2*x[1]-1)-x[4]^(2*x[1]-1)),
    F2 = 2^(-x[1]-3/2) * sqrt(pi) * (factorial(2)/factorial(2/2)) * x[2] *
      gamma((2/2)-x[1]) * (x[3]^(2*x[1]-2)+x[4]^(2*x[1]-2)),
    F3 = 2^(3-x[1]-3/2) * factorial((3-1)/2) * x[2] * gamma(3/2-x[1]) *
      (x[3]^(2*x[1]-3)+x[4]^(2*x[1]-3)),
    F4 = 2^(-x[1]-3/2) * sqrt(pi) * (factorial(4)/factorial(4/2)) * x[2] *
      gamma((4/2)-x[1]) * (x[3]^(2*x[1]-4)+x[4]^(2*x[1]-4)),
    F5 = 2^(5-x[1]-3/2) * factorial((5-1)/2) * x[2] * gamma(5/2-x[1]) *
      (x[3]^(2*x[1]-5)+x[4]^(2*x[1]-5)),
    )
}

# Function origin is kim09. See charMTS().
#' @importFrom rootSolve multiroot
MoC_MTS <- function(x, theta0 = c(0.6, 1, 1, 1, 0), eps = 1e-06) {
    cumulants <- CumFinder_MTS(x)
    utils::capture.output(parroot <-
        rootSolve::multiroot(MoCObjective_MTS, theta0, parms = cumulants))
    theta <- parroot$root

    if (theta[1] > 1) {
        theta <- theta0
    }
    return(theta)
}
