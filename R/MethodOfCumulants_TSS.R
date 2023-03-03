#' @importFrom moments all.moments
CumFinder_TSS <- function(x, jmax = 3) {
    raw.moments <- moments::all.moments(x, order.max = jmax)
    cumfrommom(raw.moments)[-1]
}

#' @importFrom moments raw2central
cumfrommom <- function(moms) {
    mu.raw <- moms
    mu.central <- moments::raw2central(mu.raw)
    order.max <- length(mu.raw) - 1
    kappa <- rep(0, order.max + 1)
    kappa[2] <- mu.raw[2]
    kappa[3] <- mu.central[3]
    n <- 3
    while (n <= order.max) {
        np1 <- n + 1
        nm1 <- n - 1
        total <- mu.raw[np1]
        for (k in 1:nm1) {
            km1 <- k - 1
            total <- total - choose(nm1, km1) * kappa[k + 1] * mu.raw[n - k + 1]
        }
        kappa[np1] <- total
        n <- n + 1
    }
    return(kappa)
}

MoCObjective_TSS <- function(x, parms) {
    c(F1 = gamma(1 - x[1]) * x[2]/x[3]^(1 - x[1]) - parms[1],
      F2 = gamma(2 - x[1]) * x[2]/x[3]^(2 - x[1]) - parms[2],
      F3 = gamma(3 - x[1]) * x[2]/x[3]^(3 - x[1]) - parms[3])
}

#' @importFrom rootSolve multiroot
MoC_TSS <- function(x, theta0 = c(0.5, 1, 1), eps = 1e-06) {
    cumulants <- CumFinder_TSS(x)
    utils::capture.output(
      parroot <-
        rootSolve::multiroot(MoCObjective_TSS, theta0, parms = cumulants)
    )
    theta <- parroot$root
    if (theta[1] < 0) {
        theta <- theta0
    }
    if (theta[1] > 1) {
        theta <- theta0
    }
    return(theta)
}
