MLParametersEstim_CTS <- function(x, theta0 = NULL, eps = 1e-06, PrintTime = FALSE, ...) {
    t_init <- getTime_()
    if (is.null(theta0))
        theta0 <- MoC_CTS(x, c(1.5, 1, 1, 1, 1, 0), eps = eps)
    method <- .methodDesML_CTS()
    dots <- list(...)
    if (is.null(dots$control)) {
        control <- list(factr = 1e+05, pgtol = 1e-05)
        Estim <- optim(par = theta0, fn = SumLogDensity_CTS, gr = NULL, x = x, control = control, ..., method = "L-BFGS-B",
            lower = c(eps, eps, eps, eps, eps, -Inf), upper = c(2 - eps, Inf, Inf, Inf, Inf, Inf))
    } else {
        Estim <- optim(par = theta0, fn = SumLogDensity_CTS, gr = NULL, x = x, ..., method = "L-BFGS-B", lower = c(eps,
            eps, eps, eps, eps, -Inf), upper = c(2 - eps, Inf, Inf, Inf, Inf, Inf))
    }
    if (PrintTime)
        PrintDuration(ComputeDuration(t_init, t_final <- getTime_()), "Classic_MLParametersEstim_CTS")
    list(Estim = Estim, duration = ComputeDuration(t_init, getTime_(), TRUE), method = method)
}

SumLogDensity_CTS <- function(theta, x, sign = -1, ...) {
    densis <- sapply(X = x, FUN = dCTS, alpha = theta[1], deltap = theta[2], deltam = theta[3], lambdap = theta[4], lambdam = theta[5],
        mu = theta[6])
    sign * sum(log(densis))
}

.methodDesML_CTS <- function(...) {
    l <- list(...)
    paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}

