MLParametersEstim_NTS <- function(x, theta0 = NULL, eps = 1e-06, PrintTime = FALSE, ...) {
    t_init <- getTime_()
    if (is.null(theta0))
        theta0 <- MoC_NTS(x, c(0.5, 0, 1, 1, 0), eps = eps)
    method <- .methodDesML_NTS()
    dots <- list(...)
    if (is.null(dots$control)) {
        control <- list(factr = 1e+05, pgtol = 1e-05)
        Estim <- optim(par = theta0, fn = SumLogDensity_NTS, gr = NULL, x = x, control = control, ..., method = "L-BFGS-B",
            lower = c(eps, -Inf, eps, eps, -Inf), upper = c(1 - eps, Inf, Inf, Inf, Inf))
    } else {
        Estim <- optim(par = theta0, fn = SumLogDensity_NTS, gr = NULL, x = x, ..., method = "L-BFGS-B", lower = c(eps,
            -Inf, eps, eps, -Inf), upper = c(1 - eps, Inf, Inf, Inf, Inf))
    }
    if (PrintTime)
        PrintDuration(ComputeDuration(t_init, t_final <- getTime_()), "Normal_MLParametersEstim_NTS")
    list(Estim = Estim, duration = ComputeDuration(t_init, getTime_(), TRUE), method = method)
}

SumLogDensity_NTS <- function(theta, x, sign = -1, ...) {
    densis <- sapply(X = x, FUN = dNTS, alpha = theta[1], beta = theta[2], delta = theta[3], lambda = theta[4], mu = theta[5])
    sign * sum(log(densis))
}


.methodDesML_NTS <- function(...) {
    l <- list(...)
    paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}

