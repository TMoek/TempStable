MLParametersEstim_GTS <-
  function(x,
           theta0 = NULL,
           eps = 1e-06,
           PrintTime = FALSE,
           ...) {
    t_init <- StableEstim::getTime_()
    if (is.null(theta0))
      theta0 <- MoC_GTS(x, c(1.5, 1.5, 1, 1, 1, 1, 0), eps = eps)
    method <- .methodDesML_GTS()
    dots <- list(...)
    if (is.null(dots$control)) {
      control <- list(factr = 1e+05, pgtol = 1e-05)
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_GTS,
          gr = NULL,
          x = x,
          control = control,
          ...,
          method = "L-BFGS-B",
          lower = c(eps, eps, eps, eps, eps, eps,
                    -Inf),
          upper = c(2 - eps, 2 - eps, Inf, Inf, Inf,
                    Inf, Inf)
        )
    } else {
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_GTS,
          gr = NULL,
          x = x,
          ...,
          method = "L-BFGS-B",
          lower = c(eps, eps, eps, eps, eps, eps,
                    -Inf),
          upper = c(2 - eps, 2 - eps, Inf, Inf, Inf,
                    Inf, Inf)
        )
    }
    if (PrintTime)
      StableEstim::PrintDuration(
        StableEstim::ComputeDuration(t_init, t_final <-
                                       StableEstim::getTime_()),
        "Classic_MLParametersEstim_GTS")
    list(
      Estim = Estim,
      duration = StableEstim::ComputeDuration(t_init, StableEstim::getTime_(), TRUE),
      method = method
    )
  }

SumLogDensity_GTS <- function(theta, x, sign = -1, ...) {
  densis <-
    sapply(
      X = x,
      FUN = dGTS,
      alphap = theta[1],
      alpham = theta[2],
      deltap = theta[3],
      deltam = theta[4],
      lambdap = theta[5],
      lambdam = theta[6],
      mu = theta[7]
    )
  sign * sum(log(densis))
}

.methodDesML_GTS <- function(...) {
  l <- list(...)
  paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}
