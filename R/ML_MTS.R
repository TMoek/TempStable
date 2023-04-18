MLParametersEstim_MTS <-
  function(x,
           theta0 = NULL,
           eps = 1e-06,
           PrintTime = FALSE,
           ...) {
    t_init <- StableEstim::getTime_()
    if (is.null(theta0))
      theta0 <- MoC_MTS(x, c(0.6, 1, 1, 1, 0), eps = eps)
    method <- .methodDesML_MTS()
    dots <- list(...)
    if (is.null(dots$control)) {
      control <- list(factr = 1e+05, pgtol = 1e-05)
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_MTS,
          gr = NULL,
          x = x,
          control = control,
          ...,
          method = "L-BFGS-B",
          lower = c(-Inf, eps, eps, eps, -Inf),
          upper =
            c(1 - eps, Inf, Inf, Inf, Inf)
        )
    } else {
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_MTS,
          gr = NULL,
          x = x,
          ...,
          method = "L-BFGS-B",
          lower = c(-Inf, eps, eps, eps, -Inf),
          upper =
            c(1 - eps, Inf, Inf, Inf, Inf)
        )
    }
    if (PrintTime)
      StableEstim::PrintDuration(
        StableEstim::ComputeDuration(t_init, t_final <-
                                       StableEstim::getTime_()),
        "Normal_MLParametersEstim_MTS")
    list(
      Estim = Estim,
      duration = StableEstim::ComputeDuration(t_init, StableEstim::getTime_(),
                                              TRUE),
      method = method
    )
  }

SumLogDensity_MTS <- function(theta, x, sign = -1, ...) {
  densis <-
    sapply(
      X = x,
      FUN = dMTS,
      alpha = theta[1],
      delta = theta[2],
      lambdap = theta[3],
      lambdam = theta[4],
      mu = theta[5]
    )
  sign * sum(log(densis))
}


.methodDesML_MTS <- function(...) {
  l <- list(...)
  paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}
