MLParametersEstim_STS <-
  function(x,
           theta0 = NULL,
           eps = 1e-06,
           PrintTime = FALSE,
           ...) {
    t_init <- StableEstim::getTime_()
    if (is.null(theta0))
      theta0 <- MoC_STS(x, c(0.5, 1, 1), eps = eps)
    method <- .methodDesML_STS()
    dots <- list(...)
    if (is.null(dots$control)) {
      control <- list(factr = 1e+05, pgtol = 1e-05)
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_STS,
          gr = NULL,
          x = x,
          control = control,
          ...,
          method = "L-BFGS-B",
          lower = c(eps, eps, eps),
          upper = c(1 - eps, Inf, Inf)
        )
    } else {
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_STS,
          gr = NULL,
          x = x,
          ...,
          method = "L-BFGS-B",
          lower = c(eps,
                    eps, eps),
          upper = c(1 - eps, Inf, Inf)
        )
    }
    if (PrintTime)
      StableEstim::PrintDuration(
        StableEstim::ComputeDuration(t_init, t_final <-
                                       StableEstim::getTime_()),
        "Subordinator_MLParametersEstim_STS"
      )
    list(
      Estim = Estim,
      duration = StableEstim::ComputeDuration(t_init, StableEstim::getTime_(),
                                              TRUE),
      method = method
    )
  }

SumLogDensity_STS <- function(theta, x, sign = -1, ...) {
  sign * sum(log(dSTS(x, theta[1], theta[2], theta[3])))
}

.methodDesML_STS <- function(...) {
  l <- list(...)
  paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}
