MLParametersEstim_KRTS <-
  function(x,
           theta0 = NULL,
           eps = 1e-06,
           PrintTime = FALSE,
           ...) {
    t_init <- StableEstim::getTime_()
    if (is.null(theta0))
      theta0 <- MoC_KRTS(x, c(1.5, 1, 1, 1, 1, 1, 1, 0), eps = eps)
    method <- .methodDesML_KRTS()
    dots <- list(...)
    if (is.null(dots$control)) {
      control <- list(factr = 1e+05, pgtol = 1e-05)
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_KRTS,
          gr = NULL,
          x = x,
          control = control,
          ...,
          method = "L-BFGS-B",
          lower = c(eps, eps, eps, eps, eps,
                    -theta0[1] -eps,
                    -theta0[1] -eps, -Inf),
          upper = c(2 - eps, Inf, Inf, Inf, Inf,
                    Inf, Inf, Inf)
        )
    } else {
      Estim <-
        stats::optim(
          par = theta0,
          fn = SumLogDensity_KRTS,
          gr = NULL,
          x = x,
          ...,
          method = "L-BFGS-B",
          lower = c(eps, eps, eps, eps, eps,
                    -theta0[1] -eps,
                    -theta0[1] -eps, -Inf),
          upper = c(2 - eps, Inf, Inf, Inf, Inf,
                    Inf, Inf, Inf)
        )
    }
    if (PrintTime)
      StableEstim::PrintDuration(
        StableEstim::ComputeDuration(t_init, t_final <-
                                       StableEstim::getTime_()),
        "Classic_MLParametersEstim_KRTS")
    list(
      Estim = Estim,
      duration = StableEstim::ComputeDuration(t_init, StableEstim::getTime_(), TRUE),
      method = method
    )
  }

SumLogDensity_KRTS <- function(theta, x, sign = -1, ...) {
  densis <-
    sapply(
      X = x,
      FUN = dKRTS,
      alpha <- theta[1],
      kp <- theta[2],
      km <- theta[3],
      rp <- theta[4],
      rm <- theta[5],
      pp <- theta[6],
      pm <- theta[7],
      mu <- theta[8]
    )
  sign * sum(log(densis))
}

.methodDesML_KRTS <- function(...) {
  l <- list(...)
  paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}
