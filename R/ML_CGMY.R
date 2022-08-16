MLParametersEstim_CGMY <- function (x,
                                    theta0 = NULL,
                                    eps = 1e-6,
                                    PrintTime = FALSE,
                                    ...) {
  t_init <- StableEstim::getTime_()
  if (is.null(theta0))
    theta0 <- MoC_CGMY(x, c(1, 1, 1, 1.5), eps = eps)
  method <- .methodDesML_CGMY()
  dots <- list(...)
  if (is.null(dots$control)) {
    control <- list(factr = 1e+05, pgtol = 1e-05)
    Estim <- stats::optim(
      par = theta0,
      fn = SumLogDensity_CGMY,
      gr = NULL,
      x = x,
      control = control,
      ...,
      method = "L-BFGS-B",
      lower = c(eps, eps, eps, eps),
      upper = c(Inf, Inf, Inf, 2 - eps)
    )
  }
  else {
    Estim <- stats::optim(
      par = theta0,
      fn = SumLogDensity_CGMY,
      gr = NULL,
      x = x,
      ...,
      method = "L-BFGS-B",
      lower = c(eps, eps, eps, eps),
      upper = c(Inf, Inf, Inf, 2 - eps)
    )
  }
  if (PrintTime)
    StableEstim::PrintDuration(
      StableEstim::ComputeDuration(t_init, t_final <- StableEstim::getTime_()),
      "CGMY_MLParametersEstim")
  list(
    Estim = Estim,
    duration = StableEstim::ComputeDuration(t_init, StableEstim::getTime_(),
                               TRUE),
    method = method
  )
}

SumLogDensity_CGMY <- function (theta, x, sign = -1, ...)
{
  densis <-
    sapply(
      X = x,
      FUN = dCGMY,
      C = theta[1],
      G = theta[2],
      M = theta[3],
      Y = theta[4]
    )
  sign * sum(log(densis))
}

.methodDesML_CGMY <- function (...) {
  l <- list(...)
  paste("ML", paste("OptimAlgo=", "L-BFGS-B", sep = ""), sep = "_")
}
