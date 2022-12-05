##### main function#####
GMMParametersEstim_NTS <-
  function(x,
           algo = c("2SGMM", "ITGMM", "CueGMM"),
           alphaReg = 0.01,
           regularization = c("Tikhonov",
                              "LF", "cut-off"),
           WeightingMatrix = c("OptAsym", "DataVar", "Id"),
           t_scheme = c("equally",
                        "NonOptAr",
                        "uniformOpt",
                        "ArithOpt",
                        "VarOpt",
                        "free"),
           theta0 = NULL,
           IterationControl = list(),
           eps = 1e-06,
           PrintTime = FALSE,
           ...) {
    if (is.null(theta0))
      theta0 <- MoC_NTS(x, c(0.5, 0, 1, 1, 0), eps = eps)
    algo <- match.arg(algo)
    regularization <- match.arg(regularization)
    WeightingMatrix <- match.arg(WeightingMatrix)
    t_scheme <- match.arg(t_scheme)
    t_init <- StableEstim::getTime_()
    method <-
      getGMMmethodName_NTS(
        algo = algo,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = WeightingMatrix,
        t_scheme = t_scheme,
        ...
      )
    Estim <- switch(algo,
                    `2SGMM` = {
                      Compute2SGMMParametersEstim_NTS(
                        x = x,
                        theta0 = theta0,
                        alphaReg = alphaReg,
                        regularization = regularization,
                        WeightingMatrix = WeightingMatrix,
                        t_scheme = t_scheme,
                        eps = eps,
                        ...
                      )
                    },
                    ITGMM = {
                      ComputeITGMMParametersEstim_NTS(
                        x = x,
                        theta0 = theta0,
                        alphaReg = alphaReg,
                        regularization = regularization,
                        WeightingMatrix = WeightingMatrix,
                        t_scheme = t_scheme,
                        IterationControl = IterationControl,
                        eps = eps,
                        ...
                      )
                    },
                    CueGMM = {
                      ComputeCueGMMParametersEstim_NTS(
                        x = x,
                        theta0 = theta0,
                        alphaReg = alphaReg,
                        regularization = regularization,
                        WeightingMatrix = WeightingMatrix,
                        t_scheme = t_scheme,
                        IterationControl = IterationControl,
                        eps = eps,
                        ...
                      )
                    },
                    stop(paste(algo, " not taken into account !")))
    if (PrintTime) {
      CallingFct <-
        paste("Normal", "GMMParametersEstim", algo, t_scheme, sep = "_")
      StableEstim::PrintDuration(
        StableEstim::ComputeDuration(t_init, StableEstim::getTime_()),
        CallingFct)
    }
    list(
      Estim = Estim$Estim,
      duration = StableEstim::ComputeDuration(t_init, StableEstim::getTime_(),
                                              TRUE),
      method = method,
      tEstim = Estim$tEstim
    )
  }

##### auxiliaries#####
getGMMmethodName_NTS <-
  function(algo,
           alphaReg,
           regularization,
           WeightingMatrix,
           t_scheme,
           ...) {
    args <- list(...)
    if (!is.null(args$nb_t))
      nt <- args$nb_t
    else if (!is.null(args$t_free))
      nt <- length(args$t_free)
    paste(
      algo,
      paste("nb_t=", nt, sep = ""),
      paste("alphaReg=", alphaReg, sep = ""),
      paste("regularization=", regularization,
            sep = ""),
      paste("WeightingMatrix=", WeightingMatrix, sep = ""),
      paste("t_scheme=", t_scheme, sep = ""),
      paste("OptimAlgo=",
            "nlminb", sep = ""),
      sep = "_"
    )
  }



##### GMM methods#####

Compute2SGMMParametersEstim_NTS <-
  function(x,
           theta0,
           alphaReg,
           regularization,
           WeightingMatrix,
           t_scheme,
           eps,
           ...) {
    iter = 0
    AllCurrentEstim <-
      ComputeCurrentEstim_NTS(
        t_scheme = t_scheme,
        theta0 = theta0,
        x = x,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = "Id",
        eps = eps,
        ...
      )
    theta1 <- (AllCurrentEstim$OptInfo)$par
    t <- AllCurrentEstim$t
    ProvidedWeightingMatrix <-
      ComputeWeightingMatrix_NTS(t, theta1, x, WeightingMatrix, ...)
    CurrentEstimOptInfo <-
      ComputeCurrentEstim_NTS(
        t_scheme = t_scheme,
        theta0 = theta1,
        x = x,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = "provided",
        eps = eps,
        ...,
        ProvidedWeightingMatrix = ProvidedWeightingMatrix
      )$OptInfo
    list(Estim = CurrentEstimOptInfo, tEstim = t)
  }

ComputeITGMMParametersEstim_NTS <-
  function(x,
           theta0,
           alphaReg,
           regularization,
           WeightingMatrix,
           t_scheme,
           IterationControl,
           eps,
           ...) {
    iter = 0
    Control <- checkIterationControl(IterationControl)
    AllCurrentEstim <-
      ComputeCurrentEstim_NTS(
        t_scheme = t_scheme,
        theta0 = theta0,
        x = x,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = "Id",
        eps = eps,
        ...
      )
    theta1 <- (AllCurrentEstim$OptInfo)$par
    t <- AllCurrentEstim$t
    PrevEstimParVal <- theta1
    RelativeErr <- rep(Control$RelativeErrMax + 5, times = length(theta0))
    RelativeErrMaxArray <- rep(Control$RelativeErrMax, times = length(theta0))
    while ((iter < Control$NbIter) &&
           ((RelativeErr[1] > RelativeErrMaxArray[1]) ||
            (RelativeErr[2] > RelativeErrMaxArray[2]) ||
            (RelativeErr[3] > RelativeErrMaxArray[3]) ||
            (RelativeErr[4] > RelativeErrMaxArray[4]) ||
            (RelativeErr[5] > RelativeErrMaxArray[5])
           )) {
      ProvidedWeightingMatrix <-
        ComputeWeightingMatrix_NTS(
          t = t,
          theta = PrevEstimParVal,
          x = x,
          WeightingMatrix = WeightingMatrix,
          ...
        )
      AllCurrentEstim <-
        ComputeCurrentEstim_NTS(
          t_scheme = t_scheme,
          theta0 = PrevEstimParVal,
          x = x,
          alphaReg = alphaReg,
          regularization = regularization,
          WeightingMatrix = "provided",
          eps = eps,
          ...,
          ProvidedWeightingMatrix = ProvidedWeightingMatrix
        )
      CurrentEstimOptInfo <- AllCurrentEstim$OptInfo
      t <- AllCurrentEstim$t
      CurrentEstimParVal <- CurrentEstimOptInfo$par
      if (Control$PrintIter)
        PrintIteration(CurrentEstimParVal, iter, Control$NbIter)
      RelativeErr <- abs(CurrentEstimParVal - PrevEstimParVal)
      PrevEstimParVal <- CurrentEstimParVal
      iter = iter + 1
    }
    list(Estim = CurrentEstimOptInfo, tEstim = t)
  }

ComputeCueGMMParametersEstim_NTS <-
  function(x,
           theta0,
           alphaReg,
           regularization,
           WeightingMatrix,
           t_scheme,
           IterationControl,
           eps,
           ...) {
    iter = 0
    Control <- checkIterationControl(IterationControl)
    PrevEstimParVal <- theta0
    RelativeErr <- rep(Control$RelativeErrMax + 5, times = length(theta0))
    RelativeErrMaxArray <- rep(Control$RelativeErrMax, times = length(theta0))
    while ((iter < Control$NbIter) &&
           ((RelativeErr[1] > RelativeErrMaxArray[1]) ||
            (RelativeErr[2] > RelativeErrMaxArray[2]) ||
            (RelativeErr[3] > RelativeErrMaxArray[3]) ||
            (RelativeErr[4] > RelativeErrMaxArray[4]) ||
            (RelativeErr[5] > RelativeErrMaxArray[5])
           )) {
      AllCurrentEstim <-
        ComputeCurrentEstim_NTS(
          t_scheme = t_scheme,
          theta0 = PrevEstimParVal,
          x = x,
          alphaReg = alphaReg,
          regularization = regularization,
          WeightingMatrix = WeightingMatrix,
          eps = eps,
          ...
        )
      CurrentEstimOptInfo <- AllCurrentEstim$OptInfo
      t <- AllCurrentEstim$t
      CurrentEstimParVal <- CurrentEstimOptInfo$par
      if (Control$PrintIter)
        PrintIteration(CurrentEstimParVal, iter, Control$NbIter)
      RelativeErr <- abs(CurrentEstimParVal - PrevEstimParVal)
      PrevEstimParVal <- CurrentEstimParVal
      iter = iter + 1
    }
    list(Estim = CurrentEstimOptInfo, tEstim = t)
  }


##### current step####
ComputeCurrentEstim_NTS <-
  function(t_scheme,
           theta0,
           x,
           alphaReg,
           regularization,
           WeightingMatrix,
           eps,
           ...) {
    t <-
      ComputeT_NTS(
        tScheme = t_scheme,
        theta = theta0,
        x = x,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = WeightingMatrix,
        eps = eps,
        ...
      )
    optOutput <-
      stats::nlminb(
        start = theta0,
        objective = ComputeObjective_NTS,
        gradient = NULL,
        hessian = NULL,
        t = t,
        x = x,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = WeightingMatrix,
        t_scheme = t_scheme,
        eps = eps,
        ...,
        lower = c(eps,-Inf, eps, eps,-Inf),
        upper = c(1 - eps, Inf, Inf, Inf, Inf)
      )
    list(OptInfo = optOutput, t = t)
  }


ComputeObjective_NTS <-
  function(theta,
           t,
           x,
           alphaReg,
           regularization,
           WeightingMatrix,
           t_scheme,
           eps,
           ...) {
    K <-
      ComputeWeightingMatrix_NTS(
        t = t,
        theta = theta,
        x = x,
        WeightingMatrix = WeightingMatrix,
        ...
      )
    gbar <- sampleRealCFMoment_NTS(x = x, t = t, theta = theta)
    K1gbar <-
      ComputeInvKbyG_NTS(
        K = K,
        G = gbar,
        alphaReg = alphaReg,
        regularization = regularization,
        eps = eps
      )
    obj <- crossprod(gbar, K1gbar)
    as.numeric(obj)
  }


##### Sample ECF related functions#####
sampleRealCFMoment_NTS <- function(x, t, theta) {
  ComplexRes <- sampleComplexCFMoment_NTS(x, t, theta)
  return(c(Re(ComplexRes), Im(ComplexRes)))
}

sampleComplexCFMoment_NTS <- function(x, t, theta) {
  ecf <- function(tt)
    mean(exp(complex(imaginary = tt * x)))
  phiXn <- sapply(t, ecf)
  phiTheta <- ComplexCF_NTS(t, theta)
  return(phiXn - phiTheta)
}

ComplexCF_NTS <- function(t, theta) {
  alpha <- theta[1]
  beta <- theta[2]
  delta <- theta[3]
  lambda <- theta[4]
  mu <- theta[5]
  CheckParametersRange_NTS(c(alpha, beta, delta, lambda, mu))
  charNTS(t, alpha, beta, delta, lambda, mu)
}



##### Weighting Matrix related functions#####

ComputeWeightingMatrix_NTS <-
  function(t, theta, x, WeightingMatrix, ...) {
    switch(
      WeightingMatrix,
      OptAsym = {
        K <- asymVarRealCFMoment_NTS(t = t, theta = theta)
      },
      DataVar = {
        K <- DataVarRealCFMoment_NTS(t = t, theta = theta, x = x)
      },
      Id = {
        K <- diag(nrow = 2 * length(t), ncol = 2 * length(t))
      },
      provided = {
        if (!is.null(Kt <-
                     list(...)$ProvidedWeightingMatrix))
          K <- Kt
        else
          stop("You need to provide a Weighting matrix")
      }
    )
    K
  }

asymVarRealCFMoment_NTS <- function(t, theta) {
  m <- length(t)
  res <- matrix(0, ncol = 2 * m, nrow = 2 * m)
  tiPlustj <- crossSum(t, t)
  tiMenostj <- crossSum(t,-t)
  phi <- ComplexCF_NTS(t, theta)
  phi_tiPlus_tj <- sapply(tiPlustj, ComplexCF_NTS, theta)
  phi_tiMenos_tj <- sapply(tiMenostj, ComplexCF_NTS, theta)
  res[1:m, 1:m] <-
    (0.5 * (Re(phi_tiPlus_tj) + Re(phi_tiMenos_tj)) - Re(phi) %*% t(Re(phi)))
  res[1:m, (m + 1):(2 * m)] <-
    (0.5 * (Im(phi_tiPlus_tj) - Im(phi_tiMenos_tj)) - Re(phi) %*% t(Im(phi)))
  res[(m + 1):(2 * m), 1:m] <-
    (0.5 * (Im(phi_tiPlus_tj) + Im(phi_tiMenos_tj)) - Im(phi) %*% t(Re(phi)))
  res[(m + 1):(2 * m), (m + 1):(2 * m)] <-
    (-0.5 * (Re(phi_tiPlus_tj) - Re(phi_tiMenos_tj)) - Im(phi) %*% t(Im(phi)))
  res
}


DataVarRealCFMoment_NTS <- function(t, theta, x) {
  gt <-
    DataMatrixRealCFMomentCondition_NTS(t = t, theta = theta, x = x)
  stats::var(gt)
}

DataMatrixRealCFMomentCondition_NTS <- function(t, theta, x) {
  x <- matrix(c(x), ncol = 1)
  x_comp <- x %*% matrix(t, nrow = 1)
  x_comp <- matrix(complex(imaginary = x_comp), ncol = length(t))
  emp_car <- exp(x_comp)
  the_car <- ComplexCF_NTS(t, theta)
  gt <- t(t(emp_car) - the_car)
  gt <- cbind(Re(gt), Im(gt))
  return(gt)
}

##### compute Inverse K#####
ComputeInvKbyG_NTS <-
  function(K, G, alphaReg, regularization, eps) {
    ComputeRegularized <- FALSE
    eigenAnalysis <- getSingularValueDecomposition_NTS(K)
    if (any(abs(eigenAnalysis$lambda) < alphaReg))
      ComputeRegularized <- TRUE
    else {
      errStatut <- tryCatch(
        solve(K, G),
        error = function(e)
          e
      )
      err <- inherits(errStatut, "error")
      if (!err)
        K1G <- errStatut
      else
        ComputeRegularized <- TRUE
    }
    if (ComputeRegularized)
      K1G <-
      ComputeRegularizedK1G_NTS(
        K = K,
        G = G,
        alphaReg = alphaReg,
        regularization = regularization,
        eps = eps
      )
    K1G
  }

getSingularValueDecomposition_NTS <- function(Kn) {
  if (isSymmetric(Kn)) {
    SingularValuesDecomposition <- eigen(x = Kn, symmetric = TRUE)
    phi <- SingularValuesDecomposition$vectors
    ksi <- SingularValuesDecomposition$vectors
    lambda <- SingularValuesDecomposition$values
  } else {
    SingularValuesDecomposition <- svd(Kn)
    phi <- SingularValuesDecomposition$v
    ksi <- SingularValuesDecomposition$u
    lambda <- SingularValuesDecomposition$d
  }
  return(list(
    lambda = lambda,
    phi = phi,
    ksi = ksi
  ))
}

ComputeRegularizedK1G_NTS <-
  function(K, G, alphaReg, regularization, eps) {
    RegularisedSol_NTS(
      Kn = K,
      alphaReg = alphaReg,
      r = G,
      regularization = regularization,
      eps = eps
    )
  }

RegularisedSol_NTS <-
  function(Kn,
           alphaReg,
           r,
           regularization = c("Tikhonov", "LF", "cut-off"),
           eps,
           ...) {
    regularization <- match.arg(regularization)
    return(
      ComputeRegularizedSol_NTS(
        Kn = Kn,
        alpha = alphaReg,
        r = r,
        regularization = regularization,
        eps = eps,
        ...
      )
    )
  }

ComputeRegularizedSol_NTS <-
  function(Kn,
           alpha,
           r,
           regularization = c("Tikhonov", "LF", "cut-off"),
           eps,
           ...) {
    regularization <- match.arg(regularization)
    singularValuesSystem <- getSingularValueDecomposition_NTS(Kn)
    lambda <- singularValuesSystem$lambda
    phi <- singularValuesSystem$phi
    ksi <- singularValuesSystem$ksi
    qAlphaLambda <-
      ComputeqAlphaLambda_NTS(
        lambda = lambda,
        alpha = alpha,
        regularization = regularization,
        eps = eps,
        ...,
        Kn = Kn
      )
    SP_rbyKsi <- crossprod(r, ksi)
    qAlphaSP <- t(qAlphaLambda / lambda) * SP_rbyKsi
    return(rowSums(
      matrix(
        data = qAlphaSP,
        ncol = ncol(phi),
        nrow = nrow(phi),
        byrow = TRUE
      ) * phi
    ))
  }

ComputeqAlphaLambda_NTS <-
  function(lambda, alpha, regularization, eps, ...) {
    switch(
      regularization,
      Tikhonov = {
        lambda2 <- lambda ^ 2
        return(lambda2 / (lambda2 + alpha))
      },
      LF = {
        c <- checkValidConstantC(eps = eps, ...)
        invAlpha <- floor((1 / alpha))
        return(1 - (1 - c * (lambda ^ 2)) ^ invAlpha)
      },
      `cut-off` = {
        sqrtAlpha <- sqrt(alpha)
        TestCutOff <- function(x, a) {
          ifelse(x < sqrt(a), x ^ 2 / a, 1)
        }
        return(sapply(lambda, TestCutOff, a = alpha))
      },
      stop(paste(
        regularization, " method not taken into account"
      ))
    )
  }

##### t-Scheme related functions#####

ComputeT_NTS <- function(tScheme, eps, ...) {
  args <- list(...)
  if (tScheme == "free") {
    checkFreeArgs(args)
    t <- args$t_free
  } else if ((tScheme == "equally") || (tScheme == "NonOptAr")) {
    checkMainArgs(args)
    t <-
      ComputeEquallySpacedPoints_NTS(tScheme = tScheme, eps = eps, ...)
  } else if ((tScheme == "uniformOpt") ||
             (tScheme == "ArithOpt") || (tScheme == "VarOpt")) {
    nb_t <- checkMainArgs(args)
    t <-
      ComputeOptimisedPoints_NTS(tScheme = tScheme, eps = eps, ...)
  }
  t
}



ComputeEquallySpacedPoints_NTS <-
  function(..., tScheme, eps, x, nb_t, Constrained = TRUE) {
    Bands <-
      Compute_tBands_NTS(x = x,
                         Constrained = Constrained,
                         eps = eps,
                         ...)
    if (tScheme == "equally") {
      t <- seq(Bands$lower, Bands$upper, length.out = nb_t)
    } else if (tScheme == "NonOptAr") {
      t <- (((1:nb_t) * (2:(nb_t + 1))) / (nb_t * (nb_t + 1))) * Bands$upper
    }
    t
  }

Compute_tBands_NTS <- function(x, Constrained, eps, ...) {
  args <- list(...)
  if (!Constrained) {
    lower <- ifelse(is.null(args$min_t), eps, args$min_t)
    upper <-
      ifelse(is.null(args$max_t),
             ComputeFirstRootRealeCF(x = x, ...) - eps,
             args$max_t)
  } else {
    An <- ComputeFirstRootRealeCF(x = x, ...)
    lower = eps
    upper = An - eps
  }
  list(lower = lower, upper = upper)
}



ComputeOptimisedPoints_NTS <-
  function(...,
           tScheme,
           eps,
           nb_t = 40,
           Constrained = TRUE,
           FastOptim = TRUE) {
    An <- StableEstim::ComputeFirstRootRealeCF(...)
    t0 <- seq(eps, An - eps, length.out = nb_t)
    if ((tScheme == "uniformOpt") || (tScheme == "ArithOpt")) {
      t <-
        ComputeApproxOptSpacing_NTS(
          ...,
          tScheme = tScheme,
          eps = eps,
          nb_t = nb_t,
          Constrained = Constrained,
          An = An
        )
    } else if (tScheme == "VarOpt") {
      # not done yet!
      if (FastOptim)
        t <-
          stats::optim(
            par = t0,
            fn = ObjectiveFctToMinIn_t_NTS,
            gr = NULL,
            ...,
            method = "Nelder-Mead"
          )$par
      else
        t <- stats::nlminb(start = t0, objective = ObjectiveFctToMinIn_t_NTS, ...)
    }
    t
  }

ComputeApproxOptSpacing_NTS <-
  function(..., tScheme, eps, nb_t, Constrained, An) {
    tau0 <-
      ComputeTau0(
        An = An,
        tScheme = tScheme,
        eps = eps,
        nb_t = nb_t
      )
    if (Constrained) {
      Bands <- Compute_tauBands_NTS(tScheme, eps, nb_t, An)
      tauInfo <-
        stats::nlminb(
          start = tau0,
          objective = ObjectiveFctToMinIn_tau_NTS,
          gradient = NULL,
          hessian = NULL,
          tScheme = tScheme,
          eps = eps,
          nb_t = nb_t,
          ...,
          lower = Bands$lower,
          upper = Bands$upper
        )
    } else {
      tauInfo <-
        stats::nlminb(
          start = tau0,
          objective = ObjectiveFctToMinIn_tau_NTS,
          gradient = NULL,
          hessian = NULL,
          tScheme = tScheme,
          eps = eps,
          nb_t = nb_t,
          ...
        )
    }
    tau <- tauInfo$par
    if (tScheme == "uniformOpt") {
      t <- (1:nb_t) * tau
    } else if (tScheme == "ArithOpt") {
      t <- ((1:nb_t) * (2:(nb_t + 1))) * tau
    }
    t
  }



Compute_tauBands_NTS <- function(tScheme, eps, nb_t, An) {
  if (tScheme == "uniformOpt") {
    lower = eps
    upper = An / nb_t
  } else if (tScheme == "ArithOpt") {
    lower = min(eps / 2, An / (nb_t * (nb_t ^ 2 + 1)))
    upper = (An - eps) / (nb_t * (nb_t + 1))
  }
  Bands <- list(lower = lower, upper = upper)
}


ObjectiveFctToMinIn_tau_NTS <-
  function(tau,
           ...,
           tScheme,
           nb_t,
           theta,
           x,
           WeightingMatrix,
           eps,
           alphaReg = 0.01,
           regularization = "Tikhonov",
           fctToApply = InvDet) {
    if (tScheme == "uniformOpt") {
      t <- (1:nb_t) * tau
    } else if (tScheme == "ArithOpt") {
      t <- ((1:nb_t) * (2:(nb_t + 1))) * tau
    }
    ObjectiveFctToMinIn_t_NTS(
      t,
      ...,
      theta = theta,
      x = x,
      WeightingMatrix = WeightingMatrix,
      eps = eps,
      alphaReg = alphaReg,
      regularization = regularization,
      fctToApply = fctToApply
    )
  }


ObjectiveFctToMinIn_t_NTS <-
  function(t,
           ...,
           theta,
           x,
           WeightingMatrix,
           eps,
           alphaReg = 0.01,
           regularization = "Tikhonov",
           fctToApply = InvDet) {
    W <-
      GMMasymptoticVarianceEstim_NTS(
        ...,
        t = t,
        theta = theta,
        x = x,
        WeightingMatrix = WeightingMatrix,
        eps = eps,
        alphaReg = alphaReg,
        regularization = regularization
      )
    fctToApply(W)
  }


GMMasymptoticVarianceEstim_NTS <-
  function(...,
           t,
           theta,
           x,
           WeightingMatrix,
           eps,
           alphaReg = 0.01,
           regularization = "Tikhonov") {
    K <-
      ComputeWeightingMatrix_NTS(
        t = t,
        theta = theta,
        x = x,
        WeightingMatrix = WeightingMatrix,
        ...
      )
    B <- jacobianSampleRealCFMoment_NTS(t, theta)
    fct <-
      function(G)
        ComputeInvKbyG_NTS(
          K = K,
          G = G,
          alphaReg = alphaReg,
          regularization = regularization,
          eps = eps
        )
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
  }

jacobianSampleRealCFMoment_NTS <- function(t, theta) {
  jac <- jacobianSampleComplexCFMoment_NTS(t, theta)
  apply(jac, 2, function(x)
    c(Re(x), Im(x)))
}

jacobianSampleComplexCFMoment_NTS <- function(t, theta) {
  -jacobianComplexCF_NTS(t, theta)
}

jacobianComplexCF_NTS <- function(t, theta) {
  ComplexCFtoDeriv <- function(th)
    ComplexCF_NTS(t, th)
  NumDeriv_jacobian_NTS(ComplexCFtoDeriv, theta)
}

NumDeriv_jacobian_NTS <-
  function(fctToDeriv, WhereFctIsEvaluated, ...) {
    numDeriv::jacobian(
      fctToDeriv,
      WhereFctIsEvaluated,
      method = "Richardson",
      method.args = list(),
      ...
    )
  }

