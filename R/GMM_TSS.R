##### main function#####
GMMParametersEstim_TSS <-
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
      theta0 <- MoC_TSS(x, c(0.5, 1, 1), eps = eps)
    algo <- match.arg(algo)
    regularization <- match.arg(regularization)
    WeightingMatrix <- match.arg(WeightingMatrix)
    t_scheme <- match.arg(t_scheme)
    t_init <- StableEstim::getTime_()

    method <-
      getGMMmethodName_TSS(
        algo = algo,
        alphaReg = alphaReg,
        regularization = regularization,
        WeightingMatrix = WeightingMatrix,
        t_scheme = t_scheme,
        ...
      )
    Estim <- switch(algo,
                    `2SGMM` = {
                      Compute2SGMMParametersEstim_TSS(
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
                      ComputeITGMMParametersEstim_TSS(
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
                      ComputeCueGMMParametersEstim_TSS(
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
        paste("TSS", "GMMParametersEstim", algo, t_scheme, sep = "_")
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
getGMMmethodName_TSS <-
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

checkIterationControl <- function(IterationControl) {
  NbIter <-
    ifelse(is.null(IterationControl$NbIter),
           10,
           IterationControl$NbIter)
  PrintIter <-
    ifelse(is.null(IterationControl$PrintIter),
           TRUE,
           IterationControl$PrintIter)
  RelativeErrMax <-
    ifelse(
      is.null(IterationControl$RelativeErrMax),
      0.001,
      IterationControl$RelativeErrMax
    )
  list(NbIter = NbIter,
       #Change to TRUE, to see results for each Iteration in console
       PrintIter = FALSE,
       RelativeErrMax = RelativeErrMax)
}

PrintIteration <- function(theta, iter, nbIterMax) {

  #Subordinator
  if(length(theta)==3){
    print(
      paste(
        "---------------iteration ",
        iter,
        "/",
        nbIterMax,
        "------------------------ (",
        theta[1],
        ",",
        theta[2],
        ",",
        theta[3],
        ")"
      )
    )
  }

  #Normal
  if(length(theta)==5){
    print(
      paste(
        "---------------iteration ",
        iter,
        "/",
        nbIterMax,
        "------------------------ (",
        theta[1],
        ",",
        theta[2],
        ",",
        theta[3],
        ",",
        theta[4],
        ",",
        theta[5],
        ")"
      )
    )
  }

  #Classic
  if (length(theta)==6){
    print(
      paste(
        "---------------iteration ",
        iter,
        "/",
        nbIterMax,
        "------------------------ (",
        theta[1],
        ",",
        theta[2],
        ",",
        theta[3],
        ",",
        theta[4],
        ",",
        theta[5],
        ",",
        theta[6],
        ")"
      )
    )
  }

  cat(" \n")
}



##### GMM methods#####

Compute2SGMMParametersEstim_TSS <-
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
      ComputeCurrentEstim_TSS(
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
      ComputeWeightingMatrix_TSS(t, theta1, x, WeightingMatrix, ...)
    CurrentEstimOptInfo <-
      ComputeCurrentEstim_TSS(
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

ComputeITGMMParametersEstim_TSS <-
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
      ComputeCurrentEstim_TSS(
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
            (RelativeErr[3] > RelativeErrMaxArray[3])
           )) {
      ProvidedWeightingMatrix <-
        ComputeWeightingMatrix_TSS(
          t = t,
          theta = PrevEstimParVal,
          x = x,
          WeightingMatrix = WeightingMatrix,
          ...
        )
      AllCurrentEstim <-
        ComputeCurrentEstim_TSS(
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

ComputeCueGMMParametersEstim_TSS <-
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
            (RelativeErr[3] > RelativeErrMaxArray[3])
           )) {
      AllCurrentEstim <-
        ComputeCurrentEstim_TSS(
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
ComputeCurrentEstim_TSS <-
  function(t_scheme,
           theta0,
           x,
           alphaReg,
           regularization,
           WeightingMatrix,
           eps,
           ...) {
    t <-
      ComputeT_TSS(
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
        objective = ComputeObjective_TSS,
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
        lower = c(eps, eps, eps),
        upper = c(1 - eps, Inf, Inf)
      )
    list(OptInfo = optOutput, t = t)
  }


ComputeObjective_TSS <-
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
      ComputeWeightingMatrix_TSS(
        t = t,
        theta = theta,
        x = x,
        WeightingMatrix = WeightingMatrix,
        ...
      )
    gbar <- sampleRealCFMoment_TSS(x = x, t = t, theta = theta)
    K1gbar <-
      ComputeInvKbyG_TSS(
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
sampleRealCFMoment_TSS <- function(x, t, theta) {
  ComplexRes <- sampleComplexCFMoment_TSS(x, t, theta)
  return(c(Re(ComplexRes), Im(ComplexRes)))
}

sampleComplexCFMoment_TSS <- function(x, t, theta) {
  ecf <- function(tt)
    mean(exp(complex(imaginary = tt * x)))
  phiXn <- sapply(t, ecf)
  phiTheta <- ComplexCF_TSS(t, theta)
  return(phiXn - phiTheta)
}

ComplexCF_TSS <- function(t, theta) {
  alpha <- theta[1]
  delta <- theta[2]
  lambda <- theta[3]
  CheckParametersRange_TSS(c(alpha, delta, lambda))
  charTSS(t, alpha, delta, lambda)
}




##### Weighting Matrix related functions#####

ComputeWeightingMatrix_TSS <-
  function(t, theta, x, WeightingMatrix, ...) {
    switch(
      WeightingMatrix,
      OptAsym = {
        K <- asymVarRealCFMoment_TSS(t = t, theta = theta)
      },
      DataVar = {
        K <- DataVarRealCFMoment_TSS(t = t, theta = theta, x = x)
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

asymVarRealCFMoment_TSS <- function(t, theta) {
  m <- length(t)
  res <- matrix(0, ncol = 2 * m, nrow = 2 * m)
  tiPlustj <- crossSum(t, t)
  tiMenostj <- crossSum(t,-t)
  phi <- ComplexCF_TSS(t, theta)
  phi_tiPlus_tj <- sapply(tiPlustj, ComplexCF_TSS, theta)
  phi_tiMenos_tj <- sapply(tiMenostj, ComplexCF_TSS, theta)
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

crossSum <- function(X, Y) {
  fct <- function(x)
    X + x
  res <- sapply(Y, fct)
  res
}

DataVarRealCFMoment_TSS <- function(t, theta, x) {
  gt <-
    DataMatrixRealCFMomentCondition_TSS(t = t, theta = theta, x = x)
  stats::var(gt)
}

DataMatrixRealCFMomentCondition_TSS <- function(t, theta, x) {
  x <- matrix(c(x), ncol = 1)
  x_comp <- x %*% matrix(t, nrow = 1)
  x_comp <- matrix(complex(imaginary = x_comp), ncol = length(t))
  emp_car <- exp(x_comp)
  the_car <- ComplexCF_TSS(t, theta)
  gt <- t(t(emp_car) - the_car)
  gt <- cbind(Re(gt), Im(gt))
  return(gt)
}

##### compute Inverse K#####
ComputeInvKbyG_TSS <-
  function(K, G, alphaReg, regularization, eps) {
    ComputeRegularized <- FALSE
    eigenAnalysis <- getSingularValueDecomposition_TSS(K)
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
      ComputeRegularizedK1G_TSS(
        K = K,
        G = G,
        alphaReg = alphaReg,
        regularization = regularization,
        eps = eps
      )
    K1G
  }

getSingularValueDecomposition_TSS <- function(Kn) {
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

ComputeRegularizedK1G_TSS <-
  function(K, G, alphaReg, regularization, eps) {
    RegularisedSol_TSS(
      Kn = K,
      alphaReg = alphaReg,
      r = G,
      regularization = regularization,
      eps = eps
    )
  }

RegularisedSol_TSS <-
  function(Kn,
           alphaReg,
           r,
           regularization = c("Tikhonov", "LF", "cut-off"),
           eps,
           ...) {
    regularization <- match.arg(regularization)
    return(
      ComputeRegularizedSol_TSS(
        Kn = Kn,
        alpha = alphaReg,
        r = r,
        regularization = regularization,
        eps = eps,
        ...
      )
    )
  }

ComputeRegularizedSol_TSS <-
  function(Kn,
           alpha,
           r,
           regularization = c("Tikhonov", "LF", "cut-off"),
           eps,
           ...) {
    regularization <- match.arg(regularization)
    singularValuesSystem <- getSingularValueDecomposition_TSS(Kn)
    lambda <- singularValuesSystem$lambda
    phi <- singularValuesSystem$phi
    ksi <- singularValuesSystem$ksi
    qAlphaLambda <-
      ComputeqAlphaLambda_TSS(
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

ComputeqAlphaLambda_TSS <-
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

checkValidConstantC <- function(eps, ...) {
  args <- list(...)
  c <- args$c
  K <- args$Kn
  if (is.null(K))
    stop("you need to provide Kn for LF Regularization")
  InvNormK2 <- 1 / (norm(K) ^ 2)
  if ((is.null(c)) || (c < 0) || (c > InvNormK2)) {
    c <- max(eps, InvNormK2 - eps)
  } else
    c <- c
  return(c)
}


##### t-Scheme related functions#####

ComputeT_TSS <- function(tScheme, eps, ...) {
  args <- list(...)
  if (tScheme == "free") {
    checkFreeArgs(args)
    t <- args$t_free
  } else if ((tScheme == "equally") || (tScheme == "NonOptAr")) {
    checkMainArgs(args)
    t <-
      ComputeEquallySpacedPoints_TSS(tScheme = tScheme, eps = eps, ...)
  } else if ((tScheme == "uniformOpt") ||
             (tScheme == "ArithOpt") || (tScheme == "VarOpt")) {
    nb_t <- checkMainArgs(args)
    t <-
      ComputeOptimisedPoints_TSS(tScheme = tScheme, eps = eps, ...)
  }
  t
}

checkFreeArgs <- function(args) {
  if (is.null(args$t_free))
    stop("You need to provide t when you choose free as tScheme")
}

checkMainArgs <- function(args) {
  if (is.null(args$x))
    stop("you didn't provide x to compute t")
  if (is.null(args$nb_t))
    stop("You didn't provide nb_t, a number of points to compute t")
}


ComputeEquallySpacedPoints_TSS <-
  function(..., tScheme, eps, x, nb_t, Constrained = TRUE) {
    Bands <-
      Compute_tBands_TSS(x = x,
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

Compute_tBands_TSS <- function(x, Constrained, eps, ...) {
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

# ComputeFirstRootRealeCF_TSS <- function (x, ..., tol = 0.001, maxIter = 100, lowerBand = 1e-04, upperBand = 30) {
# WelshSol <- WelshFirstRootRealeCF_TSS(x, tol, maxIter) if (WelshSol$phinR < tol) return(WelshSol$t) else
# return(numFirstRootRealeCF_TSS(x, tol, lowerBand, upperBand, ...)$t) } WelshFirstRootRealeCF_TSS <- function (x,
# tol = 0.001, maxIter = 100){ A = 0 iter = 0 m = mean(abs(x)) #val = phinR(A, x) val = mean(cos(A * x)) # this has
# to be double checked!  while ((abs(val) > tol) && (iter < maxIter)) { A = A + val/m #val = phinR(A, x) val =
# mean(cos(A * x)) iter = iter + 1 } list(t = A, phinR = val) } numFirstRootRealeCF_TSS <- function (x, tol = 0.001,
# lowerBand = 1e-04, upperBand = 30, ...) { t_init <- graphFirstRootRealeCF_TSS(x, tol = tol, lowerBand = lowerBand,
# upperBand = upperBand)$t if (is.na(t_init)) t_init <- upperBand objectiveFct <- function(t) abs(mean(cos(t * x)))
# optInfo <- nlminb(start = t_init, objective = objectiveFct, lower = lowerBand, upper = upperBand) list(t =
# as.numeric(optInfo$par), phinR = optInfo$objective) } graphFirstRootRealeCF_TSS <- function (x, tol = 0.001,
# lowerBand = 1e-04, upperBand = 30) { t_seq <- seq(lowerBand, upperBand, tol) phinR <- function (t, x) mean(cos(t *
# x)) phiVal <- sapply(t_seq, phinR, x = x) t <- t_seq[abs(phiVal) < tol][1] list(t = t, phinR = phinR(t, x)) }

ComputeOptimisedPoints_TSS <-
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
        ComputeApproxOptSpacing_TSS(
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
            fn = ObjectiveFctToMinIn_t_TSS,
            gr = NULL,
            ...,
            method = "Nelder-Mead"
          )$par
      else
        t <- stats::nlminb(start = t0, objective = ObjectiveFctToMinIn_t_TSS, ...)
    }
    t
  }

ComputeApproxOptSpacing_TSS <-
  function(..., tScheme, eps, nb_t, Constrained, An) {
    tau0 <-
      ComputeTau0(
        An = An,
        tScheme = tScheme,
        eps = eps,
        nb_t = nb_t
      )
    if (Constrained) {
      Bands <- Compute_tauBands_TSS(tScheme, eps, nb_t, An)
      tauInfo <-
        stats::nlminb(
          start = tau0,
          objective = ObjectiveFctToMinIn_tau_TSS,
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
          objective = ObjectiveFctToMinIn_tau_TSS,
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

ComputeTau0 <- function(An, tScheme, eps, nb_t) {
  if (tScheme == "uniformOpt")
    tau0 <- An / (2 * nb_t)
  else if (tScheme == "ArithOpt")
    tau0 <- (An - eps) / (nb_t * (2 * nb_t + 1))
}

Compute_tauBands_TSS <- function(tScheme, eps, nb_t, An) {
  if (tScheme == "uniformOpt") {
    lower = eps
    upper = An / nb_t
  } else if (tScheme == "ArithOpt") {
    lower = min(eps / 2, An / (nb_t * (nb_t ^ 2 + 1)))
    upper = (An - eps) / (nb_t * (nb_t + 1))
  }
  Bands <- list(lower = lower, upper = upper)
}


ObjectiveFctToMinIn_tau_TSS <-
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
    ObjectiveFctToMinIn_t_TSS(
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

InvDet <- function(Mat)
  1 / abs(det(Mat))

ObjectiveFctToMinIn_t_TSS <-
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
      GMMasymptoticVarianceEstim_TSS(
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


GMMasymptoticVarianceEstim_TSS <-
  function(...,
           t,
           theta,
           x,
           WeightingMatrix,
           eps,
           alphaReg = 0.01,
           regularization = "Tikhonov") {
    K <-
      ComputeWeightingMatrix_TSS(
        t = t,
        theta = theta,
        x = x,
        WeightingMatrix = WeightingMatrix,
        ...
      )
    B <- jacobianSampleRealCFMoment_TSS(t, theta)
    fct <-
      function(G)
        ComputeInvKbyG_TSS(
          K = K,
          G = G,
          alphaReg = alphaReg,
          regularization = regularization,
          eps = eps
        )
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
  }

jacobianSampleRealCFMoment_TSS <- function(t, theta) {
  jac <- jacobianSampleComplexCFMoment_TSS(t, theta)
  apply(jac, 2, function(x)
    c(Re(x), Im(x)))
}

jacobianSampleComplexCFMoment_TSS <- function(t, theta) {
  -jacobianComplexCF_TSS(t, theta)
}

jacobianComplexCF_TSS <- function(t, theta) {
  ComplexCFtoDeriv <- function(th)
    ComplexCF_TSS(t, th)
  NumDeriv_jacobian_TSS(ComplexCFtoDeriv, theta)
}

NumDeriv_jacobian_TSS <-
  function(fctToDeriv, WhereFctIsEvaluated, ...) {
    numDeriv::jacobian(
      fctToDeriv,
      WhereFctIsEvaluated,
      method = "Richardson",
      method.args = list(),
      ...
    )
  }



