##### main function#####
CgmmParametersEstim_CGMY <- function(x, algo = c("2SCgmm", "ITCgmm", "CueCgmm"),
                                     alphaReg = 0.01, subdivisions = 50,
                                     IntegrationMethod = c("Uniform",
                                                           "Simpson"),
                                     randomIntegrationLaw = c("unif", "norm"),
                                     s_min = 0, s_max = 1, theta0 = NULL,
                                     IterationControl = list(), eps = 1e-06,
                                     PrintTime = FALSE, ...) {
    if (is.null(theta0))
        theta0 <- MoC_CGMY(x, c(1, 1, 1, 1.5), eps = eps)
    algo <- match.arg(algo)
    t_init <- StableEstim::getTime_()
    method <- getCgmmMethodName_CGMY(algo = algo, alphaReg = alphaReg,
                                     subdivisions = subdivisions,
                                     IntegrationMethod = IntegrationMethod,
                                     randomIntegrationLaw =
                                       randomIntegrationLaw,
                                     s_min = s_min, s_max = s_max)
    Estim <- switch(algo, `2SCgmm` = {
        Compute2SCgmmParametersEstim_CGMY(x = x, theta0 = theta0,
                                          alphaReg = alphaReg, eps = eps,
                                          s_min = s_min, s_max = s_max,
                                          IntegrationMethod = IntegrationMethod,
                                          randomIntegrationLaw =
                                            randomIntegrationLaw,
                                          subdivisions = subdivisions, ...)
    }, ITCgmm = {
        ComputeITCgmmParametersEstim_CGMY(x = x, theta0 = theta0,
                                          alphaReg = alphaReg, eps = eps,
                                          s_min = s_min, s_max = s_max,
                                          IntegrationMethod = IntegrationMethod,
                                          randomIntegrationLaw =
                                            randomIntegrationLaw,
                                          subdivisions = subdivisions,
                                          IterationControl = IterationControl,
                                          ...)
    }, CueCgmm = {
        ComputeCueCgmmParametersEstim_CGMY(x = x, theta0 = theta0,
                                           alphaReg = alphaReg, eps = eps,
                                           s_min = s_min, s_max = s_max,
                                           IntegrationMethod =
                                             IntegrationMethod,
                                           randomIntegrationLaw =
                                             randomIntegrationLaw,
                                           subdivisions = subdivisions,
                                           IterationControl = IterationControl,
                                           ...)
    }, stop(paste(algo, " not taken into account for Cgmm procedure")))
    if (PrintTime) {
        CallingFct <- paste("CGMY", "CgmmParametersEstim", algo, sep = "_")
        StableEstim::PrintDuration(
          StableEstim::ComputeDuration(t_init, StableEstim::getTime_()),
          CallingFct)
    }
    list(Estim = Estim, duration = as.numeric(
      StableEstim::ComputeDuration (t_init, StableEstim::getTime_(), TRUE)),
      method = method)
}

##### auxiliaries#####
getCgmmMethodName_CGMY <- function(algo, alphaReg, subdivisions,
                                   IntegrationMethod, randomIntegrationLaw,
                                   s_min, s_max, ...) {
    args <- list(...)
    paste("Cgmm", paste("algo=", algo, sep = ""),
          paste("alphaReg=", alphaReg, sep = ""),
          paste("OptimAlgo=", "nlminb", sep = ""),
          paste("sd=", subdivisions, sep = ""),
          paste("IM=", IntegrationMethod, sep = ""),
          paste("RIL=", randomIntegrationLaw, sep = ""),
          paste("s_min=", s_min, sep = ""),
          paste("s_max=",s_max, sep = ""), sep = "_")
}

##### CGMM methods#####
Compute2SCgmmParametersEstim_CGMY <- function(x, theta0, alphaReg, eps, s_min,
                                              s_max, IntegrationMethod,
                                              randomIntegrationLaw,
                                              subdivisions, ...) {
    dots <- list(...)
    if (is.null(dots$control)) {
        control <- list(abs.tol = 1e-15, rel.tol = 1e-07,
                        x.tol = 1.5e-05, xf.tol = 2.2e-10)
        thetaHat <- as.numeric(stats::nlminb(start = theta0,
                                      objective = ComputeObjectiveCgmm_CGMY,
                                      gradient = NULL, hessian = NULL,
                                      Weighting = "Id", x = x,
                                      alphaReg = alphaReg, thetaHat = NULL,
                                      s_min = s_min, s_max = s_max,
                                      IntegrationMethod = IntegrationMethod,
                                      randomIntegrationLaw =
                                        randomIntegrationLaw,
                                      subdivisions = subdivisions, ...,
                                      control = control,
                                      lower = c(eps, eps, eps, eps),
                                      upper = c(Inf, Inf, Inf, 2 - eps))$par)
    } else {
        thetaHat <- as.numeric(stats::nlminb(start = theta0,
                                      objective = ComputeObjectiveCgmm_CGMY,
                                      gradient = NULL, hessian = NULL,
                                      Weighting = "Id", x = x,
                                      alphaReg = alphaReg, thetaHat = NULL,
                                      s_min = s_min, s_max = s_max,
                                      IntegrationMethod = IntegrationMethod,
                                      randomIntegrationLaw =
                                        randomIntegrationLaw,
                                      subdivisions = subdivisions, ...,
                                      lower = c(eps, eps, eps, eps),
                                      upper = c(Inf, Inf, Inf, 2 - eps))$par)
    }
    Cmat <- ComputeCmat_CGMY(x = x, thetaHat = thetaHat, s_min = s_min,
                             s_max = s_max,
                             IntegrationMethod = IntegrationMethod,
                             randomIntegrationLaw = randomIntegrationLaw,
                             subdivisions = subdivisions)
    if (is.null(dots$control)) {
        control <- list(abs.tol = 1e-15, rel.tol = 1e-07, x.tol = 1.5e-05,
                        xf.tol = 2.2e-10)
        res <- stats::nlminb(start = theta0, objective = ComputeObjectiveCgmm_CGMY,
                      Weighting = "optimal", Cmat = Cmat, x = x,
                      alphaReg = alphaReg, thetaHat = thetaHat, s_min = s_min,
                      s_max = s_max, IntegrationMethod = IntegrationMethod,
                      randomIntegrationLaw = randomIntegrationLaw,
                      subdivisions = subdivisions, ..., control = control,
                      lower = c(eps, eps, eps, eps),
                      upper = c(Inf, Inf, Inf, 2 - eps))
    } else {
        res <- stats::nlminb(start = theta0, objective = ComputeObjectiveCgmm_CGMY,
                      Weighting = "optimal", Cmat = Cmat, x = x,
                      alphaReg = alphaReg, thetaHat = thetaHat, s_min = s_min,
                      s_max = s_max, IntegrationMethod = IntegrationMethod,
                      randomIntegrationLaw = randomIntegrationLaw,
                      subdivisions = subdivisions, ...,
                      lower = c(eps, eps, eps, eps),
                      upper = c(Inf, Inf, Inf, 2 - eps))
    }
    list(par = as.numeric(res$par), all = res)
}


ComputeITCgmmParametersEstim_CGMY <- function(x, theta0, alphaReg, eps, s_min,
                                              s_max, IntegrationMethod,
                                              randomIntegrationLaw,
                                              subdivisions, IterationControl,
                                              ...) {
    iter = 0
    IterationControl <- checkIterationControl(IterationControl)
    theta1 <- as.numeric(stats::nlminb(start = theta0,
                                objective = ComputeObjectiveCgmm_CGMY,
                                Weighting = "Id", x = x, alphaReg = alphaReg,
                                thetaHat = NULL, s_min = s_min, s_max = s_max,
                                IntegrationMethod = IntegrationMethod,
                                randomIntegrationLaw = randomIntegrationLaw,
                                subdivisions = subdivisions, ...,
                                lower = c(eps, eps, eps, eps),
                                upper = c(Inf, Inf, Inf, 2 - eps))$par)
    PrevEstimParVal <- theta1
    RelativeErr = IterationControl$RelativeErrMax + 5
    while ((iter < IterationControl$NbIter) &&
           (RelativeErr > IterationControl$RelativeErrMax)) {
        Cmat <- ComputeCmat_CGMY(x = x, thetaHat = PrevEstimParVal,
                                 s_min = s_min, s_max = s_max,
                                 IntegrationMethod = IntegrationMethod,
                                 randomIntegrationLaw = randomIntegrationLaw,
                                 subdivisions = subdivisions)
        dots <- list(...)
        if (is.null(dots$control)) {
            control <- list(abs.tol = 1e-15, rel.tol = 1e-07, x.tol = 1.5e-05,
                            xf.tol = 2.2e-10)
            CurrentEstimAllInfo <- stats::nlminb(start = PrevEstimParVal,
                                          objective = ComputeObjectiveCgmm_CGMY,
                                          Weighting = "optimal", Cmat = Cmat,
                                          alphaReg = alphaReg, x = x,
                                          thetaHat = PrevEstimParVal,
                                          s_min = s_min, s_max = s_max,
                                          IntegrationMethod = IntegrationMethod,
                                          randomIntegrationLaw =
                                            randomIntegrationLaw,
                                          subdivisions = subdivisions, ...,
                                          control = control,
                                          lower = c(eps, eps, eps, eps),
                                          upper = c(Inf, Inf, Inf, 2 - eps))
        } else {
            CurrentEstimAllInfo <- stats::nlminb(start = PrevEstimParVal,
                                          objective = ComputeObjectiveCgmm_CGMY,
                                          Weighting = "optimal", Cmat = Cmat,
                                          alphaReg = alphaReg, x = x,
                                          thetaHat = PrevEstimParVal,
                                          s_min = s_min, s_max = s_max,
                                          IntegrationMethod = IntegrationMethod,
                                          randomIntegrationLaw =
                                            randomIntegrationLaw,
                                          subdivisions = subdivisions, ...,
                                          lower = c(eps, eps, eps, eps),
                                          upper = c(Inf, Inf, Inf, 2 - eps))
        }
        CurrentEstimParVal <- CurrentEstimAllInfo$par
        if (IterationControl$PrintIter)
            PrintIteration(CurrentEstimParVal, iter, IterationControl$NbIter)
        RelativeErr <- abs(CurrentEstimParVal - PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter = iter + 1
    }
    list(par = as.numeric(CurrentEstimParVal), all = CurrentEstimAllInfo)
}

ComputeCueCgmmParametersEstim_CGMY <- function(x, theta0, alphaReg, eps,
                                               s_min, s_max, IntegrationMethod,
                                               randomIntegrationLaw,
                                               subdivisions, IterationControl,
                                               ...) {
    iter = 0
    IterationControl <- checkIterationControl(IterationControl)
    theta1 <- as.numeric(stats::nlminb(start = theta0,
                                objective = ComputeObjectiveCgmm_CGMY,
                                Weighting = "Id", x = x, alphaReg = alphaReg,
                                thetaHat = NULL, s_min = s_min, s_max = s_max,
                                IntegrationMethod = IntegrationMethod,
                                randomIntegrationLaw = randomIntegrationLaw,
                                subdivisions = subdivisions, ...,
                                lower = c(eps, eps, eps, eps),
                                upper = c(Inf, Inf, Inf, 2 - eps))$par)
    PrevEstimParVal <- theta1
    RelativeErr = IterationControl$RelativeErrMax + 5
    while ((iter < IterationControl$NbIter) &&
           (RelativeErr > IterationControl$RelativeErrMax)) {
        dots <- list(...)
        if (is.null(dots$control)) {
            control <- list(abs.tol = 1e-15, rel.tol = 1e-07, x.tol = 1.5e-05,
                            xf.tol = 2.2e-10)
            CurrentEstimAllInfo <- stats::nlminb(start = PrevEstimParVal,
                                          objective = ComputeObjectiveCgmm_CGMY,
                                          Cmat = NULL, Weighting = "optimal",
                                          alphaReg = alphaReg, x = x,
                                          thetaHat = PrevEstimParVal,
                                          s_min = s_min, s_max = s_max,
                                          IntegrationMethod = IntegrationMethod,
                                          randomIntegrationLaw =
                                            randomIntegrationLaw,
                                          subdivisions = subdivisions, ...,
                                          control = control,
                                          lower = c(eps, eps, eps, eps),
                                          upper = c(Inf, Inf, Inf, 2 - eps))
        } else {
            CurrentEstimAllInfo <- stats::nlminb(start = PrevEstimParVal,
                                          objective = ComputeObjectiveCgmm_CGMY,
                                          Cmat = NULL, Weighting = "optimal",
                                          alphaReg = alphaReg, x = x,
                                          thetaHat = PrevEstimParVal,
                                          s_min = s_min, s_max = s_max,
                                          IntegrationMethod = IntegrationMethod,
                                          randomIntegrationLaw =
                                            randomIntegrationLaw,
                                          subdivisions = subdivisions, ...,
                                          lower = c(eps, eps, eps, eps),
                                          upper = c(Inf, Inf, Inf, 2 - eps))
        }
        CurrentEstimParVal <- as.numeric(CurrentEstimAllInfo$par)
        if (IterationControl$PrintIter)
            PrintIteration(CurrentEstimParVal, iter, IterationControl$NbIter)
        RelativeErr <- abs(CurrentEstimParVal - PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter = iter + 1
    }
    list(par = as.numeric(CurrentEstimParVal), all = CurrentEstimAllInfo)
}

##### Objective related functions#####
ComputeObjectiveCgmm_CGMY <- function(theta, Cmat = NULL, x,
                                      Weighting = c("optimal", "Id"), alphaReg,
                                      thetaHat, s_min, s_max, subdivisions = 50,
                                      IntegrationMethod =
                                        c("Uniform", "Simpson"),
                                      randomIntegrationLaw = c("norm", "unif"),
                                      ...) {
    n <- length(x)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    Weighting <- match.arg(Weighting)
    ObjectiveVal <- ComputeCgmmFcts_CGMY(Fct = "Objective", theta = theta,
                                         Cmat = Cmat, x = x,
                                         Weighting = Weighting,
                                         alphaReg = alphaReg,
                                         thetaHat = thetaHat, s_min = s_min,
                                         s_max = s_max,
                                         subdivisions = subdivisions,
                                         IntegrationMethod = IntegrationMethod,
                                         randomIntegrationLaw =
                                           randomIntegrationLaw,
                                         ...)
    as.numeric(Mod(ObjectiveVal))
}


ComputeCgmmFcts_CGMY <- function(Fct = c("Objective", "Covariance"), theta,
                                 Cmat = NULL, x, Weighting = c("optimal", "Id"),
                                 alphaReg, thetaHat, s_min, s_max,
                                 subdivisions = 50,
                                 IntegrationMethod = c("Uniform", "Simpson"),
                                 randomIntegrationLaw = c("norm", "unif"), ...){
    n <- length(x)
    Fct <- match.arg(Fct)
    IntegrationMethod <- match.arg(IntegrationMethod)
    randomIntegrationLaw <- match.arg(randomIntegrationLaw)
    Weighting <- match.arg(Weighting)
    ghatBarFctOft <- function(t_var, X){
      Conj(sampleComplexCFMoment_CGMY(x = X, t = t_var, theta = theta))
    }
    ghatFctOft <- function(t_var, X) Conj(ghatBarFctOft(t_var, X))
    if (Weighting == "Id") {
      ObjectiveVal <-
        StableEstim::IntegrateRandomVectorsProduct(
          f_fct = ghatFctOft,
          X = x,
          g_fct = ghatBarFctOft,
          Y = x,
          s_min = s_min,
          s_max = s_max,
          subdivisions =
            subdivisions,
          IntegrationMethod =
            IntegrationMethod,
          randomIntegrationLaw =
            randomIntegrationLaw,
          ...
        )
    } else {
        V <- ComputeV_CGMY(Fct = Fct, theta = theta, thetaHat = thetaHat, X = x,
                           IntegrationMethod = IntegrationMethod, s_min = s_min,
                           randomIntegrationLaw = randomIntegrationLaw,
                           s_max = s_max, subdivisions = subdivisions, ...)
        if (is.null(Cmat)) {
            if (Fct == "Covariance")
                thetaToUse <- thetaHat else thetaToUse <- theta
            Cmat <- ComputeCmat_CGMY(x = x, thetaHat = thetaToUse,
                                     s_min = s_min, s_max = s_max,
                                     IntegrationMethod = IntegrationMethod,
                                     randomIntegrationLaw =
                                       randomIntegrationLaw,
                                     subdivisions = subdivisions, ...)/(n - 6)
        } else Cmat <- Cmat/(n - 6)
        In <- diag(nrow = n, ncol = n)
        Cmat2 <- Cmat %*% Cmat
        matrixToInverse <- alphaReg * In + Cmat2
        ObjectiveVal <- crossprod(Conj(V), solve(a = matrixToInverse, b = V))
    }
    ObjectiveVal
}

ComputeV_CGMY <- function(Fct = c("Objective", "Covariance"), theta, thetaHat,
                          X, s_min, s_max, IntegrationMethod,
                          randomIntegrationLaw, subdivisions, ...) {
    Fct <- match.arg(Fct)
    g_hat_fct <- function(s, x) {
        sampleComplexCFMoment_CGMY(x = x, t = s, theta = theta)
    }
    g_bar_fct <- function(s, x) {
        Conj(sapply(X = x, FUN = sampleComplexCFMoment_CGMY, t = s,
                    theta = thetaHat))
    }
    Jac_g_hat_fct <- function(s, x) {
        jacobianSampleComplexCFMoment_CGMY(t = s, theta = theta)
    }
    if (Fct == "Covariance") {
      res <-
        StableEstim::IntegrateRandomVectorsProduct(
          f_fct = g_bar_fct,
          X = X,
          g_fct = Jac_g_hat_fct,
          Y = X,
          s_min = s_min,
          s_max = s_max,
          subdivisions = subdivisions,
          IntegrationMethod =
            IntegrationMethod,
          randomIntegrationLaw =
            randomIntegrationLaw,
          ...
        )
    } else if (Fct == "Objective") {
      res <-
        StableEstim::IntegrateRandomVectorsProduct(
          f_fct = g_bar_fct,
          X = X,
          g_fct = g_hat_fct,
          Y = X,
          s_min = s_min,
          s_max = s_max,
          subdivisions = subdivisions,
          IntegrationMethod =
            IntegrationMethod,
          randomIntegrationLaw =
            randomIntegrationLaw,
          ...
        )
    }
    res
}


ComputeCmat_CGMY <- function(x, thetaHat, s_min, s_max, IntegrationMethod,
                             randomIntegrationLaw, subdivisions, ...) {
    f_fct <- function(s, x) {
        sapply(X = x, FUN = sampleComplexCFMoment_CGMY, t = s, theta = thetaHat)
    }
    f_bar_fct <- function(s, x) {
        Conj(f_fct(s, x))
    }
    StableEstim::IntegrateRandomVectorsProduct(
      f_fct = f_bar_fct,
      X = x,
      g_fct = f_fct,
      Y = x,
      s_min = s_min,
      s_max = s_max,
      subdivisions = subdivisions,
      IntegrationMethod = IntegrationMethod,
      randomIntegrationLaw = randomIntegrationLaw,
      ...
    )
}

