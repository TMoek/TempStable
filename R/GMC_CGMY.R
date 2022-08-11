##### main function#####
GMCParametersEstim_CGMY <- function(x, algo = c("2SGMC", "ITGMC", "CueGMC"),
                                    ncond, alphaReg = 0.01,
                                    regularization =
                                      c("Tikhonov", "LF", "cut-off"),
                                    WeightingMatrix = c("OptAsym", "Id"),
                                    theta0 = NULL, IterationControl = list(),
                                    eps = 1e-06, PrintTime = FALSE, ...) {
    if (ncond < 4) {
        stop("You need at least 4 conditions")
    } else {
        ncondfl <- floor(ncond)
    }
    if (is.null(theta0))
        theta0 <- MoC_CGMY(x, c(1, 1, 1, 1.5), eps = eps)
    algo <- match.arg(algo)
    regularization <- match.arg(regularization)
    WeightingMatrix <- match.arg(WeightingMatrix)
    t_init <- StableEstim::getTime_()
    method <- getGMCmethodName_CGMY(algo = algo, ncond = ncondfl,
                                    alphaReg = alphaReg,
                                    regularization = regularization,
                                    WeightingMatrix = WeightingMatrix, ...)
    Estim <- switch(algo, `2SGMC` = {
        Compute2SGMCParametersEstim_CGMY(x = x, ncond = ncondfl,
                                         theta0 = theta0, alphaReg = alphaReg,
                                         regularization = regularization,
                                         WeightingMatrix = WeightingMatrix,
                                         eps = eps, ...)
    }, ITGMM = {
        ComputeITGMCParametersEstim_CGMY(x = x, ncond = ncondfl,
                                         theta0 = theta0, alphaReg = alphaReg,
                                         regularization = regularization,
                                         WeightingMatrix = WeightingMatrix,
                                         IterationControl = IterationControl,
                                         eps = eps, ...)
    }, CueGMM = {
        ComputeCueGMCParametersEstim_CGMY(x = x, ncond = ncondfl,
                                          theta0 = theta0, alphaReg = alphaReg,
                                          regularization = regularization,
                                          WeightingMatrix = WeightingMatrix,
                                          IterationControl = IterationControl,
                                          eps = eps, ...)
    }, stop(paste(algo, " not taken into account !")))
    if (PrintTime) {
        CallingFct <- paste("CGMY", "GMCParametersEstim", algo, "ncond=",
                            ncondfl, sep = "_")
        StableEstim::PrintDuration(
          StableEstim::ComputeDuration(t_init, StableEstim::getTime_()),
          CallingFct)
    }
    list(Estim = Estim$Estim,
         duration = StableEstim::ComputeDuration(t_init,
                                                 StableEstim::getTime_(), TRUE),
         method = method,
         ncond = ncondfl)
}

getGMCmethodName_CGMY <- function(algo, ncond, alphaReg, regularization,
                                  WeightingMatrix, ...) {
    args <- list(...)
    paste(algo,
          paste("ncond=", ncond, sep = ""),
          paste("alphaReg=", alphaReg, sep = ""),
          paste("regularization=", regularization, sep = ""),
          paste("WeightingMatrix=", WeightingMatrix, sep = ""),
          paste("OptimAlgo=", "nlminb", sep = ""), sep = "_")
}


##### GMC methods#####
Compute2SGMCParametersEstim_CGMY <- function(x, ncond, theta0, alphaReg,
                                             regularization, WeightingMatrix,
                                             eps, ...) {
    iter = 0
    AllCurrentEstim <- ComputeCurrentGMC_CGMY(theta0 = theta0, x = x,
                                              ncond = ncond,
                                              alphaReg = alphaReg,
                                              regularization = regularization,
                                              WeightingMatrix = "Id",
                                              eps = eps, ...)
    theta1 <- (AllCurrentEstim$OptInfo)$par
    ProvidedWeightingMatrix <- ComputeGMCWeightingMatrix_CGMY(theta1, x, ncond,
                                                              WeightingMatrix,
                                                              ...)
    CurrentEstimOptInfo <-
      ComputeCurrentGMC_CGMY(theta0 = theta1, x = x, ncond = ncond,
                             alphaReg = alphaReg,
                             regularization = regularization,
                             WeightingMatrix = "provided", eps = eps, ...,
                             ProvidedWeightingMatrix =
                               ProvidedWeightingMatrix)$OptInfo
    list(Estim = CurrentEstimOptInfo)
}


ComputeITGMCParametersEstim_CGMY <- function(x, ncond, theta0, alphaReg,
                                             regularization, WeightingMatrix,
                                             IterationControl, eps, ...) {
    iter = 0
    Control <- checkIterationControl(IterationControl)
    AllCurrentEstim <- ComputeCurrentGMC_CGMY(theta0 = theta0, x = x,
                                              ncond = ncond,
                                              alphaReg = alphaReg,
                                              regularization = regularization,
                                              WeightingMatrix = "Id", eps = eps,
                                              ...)
    theta1 <- (AllCurrentEstim$OptInfo)$par
    PrevEstimParVal <- theta1
    RelativeErr = Control$RelativeErrMax + 5
    while ((iter < Control$NbIter) && (RelativeErr > Control$RelativeErrMax)) {
        ProvidedWeightingMatrix <-
          ComputeWeightingMatrix_CGMY(theta = PrevEstimParVal, x = x,
                                      ncond = ncond,
                                      WeightingMatrix = WeightingMatrix, ...)
        AllCurrentEstim <- ComputeCurrentGMC_CGMY(theta0 = PrevEstimParVal,
                                                  x = x, ncond = ncond,
                                                  alphaReg = alphaReg,
                                                  regularization =
                                                    regularization,
                                                  WeightingMatrix = "provided",
                                                  eps = eps, ...,
                                                  ProvidedWeightingMatrix =
                                                    ProvidedWeightingMatrix)
        CurrentEstimOptInfo <- AllCurrentEstim$OptInfo
        CurrentEstimParVal <- CurrentEstimOptInfo$par
        if (Control$PrintIter)
            PrintIteration(CurrentEstimParVal, iter, Control$NbIter)
        RelativeErr <- abs(CurrentEstimParVal - PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter = iter + 1
    }
    list(Estim = CurrentEstimOptInfo)
}

ComputeCueGMCParametersEstim_CGMY <- function(x, ncond, theta0, alphaReg,
                                              regularization, WeightingMatrix,
                                              IterationControl, eps, ...) {
    iter = 0
    Control <- checkIterationControl(IterationControl)
    PrevEstimParVal <- theta0
    RelativeErr = Control$RelativeErrMax + 5
    while ((iter < Control$NbIter) && (RelativeErr > Control$RelativeErrMax)) {
        AllCurrentEstim <- ComputeCurrentGMC_CGMY(theta0 = PrevEstimParVal,
                                                  x = x, ncond = ncond,
                                                  alphaReg = alphaReg,
                                                  regularization =
                                                    regularization,
                                                  WeightingMatrix =
                                                    WeightingMatrix, eps = eps,
                                                  ...)
        CurrentEstimOptInfo <- AllCurrentEstim$OptInfo
        CurrentEstimParVal <- CurrentEstimOptInfo$par
        if (Control$PrintIter)
            PrintIteration(CurrentEstimParVal, iter, Control$NbIter)
        RelativeErr <- abs(CurrentEstimParVal - PrevEstimParVal)
        PrevEstimParVal <- CurrentEstimParVal
        iter = iter + 1
    }
    list(Estim = CurrentEstimOptInfo)
}

##### current step####
ComputeCurrentGMC_CGMY <- function(theta0, x, ncond, alphaReg, regularization,
                                   WeightingMatrix, eps, ...) {
    optOutput <- stats::nlminb(start = theta0,
                               objective = ComputeGMCObjective_CGMY,
                        gradient = NULL, hessian = NULL, x = x, ncond = ncond,
                        alphaReg = alphaReg, regularization = regularization,
                        WeightingMatrix = WeightingMatrix, eps = eps, ...,
                        lower = c(eps, eps, eps, eps),
                        upper = c(Inf, Inf, Inf, 2 - eps))
    list(OptInfo = optOutput)
}


ComputeGMCObjective_CGMY <- function(theta, x, ncond, alphaReg, regularization,
                                     WeightingMatrix, eps, ...) {
    W <- ComputeGMCWeightingMatrix_CGMY(theta = theta, x = x, ncond = ncond,
                                        WeightingMatrix = WeightingMatrix, ...)
    gbar <- sampleCumulants_CGMY(x = x, ncond = ncond, theta = theta)
    W1gbar <- ComputeInvKbyG_CGMY(K = W, G = gbar, alphaReg = alphaReg,
                                  regularization = regularization, eps = eps)
    obj <- crossprod(gbar, W1gbar)
    as.numeric(obj)
}



##### sample moment condition functions#####
sampleCumulants_CGMY <- function(x, ncond, theta) {
    sampleCumulants <- empiricalCumulants_CGMY(x = x, ncond = ncond)
    -theoreticalCumulants_CGMY(ncond = ncond, theta = theta)
    return(sampleCumulants)
}

empiricalCumulants_CGMY <- function(x, ncond) {
    CumFinder_CGMY(x = x, jmax = ncond)
}

theoreticalCumulants_CGMY <- function(ncond, theta) {
    C <- theta[1]
    G <- theta[2]
    M <- theta[3]
    Y <- theta[4]
    CheckParametersRange_CGMY(c(C, G, M, Y))
    theoretical <- numeric(ncond)
    X <- 2:(ncond + 1)
    theoretical <- sapply(X = X, jththeoreticalCumulant_CGMY, theta = theta)
    return(theoretical)
}

jththeoreticalCumulant_CGMY <- function(j, theta) {
    gamma(j - theta[4]) * (theta[1]/theta[2]^(j - theta[4]) + (-1)^j *
                             theta[1]/theta[3]^(j - theta[4]))
}


##### Weighting Matrix related functions#####

ComputeGMCWeightingMatrix_CGMY <- function(theta, x, ncond, WeightingMatrix,
                                           ...) {
    switch(WeightingMatrix, OptAsym = {
        W <- asymVarCumulants_CGMY(theta = theta, x = x, ncond = ncond)
    }, Id = {
        W <- diag(nrow = ncond, ncol = ncond)
    }, provided = {
        if (!is.null(Wt <- list(...)$ProvidedWeightingMatrix)) W <- Wt
        else stop("You need to provide a Weighting matrix")
    })
    W
}


asymVarCumulants_CGMY <- function(theta, x, ncond) {
    k.pow <- (1:(ncond - 1))
    cumus <- CumFinder_CGMY(x = x, jmax = ncond)[1:(ncond - 1)]
    g <- function(n, z) {
        z^n - ifelse(n > 1, 1, 0) *
        sum(choose(n - 1, (0:(n - 2))) * cumus[1:(n - 1)] * z^(n - (1:(n - 1))))
    }
    gvec <- function(z) {
        sapply(X = (1:ncond), FUN = g, z = z) -
        theoreticalCumulants_CGMY(ncond = ncond, theta = theta)
    }
    gmat <- function(z) {
        matrix(gvec(z)) %*% gvec(z)
    }
    gmatlist <- sapply(X = x, FUN = gmat, simplify = FALSE)
    W <- Reduce("+", gmatlist)/length(x)
    return(W)
}

