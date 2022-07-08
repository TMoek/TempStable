##### main function#####
GMMParametersEstim_CTS <- function(x, algo = c("2SGMM", "ITGMM", "CueGMM"), alphaReg = 0.01, regularization = c("Tikhonov",
    "LF", "cut-off"), WeightingMatrix = c("OptAsym", "DataVar", "Id"), t_scheme = c("equally", "NonOptAr", "uniformOpt",
    "ArithOpt", "VarOpt", "free"), theta0 = NULL, IterationControl = list(), eps = 1e-06, PrintTime = FALSE, ...) {
    if (is.null(theta0))
        theta0 <- MoC_CTS(x, c(1.5, 1, 1, 1, 1, 0), eps = eps)
    algo <- match.arg(algo)
    regularization <- match.arg(regularization)
    WeightingMatrix <- match.arg(WeightingMatrix)
    t_scheme <- match.arg(t_scheme)
    t_init <- getTime_()
    method <- getGMMmethodName_CTS(algo = algo, alphaReg = alphaReg, regularization = regularization, WeightingMatrix = WeightingMatrix,
        t_scheme = t_scheme, ...)
    Estim <- switch(algo, `2SGMM` = {
        Compute2SGMMParametersEstim_CTS(x = x, theta0 = theta0, alphaReg = alphaReg, regularization = regularization, WeightingMatrix = WeightingMatrix,
            t_scheme = t_scheme, eps = eps, ...)
    }, ITGMM = {
        ComputeITGMMParametersEstim_CTS(x = x, theta0 = theta0, alphaReg = alphaReg, regularization = regularization, WeightingMatrix = WeightingMatrix,
            t_scheme = t_scheme, IterationControl = IterationControl, eps = eps, ...)
    }, CueGMM = {
        ComputeCueGMMParametersEstim_CTS(x = x, theta0 = theta0, alphaReg = alphaReg, regularization = regularization,
            WeightingMatrix = WeightingMatrix, t_scheme = t_scheme, IterationControl = IterationControl, eps = eps, ...)
    }, stop(paste(algo, " not taken into account !")))
    if (PrintTime) {
        CallingFct <- paste("Classic", "GMMParametersEstim", algo, t_scheme, sep = "_")
        PrintDuration(ComputeDuration(t_init, getTime_()), CallingFct)
    }
    list(Estim = Estim$Estim, duration = ComputeDuration(t_init, getTime_(), TRUE), method = method, tEstim = Estim$tEstim)
}

##### auxiliaries#####
getGMMmethodName_CTS <- function(algo, alphaReg, regularization, WeightingMatrix, t_scheme, ...) {
    args <- list(...)
    if (!is.null(args$nb_t))
        nt <- args$nb_t else if (!is.null(args$t_free))
        nt <- length(args$t_free)
    paste(algo, paste("nb_t=", nt, sep = ""), paste("alphaReg=", alphaReg, sep = ""), paste("regularization=", regularization,
        sep = ""), paste("WeightingMatrix=", WeightingMatrix, sep = ""), paste("t_scheme=", t_scheme, sep = ""), paste("OptimAlgo=",
        "nlminb", sep = ""), sep = "_")
}



##### GMM methods#####

Compute2SGMMParametersEstim_CTS <- function(x, theta0, alphaReg, regularization, WeightingMatrix, t_scheme, eps, ...) {
    iter = 0
    AllCurrentEstim <- ComputeCurrentEstim_CTS(t_scheme = t_scheme, theta0 = theta0, x = x, alphaReg = alphaReg, regularization = regularization,
        WeightingMatrix = "Id", eps = eps, ...)
    theta1 <- (AllCurrentEstim$OptInfo)$par
    t <- AllCurrentEstim$t
    ProvidedWeightingMatrix <- ComputeWeightingMatrix_CTS(t, theta1, x, WeightingMatrix, ...)
    CurrentEstimOptInfo <- ComputeCurrentEstim_CTS(t_scheme = t_scheme, theta0 = theta1, x = x, alphaReg = alphaReg, regularization = regularization,
        WeightingMatrix = "provided", eps = eps, ..., ProvidedWeightingMatrix = ProvidedWeightingMatrix)$OptInfo
    list(Estim = CurrentEstimOptInfo, tEstim = t)
}

ComputeITGMMParametersEstim_CTS <- function(x, theta0, alphaReg, regularization, WeightingMatrix, t_scheme, IterationControl,
    eps, ...) {
    iter = 0
    Control <- checkIterationControl(IterationControl)
    AllCurrentEstim <- ComputeCurrentEstim_CTS(t_scheme = t_scheme, theta0 = theta0, x = x, alphaReg = alphaReg, regularization = regularization,
        WeightingMatrix = "Id", eps = eps, ...)
    theta1 <- (AllCurrentEstim$OptInfo)$par
    t <- AllCurrentEstim$t
    PrevEstimParVal <- theta1
    RelativeErr = Control$RelativeErrMax + 5
    while ((iter < Control$NbIter) && (RelativeErr > Control$RelativeErrMax)) {
        ProvidedWeightingMatrix <- ComputeWeightingMatrix_CTS(t = t, theta = PrevEstimParVal, x = x, WeightingMatrix = WeightingMatrix,
            ...)
        AllCurrentEstim <- ComputeCurrentEstim_CTS(t_scheme = t_scheme, theta0 = PrevEstimParVal, x = x, alphaReg = alphaReg,
            regularization = regularization, WeightingMatrix = "provided", eps = eps, ..., ProvidedWeightingMatrix = ProvidedWeightingMatrix)
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

ComputeCueGMMParametersEstim_CTS <- function(x, theta0, alphaReg, regularization, WeightingMatrix, t_scheme, IterationControl,
    eps, ...) {
    iter = 0
    Control <- checkIterationControl(IterationControl)
    PrevEstimParVal <- theta0
    RelativeErr = Control$RelativeErrMax + 5
    while ((iter < Control$NbIter) && (RelativeErr > Control$RelativeErrMax)) {
        AllCurrentEstim <- ComputeCurrentEstim_CTS(t_scheme = t_scheme, theta0 = PrevEstimParVal, x = x, alphaReg = alphaReg,
            regularization = regularization, WeightingMatrix = WeightingMatrix, eps = eps, ...)
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
ComputeCurrentEstim_CTS <- function(t_scheme, theta0, x, alphaReg, regularization, WeightingMatrix, eps, ...) {
    t <- ComputeT_CTS(tScheme = t_scheme, theta = theta0, x = x, alphaReg = alphaReg, regularization = regularization,
        WeightingMatrix = WeightingMatrix, eps = eps, ...)
    optOutput <- nlminb(start = theta0, objective = ComputeObjective_CTS, gradient = NULL, hessian = NULL, t = t, x = x,
        alphaReg = alphaReg, regularization = regularization, WeightingMatrix = WeightingMatrix, t_scheme = t_scheme, eps = eps,
        ..., lower = c(eps, eps, eps, eps, eps, -Inf), upper = c(2 - eps, Inf, Inf, Inf, Inf, Inf))
    list(OptInfo = optOutput, t = t)
}


ComputeObjective_CTS <- function(theta, t, x, alphaReg, regularization, WeightingMatrix, t_scheme, eps, ...) {
    K <- ComputeWeightingMatrix_CTS(t = t, theta = theta, x = x, WeightingMatrix = WeightingMatrix, ...)
    gbar <- sampleRealCFMoment_CTS(x = x, t = t, theta = theta)
    K1gbar <- ComputeInvKbyG_CTS(K = K, G = gbar, alphaReg = alphaReg, regularization = regularization, eps = eps)
    obj <- crossprod(gbar, K1gbar)
    as.numeric(obj)
}


##### Sample ECF related functions#####
sampleRealCFMoment_CTS <- function(x, t, theta) {
    ComplexRes <- sampleComplexCFMoment_CTS(x, t, theta)
    return(c(Re(ComplexRes), Im(ComplexRes)))
}

sampleComplexCFMoment_CTS <- function(x, t, theta) {
    ecf <- function(tt) mean(exp(complex(imaginary = tt * x)))
    phiXn <- sapply(t, ecf)
    phiTheta <- ComplexCF_CTS(t, theta)
    return(phiXn - phiTheta)
}

ComplexCF_CTS <- function(t, theta) {
    alpha <- theta[1]
    deltap <- theta[2]
    deltam <- theta[3]
    lambdap <- theta[4]
    lambdam <- theta[5]
    mu <- theta[6]
    CheckParametersRange_CTS(c(alpha, deltap, deltam, lambdap, lambdam, mu))
    charCTS(t, alpha, deltap, deltam, lambdap, lambdam, mu)
}



##### Weighting Matrix related functions#####

ComputeWeightingMatrix_CTS <- function(t, theta, x, WeightingMatrix, ...) {
    switch(WeightingMatrix, OptAsym = {
        K <- asymVarRealCFMoment_CTS(t = t, theta = theta)
    }, DataVar = {
        K <- DataVarRealCFMoment_CTS(t = t, theta = theta, x = x)
    }, Id = {
        K <- diag(nrow = 2 * length(t), ncol = 2 * length(t))
    }, provided = {
        if (!is.null(Kt <- list(...)$ProvidedWeightingMatrix)) K <- Kt else stop("You need to provide a Weighting matrix")
    })
    K
}

asymVarRealCFMoment_CTS <- function(t, theta) {
    m <- length(t)
    res <- matrix(0, ncol = 2 * m, nrow = 2 * m)
    tiPlustj <- crossSum(t, t)
    tiMenostj <- crossSum(t, -t)
    phi <- ComplexCF_CTS(t, theta)
    phi_tiPlus_tj <- sapply(tiPlustj, ComplexCF_CTS, theta)
    phi_tiMenos_tj <- sapply(tiMenostj, ComplexCF_CTS, theta)
    res[1:m, 1:m] <- (0.5 * (Re(phi_tiPlus_tj) + Re(phi_tiMenos_tj)) - Re(phi) %*% t(Re(phi)))
    res[1:m, (m + 1):(2 * m)] <- (0.5 * (Im(phi_tiPlus_tj) - Im(phi_tiMenos_tj)) - Re(phi) %*% t(Im(phi)))
    res[(m + 1):(2 * m), 1:m] <- (0.5 * (Im(phi_tiPlus_tj) + Im(phi_tiMenos_tj)) - Im(phi) %*% t(Re(phi)))
    res[(m + 1):(2 * m), (m + 1):(2 * m)] <- (-0.5 * (Re(phi_tiPlus_tj) - Re(phi_tiMenos_tj)) - Im(phi) %*% t(Im(phi)))
    res
}


DataVarRealCFMoment_CTS <- function(t, theta, x) {
    gt <- DataMatrixRealCFMomentCondition_CTS(t = t, theta = theta, x = x)
    var(gt)
}

DataMatrixRealCFMomentCondition_CTS <- function(t, theta, x) {
    x <- matrix(c(x), ncol = 1)
    x_comp <- x %*% matrix(t, nrow = 1)
    x_comp <- matrix(complex(imaginary = x_comp), ncol = length(t))
    emp_car <- exp(x_comp)
    the_car <- ComplexCF_CTS(t, theta)
    gt <- t(t(emp_car) - the_car)
    gt <- cbind(Re(gt), Im(gt))
    return(gt)
}

##### compute Inverse K#####
ComputeInvKbyG_CTS <- function(K, G, alphaReg, regularization, eps) {
    ComputeRegularized <- FALSE
    eigenAnalysis <- getSingularValueDecomposition_CTS(K)
    if (any(abs(eigenAnalysis$lambda) < alphaReg))
        ComputeRegularized <- TRUE else {
        errStatut <- tryCatch(solve(K, G), error = function(e) e)
        err <- inherits(errStatut, "error")
        if (!err)
            K1G <- errStatut else ComputeRegularized <- TRUE
    }
    if (ComputeRegularized)
        K1G <- ComputeRegularizedK1G_CTS(K = K, G = G, alphaReg = alphaReg, regularization = regularization, eps = eps)
    K1G
}

getSingularValueDecomposition_CTS <- function(Kn) {
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
    return(list(lambda = lambda, phi = phi, ksi = ksi))
}

ComputeRegularizedK1G_CTS <- function(K, G, alphaReg, regularization, eps) {
    RegularisedSol_CTS(Kn = K, alphaReg = alphaReg, r = G, regularization = regularization, eps = eps)
}

RegularisedSol_CTS <- function(Kn, alphaReg, r, regularization = c("Tikhonov", "LF", "cut-off"), eps, ...) {
    regularization <- match.arg(regularization)
    return(ComputeRegularizedSol_CTS(Kn = Kn, alpha = alphaReg, r = r, regularization = regularization, eps = eps, ...))
}

ComputeRegularizedSol_CTS <- function(Kn, alpha, r, regularization = c("Tikhonov", "LF", "cut-off"), eps, ...) {
    regularization <- match.arg(regularization)
    singularValuesSystem <- getSingularValueDecomposition_CTS(Kn)
    lambda <- singularValuesSystem$lambda
    phi <- singularValuesSystem$phi
    ksi <- singularValuesSystem$ksi
    qAlphaLambda <- ComputeqAlphaLambda_CTS(lambda = lambda, alpha = alpha, regularization = regularization, eps = eps,
        ..., Kn = Kn)
    SP_rbyKsi <- crossprod(r, ksi)
    qAlphaSP <- t(qAlphaLambda/lambda) * SP_rbyKsi
    return(rowSums(matrix(data = qAlphaSP, ncol = ncol(phi), nrow = nrow(phi), byrow = TRUE) * phi))
}

ComputeqAlphaLambda_CTS <- function(lambda, alpha, regularization, eps, ...) {
    switch(regularization, Tikhonov = {
        lambda2 <- lambda^2
        return(lambda2/(lambda2 + alpha))
    }, LF = {
        c <- checkValidConstantC(eps = eps, ...)
        invAlpha <- floor((1/alpha))
        return(1 - (1 - c * (lambda^2))^invAlpha)
    }, `cut-off` = {
        sqrtAlpha <- sqrt(alpha)
        TestCutOff <- function(x, a) {
            ifelse(x < sqrt(a), x^2/a, 1)
        }
        return(sapply(lambda, TestCutOff, a = alpha))
    }, stop(paste(regularization, " method not taken into account")))
}

##### t-Scheme related functions#####

ComputeT_CTS <- function(tScheme, eps, ...) {
    args <- list(...)
    if (tScheme == "free") {
        checkFreeArgs(args)
        t <- args$t_free
    } else if ((tScheme == "equally") || (tScheme == "NonOptAr")) {
        checkMainArgs(args)
        t <- ComputeEquallySpacedPoints_CTS(tScheme = tScheme, eps = eps, ...)
    } else if ((tScheme == "uniformOpt") || (tScheme == "ArithOpt") || (tScheme == "VarOpt")) {
        nb_t <- checkMainArgs(args)
        t <- ComputeOptimisedPoints_CTS(tScheme = tScheme, eps = eps, ...)
    }
    t
}



ComputeEquallySpacedPoints_CTS <- function(..., tScheme, eps, x, nb_t, Constrained = TRUE) {
    Bands <- Compute_tBands_CTS(x = x, Constrained = Constrained, eps = eps, ...)
    if (tScheme == "equally") {
        t <- seq(Bands$lower, Bands$upper, length.out = nb_t)
    } else if (tScheme == "NonOptAr") {
        t <- (((1:nb_t) * (2:(nb_t + 1)))/(nb_t * (nb_t + 1))) * Bands$upper
    }
    t
}

Compute_tBands_CTS <- function(x, Constrained, eps, ...) {
    args <- list(...)
    if (!Constrained) {
        lower <- ifelse(is.null(args$min_t), eps, args$min_t)
        upper <- ifelse(is.null(args$max_t), ComputeFirstRootRealeCF(x = x, ...) - eps, args$max_t)
    } else {
        An <- ComputeFirstRootRealeCF(x = x, ...)
        lower = eps
        upper = An - eps
    }
    list(lower = lower, upper = upper)
}



ComputeOptimisedPoints_CTS <- function(..., tScheme, eps, nb_t = 40, Constrained = TRUE, FastOptim = TRUE) {
    An <- ComputeFirstRootRealeCF(...)
    t0 <- seq(eps, An - eps, length.out = nb_t)
    if ((tScheme == "uniformOpt") || (tScheme == "ArithOpt")) {
        t <- ComputeApproxOptSpacing_CTS(..., tScheme = tScheme, eps = eps, nb_t = nb_t, Constrained = Constrained, An = An)
    } else if (tScheme == "VarOpt") {
        # not done yet!
        if (FastOptim)
            t <- optim(par = t0, fn = ObjectiveFctToMinIn_t_CTS, gr = NULL, ..., method = "Nelder-Mead")$par else t <- nlminb(start = t0, objective = ObjectiveFctToMinIn_t_CTS, ...)
    }
    t
}

ComputeApproxOptSpacing_CTS <- function(..., tScheme, eps, nb_t, Constrained, An) {
    tau0 <- ComputeTau0(An = An, tScheme = tScheme, eps = eps, nb_t = nb_t)
    if (Constrained) {
        Bands <- Compute_tauBands_CTS(tScheme, eps, nb_t, An)
        tauInfo <- nlminb(start = tau0, objective = ObjectiveFctToMinIn_tau_CTS, gradient = NULL, hessian = NULL, tScheme = tScheme,
            eps = eps, nb_t = nb_t, ..., lower = Bands$lower, upper = Bands$upper)
    } else {
        tauInfo <- nlminb(start = tau0, objective = ObjectiveFctToMinIn_tau_CTS, gradient = NULL, hessian = NULL, tScheme = tScheme,
            eps = eps, nb_t = nb_t, ...)
    }
    tau <- tauInfo$par
    if (tScheme == "uniformOpt") {
        t <- (1:nb_t) * tau
    } else if (tScheme == "ArithOpt") {
        t <- ((1:nb_t) * (2:(nb_t + 1))) * tau
    }
    t
}



Compute_tauBands_CTS <- function(tScheme, eps, nb_t, An) {
    if (tScheme == "uniformOpt") {
        lower = eps
        upper = An/nb_t
    } else if (tScheme == "ArithOpt") {
        lower = min(eps/2, An/(nb_t * (nb_t^2 + 1)))
        upper = (An - eps)/(nb_t * (nb_t + 1))
    }
    Bands <- list(lower = lower, upper = upper)
}


ObjectiveFctToMinIn_tau_CTS <- function(tau, ..., tScheme, nb_t, theta, x, WeightingMatrix, eps, alphaReg = 0.01, regularization = "Tikhonov",
    fctToApply = InvDet) {
    if (tScheme == "uniformOpt") {
        t <- (1:nb_t) * tau
    } else if (tScheme == "ArithOpt") {
        t <- ((1:nb_t) * (2:(nb_t + 1))) * tau
    }
    ObjectiveFctToMinIn_t_CTS(t, ..., theta = theta, x = x, WeightingMatrix = WeightingMatrix, eps = eps, alphaReg = alphaReg,
        regularization = regularization, fctToApply = fctToApply)
}


ObjectiveFctToMinIn_t_CTS <- function(t, ..., theta, x, WeightingMatrix, eps, alphaReg = 0.01, regularization = "Tikhonov",
    fctToApply = InvDet) {
    W <- GMMasymptoticVarianceEstim_CTS(..., t = t, theta = theta, x = x, WeightingMatrix = WeightingMatrix, eps = eps,
        alphaReg = alphaReg, regularization = regularization)
    fctToApply(W)
}


GMMasymptoticVarianceEstim_CTS <- function(..., t, theta, x, WeightingMatrix, eps, alphaReg = 0.01, regularization = "Tikhonov") {
    K <- ComputeWeightingMatrix_CTS(t = t, theta = theta, x = x, WeightingMatrix = WeightingMatrix, ...)
    B <- jacobianSampleRealCFMoment_CTS(t, theta)
    fct <- function(G) ComputeInvKbyG_CTS(K = K, G = G, alphaReg = alphaReg, regularization = regularization, eps = eps)
    invKcrossB <- apply(X = B, MARGIN = 2, FUN = fct)
    crossprod(B, invKcrossB)
}

jacobianSampleRealCFMoment_CTS <- function(t, theta) {
    jac <- jacobianSampleComplexCFMoment_CTS(t, theta)
    apply(jac, 2, function(x) c(Re(x), Im(x)))
}

jacobianSampleComplexCFMoment_CTS <- function(t, theta) {
    -jacobianComplexCF_CTS(t, theta)
}

jacobianComplexCF_CTS <- function(t, theta) {
    ComplexCFtoDeriv <- function(th) ComplexCF_CTS(t, th)
    NumDeriv_jacobian_CTS(ComplexCFtoDeriv, theta)
}

NumDeriv_jacobian_CTS <- function(fctToDeriv, WhereFctIsEvaluated, ...) {
    library(numDeriv)
    jacobian(fctToDeriv, WhereFctIsEvaluated, method = "Richardson", method.args = list(), ...)
}

